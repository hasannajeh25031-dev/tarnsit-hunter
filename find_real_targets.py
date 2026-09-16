"""Find the REAL pixel positions of the target stars in every observation frame.

Unlike radec_to_xy.py (which only predicts a position from the header pointing,
assuming the target sits exactly where the telescope reports), this script:

  1. Detects the stars actually present in each dark-subtracted frame.
  2. Projects the Gaia DR3 catalog of the field through a header-built WCS
     for all four flip conventions (orientation was unverified).
  3. Pattern-matches catalog vs detected stars to find the true orientation
     and the frame's pointing offset (the telescope drifts / mispoints).
  4. Reports the real target position: the detected star that coincides with
     the catalog position of the target after the field is solved.

Outputs:
  results/real_target_positions.csv   per-frame real target X/Y
  results/real_target_summary.csv     per-session summary
  results/real_targets_map.png        one solved frame per target, target circled
"""

import csv
import glob
import os
import warnings

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from photutils.detection import DAOStarFinder
from scipy.spatial.distance import cdist
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

BASE = os.path.dirname(os.path.abspath(__file__))
DB = os.path.join(BASE, "database")
OUTDIR = os.path.join(BASE, "results")
MATCH_TOL = 5.0        # px: catalog-to-detected match tolerance
TARGET_TOL = 6.0       # px: target must have a detected star this close
MIN_INLIERS = 5        # stars that must match for a field to count as solved
MAX_OFFSET = 160.0     # px: max plausible pointing error + drift

CONVENTIONS = [
    ("N-up_E-left", False, False),
    ("flip-x", True, False),
    ("flip-y", False, True),
    ("flip-x_flip-y", True, True),
]


def build_wcs(header, flip_x, flip_y):
    scale = header["IM_SCALE"] / 3600.0
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [header["RA"], header["DEC"]]
    w.wcs.crpix = [(header["NAXIS1"] + 1) / 2.0, (header["NAXIS2"] + 1) / 2.0]
    w.wcs.cdelt = [(1.0 if flip_x else -1.0) * scale,
                   (-1.0 if flip_y else 1.0) * scale]
    w.wcs.equinox = header.get("EQUINOX", 2000.0)
    return w


def detect_stars(img, threshold_sigma=7, max_stars=45):
    mean, median, std = sigma_clipped_stats(img, sigma=3)
    if std <= 0:
        return np.empty((0, 2)), np.empty(0), median, std
    finder = DAOStarFinder(fwhm=3.5, threshold=threshold_sigma * std)
    tbl = finder(img - median)
    if tbl is None or len(tbl) == 0:
        return np.empty((0, 2)), np.empty(0), median, std
    tbl.sort("flux")
    tbl.reverse()
    pos = np.array([tbl["x_centroid"], tbl["y_centroid"]]).T[:max_stars]
    flux = np.array(tbl["flux"])[:max_stars]
    return pos, flux, median, std


def solve_field(detected, cat_xy):
    """Best translation aligning predicted catalog positions onto detections.

    Tries offsets implied by pairing bright catalog stars with bright
    detections; scores each by inlier count within MATCH_TOL.
    Returns (offset, n_inliers) or (None, 0).
    """
    if len(detected) < 3 or len(cat_xy) < 3:
        return None, 0
    best_off, best_n, best_res = None, 0, 1e9
    for q in cat_xy[:15]:
        for p in detected[:15]:
            off = p - q
            if np.hypot(*off) > MAX_OFFSET:
                continue
            d = cdist(cat_xy + off, detected).min(axis=1)
            inl = d < MATCH_TOL
            n = int(inl.sum())
            res = d[inl].mean() if n else 1e9
            if n > best_n or (n == best_n and res < best_res):
                best_off, best_n, best_res = off, n, res
    if best_n >= 3 and best_off is not None:
        # refine: recompute offset from the matched pairs
        d = cdist(cat_xy + best_off, detected)
        j = d.argmin(axis=1)
        inl = d.min(axis=1) < MATCH_TOL
        if inl.sum() >= 3:
            best_off = (detected[j[inl]] - cat_xy[inl]).mean(axis=0)
            d = cdist(cat_xy + best_off, detected).min(axis=1)
            best_n = int((d < MATCH_TOL).sum())
    return best_off, best_n


def load_catalogs():
    coords = {r["folder"]: (float(r["ra_deg"]), float(r["dec_deg"]))
              for r in csv.DictReader(open(os.path.join(BASE, "target_coords.csv")))}
    fields = {}
    for r in csv.DictReader(open(os.path.join(BASE, "field_stars.csv"))):
        fields.setdefault(r["folder"], []).append(
            (float(r["ra_deg"]), float(r["dec_deg"]), float(r["gmag"])))
    return coords, {k: np.array(v) for k, v in fields.items()}


def master_dark(night):
    for d in (os.path.join(DB, "calibration", night),
              os.path.join(BASE, "calibration", night)):
        fs = sorted(glob.glob(os.path.join(d, "*.fits")))
        if fs:
            return np.median([fits.getdata(f).astype(float) for f in fs], axis=0)
    return 0.0


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    coords, fields = load_catalogs()

    sessions = sorted(glob.glob(os.path.join(DB, "observations", "*", "*", "session_*")))
    print(f"{len(sessions)} sessions found\n")

    rows = []
    summary = []
    conv_tally = {name: 0 for name, _, _ in CONVENTIONS}
    best_frame_per_target = {}   # target -> (n_inliers, img, tx, ty, night, file)

    for sess in sessions:
        parts = sess.split(os.sep)
        night, target, sess_id = parts[-3], parts[-2], parts[-1]
        if target not in coords:
            print(f"[skip] unknown target {target}")
            continue
        files = sorted(glob.glob(os.path.join(sess, "*.fits")))
        if not files:
            continue
        dark = master_dark(night)
        t_ra, t_dec = coords[target]
        cat = fields[target]

        n_solved = 0
        sess_conv = {}
        xs, ys = [], []
        for f in files:
            with fits.open(f) as hdul:
                h = hdul[0].header
                img = hdul[0].data.astype(float) - dark
            det, dflux, sky, std = detect_stars(img)
            H, W = img.shape

            best = None  # (n_inliers, conv_name, off, cat_xy, tx_pred, ty_pred)
            for name, fx, fy in CONVENTIONS:
                w = build_wcs(h, fx, fy)
                cxy = np.array(w.wcs_world2pix(cat[:, 0], cat[:, 1], 0)).T
                inside = ((cxy[:, 0] > -40) & (cxy[:, 0] < W + 40) &
                          (cxy[:, 1] > -40) & (cxy[:, 1] < H + 40))
                cxy_in = cxy[inside]
                off, n = solve_field(det, cxy_in)
                if off is not None and (best is None or n > best[0]):
                    tx, ty = w.wcs_world2pix(t_ra, t_dec, 0)
                    best = (n, name, off, float(tx), float(ty))

            row = {
                "night": night, "target": target, "session": sess_id,
                "file": os.path.basename(f),
                "mjd_obs": h["MJD-OBS"], "ut_obs": h["UT-OBS"],
                "weather": h.get("WEATHER", -1), "n_detected": len(det),
                "solved": "no", "convention": "", "n_matched": 0,
                "x_real": "", "y_real": "", "sep_pix": "", "target_detected": "no",
            }
            if best is not None and best[0] >= MIN_INLIERS:
                n, name, off, tx, ty = best
                x_solved, y_solved = tx + off[0], ty + off[1]
                row.update(solved="yes", convention=name, n_matched=n)
                conv_tally[name] += 1
                sess_conv[name] = sess_conv.get(name, 0) + 1
                n_solved += 1
                # is there a real detected star at the solved target position?
                if len(det):
                    d = np.hypot(det[:, 0] - x_solved, det[:, 1] - y_solved)
                    i = d.argmin()
                    if d[i] < TARGET_TOL:
                        row.update(x_real=round(float(det[i, 0]), 2),
                                   y_real=round(float(det[i, 1]), 2),
                                   sep_pix=round(float(d[i]), 2),
                                   target_detected="yes")
                        xs.append(det[i, 0]); ys.append(det[i, 1])
                        prev = best_frame_per_target.get(target)
                        if prev is None or n > prev[0]:
                            best_frame_per_target[target] = (
                                n, img, det[i, 0], det[i, 1], night,
                                os.path.basename(f))
                if row["target_detected"] == "no":
                    # field solved but no star at target loc (faint/clouded)
                    row.update(x_real=round(x_solved, 2), y_real=round(y_solved, 2))
            rows.append(row)

        conv_str = max(sess_conv, key=sess_conv.get) if sess_conv else "-"
        med_x = round(float(np.median(xs)), 1) if xs else ""
        med_y = round(float(np.median(ys)), 1) if ys else ""
        summary.append({
            "night": night, "target": target, "session": sess_id,
            "n_frames": len(files), "n_solved": n_solved,
            "n_target_detected": len(xs), "convention": conv_str,
            "x_median": med_x, "y_median": med_y,
        })
        print(f"{night} {target:9s} {sess_id}: {len(files):3d} frames, "
              f"{n_solved:3d} solved, target detected in {len(xs):3d}, "
              f"conv={conv_str}, median xy=({med_x}, {med_y})", flush=True)

    with open(os.path.join(OUTDIR, "real_target_positions.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader(); w.writerows(rows)
    with open(os.path.join(OUTDIR, "real_target_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=summary[0].keys())
        w.writeheader(); w.writerows(summary)

    print("\nOrientation vote across all solved frames:")
    for k, v in sorted(conv_tally.items(), key=lambda kv: -kv[1]):
        print(f"  {k:15s} {v}")

    # figure: one solved frame per target with the real target circled
    targets = sorted(best_frame_per_target)
    if targets:
        ncol = 4
        nrow = int(np.ceil(len(targets) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.6 * nrow))
        axes = np.atleast_1d(axes).ravel()
        for ax, t in zip(axes, targets):
            n, img, tx, ty, night, fname = best_frame_per_target[t]
            lo, hi = np.percentile(img, [5, 99.7])
            ax.imshow(img, cmap="gray", vmin=lo, vmax=hi, origin="lower")
            c = plt.Circle((tx, ty), 12, ec="red", fc="none", lw=1.5)
            ax.add_patch(c)
            ax.set_title(f"{t}  ({tx:.1f}, {ty:.1f})\n{night} {fname}", fontsize=8)
            ax.set_xticks([]); ax.set_yticks([])
        for ax in axes[len(targets):]:
            ax.axis("off")
        plt.suptitle("Real target star positions (Gaia-matched, dark-subtracted)")
        plt.tight_layout()
        out = os.path.join(OUTDIR, "real_targets_map.png")
        plt.savefig(out, dpi=110)
        print(f"\nFigure: {out}")


if __name__ == "__main__":
    main()
