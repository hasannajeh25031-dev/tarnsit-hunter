"""Measure the correct (centroid-refined) target star position in each frame.

Starting from the Gaia-solved positions of find_real_targets.py, each frame is
dark-subtracted and the star position is refined to sub-pixel accuracy with a
photutils 2D centroid inside a small box. A refinement is accepted only when
the star is significantly above the local background and the centroid stays
close to the solved position; otherwise the frame keeps the solved position
and is flagged detected = no.

Outputs:
  results/star_positions.csv    the per-frame position list
  printed python list           (file, x, y) per solved session
"""

import csv
import os
from collections import defaultdict

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.nddata import Cutout2D
from photutils.centroids import centroid_com

from find_real_targets import DB, OUTDIR, master_dark

BOX = 9          # centroid box size in pixels
MAX_SHIFT = 4.0  # px: max allowed centroid shift from the solved position
MIN_SNR = 4.0    # peak significance needed to call the star detected


def refine(img, x0, y0):
    """Return (x, y, snr, ok) with a background-subtracted COM centroid."""
    cut = Cutout2D(img, (x0, y0), BOX, mode="partial", fill_value=np.nan)
    stamp = cut.data
    # local background from an annulus-like border of a wider box
    wide = Cutout2D(img, (x0, y0), BOX + 12, mode="partial", fill_value=np.nan).data
    border = wide.copy()
    b = 6
    border[b:-b, b:-b] = np.nan
    bkg = np.nanmedian(border)
    std = 1.4826 * np.nanmedian(np.abs(border - bkg))
    sub = np.nan_to_num(stamp - bkg)
    snr = float(np.nanmax(stamp) - bkg) / std if std > 0 else 0.0
    if snr < MIN_SNR:
        return x0, y0, snr, False
    sub[sub < 0] = 0.0
    cx, cy = centroid_com(sub)
    x, y = cut.to_original_position((cx, cy))
    if not np.isfinite(x) or np.hypot(x - x0, y - y0) > MAX_SHIFT:
        return x0, y0, snr, False
    return float(x), float(y), snr, True


def main():
    rows = list(csv.DictReader(open(os.path.join(OUTDIR, "real_target_positions.csv"))))
    out_rows = []
    positions = defaultdict(list)   # (night, target, session) -> [(file, x, y)]

    for r in rows:
        if r["solved"] != "yes" or not r["x_real"]:
            continue
        night, target, sess = r["night"], r["target"], r["session"]
        path = os.path.join(DB, "observations", night, target, sess, r["file"])
        img = fits.getdata(path).astype(float) - master_dark(night)
        x, y, snr, ok = refine(img, float(r["x_real"]), float(r["y_real"]))
        out_rows.append({
            "night": night, "target": target, "session": sess, "file": r["file"],
            "mjd_obs": r["mjd_obs"],
            "x": round(x, 2), "y": round(y, 2),
            "snr": round(snr, 1), "detected": "yes" if ok else "no",
        })
        if ok:
            positions[(night, target, sess)].append((r["file"], round(x, 2), round(y, 2)))

    out = os.path.join(OUTDIR, "star_positions.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_rows[0].keys())
        w.writeheader(); w.writerows(out_rows)

    n_ok = sum(r["detected"] == "yes" for r in out_rows)
    print(f"{len(out_rows)} solved frames, star measured in {n_ok}")
    print(f"Saved: {out}\n")

    for key in sorted(positions):
        night, target, sess = key
        plist = positions[key]
        print(f"# {target} {night} {sess} — {len(plist)} frames")
        print(f"{target.lower().replace('-', '_')}_{night.replace('-', '')} = [")
        for fname, x, y in plist:
            print(f"    ({fname!r}, {x}, {y}),")
        print("]\n")


if __name__ == "__main__":
    main()
