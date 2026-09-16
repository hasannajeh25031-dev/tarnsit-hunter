"""Differential aperture photometry of the target star in every solved session.

No frame is deleted: every frame with a measured target position is
photometered and classified —
  ok        clean frame, used for all statistics
  weather   header WEATHER <= MIN_WEATHER (hazy/cloudy sky report)
  position  frame's solved position is off the night's drift trend
  cloud     comparison-ensemble flux collapsed (cloud crossing the field)
  outlier   >3 sigma from the running median of clean frames
  no_flux   star or comparisons unmeasurable in this frame
All frames appear in the CSV (status column) and in the per-session graph
(distinct markers); only `ok` frames define the detrending, rms and the
binned curve used downstream.

Per-session graphs use the real UT clock time of the night on the x axis and
show every frame individually with its error bar; the predicted transit
window is shaded.

Outputs per session in results/photometry/:
  <session>_lightcurve.csv         every frame, with status
  <session>_lightcurve_binned.csv  10-min binned curve (ok frames)
  <session>_lightcurve.png         detailed graph, real UT time
Plus results/photometry/all_sessions.png (binned overview).
"""

import csv
import os
from collections import defaultdict

import numpy as np
from astropy.io import fits
from astropy.time import Time
from photutils.aperture import (CircularAperture, CircularAnnulus,
                                ApertureStats, aperture_photometry)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from find_real_targets import (DB, OUTDIR, master_dark, build_wcs,
                               load_catalogs, detect_stars)
from check_transits import EPHEMERIDES, MJD_TO_JD

PHOTDIR = os.path.join(OUTDIR, "photometry")
RADII = [2.5, 3.0, 3.5, 4.0, 5.0]
R_IN, R_OUT = 7.0, 11.0
N_COMPS = 8            # max comparison stars
COMP_DMAG = 1.5        # comps within this many mag of the target (Gaia G)
MATCH_R = 6.0          # px: comp must have a real detection this close
MIN_MATCH = 0.85       # comp must be detected in this fraction of clean frames
EDGE = 12              # px: comps must stay this far inside the detector
TRAJ_CLIP = 5.0        # px: flag frames whose position is off the trend
MIN_FRAMES = 5         # skip sessions with fewer measured frames
MIN_WEATHER = 90       # frames at or below this WEATHER value are flagged
CLOUD_CUT = 0.6        # flag frames with ensemble below this frac of median
CLIP_SIGMA = 3.0       # running-median outlier flag
BIN_MIN = 10.0         # bin width in minutes

STATUS_STYLE = {       # marker style per category in the detailed graph
    "weather":  dict(fmt="^", color="darkorange", label="weather flagged"),
    "position": dict(fmt="s", color="gray", label="bad position"),
    "cloud":    dict(fmt="x", color="steelblue", label="cloud (low ensemble)"),
    "outlier":  dict(fmt="o", color="red", mfc="none", label="outlier"),
}

CONV_FLAGS = {"N-up_E-left": (False, False), "flip-x": (True, False),
              "flip-y": (False, True), "flip-x_flip-y": (True, True)}


def clip_trajectory(mjd, x, y):
    keep = np.ones(len(mjd), bool)
    t = mjd - mjd.min()
    for _ in range(3):
        if keep.sum() < 5:
            break
        for arr in (x, y):
            c = np.polyfit(t[keep], arr[keep], 2)
            keep &= np.abs(arr - np.polyval(c, t)) < TRAJ_CLIP
    return keep


def measure_multi(img, positions):
    pos = np.atleast_2d(positions)
    ann = CircularAnnulus(pos, r_in=R_IN, r_out=R_OUT)
    bkg = ApertureStats(img, ann).median
    out = np.full((len(RADII), len(pos)), np.nan)
    for i, r in enumerate(RADII):
        ap = CircularAperture(pos, r=r)
        f = aperture_photometry(img, ap)["aperture_sum"].data - bkg * ap.area
        f[f <= 0] = np.nan
        out[i] = f
    return out


def detrend_position(rel, x, y, fit_mask):
    good = fit_mask & np.isfinite(rel)
    if good.sum() < 5:
        return rel
    A = np.column_stack([np.ones(good.sum()), x[good] - x.mean(), y[good] - y.mean()])
    coef, *_ = np.linalg.lstsq(A, rel[good], rcond=None)
    model = coef[0] + coef[1] * (x - x.mean()) + coef[2] * (y - y.mean())
    return rel / model * np.nanmedian(rel[good])


def flag_running_median(hours, v, base_mask, sigma=CLIP_SIGMA, win=7):
    out = np.zeros(len(v), bool)
    idx = np.where(base_mask & np.isfinite(v))[0]
    vv = v[idx]
    for k, i in enumerate(idx):
        lo, hi = max(0, k - win // 2), min(len(vv), k + win // 2 + 1)
        loc = np.delete(vv[lo:hi], k - lo)
        m = np.nanmedian(loc)
        s = 1.4826 * np.nanmedian(np.abs(loc - m))
        if s > 0 and abs(vv[k] - m) > sigma * s:
            out[i] = True
    return out


def bin_curve(hours, v, width_min=BIN_MIN):
    if len(hours) == 0:
        return np.array([]), np.array([]), np.array([])
    edges = np.arange(0, hours.max() + width_min / 60, width_min / 60)
    tb, vb, eb = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (hours >= lo) & (hours < hi) & np.isfinite(v)
        if m.sum() >= 2:
            tb.append(hours[m].mean())
            vb.append(v[m].mean())
            eb.append(v[m].std(ddof=1) / np.sqrt(m.sum()))
    return np.array(tb), np.array(vb), np.array(eb)


def transit_window_mjd(target, mjd_lo, mjd_hi):
    """Predicted (mid, ingress, egress) in MJD for the night, or None."""
    if target not in EPHEMERIDES:
        return None
    t0, p, dur_h, depth = EPHEMERIDES[target]
    mid_jd = (mjd_lo + mjd_hi) / 2 + MJD_TO_JD
    tc = t0 + round((mid_jd - t0) / p) * p - MJD_TO_JD
    return tc, tc - dur_h / 2 / 24, tc + dur_h / 2 / 24


def main():
    os.makedirs(PHOTDIR, exist_ok=True)
    coords, fields = load_catalogs()
    conv = {(r["night"], r["target"], r["session"]): r["convention"]
            for r in csv.DictReader(open(os.path.join(OUTDIR, "real_target_summary.csv")))}
    weather_map = {(r["night"], r["file"]): int(r["weather"])
                   for r in csv.DictReader(open(os.path.join(OUTDIR, "real_target_positions.csv")))}

    sessions = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(OUTDIR, "star_positions.csv"))):
        if r["detected"] == "yes":
            sessions[(r["night"], r["target"], r["session"])].append(r)

    all_curves = []
    for key in sorted(sessions):
        night, target, sess = key
        rows = sessions[key]
        if len(rows) < MIN_FRAMES or conv.get(key) not in CONV_FLAGS:
            print(f"[skip] {night} {target} ({len(rows)} frames)")
            continue

        mjd = np.array([float(r["mjd_obs"]) for r in rows])
        tx = np.array([float(r["x"]) for r in rows])
        ty = np.array([float(r["y"]) for r in rows])
        n = len(rows)

        status = np.array(["ok"] * n, dtype=object)
        wx = np.array([weather_map.get((night, r["file"]), 100) for r in rows])
        status[wx <= MIN_WEATHER] = "weather"
        traj_ok = clip_trajectory(mjd, tx, ty)
        status[(~traj_ok) & (status == "ok")] = "position"
        base = status == "ok"
        if base.sum() < MIN_FRAMES:
            print(f"[skip] {night} {target}: only {base.sum()} clean frames")
            continue

        sdir = os.path.join(DB, "observations", night, target, sess)
        h = fits.getheader(os.path.join(sdir, rows[0]["file"]))
        W, H = h["NAXIS1"], h["NAXIS2"]
        w = build_wcs(h, *CONV_FLAGS[conv[key]])
        t_ra, t_dec = coords[target]
        t_cat = np.array(w.wcs_world2pix(t_ra, t_dec, 0))
        cat = fields[target]
        cxy = np.array(w.wcs_world2pix(cat[:, 0], cat[:, 1], 0)).T
        deltas = cxy - t_cat
        sep = np.hypot(*deltas.T)
        t_g = cat[sep.argmin(), 2]
        ok = (sep > 4) & (cat[:, 2] <= t_g + COMP_DMAG)
        for i in np.where(ok)[0]:
            px, py = tx + deltas[i, 0], ty + deltas[i, 1]
            if (px.min() < EDGE or px.max() > W - EDGE or
                    py.min() < EDGE or py.max() > H - EDGE):
                ok[i] = False
        comp_idx = np.where(ok)[0][np.argsort(cat[ok, 2])][:N_COMPS]
        if len(comp_idx) < 2:
            print(f"[skip] {night} {target}: only {len(comp_idx)} comps")
            continue

        dark = master_dark(night)
        nR = len(RADII)
        t_flux = np.full((nR, n), np.nan)
        c_flux = np.full((nR, n, len(comp_idx)), np.nan)
        for j, r in enumerate(rows):
            img = fits.getdata(os.path.join(sdir, r["file"])).astype(float) - dark
            t_flux[:, j] = measure_multi(img, (tx[j], ty[j]))[:, 0]
            det, _, _, _ = detect_stars(img, threshold_sigma=5, max_stars=60)
            if len(det) == 0:
                continue
            for k, i in enumerate(comp_idx):
                px, py = tx[j] + deltas[i, 0], ty[j] + deltas[i, 1]
                d = np.hypot(det[:, 0] - px, det[:, 1] - py)
                m = d.argmin()
                if d[m] < MATCH_R:
                    c_flux[:, j, k] = measure_multi(img, (det[m, 0], det[m, 1]))[:, 0]

        good_c = np.isfinite(c_flux[0][base]).mean(axis=0) >= MIN_MATCH
        c_flux = c_flux[:, :, good_c]
        comp_ids = comp_idx[good_c]
        if c_flux.shape[2] < 2:
            print(f"[skip] {night} {target}: comps too unstable")
            continue

        best_r, best_rms = 0, np.inf
        for ri in range(nR):
            gf = base & np.isfinite(t_flux[ri]) & np.all(np.isfinite(c_flux[ri]), axis=1)
            if gf.sum() < MIN_FRAMES:
                continue
            chk_ = c_flux[ri, gf, 0] / c_flux[ri, gf, 1:].sum(axis=1)
            rms = np.nanstd(chk_ / np.nanmedian(chk_))
            if rms < best_rms:
                best_r, best_rms = ri, rms
        tf, cf = t_flux[best_r], c_flux[best_r]

        measurable = np.isfinite(tf) & np.all(np.isfinite(cf), axis=1)
        status[(~measurable) & (status == "ok")] = "no_flux"
        base = status == "ok"

        ens = np.where(measurable, np.nansum(cf, axis=1), np.nan)
        cloudy = measurable & (ens < CLOUD_CUT * np.nanmedian(ens[base]))
        status[cloudy & (status == "ok")] = "cloud"
        base = status == "ok"
        if base.sum() < MIN_FRAMES:
            print(f"[skip] {night} {target}: too few clear-sky frames")
            continue

        with np.errstate(invalid="ignore"):
            rel = tf / ens
            chk = cf[:, 0] / cf[:, 1:].sum(axis=1)
        rel = detrend_position(rel, tx, ty, base)
        chk = detrend_position(chk, tx, ty, base)
        rel /= np.nanmedian(rel[base])
        chk /= np.nanmedian(chk[base])

        hours = (mjd - mjd.min()) * 24
        out_mask = flag_running_median(hours, rel, base)
        status[out_mask & (status == "ok")] = "outlier"
        base = status == "ok"

        rms_t = np.nanstd(rel[base])
        rms_c = np.nanstd(chk[base])
        sigma_frame = 1.4826 * np.nanmedian(np.abs(rel[base] - np.nanmedian(rel[base])))
        hb, rb, re_ = bin_curve(hours[base], rel[base])

        name = f"{night}_{target}_{sess}"
        with open(os.path.join(PHOTDIR, f"{name}_lightcurve.csv"), "w", newline="") as fh:
            wcsv = csv.writer(fh)
            wcsv.writerow(["file", "mjd_obs", "ut_time", "x", "y", "weather",
                           "target_flux", "ensemble_flux", "rel_flux",
                           "check_flux", "status", "used"]
                          + [f"comp{cid}_flux" for cid in comp_ids])
            for j, r in enumerate(rows):
                ut = Time(mjd[j], format="mjd").isot[11:19]
                wcsv.writerow([r["file"], mjd[j], ut,
                               round(tx[j], 2), round(ty[j], 2), wx[j],
                               "" if np.isnan(tf[j]) else round(tf[j], 1),
                               "" if np.isnan(ens[j]) else round(ens[j], 1),
                               "" if np.isnan(rel[j]) else round(rel[j], 4),
                               "" if np.isnan(chk[j]) else round(chk[j], 4),
                               status[j], "yes" if status[j] == "ok" else "no"]
                              + [round(v, 1) if np.isfinite(v) else ""
                                 for v in cf[j]])
        with open(os.path.join(PHOTDIR, f"{name}_lightcurve_binned.csv"), "w", newline="") as fh:
            wcsv = csv.writer(fh)
            wcsv.writerow(["hours", "rel_flux", "rel_flux_err"])
            for a, b, e in zip(hb, rb, re_):
                wcsv.writerow([round(a, 4), round(b, 4), round(e, 4)])

        # ---- detailed graph: real UT time, every frame shown ----
        times = Time(mjd, format="mjd").datetime
        fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 7.5), sharex=True,
                                     gridspec_kw={"height_ratios": [2.4, 1]})
        tw = transit_window_mjd(target, mjd.min(), mjd.max())
        if tw:
            tc, t_in, t_eg = tw
            for ax in (a1, a2):
                ax.axvspan(Time(t_in, format="mjd").datetime,
                           Time(t_eg, format="mjd").datetime,
                           color="orange", alpha=0.13, zorder=0)
                ax.axvline(Time(tc, format="mjd").datetime,
                           color="orange", lw=0.9, ls="--", zorder=0)

        tt = np.array(times)
        a1.errorbar(tt[base], rel[base], yerr=sigma_frame, fmt="o", ms=4.5,
                    color="black", ecolor="0.6", elinewidth=0.8, capsize=2,
                    label=f"clean frame ({base.sum()})", zorder=3)
        for st, style in STATUS_STYLE.items():
            m = (status == st) & np.isfinite(rel)
            if m.any():
                kw = dict(ms=5, ls="none", zorder=2)
                if "mfc" in style:
                    kw["mfc"] = style["mfc"]
                a1.plot(tt[m], rel[m], style["fmt"], color=style["color"],
                        label=f"{style['label']} ({m.sum()})", **kw)
        m = status == "no_flux"
        if m.any():
            y0 = np.nanmin(rel) if np.isfinite(rel).any() else 0.9
            a1.plot(tt[m], np.full(m.sum(), y0), "|", color="0.4", ms=10,
                    label=f"unmeasurable ({m.sum()})", zorder=2)
        a1.axhline(1, color="gray", lw=0.5)
        a1.set_ylabel("relative flux")
        a1.legend(fontsize=7, ncol=3, loc="best", framealpha=0.9)
        a1.set_title(f"{target}  night {night}  UT times  "
                     f"r={RADII[best_r]}px  {cf.shape[1]} comps  "
                     f"per-frame precision {sigma_frame*100:.1f}%"
                     + ("  (orange = predicted transit)" if tw else ""))

        a2.errorbar(tt[base], chk[base], yerr=rms_c, fmt="o", ms=3.5,
                    color="tab:blue", ecolor="lightsteelblue",
                    elinewidth=0.8, capsize=2)
        a2.axhline(1, color="gray", lw=0.5)
        a2.set_ylabel(f"check star\nrms={rms_c*100:.1f}%")
        a2.set_xlabel("UT time")
        a2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        a2.xaxis.set_major_locator(mdates.MinuteLocator(byminute=[0, 30]))
        plt.tight_layout()
        plt.savefig(os.path.join(PHOTDIR, f"{name}_lightcurve.png"), dpi=110)
        plt.close(fig)

        all_curves.append((target, night, hb, rb, re_, np.nanstd(rb)))
        counts = {s: int((status == s).sum()) for s in
                  ["ok", "weather", "position", "cloud", "outlier", "no_flux"]
                  if (status == s).any()}
        print(f"{night} {target}: r={RADII[best_r]}px, {cf.shape[1]} comps, "
              f"frame rms {rms_t*100:.1f}%, check rms {rms_c*100:.1f}%, "
              f"frames {counts}")

    ncol = 3
    nrow = max(1, int(np.ceil(len(all_curves) / ncol)))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 2.9 * nrow),
                             squeeze=False)
    for ax, (target, night, hb, rb, re_, binrms) in zip(axes.ravel(), all_curves):
        ax.errorbar(hb, rb, yerr=re_, fmt="ko", ms=4, capsize=2)
        if len(hb) > 1:
            gap = np.where(np.diff(hb) > 2.5 * BIN_MIN / 60)[0]
            hb2 = np.insert(hb.astype(float), gap + 1, np.nan)
            rb2 = np.insert(rb.astype(float), gap + 1, np.nan)
            ax.plot(hb2, rb2, "k-", lw=0.8)
        ax.axhline(1, color="gray", lw=0.5)
        ax.set_title(f"{target} {night}  binned rms={binrms*100:.1f}%", fontsize=9)
        ax.set_xlabel("hours"); ax.set_ylabel("rel flux")
    for ax in axes.ravel()[len(all_curves):]:
        ax.axis("off")
    plt.suptitle(f"Binned ({BIN_MIN:.0f} min) differential light curves")
    plt.tight_layout()
    plt.savefig(os.path.join(PHOTDIR, "all_sessions.png"), dpi=110)
    print(f"\nOutputs in {PHOTDIR}")


if __name__ == "__main__":
    main()
