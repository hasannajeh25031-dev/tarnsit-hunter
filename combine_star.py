"""Per-star campaign report: overview across all nights + readable panels.

For every star in the database, one figure:
  - top strip: the whole campaign on a real date axis — one median point per
    night with data, "no data" marks for the others, empty where nothing
    was observed;
  - below: one panel per scheduled night, on the real UT clock time of that
    night, with every measured frame (grey), the 10-min binned curve
    (black) and the PREDICTION: the expected transit light curve (orange,
    from published depth/duration/timing) drawn whether or not the night
    produced data.

Nights share one brightness scale through the comparison stars common to
every night with data; otherwise nights are self-normalized (stated).

Outputs per star in results/photometry/:
  <target>_combined.png / <target>_combined.csv
"""

import csv
import glob
import os
import re
from collections import defaultdict

import numpy as np
from astropy.time import Time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from find_real_targets import OUTDIR
from photometry import PHOTDIR, bin_curve, BIN_MIN
from check_transits import EPHEMERIDES, MJD_TO_JD

MIN_NIGHT_FRAMES = 4
NOISY_RMS = 0.12
NCOL = 3
YLIM = (0.86, 1.14)


def mjd_to_dt(mjd):
    return np.array(Time(np.atleast_1d(mjd), format="mjd").datetime)


def transit_model_mjd(target, mlo, mhi, n=400):
    """(mjd_grid, model_flux, tc) for the night, or None. Trapezoid shape."""
    if target not in EPHEMERIDES:
        return None
    t0, p, dur_h, depth = EPHEMERIDES[target]
    mid_jd = (mlo + mhi) / 2 + MJD_TO_JD
    tc = t0 + round((mid_jd - t0) / p) * p - MJD_TO_JD
    pad = 0.25 / 24
    grid = np.linspace(min(mlo, tc - dur_h / 48) - pad,
                       max(mhi, tc + dur_h / 48) + pad, n)
    half = dur_h / 2 / 24            # T14/2 in days
    flat = 0.6 * half                # assume T23 = 0.6 T14
    x = np.abs(grid - tc)
    model = np.ones_like(grid)
    ing = (x >= flat) & (x < half)
    model[x < flat] = 1 - depth / 100
    model[ing] = 1 - depth / 100 * (half - x[ing]) / (half - flat)
    return grid, model, tc


def main():
    windows = defaultdict(lambda: [np.inf, -np.inf])
    for r in csv.DictReader(open(os.path.join(OUTDIR, "real_target_positions.csv"))):
        k = (r["target"], r["night"])
        m = float(r["mjd_obs"])
        windows[k][0] = min(windows[k][0], m)
        windows[k][1] = max(windows[k][1], m)
    targets = sorted({t for t, _ in windows})
    campaign_lo = min(w[0] for w in windows.values()) - 0.7
    campaign_hi = max(w[1] for w in windows.values()) + 0.7

    nights = defaultdict(dict)
    for path in sorted(glob.glob(os.path.join(PHOTDIR, "*_lightcurve.csv"))):
        m = re.match(r"(\d{4}-\d{2}-\d{2})_(.+)_session_\d+_lightcurve\.csv",
                     os.path.basename(path))
        if not m:
            continue
        night, target = m.groups()
        rows = [r for r in csv.DictReader(open(path)) if r["used"] == "yes"]
        if len(rows) >= MIN_NIGHT_FRAMES:
            nights[target][night] = rows

    print(f"{len(targets)} stars, {len(windows)} scheduled nights\n")

    for target in targets:
        per_night = nights.get(target, {})
        dates = sorted(per_night)
        sched = sorted(n for t_, n in windows if t_ == target)

        common = []
        if dates:
            comp_sets = [set(k for k in per_night[d][0] if k.startswith("comp"))
                         for d in dates]
            common = sorted(set.intersection(*comp_sets))
        if common:
            def ratio_of(r):
                return float(r["target_flux"]) / sum(float(r[c]) for c in common)
            scale_note = f"{len(common)} shared comps"
        elif dates:
            def ratio_of(r):
                return float(r["target_flux"]) / float(r["ensemble_flux"])
            scale_note = "nights normalized individually (no shared comps)"
        else:
            scale_note = "no usable photometry — predictions only"
        if dates:
            g_med = np.median(np.concatenate(
                [[ratio_of(r) for r in per_night[d]] for d in dates]))

        nrow_panels = int(np.ceil(len(sched) / NCOL))
        fig = plt.figure(figsize=(4.3 * NCOL, 2.6 + 2.9 * nrow_panels))
        gs = fig.add_gridspec(1 + nrow_panels, NCOL,
                              height_ratios=[1.1] + [1.6] * nrow_panels,
                              hspace=0.55, wspace=0.25)

        # ---------- top strip: whole campaign ----------
        ax0 = fig.add_subplot(gs[0, :])
        out_rows = []
        night_data = {}
        for night in sched:
            mlo, mhi = windows[(target, night)]
            mid_dt = mjd_to_dt((mlo + mhi) / 2)[0]
            if night in per_night:
                rows = per_night[night]
                mjd = np.array([float(r["mjd_obs"]) for r in rows])
                order = np.argsort(mjd)
                mjd = mjd[order]
                rows = [rows[i] for i in order]
                ratio = np.array([ratio_of(r) for r in rows]) / g_med
                if not common:
                    ratio /= np.median(ratio)
                noisy = np.std(ratio) > NOISY_RMS
                med = float(np.median(ratio))
                err = (1.4826 * np.median(np.abs(ratio - med))
                       / np.sqrt(len(ratio)))
                night_data[night] = (mjd, ratio, rows, noisy, med, err)
                mk = dict(color="0.5", mfc="white") if noisy else dict(color="k")
                ax0.errorbar([mid_dt], [med], yerr=[err], fmt="s", ms=6,
                             capsize=3, **mk)
                ax0.annotate(night[5:], (mid_dt, med),
                             textcoords="offset points", xytext=(0, 9),
                             ha="center", fontsize=7)
                for r, m_, v in zip(rows, mjd, ratio):
                    out_rows.append({"target": target, "night": night,
                                     "file": r["file"], "mjd_obs": m_,
                                     "ut_time": r.get("ut_time", ""),
                                     "rel_flux_global": round(v, 4),
                                     "noisy_night": "yes" if noisy else "no"})
            else:
                ax0.plot([mid_dt], [1.0], "x", ms=6, color="0.6")
                ax0.annotate(f"{night[5:]}\nno data", (mid_dt, 1.0),
                             textcoords="offset points", xytext=(0, 8),
                             ha="center", fontsize=6.5, color="0.5")
        ax0.axhline(1, color="gray", lw=0.5)
        ax0.set_xlim(mjd_to_dt(campaign_lo)[0], mjd_to_dt(campaign_hi)[0])
        ax0.xaxis.set_major_locator(mdates.DayLocator(interval=3))
        ax0.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        ax0.grid(axis="x", color="0.93", lw=0.6)
        ax0.set_ylabel("night median")
        ax0.set_ylim(0.85, 1.15)
        ax0.set_title(f"{target} — campaign overview ({scale_note})",
                      fontsize=11)

        # ---------- one readable panel per night ----------
        for i, night in enumerate(sched):
            ax = fig.add_subplot(gs[1 + i // NCOL, i % NCOL])
            mlo, mhi = windows[(target, night)]

            tm = transit_model_mjd(target, mlo, mhi)
            if tm:
                grid, model, tc = tm
                ax.plot(mjd_to_dt(grid), model, "-", color="darkorange",
                        lw=1.6, alpha=0.85, zorder=2,
                        label="predicted transit")
                ax.axvline(mjd_to_dt(tc)[0], color="darkorange", lw=0.7,
                           ls="--", alpha=0.6, zorder=1)

            if night in night_data:
                mjd, ratio, rows, noisy, med, err = night_data[night]
                ax.plot(mjd_to_dt(mjd), ratio, "o", ms=3,
                        color="lightgray", zorder=3)
                hours = (mjd - mjd.min()) * 24
                hb, rb, re_ = bin_curve(hours, ratio)
                if len(hb):
                    tb = mjd_to_dt(mjd.min() + hb / 24)
                    col = "0.45" if noisy else "black"
                    ax.errorbar(tb, rb, yerr=re_, fmt="o", ms=4.5, lw=1,
                                capsize=2, color=col, zorder=4)
                    gap = np.where(np.diff(hb) > 2.5 * BIN_MIN / 60)[0]
                    tb2 = np.insert(tb, gap + 1, None)
                    rb2 = np.insert(rb.astype(float), gap + 1, np.nan)
                    ax.plot(tb2, rb2, "-", lw=0.8, color=col, zorder=4)
                extra = " NOISY" if noisy else ""
                title = f"{night}   {med:.3f}  ({len(rows)} fr){extra}"
                if not (YLIM[0] < med < YLIM[1]):
                    title += "  [median off scale]"
                ax.set_title(title, fontsize=8.5)
            else:
                ax.text(0.5, 0.35, "no data\n(clouds / unsolved field)",
                        ha="center", va="center", transform=ax.transAxes,
                        fontsize=9, color="0.45")
                ax.set_title(f"{night}", fontsize=8.5)

            ax.axhline(1, color="gray", lw=0.5)
            ax.set_ylim(*YLIM)
            ax.set_xlim(mjd_to_dt(min(mlo, tm[0][0]) if tm else mlo)[0],
                        mjd_to_dt(max(mhi, tm[0][-1]) if tm else mhi)[0])
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
            ax.xaxis.set_major_locator(mdates.HourLocator())
            ax.tick_params(labelsize=7.5)
            if i % NCOL == 0:
                ax.set_ylabel("relative flux")
            if i == 0:
                ax.legend(fontsize=7, loc="lower right")
            ax.set_xlabel("UT time", fontsize=8)

        out_png = os.path.join(PHOTDIR, f"{target}_combined.png")
        fig.savefig(out_png, dpi=110, bbox_inches="tight")
        plt.close(fig)

        if out_rows:
            with open(os.path.join(PHOTDIR, f"{target}_combined.csv"),
                      "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=out_rows[0].keys())
                w.writeheader(); w.writerows(out_rows)

        got = sorted(night_data)
        print(f"{target}: {len(sched)} nights, {len(got)} with data -> "
              f"{os.path.basename(out_png)}")


if __name__ == "__main__":
    main()
