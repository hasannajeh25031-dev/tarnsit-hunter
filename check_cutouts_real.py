"""Cutouts at the REAL target positions solved by find_real_targets.py.

For every session where the field was solved against Gaia, cut a SIZE x SIZE
stamp around the per-frame real target X/Y (so telescope drift is followed
frame by frame, unlike the old fixed-position check_cutouts.py). Frames are
dark-subtracted. Per session the script:
  - shows a sample of individual-frame cutouts
  - median-stacks all cutouts
  - reports peak significance vs local background

A real star must appear as a bright core at the stamp center.
"""

import csv
import os
from collections import defaultdict

import numpy as np
from astropy.io import fits
from astropy.nddata import Cutout2D
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from find_real_targets import DB, OUTDIR, master_dark

BASE = os.path.dirname(os.path.abspath(__file__))
SIZE = 61
N_SAMPLE = 5


def main():
    rows = list(csv.DictReader(open(os.path.join(OUTDIR, "real_target_positions.csv"))))
    sessions = defaultdict(list)
    for r in rows:
        if r["solved"] == "yes" and r["x_real"]:
            sessions[(r["night"], r["target"], r["session"])].append(r)

    keys = sorted(sessions)
    print(f"{len(keys)} solved sessions")

    fig, axes = plt.subplots(len(keys), N_SAMPLE + 1,
                             figsize=(2.1 * (N_SAMPLE + 1), 2.1 * len(keys)),
                             squeeze=False)

    print(f"\n{'session':<34} {'frames':>6} {'stack peak':>10} {'bkg std':>8} {'S/N':>7}  star?")
    summary = []
    for row_i, key in enumerate(keys):
        night, target, sess = key
        frames = sessions[key]
        dark = master_dark(night)
        cuts = []
        for r in frames:
            path = os.path.join(DB, "observations", night, target, sess, r["file"])
            data = fits.getdata(path).astype(float) - dark
            cut = Cutout2D(data, position=(float(r["x_real"]), float(r["y_real"])),
                           size=SIZE, mode="partial", fill_value=np.nan)
            cuts.append(cut.data)
        cuts = np.array(cuts)

        stack = np.nanmedian(cuts, axis=0)
        bkg = np.nanmedian(stack)
        std = 1.4826 * np.nanmedian(np.abs(stack - bkg))
        core = stack[SIZE // 2 - 2:SIZE // 2 + 3, SIZE // 2 - 2:SIZE // 2 + 3]
        snr = (np.nanmax(core) - bkg) / std if std > 0 else 0.0
        detected = "YES" if snr > 5 else "no"
        label = f"{night} {target} {sess}"
        print(f"{label:<34} {len(frames):>6} {np.nanmax(core) - bkg:>10.1f} "
              f"{std:>8.2f} {snr:>7.1f}  {detected}")
        summary.append({"night": night, "target": target, "session": sess,
                        "n_frames": len(frames),
                        "stack_peak": round(float(np.nanmax(core) - bkg), 1),
                        "bkg_std": round(float(std), 2),
                        "snr": round(float(snr), 1), "star": detected})

        idx = np.linspace(0, len(cuts) - 1, N_SAMPLE).astype(int)
        for col, i in enumerate(idx):
            ax = axes[row_i, col]
            img = cuts[i]
            lo, hi = np.nanpercentile(img, [5, 99.5])
            ax.imshow(img, cmap="gray", vmin=lo, vmax=hi, origin="lower")
            ax.set_xticks([]); ax.set_yticks([])
        ax = axes[row_i, N_SAMPLE]
        lo, hi = np.nanpercentile(stack, [5, 99.5])
        ax.imshow(stack, cmap="gray", vmin=lo, vmax=hi, origin="lower")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"stack S/N={snr:.0f}", fontsize=7)
        axes[row_i, 0].set_ylabel(f"{target}\n{night}", fontsize=6)
        for col in range(N_SAMPLE + 1):
            axes[row_i, col].plot(SIZE // 2, SIZE // 2, "r+", ms=9, mew=0.7, alpha=0.7)

    axes[0, 0].set_title("sample frames", fontsize=8)
    plt.suptitle(f"{SIZE}x{SIZE} cutouts at REAL (Gaia-solved) target positions", y=1.0)
    plt.tight_layout()
    out = os.path.join(BASE, "cutouts_check_real.png")
    plt.savefig(out, dpi=110)
    with open(os.path.join(OUTDIR, "cutouts_real_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=summary[0].keys())
        w.writeheader(); w.writerows(summary)
    print(f"\nFigure : {out}")
    print(f"Summary: {os.path.join(OUTDIR, 'cutouts_real_summary.csv')}")


if __name__ == "__main__":
    main()
