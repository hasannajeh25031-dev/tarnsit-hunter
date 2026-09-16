"""Verify the predicted star position by cutting out square stamps (SIZE px).

For every FITS frame, a SIZE x SIZE Cutout2D is taken around the predicted X/Y of
TrES-1 (all four flip conventions, since orientation is unverified). Frames
are dark-subtracted first. For each candidate position the script:
  - shows a sample of individual-frame cutouts
  - median-stacks all cutouts (beats noise down by ~sqrt(N))
  - reports peak significance vs local background

A real star should appear as a persistent bright core at the stamp center.
"""

import glob
import os

import numpy as np
from astropy.io import fits
from astropy.nddata import Cutout2D
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from radec_to_xy import build_wcs, DEFAULT_RA, DEFAULT_DEC, DEFAULT_DATA_DIR

BASE = os.path.dirname(os.path.abspath(__file__))
DARK_DIR = os.path.join(BASE, "calibration", "2026-08-06")
SIZE = 61  # cutout size in pixels


def load_master_dark():
    darks = sorted(glob.glob(os.path.join(DARK_DIR, "*.fits")))
    return np.median([fits.getdata(d).astype(float) for d in darks], axis=0)


def main():
    files = sorted(glob.glob(os.path.join(DEFAULT_DATA_DIR, "*.fits")))
    dark = load_master_dark()

    conventions = [
        ("N-up E-left (default)", False, False),
        ("flip-x", True, False),
        ("flip-y", False, True),
        ("flip-x flip-y", True, True),
    ]

    # Predicted pixel position for each convention (same for all frames,
    # pointing is constant) taken from the first header.
    h0 = fits.getheader(files[0])
    positions = []
    for name, fx, fy in conventions:
        w = build_wcs(h0, flip_x=fx, flip_y=fy)
        x, y = w.wcs_world2pix(DEFAULT_RA, DEFAULT_DEC, 0)
        positions.append((name, float(x), float(y)))
        print(f"{name:<24} -> x={x:7.2f}  y={y:7.2f}")

    # Collect cutouts of every frame at every candidate position
    stacks = {name: [] for name, _, _ in positions}
    for f in files:
        data = fits.getdata(f).astype(float) - dark
        for name, x, y in positions:
            cut = Cutout2D(data, position=(x, y), size=SIZE, mode="partial")
            stacks[name].append(cut.data)

    n_sample = 6
    fig, axes = plt.subplots(len(positions), n_sample + 1,
                             figsize=(2.1 * (n_sample + 1), 2.3 * len(positions)))
    sample_idx = np.linspace(0, len(files) - 1, n_sample).astype(int)

    print(f"\n{'position':<24} {'stack peak':>10} {'bkg std':>8} {'S/N':>6}  star?")
    for row, (name, x, y) in enumerate(positions):
        cuts = np.array(stacks[name])
        for col, i in enumerate(sample_idx):
            ax = axes[row, col]
            img = cuts[i]
            lo, hi = np.percentile(img, [5, 99.5])
            ax.imshow(img, cmap="gray", vmin=lo, vmax=hi, origin="lower")
            ax.set_xticks([]); ax.set_yticks([])
            if row == 0:
                ax.set_title(f"frame {i}", fontsize=9)

        stack = np.median(cuts, axis=0)
        bkg = np.median(stack)
        std = 1.4826 * np.median(np.abs(stack - bkg))  # robust sigma
        core = stack[SIZE // 2 - 2:SIZE // 2 + 3, SIZE // 2 - 2:SIZE // 2 + 3]
        snr = (core.max() - bkg) / std if std > 0 else 0.0
        detected = "YES" if snr > 5 else "no"
        print(f"{name:<24} {core.max() - bkg:>10.2f} {std:>8.2f} {snr:>6.1f}  {detected}")

        ax = axes[row, n_sample]
        lo, hi = np.percentile(stack, [5, 99.5])
        ax.imshow(stack, cmap="gray", vmin=lo, vmax=hi, origin="lower")
        ax.set_xticks([]); ax.set_yticks([])
        if row == 0:
            ax.set_title("median stack (91)", fontsize=9)
        axes[row, 0].set_ylabel(f"{name}\n({x:.0f}, {y:.0f})", fontsize=8)

    for row in range(len(positions)):
        for col in range(n_sample + 1):
            # mark stamp center where the star should sit
            axes[row, col].plot(SIZE // 2, SIZE // 2, "r+", ms=10, mew=0.8, alpha=0.7)

    plt.suptitle(f"{SIZE}x{SIZE} cutouts at predicted TrES-1 positions (dark-subtracted)")
    plt.tight_layout()
    out = os.path.join(BASE, "cutouts_check.png")
    plt.savefig(out, dpi=110)
    print(f"\nFigure saved to {out}")


if __name__ == "__main__":
    main()
