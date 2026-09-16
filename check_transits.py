"""Check whether an exoplanet transit was predicted during each observed session.

Uses literature ephemerides (NASA Exoplanet Archive / ETD values): mid-transit
epoch T0 (BJD_TDB), orbital period P (days) and total duration T14 (hours).
Each session's frame-time span (from results/real_target_positions.csv) is
compared with the predicted transit windows.

Note: predictions propagated >10 yr from T0 carry ~10-30 min uncertainty,
and BJD_TDB vs UTC differs by up to ~8 min — fine for yes/no overlap checks.
"""

import csv
import os
from collections import defaultdict

from find_real_targets import OUTDIR

# folder -> (T0 [BJD_TDB], P [d], T14 [h], depth [%])
EPHEMERIDES = {
    "TRES-1":  (2453186.8060,  3.0300699,  2.5, 2.3),
    "TRES-3":  (2454185.9101,  1.30618581, 1.4, 2.7),
    "TRES-5":  (2455443.2515,  1.48224686, 1.8, 2.1),
    "Qatar-1": (2455518.4109,  1.42002504, 1.6, 2.1),
    "WASP-2":  (2453991.5146,  2.15222144, 1.8, 1.7),
    "WASP-10": (2454357.8581,  3.0927616,  2.2, 2.9),
    "CoRoT-2": (2454237.5356,  1.7429964,  2.3, 2.7),
    "HATP-10": (2454729.9063,  3.7224747,  2.7, 1.5),
}

MJD_TO_JD = 2400000.5


def main():
    sessions = defaultdict(list)
    for r in csv.DictReader(open(os.path.join(OUTDIR, "real_target_positions.csv"))):
        sessions[(r["night"], r["target"])].append(float(r["mjd_obs"]))

    print(f"{'night':<12} {'target':<9} {'obs window (JD frac)':<24} "
          f"{'transit mid':<13} {'sep(h)':>7}  verdict")
    for (night, target), mjds in sorted(sessions.items()):
        if target not in EPHEMERIDES:
            print(f"{night:<12} {target:<9} no ephemeris")
            continue
        t0, p, dur_h, depth = EPHEMERIDES[target]
        start = min(mjds) + MJD_TO_JD
        end = max(mjds) + MJD_TO_JD
        mid = (start + end) / 2
        # nearest transit epoch to the middle of the observation
        n = round((mid - t0) / p)
        tc = t0 + n * p
        half = dur_h / 2 / 24
        # overlap of [tc-half, tc+half] with [start, end]
        ov_lo, ov_hi = max(start, tc - half), min(end, tc + half)
        overlap_h = max(0.0, (ov_hi - ov_lo) * 24)
        sep_h = (tc - mid) * 24
        if overlap_h > 0:
            frac = overlap_h / dur_h
            verdict = (f"TRANSIT ({overlap_h:.1f}h of {dur_h:.1f}h in window, "
                       f"depth {depth}%)")
        else:
            verdict = "no transit"
        print(f"{night:<12} {target:<9} "
              f"{start:.3f} - {end:.3f}   {tc:.4f} {sep_h:>7.1f}  {verdict}")


if __name__ == "__main__":
    main()
