# Exoplanet Transit Photometry Pipeline — HACK4DEV Iraq 2026

A from-scratch pipeline that takes raw MicroObservatory FITS images of
exoplanet host stars and produces calibrated, differential light curves with
predicted transit windows overlaid — turning uncalibrated pixels into a
measurement of each star's brightness over time.

## What it does

1. **Astrometric solving** (`find_real_targets.py`) — the FITS headers carry
   no WCS, so a tangent-plane WCS is built from the pointing and plate scale,
   the Gaia DR3 field is projected through all four camera orientations, and
   pattern-matched to the stars actually detected in each frame. This yields
   the true pixel position of every target and the correct orientation
   (the data is **flip-y**).
2. **Sub-pixel centroiding** (`get_star_positions.py`) — refines each target
   position with a photutils centroid on the dark-subtracted frame.
3. **Differential photometry** (`photometry.py`) — dark-subtracts every
   frame, does aperture photometry on the target and an ensemble of Gaia
   comparison stars of similar brightness, and divides the two to cancel
   clouds, airmass and transparency. Frames are never silently dropped —
   they are classified (weather / cloud / bad-position / outlier) and only
   clean frames feed the statistics.
4. **Per-star campaign reports** (`combine_star.py`) — every observed night
   on a real time axis, every frame plotted, with the **predicted transit
   light curve** (`check_transits.py`, from published ephemerides) drawn on
   each night whether or not it produced data.

## Datasets included

| Folder                | Frames | Description                                   |
|-----------------------|--------|-----------------------------------------------|
| `database/`           | 361    | **Cleaned** data (quality-filtered) — pipeline input |
| `database_original/`  | 1741   | Original raw data, 8 stars over 22 nights (reference / reproducibility) |

The primary analysis runs on `database/`. To rerun on the raw set instead,
point the scripts' `DB` at `database_original/`.

## Run order

```bash
python3 find_real_targets.py    # solve fields vs Gaia -> results/real_target_positions.csv
python3 get_star_positions.py   # sub-pixel target positions -> results/star_positions.csv
python3 photometry.py           # per-session differential light curves
python3 combine_star.py         # per-star campaign reports with transit predictions
```

Requires: `numpy`, `astropy`, `photutils`, `scipy`, `matplotlib`.

## Outputs (in `results/`)

- `photometry/<night>_<target>_session_01_lightcurve.png` — per-night light
  curve on real UT time, every frame shown with its quality flag.
- `photometry/<target>_combined.png` — the campaign report per star: overview
  strip across all nights + one readable panel per night with the predicted
  transit overlaid.
- `photometry/*.csv` — the numbers behind every graph.
- `real_targets_map.png` — the Gaia-matched target position in each field.

## Results summary

Two credible transit signatures were recovered from the cleaned data:
**Qatar-1 (2026-08-21)** and **TRES-5 (2026-08-29)**, both showing a dip of
roughly the predicted depth at the predicted time. Field solving succeeded on
83% of the cleaned frames. Nights dominated by cloud (e.g. Qatar-1 2026-08-18)
are kept but flagged, not hidden.

## Catalogs

- `target_coords.csv` — target RA/DEC (J2000) per star.
- `field_stars.csv` — Gaia DR3 field stars (RA/DEC, G mag) used for solving
  and as photometric comparisons.
- `check_transits.py` — transit ephemerides (T0, period, duration, depth).
