"""Convert a star's RA/DEC (J2000) to pixel X/Y in MicroObservatory FITS images.

The FITS headers contain no WCS keywords, so a tangent-plane (TAN) WCS is
built manually for each frame from:
  - RA, DEC   : telescope pointing, assumed to map to the image center
  - IM_SCALE  : plate scale in arcsec/pixel (5.0"/px for these images)
  - NAXIS1/2  : image dimensions (650 x 500)

Orientation assumption: North up, East left (standard sky view), i.e.
CDELT1 = -scale, CDELT2 = +scale, with Y increasing upward (matplotlib
origin='lower'). The observing session was clouded out (WEATHER=0, no
stars detectable), so this orientation could not be verified against the
sky — use --flip-x / --flip-y to try the other conventions if needed.

Usage:
    python radec_to_xy.py                      # TrES-1, default session
    python radec_to_xy.py --ra 286.0408 --dec 36.6326
    python radec_to_xy.py --data-dir Observations/TRES-1/session_01
"""

import argparse
import csv
import glob
import os

from astropy.io import fits
from astropy.wcs import WCS

# TrES-1 (host star of exoplanet TrES-1b), J2000
DEFAULT_RA = 286.040791   # 19h04m09.79s
DEFAULT_DEC = 36.632622   # +36d37m57.4s

DEFAULT_DATA_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "Observations", "TRES-1", "session_01",
)


def build_wcs(header, flip_x=False, flip_y=False):
    """Build a TAN-projection WCS from a MicroObservatory FITS header."""
    scale = header["IM_SCALE"] / 3600.0  # deg/pixel
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [header["RA"], header["DEC"]]
    # Pointing assumed at image center (FITS pixels are 1-indexed)
    w.wcs.crpix = [(header["NAXIS1"] + 1) / 2.0, (header["NAXIS2"] + 1) / 2.0]
    sx = 1.0 if flip_x else -1.0   # default: East left
    sy = -1.0 if flip_y else 1.0   # default: North up
    w.wcs.cdelt = [sx * scale, sy * scale]
    w.wcs.equinox = header.get("EQUINOX", 2000.0)
    return w


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--ra", type=float, default=DEFAULT_RA, help="star RA in degrees (J2000)")
    p.add_argument("--dec", type=float, default=DEFAULT_DEC, help="star DEC in degrees (J2000)")
    p.add_argument("--data-dir", default=DEFAULT_DATA_DIR, help="directory with FITS frames")
    p.add_argument("--out", default="star_xy.csv", help="output CSV path")
    p.add_argument("--flip-x", action="store_true", help="flip East-West orientation")
    p.add_argument("--flip-y", action="store_true", help="flip North-South orientation")
    args = p.parse_args()

    files = sorted(glob.glob(os.path.join(args.data_dir, "*.fits")))
    if not files:
        raise SystemExit(f"No FITS files found in {args.data_dir}")

    rows = []
    for f in files:
        header = fits.getheader(f)
        w = build_wcs(header, args.flip_x, args.flip_y)
        # origin=0 -> numpy/matplotlib 0-indexed pixel coordinates
        x, y = w.wcs_world2pix(args.ra, args.dec, 0)
        rows.append({
            "file": os.path.basename(f),
            "mjd_obs": header["MJD-OBS"],
            "ut_obs": header["UT-OBS"],
            "x_pix": round(float(x), 2),
            "y_pix": round(float(y), 2),
        })

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    print(f"Star position: RA={args.ra:.6f}  DEC={args.dec:.6f}  (J2000)")
    print(f"{len(rows)} frames -> {args.out}\n")
    print(f"{'file':<28} {'MJD-OBS':>12} {'X (pix)':>9} {'Y (pix)':>9}")
    for r in rows:
        print(f"{r['file']:<28} {r['mjd_obs']:>12.3f} {r['x_pix']:>9.2f} {r['y_pix']:>9.2f}")


if __name__ == "__main__":
    main()
