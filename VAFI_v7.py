#!/usr/bin/env python3
"""Venus night-side flash search in BMP/FITS sequences or AVI videos (v7).

VAFI (Venus Atmospheric Flash Investigator) - Version 7
========================================================

Python software for the automatic detection of possible fireballs in
Venus's nightside atmosphere, inspired by the ALFI (Automatic Lunar
Flash Investigation) algorithm developed for lunar impact detection: 
3x3 local maxima above N spatial standard deviations AND above M temporal 
standard deviations, using 1x1/2x2/3x3... spatial block averaging and 1/2/3... 
frame temporal averaging. 
Different from ALFI, frames are registered on Venus's bright disk before 
temporal comparisons. This is an independent implementation of the 
publicly described ALFI method, not a reproduction of the original ALFI 
source code.

The program accepts sequences of BMP or FITS images, as well as AVI
videos, and performs the following operations:

1. VENUS PHASE AND DISK GEOMETRY
   The positions of Earth and Venus are computed from fixed Keplerian
   orbital elements at the observation date (UTC). The Sun-Venus-Earth
   phase angle, illuminated fraction, and apparent diameter are
   calculated. Using the image scale (arcsec/pixel), an analytical
   model of the illuminated Venus disk and its terminator is generated.

2. MANUAL MASK FITTING
   In manual mode, the theoretical Venus profile is superimposed on
   the first image in an interactive window. The user can translate,
   rotate, and resize the profile to obtain the best fit to the
   observed disk. The accepted geometry is saved in a JSON file and
   can be reused in subsequent analyses. Automatic intensity-based,
   geometric, and phase-based masking modes are also available.

3. IMAGE ALIGNMENT
   Since telescope tracking errors cause Venus to shift between
   consecutive frames, all images are automatically registered to
   the first frame using template matching on the bright hemisphere.
   Subpixel translations are applied, and poorly aligned frames
   are excluded from the analysis.

4. NIGHT-SIDE MASK
   The fitted theoretical profile is used to exclude the illuminated
   hemisphere, a region around the terminator, and the planetary
   limb. Additional masking reduces contamination from scattered
   light. Only the remaining nightside pixels are searched for flashes.

5. SPATIAL AND TEMPORAL FLASH DETECTION
   The algorithm searches for local brightness maxima within 3x3
   pixel regions exceeding the eight neighboring pixels by at least
   N standard deviations. Each candidate must also exceed the
   temporal background at the same registered position by at least
   M standard deviations.

   The analysis is repeated using spatial averaging over different
   pixel blocks (1x1, 2x2, 3x3, ...) and temporal averaging over
   different numbers of consecutive frames, allowing the detection
   of flashes with different angular sizes and durations.

6. CANDIDATE SELECTION AND OUTPUT
   Detections close in position and time are merged into individual
   candidate events. The program produces CSV and text tables with
   frame numbers, timestamps, coordinates, detection significance,
   and signal excess, together with diagnostic images, alignment
   information, and the adopted Venus mask.

The software is an independent ALFI-inspired implementation.
Detected flashes are candidates, not confirmed atmospheric impacts.
Visual inspection and, when possible, independent confirmation
are necessary to reject instrumental artifacts and seeing effects.

Main dependencies: NumPy, OpenCV, and Astropy (for FITS input).
Compatible with Python 3.10 or later.


DETECTION PARAMETERS
1.Parameters that directly control the sensitivity to fireballs.

Parameter					Meaning and effect

--spatial-scales 1,2,3,4	Searches the original images and images averaged into 2×2, 3×3 and 4×4 pixel blocks. 
Larger scales help detect spatially extended, faint flashes but can dilute very compact flashes.

--time-scales 1,2,3,4	Searches individual frames and averages of 2, 3 and 4 consecutive frames. 
Improves sensitivity to flashes persisting over several frames, but longer averaging dilutes shorter flashes.

--N 1.5	Minimum spatial significance relative to the eight neighboring pixels. 
Increasing N rejects more spatial candidates, including genuine PSF-like flashes.

--M 7	Minimum temporal significance relative to the same aligned location in surrounding frames. 
Increasing M reduces random temporal false detections but can miss faint flashes.

--min-psf-pixels 2	Requires at least two connected pixels with significant positive temporal residuals. 
Increasing this rejects isolated noise but can miss very compact or weak fireballs.

--pixel-support-sigma 3	Minimum temporal significance of each pixel counted toward the connected-pixel requirement. 
Increasing it makes the multipixel test stricter.

2.Input, timing and mask parameters

Parameter			Meaning and effect
--input				Folder containing the BMP/FITS sequence, or an AVI file.
--output			Destination of detection tables, masks, alignment diagnostics and candidate images.
--start-utc			UTC time of the first frame. In manual mode, also determines the theoretical Venus phase unless --observation-utc is supplied.
--cadence 0.5		Time between consecutive images, in seconds. Used for timestamps, not directly for detection significance.
--plate-scale 0.088	Angular image scale in arcsec/pixel. Determines the predicted Venus radius for fitting.
--mask-mode manual	Uses the Keplerian phase model with manually accepted disk position, radius and orientation.
--manual-config		Loads the saved geometric fit, avoiding repeated manual fitting. The file must correspond to the same input's first frame and observation parameters.

3.The multipixel filter introduced in version 7
This is the most interesting improvement, but the parameters must be carefully selected.
The program first detects candidates using N and M. Only afterwards does it apply the additional multipixel requirement.
It calculates the positive temporal residual at each original-resolution pixel:

\[\Delta I(x,y)=I_{\mathrm{event}}(x,y)-I_{\mathrm{background}}(x,y)\]

A supporting pixel must satisfy
\[\Delta I(x,y)\geq K\sigma_{\mathrm{pixel}}(x,y)\]

where \(K\) is --pixel-support-sigma.

The script then finds connected groups of supporting pixels using eight-neighbor connectivity.
Effect of the two parameters:

Configuration			Expected behavior
--min-psf-pixels 1		Disables the multipixel rejection test
--min-psf-pixels 2		Rejects isolated significant pixels
--min-psf-pixels 3		Requires a more extended signal
--min-psf-pixels 4		More conservative for well-sampled PSFs

--pixel-support-sigma 2	More sensitive to faint PSF wings
--pixel-support-sigma 3	Requires stronger supporting pixels
--pixel-support-sigma 4	Very restrictive for faint events

EXAMPLE:
    Create the manual mask for Venus:
    py VAFI_v7.py --input "D:\Venus\BMP" --output "D:\Venus\MaskPreview_BMP" --start-utc 2026-07-05T18:00:00 --plate-scale 0.1 --mask-mode manual --preview-only

    Use the manual mask:
    py VAFI_v7.py --input "D:\Venus\BMP" --output "D:\Venus\results" --start-utc 2026-07-05T18:00:00 --cadence 0.0302 --plate-scale 0.1 --mask-mode manual --manual-config "D:\Venus\MaskPreview_BMP\manual_mask_geometry.json" --spatial-scales 1,2,3 --time-scales 1,2,3 --N 1.5 --M 2 --min-psf-pixels 15

    py VAFI_v7.py --input "D:\Venus\BMP" --output "D:\Venus\results" --start-utc 2026-07-05T18:00:00 --cadence 0.0302
    py VAFI_v7.py --input "D:\Venus\FITS" --output "D:\Venus\results"
    py VAFI_v7.py --input "D:\Venus\observation.avi" --output "D:\Venus\results" --start-utc 2026-07-05T18:00:00 --spatial-scales 1,2,3,4 --time-scales 1,2,3,4
    
For BMP images, timestamps can be assigned with e.g.:
    --start-utc 2026-07-05T18:00:00 --cadence 0.030
For AVI, frame cadence is automatically read from the movie's FPS metadata;
--cadence overrides it when necessary. To obtain UTC timestamps from AVI,
pass --start-utc (AVI files generally do not store an absolute UTC start).

Albino Carbognani,
INAF-OAS, version November 02, 2026

"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import csv
import hashlib
import json
import os
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
import re
import shutil
import tempfile

import cv2
import numpy as np

SUPPORTED = {".bmp", ".fit", ".fits", ".fts"}
VIDEO_EXTENSIONS = {".avi"}


def natural_key(path: Path):
    return [int(t) if t.isdigit() else t.casefold()
            for t in re.split(r"(\d+)", path.name)]


def get_files(folder: Path) -> list[Path]:
    if not folder.is_dir():
        raise ValueError(f"Input folder does not exist: {folder}")
    files = sorted((p for p in folder.iterdir()
                    if p.is_file() and p.suffix.lower() in SUPPORTED), key=natural_key)
    if not files:
        raise ValueError("No BMP or FITS images found in the input folder")
    types = {"bmp" if f.suffix.lower() == ".bmp" else "fits" for f in files}
    if len(types) != 1:
        raise ValueError("Please use a folder containing BMP OR FITS frames, not both")
    return files


def read_image(path: Path) -> tuple[np.ndarray, str]:
    """Return floating-point 2-D image and FITS DATE-OBS string, if any."""
    date_obs = ""
    if path.suffix.lower() == ".bmp":
        img = cv2.imread(str(path), cv2.IMREAD_UNCHANGED)
        if img is None:
            raise OSError(f"Cannot read BMP: {path}")
        if img.ndim == 3:
            img = cv2.cvtColor(img, cv2.COLOR_BGRA2GRAY if img.shape[2] == 4
                               else cv2.COLOR_BGR2GRAY)
        if img.ndim != 2:
            raise ValueError(f"Expected a 2-D BMP image: {path}")
        img = img.astype(np.float32)
    else:
        try:
            from astropy.io import fits
        except ImportError as exc:
            raise ImportError("FITS files require 'pip install astropy'") from exc
        with fits.open(path, memmap=False) as hdu_list:
            valid_hdu = next((h for h in hdu_list if h.data is not None
                              and h.data.ndim == 2), None)
            if valid_hdu is None:
                raise ValueError(f"No two-dimensional image HDU in {path}")
            img = np.asarray(valid_hdu.data, dtype=np.float32)
            date_obs = str(valid_hdu.header.get(
                "DATE-OBS", hdu_list[0].header.get("DATE-OBS", "")))
            exposure = valid_hdu.header.get(
                "EXPTIME", hdu_list[0].header.get("EXPTIME", 0))
            # DATE-OBS normally marks exposure START; record center if possible.
            if date_obs:
                try:
                    dt = datetime.fromisoformat(date_obs.strip().replace("Z", "+00:00"))
                    dt += timedelta(seconds=float(exposure) / 2)
                    date_obs = as_utc(dt).isoformat(timespec="milliseconds")
                except (ValueError, TypeError, OverflowError):
                    pass
    if not np.isfinite(img).all():
        good = np.isfinite(img)
        if not good.any():
            raise ValueError(f"All pixels are invalid: {path}")
        img[~good] = np.median(img[good])
    return img, date_obs


def avi_to_image(frame: np.ndarray, path: Path) -> np.ndarray:
    """Convert a decoded AVI frame to the detector's 2-D float image.

    OpenCV normally decodes video as 8-bit BGR. We deliberately do not
    histogram-normalize each frame: temporal comparisons need a stable
    intensity scale. A compressed/8-bit AVI cannot recover original SER
    bit depth or lost faint flashes.
    """
    if frame.ndim == 3:
        if frame.shape[2] == 4:
            frame = cv2.cvtColor(frame, cv2.COLOR_BGRA2GRAY)
        elif frame.shape[2] == 3:
            frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        elif frame.shape[2] == 1:
            frame = frame[:, :, 0]
        else:
            raise ValueError(f"Unsupported AVI frame color channels in {path}")
    if frame.ndim != 2:
        raise ValueError(f"Unexpected AVI frame dimensions in {path}: {frame.shape}")
    frame = frame.astype(np.float32)
    if not np.isfinite(frame).all():
        good = np.isfinite(frame)
        if not good.any():
            raise ValueError(f"All pixels in AVI frame are invalid: {path}")
        frame[~good] = np.median(frame[good])
    return frame


def open_avi_capture(path: Path) -> cv2.VideoCapture:
    cap = cv2.VideoCapture(str(path))
    if not cap.isOpened():
        cap.release()
        raise OSError(f"Cannot open AVI video: {path}. Check the installed video codec.")
    return cap


def count_avi_frames(path: Path) -> tuple[int, float | None]:
    """Count actual AVI frames (some AVI frame-count headers are inaccurate).

    One fast demux/grab pass ensures the disk-backed data cube is allocated
    for the correct number of frames; frames are decoded on the second pass.
    """
    print(f"Counting frames in AVI: {path.name}", flush=True)
    cap = open_avi_capture(path)
    try:
        fps = float(cap.get(cv2.CAP_PROP_FPS))
        if not np.isfinite(fps) or fps <= 0:
            fps = None
        count = 0
        while cap.grab():
            count += 1
            if count % 10000 == 0:
                print(f"  Counted {count} AVI frames", flush=True)
    finally:
        cap.release()
    if count == 0:
        raise ValueError(f"No readable frames found in AVI: {path}")
    return count, fps


@dataclass
class FrameSource:
    """Common sequential interface for original image folders and AVI files."""
    input_path: Path
    files: list[Path] | None
    nframes: int
    first: np.ndarray
    fps: float | None = None

    @property
    def is_video(self) -> bool:
        return self.files is None

    def frame_name(self, i: int) -> str:
        if self.files is not None:
            return self.files[i].name
        return f"{self.input_path.name}:frame_{i:08d}"

    def iter_frames(self):
        """Yield (image, FITS exposure-center UTC if available)."""
        if self.files is not None:
            for file_path in self.files:
                yield read_image(file_path)
            return
        cap = open_avi_capture(self.input_path)
        try:
            for i in range(self.nframes):
                ok, frame = cap.read()
                if not ok or frame is None:
                    raise RuntimeError(
                        f"AVI decoder stopped at frame {i} of {self.nframes}. "
                        "Check whether the AVI file is truncated or uses an unsupported codec.")
                yield avi_to_image(frame, self.input_path), ""
            ok, _ = cap.read()
            if ok:
                raise RuntimeError(
                    "More AVI frames were decoded than counted. Please convert "
                    "the AVI to an uncompressed image sequence and rerun.")
        finally:
            cap.release()


def open_frame_source(path: Path) -> FrameSource:
    if path.is_dir():
        files = get_files(path)
        first, _ = read_image(files[0])
        return FrameSource(path, files, len(files), first)
    if path.is_file() and path.suffix.lower() in VIDEO_EXTENSIONS:
        nframes, fps = count_avi_frames(path)
        cap = open_avi_capture(path)
        try:
            ok, frame = cap.read()
            if not ok or frame is None:
                raise RuntimeError(f"Cannot decode the first AVI frame: {path}")
            first = avi_to_image(frame, path)
        finally:
            cap.release()
        return FrameSource(path, None, nframes, first, fps)
    raise ValueError(
        "--input must be a folder of BMP/FITS frames or the filename of an AVI video")


def as_utc(dt: datetime) -> datetime:
    return dt.replace(tzinfo=timezone.utc) if dt.tzinfo is None else dt.astimezone(timezone.utc)


def parse_start(s: str) -> datetime:
    return as_utc(datetime.fromisoformat(s.strip().replace("Z", "+00:00")))


def robust_sigma(a: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64).ravel()
    if a.size == 0:
        return 0.0
    med = np.median(a)
    return float(1.4826 * np.median(np.abs(a - med)))


# -----------------------------------------------------------------------------
# Phase-constrained Venus geometry. No network, no SPICE, no third-party ephemeris.
# Fixed Keplerian J2000 elements from NASA/JPL SSD Approximate Positions of the
# Planets, Table 1 (1800-2050): https://ssd.jpl.nasa.gov/planets/approx_pos.html
# We deliberately keep all orbital elements fixed and propagate the mean anomaly
# with n=sqrt(mu/a^3), as requested (strictly two-body Keplerian approximation).
# For precise absolute angles, substitute a numerical ephemeris, not this model.
# -----------------------------------------------------------------------------
J2000_JD = 2451545.0
GAUSSIAN_K = 0.01720209895                  # radians / day for 1 AU
KM_PER_AU = 149597870.700
VENUS_RADIUS_KM = 6051.8
OBLIQUITY_J2000 = np.deg2rad(23.439291111)
# (a in AU, eccentricity, inclination, mean longitude, longitude of perihelion,
#  longitude of ascending node; all angles degrees at J2000.)
KEPLER_ELEMENTS = {
    'earth': (1.00000261, 0.01671123, -0.00001531,
              100.46457166, 102.93768193, 0.0),  # Earth-Moon barycenter
    'venus': (0.72333566, 0.00677672, 3.39467605,
              181.97909950, 131.60246718, 76.67984255),
}


def kepler_position(body: str, jd: float) -> np.ndarray:
    """Heliocentric ecliptic J2000 xyz [AU] using a fixed Kepler ellipse."""
    a, e, inc, L, varpi, node = KEPLER_ELEMENTS[body]
    n = GAUSSIAN_K / a ** 1.5
    M = np.deg2rad((L - varpi) % 360.0) + n * (jd - J2000_JD)
    M = float((M + np.pi) % (2*np.pi) - np.pi)
    E = M
    for _ in range(20):
        dE = (E - e * np.sin(E) - M) / (1.0 - e*np.cos(E))
        E -= dE
        if abs(dE) < 1e-14:
            break
    xp = a*(np.cos(E) - e)
    yp = a*np.sqrt(1.0-e*e)*np.sin(E)
    w = np.deg2rad(varpi - node)
    Om = np.deg2rad(node)
    incl = np.deg2rad(inc)
    # rotate from perifocal to J2000 ecliptic
    cw, sw, cO, sO, ci, si = np.cos(w), np.sin(w), np.cos(Om), np.sin(Om), np.cos(incl), np.sin(incl)
    return np.array([
        (cw*cO-sw*sO*ci)*xp + (-sw*cO-cw*sO*ci)*yp,
        (cw*sO+sw*cO*ci)*xp + (-sw*sO+cw*cO*ci)*yp,
        (sw*si)*xp + (cw*si)*yp,
    ])


def unit(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v)


def venus_phase_geometry(utc: datetime, plate_scale: float) -> dict:
    """Keplerian phase, angular radius and Sunward sky position angle.

    The Earth is represented by the Earth-Moon barycenter. The projected Sun
    position angle is east of celestial north; actual camera roll/mirroring is
    unknown and is FIT from each image instead of assumed.
    """
    if plate_scale <= 0 or not np.isfinite(plate_scale):
        raise ValueError('The plate scale must be a positive number (arcsec/pixel)')
    jd = 2440587.5 + as_utc(utc).timestamp() / 86400.0
    re, rv = kepler_position('earth', jd), kepler_position('venus', jd)
    earth_from_v = unit(re - rv)
    sun_from_v = unit(-rv)
    alpha = float(np.arccos(np.clip(np.dot(earth_from_v, sun_from_v), -1, 1)))
    delta = float(np.linalg.norm(rv - re))
    radius_arcsec = float(np.rad2deg(np.arcsin(VENUS_RADIUS_KM / (delta*KM_PER_AU))) * 3600.0)
    # J2000 equatorial sky basis for Venus as seen from Earth
    line_ecl = rv - re
    eps = OBLIQUITY_J2000
    eq = np.array([line_ecl[0], line_ecl[1]*np.cos(eps)-line_ecl[2]*np.sin(eps),
                   line_ecl[1]*np.sin(eps)+line_ecl[2]*np.cos(eps)])
    ra, dec = np.arctan2(eq[1],eq[0]), np.arctan2(eq[2],np.hypot(eq[0],eq[1]))
    sun_eq = np.array([sun_from_v[0], sun_from_v[1]*np.cos(eps)-sun_from_v[2]*np.sin(eps),
                       sun_from_v[1]*np.sin(eps)+sun_from_v[2]*np.cos(eps)])
    east = np.array([-np.sin(ra), np.cos(ra), 0.0])
    north = np.array([-np.sin(dec)*np.cos(ra), -np.sin(dec)*np.sin(ra), np.cos(dec)])
    sun_pa = float(np.rad2deg(np.arctan2(np.dot(sun_eq,east), np.dot(sun_eq,north))) % 360.0)
    return dict(date_utc=as_utc(utc).isoformat(), jd=jd,
                earth_venus_au=delta, phase_angle_deg=float(np.rad2deg(alpha)),
                illuminated_fraction=float((1+np.cos(alpha))/2),
                lambert_integrated_brightness_fraction=lambert_phase_function(float(np.rad2deg(alpha))),
                radius_arcsec=radius_arcsec, radius_px=radius_arcsec/plate_scale,
                diameter_arcsec=2*radius_arcsec,
                sun_pa_east_of_north_deg=sun_pa,
                plate_scale_arcsec_px=plate_scale,
                ephemeris='fixed J2000 two-body Keplerian, Earth-Moon barycenter')


def phase_maps(shape: tuple[int, int], cx: float, cy: float,
               radius: float, phase_deg: float, sun_angle_deg: float):
    """Exact projected illuminated silhouette and local incidence cosine.

    Coordinates (u,v) are in Venus radii, with u toward the projected Sun and
    v orthogonal to it. The visible circle has u^2+v^2 <= 1. On that circle,
    the projected terminator is the ellipse
        u_terminator(v) = -cos(alpha)*sqrt(1-v^2).
    This is the SHAPE of the illuminated part (NOT the integrated Lambert phase
    function). Alpha is the Sun-Venus-Earth angle, in degrees.
    """
    if radius <= 0:
        raise ValueError('Projected Venus radius must be positive')
    h,w=shape
    yy,xx=np.indices(shape,dtype=np.float32)
    x=(xx-cx)/radius
    y=(yy-cy)/radius
    theta=np.deg2rad(sun_angle_deg)
    ct,st=np.cos(theta),np.sin(theta)
    u=ct*x+st*y
    v=-st*x+ct*y
    rho2=u*u+v*v
    disk=rho2<=1.0
    a=np.deg2rad(phase_deg)
    # Analytic ellipse: positive u is sunward in the camera image.
    u_terminator=-np.cos(a)*np.sqrt(np.maximum(1.0-v*v,0.0))
    day=disk & (u>=u_terminator)
    # Lambertian local surface brightness, up to an arbitrary scale; useful
    # only as a FITTING TEMPLATE, not a calibrated model for Venusian clouds.
    z=np.sqrt(np.maximum(1.0-rho2,0.0))
    mu=np.sin(a)*u+np.cos(a)*z
    return disk,day,mu


def terminator_polyline(cx: float, cy: float, radius: float,
                        phase_deg: float, image_sun_angle_deg: float,
                        n: int = 501) -> np.ndarray:
    """2-D analytical terminator ellipse, pixel coords; open arc limb-to-limb."""
    v=np.linspace(-1.0,1.0,n,dtype=np.float64)
    u=-np.cos(np.deg2rad(phase_deg))*np.sqrt(np.maximum(1.0-v*v,0.0))
    theta=np.deg2rad(image_sun_angle_deg)
    ct,st=np.cos(theta),np.sin(theta)
    x=cx+radius*(ct*u-st*v)
    y=cy+radius*(st*u+ct*v)
    return np.rint(np.stack([x,y],axis=-1)).astype(np.int32).reshape(-1,1,2)


def lambert_phase_function(alpha_deg: float) -> float:
    """Disk-INTEGRATED Lambert brightness, normalized to full phase.

    Unlike phase_maps(), this returns ONE NUMBER, not a disk silhouette.
    https://pyastronomy.readthedocs.io/en/latest/pyaslDoc/aslDoc/planetPhase.html
    """
    a=np.deg2rad(alpha_deg)
    return float((np.sin(a)+(np.pi-a)*np.cos(a))/np.pi)


def phase_template(radius: float, phase_deg: float,
                   image_sun_angle_deg: float,
                   profile: str = 'binary'):
    """Analytical bright silhouette optionally weighted by Lambert incidence.

    The binary model is appropriate for a nearly saturated Venus disk, while
    Lambert is an optional approximate radiance template for unsaturated data.
    """
    half=int(np.ceil(radius))+3
    template_shape=(2*half+1,2*half+1)
    _,day,mu=phase_maps(template_shape,half,half,radius,
                        phase_deg,image_sun_angle_deg)
    if profile=='binary':
        templ=day.astype(np.float32)
    elif profile=='lambert':
        # Observed contrast is square-root transformed during registration.
        # Approximate the same transform for the Lambertian brightness model.
        templ=np.sqrt(np.clip(mu,0,1)).astype(np.float32)
        templ[~day]=0.0
    else:
        raise ValueError('Unknown phase profile '+profile)
    templ=cv2.GaussianBlur(templ,(0,0),max(0.8,0.015*radius))
    return templ,half


def phase_fit(img: np.ndarray, radius: float, phase_deg: float,
              fixed_center: tuple[float,float]|None = None,
              fixed_angle_deg: float|None = None,
              fit_profile: str = 'auto') -> dict:
    """Fit predicted bright SHAPE to camera pixels; physical phase/radius fixed.

    Camera angle and center are unknown and searched. The profile can be a
    binary analytical silhouette or a Lambert-incidence-weighted version.
    The latter is a geometrical approximation, NOT a Venus-cloud radiative
    transfer model. Both use normalized cross-correlation to reduce reliance
    on arbitrary video gain/exposure.
    """
    h,w=img.shape
    fit_scale=min(1.0,90.0/max(radius,1.0),600.0/max(h,w))
    sz=(max(24,round(w*fit_scale)),max(24,round(h*fit_scale)))
    small=cv2.resize(img,sz,interpolation=cv2.INTER_AREA).astype(np.float32)
    bg=float(np.percentile(small,25))
    high=float(np.percentile(small,99.8))
    if high<=bg:
        raise RuntimeError('Phase fit: Venus not distinguished from sky')
    feature=np.sqrt(np.clip((small-bg)/max(high-bg,1e-6),0,1)).astype(np.float32)
    feature=cv2.GaussianBlur(feature,(0,0),0.9)
    rad=radius*fit_scale
    template_d=2*(int(np.ceil(rad))+3)+1
    if template_d>min(feature.shape):
        raise RuntimeError('Predicted Venus disk does not fit the image. '
                           'Check --plate-scale or use an image with more sky.')
    center_small=(fixed_center[0]*fit_scale,fixed_center[1]*fit_scale) if fixed_center else None
    styles=['binary','lambert'] if fit_profile=='auto' else [fit_profile]

    def evaluate(ang: float,style: str):
        templ,half=phase_template(rad,phase_deg,ang,style)
        if center_small is None:
            corr=cv2.matchTemplate(feature,templ,cv2.TM_CCOEFF_NORMED)
            _,score,_,(ix,iy)=cv2.minMaxLoc(corr)
            fx,fy=subpixel_parabola(corr,ix,iy)
            return float(score),(float(ix+fx+half),float(iy+fy+half))
        sx,sy=center_small
        patch=cv2.getRectSubPix(feature,(templ.shape[1],templ.shape[0]),
                                (float(sx),float(sy)))
        s0=patch-patch.mean()
        t0=templ-templ.mean()
        score=float(np.sum(s0*t0)/max(np.sqrt(np.sum(s0*s0)*np.sum(t0*t0)),1e-10))
        return score,(sx,sy)

    angles=[fixed_angle_deg%360.0] if fixed_angle_deg is not None else np.arange(0.,360.,15.)
    scores=[]
    for style in styles:
        for ang in angles:
            score,xy=evaluate(float(ang),style)
            scores.append((score,float(ang)%360.0,xy,style))
    best=max(scores,key=lambda q:q[0])
    if fixed_angle_deg is None:
        for step,width in [(3.,15.),(0.6,3.)]:
            for da in np.arange(-width,width+step*.1,step):
                angle=(best[1]+da)%360.0
                for style in styles:
                    score,xy=evaluate(angle,style)
                    if score>best[0]:
                        best=(score,angle,xy,style)
    score,angle,(sx,sy),style=best
    return dict(center_x=sx/fit_scale,center_y=sy/fit_scale,
                radius_px=radius,sun_angle_image_deg=angle,
                score=score,fit_downsampling=fit_scale,
                fitted_brightness_profile=style,
                center_fixed=bool(fixed_center is not None),
                sun_angle_fixed=bool(fixed_angle_deg is not None))


def phase_masks(composite: np.ndarray, cx: float, cy: float, radius: float,
                phase_geometry: dict, image_angle: float,
                bright_margin: int, limb_margin: float, glare_fraction: float):
    """Physical terminator mask plus exclusion for high-intensity light spill.

    The bright hemisphere and a dilated terminator band are excluded based on
    predicted GEOMETRY, not a low intensity threshold. This leaves a nightside
    even if the scattered-light background is not black.
    """
    disk,day,_ = phase_maps(composite.shape,cx,cy,radius,
                            phase_geometry['phase_angle_deg'],image_angle)
    yy,xx=np.indices(composite.shape,dtype=np.float32)
    rr=np.hypot(xx-cx,yy-cy)
    safe_limb=rr <= max(0.0,radius-limb_margin)
    k=2*bright_margin+1
    day_excluded=cv2.dilate(np.uint8(day),cv2.getStructuringElement(
        cv2.MORPH_ELLIPSE,(k,k)))>0 if bright_margin else day.copy()
    off_disk=composite[rr >= radius + max(3,bright_margin)]
    bg=float(np.median(off_disk)) if off_disk.size else float(np.percentile(composite,10))
    sig=np.maximum(composite-bg,0)
    peak=float(np.percentile(sig[disk],99)) if disk.any() else 0.0
    # Only veto very bright pixels on the theoretical nightside. Using an
    # intensity threshold near 3% would reject the entire nightside under glare.
    if glare_fraction>=1.0 or peak<=0:
        glare=np.zeros_like(disk)
    else:
        sm=cv2.GaussianBlur(sig,(0,0),1.2)
        glare=(sm>glare_fraction*peak)
        glare=cv2.dilate(np.uint8(glare),cv2.getStructuringElement(
            cv2.MORPH_ELLIPSE,(5,5)))>0
    bright=day_excluded|glare
    night=safe_limb & (~bright)
    night[:2,:]=False
    night[-2:,:]=False
    night[:,:2]=False
    night[:,-2:]=False
    info=dict(mode='phase',center_x=float(cx),center_y=float(cy),radius_px=float(radius),
              disk_pixels=int(disk.sum()),background=bg,
              bright_margin_px=bright_margin,limb_margin_px=limb_margin,
              phase_glare_fraction=glare_fraction,
              theoretical_illuminated_fraction=phase_geometry['illuminated_fraction'],
              image_sun_angle_deg=float(image_angle),
              bright_pixels_in_disk=int((bright&disk).sum()),
              geometric_cap_pixels_in_disk=int((~day&disk).sum()),
              safe_night_pixels=int(night.sum()),
              safe_night_percent_of_disk=round(100*night.sum()/max(disk.sum(),1),3))
    return disk,night,bright,bg,info


def save_phase_overlay(path: Path, img: np.ndarray, disk: np.ndarray,
                       night: np.ndarray, cx: float, cy: float, radius: float,
                       phase_deg: float, image_angle: float):
    """Four-panel, unambiguous QA: the bright outline is NOT a circle.

    The true disk LIMB is circular, but the illuminated portion consists of
    that limb PLUS a curved projected TERMINATOR, shown in cyan/magenta.
    """
    disk,day,_=phase_maps(img.shape,cx,cy,radius,phase_deg,image_angle)
    base=cv2.cvtColor(grayscale(img),cv2.COLOR_GRAY2BGR)
    # Pure ANALYTICAL bright silhouette; nothing circular is drawn around it.
    silhouette=np.zeros_like(base)
    silhouette[day]=(245,245,245)
    # Fitted analytic profile on top of actual Venus video image.
    fit=base.copy()
    fit[day]=(0.48*fit[day]+0.52*np.array([30,115,245])).astype(np.uint8)
    # In the phase model we emphasize the elliptical TERMINATOR, not limb.
    curve=terminator_polyline(cx,cy,radius,phase_deg,image_angle)
    cv2.polylines(fit,[curve],False,(80,255,255),2,cv2.LINE_AA)
    # Outer limb is deliberately much fainter to avoid confusion.
    cv2.circle(fit,(round(cx),round(cy)),round(radius),(125,125,125),1,cv2.LINE_AA)
    theta=np.deg2rad(image_angle)
    end=(round(cx+0.48*radius*np.cos(theta)),
         round(cy+0.48*radius*np.sin(theta)))
    cv2.arrowedLine(fit,(round(cx),round(cy)),end,(255,180,255),1,tipLength=0.22)
    roi=base.copy()
    roi[day]=(0.35*roi[day]+0.65*np.array([30,60,200])).astype(np.uint8)
    roi[night]=(0.25*roi[night]+0.75*np.array([50,210,35])).astype(np.uint8)
    cv2.polylines(roi,[curve],False,(60,255,255),2,cv2.LINE_AA)
    # Standalone mask ONLY, shape must NOT be mistaken for an entire circle.
    mask_only=np.zeros_like(base)
    mask_only[day]=(30,65,200)
    mask_only[night]=(35,210,55)
    cv2.polylines(mask_only,[curve],False,(60,255,255),1,cv2.LINE_AA)
    titles=['ORIGINAL IMAGE','ANALYTICAL BRIGHT SILHOUETTE',
            'FITTED TERMINATOR / LIMB','SAFE NIGHT (GREEN)']
    panels=[base,silhouette,fit,roi]
    visual=[]
    for title,panel in zip(titles,panels):
        strip=np.zeros((34,panel.shape[1],3),dtype=np.uint8)
        cv2.putText(strip,title,(6,22),cv2.FONT_HERSHEY_SIMPLEX,
                    0.46,(255,255,255),1,cv2.LINE_AA)
        visual.append(np.vstack([strip,panel]))
    cv2.imwrite(str(path),np.hstack(visual))
    # Save the mathematical silhouette and the ACTUAL candidate ROI as two
    # independent PNGs for quick inspection, without explanatory circles.
    cv2.imwrite(str(path.with_name('phase_bright_silhouette.png')),silhouette)
    cv2.imwrite(str(path.with_name('phase_search_mask.png')),mask_only)


def write_phase_geometry(path: Path, ephem: dict, fit: dict):
    with path.open('w',encoding='utf-8') as f:
        f.write('VENUS TWO-BODY KEPLERIAN PHASE AND IMAGE FIT\n')
        f.write('NOTE: image roll/parity are fitted, not known from the observation time.\n')
        f.write('No light-time, planetary perturbations, topocentric corrections.\n')
        f.write('Projected TERMINATOR: u = -cos(alpha)*sqrt(1-v^2).\n')
        f.write('Bright SHAPE is NOT the circular limb or integrated Lambert phase function.\n')
        for label, d in [('EPHEMERIS',ephem),('IMAGE FIT',fit)]:
            f.write('\n'+label+'\n')
            for key,value in d.items():
                f.write(f'{key}: {value}\n')
        f.write('\nJPL source: https://ssd.jpl.nasa.gov/planets/approx_pos.html\n')


# -----------------------------------------------------------------------------
# Manual fitting of the exact phase silhouette (requires a local graphical
# desktop for first use). Coordinates are always ORIGINAL camera pixels.
# Profile has only center x/y, scale and orientation DOFs: its illuminated
# fraction / terminator curvature stays fixed by the Keplerian phase angle.
# -----------------------------------------------------------------------------
def first_frame_digest(img: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(img, dtype=np.float32).tobytes()).hexdigest()


def manual_initial_angle(img: np.ndarray, cx: float, cy: float) -> float:
    """Very rough hint from bright-image centroid, never a fixed constraint."""
    sm = cv2.GaussianBlur(img.astype(np.float32), (0, 0), 3)
    lo, hi = float(np.percentile(sm, 20)), float(np.percentile(sm, 99))
    region = sm > lo + 0.65 * (hi - lo)
    if np.count_nonzero(region) < 5:
        return 0.0
    ys, xs = np.nonzero(region)
    vx, vy = float(np.median(xs) - cx), float(np.median(ys) - cy)
    return float(np.degrees(np.arctan2(vy, vx)) % 360) if np.hypot(vx, vy) >= 3 else 0.0


def manual_display_gray(img: np.ndarray, stretch: int, clip_pct: int) -> np.ndarray:
    """Only alters on-screen visual contrast, never photometric input data."""
    low = float(np.percentile(img, 10))
    top = float(np.percentile(img, 99.9))
    high = low + max(top-low, 1.0) * max(1, clip_pct) / 100.0
    x = np.clip((img-low)/max(1e-6, high-low), 0, 1)
    if stretch == 1:
        x = np.sqrt(x)
    elif stretch == 2:
        x = np.arcsinh(12*x)/np.arcsinh(12)
    return np.uint8(np.rint(x*255))


def render_manual_fit(img: np.ndarray, ephem: dict, geometry: dict,
                      bright_margin: int, limb_margin: float,
                      glare_fraction: float, stretch: int = 1,
                      clip_pct: int = 100, display_max: int = 540):
    """Two side-by-side images; return image, GUI image-panel transform and stats."""
    h,w = img.shape
    cx,cy,r,angle = (float(geometry[k]) for k in
                     ('center_x','center_y','radius_px','sun_angle_image_deg'))
    disk,day,_ = phase_maps(img.shape,cx,cy,r,ephem['phase_angle_deg'],angle)
    _,night,_,_,info = phase_masks(img,cx,cy,r,ephem,angle,
                                  bright_margin,limb_margin,glare_fraction)
    gray = manual_display_gray(img,stretch,clip_pct)
    original=cv2.cvtColor(gray,cv2.COLOR_GRAY2BGR)
    fit=original.copy()
    fit[day]=(fit[day].astype(np.float32)*0.69 +
              np.array([30,75,230],dtype=np.float32)*0.31).astype(np.uint8)
    fit[night]=(fit[night].astype(np.float32)*0.48 +
                np.array([40,220,65],dtype=np.float32)*0.52).astype(np.uint8)
    curve=terminator_polyline(cx,cy,r,ephem['phase_angle_deg'],angle)
    cv2.polylines(fit,[curve],False,(0,255,255),2,cv2.LINE_AA)
    cv2.circle(fit,(round(cx),round(cy)),round(r),(180,180,180),1,cv2.LINE_AA)
    cv2.drawMarker(fit,(round(cx),round(cy)),(255,255,255),cv2.MARKER_CROSS,12,1)
    theta=np.deg2rad(angle)
    end=(round(cx+0.50*r*np.cos(theta)),round(cy+0.50*r*np.sin(theta)))
    cv2.arrowedLine(fit,(round(cx),round(cy)),end,(255,0,255),1,tipLength=0.20)
    scale=min(1.0,display_max/max(w,h))
    vw,vh=max(1,int(round(w*scale))),max(1,int(round(h*scale)))
    left=cv2.resize(original,(vw,vh),interpolation=cv2.INTER_AREA)
    right=cv2.resize(fit,(vw,vh),interpolation=cv2.INTER_AREA)
    pad=10; top=55; footer=95
    canvas=np.zeros((top+vh+footer,2*vw+3*pad,3),dtype=np.uint8)
    canvas[:]=(22,22,22)
    lx=pad; rx=2*pad+vw
    canvas[top:top+vh,lx:lx+vw]=left
    canvas[top:top+vh,rx:rx+vw]=right
    cv2.putText(canvas,'ORIGINAL', (lx,25), cv2.FONT_HERSHEY_SIMPLEX,
                0.55,(225,225,225),1,cv2.LINE_AA)
    cv2.putText(canvas,'MANUAL FIT  red=day  green=SEARCH  yellow=terminator',
                (rx,25),cv2.FONT_HERSHEY_SIMPLEX,0.46,(225,225,225),1,cv2.LINE_AA)
    ybase=top+vh+22
    cv2.putText(canvas, 'Drag RIGHT frame to translate   |   A/D rotate   |   +/- radius   |   I/J/K/L shift 1 px',
                (pad,ybase),cv2.FONT_HERSHEY_SIMPLEX,0.46,(240,240,240),1,cv2.LINE_AA)
    cv2.putText(canvas,'S or SPACE or ENTER: accept / save     ESC: cancel     R: reset to predicted radius',
                (pad,ybase+23),cv2.FONT_HERSHEY_SIMPLEX,0.46,(240,240,240),1,cv2.LINE_AA)
    frac=100*info['safe_night_pixels']/max(1,int(disk.sum()))
    msg=(f"center=({cx:.0f},{cy:.0f})  radius={r:.0f}px  sun-angle={angle:.0f}deg  "
         f"safe night={info['safe_night_pixels']}px ({frac:.1f}% disk)  "
         f"phase={100*ephem['illuminated_fraction']:.1f}%")
    msgcol=(70,80,255) if info['safe_night_pixels']<30 else (100,220,220)
    cv2.putText(canvas,msg,(pad,ybase+49),cv2.FONT_HERSHEY_SIMPLEX,0.45,msgcol,1,cv2.LINE_AA)
    return canvas, (rx,top,vw,vh,scale), info


def manual_mask_gui(img: np.ndarray, ephem: dict, starting: dict,
                    bright_margin: int, limb_margin: float,
                    glare_fraction: float):
    """Undistorted, scrollable Tkinter GUI for the manual Venus phase fit.

    Only this function replaces the OpenCV HighGUI implementation.
    The returned values match the original function: (geometry, full_screenshot).
    """
    import base64
    import tkinter as tk
    from tkinter import messagebox

    h, w = img.shape
    maxrad = max(15, round(1.5 * max(h, w)))
    start = dict(starting)
    start['center_x'] = float(np.clip(start['center_x'], 0, w - 1))
    start['center_y'] = float(np.clip(start['center_y'], 0, h - 1))
    start['radius_px'] = float(np.clip(start['radius_px'], 5, maxrad))
    start['sun_angle_image_deg'] = float(start['sun_angle_image_deg'] % 360)

    root = tk.Tk()
    root.title('Venus - Manual phase and nightside mask')
    sw, sh = root.winfo_screenwidth(), root.winfo_screenheight()
    root.geometry(f'{min(sw - 80, 1000)}x{min(sh - 90, 760)}+20+20')
    root.minsize(460, 320)

    control_frame = tk.Frame(root)
    control_frame.pack(fill='x', side='top', padx=5, pady=3)

    vals = {
        'x': tk.IntVar(master=root, value=round(start['center_x'])),
        'y': tk.IntVar(master=root, value=round(start['center_y'])),
        'r': tk.IntVar(master=root, value=round(start['radius_px'])),
        'a': tk.IntVar(master=root, value=round(start['sun_angle_image_deg']) % 360),
        'clip': tk.IntVar(master=root, value=100),
        'stretch': tk.IntVar(master=root, value=1),
    }
    ui = {'photo': None, 'image_item': None, 'screen': None,
          'viewport': None, 'info': None, 'render_queued': False}
    result = {'accepted': False, 'geometry': None, 'screen': None}
    drag = {'active': False, 'previous': (0.0, 0.0)}

    # A Canvas scrolls the original pixel array. It never stretches it.
    main = tk.Frame(root)
    main.pack(side='top', fill='both', expand=True)
    main.rowconfigure(0, weight=1)
    main.columnconfigure(0, weight=1)
    canvas = tk.Canvas(main, highlightthickness=0, bg='#202020')
    canvas.grid(row=0, column=0, sticky='nsew')
    scroll_x = tk.Scrollbar(main, orient='horizontal', command=canvas.xview)
    scroll_y = tk.Scrollbar(main, orient='vertical', command=canvas.yview)
    scroll_x.grid(row=1, column=0, sticky='ew')
    scroll_y.grid(row=0, column=1, sticky='ns')
    canvas.configure(xscrollcommand=scroll_x.set, yscrollcommand=scroll_y.set)

    status = tk.StringVar(master=root, value='Move the phase profile to Venus, then accept the fit.')
    tk.Label(root, textvariable=status, anchor='w', padx=6).pack(fill='x')
    tk.Label(root,
             text='Mouse drag on RIGHT panel | A/D rotate | +/- radius | I/J/K/L shift | '
                  'S/Space/Enter accept | Esc cancel',
             anchor='w', padx=6).pack(fill='x')

    def geometry():
        return dict(center_x=float(vals['x'].get()),
                    center_y=float(vals['y'].get()),
                    radius_px=float(max(5, vals['r'].get())),
                    sun_angle_image_deg=float(vals['a'].get() % 360))

    def redraw():
        ui['render_queued'] = False
        g = geometry()
        complete, viewport, info = render_manual_fit(
            img, ephem, g, bright_margin, limb_margin, glare_fraction,
            stretch=vals['stretch'].get(), clip_pct=max(1, vals['clip'].get()))
        success, encoded = cv2.imencode('.png', complete)
        if not success:
            raise RuntimeError('Could not create the manual-fit preview')
        encoded64 = base64.b64encode(encoded.tobytes()).decode('ascii')
        photo = tk.PhotoImage(master=root, data=encoded64, format='png')
        # The array and widget have EXACTLY the same number of pixels.
        if (photo.width(), photo.height()) != (complete.shape[1], complete.shape[0]):
            raise RuntimeError('The GUI image dimensions differ from the preview pixels')
        if ui['image_item'] is None:
            ui['image_item'] = canvas.create_image(0, 0, anchor='nw', image=photo)
        else:
            canvas.itemconfigure(ui['image_item'], image=photo)
        ui['photo'] = photo  # Keep the Tk image alive.
        ui['screen'], ui['viewport'], ui['info'] = complete, viewport, info
        canvas.configure(scrollregion=(0, 0, complete.shape[1], complete.shape[0]))
        fraction = 100.0 * info['safe_night_pixels'] / max(1, info['disk_pixels'])
        status.set(f"Center ({g['center_x']:.0f}, {g['center_y']:.0f}), "
                   f"radius {g['radius_px']:.0f}px, angle {g['sun_angle_image_deg']:.0f} deg; "
                   f"safe nightside: {info['safe_night_pixels']}px ({fraction:.1f}% of disk)")

    def request_redraw(*_):
        if not ui['render_queued']:
            ui['render_queued'] = True
            root.after_idle(redraw)

    controls = [
        ('Center X (px)', 'x', w - 1),
        ('Center Y (px)', 'y', h - 1),
        ('Radius (px)', 'r', maxrad),
        ('Sun angle (deg)', 'a', 359),
        ('Display clip (%)', 'clip', 100),
        ('Display stretch 0/1/2', 'stretch', 2),
    ]
    for i, (label, key, maxval) in enumerate(controls):
        column = i % 2
        row = i // 2
        s = tk.Scale(control_frame, from_=(5 if key == 'r' else 0), to=maxval,
                     orient='horizontal', label=label, showvalue=True,
                     variable=vals[key], command=request_redraw,
                     length=300, resolution=1)
        s.grid(row=row, column=column, sticky='ew', padx=5)
    for column in range(2):
        control_frame.columnconfigure(column, weight=1)

    def canvas_coordinates(event):
        return float(canvas.canvasx(event.x)), float(canvas.canvasy(event.y))

    def mouse_press(event):
        x, y = canvas_coordinates(event)
        vx, vy, vw, vh, sc = ui['viewport']
        if vx <= x < vx + vw and vy <= y < vy + vh:
            drag['active'] = True
            drag['previous'] = (x, y)

    def mouse_move(event):
        if not drag['active']:
            return
        x, y = canvas_coordinates(event)
        old_x, old_y = drag['previous']
        sc = ui['viewport'][4]
        dx = int(round((x - old_x) / sc))
        dy = int(round((y - old_y) / sc))
        if dx or dy:
            vals['x'].set(int(np.clip(vals['x'].get() + dx, 0, w - 1)))
            vals['y'].set(int(np.clip(vals['y'].get() + dy, 0, h - 1)))
            drag['previous'] = (x, y)
            request_redraw()

    def mouse_release(event):
        drag['active'] = False

    canvas.bind('<ButtonPress-1>', mouse_press)
    canvas.bind('<B1-Motion>', mouse_move)
    canvas.bind('<ButtonRelease-1>', mouse_release)

    def accept():
        if ui['info'] is None:
            return
        if ui['info']['safe_night_pixels'] < 30:
            messagebox.showwarning('Venus mask',
                'Fewer than 30 safe nightside pixels. Adjust the mask before accepting.',
                parent=root)
            return
        result['accepted'] = True
        result['geometry'] = geometry()
        result['screen'] = ui['screen'].copy()
        root.quit()

    def cancel():
        root.quit()

    def on_key(event):
        k = event.keysym
        c = event.char
        if k == 'Escape':
            cancel()
        elif k in ('Return', 'KP_Enter', 'space') or c in ('s', 'S'):
            accept()
        elif c in ('a','A','d','D'):
            step = 5 if c.isupper() else 1
            change = -step if c.lower() == 'a' else step
            vals['a'].set((vals['a'].get() + change) % 360)
        elif c in ('+', '=') or k == 'KP_Add':
            vals['r'].set(min(maxrad, vals['r'].get() + 1))
        elif c in ('-', '_') or k == 'KP_Subtract':
            vals['r'].set(max(5, vals['r'].get() - 1))
        elif c.lower() in ('i', 'j', 'k', 'l'):
            key, change = {'i': ('y', -1), 'k': ('y', 1),
                           'j': ('x', -1), 'l': ('x', 1)}[c.lower()]
            limit = h - 1 if key == 'y' else w - 1
            vals[key].set(int(np.clip(vals[key].get() + change, 0, limit)))
        elif c.lower() == 'r':
            vals['r'].set(int(np.clip(round(ephem['radius_px']), 5, maxrad)))
        else:
            return
        request_redraw()

    root.bind('<KeyPress>', on_key)
    root.protocol('WM_DELETE_WINDOW', cancel)
    print('MANUAL FIT: drag right image, use sliders or keyboard; '
          'Tkinter scrollbars preserve aspect ratio at 1:1 preview pixels.')
    request_redraw()
    root.mainloop()
    root.destroy()
    if not result['accepted']:
        raise RuntimeError('Manual mask selection was cancelled.')
    return result['geometry'], result['screen']


def write_manual_config(path: Path, source: FrameSource, first: np.ndarray,
                        ephem: dict, geometry: dict, phase_time: datetime,
                        bright_margin: int, limb_margin: float,
                        glare_fraction: float) -> None:
    payload={
        'format': 'venus-alfi-manual-phase-v1',
        'source_path':str(source.input_path.resolve()),
        'first_frame_sha256_float32':first_frame_digest(first),
        'image_shape_hw':list(first.shape),
        'observation_utc':as_utc(phase_time).isoformat(),
        'plate_scale_arcsec_per_pixel':ephem['plate_scale_arcsec_px'],
        'keplerian_phase_angle_deg':ephem['phase_angle_deg'],
        'keplerian_illuminated_fraction':ephem['illuminated_fraction'],
        'predicted_radius_px':ephem['radius_px'],
        'manual_center_x_px':float(geometry['center_x']),
        'manual_center_y_px':float(geometry['center_y']),
        'manual_radius_px':float(geometry['radius_px']),
        'manual_sun_angle_image_deg':float(geometry['sun_angle_image_deg'])%360,
        'bright_margin_px':bright_margin,
        'limb_margin_px':limb_margin,
        'phase_glare_fraction':glare_fraction,
        'note':'Positive x right, positive y down; angle from +x toward +y in image. '
               'Reuse only with the same original first frame and observing parameters.'
    }
    path.write_text(json.dumps(payload,indent=2,ensure_ascii=True)+'\n',encoding='utf-8')


def read_manual_config(path: Path, source: FrameSource, first: np.ndarray,
                       ephem: dict, phase_time: datetime) -> tuple[dict,dict]:
    try:
        payload=json.loads(path.read_text(encoding='utf-8'))
    except (OSError,ValueError) as exc:
        raise ValueError(f'Unable to read manual mask configuration {path}: {exc}') from exc
    if payload.get('format')!='venus-alfi-manual-phase-v1':
        raise ValueError('Unsupported or invalid manual mask configuration format')
    if payload.get('image_shape_hw')!=list(first.shape):
        raise ValueError('Manual geometry was saved for different frame dimensions')
    saved_source = os.path.normcase(os.path.abspath(str(payload.get('source_path',''))))
    current_source = os.path.normcase(os.path.abspath(str(source.input_path)))
    if saved_source != current_source:
        raise ValueError('The saved mask refers to another input path. Refit on this '
                         'recording to avoid applying an obsolete geometry.')
    if payload.get('first_frame_sha256_float32')!=first_frame_digest(first):
        raise ValueError('Manual geometry was saved for a different first frame. '
                         'Run the interactive fit again on this recording.')
    if abs(float(payload.get('plate_scale_arcsec_per_pixel',0))-
           ephem['plate_scale_arcsec_px'])>1e-5:
        raise ValueError('Manual configuration plate scale differs from --plate-scale')
    if abs(float(payload.get('keplerian_phase_angle_deg',0))-
           ephem['phase_angle_deg'])>0.01:
        raise ValueError('Manual configuration phase angle differs from the observation time')
    if payload.get('observation_utc')!=as_utc(phase_time).isoformat():
        raise ValueError('Manual configuration UTC differs from this analysis. '
                         'Use the same --start-utc / --observation-utc.')
    geo=dict(center_x=float(payload['manual_center_x_px']),
             center_y=float(payload['manual_center_y_px']),
             radius_px=float(payload['manual_radius_px']),
             sun_angle_image_deg=float(payload['manual_sun_angle_image_deg']))
    return geo,payload


def estimate_disk(img: np.ndarray, threshold_fraction: float = 0.12) -> tuple[float, float, float]:
    """Estimate Venus's limb by enclosing its strongest extended bright region."""
    bg = float(np.median(img))
    high = float(np.percentile(img, 99.9))
    if high <= bg:
        raise RuntimeError("Could not identify a bright Venus disk")
    sm = cv2.GaussianBlur(img, (0, 0), 1.4)
    threshold = bg + threshold_fraction * (high - bg)
    bw = np.uint8(sm > threshold)
    bw = cv2.morphologyEx(bw, cv2.MORPH_CLOSE,
                          cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5)))
    contours, _ = cv2.findContours(bw, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    contours = [c for c in contours if len(c) >= 20 and cv2.contourArea(c) >= 20]
    if not contours:
        raise RuntimeError("No extended bright crescent detected; specify --cx --cy --radius")
    contour = max(contours, key=cv2.contourArea)
    (cx, cy), radius = cv2.minEnclosingCircle(contour)
    h, w = img.shape
    if not (8 <= radius < 0.5 * min(w, h)):
        raise RuntimeError("Implausible automatic Venus radius; use --cx --cy --radius")
    return float(cx), float(cy), float(radius)


def register_feature(img: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    """High-SNR crescent image for translation-only phase correlation."""
    out_w, out_h = shape
    scaled = cv2.resize(img, (out_w, out_h), interpolation=cv2.INTER_AREA)
    bg = float(np.percentile(scaled, 40))
    hi = float(np.percentile(scaled, 99.8))
    if hi <= bg:
        return np.zeros_like(scaled, dtype=np.float32)
    x = np.clip((scaled - bg) / (hi - bg), 0.0, 1.0)
    x = np.sqrt(x).astype(np.float32)
    return cv2.GaussianBlur(x, (0, 0), 0.75)


def prepare_template(ref_feature: np.ndarray, cx: float, cy: float,
                     radius: float, sx: float, sy: float):
    """Reference crescent template (resized coordinates), excluding much sky."""
    h, w = ref_feature.shape
    rcx, rcy = cx * sx, cy * sy
    half = max(12, int(round(radius * min(sx, sy) * 0.9)))
    x0 = max(0, int(round(rcx)) - half)
    x1 = min(w, int(round(rcx)) + half + 1)
    y0 = max(0, int(round(rcy)) - half)
    y1 = min(h, int(round(rcy)) + half + 1)
    template = np.ascontiguousarray(ref_feature[y0:y1, x0:x1])
    if min(template.shape) < 10 or np.std(template) < 1e-5:
        raise RuntimeError("Cannot extract a reliable registration template")
    return template, x0, y0


def subpixel_parabola(corr: np.ndarray, x: int, y: int) -> tuple[float, float]:
    """Subpixel parabola refinement around a correlation maximum."""
    h, w = corr.shape
    vals = []
    for dx, dy in [(1, 0), (0, 1)]:
        if not (0 <= x-dx < w and 0 <= x+dx < w and
                0 <= y-dy < h and 0 <= y+dy < h):
            vals.append(0.0)
            continue
        left = float(corr[y-dy, x-dx])
        mid = float(corr[y, x])
        right = float(corr[y+dy, x+dx])
        den = left - 2*mid + right
        offset = 0.5 * (left - right) / den if abs(den) > 1e-12 else 0.0
        vals.append(float(np.clip(offset, -0.7, 0.7)))
    return vals[0], vals[1]


def template_shift(feature: np.ndarray, template: np.ndarray,
                   ref_x0: int, ref_y0: int, sx: float, sy: float):
    """Locate the bright crescent, robust to shifts and strong contrast."""
    corr = cv2.matchTemplate(feature, template, cv2.TM_CCOEFF_NORMED)
    _, response, _, (x, y) = cv2.minMaxLoc(corr)
    fx, fy = subpixel_parabola(corr, x, y)
    return (x + fx - ref_x0) / sx, (y + fy - ref_y0) / sy, float(response)

def crop_registered(img: np.ndarray, dx: float, dy: float,
                    crop: tuple[int, int, int, int]) -> np.ndarray:
    x0, y0, x1, y1 = crop
    matrix = np.array([[1.0, 0.0, -dx - x0],
                       [0.0, 1.0, -dy - y0]], dtype=np.float32)
    return cv2.warpAffine(img, matrix, (x1 - x0, y1 - y0),
                          flags=cv2.INTER_LINEAR,
                          borderMode=cv2.BORDER_CONSTANT, borderValue=0)


def build_masks(composite: np.ndarray, cx: float, cy: float, radius: float,
                bright_margin: int, limb_margin: float, bright_fraction: float,
                mask_mode: str = "intensity", night_cap: float = 0.25):
    """Produce a CONSERVATIVE nightside ROI and diagnostics.

    intensity: threshold image by peak brightness (v2-compatible).
    geometric: use direction opposite brightest part of disk, excluding
    its terminator and retaining only the deep nightside cap. This is an
    optional fallback for intense scattered-light halos, not a precise
    photometric fit of the Venusian terminator. Visual QA is essential.
    """
    h, w = composite.shape
    yy, xx = np.indices((h, w), dtype=np.float32)
    rr = np.hypot(xx - cx, yy - cy)
    in_disk = rr <= radius
    safe_limb = rr < (radius - limb_margin)
    bg_pixels = composite[rr >= radius + 2]
    background = float(np.median(bg_pixels)) if bg_pixels.size else float(np.percentile(composite, 15))
    signal = np.clip(composite - background, 0.0, None)
    disk_signal = signal[in_disk]
    if disk_signal.size == 0:
        raise RuntimeError("Empty planetary disk. Check the manual disk parameters")
    peak = float(np.percentile(disk_signal, 99))
    if peak <= 0:
        raise RuntimeError("No bright-side signal found. Check disk geometry")
    smooth_signal = cv2.GaussianBlur(signal, (0, 0), 1.4)
    info = {
        "mode": mask_mode, "center_x": float(cx), "center_y": float(cy),
        "radius_px": float(radius), "disk_pixels": int(in_disk.sum()),
        "background": background, "peak_99_minus_bg": peak,
        "bright_margin_px": bright_margin, "limb_margin_px": limb_margin,
        "bright_fraction": bright_fraction, "night_cap": night_cap,
        "bright_centroid_x": None, "bright_centroid_y": None,
    }

    # The intensity cutoff may be too low when bright Venus saturates the
    # sensor and scattered light extends onto the nightside. Unlike v2,
    # preserve an empty result long enough to write disk_and_mask.png.
    if mask_mode == "intensity":
        exclude_threshold = bright_fraction * peak
        bright = smooth_signal > exclude_threshold
        geometric_cap = np.ones_like(in_disk, dtype=bool)
    elif mask_mode == "geometric":
        # Locate the sunward side using a weighted centroid of the brightest
        # pixels, restricted to the disk interior to reduce stray limb light.
        centroid_region = rr <= 0.90 * radius
        sample = smooth_signal[centroid_region]
        if sample.size < 40:
            raise RuntimeError("Not enough disk pixels for geometric masking")
        high_cut = float(np.percentile(sample, 85))
        tops = centroid_region & (smooth_signal >= high_cut) & (smooth_signal > 0)
        if np.count_nonzero(tops) < 10:
            raise RuntimeError("Cannot determine bright-side direction")
        weights = np.maximum(smooth_signal[tops] - high_cut * 0.7, 1e-6)
        bcx = float(np.average(xx[tops], weights=weights))
        bcy = float(np.average(yy[tops], weights=weights))
        vecx, vecy = bcx-cx, bcy-cy
        length = float(np.hypot(vecx, vecy))
        info["bright_centroid_x"], info["bright_centroid_y"] = bcx, bcy
        info["bright_centroid_offset_px"] = length
        if length < radius * 0.035:
            raise RuntimeError(
                "Bright-side direction is ambiguous. Provide better disk geometry "
                "or use --mask-mode intensity")
        ux, uy = vecx/length, vecy/length
        # Exclude the center and the terminator by retaining only a far cap
        # opposite the bright centroid. A sufficiently distant cap works for
        # crescent and gibbous views without modelling the exact terminator.
        projection = (xx-cx)*ux + (yy-cy)*uy
        geometric_cap = projection <= -night_cap * radius
        # Use a less aggressive photometric limit than intensity mode, because
        # scattered light may exceed 3% of the peak across the entire disk.
        # We STILL mask the brightest (>=20%-peak) pixels as a safety guard.
        exclude_threshold = max(bright_fraction, 0.20) * peak
        bright = smooth_signal > exclude_threshold
    else:
        raise ValueError(f"Unrecognized mask_mode: {mask_mode}")

    bright = np.uint8(bright)
    if bright_margin:
        k = 2 * bright_margin + 1
        bright = cv2.dilate(bright, cv2.getStructuringElement(
            cv2.MORPH_ELLIPSE, (k, k)))
    bright = bright.astype(bool)
    night = safe_limb & geometric_cap & (~bright)
    # Reject border pixels where 3x3 statistics would be incomplete.
    night[:2, :] = False
    night[-2:, :] = False
    night[:, :2] = False
    night[:, -2:] = False
    info["exclusion_signal_threshold"] = float(exclude_threshold)
    info["bright_pixels_in_disk"] = int((bright & in_disk).sum())
    info["geometric_cap_pixels_in_disk"] = int((geometric_cap & in_disk).sum())
    info["safe_night_pixels"] = int(night.sum())
    info["safe_night_percent_of_disk"] = round(100.0 * int(night.sum()) / max(int(in_disk.sum()), 1), 3)
    return in_disk, night, bright, background, info


def save_mask_diagnostics(path: Path, info: dict, preview: bool) -> None:
    with path.open("w", encoding="utf-8") as fp:
        fp.write("VENUS MASK DIAGNOSTICS\n")
        fp.write(f"Mode: {'FIRST FRAME (preview only)' if preview else 'MEDIAN OF ALIGNED FRAMES'}\n")
        for key, value in info.items():
            fp.write(f"{key}: {value}\n")
        fp.write("\nIMPORTANT: Check the green pixels in disk_and_mask.png.\n")
        fp.write("They must lie INSIDE the real Venus nightside, away from the terminator.\n")
        fp.write("An incorrectly measured circle or bright-side direction invalidates the detection.\n")
        fp.write("Inspect phase_model_overlay.png for phase mode; review the image angle, "
                 "--plate-scale, --cx, --cy and --phase-glare-fraction.\n")


def print_mask_stats(info: dict, folder: Path) -> None:
    print(f"Mask [{info['mode']}]: {info['safe_night_pixels']} safe night-side pixels; "
          f"{info['safe_night_percent_of_disk']:.1f}% of fitted disk")
    print(f"Saved mask preview: {folder / 'disk_and_mask.png'}")
    print(f"Saved mask statistics: {folder / 'mask_diagnostics.txt'}")
    if info['safe_night_pixels'] < 30:
        print("WARNING: Mask has fewer than 30 safe nightside pixels. "
              "Do not run candidate detection with this mask.")

def grayscale(img: np.ndarray, stretch: str = "asinh") -> np.ndarray:
    lower = float(np.percentile(img, 1))
    high = float(np.percentile(img, 99.9))
    if high <= lower:
        high = lower + 1
    z = np.clip((img - lower) / (high - lower), 0.0, 1.0)
    if stretch == "asinh":
        z = np.arcsinh(15 * z) / np.arcsinh(15)
    return np.uint8(255 * z)


def save_mask_preview(path: Path, composite: np.ndarray, in_disk: np.ndarray,
                      night: np.ndarray, bright: np.ndarray,
                      cx: float, cy: float, radius: float,
                      draw_limb: bool = True):
    g = grayscale(composite)
    base = cv2.cvtColor(g, cv2.COLOR_GRAY2BGR)
    overlay = base.copy()
    overlay[night] = (40, 230, 40)    # GREEN = accepted search region
    overlay[bright & in_disk] = (40, 40, 250)  # RED = excluded bright region
    cv2.addWeighted(base, 0.48, overlay, 0.52, 0, dst=overlay)
    if draw_limb:
        cv2.circle(overlay, (round(cx), round(cy)), round(radius), (255, 240, 0), 1)
    cv2.putText(overlay, "GREEN: searched nightside; RED: excluded bright side",
                (8, min(24, overlay.shape[0]-3)), cv2.FONT_HERSHEY_SIMPLEX,
                0.43, (255, 255, 255), 1, cv2.LINE_AA)
    cv2.imwrite(str(path), np.hstack([base, overlay]))


def block_mean(image: np.ndarray, size: int) -> np.ndarray:
    h, w = image.shape
    hh, ww = (h // size) * size, (w // size) * size
    small = image[:hh, :ww]
    if size == 1:
        return small
    return small.reshape(hh // size, size, ww // size, size).mean(axis=(1, 3))


def block_all(mask: np.ndarray, size: int) -> np.ndarray:
    h, w = mask.shape
    hh, ww = (h // size) * size, (w // size) * size
    small = mask[:hh, :ww]
    if size == 1:
        return small
    return small.reshape(hh // size, size, ww // size, size).all(axis=(1, 3))


def spatial_candidates(img: np.ndarray, allowed: np.ndarray,
                       threshold_n: float, sigma_floor: float,
                       min_excess: float, max_candidates: int):
    """3x3 local maxima, excluding center from the 8-pixel statistics."""
    if min(img.shape) < 5:
        return (np.empty(0, np.int32),) * 4
    ker = np.ones((3, 3), np.uint8)
    ker[1, 1] = 0
    # Max of the eight neighboring pixels only.
    max8 = cv2.dilate(img, ker, borderType=cv2.BORDER_CONSTANT,
                      borderValue=float('-inf'))
    sums = cv2.boxFilter(img, -1, (3, 3), normalize=False)
    sums2 = cv2.boxFilter(img * img, -1, (3, 3), normalize=False)
    mean8 = (sums - img) / 8.0
    # Sample variance of the eight surrounding pixels.
    var8 = np.maximum((sums2 - img * img - 8.0 * mean8 * mean8) / 7.0, 0.0)
    sigma8 = np.maximum(np.sqrt(var8), sigma_floor)
    excess = img - mean8
    is_peak = (img > max8) & (excess >= min_excess) & (excess >= threshold_n * sigma8)
    is_peak &= allowed
    y, x = np.nonzero(is_peak)
    if len(y) == 0:
        return y, x, np.empty(0, np.float32), np.empty(0, np.float32)
    snr = excess[y, x] / sigma8[y, x]
    if len(y) > max_candidates:
        top = np.argpartition(snr, -max_candidates)[-max_candidates:]
        y, x, snr = y[top], x[top], snr[top]
    return y, x, snr.astype(np.float32), excess[y, x].astype(np.float32)


def frame_range_good(valid: np.ndarray, j: int, t: int) -> bool:
    return j >= 0 and j + t <= len(valid) and bool(valid[j:j + t].all())


def sample_group(cube: np.memmap, j: int, t: int,
                 yy: np.ndarray, xx: np.ndarray, s: int) -> np.ndarray:
    """Mean at C binned target pixels over s*s pixels and t frames."""
    value = np.zeros(len(yy), dtype=np.float64)
    offsets = np.arange(s)
    yg = yy[:, None, None] + offsets[None, :, None]
    xg = xx[:, None, None] + offsets[None, None, :]
    for frame in range(j, j + t):
        value += cube[frame][yg, xg].mean(axis=(1, 2))
    return value / t


def transient_spatial_support(cube: np.memmap, averaged: np.ndarray,
                              ref_starts: list[int], j: int, t: int,
                              x_bin: int, y_bin: int, s: int,
                              night: np.ndarray, noise_floor: float,
                              support_sigma: float, min_pixels: int) -> tuple[bool, int]:
    """Require an 8-connected multi-pixel POSITIVE TEMPORAL residual.

    Work at full original pixel resolution, even when the initial spatial
    detection was made on s x s averages. Avoid counting unrelated positive
    pixels in the box: only the connected component containing the brightest
    temporally variable pixel in the candidate block is considered.

    Background is the per-pixel median of the same disjoint temporal reference
    windows already used by the main detector. This suppresses stationary hot
    pixels and tests the true shape of the transient rather than the shape of
    Venus's static scattered-light background.
    """
    if min_pixels <= 1:
        return True, 1

    # Original-resolution region: the binned candidate footprint plus 2 pixels
    # on each side. Include surrounding PSF pixels when the peak is near a
    # block boundary.
    h, w = averaged.shape
    bx, by = int(x_bin) * s, int(y_bin) * s
    x0, x1 = max(0, bx-2), min(w, bx+s+2)
    y0, y1 = max(0, by-2), min(h, by+s+2)
    if x1 <= x0 or y1 <= y0:
        return False, 0

    target = averaged[y0:y1, x0:x1]
    # Each reference value is a temporal average over t frames, matching target.
    refs = np.stack([
        np.asarray(cube[r:r+t, y0:y1, x0:x1]).mean(axis=0)
        for r in ref_starts
    ], axis=0)
    baseline = np.median(refs, axis=0)
    # Use local temporal scatter where available, but never below the global
    # per-pixel noise floor scaled for a t-frame average.
    sigma_ref = np.std(refs, axis=0, ddof=1)
    sigma_px = np.maximum(0.8 * sigma_ref, noise_floor / np.sqrt(t))
    residual = target - baseline
    above = (residual >= support_sigma * sigma_px) & night[y0:y1, x0:x1]

    # The triggering spatial block must contain at least one qualifying pixel.
    core_x0, core_x1 = bx-x0, min(bx+s, x1)-x0
    core_y0, core_y1 = by-y0, min(by+s, y1)-y0
    core = residual[core_y0:core_y1, core_x0:core_x1]
    core_allowed = above[core_y0:core_y1, core_x0:core_x1]
    if not core_allowed.any():
        return False, 0
    # Seed at the brightest qualifying residual within the trigger block.
    ranked_core = np.where(core_allowed, core, -np.inf)
    iy, ix = np.unravel_index(np.argmax(ranked_core), core.shape)
    seed_y, seed_x = core_y0 + int(iy), core_x0 + int(ix)

    num, labels, stats, _ = cv2.connectedComponentsWithStats(
        np.uint8(above), connectivity=8)
    label = int(labels[seed_y, seed_x])
    area = int(stats[label, cv2.CC_STAT_AREA]) if label > 0 else 0
    return area >= min_pixels, area


@dataclass
class Detection:
    start: int
    end: int
    x: float
    y: float
    spatial_scale: int
    time_scale: int
    spatial_snr: float
    temporal_snr: float
    excess: float
    value: float
    background: float
    support_pixels: int = 1

    @property
    def score(self):
        return min(self.spatial_snr, self.temporal_snr)


def search(cube: np.memmap, valid: np.ndarray, night: np.ndarray,
           spatial_scales: list[int], time_scales: list[int],
           n_sigma: float, m_sigma: float, temporal_neighbors: int,
           min_temporal_refs: int, noise_floor: float,
           min_excess: float, max_spatial_candidates: int,
           crop_x0: int, crop_y0: int,
           min_psf_pixels: int = 1,
           pixel_support_sigma: float = 2.5) -> list[Detection]:
    n, h, w = cube.shape
    results: list[Detection] = []
    allowed_scales = {}
    for s in spatial_scales:
        allowed = np.uint8(block_all(night, s))
        # The candidate and its eight neighbors must be inside the nightside.
        allowed_scales[s] = cv2.erode(allowed, np.ones((3, 3), np.uint8)) > 0
    max_neighbor = temporal_neighbors
    for t in time_scales:
        if n < t + min_temporal_refs * t:
            print(f"Skipping t={t}: insufficient frames for reference windows")
            continue
        print(f"Temporal average = {t} frame(s), spatial scales = {spatial_scales}", flush=True)
        rolling = np.zeros((h, w), dtype=np.float64)
        for q in range(t):
            rolling += cube[q]
        for j in range(n - t + 1):
            if j > 0:
                rolling += cube[j + t - 1]
                rolling -= cube[j - 1]
            if not frame_range_good(valid, j, t):
                continue
            ref_starts = [j + sign * k * t for k in range(1, max_neighbor + 1)
                          for sign in (-1, +1)]
            ref_starts = [r for r in ref_starts if frame_range_good(valid, r, t)]
            if len(ref_starts) < min_temporal_refs:
                continue
            averaged = rolling / t
            for s in spatial_scales:
                image_s = block_mean(averaged, s)
                floor_s = noise_floor / (s * np.sqrt(t))
                y, x, spatial_snr, spatial_excess = spatial_candidates(
                    image_s, allowed_scales[s], n_sigma, floor_s,
                    min_excess, max_spatial_candidates)
                if not len(y):
                    continue
                yy, xx = y * s, x * s
                ref_values = np.stack(
                    [sample_group(cube, r, t, yy, xx, s) for r in ref_starts], axis=0)
                med = np.median(ref_values, axis=0)
                mad = 1.4826 * np.median(np.abs(ref_values - med), axis=0)
                # MAD is robust against background flashes; a standard deviation
                # guard reduces false alarms from an accidentally tiny MAD.
                sample_std = np.std(ref_values, axis=0, ddof=1)
                temporal_sigma = np.maximum.reduce((
                    mad, 0.8 * sample_std, np.full_like(mad, floor_s)))
                observed = image_s[y, x]
                temporal_excess = observed - med
                temporal_snr = temporal_excess / temporal_sigma
                hits = np.nonzero((temporal_snr >= m_sigma) &
                                  (temporal_excess > min_excess))[0]
                for k in hits:
                    # Reject spatially isolated temporal outliers (hot pixels,
                    # cosmic rays, random single-pixel noise). All comparisons
                    # are made in the ORIGINAL pixel grid, not the binned one.
                    keep, support_area = transient_spatial_support(
                        cube, averaged, ref_starts, j, t,
                        int(x[k]), int(y[k]), s, night, noise_floor,
                        pixel_support_sigma, min_psf_pixels)
                    if not keep:
                        continue
                    results.append(Detection(
                        start=j, end=j + t - 1,
                        x=float(crop_x0 + (x[k] + 0.5)*s - 0.5),
                        y=float(crop_y0 + (y[k] + 0.5)*s - 0.5),
                        spatial_scale=s, time_scale=t,
                        spatial_snr=float(spatial_snr[k]),
                        temporal_snr=float(temporal_snr[k]),
                        excess=float(temporal_excess[k]),
                        value=float(observed[k]), background=float(med[k]),
                        support_pixels=int(support_area)))
            if (j + 1) % 250 == 0:
                print(f"  {j+1}/{n-t+1} starting frames, raw detections {len(results)}",
                      flush=True)
    return results


def consolidate(detections: list[Detection], frame_gap: int,
                pixel_radius: float):
    """Group multiple scales/overlapping time averages of the same flash."""
    ordered = sorted(detections, key=lambda d: (d.start, -d.score))
    clusters = []
    for d in ordered:
        match = None
        for cl in reversed(clusters):
            if cl["start"] > d.end + frame_gap:
                continue
            if cl["end"] + frame_gap < d.start:
                # Chronological order: earlier clusters cannot match either.
                break
            anchor = cl["best"]
            if np.hypot(d.x - anchor.x, d.y - anchor.y) <= pixel_radius:
                match = cl
                break
        if match is None:
            clusters.append({"start": d.start, "end": d.end,
                             "best": d, "n": 1})
        else:
            match["start"] = min(match["start"], d.start)
            match["end"] = max(match["end"], d.end)
            match["n"] += 1
            if d.score > match["best"].score:
                match["best"] = d
    return sorted(clusters, key=lambda c: c["start"])


def get_utc_times(nframes: int, file_times: list[str],
                  start: str, cadence: float | None):
    if start:
        if cadence is None or cadence <= 0:
            raise ValueError("--start-utc requires positive --cadence in seconds")
        t0 = parse_start(start)
        return [(t0 + timedelta(seconds=i*cadence)).isoformat(timespec="milliseconds")
                for i in range(nframes)]
    return file_times


def refine_peak_frame(cube: np.memmap, cluster: dict,
                      crop_x0: int, crop_y0: int) -> int:
    """Locate actual peak frame inside event interval (time averaging can blur it)."""
    d = cluster["best"]
    lo, hi = cluster["start"], cluster["end"]
    s = d.spatial_scale
    ix = int(round(d.x - crop_x0 - (s-1)/2.0))
    iy = int(round(d.y - crop_y0 - (s-1)/2.0))
    ix = int(np.clip(ix, 0, cube.shape[2] - s))
    iy = int(np.clip(iy, 0, cube.shape[1] - s))
    brightness = [float(np.mean(cube[k, iy:iy+s, ix:ix+s]))
                  for k in range(lo, hi+1)]
    return lo + int(np.argmax(brightness))


def describe_peak_time(times: list[str], start: int, end: int) -> str:
    available = [times[j] for j in range(start, end + 1) if times[j]]
    if not available:
        return ""
    try:
        dts = [parse_start(v) for v in available]
        return (dts[0] + (dts[-1] - dts[0]) / 2).isoformat(timespec="milliseconds")
    except ValueError:
        return available[len(available)//2]


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]):
    with path.open("w", newline="", encoding="utf-8") as fp:
        wr = csv.DictWriter(fp, fieldnames=fieldnames)
        wr.writeheader()
        wr.writerows(rows)


def save_event_cutout(path: Path, cube: np.memmap, valid: np.ndarray,
                      cluster: dict, crop_x0: int, crop_y0: int,
                      radius: int, night: np.ndarray,
                      in_disk: np.ndarray, bright: np.ndarray):
    d = cluster["best"]
    frame = refine_peak_frame(cube, cluster, crop_x0, crop_y0)
    current = np.asarray(cube[frame], dtype=np.float32)
    background_indices = [k for k in
                          (frame - 12, frame - 9, frame - 6,
                           frame + 6, frame + 9, frame + 12)
                          if 0 <= k < len(valid) and valid[k]
                          and not (d.start <= k <= d.end)]
    background = np.median(cube[background_indices], axis=0) if background_indices else np.zeros_like(current)
    diff = current - background
    cx = int(round(d.x - crop_x0))
    cy = int(round(d.y - crop_y0))
    h, w = current.shape
    x0, x1 = max(0, cx - radius), min(w, cx + radius + 1)
    y0, y1 = max(0, cy - radius), min(h, cy + radius + 1)
    before = current[y0:y1, x0:x1]
    delta = diff[y0:y1, x0:x1]
    image_a = cv2.cvtColor(grayscale(before), cv2.COLOR_GRAY2BGR)
    max_abs = max(float(np.percentile(np.abs(delta), 99.5)), 1e-6)
    image_b = np.uint8(np.clip(delta / max_abs * 125 + 127, 0, 255))
    image_b = cv2.applyColorMap(image_b, cv2.COLORMAP_TURBO)

    # Third panel: actual search mask. BGR colors:
    # green = searched nightside, red = excluded bright/glare,
    # orange = excluded region inside disk, gray = outside disk.
    colored = np.empty_like(image_a)
    colored[:] = (100, 100, 100)
    local_disk = in_disk[y0:y1, x0:x1]
    local_bright = bright[y0:y1, x0:x1]
    local_night = night[y0:y1, x0:x1]
    colored[local_disk] = (0, 165, 255)
    colored[local_bright & local_disk] = (0, 0, 220)
    colored[local_night] = (0, 200, 0)
    image_c = cv2.addWeighted(image_a, 0.5, colored, 0.5, 0)

    loc = (cx - x0, cy - y0)
    for panel in (image_a, image_b, image_c):
        cv2.drawMarker(panel, loc, (255, 255, 255), cv2.MARKER_CROSS, 12, 1)
        cv2.circle(panel, loc, 5, (255, 255, 255), 1)
    # 3x nearest-neighbor expansion keeps individual pixels visible.
    panel = np.hstack([image_a, image_b, image_c])
    panel = cv2.resize(panel, (panel.shape[1]*3, panel.shape[0]*3),
                       interpolation=cv2.INTER_NEAREST)
    cv2.imwrite(str(path), panel)


@contextmanager
def managed_memmap(path: Path, shape: tuple[int, int, int]):
    """Close mmap BEFORE TemporaryDirectory tries to unlink it on Windows.

    This also runs after exceptions in alignment or masking.
    """
    arr = np.memmap(path, dtype="float32", mode="w+", shape=shape)
    try:
        yield arr
    finally:
        try:
            arr.flush()
        finally:
            # Relying on 'del arr' alone does not release Windows file locks
            # while the outer main() function still holds 'cube'.
            if getattr(arr, "_mmap", None) is not None:
                arr._mmap.close()


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, type=Path,
                   help="Folder of BMP/FITS frames OR an individual .avi video file")
    p.add_argument("--output", required=True, type=Path, help="Output folder")
    p.add_argument("--spatial-scales", default="1,2,3", help="Square block sizes, e.g. 1,2,3,4")
    p.add_argument("--time-scales", default="1,2,4,16,32",
                   help="Numbers of frames averaged [1,2,4,16,32]; for 30 ms video also try 64")
    p.add_argument("--N", type=float, default=1.5,
                   help="Spatial sigma threshold [1.5]; low to admit seeing-blurred PSFs")
    p.add_argument("--M", type=float, default=7.0,
                   help="Temporal sigma threshold [7.0]; the primary false-alarm rejection")
    p.add_argument("--neighbors", type=int, default=6,
                   help="Number of disjoint background windows on each side [6]")
    p.add_argument("--min-references", type=int, default=6,
                   help="Minimum usable temporal comparison windows [6]")
    p.add_argument("--preview-only", action="store_true",
                   help="Estimate disk/mask from FIRST frame; save PNG and exit without processing video")
    p.add_argument("--disk-threshold-fraction", type=float, default=0.12,
                   help="Detection level for automatic disk geometry [0.12]; try 0.25-0.40 for heavy glare")
    p.add_argument("--mask-mode", choices=["intensity", "geometric", "phase", "manual"], default="intensity",
                   help="Mask [intensity]; phase=automatic Keplerian fit; manual=interactive Keplerian fit")
    p.add_argument("--plate-scale", type=float, default=None,
                   help="Image scale in arcsec/pixel, REQUIRED for --mask-mode phase or manual")
    p.add_argument("--manual-config", type=Path, default=None,
                   help="Load an accepted JSON manual fit; skips GUI and requires the same first frame")
    p.add_argument("--observation-utc", default="",
                   help="Override date/time used to calculate Venus phase (UTC ISO timestamp)")
    p.add_argument("--phase-image-angle", type=float, default=None,
                   help="Manual Sunward image direction in degrees, clockwise from +x to +y; "
                        "otherwise fit crescent orientation automatically")
    p.add_argument("--phase-fit-profile", choices=["auto", "binary", "lambert"], default="auto",
                   help="Analytical silhouette fit: auto selects binary or Lambertian photometric profile [auto]")
    p.add_argument("--phase-glare-fraction", type=float, default=0.65,
                   help="Additional intensity veto above this fraction of the bright peak "
                        "in PHASE mode [0.65], 1 disables veto")
    p.add_argument("--phase-min-fit-score", type=float, default=0.30,
                   help="Minimum template-correlation quality required for full PHASE search [0.30]")
    p.add_argument("--night-cap", type=float, default=0.25,
                   help="In geometric mode, require projection >= this radius fraction into dark hemisphere [0.25]")
    p.add_argument("--cx", type=float, default=None, help="Disk center x in original pixels")
    p.add_argument("--cy", type=float, default=None, help="Disk center y in original pixels")
    p.add_argument("--radius", type=float, default=None, help="Disk radius in original pixels")
    p.add_argument("--crop-pad", type=int, default=24, help="Pixels around disk to retain [24]")
    p.add_argument("--bright-margin", type=int, default=10,
                   help="Dilation of bright-side/terminator mask in pixels [10]")
    p.add_argument("--limb-margin", type=float, default=5.0,
                   help="Exclude outer this many pixels of planet [5]")
    p.add_argument("--bright-fraction", type=float, default=0.03,
                   help="Bright mask cutoff fraction of disk peak [0.03]")
    p.add_argument("--align-size", type=int, default=512,
                   help="Longest side of downsampled frame used for alignment [512]")
    p.add_argument("--registration-method", choices=["template", "phase"],
                   default="template", help="Bright-crescent registration method [template]")
    p.add_argument("--min-registration-response", type=float, default=None,
                   help="Skip uncertain registrations (default: 0.7 template, 0.03 phase)")
    p.add_argument("--noise-floor", type=float, default=None,
                   help="Minimum per-pixel rms noise in input units (auto by default)")
    p.add_argument("--min-excess", type=float, default=0.0,
                   help="Minimum positive excess in image units [0]")
    p.add_argument("--min-psf-pixels", type=int, default=2,
                   help="Require this many 8-connected positive temporal-residual pixels [2]; 1 disables")
    p.add_argument("--pixel-support-sigma", type=float, default=2.5,
                   help="Per-pixel temporal-residual significance used to count PSF pixels [2.5]")
    p.add_argument("--max-spatial-candidates", type=int, default=300,
                   help="Keep top this many peaks per frame/scale [300]")
    p.add_argument("--merge-frames", type=int, default=3,
                   help="Merge events separated by up to this many frames [3]")
    p.add_argument("--merge-radius", type=float, default=8.0,
                   help="Merge detections this many pixels apart [8]")
    p.add_argument("--start-utc", default="",
                   help="First frame time in ISO UTC (e.g. 2026-07-05T18:00:00)")
    p.add_argument("--cadence", type=float, default=None,
                   help="Frame interval in seconds (AVI: overrides metadata FPS; "
                        "BMP: required with --start-utc)")
    p.add_argument("--max-cutouts", type=int, default=60,
                   help="Max number of event diagnostic PNGs [60]")
    args = p.parse_args()

    s_scales = sorted(set(int(x) for x in args.spatial_scales.split(",")))
    t_scales = sorted(set(int(x) for x in args.time_scales.split(",")))
    if min(s_scales + t_scales) < 1:
        p.error("All spatial and temporal scales must be >=1")
    if not 1 <= args.min_references <= 2*args.neighbors:
        p.error("--min-references must be between 1 and twice --neighbors")
    if args.N <= 0 or args.M <= 0 or args.crop_pad < 0 or args.bright_margin < 0:
        p.error("N, M must be positive; margins must be non-negative")
    if args.min_psf_pixels < 1 or args.pixel_support_sigma <= 0:
        p.error("--min-psf-pixels must be >=1; --pixel-support-sigma must be >0")
    if not (0.01 < args.disk_threshold_fraction < 0.95):
        p.error("--disk-threshold-fraction must be between 0.01 and 0.95")
    if not (0.0 < args.bright_fraction < 1.0):
        p.error("--bright-fraction must be between 0 and 1")
    if not (0.0 <= args.night_cap < 0.95):
        p.error("--night-cap must be between 0 and 0.95")
    if args.limb_margin < 0:
        p.error("--limb-margin must be nonnegative")
    if not 0 < args.phase_glare_fraction <= 1:
        p.error("--phase-glare-fraction must be in (0,1]")
    if not -1 <= args.phase_min_fit_score <= 1:
        p.error("--phase-min-fit-score must lie between -1 and 1")
    if args.mask_mode in ('phase','manual') and (args.plate_scale is None or args.plate_scale <= 0):
        p.error("--mask-mode phase/manual requires --plate-scale in arcsec/pixel")
    if args.manual_config is not None and args.mask_mode != 'manual':
        p.error("--manual-config can only be used with --mask-mode manual")
    if args.cadence is not None and (not np.isfinite(args.cadence) or args.cadence <= 0):
        p.error("--cadence must be a positive number of seconds")
    manual = (args.cx, args.cy, args.radius)
    if args.mask_mode in ('phase','manual'):
        if (args.cx is None) != (args.cy is None):
            p.error("In phase mode provide BOTH --cx and --cy, or neither")
    elif sum(x is not None for x in manual) not in (0, 3):
        p.error("Provide all of --cx --cy --radius, or none of them")

    source = open_frame_source(args.input)
    frame_cadence = args.cadence
    if source.is_video:
        if frame_cadence is None and source.fps is not None:
            frame_cadence = 1.0 / source.fps
        fps_text = f"{source.fps:.6g}" if source.fps is not None else "unavailable"
        dt_text = f"{frame_cadence:.6f} s" if frame_cadence else "unknown"
        print(f"AVI input: {source.nframes} frames; AVI metadata FPS={fps_text}; "
              f"adopted frame cadence={dt_text}")
        if args.cadence is not None:
            print("AVI frame timing overridden by --cadence")
        print("NOTICE: compressed AVI may introduce false flashes; "
              "OpenCV may output 8-bit images even from higher-bit-depth sources")
    if args.start_utc and frame_cadence is None:
        p.error("Cannot calculate UTC times: supply --cadence, "
                "or use an AVI with valid FPS metadata")
    args.output.mkdir(parents=True, exist_ok=True)
    n = source.nframes
    first = source.first
    H, W = first.shape
    phase_ephem=None
    phase_fit_params=None
    if args.mask_mode in ('phase','manual'):
        if args.observation_utc:
            phase_time=parse_start(args.observation_utc)
        elif args.start_utc:
            phase_time=parse_start(args.start_utc)
        elif source.files is not None:
            _, first_date_obs = read_image(source.files[0])
            if not first_date_obs:
                p.error('For phase/manual mode supply --observation-utc or '
                        '--start-utc (no valid DATE-OBS in the first image)')
            phase_time=parse_start(first_date_obs)
        else:
            p.error('For AVI phase/manual mode supply --start-utc or --observation-utc')
        phase_ephem=venus_phase_geometry(phase_time,args.plate_scale)
        radius=float(args.radius) if args.radius is not None else phase_ephem['radius_px']
        print(f"Keplerian Venus phase: {phase_ephem['illuminated_fraction']*100:.1f}% illuminated; "
              f"phase angle {phase_ephem['phase_angle_deg']:.2f} deg; "
              f"diameter {phase_ephem['diameter_arcsec']:.2f} arcsec; "
              f"predicted radius {phase_ephem['radius_px']:.2f} px")
        if args.mask_mode=='manual':
            if args.manual_config is not None:
                geometry,payload=read_manual_config(args.manual_config,source,first,
                                                    phase_ephem,phase_time)
                print(f"Using saved manual geometry: {args.manual_config}")
                # Keep the safety exclusions identical to those at manual approval.
                if (args.bright_margin != payload['bright_margin_px'] or
                    abs(args.limb_margin-payload['limb_margin_px'])>1e-8 or
                    abs(args.phase_glare_fraction-payload['phase_glare_fraction'])>1e-8):
                    raise ValueError('Mask margins or glare fraction differ from the saved '
                                     'geometry; rerun with matching parameters.')
                screen,_,_ = render_manual_fit(first,phase_ephem,geometry,
                    args.bright_margin,args.limb_margin,args.phase_glare_fraction)
            else:
                initial_cx=float(args.cx) if args.cx is not None else (W-1)/2
                initial_cy=float(args.cy) if args.cy is not None else (H-1)/2
                initial_angle=(args.phase_image_angle if args.phase_image_angle is not None
                               else manual_initial_angle(first,initial_cx,initial_cy))
                initial=dict(center_x=initial_cx,center_y=initial_cy,radius_px=radius,
                             sun_angle_image_deg=initial_angle)
                print(f'Manual mask initial center=({initial_cx:.1f}, {initial_cy:.1f}), '
                      f'radius={radius:.1f} px, approximate sun angle={initial_angle:.1f} deg')
                geometry,screen=manual_mask_gui(first,phase_ephem,initial,
                    args.bright_margin,args.limb_margin,args.phase_glare_fraction)
            cx=float(geometry['center_x']); cy=float(geometry['center_y'])
            radius=float(geometry['radius_px'])
            phase_image_angle=float(geometry['sun_angle_image_deg'])%360.0
            phase_fit_params={**geometry,'score':None,
                'fitted_brightness_profile':'manual-visual-fit',
                'center_fixed':True,'sun_angle_fixed':True,
                'fitting_method':'human visual alignment to analytical Keplerian phase'}
            write_manual_config(args.output/'manual_mask_geometry.json',source,first,
                phase_ephem,geometry,phase_time,args.bright_margin,
                args.limb_margin,args.phase_glare_fraction)
            cv2.imwrite(str(args.output/'manual_fit_preview.png'),screen)
            print(f'Manual fit accepted: center=({cx:.1f},{cy:.1f}), '
                  f'radius={radius:.1f} px, sun angle={phase_image_angle:.1f} deg')
            if abs(radius/phase_ephem['radius_px']-1) > 0.20:
                print('WARNING: manual radius differs from Keplerian predicted radius '
                      'by more than 20%. Bright optical glare is NOT the physical limb. '
                      'Check the actual AVI pixel scale and geometry.')
            print(f"Saved geometry: {args.output / 'manual_mask_geometry.json'}")
        else:
            if args.radius is not None:
                print('NOTICE: --radius overrides predicted physical radius in phase mode.')
            fixed_center=None if args.cx is None else (float(args.cx),float(args.cy))
            phase_fit_params=phase_fit(first,radius,phase_ephem['phase_angle_deg'],
                                       fixed_center,args.phase_image_angle,
                                       args.phase_fit_profile)
            cx=float(phase_fit_params['center_x'])
            cy=float(phase_fit_params['center_y'])
            phase_image_angle=float(phase_fit_params['sun_angle_image_deg'])
            print(f"Fitted image center=({cx:.1f},{cy:.1f}) px; bright-side image angle="
                  f"{phase_image_angle:.1f} deg; fit correlation={phase_fit_params['score']:.3f}; "
                  f"profile={phase_fit_params['fitted_brightness_profile']}")
        write_phase_geometry(args.output / 'phase_geometry.txt',phase_ephem,phase_fit_params)
    elif manual[0] is None:
        cx, cy, radius = estimate_disk(first, args.disk_threshold_fraction)
        print(f"Estimated Venus disk: center=({cx:.1f},{cy:.1f}), radius={radius:.1f} px")
    else:
        cx, cy, radius = manual
        print(f"User-defined Venus disk: center=({cx:.1f},{cy:.1f}), radius={radius:.1f} px")
    if radius <= 0 or not (0 <= cx < W and 0 <= cy < H):
        p.error("Invalid disk center or radius")
    def make_masks(image: np.ndarray, center_x: float, center_y: float):
        if args.mask_mode in ('phase','manual'):
            assert phase_ephem is not None
            result=phase_masks(image,center_x,center_y,radius,phase_ephem,
                               phase_image_angle,args.bright_margin,
                               args.limb_margin,args.phase_glare_fraction)
            if args.mask_mode=='manual':
                result[-1]['mode']='manual'
            return result
        return build_masks(image,center_x,center_y,radius,args.bright_margin,
                           args.limb_margin,args.bright_fraction,
                           args.mask_mode,args.night_cap)

    def emit_mask_diagnostics(image, inside, nightside, excluded,
                              center_x: float, center_y: float, info: dict,
                              preview: bool):
        save_mask_preview(args.output / 'disk_and_mask.png',image,
                          inside,nightside,excluded,center_x,center_y,radius,
                          draw_limb=(args.mask_mode not in ('phase','manual')))
        if args.mask_mode in ('phase','manual'):
            save_phase_overlay(args.output / 'phase_model_overlay.png', image,
                               inside,nightside,center_x,center_y,radius,
                               phase_ephem['phase_angle_deg'],phase_image_angle)
            info.update(phase_ephem)
            info.update({'phase_fit_score':phase_fit_params['score'],
                         'phase_fit_image_angle_deg':phase_image_angle,
                         'phase_fit_center_x':cx,'phase_fit_center_y':cy,
                         'phase_fit_profile':phase_fit_params['fitted_brightness_profile']})
        save_mask_diagnostics(args.output / 'mask_diagnostics.txt', info, preview=preview)
        print_mask_stats(info,args.output)
        if args.mask_mode == 'phase' and phase_fit_params['score'] < args.phase_min_fit_score:
            print('WARNING: Low phase-fit correlation. The orientation/center may be unreliable.')
            if not preview:
                raise RuntimeError('Phase-profile fit score is below --phase-min-fit-score. '
                                   'Inspect phase_model_overlay.png and set --cx --cy '
                                   'and/or --phase-image-angle; you may lower '
                                   '--phase-min-fit-score after verifying the geometry.')

    # Fast geometry QA: do not allocate the huge aligned-frame data file when
    # the first-frame disk/mask already looks implausible.
    if args.preview_only:
        in_disk, night, bright, bg, info = make_masks(first,cx,cy)
        emit_mask_diagnostics(first,in_disk,night,bright,cx,cy,info,True)
        print("PREVIEW ONLY: no video frames were aligned and no flashes were searched.")
        return
    x0 = max(0, int(np.floor(cx - radius - args.crop_pad)))
    x1 = min(W, int(np.ceil(cx + radius + args.crop_pad)) + 1)
    y0 = max(0, int(np.floor(cy - radius - args.crop_pad)))
    y1 = min(H, int(np.ceil(cy + radius + args.crop_pad)) + 1)
    crop = (x0, y0, x1, y1)
    hc, wc = y1 - y0, x1 - x0
    if min(hc, wc) < 10:
        p.error("Empty or too-small crop. Check disk geometry")
    if min(hc, wc) // max(s_scales) < 5:
        p.error("Some spatial scales are too large for the cropped image")
    need_bytes = n * hc * wc * np.dtype("float32").itemsize
    free_bytes = shutil.disk_usage(args.output).free
    print(f"{n} frames, original {W}x{H}, stored Venus crop {wc}x{hc}; "
          f"temporary disk use ~{need_bytes / 1024**3:.2f} GiB")
    if need_bytes > 0.9 * free_bytes:
        raise RuntimeError("Insufficient free space for the aligned image sequence; "
                           "free space in the output drive or process shorter batches")

    scale = min(1.0, args.align_size / max(W, H))
    smallsize = (max(16, int(round(W * scale))), max(16, int(round(H * scale))))
    ref_feature = register_feature(first, smallsize)
    hann = cv2.createHanningWindow(smallsize, cv2.CV_32F)
    template, template_x0, template_y0 = prepare_template(
        ref_feature, cx, cy, radius, smallsize[0]/W, smallsize[1]/H)
    if args.min_registration_response is None:
        args.min_registration_response = 0.7 if args.registration_method == "template" else 0.03
    print(f"Registering with {args.registration_method} method "
          f"(minimum score {args.min_registration_response:g})")
    valid = np.zeros(n, dtype=bool)
    align_rows: list[dict] = []
    file_times: list[str] = []

    # An on-disk array avoids keeping every video frame in RAM.
    with tempfile.TemporaryDirectory(prefix="venus_alfi_", dir=str(args.output)) as tmp:
        path_cube = Path(tmp) / "aligned.dat"
        with managed_memmap(path_cube, shape=(n, hc, wc)) as cube:
            for i, (raw, date_obs) in enumerate(source.iter_frames()):
                if raw.shape != (H, W):
                    raise ValueError(f"Image shape changed in frame {i}: "
                                     f"{raw.shape}, expected {(H,W)}")
                file_times.append(date_obs)
                if i == 0:
                    dx, dy, response = 0.0, 0.0, 1.0
                else:
                    feat = register_feature(raw, smallsize)
                    if args.registration_method == "template":
                        dx, dy, response = template_shift(
                            feat, template, template_x0, template_y0,
                            smallsize[0]/W, smallsize[1]/H)
                    else:
                        (shift_x, shift_y), response = cv2.phaseCorrelate(
                            ref_feature, feat, hann)
                        dx = shift_x / (smallsize[0]/W)
                        dy = shift_y / (smallsize[1]/H)
                # If the Venus disk falls outside the original frame in this image,
                # do not analyze: interpolation would create artificial flashes.
                disk_fits = (cx + dx - radius + args.limb_margin >= 0 and
                             cx + dx + radius - args.limb_margin < W and
                             cy + dy - radius + args.limb_margin >= 0 and
                             cy + dy + radius - args.limb_margin < H)
                valid[i] = bool((i == 0 or response >= args.min_registration_response)
                                and disk_fits)
                if valid[i]:
                    cube[i] = crop_registered(raw, dx, dy, crop)
                else:
                    cube[i] = 0.0
                align_rows.append(dict(frame=i, filename=source.frame_name(i),
                                       dx_pixels=f"{dx:.3f}", dy_pixels=f"{dy:.3f}",
                                       registration_response=f"{response:.4f}",
                                       usable=int(valid[i]), utc=date_obs))
                if (i + 1) % 250 == 0 or i == n - 1:
                    print(f"Aligned {i+1}/{n} frames; {valid[:i+1].sum()} usable", flush=True)
            cube.flush()
            times = get_utc_times(n, file_times, args.start_utc, frame_cadence)
            for i in range(n):
                align_rows[i]["utc"] = times[i]
                if source.is_video:
                    align_rows[i]["video_time_s"] = (
                        f"{i * frame_cadence:.6f}" if frame_cadence else "")
            alignment_fields = ["frame", "filename", "dx_pixels", "dy_pixels",
                                "registration_response", "usable", "utc"]
            if source.is_video:
                alignment_fields.append("video_time_s")
            write_csv(args.output / "alignment.csv", align_rows,
                      alignment_fields)
            if np.count_nonzero(valid) < max(12, 2 * args.neighbors + 1):
                raise RuntimeError("Too few usable registered frames. Check alignment.csv, "
                                   "disk estimate, or lower --min-registration-response")
            sample_indices = np.flatnonzero(valid)
            sample_indices = sample_indices[np.unique(
                np.linspace(0, len(sample_indices)-1, min(31, len(sample_indices))).astype(int))]
            composite = np.median(cube[sample_indices], axis=0)
            in_disk, night, bright, bg, info = make_masks(composite,cx-x0,cy-y0)
            emit_mask_diagnostics(composite,in_disk,night,bright,
                                  cx-x0,cy-y0,info,False)
            if info["safe_night_pixels"] < 30:
                raise RuntimeError(
                    "Too few safe night-side pixels. The PREVIEW and diagnostics were saved. "
                    "First inspect disk_and_mask.png: if the circle is wrong, use "
                    "--cx --cy --radius or --disk-threshold-fraction. "
                    "If the circle is right but scattered light covers the nightside, "
                    "try --mask-mode phase --plate-scale ... (with UTC time), "
                    "or review --phase-glare-fraction in preview mode. "
                    "Do not blindly lower the masking safety margins.")

            if args.noise_floor is not None:
                noise_floor = args.noise_floor
            else:
                # Measure the noise in pairwise registered-frame differences in the
                # searched region. Fixed pattern and planetary texture subtract out.
                measures = []
                for a, b in zip(sample_indices[:-1], sample_indices[1:]):
                    d = cube[a] - cube[b]
                    sig = robust_sigma(d[night]) / np.sqrt(2.0)
                    if sig > 0:
                        measures.append(sig)
                noise_floor = float(np.median(measures)) if measures else 0.5
                noise_floor = max(noise_floor, 0.5)  # minimum quantization/read-noise guard
            if noise_floor <= 0:
                p.error("--noise-floor must be positive")
            print(f"Adopted per-pixel RMS noise floor: {noise_floor:.3f} image units")
            print(f"Minimum PSF connected pixels: {args.min_psf_pixels}; "
                  f"support threshold: {args.pixel_support_sigma:.2f} sigma")
            print("Inspect disk_and_mask.png before interpreting the detections!")

            detections = search(cube, valid, night, s_scales, t_scales,
                                args.N, args.M, args.neighbors,
                                args.min_references, noise_floor,
                                args.min_excess, args.max_spatial_candidates,
                                x0, y0, args.min_psf_pixels,
                                args.pixel_support_sigma)
            clusters = consolidate(detections, args.merge_frames, args.merge_radius)

            rows = []
            for idx, cl in enumerate(clusters, 1):
                d = cl["best"]
                peak_idx = refine_peak_frame(cube, cl, x0, y0)
                row = dict(event_id=f"V{idx:04d}",
                           first_frame=cl["start"], last_frame=cl["end"],
                           peak_frame=peak_idx,
                           peak_filename=source.frame_name(peak_idx),
                           utc=times[peak_idx],
                           x_reference_pixels=f"{d.x:.2f}",
                           y_reference_pixels=f"{d.y:.2f}",
                           spatial_bin=d.spatial_scale,
                           temporal_average_frames=d.time_scale,
                           spatial_sigma=f"{d.spatial_snr:.2f}",
                           temporal_sigma=f"{d.temporal_snr:.2f}",
                           signal_excess=f"{d.excess:.3f}",
                           original_pixel_value=f"{d.value:.3f}",
                           local_temporal_baseline=f"{d.background:.3f}",
                           combined_score=f"{d.score:.2f}",
                           raw_detections_merged=cl["n"],
                           support_pixels=d.support_pixels)
                if source.is_video:
                    row["video_time_s"] = (
                        f"{peak_idx * frame_cadence:.6f}" if frame_cadence else "")
                rows.append(row)
            fields = ["event_id", "first_frame", "last_frame", "peak_frame",
                      "peak_filename", "utc", "x_reference_pixels", "y_reference_pixels",
                      "spatial_bin", "temporal_average_frames", "spatial_sigma",
                      "temporal_sigma", "signal_excess", "original_pixel_value",
                      "local_temporal_baseline", "combined_score", "raw_detections_merged",
                      "support_pixels"]
            if source.is_video:
                fields.append("video_time_s")
            write_csv(args.output / "venus_candidates.csv", rows, fields)
            with (args.output / "venus_candidates.txt").open("w", encoding="utf-8") as fp:
                fp.write("VENUS NIGHT-SIDE FLASH CANDIDATES (not confirmed impacts)\n")
                fp.write(f"Spatial N={args.N:g} sigma; temporal M={args.M:g} sigma; "
                         f"spatial scales={s_scales}; frame averages={t_scales}\n")
                fp.write(f"Registered frames: {valid.sum()} of {n}; "
                         f"local noise floor: {noise_floor:.3f} image units\n\n")
                fp.write(f"{'ID':<7s} {'UTC':<27s} {'Frame':>7s} "
                         f"{'X(px)':>8s} {'Y(px)':>8s} {'Bin':>5s} "
                         f"{'t':>4s} {'S-sig':>7s} {'T-sig':>7s} {'Excess':>10s}\n")
                for row in rows:
                    fp.write(f"{row['event_id']:<7s} {row['utc']:<27s} "
                             f"{row['peak_frame']:>7d} {row['x_reference_pixels']:>8s} "
                             f"{row['y_reference_pixels']:>8s} {row['spatial_bin']:>5d} "
                             f"{row['temporal_average_frames']:>4d} "
                             f"{row['spatial_sigma']:>7s} {row['temporal_sigma']:>7s} "
                             f"{row['signal_excess']:>10s}\n")
            if rows:
                cutout_folder = args.output / "event_cutouts"
                cutout_folder.mkdir(exist_ok=True)
                top = sorted(range(len(clusters)),
                             key=lambda k: clusters[k]["best"].score, reverse=True)
                for idx in top[:max(0, args.max_cutouts)]:
                    try:
                        save_event_cutout(cutout_folder / f"V{idx+1:04d}.png",
                                          cube, valid, clusters[idx], x0, y0, 30,
                                          night, in_disk, bright)
                    except Exception as exc:
                        print(f"Could not make event cutout for V{idx+1:04d}: {exc}")
            print(f"Finished: {len(detections)} raw triggers consolidated into "
                  f"{len(rows)} candidates")
            print(f"Results: {args.output / 'venus_candidates.csv'}")
            print(f"Diagnostics: {args.output / 'disk_and_mask.png'}")
            print("CAUTION: these are statistical candidates, NOT confirmed impacts. "
                  "Check cutouts and independent observations.")


if __name__ == "__main__":
    main()
