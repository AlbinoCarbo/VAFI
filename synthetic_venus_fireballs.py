#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate a synthetic BMP sequence of Venus as seen from Earth, with:
- analytically computed illuminated phase (simple Keplerian orbits),
- Venus bright silhouette on black background,
- synthetic fireballs injected on the dark side,
       Generate a normalized temporal brightness profile for a fireball lasting n frames.
       Supported profiles:
           triangular : brightness increases and then decreases
           gaussian   : Gaussian brightness evolution
           flat       : constant brightness
       The maximum brightness is always normalized to 1.0.
       Works correctly for any duration >= 1 frame.
- optional frame-to-frame jitter and Gaussian noise.

The parameter --noise-sigma represents the standard deviation of random fluctuations in pixel intensity, 
simulating the noise present in astronomical images acquired with a CMOS or CCD camera.
The parameter --noise-sigma is the standard deviation of the random Gaussian noise added to each pixel, expressed in image 
intensity units (gray levels for 8-bit BMP images). It does not describe the shape of the fireball.

Consider a synthetic fireball with a Gaussian PSF and a peak amplitude of 100 intensity units.
Its approximate peak-pixel signal-to-noise ratio, considering only the added Gaussian noise, is
\[\mathrm{SNR}_{\mathrm{peak}}\approx\frac{A_{\mathrm{peak}}}{\sigma_{\mathrm{noise}}}\]

Noise σ		Approximate peak SNR	Expected difficulty
1			100						Very easy
3			33						Easy
5			20						Moderate
10			10						More challenging
20			5						Difficult

Parameter for the PSF of the fireball:
--min-sigma 1.0 and --max-sigma 2.5: These define the range of standard deviations of the spatial Gaussian PSF 
of the synthetic fireballs, expressed in pixels. Each fireball is assigned a random σ within this range.

Outputs:
- BMP frames
- truth_events.csv
- sequence_info.txt

Examples:
py synthetic_venus_fireballs.py --output .\Synthetic_Venus_Test --date 2026-07-05T17:34:00 --diameter 188 --width 448 --height 448 --sun-angle 315 --n-frames 300 --n-events 8 --noise-sigma 2 --jitter-sigma 0.8

Meaning of the main parameters
- --date  UTC observing date. The script computes the Venus phase for that date.
- --diameter  Diameter of Venus in pixels.
- --sun-angle  Direction of the illuminated side in the image:
  - 0° = Sun toward the right
  - 90° = Sun toward the bottom
  - 180° = Sun toward the left
  - 270° = Sun toward the top
- --jitter-sigma  Adds frame-to-frame shifts of Venus, useful to test the alignment part of your detector.
- --n-events  Number of synthetic fireballs injected randomly on the dark side.

If you want to define the fireballs yourself, create a CSV like this:

event_id,start_frame,duration,peak_amplitude,sigma_px,x_rel,y_rel,auto_place
1,40,2,90,1.5,, ,True
2,110,3,120,2.0,-40,25,False
3,180,1,70,1.2,30,50,False

And run:

py synthetic_venus_fireballs.py --output .\Synthetic_Venus_Test --date 2026-07-05T17:34:00 --diameter 188 --width 448 --height 448 --events-file .\my_events.csv

With:
x_rel, y_rel are in pixels relative to the Venus center.
If auto_place = True, the script automatically places the event in the dark hemisphere.

Albino Carbognani, INAF-OAS
October 29, 2026
"""

from __future__ import annotations
import os
import csv
import math
import argparse
from dataclasses import dataclass
from datetime import datetime, timezone, timedelta
from pathlib import Path

import numpy as np
import cv2


# ============================================================
# BASIC ASTRONOMY
# ============================================================

DEG = np.pi / 180.0
KGAUSS = 0.01720209895          # au^(3/2)/day
EPS_J2000 = 23.439291111 * DEG  # obliquity


ORBITS = {
    "earth": dict(
        a=1.00000011,
        e=0.01671022,
        i=0.00005,
        Omega=-11.26064,
        omega=102.94719 - (-11.26064),
        M=100.46435 - 102.94719,
        epoch=2451545.0
    ),
    "venus": dict(
        a=0.72333199,
        e=0.00677323,
        i=3.39471,
        Omega=76.68069,
        omega=131.53298 - 76.68069,
        M=181.97973 - 131.53298,
        epoch=2451545.0
    )
}


def jd_from_datetime(dt: datetime) -> float:
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    y, m = dt.year, dt.month
    d = dt.day + (dt.hour + dt.minute / 60.0 + dt.second / 3600.0 + dt.microsecond / 3.6e9) / 24.0
    if m <= 2:
        y -= 1
        m += 12
    A = y // 100
    B = 2 - A + A // 4
    jd = int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + B - 1524.5
    return jd


def solve_kepler(M, e, tol=1e-12):
    M = np.mod(M, 2*np.pi)
    E = M if e < 0.8 else np.pi
    for _ in range(100):
        dE = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        E -= dE
        if abs(dE) < tol:
            break
    return E


def rot_z(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0],
                     [s,  c, 0],
                     [0,  0, 1]], dtype=float)


def rot_x(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[1, 0, 0],
                     [0, c, -s],
                     [0, s,  c]], dtype=float)


def heliocentric_position(el, jd):
    a, e = el["a"], el["e"]
    inc = el["i"] * DEG
    Om = el["Omega"] * DEG
    om = el["omega"] * DEG
    M0 = el["M"] * DEG

    n = KGAUSS / (a ** 1.5)
    M = M0 + n * (jd - el["epoch"])
    E = solve_kepler(M, e)

    x = a * (np.cos(E) - e)
    y = a * np.sqrt(1 - e * e) * np.sin(E)
    r_pf = np.array([x, y, 0.0])

    R = rot_z(Om) @ rot_x(inc) @ rot_z(om)
    return R @ r_pf


def venus_phase_and_distance(obs_dt: datetime):
    """
    Compute Sun-Venus-Earth phase angle and Earth-Venus distance
    using simple Keplerian heliocentric orbits.
    """
    jd = jd_from_datetime(obs_dt)
    rE = heliocentric_position(ORBITS["earth"], jd)
    rV = heliocentric_position(ORBITS["venus"], jd)

    sun_from_venus = -rV
    earth_from_venus = rE - rV

    cang = np.dot(sun_from_venus, earth_from_venus) / (
        np.linalg.norm(sun_from_venus) * np.linalg.norm(earth_from_venus)
    )
    cang = np.clip(cang, -1.0, 1.0)
    alpha = np.degrees(np.arccos(cang))  # Sun-Venus-Earth phase angle

    illum_frac = 0.5 * (1.0 + np.cos(np.radians(alpha)))
    delta = np.linalg.norm(rV - rE)  # Earth-Venus distance in au

    return alpha, illum_frac, delta


# ============================================================
# IMAGE SYNTHESIS
# ============================================================

@dataclass
class FireballEvent:
    event_id: int
    start_frame: int
    duration: int
    peak_amplitude: float
    sigma_px: float
    x_rel: float | None = None  # relative to center, pixels
    y_rel: float | None = None
    auto_place: bool = True


def analytical_masks(
    width: int,
    height: int,
    cx: float,
    cy: float,
    radius: float,
    phase_angle_deg: float,
    sun_angle_img_deg: float
):
    """
    Return:
      disk_mask          full Venus disk
      bright_mask        illuminated region
      dark_mask          nightside region
    """
    yy, xx = np.indices((height, width), dtype=float)
    x = (xx - cx) / radius
    y = (yy - cy) / radius

    disk = (x*x + y*y) <= 1.0

    # Rotate coordinates so +u points toward the Sun in the image plane
    th = np.radians(sun_angle_img_deg)
    u = x * np.cos(th) + y * np.sin(th)
    v = -x * np.sin(th) + y * np.cos(th)

    # Terminator:
    # u_term = -cos(alpha) * sqrt(1 - v^2)
    alpha = np.radians(phase_angle_deg)

    inside_strip = np.abs(v) <= 1.0
    u_term = np.zeros_like(u)
    u_term[inside_strip] = -np.cos(alpha) * np.sqrt(np.maximum(0.0, 1.0 - v[inside_strip]**2))

    bright = disk & inside_strip & (u >= u_term)
    dark = disk & (~bright)

    return disk, bright, dark


def render_venus_frame(
    width: int,
    height: int,
    cx: float,
    cy: float,
    radius: float,
    phase_angle_deg: float,
    sun_angle_img_deg: float,
    bright_level: float = 230.0,
    blur_sigma: float = 1.2,
    background_level: float = 0.0
):
    """
    Render Venus bright silhouette on black background.
    """
    img = np.full((height, width), background_level, dtype=np.float32)

    _, bright, _ = analytical_masks(
        width, height, cx, cy, radius, phase_angle_deg, sun_angle_img_deg
    )

    img[bright] = bright_level

    if blur_sigma > 0:
        ksize = int(max(3, 2 * round(3 * blur_sigma) + 1))
        img = cv2.GaussianBlur(img, (ksize, ksize), blur_sigma)

    return img


def gaussian_spot(height, width, x0, y0, amp, sigma):
    yy, xx = np.indices((height, width), dtype=np.float32)
    rr2 = (xx - x0)**2 + (yy - y0)**2
    spot = amp * np.exp(-0.5 * rr2 / (sigma * sigma))
    return spot


def temporal_profile(n, kind="triangular"):
    """
    Generate a normalized temporal brightness profile for a fireball
    lasting n frames.

    Supported profiles:
        triangular : brightness increases and then decreases
        gaussian   : Gaussian brightness evolution
        flat       : constant brightness

    The maximum brightness is always normalized to 1.0.
    Works correctly for any duration >= 1 frame.
    """

    if n < 1:
        raise ValueError("Fireball duration must be >= 1 frame.")

    if n == 1:
        return np.array([1.0], dtype=np.float32)

    t = np.arange(n, dtype=np.float32)
    center = 0.5 * (n - 1)

    if kind == "triangular":
        # Avoid zero brightness for two-frame events
        prof = 1.0 - np.abs(t - center) / (center + 1.0)

    elif kind == "gaussian":
        sigma = max(0.35 * n, 0.8)
        prof = np.exp(-0.5 * ((t - center) / sigma)**2)

    elif kind == "flat":
        prof = np.ones(n, dtype=np.float32)

    else:
        raise ValueError(f"Unknown temporal profile: {kind}")

    # Normalize safely
    peak = float(np.max(prof))

    if not np.isfinite(peak) or peak <= 0:
        raise ValueError("Invalid fireball temporal profile.")

    prof /= peak

    return prof.astype(np.float32)


def choose_random_dark_position(dark_mask, rng, margin=3):
    ys, xs = np.where(dark_mask)
    if len(xs) == 0:
        raise RuntimeError("No dark-side pixels available for fireball placement.")
    idx = rng.integers(0, len(xs))
    return float(xs[idx]), float(ys[idx])


def translate_image(img, dx, dy):
    M = np.array([[1, 0, dx],
                  [0, 1, dy]], dtype=np.float32)
    out = cv2.warpAffine(
        img, M, (img.shape[1], img.shape[0]),
        flags=cv2.INTER_LINEAR,
        borderMode=cv2.BORDER_CONSTANT,
        borderValue=0
    )
    return out


# ============================================================
# MAIN SEQUENCE GENERATION
# ============================================================

def generate_sequence(args):
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    obs_dt = datetime.fromisoformat(args.date.replace("Z", "")).replace(tzinfo=timezone.utc)
    phase_angle_deg, illum_frac, delta_au = venus_phase_and_distance(obs_dt)

    rng = np.random.default_rng(args.seed)

    width = args.width
    height = args.height
    radius = args.diameter / 2.0

    cx0 = width / 2.0
    cy0 = height / 2.0

    # Build nominal dark-side mask at reference position
    _, _, dark_nominal = analytical_masks(
        width, height, cx0, cy0, radius, phase_angle_deg, args.sun_angle
    )

    # Parse fireball events
    events = []
    if args.events_file:
        # CSV format:
        # event_id,start_frame,duration,peak_amplitude,sigma_px,x_rel,y_rel,auto_place
        with open(args.events_file, "r", newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                events.append(FireballEvent(
                    event_id=int(row["event_id"]),
                    start_frame=int(row["start_frame"]),
                    duration=int(row["duration"]),
                    peak_amplitude=float(row["peak_amplitude"]),
                    sigma_px=float(row["sigma_px"]),
                    x_rel=None if row["x_rel"] == "" else float(row["x_rel"]),
                    y_rel=None if row["y_rel"] == "" else float(row["y_rel"]),
                    auto_place=str(row.get("auto_place", "True")).lower() in ("1", "true", "yes", "y")
                ))
    else:
        # Automatic random events
        for i in range(args.n_events):
            duration = int(rng.integers(args.min_duration, args.max_duration + 1))
            start_frame = int(rng.integers(0, max(1, args.n_frames - duration)))
            amp = float(rng.uniform(args.min_amp, args.max_amp))
            sigma = float(rng.uniform(args.min_sigma, args.max_sigma))
            events.append(FireballEvent(
                event_id=i + 1,
                start_frame=start_frame,
                duration=duration,
                peak_amplitude=amp,
                sigma_px=sigma,
                auto_place=True
            ))

    # Assign event positions on the dark side
    placed_events = []
    for ev in events:
        if ev.auto_place or ev.x_rel is None or ev.y_rel is None:
            x_abs, y_abs = choose_random_dark_position(dark_nominal, rng)
            x_rel = x_abs - cx0
            y_rel = y_abs - cy0
        else:
            x_rel = ev.x_rel
            y_rel = ev.y_rel

        placed_events.append(dict(
            event_id=ev.event_id,
            start_frame=ev.start_frame,
            duration=ev.duration,
            peak_amplitude=ev.peak_amplitude,
            sigma_px=ev.sigma_px,
            x_rel=x_rel,
            y_rel=y_rel
        ))

    # Save truth table
    truth_csv = outdir / "truth_events.csv"
    with open(truth_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "event_id", "start_frame", "duration",
            "peak_amplitude", "sigma_px", "x_rel_px", "y_rel_px"
        ])
        for ev in placed_events:
            writer.writerow([
                ev["event_id"], ev["start_frame"], ev["duration"],
                ev["peak_amplitude"], ev["sigma_px"],
                f"{ev['x_rel']:.3f}", f"{ev['y_rel']:.3f}"
            ])

    # Info file
    with open(outdir / "sequence_info.txt", "w", encoding="utf-8") as f:
        f.write("Synthetic Venus sequence\n")
        f.write(f"UTC date: {obs_dt.isoformat()}\n")
        f.write(f"Phase angle Sun-Venus-Earth (deg): {phase_angle_deg:.4f}\n")
        f.write(f"Illuminated fraction: {illum_frac:.6f}\n")
        f.write(f"Earth-Venus distance (au): {delta_au:.6f}\n")
        f.write(f"Image size: {width} x {height} px\n")
        f.write(f"Venus diameter: {args.diameter:.2f} px\n")
        f.write(f"Sun angle in image: {args.sun_angle:.2f} deg\n")
        f.write(f"Frames: {args.n_frames}\n")
        f.write(f"Noise sigma: {args.noise_sigma:.3f}\n")
        f.write(f"Jitter sigma (px): {args.jitter_sigma:.3f}\n")

    # Optional reference masks
    disk_ref, bright_ref, dark_ref = analytical_masks(
        width, height, cx0, cy0, radius, phase_angle_deg, args.sun_angle
    )

    preview = np.zeros((height, width, 3), dtype=np.uint8)
    preview[bright_ref] = (220, 220, 220)
    preview[dark_ref] = (0, 140, 0)
    cv2.imwrite(str(outdir / "reference_masks.png"), preview)

    # Precompute temporal profiles
    event_profiles = {}
    for ev in placed_events:
        event_profiles[ev["event_id"]] = temporal_profile(ev["duration"], kind=args.temporal_profile)

    # Generate frames
    for i in range(args.n_frames):
        # Random frame-to-frame shift
        dx = rng.normal(0.0, args.jitter_sigma) if args.jitter_sigma > 0 else 0.0
        dy = rng.normal(0.0, args.jitter_sigma) if args.jitter_sigma > 0 else 0.0

        cx = cx0 + dx
        cy = cy0 + dy

        frame = render_venus_frame(
            width=width,
            height=height,
            cx=cx,
            cy=cy,
            radius=radius,
            phase_angle_deg=phase_angle_deg,
            sun_angle_img_deg=args.sun_angle,
            bright_level=args.bright_level,
            blur_sigma=args.blur_sigma,
            background_level=0.0
        )

        # Add fireballs
        for ev in placed_events:
            start = ev["start_frame"]
            stop = start + ev["duration"]
            if start <= i < stop:
                k = i - start
                prof = event_profiles[ev["event_id"]][k]

                x_ev = cx + ev["x_rel"]
                y_ev = cy + ev["y_rel"]

                spot = gaussian_spot(
                    height, width,
                    x_ev, y_ev,
                    amp=ev["peak_amplitude"] * prof,
                    sigma=ev["sigma_px"]
                )
                frame += spot

        # Add noise
        if args.noise_sigma > 0:
            frame += rng.normal(0.0, args.noise_sigma, size=frame.shape).astype(np.float32)

        # Clip and save
        frame = np.clip(frame, 0, 255).astype(np.uint8)

        # If you want extra optical smear after adding fireball + noise
        if args.final_blur_sigma > 0:
            ksize = int(max(3, 2 * round(3 * args.final_blur_sigma) + 1))
            frame = cv2.GaussianBlur(frame, (ksize, ksize), args.final_blur_sigma)

        outname = outdir / f"venus_{i:04d}.bmp"
        cv2.imwrite(str(outname), frame)

    print("Done.")
    print(f"Output folder: {outdir}")
    print(f"Phase angle (deg): {phase_angle_deg:.3f}")
    print(f"Illuminated fraction: {illum_frac:.5f}")
    print(f"Truth table: {truth_csv}")


# ============================================================
# CLI
# ============================================================

def build_argparser():
    p = argparse.ArgumentParser(description="Generate synthetic BMP sequence of Venus with fireballs.")
    p.add_argument("--output", required=True, help="Output folder")
    p.add_argument("--date", required=True,
                   help="UTC date/time, e.g. 2026-07-05T17:34:00")
    p.add_argument("--width", type=int, default=512, help="Image width in pixels")
    p.add_argument("--height", type=int, default=512, help="Image height in pixels")
    p.add_argument("--diameter", type=float, required=True,
                   help="Venus disk diameter in pixels")
    p.add_argument("--sun-angle", type=float, default=0.0,
                   help="Direction toward the Sun in image coordinates (deg). "
                        "0=right, 90=down, 180=left, 270=up")
    p.add_argument("--n-frames", type=int, default=200, help="Number of frames")
    p.add_argument("--noise-sigma", type=float, default=4.0, help="Gaussian noise sigma")
    p.add_argument("--jitter-sigma", type=float, default=0.8,
                   help="Frame-to-frame Venus centroid jitter sigma in pixels")
    p.add_argument("--bright-level", type=float, default=230.0,
                   help="Brightness level of illuminated Venus")
    p.add_argument("--blur-sigma", type=float, default=1.2,
                   help="Blur sigma applied to Venus silhouette")
    p.add_argument("--final-blur-sigma", type=float, default=0.0,
                   help="Optional extra blur applied to final frame")
    p.add_argument("--seed", type=int, default=42, help="Random seed")

    # Fireball generation
    p.add_argument("--n-events", type=int, default=5,
                   help="Number of random events if --events-file is not provided")
    p.add_argument("--min-duration", type=int, default=1, help="Minimum event duration in frames")
    p.add_argument("--max-duration", type=int, default=4, help="Maximum event duration in frames")
    p.add_argument("--min-amp", type=float, default=40.0, help="Minimum peak amplitude")
    p.add_argument("--max-amp", type=float, default=120.0, help="Maximum peak amplitude")
    p.add_argument("--min-sigma", type=float, default=1.0, help="Minimum fireball sigma in pixels")
    p.add_argument("--max-sigma", type=float, default=2.5, help="Maximum fireball sigma in pixels")
    p.add_argument("--temporal-profile", choices=["triangular", "gaussian", "flat"],
                   default="triangular", help="Temporal profile of fireball brightness")
    p.add_argument("--events-file", default="",
                   help="Optional CSV file specifying custom events")

    return p


def main():
    parser = build_argparser()
    args = parser.parse_args()
    generate_sequence(args)


if __name__ == "__main__":
    main()