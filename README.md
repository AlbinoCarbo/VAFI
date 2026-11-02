# VAFI - Venus Atmospheric Flash Investigator

VAFI is a Python program for detecting candidate luminous events
in the nightside atmosphere of Venus.

The software accepts BMP and FITS image sequences and AVI videos.
It uses an ALFI-inspired spatial and temporal detection algorithm,
automatic image registration, and an interactive Venus phase mask.

The repository also contains a synthetic Venus fireball generator
for injection-and-recovery tests.

## Installation

Python 3.10 or later is required.

```bash
py -m pip install -r requirements.txt
```

## Fireball detection

```bash
py VAFI_v7.py --help
```

## Synthetic fireball generator

```bash
py synthetic_venus_fireballs.py --help
```

## Documentation

See docs/VAFI_Observer_Quick_Manual_v7.pdf for detailed
instructions, parameters, examples, and output descriptions.

## Scientific use

VAFI identifies candidate flashes, not confirmed atmospheric
impacts. Candidate events require independent verification.
