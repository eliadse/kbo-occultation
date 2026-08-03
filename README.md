# kbo-occultation

Simulation and detection of stellar occultation diffraction patterns by sub-km
Kuiper Belt Objects, for fast photometry with the MAGIC + LST-1 array.
Authors: E. do Souto Espiñeira, T. Hassan, V. Pascual

## Features

- Fresnel diffraction (Lommel–Hankel formulation)
- Polychromatic integration (Planck + filter response)
- Configurable astrophysical parameters
- Blind matched-filter search with χ² shape veto
- Signal injection / recovery Monte Carlo (single telescope and full array)
- 3-telescope coincidence matching (MAGIC-1 / MAGIC-2 / LST-1)
- Upper limits on the cumulative KBO surface density

## Installation

```bash
pip install -e .
```

## Example

```python
from kbo_occultation import *

kbo = KBOConfig(radius_m=500.0, distance_au=40.0)
star = StarConfig(temperature_K=20000, angular_radius_mas=0.03)
band = BandpassConfig(400, 430, 25)
grid = GridConfig(8.0, 1000)
num = NumericalConfig()

x, I = compute_lightcurve(kbo, star, band, grid, num)
```

## Version

```toml
version = "0.2.0"
```
