import matplotlib.pyplot as plt

from kbo_occultation import (
    compute_lightcurve,
    KBOConfig,
    StarConfig,
    BandpassConfig,
    GridConfig,
    NumericalConfig,
)

kbo = KBOConfig(radius_km=0.5, distance_au=40.0)
star = StarConfig(temperature_K=20000, angular_radius_mas=0.03)
band = BandpassConfig(400, 430, 25)
grid = GridConfig(8.0, 1000)
num = NumericalConfig()

x, I = compute_lightcurve(kbo, star, band, grid, num)

plt.plot(x, I)
plt.xlabel("Position (km)")
plt.ylabel("Normalized intensity")
plt.show()
