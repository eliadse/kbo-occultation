import matplotlib.pyplot as plt

from kbo_occultation import (
    compute_lightcurve,
    KBOConfig,
    StarConfig,
    BandpassConfig,
    GridConfig,
    NumericalConfig,
    simulate_poly_point
)

kbo = KBOConfig(radius_m=500, distance_au=40.0)
star = StarConfig(temperature_K=20000, angular_radius_mas=0.03)
band = BandpassConfig(300, 500, 25)
grid = GridConfig(5000, 500)
num = NumericalConfig()

x, Ipoly = simulate_poly_point(kbo, star, band, grid, num)
x, Istar = compute_lightcurve(kbo, star, band, grid, num)

band = BandpassConfig(400, 430, 5)
x, Imono = simulate_poly_point(kbo, star, band, grid, num)

plt.plot(x, Imono, label='Monochromatic point source')
plt.plot(x, Ipoly, label='Polychromatic point source')
plt.plot(x, Istar, label='0.03mas radius star')
plt.grid()
plt.legend(loc="upper left")
plt.xlabel("Position (m)")
plt.ylabel("Normalized intensity")
plt.show()
