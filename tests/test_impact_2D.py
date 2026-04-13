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

grid = GridConfig(5000, 500)
num = NumericalConfig()
star = StarConfig(temperature_K=20000, angular_radius_mas=0.03)
kbo = KBOConfig(radius_m=500, distance_au=40.0)

# Simulate mono point
band = BandpassConfig(400, 430, 5)
x, Imono = simulate_poly_point(kbo, star, band, grid, num)

# Simulate poly point
band = BandpassConfig(300, 500, 25)
#x, Ipoly = simulate_poly_point(kbo, star, band, grid, num)

# Simulate poly star, no offset
x, Istar = compute_lightcurve(kbo, star, band, grid, num)

kbo = KBOConfig(radius_m=500, distance_au=40.0, impact_parameter_m = 100)
x, Ioff = compute_lightcurve(kbo, star, band, grid, num)

plt.plot(x, Imono, label='Monochromatic point source')
#plt.plot(x, Ipoly, label='Polychromatic point source')
plt.plot(x, Istar, label='0.03mas radius star')
plt.plot(x, Ioff, label='Monochromatic offset source')
plt.grid()
plt.legend(loc="upper left")
plt.xlabel("Position (m)")
plt.ylabel("Normalized intensity")
plt.show()