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
star = StarConfig(temperature_K=30000, angular_radius_mas=0.03)
kbo = KBOConfig(radius_m=50, distance_au=40.0)

# Simulate mono point
band = BandpassConfig(400, 430, 5)
x, Imono = simulate_poly_point(kbo, star, band, grid, num)

# Simulate poly point
band = BandpassConfig(100, 800, 25)
#x, Ipoly = simulate_poly_point(kbo, star, band, grid, num)

# Simulate poly star, no offset
x, Istar = compute_lightcurve(kbo, star, band, grid, num)
x, Istar_SII = compute_lightcurve(kbo, star, band, grid, num, True)

#kbo = KBOConfig(radius_m=500, distance_au=40.0, impact_parameter_m = 100)
#x, Ioff = compute_lightcurve(kbo, star, band, grid, num)

plt.plot(x, Imono, label='Monochromatic point source')
#plt.plot(x, Ipoly, label='Polychromatic point source')
plt.plot(x, Istar, label=f'{star.angular_radius_mas}mas radius star')
plt.plot(x, Istar_SII, label=f'{star.angular_radius_mas}mas radius star SII')
#plt.plot(x, Ioff, label=f'{star.angular_radius_mas}mas radius star, {kbo.impact_parameter_m}m offset')
plt.grid()
plt.legend(loc="upper left")
plt.xlabel("Position (m)")
plt.ylabel("Normalized intensity")
plt.savefig(f"./{kbo.radius_m}mKBO_{kbo.impact_parameter_m}mOffset_{star.angular_radius_mas}mas.png")
plt.show()
