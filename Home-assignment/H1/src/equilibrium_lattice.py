import os

import numpy as np
import matplotlib.pyplot as plt

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
fig_dir = os.path.join(base_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)
data_file = os.path.join(base_dir, "data", "lattice_potential_energy.csv")
a, E = np.genfromtxt(data_file, delimiter=",", usecols=(0, 1), unpack=True)

n_cells = 64
mask = (a >= 60**(1/3)) & (a <= 70**(1/3))
a_fit, E_fit = a[mask], E[mask]

p = np.polyfit(a_fit, E_fit, 2)

a0_quad = -p[1] / (2 * p[0]) # analytical minimum of quadratic fit
xs = np.linspace(a_fit.min(), a_fit.max(), 1000)
print(f"quadratic a0 = {a0_quad:.8f} Å")
plt.figure(figsize=(10, 6))
plt.scatter(a_fit**3, E_fit / n_cells, label=r"Calculated energies")
plt.plot(xs**3, np.polyval(p, xs) / n_cells, "r", label=r"Quadratic fit")
plt.axvline(a0_quad**3, color="r", ls=":", lw=0.8, label=rf"Fit minimum: $a_0 = {a0_quad:.6f}$ Å")
plt.xlabel(r'Unit cell volume (Å$^3$)')
plt.ylabel(r'Energy (eV/unit cell)')

plt.xlim((63,68))
plt.ylim((-10.50, -10.46))
plt.legend()


plt.savefig(os.path.join(fig_dir, "equilibrium_lattice_constant.png"), dpi=800)
plt.show()

#[4.03230237]