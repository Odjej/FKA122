from pathlib import Path

from scipy.optimize import minimize

import numpy as np
import matplotlib.pyplot as plt
data_file = Path(__file__).with_name("lattice_potential_energy.csv")
lattice_constant, potential_energy = np.genfromtxt(
    data_file, delimiter=",", usecols=(0, 1), unpack=True
)

n_atoms = 64


mask1 = lattice_constant >= 60**(1/3)

cropped_lattice_constant = lattice_constant[mask1]
cropped_potential_energy = potential_energy[mask1]

mask2 = cropped_lattice_constant <= 70**(1/3)

cropped_lattice_constant = cropped_lattice_constant[mask2]
cropped_potential_energy = cropped_potential_energy[mask2]

p = np.polyfit(cropped_lattice_constant, cropped_potential_energy, 4,)
x = np.linspace(3,6, 1000)

def f(x):
    return (p[4] + p[3] * x + p[2] * x**2 + p[1] * x**3 + p[0] * x**4)


plt.scatter(cropped_lattice_constant**3, cropped_potential_energy / n_atoms)
plt.plot(x**3, (p[4] + p[3] * x + p[2] * x**2 + p[1] * x**3 + p[0] * x**4) / n_atoms)

plt.xlabel('Unit cell volume (Å$^3$)')
plt.ylabel('Energy (eV/unit cell)')

plt.xlim((64,68))
plt.ylim((-10.50, -10.46))

fit = minimize(f, x0 = 4, bounds=((3.5,4.5),))
print(fit.x)
plt.show()

#[4.02977927]