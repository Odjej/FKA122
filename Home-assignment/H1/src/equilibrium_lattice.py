from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

data_file = Path(__file__).with_name("lattice_potential_energy.csv")
lattice_constant, potential_energy = np.genfromtxt(
    data_file, delimiter=",", usecols=(0, 1), unpack=True
)

p = np.polyfit(lattice_constant, potential_energy, 4,)
x = np.linspace(3,6, 1000)

plt.scatter(lattice_constant, potential_energy)
plt.plot(x, p[4] + p[3] * x + p[2] * x**2 + p[1] * x**3 + p[0] * x**4 )
plt.show()