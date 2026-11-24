from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
data_file = Path(__file__).parent.parent / "energies_over_time.csv"
time, potential_energy, kinetic_energy, total_energy = np.genfromtxt(
    data_file,
    delimiter=",",
    skip_header=1,
    usecols=(0, 1, 2, 3),
    unpack=True,
)

boltzmann = scipy.constants.Boltzmann
number_of_atoms = 4**4
mask = time >= 0.5
average_kinetic_energy = np.mean(kinetic_energy[mask]) * scipy.constants.electron_volt
average_temperature = average_kinetic_energy * 2 / (3* boltzmann * number_of_atoms)
print(f"Average kinetic energy: {average_kinetic_energy:.12f}")
print(f"Average temperature: {average_temperature:.12f}")
plt.plot(time, kinetic_energy)
plt.show()
