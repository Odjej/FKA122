from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

DATA_DIR = Path(__file__).resolve().parent

positions = np.loadtxt(DATA_DIR / "positions.csv", delimiter=",")
potential = np.loadtxt(DATA_DIR / "potential_energy.csv")
kinetic = np.loadtxt(DATA_DIR / "kinetic_energy.csv")
frequencies = np.loadtxt(DATA_DIR / "fft_frqs.csv")
spectrum = np.loadtxt(DATA_DIR / "powerspectrum.csv")

dt = 2 * np.pi / (len(frequencies) * (frequencies[1] - frequencies[0]))
time = np.arange(len(positions)) * dt

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

for i in range(3):
    ax1.plot(time, positions[:, i], label=f"Particle {i+1}")

ax1.set_ylabel("Position (Å)")
ax1.legend()

ax2.plot(time, potential, label="Potential")
ax2.plot(time, kinetic, label="Kinetic")
ax2.plot(time, potential + kinetic, label="Total")

ax2.set_xlabel("Time (ps)")
ax2.set_ylabel("Energy (eV)")
ax2.legend()
ax2.set_xlim(0, 0.25)

plt.tight_layout()
plt.savefig(DATA_DIR / "simulation.png", dpi=150)

plt.figure(figsize=(9, 5))

plt.plot(frequencies / (2 * np.pi), spectrum)

plt.xlabel("Frequency (ps$^{-1}$)")
plt.ylabel("Power spectrum (eV^2 * ps)")
plt.yscale("log")
#plt.xlim(left=0)

plt.tight_layout()
plt.savefig(DATA_DIR / "spectrum.png", dpi=150)

plt.show()