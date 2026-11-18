import os
import numpy as np
import matplotlib.pyplot as plt

data_dir = "data"
fig_dir = "figures"
os.makedirs(fig_dir, exist_ok=True)

# ------------------------------------------------------------------
# Trajectory and energies
# ------------------------------------------------------------------
data = np.genfromtxt(os.path.join(data_dir, "verlet_output.csv"),
                     delimiter=",", names=True)

t = data["time"]                 # ps
q1, q2, q3 = data["q1"], data["q2"], data["q3"]   # Å
ekin = data["Ekin"] * 1000       # eV -> meV
epot = data["Epot"] * 1000
etot = data["Etot"] * 1000

t_max = 0.25                     # ps, plotting window
mask = t <= t_max

# displacements
fig1, ax1 = plt.subplots(figsize=(6, 4))
ax1.plot(t[mask], q1[mask], label=r"$q_1$ (O)")
ax1.plot(t[mask], q2[mask], label=r"$q_2$ (C)")
ax1.plot(t[mask], q3[mask], label=r"$q_3$ (O)")
ax1.set_xlabel("Time (ps)")
ax1.set_ylabel(r"Displacement ($\mathrm{\AA}$)")
ax1.set_xlim(0, t_max)
ax1.legend()
fig1.tight_layout()
fig1.savefig(os.path.join(fig_dir, "fig2_displacements.png"), dpi=300)

# energies
fig2, ax2 = plt.subplots(figsize=(6, 4))
ax2.plot(t[mask], etot[mask], label="Total energy")
ax2.plot(t[mask], epot[mask], label="Pot. energy")
ax2.plot(t[mask], ekin[mask], label="Kin. energy")
ax2.set_xlabel("Time (ps)")
ax2.set_ylabel("Energy (meV)")
ax2.set_xlim(0, t_max)
ax2.set_ylim(0, 8)
ax2.legend()
fig2.tight_layout()
fig2.savefig(os.path.join(fig_dir, "fig3_energies.png"), dpi=300)

# ------------------------------------------------------------------
# Energy conservation vs time step
# ------------------------------------------------------------------
scan = np.genfromtxt(os.path.join(data_dir, "energy_vs_dt.csv"),
                     delimiter=",", names=True)
dt_fs = scan["dt_ps"] * 1000     # ps -> fs
dev = scan["max_rel_dev"]

fig3, ax3 = plt.subplots(figsize=(6, 4))
ax3.loglog(dt_fs, dev, "o-", label="Simulation")
# reference line with slope 2, anchored at the largest time step
ax3.loglog(dt_fs, dev[0] * (dt_fs / dt_fs[0])**2, "--", color="gray",
           label=r"$\propto \Delta t^2$")
ax3.set_xlabel(r"$\Delta t$ (fs)")
ax3.set_ylabel(r"max $|E(t)-E(0)|/E(0)$")
ax3.legend()
fig3.tight_layout()
fig3.savefig(os.path.join(fig_dir, "fig4_energy_vs_dt.png"), dpi=300)

# ------------------------------------------------------------------
# Power spectrum
# ------------------------------------------------------------------
spec = np.genfromtxt(os.path.join(data_dir, "spectrum.csv"),
                     delimiter=",", names=True)
f_thz = spec["freq_THz"]
power = spec["power"]

fig4, ax4 = plt.subplots(figsize=(6, 4))
ax4.semilogy(f_thz, power)
ax4.set_xlabel("Frequency (THz)")
ax4.set_ylabel(r"Power spectrum of $q_1$ (arb. units)")
ax4.set_xlim(-5, 120)
fig4.tight_layout()
fig4.savefig(os.path.join(fig_dir, "fig5_spectrum.png"), dpi=300)

plt.show()