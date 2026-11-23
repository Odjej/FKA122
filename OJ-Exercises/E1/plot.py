import os
import numpy as np
import matplotlib.pyplot as plt

data_dir = "data"
fig_dir = "figures"
os.makedirs(fig_dir, exist_ok=True)


def load(name):
    """Load a CSV from data_dir, or return None (with a message) if missing."""
    path = os.path.join(data_dir, name)
    if not os.path.exists(path):
        print(f"Skipping: {path} not found (run the matching ./program task first)")
        return None
    return np.genfromtxt(path, delimiter=",", names=True)


# ------------------------------------------------------------------
# Trajectory and energies (./program trajectory)
# ------------------------------------------------------------------
data = load("verlet_output.csv")
if data is not None:
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
# Energy conservation vs time step (./program scan)
# ------------------------------------------------------------------
scan = load("energy_vs_dt.csv")
if scan is not None:
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
# Power spectrum (./program spectrum)
# ------------------------------------------------------------------
spec = load("spectrum.csv")
if spec is not None:
    f_thz = spec["freq_THz"]
    power = spec["power"]

    fig4, ax4 = plt.subplots(figsize=(6, 4))
    ax4.semilogy(f_thz, power)
    ax4.set_xlabel("Frequency (THz)")
    ax4.set_ylabel(r"Power spectrum of $q_1$ (arb. units)")
    ax4.set_xlim(-5, 120)
    fig4.tight_layout()
    fig4.savefig(os.path.join(fig_dir, "fig5_spectrum.png"), dpi=300)

# ------------------------------------------------------------------
# Harmonic FPUT chain, normal-mode energies (./program modes)
# ------------------------------------------------------------------
modes = load("normal_mode_energies.csv")
if modes is not None:
    tm = modes["time"]

    # E_k(t), k = 1..5 on a linear axis (the plot the assignment asks for)
    fig5, ax5 = plt.subplots(figsize=(6, 4))
    for k in range(1, 6):
        ax5.plot(tm, modes[f"E{k}"], label=rf"$E_{k}$")
    ax5.set_xlabel("Time")
    ax5.set_ylabel(r"$E_k$")
    ax5.set_xlim(0, tm.max())
    ax5.legend()
    fig5.tight_layout()
    fig5.savefig(os.path.join(fig_dir, "fig6_harmonic_mode_energies.png"), dpi=300)

    # E_2..E_5 on a log axis, since E_1 hides them on the linear plot
    fig6, ax6 = plt.subplots(figsize=(6, 4))
    for k in range(2, 6):
        ax6.semilogy(tm, np.maximum(modes[f"E{k}"], 1e-30), label=rf"$E_{k}$")
    ax6.set_xlabel("Time")
    ax6.set_ylabel(r"$E_k$")
    ax6.set_xlim(0, tm.max())
    ax6.legend()
    fig6.tight_layout()
    fig6.savefig(os.path.join(fig_dir, "fig7_harmonic_modes_log.png"), dpi=300)

    # Total energy: direct (kinetic + potential) vs sum of mode energies
    fig7, ax7 = plt.subplots(figsize=(6, 4))
    e0 = modes["Etot_direct"][0]
    ax7.plot(tm, (modes["Etot_direct"] - e0) / e0, ls="-", label="direct (kin + pot)")
    ax7.plot(tm, (modes["Etot_modes"] - e0) / e0, ls="--", label=r"$\sum_k E_k$")
    ax7.set_xlabel("Time")
    ax7.set_ylabel(r"$(E - E_0)/E_0$")
    ax7.set_xlim(0, tm.max())
    ax7.legend()
    fig7.tight_layout()
    fig7.savefig(os.path.join(fig_dir, "fig8_harmonic_total_energy.png"), dpi=300)

plt.show()