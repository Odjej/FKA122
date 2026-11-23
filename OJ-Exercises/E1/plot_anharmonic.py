import os
import numpy as np
import matplotlib.pyplot as plt

data_dir, fig_dir = "data", "figures"
os.makedirs(fig_dir, exist_ok=True)
N = 32

def load(name):
    return np.genfromtxt(os.path.join(data_dir, name), delimiter=",", names=True)

# Task 9: E_k(t), k = 1..5, for each alpha
for alpha in ["0.01", "0.1"]:
    d = load(f"modes_alpha{alpha}.csv")
    fig, ax = plt.subplots(figsize=(6, 4))
    for k in range(1, 6):
        ax.plot(d["time"], d[f"E{k}"], label=rf"$E_{k}$")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$E_k$")
    ax.set_title(rf"$\alpha = {alpha}$")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, f"modes_alpha{alpha}.png"), dpi=300)

# Task 10: time-averaged energies of all modes, log-log
for alpha in ["0.01", "0.1"]:
    d = load(f"long_alpha{alpha}.csv")
    t = d["time"][1:]                      # skip t = 0 on a log axis
    fig, ax = plt.subplots(figsize=(6, 4))
    colors = plt.cm.viridis(np.linspace(0, 1, N))
    for k in range(1, N + 1):
        ax.loglog(t, d[f"A{k}"][1:], color=colors[k - 1], lw=0.8,
                  label=rf"$k={k}$" if k in (1, 2, 3, N) else None)
    ax.axhline(1.0, color="k", ls="--", lw=0.8, label=r"$E_0/N$")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"$\langle E_k\rangle_{\mathrm{time}}$")
    ax.set_title(rf"$\alpha = {alpha}$")
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, f"time_avg_alpha{alpha}.png"), dpi=300)

plt.show()

'''
0.1
The time averages A₁…A₃₂ all lie between about 0.76 and 1.29, close to E₀/N = 32/32 = 1. 
They don't all sit exactly at 1, and the low modes are slightly above and the middle and high modes slightly below. 
That is the expected finite-time scatter.
The instantaneous Eₖ still fluctuate strongly (0.03 to 3.9). 
That is normal, because equipartition applies to the time average and not to a single snapshot.
Etot = 31.97 against 32, a deviation of about 0.1%.

0.01
The time averages A₁…A₃₂ are very far from equipartition.
The averages fall off roughly exponentially with mode number: A₁ ≈ 11.9, A₂ ≈ 8.8, A₃ ≈ 6.0, A₄ ≈ 3.3, A₅ ≈ 1.3, A₆ ≈ 0.42, 
and so on down to about 10⁻¹³ for the highest modes.
Almost all the energy (the A values sum to about 32) sits in the lowest ~5 modes. 
Modes above ~10 hold almost nothing, so the high modes are many orders of magnitude from the equipartition value of 1.
Energy has left mode 1 (A₁ ≈ 12 compared with the initial 32), so there is exchange.
But it is confined to the lowest modes, which is the original FPUT result. 
The system sits in a long-lived, non-equipartitioned state.
Etot = 31.9985, so the energy is conserved very well
'''
