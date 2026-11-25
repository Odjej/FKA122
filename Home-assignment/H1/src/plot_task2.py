
import os
import re
import glob

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------- paths and settings

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
fig_dir = os.path.join(base_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)

n_atoms = 256
k_B = 8.617333e-5       # eV/K
t_transient = 1000.0    # fs, discarded when averaging T
t_total = 5000.0        # fs, run length

# Consistent plotting settings
plt.rcParams.update({
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "legend.fontsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 600,
    "lines.linewidth": 1.0,
})


# -------------------------------------------------------------- load data

runs = []

for path in glob.glob(os.path.join(base_dir, "data", "energies_dt_*fs.csv")):

    match = re.search(r"energies_dt_([\d.]+)fs\.csv", path)

    if match is None:
        continue

    dt = float(match.group(1))

    t, e_pot, e_kin, e_tot = np.genfromtxt(
        path, delimiter=",", skip_header=1, unpack=True
    )

    ok = np.isfinite(e_tot)

    t = t[ok]
    e_pot = e_pot[ok]
    e_kin = e_kin[ok]
    e_tot = e_tot[ok]

    if len(t) == 0:
        continue

    runs.append({
        "dt": dt,
        "t": t,
        "e_kin": e_kin,
        "e_tot": e_tot,
        "completed": t[-1] >= 0.99 * t_total
    })

runs.sort(key=lambda r: r["dt"])
by_dt = {r["dt"]: r for r in runs}


# ---------------------------------------------------------- statistics

for r in runs:

    r["dE"] = r["e_tot"] - r["e_tot"][0]

    r["std"] = np.std(r["e_tot"])

    # Linear regression to estimate systematic energy drift
    slope = np.polyfit(r["t"], r["e_tot"], 1)[0]

    r["drift"] = abs(slope) * t_total

    # Temperature
    r["T"] = 2 * r["e_kin"] / (3 * n_atoms * k_B)

    valid_T = r["t"] >= t_transient

    r["T_mean"] = (
        np.mean(r["T"][valid_T]) if np.any(valid_T) else np.nan
    )

    # Flag unstable simulations
    r["exploded"] = (not r["completed"]) or r["std"] > 5.0


# ---------------------------------------------------------- summary table

print(
    f"{'dt (fs)':>8} {'std(Etot) eV':>13} "
    f"{'drift (eV)':>12} {'<T> K':>8}  note"
)

for r in runs:

    note = "exploded" if r["exploded"] else ""

    T_str = "-" if r["exploded"] else f"{r['T_mean']:.1f}"

    print(
        f"{r['dt']:8.1f} {r['std']:13.3e} "
        f"{r['drift']:12.3e} {T_str:>8}  {note}"
    )


# ==========================================================
# FIGURE 1-3: ENERGY CONSERVATION REGIMES
# ==========================================================

regimes = [
    (
        "small_enough",
        "Small enough",
        [0.5, 1.0, 2.0, 4.0],
        (-0.5, 0.10)
    ),
    (
        "a_bit_too_large",
        "A bit too large",
        [6.0, 8.0, 12.0, 16.0, 18.0],
        (-8, 2)
    ),
    (
        "far_too_large",
        "Far too large",
        [20.0, 24.0],
        None
    ),
]

for filename, title, dts_here, ylim in regimes:

    fig, ax = plt.subplots(figsize=(7, 4.5))

    cols = plt.cm.viridis(
        np.linspace(0, 0.9, len(dts_here))
    )

    for dt, c in zip(dts_here, cols):

        if dt not in by_dt:
            continue

        r = by_dt[dt]

        ax.plot(
            r["t"] / 1000,
            r["dE"],
            color=c,
            label=rf"$\Delta t = {dt:g}$ fs"
        )

    if ylim is None:
        ax.set_yscale("symlog", linthresh=1.0)
    else:
        ax.set_ylim(ylim)

    ax.set_xlabel("Time (ps)")
    ax.set_ylabel(
        r"$E_{\mathrm{tot}}(t)-E_{\mathrm{tot}}(0)$ (eV)"
    )

    ax.legend(ncol=2, loc="best")
    ax.grid(alpha=0.2)

    fig.tight_layout()

    fig.savefig(
        os.path.join(fig_dir, f"energy_{filename}.png"),
        bbox_inches="tight",
        dpi=600
    )

    plt.close(fig)


# ==========================================================
# FIGURE 4: SCALING OF std(E_tot) WITH TIMESTEP
# ==========================================================

good = [r for r in runs if not r["exploded"]]

dts = np.array([r["dt"] for r in good])
stds = np.array([r["std"] for r in good])

dt_fit_max = 8.0

m = dts <= dt_fit_max

slope, icpt = np.polyfit(
    np.log(dts[m]),
    np.log(stds[m]),
    1
)

print(
    f"\nstd(E_tot) ~ dt^{slope:.2f} "
    f"(fit for dt <= {dt_fit_max:g} fs)"
)

fig, ax = plt.subplots(figsize=(7, 4.5))

ax.loglog(
    dts[m],
    stds[m],
    "o",
    color="C0",
    label="Used in fit"
)

ax.loglog(
    dts[~m],
    stds[~m],
    "o",
    color="C1",
    mfc="none",
    label="Not used in fit"
)

xx = np.logspace(
    np.log10(dts.min()),
    np.log10(dts.max()),
    100
)

ax.loglog(
    xx,
    np.exp(icpt) * xx**slope,
    "r--",
    label=rf"Fit: $\propto \Delta t^{{{slope:.2f}}}$"
)

ax.loglog(
    xx,
    stds[0] * (xx / dts[0])**2,
    "k:",
    label=r"$\propto \Delta t^2$ (reference)"
)

ax.set_xlabel(r"Timestep $\Delta t$ (fs)")
ax.set_ylabel(r"$\mathrm{std}(E_{\mathrm{tot}})$ (eV)")

ax.legend()
ax.grid(alpha=0.2, which="both")

fig.tight_layout()

fig.savefig(
    os.path.join(fig_dir, "dt_scaling.png"),
    bbox_inches="tight",
    dpi=600
)

plt.close(fig)


# ==========================================================
# FIGURE 5: MEAN TEMPERATURE VERSUS TIMESTEP
# ==========================================================

fig, ax = plt.subplots(figsize=(7, 4.5))

ax.plot(
    dts,
    [r["T_mean"] for r in good],
    "o-",
    label="Mean temperature"
)

ax.axhspan(
    600,
    850,
    color="gray",
    alpha=0.15,
    label="Expected range"
)

ax.set_xlabel(r"Timestep $\Delta t$ (fs)")
ax.set_ylabel(r"Mean temperature $\langle T \rangle$ (K)")

ax.legend()
ax.grid(alpha=0.2)

fig.tight_layout()

fig.savefig(
    os.path.join(fig_dir, "temperature_vs_dt.png"),
    bbox_inches="tight",
    dpi=600
)

plt.close(fig)


print(f"\nFigures saved in: {fig_dir}")
