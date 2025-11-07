"""Plot Mikheyev–Smirnov–Wolfenstein (MSW) resonance curves.

This script computes effective matter-modified neutrino mass-squared
eigenvalues as a function of electron density and produces a publication
quality plot (`ThesisPlots/MSW/MSW_Resonance_<E>MeV.pdf`). The code uses
standard units (eV, eV^2) and simple illustrative parameters for plotting.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

Gf = 1.166e-23  # Fermi constant in 1/eV²
m1 = 0.0  # mass eigenstate in eV²
m2 = 7.5e-5  # mass eigenstate in eV² (corrected value)
dm2 = m2 - m1  # mass squared difference in eV²
E = 1e6  # neutrino energy in eV
theta = 33.68 / 180 * np.pi  # mixing angle in radians
NA = 6.022e23 * (1.97327e-5) ** 3  # Avogadro's number in ev^3

B = 2 * np.sqrt(2) * E * Gf  # matter potential in eV²


def lambda_m(n_e, sign):
    Acc = 2 * np.sqrt(2) * E * Gf * n_e * NA
    Deltam2_M = np.sqrt(
        (dm2 * np.cos(2 * theta) - Acc) ** 2 + (dm2 * np.sin(2 * theta)) ** 2
    )
    return 0.5 * (m1 + m2 + Acc + sign * Deltam2_M)


def extrapolate_p(n_e):
    return 0.5 * (m2 + 2 * (B * n_e * NA) - dm2 * np.cos(2 * theta))


n_e_Reso = dm2 * np.cos(2 * theta) / (2 * np.sqrt(2) * E * Gf * NA)

density_max = 1.1e3
ne_vals = np.linspace(0, density_max, 1000)  # electron density in g/cm⁻³
lambda_plus = lambda_m(ne_vals, 1)
lambda_minus = lambda_m(ne_vals, -1)
extrapolated_p_vals = extrapolate_p(ne_vals)
extrapolated_m_vals = 0.5 * (m2 + dm2 * np.cos(2 * theta)) * np.ones_like(ne_vals)

plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "font.size": 30,
    }
)

plt.figure(figsize=(12, 9))

plt.plot(ne_vals, lambda_plus, label=r"$m^2_M(\nu_2)$", color="red", linewidth=4)
plt.plot(ne_vals, lambda_minus, label=r"$m^2_M(\nu_1)$", color="blue", linewidth=4)

plt.plot(
    ne_vals,
    extrapolated_p_vals,
    "k",
    dashes=(10, 5),
    linewidth=2,
    label=r"$\lim\limits_{N \to \infty} m^2_M(\nu_{1,2})$",
)
plt.plot(ne_vals, extrapolated_m_vals, "k", dashes=(10, 5), linewidth=2)

plt.axvline(x=n_e_Reso, color="gray", linestyle="dotted", linewidth=3, zorder=0)
plt.text(
    n_e_Reso + plt.xlim()[1] * 0.01,
    plt.ylim()[1] * 0.7,
    r"$N^{\mathrm{resonance}}_\odot$",
    fontsize=35,
    color="gray",
    ha="left",
    va="top",
)

plt.axvline(x=160, color="gray", linestyle="dotted", linewidth=3, zorder=0)
plt.text(
    160 - plt.xlim()[1] * 0.11,
    plt.ylim()[1] * 0.56,
    r"$N^{\mathrm{core}}_\odot$",
    fontsize=35,
    color="gray",
    ha="left",
    va="top",
)

ax = plt.gca()
formatter = ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-2, 2))  # force scientific notation
ax.yaxis.set_major_formatter(formatter)
ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))

plt.xlabel(r"Density [g/cm$^3$]", fontsize=35)
plt.ylabel(r"$m^2_M$ [eV$^2$]", fontsize=35, labelpad=20)
plt.grid(True)
plt.legend(handlelength=1.2, loc="upper center", bbox_to_anchor=(0.8, 0.66))
plt.xlim(0, density_max)
# plt.ylim(0, 3.5e-5)
plt.savefig(f"ThesisPlots/MSW/MSW_Resonance_{E/1e6:.1f}MeV.pdf", bbox_inches="tight")
