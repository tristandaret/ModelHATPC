"""Interactive minimal oscillation probability explorer.

This script provides an interactive matplotlib visualization of neutrino
oscillation probabilities P(ν_alpha -> ν_beta) as a function of baseline,
for configurable neutrino energy and mixing parameters. It was used to
produce illustrative plots for thesis figures and teaching demonstrations.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Button

# VARIABLES -----------------------------------------------------------------------------
# Flavor to index mapping
flavors = ["e", "mu", "tau"]
iflavor = {"e": 0, "mu": 1, "tau": 2}
ilatex = {0: r"\nu_e", 1: r"\nu_\mu", 2: r"\nu_\tau"}
alpha = 1
#  Cases: reactor: 0, accelerator: 1
case = 1
decimal = [1, 0]
Emin = [0.1, 100]
Enom = [4, 600]
Emax = [10, 4e3]
Lnom = [300, 2e3]
Lmax = [1e3, 5e3]
# NuFIT 6.0 (Normal Ordering) best-fit values
theta_12_0 = np.deg2rad(33.68)
theta_23_0 = np.deg2rad(43.3)
theta_13_0 = np.deg2rad(8.56)
delta_cp_0 = np.deg2rad(122)
# Mass-squared differences in eV^2
dm12_0 = 7.49e-5
dm23_0 = 2.513e-3
dm13_0 = dm12_0 + dm23_0


# FUNCTIONS -----------------------------------------------------------------------------
# Build PMNS matrix using standard parametrization
def pmns_matrix(
    theta_12=theta_12_0, theta_13=theta_13_0, theta_23=theta_23_0, delta_cp=delta_cp_0
):
    """Build the PMNS mixing matrix using standard parametrization.

    Parameters
    ----------
    theta_12, theta_13, theta_23 : float
        Mixing angles in radians.
    delta_cp : float
        CP-violating phase in radians.

    Returns
    -------
    ndarray
        Complex 3x3 PMNS matrix.
    """
    c12, s12 = np.cos(theta_12), np.sin(theta_12)
    c23, s23 = np.cos(theta_23), np.sin(theta_23)
    c13, s13 = np.cos(theta_13), np.sin(theta_13)

    U = np.zeros((3, 3), dtype=complex)
    U[0, 0] = c12 * c13
    U[0, 1] = s12 * c13
    U[0, 2] = s13 * np.exp(-1j * delta_cp)

    U[1, 0] = -s12 * c23 - c12 * s23 * s13 * np.exp(1j * delta_cp)
    U[1, 1] = c12 * c23 - s12 * s23 * s13 * np.exp(1j * delta_cp)
    U[1, 2] = s23 * c13

    U[2, 0] = s12 * s23 - c12 * c23 * s13 * np.exp(1j * delta_cp)
    U[2, 1] = -c12 * s23 - s12 * c23 * s13 * np.exp(1j * delta_cp)
    U[2, 2] = c23 * c13

    return U


# Compute oscillation probability P(ν_alpha -> ν_beta)
def oscillation_probability(
    alpha,
    beta,
    E_MeV,
    L_km,
    theta_12=theta_12_0,
    theta_13=theta_13_0,
    theta_23=theta_23_0,
    delta_cp=delta_cp_0,
    dm12=dm12_0,
    dm23=dm23_0,
):
    """Compute oscillation probability P(nu_alpha -> nu_beta).

    This function computes the oscillation probability as a function of
    baseline using standard three-flavour vacuum oscillation formulae.

    Parameters
    ----------
    alpha, beta : str
        Flavor labels ('e', 'mu', 'tau').
    E_MeV : float
        Neutrino energy in MeV.
    L_km : float
        Baseline in km.
    theta_12, theta_13, theta_23, delta_cp : float
        Mixing parameters in radians.
    dm12, dm23 : float
        Mass-squared differences in eV^2.

    Returns
    -------
    ndarray
        Oscillation probability sampled on a geometric baseline grid.
    """
    U = pmns_matrix(theta_12, theta_13, theta_23, delta_cp)
    i = iflavor[alpha]
    j = iflavor[beta]

    L_vals = np.geomspace(1, L_km, 1000)
    prob = np.zeros_like(L_vals)
    delta_m2 = [dm12, dm12 + dm23, dm23]
    mass_pairs = [(0, 1), (0, 2), (1, 2)]

    for (k, l), dm2 in zip(mass_pairs, delta_m2):
        delta = 1.267e3 * dm2 * L_vals / E_MeV
        Uik_Ujk = U[i, k] * np.conj(U[j, k])
        Uil_Ujl = U[i, l] * np.conj(U[j, l])
        re_part = np.real(Uik_Ujk * np.conj(Uil_Ujl))
        im_part = np.imag(Uik_Ujk * np.conj(Uil_Ujl))

        prob -= 4 * re_part * np.sin(delta) ** 2
        prob += 2 * im_part * np.sin(2 * delta)

    if alpha == beta:
        prob += 1.0

    return prob


# PLOTTING SETUP ------------------------------------------------------------------------
# Set up overall figure with 2 columns: left for sliders, right for plot
fig = plt.figure(figsize=(12, 6))
gs = gridspec.GridSpec(1, 2, width_ratios=[0.2, 0.8], wspace=0.1)
# left part: sliders and reset button
left_gs = gridspec.GridSpecFromSubplotSpec(
    2, 1, subplot_spec=gs[0], height_ratios=[0.9, 0.1], hspace=0.2
)
# Left top grid: 3x3 sliders
lefttop_gs = gridspec.GridSpecFromSubplotSpec(
    3, 3, subplot_spec=left_gs[0], wspace=0, hspace=0.2
)
# Create axes
ax_E = fig.add_subplot(lefttop_gs[0, 0])
# Set sliders
sE = Slider(
    ax_E, "E [MeV]", Emin[case], Emax[case], valinit=Enom[case], orientation="vertical"
)
# Left bottom part: reset button
leftbottom_gs = gridspec.GridSpecFromSubplotSpec(
    1, 1, subplot_spec=left_gs[1], wspace=0, hspace=0.2
)
ax_reset = plt.subplot(leftbottom_gs[0, 0])
reset_button = Button(ax_reset, "RESET", color="lightgray", hovercolor="0.85")
# Main plot on the right
ax_plot = fig.add_subplot(gs[1])
fig.subplots_adjust(left=0.02, right=0.98, top=0.92, bottom=0.1)


# INITIAL PLOTTING ----------------------------------------------------------------------
(line_Xe,) = ax_plot.plot(
    [], [], label=rf"$P({ilatex[alpha]} \to {ilatex[0]})$", color="blue"
)
(line_Xmu,) = ax_plot.plot(
    [], [], label=rf"$P({ilatex[alpha]} \to {ilatex[1]})$", color="orange"
)
(line_Xtau,) = ax_plot.plot(
    [], [], label=rf"$P({ilatex[alpha]} \to {ilatex[2]})$", color="green"
)
(line_all,) = ax_plot.plot(
    [], [], label=rf"$P({ilatex[alpha]} \to \nu_\ell)$", color="red", linestyle="--"
)

# Make line invisible
line_all.set_visible(False)

ax_plot.set_xscale("log")
ax_plot.set_xlim(1, Lnom[case])
ax_plot.set_ylim(0, 1.01)
ax_plot.set_xlabel("Baseline [km]")
ax_plot.set_ylabel("Oscillation Probability")
ax_plot.set_title(
    f"Oscillation Probabilities for ${ilatex[alpha]}$ at $E = {sE.val:.{decimal[case]}f}$ MeV"
)
ax_plot.legend()
ax_plot.grid()


# UPDATE FUNCTION FOR SLIDER ------------------------------------------------------------
def reset(event):
    """Reset interactive sliders to their default values.

    Parameters
    ----------
    event : matplotlib event
        Event object passed by Matplotlib widgets (ignored).
    """
    sE.reset()


def update(val):
    """Update plot lines when a slider changes.

    Parameters
    ----------
    val : float
        Slider value passed by Matplotlib (ignored here); current slider
        values are read directly.
    """
    X_e_vals = oscillation_probability(flavors[alpha], "e", sE.val, Lnom[case])
    X_mu_vals = oscillation_probability(flavors[alpha], "mu", sE.val, Lnom[case])
    X_tau_vals = oscillation_probability(flavors[alpha], "tau", sE.val, Lnom[case])
    X_all_vals = X_e_vals + X_mu_vals + X_tau_vals

    L_vals = np.geomspace(1, Lnom[case], 1000)
    line_Xe.set_data(L_vals, X_e_vals)
    line_Xmu.set_data(L_vals, X_mu_vals)
    line_Xtau.set_data(L_vals, X_tau_vals)
    line_all.set_data(L_vals, X_all_vals)

    ax_plot.set_xlim(1, Lnom[case])
    ax_plot.set_title(
        f"Oscillation Probabilities for ${ilatex[alpha]}$ at $E = {sE.val:.{decimal[case]}f}$ MeV"
    )
    fig.canvas.draw_idle()


update(sE.val)
sE.on_changed(update)
reset_button.on_clicked(reset)

from matplotlib.widgets import CheckButtons

# Create axes for check buttons (adjust position to fit your layout)
ax_check = fig.add_subplot(lefttop_gs[2, 0])
check = CheckButtons(
    ax_check,
    [r"$P(\nu_e)$", r"$P(\nu_\mu)$", r"$P(\nu_\tau)$", r"$P(\nu_\ell)$"],
    [True, True, True, False],
)


# Define visibility toggle function
def toggle_visibility(label):
    """Toggle visibility of a plotted line associated with the label.

    Parameters
    ----------
    label : str
        Label string used in the checkbox widget.
    """
    label_map = {
        r"$P(\nu_e)$": line_Xe,
        r"$P(\nu_\mu)$": line_Xmu,
        r"$P(\nu_\tau)$": line_Xtau,
        r"$P(\nu_\ell)$": line_all,
    }
    line = label_map[label]
    line.set_visible(not line.get_visible())
    fig.canvas.draw_idle()


check.on_clicked(toggle_visibility)

plt.show()
