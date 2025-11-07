"""Oscillation probability visualizer with checkboxes.

This variant of the oscillation explorer adds checkboxes to toggle which
flavor probabilities are visible, and provides a compact layout useful for
interactive demonstrations or classroom use. The script includes the usual
mixing-parameter sliders and baseline/energy controls.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Button
from matplotlib.widgets import CheckButtons

# LATEX FONT RENDERING ------------------------------------------------------------------
plt.rcParams.update(
    {
        #  "text.usetex": True,
        #  "font.family": "serif",
        #  "font.serif": ["Computer Modern Roman"],
        "font.size": 30
    }
)


# PLOT AESTHETICS -----------------------------------------------------------------------
def set_slider_fontsize(slider, fontsize=16):
    slider.label.set_fontsize(fontsize)
    slider.valtext.set_fontsize(fontsize)


# VARIABLES -----------------------------------------------------------------------------
# Flavor to index mapping
flavors = ["e", "mu", "tau"]
iflavor = {"e": 0, "mu": 1, "tau": 2}
ilatex = {
    0: r"\nu_e",
    1: r"\nu_\mu",
    2: r"\nu_\tau",
    3: r"\bar{\nu}_e",
    4: r"\bar{\nu}_\mu",
    5: r"\bar{\nu}_\tau",
}
# Flavor type: 0: electron, 1: muon, 2: tau
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
    sign=1,
):
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
        prob += sign * 2 * im_part * np.sin(2 * delta)

    if alpha == beta:
        prob += 1.0

    return prob


# PLOTTING SETUP ------------------------------------------------------------------------
# Set up overall figure with 2 columns: left for sliders, right for plot
fig = plt.figure(figsize=(12, 9))
gs = gridspec.GridSpec(1, 2, width_ratios=[0.2, 0.8], wspace=0.2)
# left part: sliders and reset button
left_gs = gridspec.GridSpecFromSubplotSpec(
    2, 1, subplot_spec=gs[0], height_ratios=[0.9, 0.1], hspace=0.2
)
# Left top grid: 3x3 sliders
lefttop_gs = gridspec.GridSpecFromSubplotSpec(3, 3, subplot_spec=left_gs[0], hspace=0.5)
# Create axes
ax_E = fig.add_subplot(lefttop_gs[0, 0])
ax_L = fig.add_subplot(lefttop_gs[0, 1])
ax_t13 = fig.add_subplot(lefttop_gs[1, 0])
ax_t12 = fig.add_subplot(lefttop_gs[1, 1])
ax_t23 = fig.add_subplot(lefttop_gs[1, 2])
ax_dcp = fig.add_subplot(lefttop_gs[2, 0])
axdm12 = fig.add_subplot(lefttop_gs[2, 1])
axdm23 = fig.add_subplot(lefttop_gs[2, 2])
# Set sliders
sE = Slider(
    ax_E, "E [MeV]", Emin[case], Emax[case], valinit=Enom[case], orientation="vertical"
)
sL = Slider(ax_L, "L [km]", 2, Lmax[case], valinit=Lnom[case], orientation="vertical")
st13deg = Slider(
    ax_t13,
    r"$\theta_{13}$ [deg]",
    0,
    90,
    valinit=np.rad2deg(theta_13_0),
    orientation="vertical",
)
st12deg = Slider(
    ax_t12,
    r"$\theta_{12}$ [deg]",
    0,
    90,
    valinit=np.rad2deg(theta_12_0),
    orientation="vertical",
)
st23deg = Slider(
    ax_t23,
    r"$\theta_{23}$ [deg]",
    0,
    90,
    valinit=np.rad2deg(theta_23_0),
    orientation="vertical",
)
sdcpdeg = Slider(
    ax_dcp,
    r"$\delta_{CP}$" f"\n[deg]",
    0,
    360,
    valinit=np.rad2deg(delta_cp_0),
    orientation="vertical",
)
sdm12e5 = Slider(
    axdm12,
    r"$\Delta m^2_{12}$" f"\n" r"[$10^{-5}$eV$^2$]",
    5,
    10,
    valinit=dm12_0 * 1e5,
    orientation="vertical",
)
sdm23e3 = Slider(
    axdm23,
    r"$\Delta m^2_{23}$" f"\n" r"[$10^{-3}$eV$^2$]",
    2,
    3,
    valinit=dm23_0 * 1e3,
    orientation="vertical",
)
sliders = [sE, sL, sdcpdeg, st13deg, st12deg, st23deg, sdm12e5, sdm23e3]
for s in sliders:
    set_slider_fontsize(s)
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
# Initial values for the plot
(line_Xe,) = ax_plot.plot(
    [], [], label=rf"$P({ilatex[alpha]} \to {ilatex[0]})$", color="orange"
)
(line_Xmu,) = ax_plot.plot(
    [], [], label=rf"$P({ilatex[alpha]} \to {ilatex[1]})$", color="blue"
)
(line_Xtau,) = ax_plot.plot(
    [], [], label=rf"$P({ilatex[alpha]} \to {ilatex[2]})$", color="green"
)
(line_all,) = ax_plot.plot(
    [], [], label=rf"$P({ilatex[alpha]} \to \nu_\ell)$", color="red", linestyle="--"
)
(line_antiXe,) = ax_plot.plot(
    [],
    [],
    label=rf"$P({ilatex[alpha+3]} \to {ilatex[3]})$",
    color="orange",
    linestyle="--",
)
(line_antiXmu,) = ax_plot.plot(
    [],
    [],
    label=rf"$P({ilatex[alpha+3]} \to {ilatex[4]})$",
    color="blue",
    linestyle="--",
)
(line_antiXtau,) = ax_plot.plot(
    [],
    [],
    label=rf"$P({ilatex[alpha+3]} \to {ilatex[5]})$",
    color="green",
    linestyle="--",
)

# Make line invisible
line_all.set_visible(False)
line_Xtau.set_visible(False)
line_antiXtau.set_visible(False)

# List of lines
lines = [
    line_Xe,
    line_Xmu,
    line_Xtau,
    line_all,
    line_antiXe,
    line_antiXmu,
    line_antiXtau,
]

# Initial setup
ax_plot.set_xscale("log")
ax_plot.set_xlim(1, sL.val)
ax_plot.set_ylim(0, 1.01)
ax_plot.set_xlabel("Baseline [km]")
ax_plot.set_ylabel("Oscillation Probability")
ax_plot.set_title(
    f"Oscillation Probabilities for ${ilatex[alpha]}$ at $E = {sE.val:.{decimal[case]}f}$ MeV"
)
ax_plot.grid()
# Only add visible lines to the initial legend
visible_lines = [
    l
    for l in [
        line_Xe,
        line_Xmu,
        line_Xtau,
        line_all,
        line_antiXe,
        line_antiXmu,
        line_antiXtau,
    ]
    if l.get_visible()
]
visible_labels = [l.get_label() for l in visible_lines]
ax_plot.legend(visible_lines, visible_labels)


# RESET SLIDERS -------------------------------------------------------------------------
def reset(event):
    sE.reset()
    sL.reset()
    st13deg.reset()
    st12deg.reset()
    st23deg.reset()
    sdcpdeg.reset()
    sdm12e5.reset()
    sdm23e3.reset()


reset_button.on_clicked(reset)


# UPDATE FUNCTION FOR SLIDER ------------------------------------------------------------
# Update function for sliders
def update(val):
    st12 = np.deg2rad(st12deg.val)
    st13 = np.deg2rad(st13deg.val)
    st23 = np.deg2rad(st23deg.val)
    sdcp = np.deg2rad(sdcpdeg.val)
    sdm12 = sdm12e5.val * 1e-5
    sdm23 = sdm23e3.val * 1e-3

    X_e_vals = oscillation_probability(
        flavors[alpha], "e", sE.val, sL.val, st12, st13, st23, sdcp, sdm12, sdm23, +1
    )
    X_mu_vals = oscillation_probability(
        flavors[alpha], "mu", sE.val, sL.val, st12, st13, st23, sdcp, sdm12, sdm23, +1
    )
    X_tau_vals = oscillation_probability(
        flavors[alpha], "tau", sE.val, sL.val, st12, st13, st23, sdcp, sdm12, sdm23, +1
    )
    X_all_vals = X_e_vals + X_mu_vals + X_tau_vals
    antiX_e_vals = oscillation_probability(
        flavors[alpha], "e", sE.val, sL.val, st12, st13, st23, sdcp, sdm12, sdm23, -1
    )
    antiX_mu_vals = oscillation_probability(
        flavors[alpha], "mu", sE.val, sL.val, st12, st13, st23, sdcp, sdm12, sdm23, -1
    )
    antiX_tau_vals = oscillation_probability(
        flavors[alpha], "tau", sE.val, sL.val, st12, st13, st23, sdcp, sdm12, sdm23, -1
    )

    L_vals = np.geomspace(1, sL.val, 1000)
    line_Xe.set_data(L_vals, X_e_vals)
    line_Xmu.set_data(L_vals, X_mu_vals)
    line_Xtau.set_data(L_vals, X_tau_vals)
    line_all.set_data(L_vals, X_all_vals)
    line_antiXe.set_data(L_vals, antiX_e_vals)
    line_antiXmu.set_data(L_vals, antiX_mu_vals)
    line_antiXtau.set_data(L_vals, antiX_tau_vals)

    ax_plot.set_xlim(1, sL.val)
    ax_plot.set_title(
        f"Oscillation probabilities for ${ilatex[alpha]}$ at $E = {sE.val:.{decimal[case]}f}$ MeV"
    )
    fig.canvas.draw_idle()


# Connect sliders to update function
update(sE.val)
for slider in sliders:
    slider.on_changed(update)

# HIDE/SHOW LINES -----------------------------------------------------------------------
ax_check = fig.add_subplot(lefttop_gs[0, 2])
check = CheckButtons(
    ax=ax_check,
    labels=[
        r"$\nu_e$",
        r"$\nu_\mu$",
        r"$\nu_\tau$",
        r"$\nu_\ell$",
        r"$\bar{\nu}_e$",
        r"$\bar{\nu}_\mu$",
        r"$\bar{\nu}_\tau$",
    ],
    actives=[l.get_visible() for l in lines],
    frame_props={"facecolor": "lightgray"},
    check_props={"facecolor": [l.get_color() for l in lines]},
)

for label in check.labels:
    label.set_fontsize(16)


# Define visibility toggle function
def toggle_visibility(label):
    label_map = {
        r"$\nu_e$": line_Xe,
        r"$\nu_\mu$": line_Xmu,
        r"$\nu_\tau$": line_Xtau,
        r"$\nu_\ell$": line_all,
        r"$\bar{\nu}_e$": line_antiXe,
        r"$\bar{\nu}_\mu$": line_antiXmu,
        r"$\bar{\nu}_\tau$": line_antiXtau,
    }
    line = label_map[label]
    line.set_visible(not line.get_visible())
    # Update legend to show only visible lines
    visible_lines = [
        l
        for l in [
            line_Xe,
            line_Xmu,
            line_Xtau,
            line_all,
            line_antiXe,
            line_antiXmu,
            line_antiXtau,
        ]
        if l.get_visible()
    ]
    visible_labels = [l.get_label() for l in visible_lines]
    ax_plot.legend(visible_lines, visible_labels)
    # Redraw
    fig.canvas.draw_idle()


check.on_clicked(toggle_visibility)

# SHOW PLOT -----------------------------------------------------------------------------
plt.show()
