"""Interactive reactor-focused oscillation probability explorer.

This script provides interactive visualisations tailored for reactor and
short-baseline experiments. It includes sliders for energy, baseline and
mixing parameters, and draws markers for representative experiments like
JUNO and KamLAND. Use it to generate thesis figures or pedagogical plots.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Button

# LATEX FONT RENDERING ------------------------------------------------------------------
plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "font.size": 35,
    }
)


# PLOT AESTHETICS -----------------------------------------------------------------------
def set_slider_fontsize(slider, fontsize=25):
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
alpha = 0
#  Cases: reactor: 0, accelerator: 1
case = 0
decimal = [1, 0]
Emin = [0.1, 100]
Enom = [4, 600]
Emax = [10, 4e3]
nLpoints = 2000
Lmin = [0.2, 1]
Lnom = [300, 2e3]
Lmax = [1e3, 5e3]
# NuFIT 6.0 (Normal Ordering) best-fit values
theta_12_0 = np.deg2rad(33.68)
theta_23_0 = np.deg2rad(43.3)
theta_13_0 = np.deg2rad(8.56)
delta_cp_0 = np.deg2rad(212)
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
    Lmin_km,
    Lmax_km,
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

    L_vals = np.geomspace(Lmin_km, Lmax_km, nLpoints)
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
fig.subplots_adjust(left=0.02, right=0.98, top=0.97, bottom=0.03)

# Create gridspec for layout
gs = gridspec.GridSpec(1, 2, width_ratios=[0.2, 0.8], wspace=0.2)
# Left part: sliders and reset button
buttons_gs = gridspec.GridSpecFromSubplotSpec(
    2, 1, subplot_spec=gs[0], height_ratios=[0.9, 0.1], hspace=0.2
)
# Sliders: E and L, and parameters
sliders_gs = gridspec.GridSpecFromSubplotSpec(
    2, 1, subplot_spec=buttons_gs[0], height_ratios=[0.3, 0.7], hspace=0.4
)
# Top grid: E and L sliders
EandL_gs = gridspec.GridSpecFromSubplotSpec(
    1, 2, subplot_spec=sliders_gs[0], wspace=0.2
)
ax_E = fig.add_subplot(EandL_gs[0, 0])
ax_Epos = ax_E.get_position()
ax_E.set_position([ax_Epos.x0, ax_Epos.y0, ax_Epos.width, ax_Epos.height * 0.9])
ax_L = fig.add_subplot(EandL_gs[0, 1])
ax_Lpos = ax_L.get_position()
ax_L.set_position([ax_Lpos.x0, ax_Lpos.y0, ax_Lpos.width, ax_Lpos.height * 0.9])
sE = Slider(
    ax_E, "E [MeV]", Emin[case], Emax[case], valinit=Enom[case], orientation="vertical"
)
sL = Slider(
    ax_L, "L [km]", Lmin[case], Lmax[case], valinit=Lnom[case], orientation="vertical"
)

# Parameters sliders
params_gs = gridspec.GridSpecFromSubplotSpec(
    2, 3, subplot_spec=sliders_gs[1], wspace=0.2, hspace=1.0
)
ax_t13 = fig.add_subplot(params_gs[0, 0])
ax_t12 = fig.add_subplot(params_gs[0, 1])
ax_t23 = fig.add_subplot(params_gs[0, 2])
axdm12 = fig.add_subplot(params_gs[1, 0])
ax_dcp = fig.add_subplot(params_gs[1, 1])
axdm23 = fig.add_subplot(params_gs[1, 2])
# Set sliders
st13deg = Slider(
    ax_t13,
    r"$\theta_{13}$ [°]",
    0,
    90,
    valinit=np.rad2deg(theta_13_0),
    orientation="vertical",
)
st12deg = Slider(
    ax_t12,
    r"$\theta_{12}$ [°]",
    0,
    90,
    valinit=np.rad2deg(theta_12_0),
    orientation="vertical",
)
st23deg = Slider(
    ax_t23,
    r"$\theta_{23}$ [°]",
    0,
    90,
    valinit=np.rad2deg(theta_23_0),
    orientation="vertical",
)
sdcpdeg = Slider(
    ax_dcp,
    r"$\delta_{CP}$ [°]",
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
    1, 1, subplot_spec=buttons_gs[1], wspace=0, hspace=0.2
)
ax_reset = plt.subplot(leftbottom_gs[0, 0])
reset_button = Button(ax_reset, "RESET", color="lightgray", hovercolor="0.85")
# Main plot on the right
ax_plot = fig.add_subplot(gs[1])
ax_pos = ax_plot.get_position()
ax_plot.set_position(
    [ax_pos.x0, ax_pos.y0 + ax_pos.height * 0.1, ax_pos.width, ax_pos.height * 0.9]
)


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
    [], [], label=rf"$P({ilatex[alpha+3]} \to {ilatex[3]})$", color="orange"
)
(line_antiXmu,) = ax_plot.plot(
    [], [], label=rf"$P({ilatex[alpha+3]} \to {ilatex[4]})$", color="blue"
)
(line_antiXtau,) = ax_plot.plot(
    [], [], label=rf"$P({ilatex[alpha+3]} \to {ilatex[5]})$", color="green"
)

reactors_L = 1.8  # km
(reactors_marker,) = ax_plot.plot(
    [],
    [],
    marker="^",
    color="red",
    markersize=20,
    markeredgecolor="black",
    markeredgewidth=2,
    label="Reactor exps.",
    linestyle="None",
)

juno_L = 65  # km
(juno_marker,) = ax_plot.plot(
    [],
    [],
    marker="s",
    color="red",
    markersize=20,
    markeredgecolor="black",
    markeredgewidth=2,
    label="JUNO",
    linestyle="None",
)

# Experiment markers
kamland_L = 180  # km
(kamland_marker,) = ax_plot.plot(
    [],
    [],
    marker="o",
    color="red",
    markersize=20,
    markeredgecolor="black",
    markeredgewidth=2,
    label="KamLAND",
    linestyle="None",
)

# Make line invisible
line_all.set_visible(False)
line_Xe.set_visible(False)
# line_antiXe.set_visible(False)
line_Xmu.set_visible(False)
# line_antiXmu.set_visible(False)
line_Xtau.set_visible(False)
# line_antiXtau.set_visible(False)

# List of lines
lines = [
    line_Xe,
    line_Xmu,
    line_Xtau,
    line_all,
    line_antiXe,
    line_antiXmu,
    line_antiXtau,
    reactors_marker,
    juno_marker,
    kamland_marker,
]

# Initial setup
ax_plot.set_xscale("log")
ax_plot.set_ylim(0, 1.07)
ax_plot.set_xlabel("Baseline [km]")
ax_plot.set_ylabel("Oscillation Probability")
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
        flavors[alpha],
        "e",
        sE.val,
        Lmin[case],
        sL.val,
        st12,
        st13,
        st23,
        sdcp,
        sdm12,
        sdm23,
        +1,
    )
    X_mu_vals = oscillation_probability(
        flavors[alpha],
        "mu",
        sE.val,
        Lmin[case],
        sL.val,
        st12,
        st13,
        st23,
        sdcp,
        sdm12,
        sdm23,
        +1,
    )
    X_tau_vals = oscillation_probability(
        flavors[alpha],
        "tau",
        sE.val,
        Lmin[case],
        sL.val,
        st12,
        st13,
        st23,
        sdcp,
        sdm12,
        sdm23,
        +1,
    )
    X_all_vals = X_e_vals + X_mu_vals + X_tau_vals
    antiX_e_vals = oscillation_probability(
        flavors[alpha],
        "e",
        sE.val,
        Lmin[case],
        sL.val,
        st12,
        st13,
        st23,
        sdcp,
        sdm12,
        sdm23,
        -1,
    )
    antiX_mu_vals = oscillation_probability(
        flavors[alpha],
        "mu",
        sE.val,
        Lmin[case],
        sL.val,
        st12,
        st13,
        st23,
        sdcp,
        sdm12,
        sdm23,
        -1,
    )
    antiX_tau_vals = oscillation_probability(
        flavors[alpha],
        "tau",
        sE.val,
        Lmin[case],
        sL.val,
        st12,
        st13,
        st23,
        sdcp,
        sdm12,
        sdm23,
        -1,
    )

    L_vals = np.geomspace(Lmin[case], sL.val, nLpoints)
    line_Xe.set_data(L_vals, X_e_vals)
    line_Xmu.set_data(L_vals, X_mu_vals)
    line_Xtau.set_data(L_vals, X_tau_vals)
    line_all.set_data(L_vals, X_all_vals)
    line_antiXe.set_data(L_vals, antiX_e_vals)
    line_antiXmu.set_data(L_vals, antiX_mu_vals)
    line_antiXtau.set_data(L_vals, antiX_tau_vals)

    # Update experiment markers
    kamland_y = oscillation_probability(
        flavors[alpha],
        "e",
        Enom[case],
        Lmin[case],
        kamland_L,
        st12,
        st13,
        st23,
        sdcp,
        sdm12,
        sdm23,
        +1,
    )[-1]
    kamland_marker.set_data([kamland_L], [kamland_y])

    juno_y = oscillation_probability(
        flavors[alpha],
        "e",
        Enom[case],
        Lmin[case],
        juno_L,
        st12,
        st13,
        st23,
        sdcp,
        sdm12,
        sdm23,
        -1,
    )[-1]
    juno_marker.set_data([juno_L], [juno_y])

    reactors_y = oscillation_probability(
        flavors[alpha],
        "e",
        Enom[case],
        Lmin[case],
        reactors_L,
        st12,
        st13,
        st23,
        sdcp,
        sdm12,
        sdm23,
        -1,
    )[-1]
    reactors_marker.set_data([reactors_L], [reactors_y])

    ax_plot.set_xlim(Lmin[case], sL.val)
    # Update legend to show only visible lines and experiment markers
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
            reactors_marker,
            juno_marker,
            kamland_marker,
        ]
        if l.get_visible()
    ]
    visible_labels = [l.get_label() for l in visible_lines]
    ax_plot.legend(visible_lines, visible_labels)

    # Remove previous arrow if it exists
    if hasattr(ax_plot, "_juno_arrow_hori"):
        ax_plot._juno_arrow_hori.remove()
        ax_plot._juno_arrow_hori_label.remove()
    if hasattr(ax_plot, "_juno_arrow_vert"):
        ax_plot._juno_arrow_vert.remove()
        ax_plot._juno_arrow_vert_label.remove()
    if hasattr(ax_plot, "_reactors_arrow_hori"):
        ax_plot._reactors_arrow_hori.remove()
        ax_plot._reactors_arrow_hori_label.remove()
    if hasattr(ax_plot, "_reactors_arrow_vert"):
        ax_plot._reactors_arrow_vert.remove()
        ax_plot._reactors_arrow_vert_label.remove()

    # parameters 21
    # Draw horizontal arrow between Lmin[case] and juno_L at y-level of JUNO marker
    ax_plot._juno_arrow_hori = ax_plot.annotate(
        "",
        xy=(juno_L * 0.9, juno_y),
        xytext=(Lmin[case], juno_y),
        arrowprops=dict(arrowstyle="<->", color="black", lw=2),
        annotation_clip=False,
    )
    ax_plot._juno_arrow_hori_label = ax_plot.text(
        np.sqrt(Lmin[case] * juno_L),
        juno_y - 0.07,
        r"$\Delta m^2_{21}$",
        ha="center",
        va="bottom",
        fontsize=30,
        color="black",
    )
    # Draw vertical line between juno_y and 1 at juno_L
    ax_plot._juno_arrow_vert = ax_plot.annotate(
        "",
        xy=(juno_L, juno_y + 0.02),
        xytext=(juno_L, 1),
        arrowprops=dict(arrowstyle="<->", color="black", lw=2),
        annotation_clip=False,
    )
    ax_plot._juno_arrow_vert_label = ax_plot.text(
        juno_L * 0.7,
        (juno_y + 1) / 2 + 0.12,
        r"$\sin^2\ 2\theta_{12}$",
        ha="left",
        va="center",
        fontsize=30,
        color="black",
        rotation=90,
    )

    # parameters 31
    # Draw horizontal arrow between Lmin[case] and reactors_L at y-level of reactors marker
    juno_y = reactors_y
    ax_plot._reactors_arrow_hori = ax_plot.annotate(
        "",
        xy=(reactors_L * 0.9, reactors_y),
        xytext=(Lmin[case], reactors_y),
        arrowprops=dict(arrowstyle="<->", color="black", lw=2),
        annotation_clip=False,
    )
    ax_plot._reactors_arrow_hori_label = ax_plot.text(
        np.sqrt(Lmin[case] * reactors_L),
        reactors_y - 0.07,
        r"$\Delta m^2_{31}$",
        ha="center",
        va="bottom",
        fontsize=30,
        color="black",
    )
    # Draw vertical line between reactors_y and 1 at reactors_L
    ax_plot._reactors_arrow_vert = ax_plot.annotate(
        "",
        xy=(reactors_L, reactors_y + 0.015),
        xytext=(reactors_L, 1),
        arrowprops=dict(arrowstyle="<->", color="black", lw=2),
        annotation_clip=False,
    )
    ax_plot._reactors_arrow_vert_label = ax_plot.text(
        reactors_L * 0.7,
        1.03,
        r"$\sin^2\ 2\theta_{13}$",
        ha="left",
        va="center",
        fontsize=30,
        color="black",
    )

    fig.canvas.draw_idle()


# Connect sliders to update function
update(sE.val)
for slider in sliders:
    slider.on_changed(update)


# HIDE/SHOW LINES -----------------------------------------------------------------------
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


plt.show()
