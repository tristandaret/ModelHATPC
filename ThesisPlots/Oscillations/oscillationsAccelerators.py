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
alpha = 1
LEmin = 20
LEmax = 2000
nLpoints = 2000
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
    LE_kmGeV,
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

    prob = 0.0
    delta_m2 = [dm12, dm12 + dm23, dm23]
    mass_pairs = [(0, 1), (0, 2), (1, 2)]

    for (k, l), dm2 in zip(mass_pairs, delta_m2):
        delta = 1.267 * dm2 * LE_kmGeV
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

# Parameters sliders
params_gs = gridspec.GridSpecFromSubplotSpec(
    2, 3, subplot_spec=buttons_gs[0], wspace=0.2, hspace=0.5
)
ax_t13 = fig.add_subplot(params_gs[0, 0])
ax_t12 = fig.add_subplot(params_gs[0, 1])
ax_t23 = fig.add_subplot(params_gs[0, 2])
axdm12 = fig.add_subplot(params_gs[1, 0])
ax_dcp = fig.add_subplot(params_gs[1, 1])
axdm23 = fig.add_subplot(params_gs[1, 2])
for ax in [ax_t13, ax_t12, ax_t23]:
    axpos = ax.get_position()
    ax.set_position([axpos.x0, axpos.y0, axpos.width, axpos.height * 0.9])

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
sliders = [sdcpdeg, st13deg, st12deg, st23deg, sdm12e5, sdm23e3]
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
    [ax_pos.x0, ax_pos.y0 + ax_pos.height * 0.1, ax_pos.width, ax_pos.height * 0.87]
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

# Experiment markers
T2K_LE = 295 / 0.6  # km/GeV
ax_plot.text(
    T2K_LE,
    1.05,
    "T2K",
    color="red",
    fontsize=30,
    ha="center",
    va="top",
    transform=ax_plot.get_xaxis_transform(),
)

# Make line invisible
# line_antiXe.set_visible(False)
# line_antiXmu.set_visible(False)
# line_Xtau.set_visible(False)
# line_antiXtau.set_visible(False)

# List of lines
lines = [line_Xe, line_Xmu, line_Xtau, line_antiXe, line_antiXmu, line_antiXtau]

# Initial setup
ax_plot.set_xscale("log")
ax_plot.set_ylim(0, 1)
ax_plot.set_xlabel("Baseline/Energy [km/GeV]")
ax_plot.set_ylabel("Oscillation Probability")
ax_plot.grid()
# Only add visible lines to the initial legend
visible_lines = [
    l
    for l in [
        line_Xe,
        line_Xmu,
        line_Xtau,
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
    st13deg.reset()
    st12deg.reset()
    st23deg.reset()
    sdcpdeg.reset()
    sdm12e5.reset()
    sdm23e3.reset()


reset_button.on_clicked(reset)


# UPDATE FUNCTION FOR SLIDER ------------------------------------------------------------
# Update function for sliders
LE_vals = np.linspace(LEmin, LEmax, nLpoints)


def update(val):
    st12 = np.deg2rad(st12deg.val)
    st13 = np.deg2rad(st13deg.val)
    st23 = np.deg2rad(st23deg.val)
    sdcp = np.deg2rad(sdcpdeg.val)
    sdm12 = sdm12e5.val * 1e-5
    sdm23 = sdm23e3.val * 1e-3

    line_Xe.set_data(
        LE_vals,
        oscillation_probability(
            flavors[alpha],
            "e",
            LE_vals,
            st12,
            st13,
            st23,
            sdcp,
            sdm12,
            sdm23,
            +1,
        ),
    )
    line_Xmu.set_data(
        LE_vals,
        oscillation_probability(
            flavors[alpha],
            "mu",
            LE_vals,
            st12,
            st13,
            st23,
            sdcp,
            sdm12,
            sdm23,
            +1,
        ),
    )
    line_Xtau.set_data(
        LE_vals,
        oscillation_probability(
            flavors[alpha],
            "tau",
            LE_vals,
            st12,
            st13,
            st23,
            sdcp,
            sdm12,
            sdm23,
            +1,
        ),
    )
    line_antiXe.set_data(
        LE_vals,
        oscillation_probability(
            flavors[alpha],
            "e",
            LE_vals,
            st12,
            st13,
            st23,
            sdcp,
            sdm12,
            sdm23,
            -1,
        ),
    )
    line_antiXmu.set_data(
        LE_vals,
        oscillation_probability(
            flavors[alpha],
            "mu",
            LE_vals,
            st12,
            st13,
            st23,
            sdcp,
            sdm12,
            sdm23,
            -1,
        ),
    )
    line_antiXtau.set_data(
        LE_vals,
        oscillation_probability(
            flavors[alpha],
            "tau",
            LE_vals,
            st12,
            st13,
            st23,
            sdcp,
            sdm12,
            sdm23,
            -1,
        ),
    )

    ax_plot.set_xlim(LEmin, LEmax)
    # Update legend to show only visible lines and experiment markers
    visible_lines = [
        l
        for l in [
            line_Xe,
            line_Xmu,
            line_Xtau,
            line_antiXe,
            line_antiXmu,
            line_antiXtau,
        ]
        if l.get_visible()
    ]
    visible_labels = [l.get_label() for l in visible_lines]
    ax_plot.legend(visible_lines, visible_labels)

    # Draw parameter arrows
    # Remove old arrows
    if hasattr(ax_plot, "theta13_arrow"):
        ax_plot.theta13_arrow.remove()
        ax_plot.theta13_label.remove()
        ax_plot.theta23_arrow.remove()
        ax_plot.theta23_label.remove()
        ax_plot.dm223_arrow.remove()
        ax_plot.dm223_label.remove()

    # Add new arrows
    T2K_ymu = oscillation_probability(
        flavors[alpha],
        "mu",
        T2K_LE,
        st12,
        st13,
        st23,
        sdcp,
        sdm12,
        sdm23,
        +1,
    )
    T2K_ye = oscillation_probability(
        flavors[alpha],
        "e",
        T2K_LE,
        st12,
        st13,
        st23,
        sdcp,
        sdm12,
        sdm23,
        +1,
    )
    ax_plot.theta13_arrow = ax_plot.annotate(
        "",
        xy=(T2K_LE * 1.08, T2K_ye + 0.007),
        xytext=(T2K_LE * 1.08, 0 - 0.005),
        arrowprops=dict(arrowstyle="<->", color="black", lw=2),
        annotation_clip=False,
    )
    ax_plot.theta13_label = ax_plot.text(
        T2K_LE * 1.08,
        T2K_ye - 0.08,
        r"$\sin^2 2 \theta_{13}$",
        fontsize=30,
        ha="center",
        va="center",
    )

    ax_plot.theta23_arrow = ax_plot.annotate(
        "",
        xy=(T2K_LE, 1 + 0.007),
        xytext=(T2K_LE, T2K_ymu - 0.005),
        arrowprops=dict(arrowstyle="<->", color="black", lw=2),
        annotation_clip=False,
    )
    ax_plot.theta23_label = ax_plot.text(
        T2K_LE * 0.9,
        (1 + T2K_ymu) / 2,
        r"$\sin^2 2 \theta_{23}$",
        fontsize=30,
        ha="center",
        va="center",
        rotation=90,
    )

    ax_plot.dm223_arrow = ax_plot.annotate(
        "",
        xy=(LEmin, T2K_ye),
        xytext=(T2K_LE, T2K_ye),
        arrowprops=dict(arrowstyle="<->", color="black", lw=2),
        annotation_clip=False,
    )
    ax_plot.dm223_label = ax_plot.text(
        2 * np.sqrt(LEmin * T2K_LE),
        T2K_ye + 0.05,
        r"$\Delta m^2_{23}$",
        fontsize=30,
        ha="center",
        va="center",
    )

    fig.canvas.draw_idle()


# Connect sliders to update function
update(sliders[0])
for slider in sliders:
    slider.on_changed(update)


# HIDE/SHOW LINES -----------------------------------------------------------------------
# Define visibility toggle function
def toggle_visibility(label):
    label_map = {
        r"$\nu_e$": line_Xe,
        r"$\nu_\mu$": line_Xmu,
        r"$\nu_\tau$": line_Xtau,
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
