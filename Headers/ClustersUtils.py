from sys import path

path.append("Headers/")
from ModelUtils import *
from GeometryUtils import *
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

plt.rcParams.update(
    {
        "font.size": 25,
        "axes.labelpad": 10,
    }
)

# ========================= Plot Setup =========================
fig, axs = plt.subplots(figsize=(12, 9))
fig.subplots_adjust(left=0.35, right=0.96, bottom=0.12, top=0.98)

# ========================= Sliders =========================
dySlider = 0.22
left1 = 0.05
left2 = left1 + 0.11
low1 = 0.16
low2 = low1 + dySlider + 0.13

slider_defs = {
    "RC": dict(
        ax=[left1, low1, 0.025, dySlider],
        label="RC\n[ns/mm$^{2}$]",
        valmin=1,
        valmax=250,
        valinit=RC,
        step=1,
        color="C4",
    ),
    "z": dict(
        ax=[left2, low1, 0.025, dySlider],
        label="Drift\n[mm]",
        valmin=0,
        valmax=1000,
        valinit=z,
        step=1,
        color="C4",
    ),
    "phi": dict(
        ax=[left1, low2, 0.025, dySlider],
        label="angle\n[°]",
        valmin=1e-6,
        valmax=90 - 1e-6,
        valinit=42,
        step=0.1,
        color="C1",
    ),
    "d": dict(
        ax=[left2, low2, 0.025, dySlider],
        label="impact\n[mm]",
        valmin=-diag / 2,
        valmax=diag / 2,
        valinit=0,
        step=0.1,
        color="C1",
    ),
}

sliders = {}
for key, cfg in slider_defs.items():
    ax = plt.axes(cfg["ax"])
    sliders[key] = Slider(
        ax=ax,
        label=cfg["label"],
        valmin=cfg["valmin"],
        valmax=cfg["valmax"],
        valinit=cfg["valinit"],
        valstep=cfg["step"],
        orientation="vertical",
        color=cfg["color"],
    )
    sliders[key].label.set_fontsize(20)
    sliders[key].valtext.set_fontsize(20)

# ========================= Reset Button =========================
resetax = fig.add_axes([0.06, 0.05, 0.11, 0.06])
button = Button(resetax, "Reset", hovercolor="0.975")
button.on_clicked(lambda event: [s.reset() for s in sliders.values()])

# ========================= Map =========================
m, q, phi_rad = compute_line_params(sliders["phi"].val, sliders["d"].val)

axMap = fig.add_axes([0.06, 0.82, 0.11, 0.12])
axMap.set_xlim(xleft, xleft + nX * xwidth)
axMap.set_ylim(ylow, ylow + nY * ywidth)
axMap.set_title("Track position", fontsize=20)
axMap.set_xticks(np.arange(xleft, xleft + nX * xwidth, xwidth))
axMap.set_yticks(np.arange(ylow, ylow + nY * ywidth, ywidth))
axMap.tick_params(
    which="both",
    labelleft=False,
    labelbottom=False,
    bottom=False,
    top=False,
    left=False,
    right=False,
)
axMap.grid()

v_points = np.linspace(-6, 10, 50)
x_line = (
    np.cos(phi_rad) ** 2 * xc + np.cos(phi_rad) * np.sin(phi_rad) * yc
) * v_points + xc
y_line = (
    np.sin(phi_rad) ** 2 * yc + np.cos(phi_rad) * np.sin(phi_rad) * xc
) * v_points + yc
(map_line,) = axMap.plot(x_line, y_line, "red")
