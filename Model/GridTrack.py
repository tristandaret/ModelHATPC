"""Interactive track (lineic deposit) visualization for the HATPC model.

This module creates an interactive matplotlib figure to visualize the
time-dependent response of grid cells to a lineic (track-like) charge
deposition. It uses helper functions from the `Headers` package:
`compute_line_params` to convert angle/impact to line parameters and
`Compute1D` to get the cell-integrated time response for a line deposit.

Controls
--------
- Vertical sliders adjust RC, drift length `z`, track angle `phi` and impact
    parameter `d`. Adjusting the sliders updates the map and recomputes every
    cell trace.

Run
---
Execute from the repository root::

        python Model/GridTrack.py

See `Headers/ModelUtils.py` and `Headers/GridUtils.py` for model details and
shared plotting defaults.
"""

from sys import path

path.append("Headers/")
import GeometryUtils as geo
import GridUtils as gutils
import matplotlib.pyplot as plt
import ModelUtils as mutils
from matplotlib.widgets import Slider

# ========================= Sliders Setup =========================
dySlider = 0.22
left1 = 0.03
left2 = left1 + 0.07
low1 = 0.16
low2 = low1 + dySlider + 0.13

slider_defs = {
    "RC": dict(
        ax=[left1, low1, 0.025, dySlider],
        label="RC\n[ns/mm$^{2}$]",
        valmin=1,
        valmax=250,
        valinit=mutils.RC,
        step=1,
        color="C4",
    ),
    "z": dict(
        ax=[left2, low1, 0.025, dySlider],
        label="Drift\n[mm]",
        valmin=0,
        valmax=1000,
        valinit=mutils.z,
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
        color=gutils.varcolor,
    ),
    "d": dict(
        ax=[left2, low2, 0.025, dySlider],
        label="impact\n[mm]",
        valmin=-geo.diag / 2,
        valmax=geo.diag / 2,
        valinit=0,
        step=0.1,
        color=gutils.varcolor,
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

# Reset button
gutils.button.on_clicked(lambda event: [s.reset() for s in sliders.values()])


# ========================= Map Plot Setup =========================
m, q, phi_rad = geo.compute_line_params(sliders["phi"].val, sliders["d"].val)

axMap = gutils.fig.add_axes((0.035, 0.82, 0.085, 0.12))
axMap.set_xlim(geo.xleft, geo.xleft + geo.nX * geo.xwidth)
axMap.set_ylim(geo.ylow, geo.ylow + geo.nY * geo.ywidth)
axMap.set_title("Track position", fontsize=20)
axMap.set_xticks(geo.np.arange(geo.xleft, geo.xleft + geo.nX * geo.xwidth, geo.xwidth))
axMap.set_yticks(geo.np.arange(geo.ylow, geo.ylow + geo.nY * geo.ywidth, geo.ywidth))
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

v_points = geo.np.linspace(-6, 10, 50)
x_line = (
    geo.np.cos(phi_rad) ** 2 * geo.xc
    + geo.np.cos(phi_rad) * geo.np.sin(phi_rad) * geo.yc
) * v_points + geo.xc
y_line = (
    geo.np.sin(phi_rad) ** 2 * geo.yc
    + geo.np.cos(phi_rad) * geo.np.sin(phi_rad) * geo.xc
) * v_points + geo.yc
(map_line,) = axMap.plot(x_line, y_line, "red")


# ========================= Subplots Initialization =========================
for iX in range(geo.nX):
    for iY in range(geo.nY):
        ax = gutils.axs[geo.nY - 1 - iY, iX]
        x0 = geo.xleft + iX * geo.xwidth
        x1 = x0 + geo.xwidth
        y0 = geo.ylow + iY * geo.ywidth
        y1 = y0 + geo.ywidth
        res = mutils.Compute1D(
            mutils.t, m, q, x0, x1, y0, y1, mutils.RC, mutils.z, gutils.vartype
        )
        if res is None:
            res = geo.np.zeros_like(mutils.t)
        sig = gutils.scalefactor * res[: mutils.t.size]
        line = ax.plot(
            mutils.t / gutils.timescale, sig, lw=8 - geo.nY, color=gutils.varcolor
        )
        txt = ax.text(
            0.96,
            0.93,
            f"{gutils.dim} = {geo.np.max(sig):.0f} {gutils.unit}\n T$_{{max}}$ = {mutils.t[geo.np.argmax(sig)]/gutils.timescale:.0f} {gutils.timeunit}",
            ha="right",
            va="top",
            transform=ax.transAxes,
            fontsize=25 - 2 * geo.nY,
            bbox=dict(boxstyle="round", facecolor=gutils.varcolor),
            color=gutils.varforeground,
        )
        gutils.lines.append(line)
        gutils.texts.append(txt)

        if iY == 0:
            ax.set_xlabel(f"Time ({gutils.timeunit})", fontsize=25 - geo.nY)
            ax.tick_params(axis="x", labelsize=20 - geo.nY)
        if iX == 0:
            ax.set_ylabel(f"{gutils.ylabel}", fontsize=25 - geo.nX)
            ax.tick_params(axis="y", labelsize=20 - geo.nX)
        ax.grid()
        ax.set_xlim(gutils.varxminplot, gutils.varxmaxplot / gutils.timescale)
        ax.set_ylim(gutils.varyminplot, gutils.varymaxplot)


# ========================= Slider Callback Update =========================
def update(val):
    """Matplotlib slider callback to update the track display.

    Parameters
    ----------
    val : float
        Callback argument passed by matplotlib (ignored). The function reads
        values from the global `sliders` dictionary.

    Effects
    -------
    - Recomputes line parameters (m, q) from `phi` and `d`.
    - Updates the mapping line shown on the map axes.
    - Recomputes time traces for each cell using `Compute1D` and updates the
      plotted lines and text summaries.
    """
    phi = sliders["phi"].val
    if abs(phi) < 1e-6:
        phi = 1e-6 * geo.np.sign(phi or 1)

    m, q, phi_rad = geo.compute_line_params(phi, sliders["d"].val)

    v_x = geo.np.linspace(
        geo.xleft - geo.xwidth, geo.xleft + (geo.nX + 1) * geo.xwidth, 10 * geo.nX
    )
    v_y = geo.np.tan(phi_rad) * v_x + q
    map_line.set_data(v_x, v_y)

    for iX in range(geo.nX):
        for iY in range(geo.nY):
            idx = iX * geo.nY + iY
            res = mutils.Compute1D(
                mutils.t,
                m,
                q,
                geo.xleft + iX * geo.xwidth,
                geo.xleft + (iX + 1) * geo.xwidth,
                geo.ylow + iY * geo.ywidth,
                geo.ylow + (iY + 1) * geo.ywidth,
                sliders["RC"].val,
                sliders["z"].val,
                gutils.vartype,
            )
            if res is None:
                res = geo.np.zeros_like(mutils.t)
            sig = gutils.scalefactor * res[: mutils.t.size]
            gutils.lines[idx][0].set_ydata(sig)
            gutils.texts[idx].set_text(
                f"{gutils.dim} = {geo.np.max(sig):.0f} {gutils.unit}\n T$_{{max}}$ = {mutils.t[geo.np.argmax(sig)]/gutils.timescale:.0f} {gutils.timeunit}"
            )
    gutils.fig.canvas.draw_idle()


for s in sliders.values():
    s.on_changed(update)

# ========================= Show Plot =========================
plt.show()
