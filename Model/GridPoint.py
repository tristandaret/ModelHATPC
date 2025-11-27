"""Interactive point-deposit visualization for the HATPC lineic-charge model.

This module builds an interactive matplotlib figure which shows the time
response (signal/charge/current) in each cell of a grid for a point-like
charge deposit. It was originally created while solving the telegrapher's
equation for lineic deposits and extending the point-like solution.

The script expects helper functions and global parameters to be provided by
the `Headers` package (notably `ModelUtils` and `GridUtils`). Typical globals
used here include `t`, `timescale`, `nX`, `nY`, `xleft`, `ylow`, `xwidth`,
`ywidth`, `xc`, `yc`, `RC`, `z`, `Dt`, `scalefactor`, and plotting defaults.

Usage
-----
Run from the repository root::

        python Model/GridPoint.py

This opens an interactive window with vertical sliders for RC, drift (z), and
the drop x/y position. Moving sliders re-computes the per-cell traces using
`Compute0D` from `ModelUtils` and updates the plots.

Notes
-----
- Docstrings are intentionally concise; see `Headers/ModelUtils.py` for the
    physical model (Compute0D) and parameter definitions.
"""

from sys import path

path.append("Headers/")
import GeometryUtils as geo
import GridUtils as gu
import matplotlib.pyplot as plt
import ModelUtils as mu
from matplotlib.widgets import Slider

# Update axes range values
if gu.vartype == "Charge":
    gu.varymaxplot = 60
if gu.vartype == "Signal":
    gu.varymaxplot = 2000
    gu.varyminplot = -500


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
        valinit=mu.RC,
        step=1,
        color="C4",
    ),
    "z": dict(
        ax=[left2, low1, 0.025, dySlider],
        label="Drift\n[mm]",
        valmin=0,
        valmax=1000,
        valinit=30,
        step=1,
        color="C4",
    ),
    "x": dict(
        ax=[left1, low2, 0.025, dySlider],
        label="x [mm]",
        valmin=-geo.xwidth / 2,
        valmax=geo.xwidth / 2,
        valinit=0,
        step=0.01,
        color=gu.varcolor,
    ),
    "y": dict(
        ax=[left2, low2, 0.025, dySlider],
        label="y [mm]",
        valmin=-geo.ywidth / 2,
        valmax=geo.ywidth / 2,
        valinit=0,
        step=0.01,
        color=gu.varcolor,
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
gu.button.on_clicked(lambda event: [s.reset() for s in sliders.values()])


# ========================= Map Plot Setup =========================

axMap = gu.fig.add_axes((0.035, 0.82, 0.085, 0.12))
axMap.set_xlim(geo.xleft, geo.xleft + geo.nX * geo.xwidth)
axMap.set_ylim(geo.ylow, geo.ylow + geo.nY * geo.ywidth)
axMap.set_title("Drop position", fontsize=15)
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

# First values
x0 = sliders["x"].val + geo.xc
y0 = sliders["y"].val + geo.yc

dropmarker = mu.Dt * sliders["z"].val / 100 * 15 + 10  # marker size for drop point
(drop_point,) = axMap.plot(x0, y0, "o", color="red", markersize=dropmarker)

# ========================= Subplots Initialization =========================
for iX in range(geo.nX):
    for iY in range(geo.nY):
        ax = gu.axs[geo.nY - 1 - iY, iX]
        xL = geo.xleft + iX * geo.xwidth
        xR = xL + geo.xwidth
        yB = geo.ylow + iY * geo.ywidth
        yT = yB + geo.ywidth
        res = mu.Compute0D(
            mu.t,
            x0,
            y0,
            xL,
            xR,
            yB,
            yT,
            sliders["RC"].val,
            sliders["z"].val,
            gu.vartype,
        )
        if res is None:
            res = geo.np.zeros_like(mu.t)
        sig = gu.scalefactor * res[: mu.t.size]
        line = ax.plot(mu.t / gu.timescale, sig, lw=8 - geo.nY, color=gu.varcolor)
        txt = ax.text(
            0.96,
            0.93,
            f"{gu.dim} = {geo.np.max(sig):.0f} {gu.unit}\n T$_{{max}}$ = {mu.t[geo.np.argmax(sig)]/gu.timescale:.0f} {gu.timeunit}",
            ha="right",
            va="top",
            transform=ax.transAxes,
            fontsize=23 - 2 * geo.nY,
            bbox=dict(boxstyle="round", facecolor=gu.varcolor),
            color=gu.varforeground,
        )
        gu.lines.append(line)
        gu.texts.append(txt)

        if iY == 0:
            ax.set_xlabel(f"Time ({gu.timeunit})", fontsize=25 - geo.nY)
            ax.tick_params(axis="x", labelsize=20 - geo.nY)
        if iX == 0:
            ax.set_ylabel(f"{gu.ylabel}", fontsize=25 - geo.nX)
            ax.tick_params(axis="y", labelsize=20 - geo.nX)
        ax.grid()
        ax.set_xlim(gu.varxminplot, gu.varxmaxplot / gu.timescale)
        ax.set_ylim(gu.varyminplot, gu.varymaxplot)


# ========================= Slider Callback Update =========================
def update(val):
    """Update callback for sliders.

    Parameters
    ----------
    val : float
        Ignored by the function (required by matplotlib slider callback
        signature). The function reads current slider values directly and
        recomputes the per-cell traces.

    Side effects
    -----------
    - Updates the map marker for the drop point.
    - Recomputes time traces for every cell using Compute0D and sets the
      y-data of the plotted lines and summary text boxes.
    - Forces the figure canvas to redraw.
    """
    x0 = sliders["x"].val + geo.xc
    y0 = sliders["y"].val + geo.yc
    drop_point.set_data([x0], [y0])
    drop_point.set_markersize(mu.Dt * sliders["z"].val / 5 + 10)  # Update marker size

    for iX in range(geo.nX):
        for iY in range(geo.nY):
            idx = iX * geo.nY + iY
            xL = geo.xleft + iX * geo.xwidth
            xR = xL + geo.xwidth
            yB = geo.ylow + iY * geo.ywidth
            yT = yB + geo.ywidth
            res = mu.Compute0D(
                mu.t,
                x0,
                y0,
                xL,
                xR,
                yB,
                yT,
                sliders["RC"].val,
                sliders["z"].val,
                gu.vartype,
            )
            if res is None:
                res = geo.np.zeros_like(mu.t)
            sig = gu.scalefactor * res[: mu.t.size]
            gu.lines[idx][0].set_ydata(sig)
            gu.texts[idx].set_text(
                f"{gu.dim} = {geo.np.max(sig):.0f} {gu.unit}\n T$_{{max}}$ = {mu.t[geo.np.argmax(sig)]/gu.timescale:.0f} {gu.timeunit}"
            )
    gu.fig.canvas.draw_idle()


for s in sliders.values():
    s.on_changed(update)

# ========================= Show Plot =========================
plt.show()
