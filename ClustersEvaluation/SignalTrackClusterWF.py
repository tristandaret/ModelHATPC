"""Interactive cluster waveform display for a projected track.

This script uses the shared `ClustersUtils` interactive figure to show
the per-cluster and per-pad contributions from a track described by angle
and impact parameter. The plotted traces are normalized by the appropriate
cluster length where relevant.
"""

# Explicit imports from the shared Headers package
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.append("Headers/")
import ClustersUtils as clu
import GeometryUtils as geo
import ModelUtils as mu

# Local aliases for objects created by the interactive helper module
fig = clu.fig
axs = clu.axs
sliders = clu.sliders
map_line = clu.map_line
# Common model/geometry aliases used throughout the script
t = mu.t
nX = geo.nX
nY = geo.nY
xleft = geo.xleft
xwidth = geo.xwidth
ylow = geo.ylow
ywidth = geo.ywidth

# ========================= Lines =========================
labels = [
    "Central pad",
    "Vertical cluster",
    "Diagonal cluster",
    "Half-cross cluster",
    "True signal",
]
colors = ["C0", "C1", "C2", "C3", "Black"]

lines = [
    axs.plot(t, np.zeros_like(t), lw=3, label=label, color=color)[0]
    for label, color in zip(labels, colors)
]
lines[-1].set_linestyle("dotted")

# Visibility
visibility = [True, True, True, True, True]
for line, vis in zip(lines, visibility):
    line.set_visible(vis)

# Axes setup
axs.set_xlabel("Time [ns]")
axs.set_xlim(t[0], t[-1])
axs.set_ylabel("Signal [ADC counts/cm]")
axs.grid()

# Initial legend (only for visible lines)
axs.legend(handles=[line for line in lines if line.get_visible()])

# ========================= Update Function =========================
initial_ylim = None


def update(val):
    """Update the plotted traces when sliders change.

    Parameters
    ----------
    val : float
        Slider value passed by Matplotlib (ignored). The function reads
        current slider values from the global `sliders` dictionary and
        recomputes per-pad/cluster traces for the projected track.

    Side effects
    -----------
    - Updates the map line, per-line `y` data and legend.
    - Sets the initial y-limits on first call and triggers a redraw.
    """
    phi_val = max(sliders["phi"].val, 1e-6)
    m, q, phi_rad = geo.compute_line_params(phi_val, sliders["d"].val)
    Ldiag, Lvert, Lcros, Ltotal = geo.ClusterLengths(phi_rad, sliders["d"].val)

    v_x = np.linspace(
        geo.xleft - geo.xwidth, geo.xleft + (geo.nX + 1) * geo.xwidth, 10 * geo.nX
    )
    v_y = np.tan(phi_rad) * v_x + q
    map_line.set_data(v_x, v_y)

    buffers = {
        name: np.zeros_like(mu.t)
        for name in ["centralPad", "vertClus", "diagClus", "halfClus"]
    }

    for iX in range(geo.nX):
        for iY in range(geo.nY):
            sig = mu.Signal1D(
                mu.t,
                m,
                q,
                geo.xleft + iX * geo.xwidth,
                geo.xleft + (iX + 1) * geo.xwidth,
                geo.ylow + iY * geo.ywidth,
                geo.ylow + (iY + 1) * geo.ywidth,
                sliders["RC"].val,
                sliders["z"].val,
            )[: len(mu.t)]
            # Leading pad
            if (iX, iY) == (nX // 2, nY // 2):
                buffers["centralPad"] += sig
            # Vertical cluster
            if iX == nX // 2:
                buffers["vertClus"] += sig
            # Diagonal cluster
            if iX + iY == nX // 2 + nY // 2:
                buffers["diagClus"] += sig
            # Half-cross cluster
            if (iX >= nX // 2 and iY == nY // 2) or (iX == nX // 2 and iY >= nY // 2):
                buffers["halfClus"] += sig

    if Ldiag != 0:
        lines[0].set_ydata(buffers["centralPad"] / Ldiag * 10)
        lines[2].set_ydata(buffers["diagClus"] / Ldiag * 10)
    if Lvert != 0:
        lines[1].set_ydata(buffers["vertClus"] / Lvert * 10)
    if Lcros != 0:
        lines[3].set_ydata(buffers["halfClus"] / Lcros * 10)

    lines[-1].set_ydata(mu.lambdaG * 10 * mu.ETF(mu.t))

    axs.legend(handles=[line for line in lines if line.get_visible()], fontsize=20)

    global initial_ylim
    if initial_ylim is None:
        all_ydata = [line.get_ydata() for line in lines if line.get_visible()]
        ymin = min(np.min(y) for y in all_ydata if np.any(np.isfinite(y)))
        ymax = max(np.max(y) for y in all_ydata if np.any(np.isfinite(y)))
        if ymax > ymin:
            margin = 0.1 * (ymax - ymin)
            initial_ylim = (ymin - margin, ymax + margin)
        else:
            initial_ylim = (-1, 1)
        axs.set_ylim(initial_ylim)

    fig.canvas.draw_idle()


# Update trigger
for s in sliders.values():
    s.on_changed(update)

update(None)
plt.show()
