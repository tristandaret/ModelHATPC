"""Interactive cluster waveform display for a projected track.

This script uses the shared `ClustersUtils` interactive figure to show
the per-cluster and per-pad contributions from a track described by angle
and impact parameter. The plotted traces are normalized by the appropriate
cluster length where relevant.
"""

from sys import path

path.append("Headers/")
from ClustersUtils import *

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
    m, q, phi_rad = compute_line_params(phi_val, sliders["d"].val)
    Ldiag, Lvert, Lcros, Ltotal = ClusterLengths(phi_rad, sliders["d"].val)

    v_x = np.linspace(xleft - xwidth, xleft + (nX + 1) * xwidth, 10 * nX)
    v_y = np.tan(phi_rad) * v_x + q
    map_line.set_data(v_x, v_y)

    buffers = {
        name: np.zeros_like(t)
        for name in ["centralPad", "vertClus", "diagClus", "halfClus"]
    }

    for iX in range(nX):
        for iY in range(nY):
            sig = Signal1D(
                t,
                m,
                q,
                xleft + iX * xwidth,
                xleft + (iX + 1) * xwidth,
                ylow + iY * ywidth,
                ylow + (iY + 1) * ywidth,
                sliders["RC"].val,
                sliders["z"].val,
            )[: len(t)]
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

    lines[-1].set_ydata(lambdaG * 10 * ETF(t))

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
