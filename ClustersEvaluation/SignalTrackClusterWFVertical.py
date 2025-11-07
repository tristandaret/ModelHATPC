"""Vertical cluster waveform explorer.

Specialised visualisation that focuses on vertical neighbour contributions
for a projected track. Note: this view was written for `nY == 5` and
contains hard-coded index assumptions; use with care if `nY` differs.
"""

from sys import path

path.append("Headers/")
from ClustersUtils import *
from GeometryUtils import *

# DOES NOT WORK FOR nY != 5
# ========================= Lines =========================
labels = [
    "2nd top neighbor",
    "1st top neighbor",
    "Central pad",
    "1st bottom neighbor",
    "2nd bottom neighbor",
    "Vertical cluster",
    "True signal",
]
C0light = "#6baed6"
C0dark = "#14507c"

colors = [C0light, C0light, "C0", C0dark, C0dark, "C1", "Black"]

lines = [
    axs.plot(t, np.zeros_like(t), lw=3, label=label, color=color)[0]
    for label, color in zip(labels, colors)
]
lines[0].set_linestyle("dashdot")
lines[1].set_linestyle("dashed")
lines[2].set_linestyle("solid")
lines[3].set_linestyle("dashed")
lines[4].set_linestyle("dashdot")
lines[-1].set_linestyle("dotted")

axs.set_xlabel("Time [ns]")
axs.set_ylabel("Signal [ADC counts/cm]")
axs.grid()
axs.legend()


# ========================= Update Function =========================
initial_ylim = None


def update(val):
    phi_val = max(sliders["phi"].val, 1e-6)
    m, q, phi_rad = compute_line_params(phi_val, sliders["d"].val)
    r_diag, r_vert, r_cros, L = ClusterLengths(phi_rad, sliders["d"].val)

    v_x = np.linspace(xleft - xwidth, xleft + (nX + 1) * xwidth, 10 * nX)
    v_y = np.tan(phi_rad) * v_x + q
    map_line.set_data(v_x, v_y)

    buffers = {
        name: np.zeros_like(t)
        for name in [
            "topNeigh2",
            "topNeigh1",
            "LeadingPad",
            "botNeigh1",
            "botNeigh2",
            "clusVert",
        ]
    }

    for iX, iY in {(2, i) for i in range(nY)}:
        print("iY:", iY)
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
        )
        sig = sig[: len(t)]  # Prevent shape mismatch
        buffers["clusVert"] += sig
        if iY == 4:
            buffers["topNeigh2"] += sig
        if iY == 3:
            buffers["topNeigh1"] += sig
        if iY == 2:
            buffers["LeadingPad"] += sig
        if iY == 1:
            buffers["botNeigh1"] += sig
        if iY == 0:
            buffers["botNeigh2"] += sig

    if r_vert != 0:
        lines[0].set_ydata(buffers["topNeigh2"] / r_vert * 10)
        lines[1].set_ydata(buffers["topNeigh1"] / r_vert * 10)
        lines[2].set_ydata(buffers["LeadingPad"] / r_vert * 10)
        lines[3].set_ydata(buffers["botNeigh1"] / r_vert * 10)
        lines[4].set_ydata(buffers["botNeigh2"] / r_vert * 10)
        lines[5].set_ydata(buffers["clusVert"] / r_vert * 10)

    lines[6].set_ydata(lambdaG * 10 * ETF(t))
    axs.legend(fontsize=20)

    # Update Y-axis limits
    global initial_ylim
    # Set Y-axis limits only once, on the first update
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


for s in sliders.values():
    s.on_changed(update)

update(None)  # Initial update to set the plot
plt.show()
