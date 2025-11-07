from sys import path

path.append("Headers/")
from ModelUtils import *
from GridUtils import *
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
        valinit=RC,
        step=1,
        color="C4"
    ),
    "z": dict(
        ax=[left2, low1, 0.025, dySlider],
        label="Drift\n[mm]",
        valmin=0,
        valmax=1000,
        valinit=z,
        step=1,
        color="C4"
    ),
    "phi": dict(
        ax=[left1, low2, 0.025, dySlider],
        label="angle\n[°]",
        valmin=1e-6,
        valmax=90 - 1e-6,
        valinit=42,
        step=0.1,
        color=varcolor
    ),
    "d": dict(
        ax=[left2, low2, 0.025, dySlider],
        label="impact\n[mm]",
        valmin=-diag / 2,
        valmax=diag / 2,
        valinit=0,
        step=0.1,
        color=varcolor
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
button.on_clicked(lambda event: [s.reset() for s in sliders.values()])


# ========================= Map Plot Setup =========================
m, q, phi_rad = compute_line_params(sliders["phi"].val, sliders["d"].val)

axMap = fig.add_axes([0.035, 0.82, 0.085, 0.12])
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


# ========================= Subplots Initialization =========================
for iX in range(nX):
    for iY in range(nY):
        ax = axs[nY - 1 - iY, iX]
        x0 = xleft + iX * xwidth
        x1 = x0 + xwidth
        y0 = ylow + iY * ywidth
        y1 = y0 + ywidth
        sig = scalefactor * Compute1D(t, m, q, x0, x1, y0, y1, RC, z, vartype)[: t.size]
        l = ax.plot(t / timescale, sig, lw=8 - nY, color=varcolor)
        txt = ax.text(
            0.96,
            0.93,
            f"{dim} = {np.max(sig):.0f} {unit}\n T$_{{max}}$ = {t[np.argmax(sig)]/timescale:.0f} {timeunit}",
            ha="right",
            va="top",
            transform=ax.transAxes,
            fontsize=25 - 2*nY,
            bbox=dict(boxstyle="round", facecolor=varcolor),
            color=varforeground,
        )
        lines.append(l)
        texts.append(txt)

        if iY == 0:
            ax.set_xlabel(f"Time ({timeunit})", fontsize=25 - nY)
            ax.tick_params(axis="x", labelsize=20 - nY)
        if iX == 0:
            ax.set_ylabel(f"{ylabel}", fontsize=25 - nX)
            ax.tick_params(axis="y", labelsize=20 - nX)
        ax.grid()
        ax.set_xlim(varxminplot, varxmaxplot / timescale)
        ax.set_ylim(varyminplot, varymaxplot)


# ========================= Slider Callback Update =========================
def update(val):
    phi = sliders["phi"].val
    if abs(phi) < 1e-6:
        phi = 1e-6 * np.sign(phi or 1)

    m, q, phi_rad = compute_line_params(phi, sliders["d"].val)

    v_x = np.linspace(xleft - xwidth, xleft + (nX + 1) * xwidth, 10 * nX)
    v_y = np.tan(phi_rad) * v_x + q
    map_line.set_data(v_x, v_y)

    for iX in range(nX):
        for iY in range(nY):
            idx = iX * nY + iY
            sig = (
                scalefactor
                * Compute1D(
                    t,
                    m,
                    q,
                    xleft + iX * xwidth,
                    xleft + (iX + 1) * xwidth,
                    ylow + iY * ywidth,
                    ylow + (iY + 1) * ywidth,
                    sliders["RC"].val,
                    sliders["z"].val,
                    vartype,
                )[: t.size]
            )
            lines[idx][0].set_ydata(sig)
            texts[idx].set_text(
                f"{dim} = {np.max(sig):.0f} {unit}\n T$_{{max}}$ = {t[np.argmax(sig)]/timescale:.0f} {timeunit}"
            )
    fig.canvas.draw_idle()


for s in sliders.values():
    s.on_changed(update)

# ========================= Show Plot =========================
plt.show()
