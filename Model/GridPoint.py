from sys import path

path.append("Headers/")
from ModelUtils import *
from GridUtils import *
from matplotlib.widgets import Slider

# Update axes range values
if vartype == "Charge":
    varymaxplot = 60
if vartype == "Signal":
    varymaxplot = 2000
    varyminplot = -500


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
        valmin=-xwidth / 2,
        valmax=xwidth / 2,
        valinit=0,
        step=0.01,
        color=varcolor,
    ),
    "y": dict(
        ax=[left2, low2, 0.025, dySlider],
        label="y [mm]",
        valmin=-ywidth / 2,
        valmax=ywidth / 2,
        valinit=0,
        step=0.01,
        color=varcolor,
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

axMap = fig.add_axes([0.035, 0.82, 0.085, 0.12])
axMap.set_xlim(xleft, xleft + nX * xwidth)
axMap.set_ylim(ylow, ylow + nY * ywidth)
axMap.set_title("Drop position", fontsize=15)
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

# First values
x0 = sliders["x"].val + xc
y0 = sliders["y"].val + yc

dropmarker = Dt * sliders["z"].val / 100 * 15 + 10  # marker size for drop point
(drop_point,) = axMap.plot(x0, y0, "o", color="red", markersize=dropmarker)

# ========================= Subplots Initialization =========================
for iX in range(nX):
    for iY in range(nY):
        ax = axs[nY - 1 - iY, iX]
        xL = xleft + iX * xwidth
        xR = xL + xwidth
        yB = ylow + iY * ywidth
        yT = yB + ywidth
        sig = (
            scalefactor
            * Compute0D(
                t, x0, y0, xL, xR, yB, yT, sliders["RC"].val, sliders["z"].val, vartype
            )[: t.size]
        )
        l = ax.plot(t / timescale, sig, lw=8 - nY, color=varcolor)
        txt = ax.text(
            0.96,
            0.93,
            f"{dim} = {np.max(sig):.0f} {unit}\n T$_{{max}}$ = {t[np.argmax(sig)]/timescale:.0f} {timeunit}",
            ha="right",
            va="top",
            transform=ax.transAxes,
            fontsize=23 - 2 * nY,
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
    x0 = sliders["x"].val + xc
    y0 = sliders["y"].val + yc
    drop_point.set_data([x0], [y0])
    drop_point.set_markersize(Dt * sliders["z"].val / 5 + 10)  # Update marker size

    for iX in range(nX):
        for iY in range(nY):
            idx = iX * nY + iY
            xL = xleft + iX * xwidth
            xR = xL + xwidth
            yB = ylow + iY * ywidth
            yT = yB + ywidth
            sig = (
                scalefactor
                * Compute0D(
                    t,
                    x0,
                    y0,
                    xL,
                    xR,
                    yB,
                    yT,
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
