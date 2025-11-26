"""Compare waveforms from point-like deposits summed over multiple columns.

This script computes the shaped signal for a line deposit (1D reference)
and compares it to the sum of individual point-like (0D) deposits placed on
different columns. The result is saved as a PDF in `Illustrations/`.

Typical usage: run interactively to reproduce figures used in evaluation of
the lineic vs punctual deposition approximations.
"""

from sys import path

import GeometryUtils as geo
import matplotlib.pyplot as plt
import ModelUtils as mu
import numpy as np

path.append("Headers/")

plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": 35})


# Computations of 0D & 1D ETFs --------------------------------------------------------------------------------------------------
# 1D reference
Convoline = mu.Signal1D(
    mu.t, 1e-6, geo.yc, geo.xmin, geo.xmax, geo.ymin, geo.ymax, mu.RC, mu.z
)[: len(mu.t)]

# 0D sum of point deposit
ncases = geo.nX // 2
v_Conv0D = [np.zeros(mu.t.size) for _ in range(ncases + 1)]
for iX in range(geo.nX):
    xdrop = geo.xleft + geo.xc + iX * geo.xwidth
    convo0D = (
        mu.Signal0D(
            mu.t, xdrop, geo.yc, geo.xmin, geo.xmax, geo.ymin, geo.ymax, mu.RC, mu.z
        )[: len(mu.t)]
        * mu.Ne
        / mu.Ne55Fe
    )
    for i in range(len(v_Conv0D)):
        if geo.nX // 2 - i <= iX <= geo.nX // 2 + i:
            v_Conv0D[i] += convo0D


# Draw plots --------------------------------------------------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(
    2, 1, figsize=(12, 9), sharex=True, gridspec_kw={"height_ratios": [2, 1]}
)
plt.subplots_adjust(left=0.08, right=0.96, bottom=0.13, top=0.98, hspace=0.08)

for i in range(ncases + 1):
    v_Conv0D[i] *= geo.xwidth
    diff = (v_Conv0D[i] - Convoline) / np.max(Convoline) * 100  # Percentage difference
    c = (i / ncases, 0, 1 - i / ncases)
    if i == 0:
        for j in range(len(v_Conv0D[i])):
            print(
                f"{j*6} {v_Conv0D[i][j]} - {Convoline[j]} = {v_Conv0D[i][j] - Convoline[j]}"
            )
        ax1.plot(mu.t, v_Conv0D[i], color=c, lw=3, label=("1 column"))
        ax2.plot(mu.t, diff, color=c, lw=3, label=("1 column"))
    else:
        ax1.plot(mu.t, v_Conv0D[i], color=c, lw=3, label=(f"{2*i+1:.0f} columns"))
        ax2.plot(mu.t, diff, color=c, lw=3, label=(f"{2*i+1:.0f} columns"))

# Main plot
ax1.plot(mu.t, Convoline, color="black", lw=3, ls="--", label="Linear deposit")
ax1.set_ylabel("Signal [ADC counts/cm]")
ax1.set_xlim(0, mu.t[-1])
ax1.tick_params(axis="x")
ax1.tick_params(axis="y")
ax1.grid()
ax1.legend()

# Difference plot
ax2.axhline(0, color="black", lw=3, ls="--")
ax2.set_xlabel("Time (ns)", labelpad=10)
ax2.set_ylabel(r"Difference [\%]", labelpad=20)
ax2.set_ylim(-15, 5)
ax2.grid()

plt.savefig(f"Illustrations/Point{geo.nX}Columns.pdf", bbox_inches="tight")
# plt.show()
