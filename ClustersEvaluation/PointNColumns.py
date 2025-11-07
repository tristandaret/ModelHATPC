import numpy as np
import matplotlib.pyplot as plt
from sys import path

path.append("Headers/")
from ModelUtils import *

plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": 35})


# Computations of 0D & 1D ETFs --------------------------------------------------------------------------------------------------
# 1D reference
Convoline = Signal1D(t, 1e-6, yc, xmin, xmax, ymin, ymax, RC, z)[: len(t)]

# 0D sum of point deposit
ncases = nX // 2
v_Conv0D = [np.zeros(t.size) for _ in range(ncases + 1)]
for iX in range(nX):
    xdrop = xleft + xc + iX * xwidth
    convo0D = (
        Signal0D(t, xdrop, yc, xmin, xmax, ymin, ymax, RC, z)[: len(t)] * Ne / Ne55Fe
    )
    for i in range(len(v_Conv0D)):
        if nX // 2 - i <= iX <= nX // 2 + i:
            v_Conv0D[i] += convo0D


# Draw plots --------------------------------------------------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(
    2, 1, figsize=(12, 9), sharex=True, gridspec_kw={"height_ratios": [2, 1]}
)
plt.subplots_adjust(left=0.08, right=0.96, bottom=0.13, top=0.98, hspace=0.08)

for i in range(ncases + 1):
    v_Conv0D[i] *= xwidth
    diff = (v_Conv0D[i] - Convoline) / np.max(Convoline) * 100  # Percentage difference
    c = (i / ncases, 0, 1 - i / ncases)
    if i == 0:
        for j in range(len(v_Conv0D[i])):
            print(
                f"{j*6} {v_Conv0D[i][j]} - {Convoline[j]} = {v_Conv0D[i][j] - Convoline[j]}"
            )
        ax1.plot(t, v_Conv0D[i], color=c, lw=3, label=("1 column"))
        ax2.plot(t, diff, color=c, lw=3, label=("1 column"))
    else:
        ax1.plot(t, v_Conv0D[i], color=c, lw=3, label=(f"{2*i+1:.0f} columns"))
        ax2.plot(t, diff, color=c, lw=3, label=(f"{2*i+1:.0f} columns"))


# Main plot
ax1.plot(t, Convoline, color="black", lw=3, ls="--", label="Linear deposit")
ax1.set_ylabel("Signal [ADC counts/cm]")
ax1.set_xlim(0, t[-1])
ax1.tick_params(axis="x")
ax1.tick_params(axis="y")
ax1.grid()
ax1.legend()

# Difference plot
ax2.axhline(0, color="black", lw=3, ls="--")
ax2.set_xlabel("Time (ns)", labelpad=10)
ax2.set_ylabel("Difference [\%]", labelpad=20)
ax2.set_ylim(-15, 5)
ax2.grid()

plt.savefig(f"Illustrations/Point{nX}Columns.pdf", bbox_inches="tight")
# plt.show()
