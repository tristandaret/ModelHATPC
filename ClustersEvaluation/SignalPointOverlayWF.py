"""Overlay of point-deposit waveforms from neighbouring pads.

This script computes the shaped signal for a chosen point-like deposit and
plots the contribution from each pad in the `nX x nY` sub-map. It is useful
to visualise how charge spreads to neighbors and to validate the summed
response against a no-diffusion reference.
"""

# Explicit Headers imports
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.append("Headers/")
import GeometryUtils as geo
import ModelUtils as mu

plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": 35})


# Computations punctual WFs --------------------------------------------------------------------------------------------------
xdrop = geo.xleft + geo.xc + geo.nX // 2 * geo.xwidth
ydrop = geo.ylow + geo.yc + geo.nY // 2 * geo.ywidth
z = 0

# Create an array to store the results
signal_results = np.zeros((geo.nX, geo.nY, len(mu.t)))

# Create the object to store the sum of all signals
SumSignal = np.zeros(len(mu.t))
SignalRCinf = (
    mu.Signal0D(
        mu.t,
        xdrop,
        ydrop,
        geo.xleft,
        geo.xleft + geo.nX * geo.xwidth,
        geo.ylow,
        geo.ylow + geo.nY * geo.ywidth,
        1e8,
        z,
    )[: len(mu.t)]
    * mu.Ne
    / mu.Ne55Fe
)

for iX in range(geo.nX):
    for iY in range(geo.nY):
        xmin = geo.xleft + iX * geo.xwidth
        xmax = xmin + geo.xwidth
        ymin = geo.ylow + iY * geo.ywidth
        ymax = ymin + geo.ywidth
        xmin = geo.xleft + iX * geo.xwidth
        xmax = xmin + geo.xwidth
        ymin = geo.ylow + iY * geo.ywidth
        ymax = ymin + geo.ywidth
        xmin = geo.xleft + iX * geo.xwidth
        xmax = xmin + geo.xwidth
        ymin = geo.ylow + iY * geo.ywidth
        ymax = ymin + geo.ywidth
        convo0D = (
            mu.Signal0D(mu.t, xdrop, ydrop, xmin, xmax, ymin, ymax, mu.RC, z)[
                : len(mu.t)
            ]
            * mu.Ne
            / mu.Ne55Fe
        )
        signal_results[iX, iY, :] = convo0D
        SumSignal += convo0D


# Draw plot --------------------------------------------------------------------------------------------------------------------
fig, ax1 = plt.subplots(figsize=(12, 9))
plt.subplots_adjust(left=0.15, right=0.96, bottom=0.13, top=0.98)

ax1.plot(mu.t, SignalRCinf, color="black", lw=3, label="No diffusion case")

ax1.set_ylabel("Signal [ADC counts]")
ax1.set_xlim(0, mu.t[-1])
ax1.set_xlabel("Time (ns)", labelpad=10)
ax1.tick_params(axis="x")
ax1.tick_params(axis="y")
ax1.grid()
ax1.legend()

# plt.savefig(f"Illustrations/SignalPointOverlayWF3x3X2mm.pdf", bbox_inches="tight")
# plt.savefig(f"Illustrations/SignalPointOverlayWF3x3NoSumX2mm.pdf", bbox_inches="tight")
plt.savefig(
    "Illustrations/SignalPointOverlayWF3x3NoDiffusionOnlyX2mm.pdf", bbox_inches="tight"
)
plt.show()
