"""Overlay of point-deposit waveforms from neighbouring pads.

This script computes the shaped signal for a chosen point-like deposit and
plots the contribution from each pad in the `nX x nY` sub-map. It is useful
to visualise how charge spreads to neighbors and to validate the summed
response against a no-diffusion reference.
"""

import numpy as np
import matplotlib.pyplot as plt
from sys import path

path.append("Headers/")
from ModelUtils import *

plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": 35})


# Computations punctual WFs --------------------------------------------------------------------------------------------------
xdrop = xleft + xc + nX // 2 * xwidth + 2
ydrop = ylow + yc + nY // 2 * ywidth
z = 0

# Create an array to store the results
signal_results = np.zeros((nX, nY, len(t)))

# Create the object to store the sum of all signals
SumSignal = np.zeros(len(t))
SignalRCinf = (
    Signal0D(t, xdrop, ydrop, xmin, xmax, ymin, ymax, 1e8, z)[: len(t)] * Ne / Ne55Fe
)

for iX in range(nX):
    for iY in range(nY):
        xmin = xleft + iX * xwidth
        xmax = xmin + xwidth
        ymin = ylow + iY * ywidth
        ymax = ymin + ywidth
        convo0D = (
            Signal0D(t, xdrop, ydrop, xmin, xmax, ymin, ymax, RC, z)[: len(t)]
            * Ne
            / Ne55Fe
        )
        signal_results[iX, iY, :] = convo0D
        SumSignal += convo0D


# Draw plot --------------------------------------------------------------------------------------------------------------------
fig, ax1 = plt.subplots(figsize=(12, 9))
plt.subplots_adjust(left=0.15, right=0.96, bottom=0.13, top=0.98)

ax1.plot(t, SignalRCinf, color="black", lw=3, label="No diffusion case")
# ax1.plot(t, SumSignal, color="C3", lw=3, ls="--", label=r"Sum of 3 $\times$ 3 signals")

# # Central pad
# ax1.plot(
#     t,
#     signal_results[nX // 2, nY // 2, :],
#     color="C0",
#     lw=2,
#     ls="--",
#     label=f"Central pad",
# )

# for iX in range(nX):
#     for iY in range(nY):
#         #  Vertical neighbors
#         if iX == 1 and iY != 1:
#             label = "Vertical neighbors" if (iX == 1 and iY == 0) else None
#             ax1.plot(
#                 t,
#                 signal_results[iX, iY, :],
#                 color="C1",
#                 lw=2,
#                 ls="--",
#                 label=label,
#             )
#         # Lateral neighbors
#         if iY == 1 and iX != 1:
#             label = "Lateral neighbors" if (iX == 0 and iY == 1) else None
#             ax1.plot(
#                 t,
#                 signal_results[iX, iY, :],
#                 color="C2",
#                 lw=2,
#                 ls="--",
#                 label=label,
#             )
#         # Corner pads
#         if (iX == 0 or iX == 2) and (iY == 0 or iY == 2):
#             label = "Corner neighbors" if (iX == 0 and iY == 0) else None
#             ax1.plot(
#                 t,
#                 signal_results[iX, iY, :],
#                 color="C4",
#                 lw=2,
#                 ls="--",
#                 label=label,
#             )

ax1.set_ylabel("Signal [ADC counts]")
ax1.set_xlim(0, t[-1])
ax1.set_xlabel("Time (ns)", labelpad=10)
ax1.tick_params(axis="x")
ax1.tick_params(axis="y")
ax1.grid()
ax1.legend()

# plt.savefig(f"Illustrations/SignalPointOverlayWF3x3X2mm.pdf", bbox_inches="tight")
# plt.savefig(f"Illustrations/SignalPointOverlayWF3x3NoSumX2mm.pdf", bbox_inches="tight")
plt.savefig(
    f"Illustrations/SignalPointOverlayWF3x3NoDiffusionOnlyX2mm.pdf", bbox_inches="tight"
)
plt.show()
