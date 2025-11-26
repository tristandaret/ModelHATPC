"""Visualize precomputed LUT-like weight arrays and length maps.

This utility builds the length map and the corresponding weight array used
to correct ADC amplitudes into dE/dx. It relies on functions in
`Headers/ModelUtils.py` and `Headers/GeometryUtils.py` and produces a heatmap
showing the weight as a function of track angle and impact parameter.
"""

import time
from typing import List

import numpy as np
from matplotlib import pyplot as plt

from Headers import GeometryUtils as geo
from Headers import ModelUtils as mu

# Computations
ETF_t = mu.lambdaG * mu.ETF(mu.t)
nphi = 250
nd = 250
RC = 120
z = 250

h2_L = np.full((nd, nphi), np.nan)
v_d = np.linspace(0, geo.diag / 2, nd)
v_phi = np.linspace(1e-6, 90 - 1e-6, nphi)
weight_array = np.full((nd, nphi), np.nan)


# Make Length map
for iphi in range(nphi):
    phi_rad = v_phi[iphi] / 180 * np.pi
    for id in range(nd):
        d = v_d[id]

        # Determine the length of the track across the central pad
        x: List[float] = []
        y: List[float] = []

        y_xmin = geo.Y(phi_rad, d, geo.xmin)
        y_xmax = geo.Y(phi_rad, d, geo.xmax)
        x_ymin = geo.X(phi_rad, d, geo.ymin)
        x_ymax = geo.X(phi_rad, d, geo.ymax)

        if geo.ymin <= y_xmin < geo.ymax:
            x.append(geo.xmin)
            y.append(y_xmin)

        if geo.ymin <= y_xmax < geo.ymax:
            x.append(geo.xmax)
            y.append(y_xmax)

        if geo.xmin <= x_ymin < geo.xmax:
            x.append(x_ymin)
            y.append(geo.ymin)

        if geo.xmin <= x_ymax < geo.xmax:
            x.append(x_ymax)
            y.append(geo.ymax)

        L = 0
        if len(x) == 2:
            L = np.sqrt((y[1] - y[0]) ** 2 + (x[1] - x[0]) ** 2)

        h2_L[id, iphi] = L


# Make LUT
start_time = time.time()
for iphi in range(nphi):
    phi_rad = v_phi[iphi] / 180 * np.pi
    print(f"phi = {v_phi[iphi]:.1f}°")

    for id in range(nd):
        # print(f"d = {d:.2f} mm")
        m = np.tan(phi_rad)
        q = (np.cos(phi_rad) * geo.yc - np.sin(phi_rad) * geo.xc + v_d[id]) / np.cos(
            phi_rad
        )  # intercept

        ETFr = h2_L[id, iphi] * np.max(ETF_t)
        ADC = np.max(
            mu.Signal1D(mu.t, m, q, geo.xmin, geo.xmax, geo.ymin, geo.ymax, RC, z)[
                : len(mu.t)
            ]
        )
        if h2_L[id, iphi] != 0:
            weight_array[id, iphi] = ETFr / ADC

print(f"LUT done in {time.time()-start_time:.1f} seconds")


# Draw the heatmap
plt.figure(figsize=(16, 9))
plt.imshow(
    weight_array,
    cmap="viridis",
    origin="lower",
    aspect="auto",
    extent=(float(v_phi[0]), float(v_phi[-1]), float(v_d[0]), float(v_d[-1])),
)
cbar = plt.colorbar(pad=0.01, fraction=0.08)
cbar.ax.tick_params(labelsize=20)

plt.xlabel(r" track angle $\varphi$ (°)", fontsize=20)
plt.xticks(fontsize=20)
plt.ylabel("impact parameter d (mm)", fontsize=20)
plt.yticks(fontsize=20)

plt.tight_layout()
plt.savefig(
    f"LUT/LUT_Dt{mu.Dt*np.power(10, 7/2):.0f}_PT{mu.PT:d}_RC{RC:d}_z{z:d}_nphi{nphi:d}_nd{nd:d}.pdf"
)
plt.savefig(
    f"LUT/LUT_Dt{mu.Dt*np.power(10, 7/2):.0f}_PT{mu.PT:d}_RC{RC:d}_z{z:d}_nphi{nphi:d}_nd{nd:d}.png"
)

plt.show()
