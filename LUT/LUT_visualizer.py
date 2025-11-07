from sys import path

path.append("Headers/")
from ModelUtils import *
import time
from matplotlib import pyplot as plt


# Computations
ETF = lambdaG * ETF(t)
nphi = 250
nd = 250
RC = 120
z = 250

h2_L = np.full((nd, nphi), np.nan)
v_d = np.linspace(0, diag / 2, nd)
v_phi = np.linspace(1e-6, 90 - 1e-6, nphi)
weight_array = np.full((nd, nphi), np.nan)


# Make Length map
for iphi in range(nphi):
    phi_rad = v_phi[iphi] / 180 * np.pi
    for id in range(nd):
        d = v_d[id]

        # Determine the length of the track across the central pad
        x = []
        y = []

        y_xmin = Y(phi_rad, d, xmin)
        y_xmax = Y(phi_rad, d, xmax)
        x_ymin = X(phi_rad, d, ymin)
        x_ymax = X(phi_rad, d, ymax)

        if ymin <= y_xmin < ymax:
            x.append(xmin)
            y.append(y_xmin)

        if ymin <= y_xmax < ymax:
            x.append(xmax)
            y.append(y_xmax)

        if xmin <= x_ymin < xmax:
            x.append(x_ymin)
            y.append(ymin)

        if xmin <= x_ymax < xmax:
            x.append(x_ymax)
            y.append(ymax)

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
        q = (np.cos(phi_rad) * yc - np.sin(phi_rad) * xc + v_d[id]) / np.cos(
            phi_rad
        )  # intercept

        ETFr = h2_L[id, iphi] * np.max(ETF)
        ADC = np.max(Signal1D(t, m, q, xmin, xmax, ymin, ymax, RC, z)[: len(t)])
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
    extent=[v_phi[0], v_phi[-1], v_d[0], v_d[-1]],
)
cbar = plt.colorbar(pad=0.01, fraction=0.08)
cbar.ax.tick_params(labelsize=20)

plt.xlabel(r" track angle $\varphi$ (°)", fontsize=20)
plt.xticks(fontsize=20)
plt.ylabel("impact parameter d (mm)", fontsize=20)
plt.yticks(fontsize=20)

plt.tight_layout()
plt.savefig(
    f"LUT/LUT_Dt{Dt*np.power(10,7/2):.0f}_PT{PT:d}_RC{RC:d}_z{z:d}_nphi{nphi:d}_nd{nd:d}.pdf"
)
plt.savefig(
    f"LUT/LUT_Dt{Dt*np.power(10,7/2):.0f}_PT{PT:d}_RC{RC:d}_z{z:d}_nphi{nphi:d}_nd{nd:d}.png"
)

plt.show()
