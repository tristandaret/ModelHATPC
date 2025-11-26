"""Estimate cluster amplitude ratios and produce diagnostic heatmaps.

This script scans track angles and impact parameters across a small
sub-map and computes the ratios between the true deposited amplitude and
different cluster estimators (diagonal, half-cross, vertical, leading pad).
It stores several diagnostic heatmaps and histograms into a multi-page PDF
in `Illustrations/`.

Note: the script can be computationally intensive depending on the grid
resolution (`nsteps`) and the number of tracks considered.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

from Headers import GeometryUtils as geo
from Headers import ModelUtils as mu

# Plotting style with LaTeX
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": 35})


# Computations
nsteps = 250
# use constants from ModelUtils when appropriate
RC = mu.RC  # Readout chip radius (from model defaults)
z = mu.z  # Distance from the readout chip to the sensor
ETF_peak = mu.lambdaG * mu.max_ETF
nd = nsteps
nphi = nsteps
v_d = np.linspace(0, geo.diag / 2, nd)
v_phi = np.linspace(1e-6, 90 - 1e-6, nphi)
print("v_phi elements:", v_phi)
arr_lead = np.full((nd, nphi), np.nan)
arr_diag = np.full((nd, nphi), np.nan)
arr_cros = np.full((nd, nphi), np.nan)
arr_vert = np.full((nd, nphi), np.nan)
list_lead = []
list_diag = []
list_cros = []
list_vert = []
list_diag_cut = []
list_cros_cut = []
list_vert_cut = []

phi_index = 0
for phi in v_phi:
    print(f"phi = {phi:.1f}°")
    d_index = 0
    for d in v_d:
        m, q, phi_rad = geo.compute_line_params(phi, d)

        # Determine the cluster lengths
        r_diag, r_vert, r_cros, L = geo.ClusterLengths(phi_rad, d)

        ETF_lead = r_diag * ETF_peak
        ETF_diag = r_diag * ETF_peak
        ETF_cros = r_cros * ETF_peak
        ETF_vert = r_vert * ETF_peak
        Lead_only = np.zeros(mu.t.size)
        Clus_Diag = np.zeros(mu.t.size)
        Clus_HC = np.zeros(mu.t.size)
        Clus_vert = np.zeros(mu.t.size)
        for iX in range(geo.nX):
            for iY in range(geo.nY):
                ADC = mu.Signal1D(
                    mu.t,
                    m,
                    q,
                    geo.xleft + iX * geo.xwidth,
                    geo.xleft + (iX + 1) * geo.xwidth,
                    geo.ylow + iY * geo.ywidth,
                    geo.ylow + (iY + 1) * geo.ywidth,
                    RC,
                    z,
                )[: len(mu.t)]
                # Leading pad only
                if iX == geo.nX // 2 and iY == geo.nY // 2:
                    Lead_only += ADC
                # Cluster distributions
                # Diagonal distribution
                if iX + iY == geo.nX // 2 + geo.nY // 2:
                    Clus_Diag += ADC
                # Half-cross distribution
                if (iX >= geo.nX // 2 and iY == geo.nY // 2) or (
                    iX == geo.nX // 2 and iY >= geo.nY // 2
                ):
                    Clus_HC += ADC
                # Vertical distribution
                if iX == geo.nX // 2:
                    Clus_vert += ADC

        # Cluster distributions (r_diag cuts on length in leading pad as a common benchmark)
        if r_diag > 0:
            arr_lead[d_index][phi_index] = ETF_diag / np.max(Lead_only)
            list_lead.append(ETF_diag / np.max(Lead_only))
            arr_diag[d_index][phi_index] = ETF_diag / np.max(Clus_Diag)
            list_diag.append(ETF_diag / np.max(Clus_Diag))
            arr_cros[d_index][phi_index] = ETF_cros / np.max(Clus_HC)
            list_cros.append(ETF_cros / np.max(Clus_HC))
            arr_vert[d_index][phi_index] = ETF_vert / np.max(Clus_vert)
            list_vert.append(ETF_vert / np.max(Clus_vert))
        if r_diag > 2:
            if 30 <= phi <= 60:
                list_diag_cut.append(ETF_diag / np.max(Clus_Diag))
                list_cros_cut.append(ETF_cros / np.max(Clus_HC))
            if phi <= 45:
                list_vert_cut.append(ETF_vert / np.max(Clus_vert))

        d_index += 1
    phi_index += 1

list_global = list_cros + list_diag + list_vert


# Draw the plots
with PdfPages(
    f"Illustrations/Ratio_Diag_HC_nphi{nphi:d}_nd{nd:d}_RC{RC:d}_z{z:d}_PT{mu.PT:d}_Dt{mu.Dt*np.power(10, 7/2):.0f}_2mm_30phi60_geo.nX{geo.nX:d}_geo.nY{geo.nY:d}.pdf"
) as pdf:

    # Diagonal heatmap #
    plt.figure(figsize=(15, 10))
    plt.imshow(
        arr_diag,
        cmap="viridis",
        origin="lower",
        aspect="auto",
        extent=(float(v_phi[0]), float(v_phi[-1]), float(v_d[0]), float(v_d[-1])),
    )
    cbar = plt.colorbar(pad=0.01, fraction=0.08)
    plt.clim(0, np.max(list_global))
    plt.xlabel(r" track angle [°]")
    plt.ylabel("impact parameter [mm]", labelpad=15)
    plt.tight_layout()
    pdf.savefig(bbox_inches="tight")
    plt.clf()

    # Distribution of cluster accuracy #
    plt.grid()
    ax = plt.gca()
    plt.xlabel(r"$A_{true}/A_{diagonal}$")
    # Transparent all data
    plt.hist(
        list_diag,
        bins=50,
        range=(0, np.max(list_global)),
        color="blue",
        histtype="stepfilled",
        alpha=0.3,
    )
    # Outline all data
    plt.hist(
        list_diag,
        bins=50,
        range=(0, np.max(list_global)),
        color="blue",
        histtype="step",
        linewidth=3,
    )
    # Transparent cut data
    plt.hist(
        list_diag_cut,
        bins=50,
        range=(0, np.max(list_global)),
        color="red",
        histtype="stepfilled",
        alpha=0.3,
    )
    # Outline cut data
    plt.hist(
        list_diag_cut,
        bins=50,
        range=(0, np.max(list_global)),
        color="red",
        histtype="step",
        linewidth=3,
    )
    # mean and std
    plt.text(
        0.95,
        0.95,
        r"$\mu$ = "
        f"{np.mean(list_diag):.2f}\n"
        r"$\sigma$ = "
        f"{np.std(list_diag):.2f}",
        ha="right",
        va="top",
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", facecolor="blue"),
        color="white",
    )
    plt.text(
        0.95,
        0.75,
        r"$\mu$ = "
        f"{np.mean(list_diag_cut):.2f}\n"
        r"$\sigma$ = "
        f"{np.std(list_diag_cut):.2f}",
        ha="right",
        va="top",
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", facecolor="red"),
        color="white",
    )
    # Cuts
    plt.text(
        0.05,
        0.95,
        r"$0^\circ < \varphi < 90^\circ$" "\n L $>$ 0 mm",
        ha="left",
        va="top",
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", facecolor="blue"),
        color="white",
    )
    plt.text(
        0.05,
        0.75,
        r"$30^\circ \leq \varphi \leq 60^\circ$" "\n L $\\geq$ 2 mm",
        ha="left",
        va="top",
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", facecolor="red"),
        color="white",
    )
    plt.tight_layout()
    pdf.savefig(bbox_inches="tight")
    plt.clf()

    # Half-cross heatmap #
    plt.imshow(
        arr_cros,
        cmap="viridis",
        origin="lower",
        aspect="auto",
        extent=(float(v_phi[0]), float(v_phi[-1]), float(v_d[0]), float(v_d[-1])),
    )
    cbar = plt.colorbar(pad=0.01, fraction=0.08)
    plt.clim(0, np.max(list_global))
    plt.xlabel(r" track angle [°]")
    plt.ylabel("impact parameter [mm]", labelpad=15)
    plt.tight_layout()
    pdf.savefig(bbox_inches="tight")
    plt.clf()

    # Distribution of cluster accuracy #
    plt.grid()
    ax = plt.gca()
    plt.xlabel(r"$A_{true}/A_{half-cross}$")
    # Transparent all data
    plt.hist(
        list_cros,
        bins=50,
        range=(0, np.max(list_global)),
        color="blue",
        histtype="stepfilled",
        alpha=0.3,
    )
    # Outline all data
    plt.hist(
        list_cros,
        bins=50,
        range=(0, np.max(list_global)),
        color="blue",
        histtype="step",
        linewidth=3,
    )
    # Transparent cut data
    plt.hist(
        list_cros_cut,
        bins=50,
        range=(0, np.max(list_global)),
        color="red",
        histtype="stepfilled",
        alpha=0.3,
    )
    # Outline cut data
    plt.hist(
        list_cros_cut,
        bins=50,
        range=(0, np.max(list_global)),
        color="red",
        histtype="step",
        linewidth=3,
    )
    # mean and std
    plt.text(
        0.95,
        0.95,
        r"$\mu$ = "
        f"{np.mean(list_cros):.2f}\n"
        r"$\sigma$ = "
        f"{np.std(list_cros):.2f}",
        ha="right",
        va="top",
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", facecolor="blue"),
        color="white",
    )
    plt.text(
        0.95,
        0.75,
        r"$\mu$ = "
        f"{np.mean(list_cros_cut):.2f}\n"
        r"$\sigma$ = "
        f"{np.std(list_cros_cut):.2f}",
        ha="right",
        va="top",
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", facecolor="red"),
        color="white",
    )
    # Cuts
    plt.text(
        0.05,
        0.95,
        r"$0^\circ < \varphi < 90^\circ$" "\n L $>$ 0 mm",
        ha="left",
        va="top",
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", facecolor="blue"),
        color="white",
    )
    plt.text(
        0.05,
        0.75,
        r"$30^\circ \leq \varphi \leq 60^\circ$" "\n L $\\geq$ 2 mm",
        ha="left",
        va="top",
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", facecolor="red"),
        color="white",
    )
    plt.tight_layout()
    pdf.savefig(bbox_inches="tight")
    plt.clf()

    # Vertical heatmap #
    plt.imshow(
        arr_vert,
        cmap="viridis",
        origin="lower",
        aspect="auto",
        extent=(float(v_phi[0]), float(v_phi[-1]), float(v_d[0]), float(v_d[-1])),
    )
    cbar = plt.colorbar(pad=0.01, fraction=0.08)
    plt.clim(0, np.max(list_global))
    plt.xlabel(r" track angle [°]")
    plt.ylabel("impact parameter [mm]", labelpad=15)
    plt.tight_layout()
    pdf.savefig(bbox_inches="tight")
    plt.clf()

    # Distribution of cluster accuracy #
    plt.grid()
    ax = plt.gca()
    plt.xlabel(r"$A_{true}/A_{vertical}$")
    # Transparent all data
    nbins = 50
    plt.hist(
        list_vert,
        bins=nbins,
        range=(0.0, 2.0),
        color="blue",
        histtype="stepfilled",
        alpha=0.3,
    )
    # Outline all data
    plt.hist(
        list_vert,
        bins=nbins,
        range=(0.0, 2.0),
        color="blue",
        histtype="step",
        linewidth=3,
    )
    # Transparent cut data
    plt.hist(
        list_vert_cut,
        bins=nbins,
        range=(0.0, 2.0),
        color="red",
        histtype="stepfilled",
        alpha=0.3,
    )
    # Outline cut data
    plt.hist(
        list_vert_cut,
        bins=nbins,
        range=(0.0, 2.0),
        color="red",
        histtype="step",
        linewidth=3,
    )
    # mean and std
    plt.text(
        0.95,
        0.95,
        r"$\mu$ = "
        f"{np.mean(list_vert):.6f}\n"
        r"$\sigma$ = "
        f"{np.std(list_vert):.6f}",
        ha="right",
        va="top",
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", facecolor="blue"),
        color="white",
    )
    plt.text(
        0.95,
        0.75,
        r"$\mu$ = "
        f"{np.mean(list_vert_cut):.6f}\n"
        r"$\sigma$ = "
        f"{np.std(list_vert_cut):.6f}",
        ha="right",
        va="top",
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", facecolor="red"),
        color="white",
    )
    # Cuts
    plt.text(
        0.05,
        0.95,
        r"$0^\circ < \varphi < 90^\circ$" "\n L $>$ 0 mm",
        ha="left",
        va="top",
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", facecolor="blue"),
        color="white",
    )
    plt.text(
        0.05,
        0.75,
        r"$\varphi \leq 45^\circ$" "\n L $\\geq$ 2 mm",
        ha="left",
        va="top",
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", facecolor="red"),
        color="white",
    )
    plt.tight_layout()
    pdf.savefig(bbox_inches="tight")
    plt.clf()

    # Leading Pad heatmap #
    plt.imshow(
        arr_lead,
        cmap="viridis",
        origin="lower",
        aspect="auto",
        extent=(float(v_phi[0]), float(v_phi[-1]), float(v_d[0]), float(v_d[-1])),
    )
    cbar = plt.colorbar(pad=0.01, fraction=0.08)
    plt.clim(0, np.max(list_global))
    plt.xlabel(r" track angle [°]")
    plt.ylabel("impact parameter [mm]", labelpad=15)
    plt.tight_layout()
    pdf.savefig(bbox_inches="tight")
    plt.clf()

plt.close()
