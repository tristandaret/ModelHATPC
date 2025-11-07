"""Estimate cluster amplitude ratios and produce diagnostic heatmaps.

This script scans track angles and impact parameters across a small
sub-map and computes the ratios between the true deposited amplitude and
different cluster estimators (diagonal, half-cross, vertical, leading pad).
It stores several diagnostic heatmaps and histograms into a multi-page PDF
in `Illustrations/`.

Note: the script can be computationally intensive depending on the grid
resolution (`nsteps`) and the number of tracks considered.
"""

from sys import path

path.append("Headers/")
from ModelUtils import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Plotting style with LaTeX
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": 35})


# Computations
nsteps = 250
RC = 120  # Readout chip radius in mm
z = 250  # Distance from the readout chip to the sensor in mm
ETF = lambdaG * ETF(t)
nd = nsteps
nphi = nsteps
v_d = np.linspace(0, diag / 2, nd)
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
        m, q, phi_rad = compute_line_params(phi, d)

        # Determine the cluster lengths
        r_diag, r_vert, r_cros, L = ClusterLengths(phi_rad, d)

        ETF_lead = r_diag * np.max(ETF)
        ETF_diag = r_diag * np.max(ETF)
        ETF_cros = r_cros * np.max(ETF)
        ETF_vert = r_vert * np.max(ETF)
        Lead_only = np.zeros(t.size)
        Clus_Diag = np.zeros(t.size)
        Clus_HC = np.zeros(t.size)
        Clus_vert = np.zeros(t.size)
        for iX in range(nX):
            for iY in range(nY):
                ADC = Signal1D(
                    t,
                    m,
                    q,
                    xleft + iX * xwidth,
                    xleft + (iX + 1) * xwidth,
                    ylow + iY * ywidth,
                    ylow + (iY + 1) * ywidth,
                    RC,
                    z,
                )[: len(t)]
                # Leading pad only
                if iX == nX // 2 and iY == nY // 2:
                    Lead_only += ADC
                # Cluster distributions
                # Diagonal distribution
                if iX + iY == nX // 2 + nY // 2:
                    Clus_Diag += ADC
                # Half-cross distribution
                if (iX >= nX // 2 and iY == nY // 2) or (
                    iX == nX // 2 and iY >= nY // 2
                ):
                    Clus_HC += ADC
                # Vertical distribution
                if iX == nX // 2:
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
    f"Illustrations/Ratio_Diag_HC_nphi{nphi:d}_nd{nd:d}_RC{RC:d}_z{z:d}_PT{PT:d}_Dt{Dt*np.power(10,7/2):.0f}_2mm_30phi60_nX{nX:d}_nY{nY:d}.pdf"
) as pdf:

    ### Diagonal heatmap ###
    plt.figure(figsize=(15, 10))
    plt.imshow(
        arr_diag,
        cmap="viridis",
        origin="lower",
        aspect="auto",
        extent=[v_phi[0], v_phi[-1], v_d[0], v_d[-1]],
    )
    cbar = plt.colorbar(pad=0.01, fraction=0.08)
    plt.clim(0, np.max(list_global))
    plt.xlabel(r" track angle [°]")
    plt.ylabel("impact parameter [mm]", labelpad=15)
    plt.tight_layout()
    pdf.savefig(bbox_inches="tight")
    plt.clf()

    ### Distribution of cluster accuracy ###
    plt.grid()
    ax = plt.gca()
    plt.xlabel(r"$A_{true}/A_{diagonal}$")
    # Transparent all data
    plt.hist(
        list_diag,
        bins=50,
        range=[0, np.max(list_global)],
        color="blue",
        histtype="stepfilled",
        alpha=0.3,
    )
    # Outline all data
    plt.hist(
        list_diag,
        bins=50,
        range=[0, np.max(list_global)],
        color="blue",
        histtype="step",
        linewidth=3,
    )
    # Transparent cut data
    plt.hist(
        list_diag_cut,
        bins=50,
        range=[0, np.max(list_global)],
        color="red",
        histtype="stepfilled",
        alpha=0.3,
    )
    # Outline cut data
    plt.hist(
        list_diag_cut,
        bins=50,
        range=[0, np.max(list_global)],
        color="red",
        histtype="step",
        linewidth=3,
    )
    # mean and std
    plt.text(
        0.95,
        0.95,
        "$\mu$ = "
        f"{np.mean(list_diag):.2f}\n"
        "$\sigma$ = "
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
        "$\mu$ = "
        f"{np.mean(list_diag_cut):.2f}\n"
        "$\sigma$ = "
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
        r"$30^\circ \leq \varphi \leq 60^\circ$" "\n L $\geq$ 2 mm",
        ha="left",
        va="top",
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", facecolor="red"),
        color="white",
    )
    plt.tight_layout()
    pdf.savefig(bbox_inches="tight")
    plt.clf()

    ####################### Half-cross heatmap ############################
    plt.imshow(
        arr_cros,
        cmap="viridis",
        origin="lower",
        aspect="auto",
        extent=[v_phi[0], v_phi[-1], v_d[0], v_d[-1]],
    )
    cbar = plt.colorbar(pad=0.01, fraction=0.08)
    plt.clim(0, np.max(list_global))
    plt.xlabel(r" track angle [°]")
    plt.ylabel("impact parameter [mm]", labelpad=15)
    plt.tight_layout()
    pdf.savefig(bbox_inches="tight")
    plt.clf()

    ### Distribution of cluster accuracy ###
    plt.grid()
    ax = plt.gca()
    plt.xlabel(r"$A_{true}/A_{half-cross}$")
    # Transparent all data
    plt.hist(
        list_cros,
        bins=50,
        range=[0, np.max(list_global)],
        color="blue",
        histtype="stepfilled",
        alpha=0.3,
    )
    # Outline all data
    plt.hist(
        list_cros,
        bins=50,
        range=[0, np.max(list_global)],
        color="blue",
        histtype="step",
        linewidth=3,
    )
    # Transparent cut data
    plt.hist(
        list_cros_cut,
        bins=50,
        range=[0, np.max(list_global)],
        color="red",
        histtype="stepfilled",
        alpha=0.3,
    )
    # Outline cut data
    plt.hist(
        list_cros_cut,
        bins=50,
        range=[0, np.max(list_global)],
        color="red",
        histtype="step",
        linewidth=3,
    )
    # mean and std
    plt.text(
        0.95,
        0.95,
        "$\mu$ = "
        f"{np.mean(list_cros):.2f}\n"
        "$\sigma$ = "
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
        "$\mu$ = "
        f"{np.mean(list_cros_cut):.2f}\n"
        "$\sigma$ = "
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
        r"$30^\circ \leq \varphi \leq 60^\circ$" "\n L $\geq$ 2 mm",
        ha="left",
        va="top",
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", facecolor="red"),
        color="white",
    )
    plt.tight_layout()
    pdf.savefig(bbox_inches="tight")
    plt.clf()

    ####################### Vertical heatmap ############################
    plt.imshow(
        arr_vert,
        cmap="viridis",
        origin="lower",
        aspect="auto",
        extent=[v_phi[0], v_phi[-1], v_d[0], v_d[-1]],
    )
    cbar = plt.colorbar(pad=0.01, fraction=0.08)
    plt.clim(0, np.max(list_global))
    plt.xlabel(r" track angle [°]")
    plt.ylabel("impact parameter [mm]", labelpad=15)
    plt.tight_layout()
    pdf.savefig(bbox_inches="tight")
    plt.clf()

    ### Distribution of cluster accuracy ###
    plt.grid()
    ax = plt.gca()
    plt.xlabel(r"$A_{true}/A_{vertical}$")
    # Transparent all data
    nbins = 50
    plt.hist(
        list_vert,
        bins=nbins,
        range=[0, 2],
        color="blue",
        histtype="stepfilled",
        alpha=0.3,
    )
    # Outline all data
    plt.hist(
        list_vert,
        bins=nbins,
        range=[0, 2],
        color="blue",
        histtype="step",
        linewidth=3,
    )
    # Transparent cut data
    plt.hist(
        list_vert_cut,
        bins=nbins,
        range=[0, 2],
        color="red",
        histtype="stepfilled",
        alpha=0.3,
    )
    # Outline cut data
    plt.hist(
        list_vert_cut,
        bins=nbins,
        range=[0, 2],
        color="red",
        histtype="step",
        linewidth=3,
    )
    # mean and std
    plt.text(
        0.95,
        0.95,
        "$\mu$ = "
        f"{np.mean(list_vert):.6f}\n"
        "$\sigma$ = "
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
        "$\mu$ = "
        f"{np.mean(list_vert_cut):.6f}\n"
        "$\sigma$ = "
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
        r"$\varphi \leq 45^\circ$" "\n L $\geq$ 2 mm",
        ha="left",
        va="top",
        transform=ax.transAxes,
        bbox=dict(boxstyle="round", facecolor="red"),
        color="white",
    )
    plt.tight_layout()
    pdf.savefig(bbox_inches="tight")
    plt.clf()

    ####################### Leading Pad heatmap ############################
    plt.imshow(
        arr_lead,
        cmap="viridis",
        origin="lower",
        aspect="auto",
        extent=[v_phi[0], v_phi[-1], v_d[0], v_d[-1]],
    )
    cbar = plt.colorbar(pad=0.01, fraction=0.08)
    plt.clim(0, np.max(list_global))
    plt.xlabel(r" track angle [°]")
    plt.ylabel("impact parameter [mm]", labelpad=15)
    plt.tight_layout()
    pdf.savefig(bbox_inches="tight")
    plt.clf()

plt.close()
