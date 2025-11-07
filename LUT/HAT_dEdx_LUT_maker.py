# Author: Tristan DARET
# Creation Date: 2024-05-16
# Last modification: 2024-11-01 by Tristan DARET

"""
Usage:
  python dEdx_XP_LUT_maker.py <transverse_diffusion_coefficient> <RC_value>

Arguments:
  <transverse_diffusion_coefficient> : The transverse diffusion coefficient (Dt) to be used in the LUT computation.
  <RC_value>                         : The RC value of the ERAMs to be used in the LUT computation.

Example:
  python dEdx_XP_LUT_maker.py 310 120

This will compute the Look Up Tables (LUTs) for the given transverse diffusion coefficient and RC value, and save the results in a ROOT file.
To compute the overall LUTs, the script should be run for Dt = (310, 350) and RC = (112, 158), then hadd the 4 files.
"""

"""
This script computes the Look Up Tables for the Crossed Pads (XP) method to get dE/dx with the ERAM modules of HATPCs.
This script is fully independent from the rest of the code and can be run on its own.
The LUTs are computed for a given set of parameters (Dt, PT, nphi, nd, nRC, nZ) and saved in a ROOT file within a TTree.
These parameters correspond to the transverse diffusion coefficient, the peaking time of the electronics, 
the track angle in a given pad, the impact parameter of the track, the RC value of the ERAMs, and the drift distance.

Computation details are in the backup slides here: https://t2k.org/nd280/physics/nd280-ccnue-em-working-group/meetings/cm_march/PiD_HATPC

With the current step sizes, the LUTs are computed in about 30 minutes in local on a recent PC.
Each value of the LUT is independent from the others, so the computation can be parallelized on a cluster if someone wishes to do so.

The LUTs were already computed and are supposed to be found in https://nd280.lancs.ac.uk/downloads/nd280files/hatTemplates/
They are downloaded automatically when sourcing for the first time. 
If for some reason, the repo disappears and every single user lost their copy of the LUTs, they can be remade with this script.
"""

import sys  # pass arguments for parallelization
import numpy as np
import scipy.special as sc  # for the error function
from ROOT import TFile, TTree
from array import array  # to store the LUT in a TTree

import time

# Output file directory (adapt it to your needs)
out_dir = "LUT/"

# Units: ns mm fC ---------------------------------------------------------------------------------------------------------------
# physics variables
t = np.linspace(1, 3000, 500)  # ns | start at 1 to avoid sigma = 0

# Electronics variables
PT = 412  # ns
ws = 2 / PT
Q = 2 / 3
A = np.sqrt((2 * Q - 1) / (2 * Q + 1))
B = ws / 2 * np.sqrt(4 - 1 / Q**2)
C = ws / (2 * Q)

# Geometry variables
nX = 5  # number of colums
xwidth = 11.28  # mm ; width of xmin pad
xc = xwidth / 2  # mm ; horizontal center of pad
xleft = (
    -(nX // 2) * xwidth
)  # Most left position of the superpad considered (wrt left of LP)
xright = (
    nX // 2 + 1
) * xwidth  # Most right position of the superpad considered (wrt left of LP)
xmin = 0  # mm ; bottom border of leading pad
xmax = xwidth  # mm ; top border of leading pad

nY = 5  # number of rows
ywidth = 10.19  # mm ; height of xmin pad
yc = ywidth / 2  # mm ; vertical center of pad
ylow = (
    -(nY // 2) * ywidth
)  # Lowest position of the superpad considered (wrt bottom of LP)
yhigh = (
    nY // 2 + 1
) * ywidth  # Highest position of the superpad considered (wrt bottom of LP)
ymin = 0  # mm ; left border of leading pad
ymax = ywidth  # mm ; right border of leading pad

diag = np.sqrt(xwidth**2 + ywidth**2)


# Functions ---------------------------------------------------------------------------------------------------------------------
# Charge function
# For linear track simulation. The charge is deposited uniformly along the track
# The track is defined as a straight line, which is correct for a local (pad-wide) approximation
def Charge(t, m, q, i, j, k, l, RC, drift, Dt):  # fC
    sigma = np.sqrt(2 * t / RC + Dt**2 * drift)  # includes transverse diffusion

    coeff1 = np.sqrt(2 * (1 + m**2) / np.pi) * sigma
    term11 = np.exp(-((-k + j * m + q) ** 2) / (2 * (1 + m**2) * sigma**2))
    term12 = np.exp(-((-k + i * m + q) ** 2) / (2 * (1 + m**2) * sigma**2))
    term13 = np.exp(-((-l + i * m + q) ** 2) / (2 * (1 + m**2) * sigma**2))
    term14 = np.exp(-((-l + j * m + q) ** 2) / (2 * (1 + m**2) * sigma**2))

    term21 = (k - i * m - q) * sc.erf(
        (-k + i * m + q) / (np.sqrt(2 * (1 + m**2)) * sigma)
    )
    term22 = (l - i * m - q) * sc.erf(
        (-l + i * m + q) / (np.sqrt(2 * (1 + m**2)) * sigma)
    )
    term23 = (k - j * m - q) * sc.erf(
        (-k + j * m + q) / (np.sqrt(2 * (1 + m**2)) * sigma)
    )
    term24 = (l - j * m - q) * sc.erf(
        (-l + j * m + q) / (np.sqrt(2 * (1 + m**2)) * sigma)
    )
    return (
        np.sqrt(1 + m**2)
        / (2 * m)
        * (
            coeff1 * (term11 - term12 + term13 - term14)
            + term21
            - term22
            - term23
            + term24
        )
    )


# Electronics transfer Function (ETF)
# Maths can be found here: https://thesis.unipd.it/handle/20.500.12608/21505
# Estimate normalization value for the transfer function of the electronics
def Get_max_ETF(t):
    ETF = np.heaviside(t, 1) * (
        np.exp(-ws * t) + np.exp(-C * t) * (A * np.sin(B * t) - np.cos(B * t))
    )
    return max(ETF)


max_ETF = Get_max_ETF(t)


def ETF(t):
    ETF = np.heaviside(t, 1) * (
        np.exp(-ws * t) + np.exp(-C * t) * (A * np.sin(B * t) - np.cos(B * t))
    )
    return 4096 / 120 * ETF / max(ETF)


# Need to convolute the transfer function with the charge derivative, but the derivative commutes with the convolution
# and it's easier to compute the derivative of the transfer function
def dETFdt(t):
    dETFdt = np.heaviside(t, 1) * (
        -ws * np.exp(-ws * t)
        + np.exp(-C * t) * ((B - A * C) * np.sin(B * t) + (A * B + C) * np.cos(B * t))
    )
    return 4096 / 120 * dETFdt / max_ETF


# Fix the value to avoid unnecessary computation
dETFdt_t = dETFdt(t)


# Convolution functions
def Signal(t, m, q, xmin, xmax, ymin, ymax, RC, drift, Dt):
    return (
        np.convolve(
            Charge(t, m, q, xmin, xmax, ymin, ymax, RC, drift, Dt),
            dETFdt(t),
            mode="full",
        )
        * np.diff(t)[0]
    )


# Geometry functions
# returns x (respectively y) for a given y (respectively x), angle and impact parameter
def X(phi_rad, d, y):
    return (
        y - (d - np.sin(phi_rad) * xc + np.cos(phi_rad) * yc) / np.cos(phi_rad)
    ) / np.tan(phi_rad)


def Y(phi_rad, d, x):
    return np.tan(phi_rad) * x + (
        d - np.sin(phi_rad) * xc + np.cos(phi_rad) * yc
    ) / np.cos(phi_rad)


start_time = time.time()

# LUT computation ---------------------------------------------------------------------------------------------------------------
# LUT parameters
ETF = ETF(t)
nphi = 25
nd = 25
nZ = 51
Dt = float(sys.argv[1])
RC = float(sys.argv[2])
print(f"Computing LUT for Dt = {Dt} and RC = {RC}")

arr_r = np.full((nd, nphi), np.nan)
v_d = np.linspace(0, diag / 2, nd)
v_phi = np.linspace(1e-6, 90 - 1e-6, nphi)
v_Z = np.linspace(0, 1000, nZ)

out_file = TFile(f"LUT/dEdx_XP_LUT_tmp_Dt{Dt:.0f}_RC{RC:.0f}.root", "RECREATE")
out_tree = TTree("outTree", "LUT")
Dt_array = array("f", [0])
RC_array = array("f", [0])
phi_array = array("f", [0])
d_array = array("f", [0])
z_array = array("f", [0])
weight_array = array("f", [0])
out_tree.Branch("transDiff", Dt_array, "transDiff/F")
out_tree.Branch("RC", RC_array, "RC/F")
out_tree.Branch("angle", phi_array, "angle/F")
out_tree.Branch("impact_param", d_array, "impact_param/F")
out_tree.Branch("drift_dist", z_array, "drift_dist/F")
out_tree.Branch("weight", weight_array, "weight/F")


# Make Length map
phi_index = 0
for phi in v_phi:
    phi_rad = phi / 180 * np.pi
    d_index = 0

    for d in v_d:
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
        arr_r[d_index, phi_index] = L
        d_index += 1

    phi_index += 1
print(f"Length map done in {time.time()-start_time:.1f} seconds")


# Make LUT for all (RC, z, phi, d)
LUT_time = time.time()
Dt_array[0] = Dt
RC_array[0] = RC
for z in v_Z:
    print(f"z = {z:.0f} mm")
    phi_index = 0
    z_array[0] = z

    for phi in v_phi:
        # print(f"phi = {phi:.1f}°")
        d_index = 0
        phi_array[0] = phi

        for d in v_d:
            # print(f"d = {d:.2f} mm")
            d_array[0] = d
            phi_rad = phi / 180 * np.pi
            m = np.tan(phi_rad)
            q = (np.cos(phi_rad) * yc - np.sin(phi_rad) * xc + d) / np.cos(
                phi_rad
            )  # intercept

            ETFr = arr_r[d_index, phi_index] * np.max(ETF)

            ADC = np.max(
                Signal(
                    t, m, q, xmin, xmax, ymin, ymax, RC, z, Dt / np.power(10, 7 / 2)
                )[: len(t)]
            )
            weight_array[0] = ETFr / ADC
            if not (np.isnan(weight_array[0])) and weight_array[0] > 0:
                out_tree.Fill()

            d_index += 1
        phi_index += 1
print(f"LUT done in {time.time()-LUT_time:.1f} seconds")

print(f"Total time: {time.time()-start_time:.1f} seconds")

out_tree.Write()
out_file.Close()
