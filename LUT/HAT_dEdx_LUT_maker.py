"""HAT dE/dx LUT maker.

This module computes Look-Up Tables (LUTs) for the Crossed Pads (XP) method
used to estimate dE/dx for HATPC ERAM modules. The LUTs are generated for a
grid of transverse diffusion coefficients, readout RC values, track angles,
impact parameters and drift distances and stored into a ROOT `TFile` with a
`TTree` for downstream analysis.

Notes
-----
- Run as a script: ``python LUT/HAT_dEdx_LUT_maker.py <Dt> <RC>`` where ``Dt`` is
    the transverse diffusion coefficient and ``RC`` the electronics RC value.
- Computation may be time-consuming depending on grid sizes; values are
    independent and can be parallelized.

Author
------
Tristan DARET
"""

import sys  # pass arguments for parallelization
import numpy as np
import scipy.special as sc  # for the error function
from ROOT import TFile, TTree
from array import array  # to store the LUT in a TTree

import time

# Output file directory (adapt it to your needs)
out_dir = "LUT/"

# ---- Units: ns mm fC ----
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


# ---- Functions ----
"""Low-level helper functions used by the LUT maker.

Functions in this section implement the analytic integrals and electronics
transfer function used when converting deposited charge into shaped ADC-like
signals. Each function follows NumPy-style docstrings.
"""


def Charge(t, m, q, i, j, k, l, RC, drift, Dt):  # fC
    """Compute integrated charge for a linear (track) deposit on a pad.

    This analytic expression integrates the Gaussian-convolved lineic
    charge density over the rectangular pad defined by the coordinates
    (i, j, k, l).

    Parameters
    ----------
    t : array_like
        Time axis in ns.
    m, q : float
        Line parameters describing the track projection (y = m*x + q).
    i, j, k, l : float
        Pad edge coordinates used in the integral.
    RC : float
        RC parameter (ns/mm) controlling transverse spread with time.
    drift : float
        Drift distance (mm).
    Dt : float
        Transverse diffusion coefficient.

    Returns
    -------
    ndarray
        Integrated charge (fC) sampled on `t`.
    """
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
    """Return the maximum of the electronics transfer function (ETF).

    Parameters
    ----------
    t : array_like
        Time axis in ns where the ETF is sampled.

    Returns
    -------
    float
        Maximum ETF value on `t`.
    """
    ETF = np.heaviside(t, 1) * (
        np.exp(-ws * t) + np.exp(-C * t) * (A * np.sin(B * t) - np.cos(B * t))
    )
    return max(ETF)


max_ETF = Get_max_ETF(t)


def ETF(t):
    """Return the normalized electronics transfer function.

    The ETF is normalized so its peak matches the ADC-like scale used in this
    code (factor 4096/120 preserved from the original implementation).

    Parameters
    ----------
    t : array_like
        Time axis in ns.

    Returns
    -------
    ndarray
        ETF sampled on `t` and normalized.
    """
    ETF = np.heaviside(t, 1) * (
        np.exp(-ws * t) + np.exp(-C * t) * (A * np.sin(B * t) - np.cos(B * t))
    )
    return 4096 / 120 * ETF / max(ETF)


# Need to convolute the transfer function with the charge derivative, but the derivative commutes with the convolution
# and it's easier to compute the derivative of the transfer function
def dETFdt(t):
    """Time derivative of the normalized ETF.

    The derivative is used when convolving the deposited charge with the
    electronics response to obtain the shaped signal.

    Parameters
    ----------
    t : array_like
        Time axis in ns.

    Returns
    -------
    ndarray
        Time derivative of the normalized ETF sampled on `t`.
    """
    dETFdt = np.heaviside(t, 1) * (
        -ws * np.exp(-ws * t)
        + np.exp(-C * t) * ((B - A * C) * np.sin(B * t) + (A * B + C) * np.cos(B * t))
    )
    return 4096 / 120 * dETFdt / max_ETF


# Fix the value to avoid unnecessary computation
dETFdt_t = dETFdt(t)


# Convolution functions
def Signal(t, m, q, xmin, xmax, ymin, ymax, RC, drift, Dt):
    """Convolve deposited charge with ETF' to obtain shaped signal.

    Parameters
    ----------
    t : array_like
        Time axis in ns.
    m, q : float
        Line parameters describing the track projection.
    xmin, xmax, ymin, ymax : float
        Pad boundaries (mm).
    RC : float
        RC parameter (ns/mm).
    drift : float
        Drift distance (mm).
    Dt : float
        Transverse diffusion coefficient.

    Returns
    -------
    ndarray
        Shaped signal (ADC-like) sampled on `t`.
    """
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
    """Return x coordinate where the projected line crosses a given y.

    Parameters
    ----------
    phi_rad : float
        Track angle in radians.
    d : float
        Impact parameter (mm).
    y : float
        y coordinate to evaluate (mm).

    Returns
    -------
    float
        x coordinate (mm) where the projected line crosses `y`.
    """
    return (
        y - (d - np.sin(phi_rad) * xc + np.cos(phi_rad) * yc) / np.cos(phi_rad)
    ) / np.tan(phi_rad)


def Y(phi_rad, d, x):
    """Return y coordinate for the projected line at a given x.

    Parameters
    ----------
    phi_rad : float
        Track angle in radians.
    d : float
        Impact parameter (mm).
    x : float
        x coordinate to evaluate (mm).

    Returns
    -------
    float
        y coordinate (mm) for the projected line at `x`.
    """
    return np.tan(phi_rad) * x + (
        d - np.sin(phi_rad) * xc + np.cos(phi_rad) * yc
    ) / np.cos(phi_rad)


start_time = time.time()

# LUT computation ---------------------------------------------------------------------------------------------------------------
# LUT parameters
ETF_samples = ETF(t)
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

            ETFr = arr_r[d_index, phi_index] * np.max(ETF_samples)

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
