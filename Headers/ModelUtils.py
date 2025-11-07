"""Physical model utilities for the HATPC lineic-charge project.

This module provides the core physics computation routines used by the
interactive visualizations. Functions compute the deposited charge (fC), the
current (fC/ns), and the shaped signal obtained by convolving with the
electronics transfer function (ETF). Units used throughout are:

- time: ns
- space: mm
- charge: fC

The module exposes helper wrappers `Compute0D` and `Compute1D` which dispatch
to the appropriate function for the chosen `vartype` ("Signal", "Charge",
"Current"). Physical and electronics constants are defined near the top of
the file.
"""

import numpy as np
import scipy.special as sc
from scipy import signal

from sys import path

path.append("Headers/")
from GeometryUtils import *

# Units: ns mm fC ---------------------------------------------------------------------------------------------------------------
# physics variables
t = np.linspace(0, 3000, 500)  # ns | start at 1 to avoid sigma = 0

G = 1493  # mean bHATPC
RC = 120  # ns/mm
PT = 412  # ns
Dt = 310 / np.power(10, 7 / 2)  # mm/sqrt(mm) CERN JF work (310 with B | 350 without B)
z = 250  # mm

# Charge variables
qe = 1.61e-4  # fC/e
# Punctual deposit case
Ne55Fe = 224 # number of electrons caused by a 55Fe decay
# Track case
Ne = 10  # e/mm
lambdaG = qe * G * Ne  # fC/mm

# Electronics variables
ws = 2 / PT
Q = 2 / 3
A = np.sqrt((2 * Q - 1) / (2 * Q + 1))
B = ws / 2 * np.sqrt(4 - 1 / Q**2)
C = ws / (2 * Q)


# Functions ---------------------------------------------------------------------------------------------------------------------
# Type of plot
def Compute0D(t, x0, y0, xmin, xmax, ymin, ymax, RC, z, type):
    """Dispatch wrapper for a point (0D) deposit.

    Parameters
    ----------
    t : array_like
        Time axis (ns).
    x0, y0 : float
        Drop coordinates (mm) relative to pad origins.
    xmin, xmax, ymin, ymax : float
        Cell boundaries (mm).
    RC : float
        RC (ns/mm) diffusion/resistance parameter.
    z : float
        Drift length (mm).
    type : str
        One of "Signal", "Charge", or "Current" determining which
        quantity is returned.

    Returns
    -------
    numpy.ndarray
        Array with the requested quantity evaluated on `t`.
    """

    if type == "Signal":
        return Signal0D(t, x0, y0, xmin, xmax, ymin, ymax, RC, z)
    elif type == "Charge":
        return Charge0D(t, x0, y0, xmin, xmax, ymin, ymax, RC, z)
    elif type == "Current":
        return Current0D(t, x0, y0, xmin, xmax, ymin, ymax, RC, z)


def Compute1D(t, m, q, xmin, xmax, ymin, ymax, RC, z, type):
    """Dispatch wrapper for a lineic (1D) deposit.

    Parameters
    ----------
    t : array_like
        Time axis (ns).
    m, q : float
        Line parameters y = m*x + q describing the track projection.
    xmin, xmax, ymin, ymax : float
        Cell boundaries (mm).
    RC : float
        RC (ns/mm) diffusion/resistance parameter.
    z : float
        Drift length (mm).
    type : str
        One of "Signal", "Charge", or "Current" determining which
        quantity is returned.

    Returns
    -------
    numpy.ndarray
        Array with the requested quantity evaluated on `t`.
    """

    if type == "Signal":
        return Signal1D(t, m, q, xmin, xmax, ymin, ymax, RC, z)
    elif type == "Charge":
        return Charge1D(t, m, q, xmin, xmax, ymin, ymax, RC, z)
    elif type == "Current":
        return Current1D(t, m, q, xmin, xmax, ymin, ymax, RC, z)


# Electronics functions
def Get_max_ETF(t):
    """Return the maximum of the electronics transfer function (ETF).

    The ETF is defined analytically below and this helper computes its
    maximum over the provided time array. The value is used to normalize the
    ETF so that signals are returned in ADC-like units.
    """

    ETF = np.heaviside(t, 1) * (
        np.exp(-ws * t) + np.exp(-C * t) * (A * np.sin(B * t) - np.cos(B * t))
    )
    return max(ETF)


max_ETF = Get_max_ETF(t)


def ETF(t):
    """Normalized electronics transfer function.

    The function returns the ETF sampled on `t` and normalized so that its
    peak corresponds to an ADC scale used elsewhere in the code (4096/120
    factor is kept from the original implementation).
    """

    ETF = np.heaviside(t, 1) * (
        np.exp(-ws * t) + np.exp(-C * t) * (A * np.sin(B * t) - np.cos(B * t))
    )
    return 4096 / 120 * ETF / max(ETF)


def dETFdt(t):
    """Time derivative of the normalized ETF.

    The derivative is used when convolving the deposited charge with the
    electronics response to obtain the shaped signal. The result is scaled
    by the same normalization factor as `ETF`.
    """

    dETFdt = np.heaviside(t, 1) * (
        -ws * np.exp(-ws * t)
        + np.exp(-C * t) * ((B - A * C) * np.sin(B * t) + (A * B + C) * np.cos(B * t))
    )
    return 4096 / 120 * dETFdt / max_ETF


# Charge functions
def Charge0D(t, x0, y0, xmin, xmax, ymin, ymax, RC, z):  # fC
    """Compute the integrated charge (fC) collected by a cell for a point
    deposit centered at (x0, y0).

    The model assumes a Gaussian transverse spread with time-dependent
    standard deviation sigma(t) which includes diffusion and an RC-dependent
    term. The integral over the rectangular pad is computed using error
    functions.
    """

    sigma = np.sqrt(2 * t / RC + Dt**2 * z)  # includes transverse diffusion
    erfx = sc.erf((xmax - x0) / (sigma * np.sqrt(2))) - sc.erf(
        (xmin - x0) / (sigma * np.sqrt(2))
    )
    erfy = sc.erf((ymax - y0) / (sigma * np.sqrt(2))) - sc.erf(
        (ymin - y0) / (sigma * np.sqrt(2))
    )
    return G * Ne55Fe * qe / 4 * erfx * erfy


def Charge1D(t, m, q, xmin, xmax, ymin, ymax, RC, z, Dt=Dt):  # fC
    """Compute the integrated charge (fC) collected by a cell for a
    lineic (track) deposit described by y = m*x + q.

    The expression is an analytic integral of the Gaussian-convolved linear
    charge density over the rectangular pad. `Dt` sets the transverse
    diffusion parametrization and can be overridden for testing.
    """

    sigma = np.sqrt(2 * t / RC + Dt**2 * z)  # includes transverse diffusion

    coeff1 = np.sqrt(2 * (1 + m**2) / np.pi) * sigma
    term11 = np.exp(-((-ymin + xmax * m + q) ** 2) / (2 * (1 + m**2) * sigma**2))
    term12 = np.exp(-((-ymin + xmin * m + q) ** 2) / (2 * (1 + m**2) * sigma**2))
    term13 = np.exp(-((-ymax + xmin * m + q) ** 2) / (2 * (1 + m**2) * sigma**2))
    term14 = np.exp(-((-ymax + xmax * m + q) ** 2) / (2 * (1 + m**2) * sigma**2))

    term21 = (ymin - xmin * m - q) * sc.erf(
        (-ymin + xmin * m + q) / (np.sqrt(2 * (1 + m**2)) * sigma)
    )
    term22 = (ymax - xmin * m - q) * sc.erf(
        (-ymax + xmin * m + q) / (np.sqrt(2 * (1 + m**2)) * sigma)
    )
    term23 = (ymin - xmax * m - q) * sc.erf(
        (-ymin + xmax * m + q) / (np.sqrt(2 * (1 + m**2)) * sigma)
    )
    term24 = (ymax - xmax * m - q) * sc.erf(
        (-ymax + xmax * m + q) / (np.sqrt(2 * (1 + m**2)) * sigma)
    )
    return (
        lambdaG
        * np.sqrt(1 + m**2)
        / (2 * m)
        * (
            coeff1 * (term11 - term12 + term13 - term14)
            + term21
            - term22
            - term23
            + term24
        )
    )


# Current functions in fC/ns = µA
def Current0D(t, x0, y0, xmin, xmax, ymin, ymax, RC, z):
    """Approximate current (fC/ns) for a point deposit via finite differences.

    The function uses a simple forward-difference approximation on the
    `t` grid. The units are fC/ns (equivalent to µA when multiplied by 1e3).
    """

    return (
        Charge0D(t + 1, x0, y0, xmin, xmax, ymin, ymax, RC, z)
        - Charge0D(t, x0, y0, xmin, xmax, ymin, ymax, RC, z)
    ) / np.diff(t)[0]


def Current1D(t, m, q, xmin, xmax, ymin, ymax, RC, z):
    """Approximate current (fC/ns) for a lineic deposit via finite differences."""

    return (
        Charge1D(t + 1, m, q, xmin, xmax, ymin, ymax, RC, z)
        - Charge1D(t, m, q, xmin, xmax, ymin, ymax, RC, z)
    ) / np.diff(t)[0]


# Signal functions
def Signal0D(t, x0, y0, xmin, xmax, ymin, ymax, RC, z):
    """Shaped signal for a point deposit: convolution of charge with ETF'"""

    return (
        np.convolve(
            Charge0D(t, x0, y0, xmin, xmax, ymin, ymax, RC, z), dETFdt(t), mode="full"
        )
        * np.diff(t)[0]
    )


def Signal1D(t, m, q, xmin, xmax, ymin, ymax, RC, z, Dt=Dt):
    """Shaped signal for a line deposit: convolution of charge with ETF'"""

    return (
        np.convolve(
            Charge1D(t, m, q, xmin, xmax, ymin, ymax, RC, z, Dt=Dt),
            dETFdt(t),
            mode="full",
        )
        * np.diff(t)[0]
    )


dETFdt_t = dETFdt(t)


def Signal1DFFT(t, m, q, xmin, xmax, ymin, ymax, RC, z):
    """An FFT-based convolution variant of Signal1D.

    This implementation uses `scipy.signal.fftconvolve` which can be faster for
    long arrays. The returned array is truncated to the original time length.
    """

    result = (
        signal.fftconvolve(
            Charge1D(t, m, q, xmin, xmax, ymin, ymax, RC, z), dETFdt(t), mode="full"
        )
        * np.diff(t)[0]
    )
    return result[: len(t)]
