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
    if type == "Signal":
        return Signal0D(t, x0, y0, xmin, xmax, ymin, ymax, RC, z)
    elif type == "Charge":
        return Charge0D(t, x0, y0, xmin, xmax, ymin, ymax, RC, z)
    elif type == "Current":
        return Current0D(t, x0, y0, xmin, xmax, ymin, ymax, RC, z)


def Compute1D(t, m, q, xmin, xmax, ymin, ymax, RC, z, type):
    if type == "Signal":
        return Signal1D(t, m, q, xmin, xmax, ymin, ymax, RC, z)
    elif type == "Charge":
        return Charge1D(t, m, q, xmin, xmax, ymin, ymax, RC, z)
    elif type == "Current":
        return Current1D(t, m, q, xmin, xmax, ymin, ymax, RC, z)


# Electronics functions
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


def dETFdt(t):
    # max_ETF = max(ETF(t))
    dETFdt = np.heaviside(t, 1) * (
        -ws * np.exp(-ws * t)
        + np.exp(-C * t) * ((B - A * C) * np.sin(B * t) + (A * B + C) * np.cos(B * t))
    )
    return 4096 / 120 * dETFdt / max_ETF


# Charge functions
def Charge0D(t, x0, y0, xmin, xmax, ymin, ymax, RC, z):  # fC
    sigma = np.sqrt(2 * t / RC + Dt**2 * z)  # includes transverse diffusion
    erfx = sc.erf((xmax - x0) / (sigma * np.sqrt(2))) - sc.erf(
        (xmin - x0) / (sigma * np.sqrt(2))
    )
    erfy = sc.erf((ymax - y0) / (sigma * np.sqrt(2))) - sc.erf(
        (ymin - y0) / (sigma * np.sqrt(2))
    )
    return G * Ne55Fe * qe / 4 * erfx * erfy


def Charge1D(t, m, q, xmin, xmax, ymin, ymax, RC, z, Dt=Dt):  # fC
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
    return (
        Charge0D(t + 1, x0, y0, xmin, xmax, ymin, ymax, RC, z)
        - Charge0D(t, x0, y0, xmin, xmax, ymin, ymax, RC, z)
    ) / np.diff(t)[0]


def Current1D(t, m, q, xmin, xmax, ymin, ymax, RC, z):
    return (
        Charge1D(t + 1, m, q, xmin, xmax, ymin, ymax, RC, z)
        - Charge1D(t, m, q, xmin, xmax, ymin, ymax, RC, z)
    ) / np.diff(t)[0]


# Signal functions
def Signal0D(t, x0, y0, xmin, xmax, ymin, ymax, RC, z):
    return (
        np.convolve(
            Charge0D(t, x0, y0, xmin, xmax, ymin, ymax, RC, z), dETFdt(t), mode="full"
        )
        * np.diff(t)[0]
    )


def Signal1D(t, m, q, xmin, xmax, ymin, ymax, RC, z, Dt=Dt):
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
    result = (
        signal.fftconvolve(
            Charge1D(t, m, q, xmin, xmax, ymin, ymax, RC, z), dETFdt(t), mode="full"
        )
        * np.diff(t)[0]
    )
    return result[: len(t)]
