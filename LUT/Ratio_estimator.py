"""Simple ratio estimator utilities used during LUT development.

Provides a function `ratio` that computes the ratio between an estimated
energy deposition (from geometry, `L`) and the maximal ADC amplitude from
the model. This helper was used to validate and debug LUT entries.
"""

from sys import path

path.append("Headers/")
from ModelUtils import *

ETF = lambdaG * ETF(t)


def ratio(RC, z, d, phi, L):
    """Estimate the ratio between geometry-based energy estimate and model ADC.

    The function computes the maximal ADC amplitude for a given projected
    line (defined by angle ``phi`` and impact parameter ``d``) and returns
    the ratio ``ETFr/ADC`` where ``ETFr`` is the geometry-derived estimate.

    Parameters
    ----------
    RC : float
        Readout RC parameter (ns/mm).
    z : float
        Drift distance (mm).
    d : float
        Impact parameter (mm).
    phi : float
        Track angle in degrees.
    L : float
        Geometry-derived length used to scale ETF.

    Returns
    -------
    float
        The ratio ETFr / ADC (unitless).
    """
    phi_rad = phi / 180 * np.pi
    m = np.tan(phi_rad)
    q = (d - xc * np.sin(phi_rad) + yc * np.cos(phi_rad)) / np.cos(phi_rad)
    ETFr = L * np.max(ETF)
    ADC = np.max(Signal1D(t, m, q, xmin, xmax, ymin, ymax, RC, z)[: len(t)])

    # Determine intersection points of the line with the central pad
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

    L = np.sqrt((y[1] - y[0]) ** 2 + (x[1] - x[0]) ** 2)
    print(f"{ETFr/ADC:.5f} = {ETFr:.2f}/{ADC:.2f} | L = {L:.3f}")
    return ETFr / ADC


def R(ra, rb, val, valmin, deltaval):
    """Linear interpolation helper.

    Parameters
    ----------
    ra, rb : float
        Values at the lower and upper bounds.
    val : float
        Query value.
    valmin : float
        Minimum of the value range.
    deltaval : float
        Range span (rb - ra corresponds to this span).

    Returns
    -------
    float
        Interpolated value corresponding to ``val``.
    """
    return ra + (val - valmin) / deltaval * (rb - ra)


# r0000 = ratio(125, 300, 5.22, 76.43, 6.16592)
# r0001 = ratio(125, 300, 5.22, 76.88, 6.16592)
# r0010 = ratio(125, 300, 5.30, 76.43, 6.16592)
# r0011 = ratio(125, 300, 5.30, 76.88, 6.16592)
# r0100 = ratio(125, 350, 5.22, 76.43, 6.16592)
# r0101 = ratio(125, 350, 5.22, 76.88, 6.16592)
# r0110 = ratio(125, 350, 5.30, 76.43, 6.16592)
# r0111 = ratio(125, 350, 5.30, 76.88, 6.16592)
# r1000 = ratio(130, 300, 5.22, 76.43, 6.16592)
# r1001 = ratio(130, 300, 5.22, 76.88, 6.16592)
# r1010 = ratio(130, 300, 5.30, 76.43, 6.16592)
# r1011 = ratio(130, 300, 5.30, 76.88, 6.16592)
# r1100 = ratio(130, 350, 5.22, 76.43, 6.16592)
# r1101 = ratio(130, 350, 5.22, 76.88, 6.16592)
# r1110 = ratio(130, 350, 5.30, 76.43, 6.16592)
# r1111 = ratio(130, 350, 5.30, 76.88, 6.16592)

# r000 = R(r0000, r0001, 76.78, 76.43, 0.45)
# r001 = R(r0010, r0011, 76.78, 76.43, 0.45)
# r010 = R(r0100, r0101, 76.78, 76.43, 0.45)
# r011 = R(r0110, r0111, 76.78, 76.43, 0.45)
# r100 = R(r1000, r1001, 76.78, 76.43, 0.45)
# r101 = R(r1010, r1011, 76.78, 76.43, 0.45)
# r110 = R(r1100, r1101, 76.78, 76.43, 0.45)
# r111 = R(r1110, r1111, 76.78, 76.43, 0.45)
# print(r0000, r0001, r000)
# print(r0010, r0011, r001)
# print(r0100, r0101, r010)
# print(r0110, r0111, r011)
# print(r1000, r1001, r100)
# print(r1010, r1011, r101)
# print(r1100, r1101, r110)
# print(r1110, r1111, r111)
# print()

# r00  = R(r000, r001, 5.28, 5.22, 0.08)
# r01  = R(r010, r011, 5.28, 5.22, 0.08)
# r10  = R(r100, r101, 5.28, 5.22, 0.08)
# r11  = R(r110, r111, 5.28, 5.22, 0.08)
# print(r000, r001, r00)
# print(r010, r011, r01)
# print(r100, r101, r10)
# print(r110, r111, r11)
# print()

# r0   = R(r00, r01, 332.8, 300, 50)
# r1   = R(r10, r11, 332.8, 300, 50)
# print(r00, r01, r0)
# print(r10, r11, r1)
# print()

# r_interpol  = R(r0, r1, 129.733, 125, 5)
# print(r0, r1, r_interpol)
# print()

# r_exact     = ratio(129.733, 332.8, 5.28314, -76.7846, 6.16592)

# print(r_exact, r_interpol)
