"""Simple ratio estimator utilities used during LUT development.

Provides a function `ratio` that computes the ratio between an estimated
energy deposition (from geometry, `L`) and the maximal ADC amplitude from
the model. This helper was used to validate and debug LUT entries.
"""

from sys import path

path.append("Headers/")
from ModelUtils import *

# compute ETF(t) once and store in a distinct variable to avoid shadowing
ETF_t = lambdaG * ETF(t)


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
    ETFr = L * np.max(ETF_t)
    ADC = np.max(Signal1D(t, m, q, xmin, xmax, ymin, ymax, RC, z)[: len(t)])
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