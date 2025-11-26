"""Simple ratio estimator utilities used during LUT development.

Provides a function `ratio` that computes the ratio between an estimated
energy deposition (from geometry, `L`) and the maximal ADC amplitude from
the model. This helper was used to validate and debug LUT entries.
"""

import numpy as np

from Headers import GeometryUtils as geo
from Headers import ModelUtils as mu

# compute ETF(t) once and store in a distinct variable to avoid shadowing
ETF_t = mu.lambdaG * mu.ETF(mu.t)


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
    q = (d - geo.xc * np.sin(phi_rad) + geo.yc * np.cos(phi_rad)) / np.cos(phi_rad)
    ETFr = L * np.max(ETF_t)
    ADC = np.max(
        mu.Signal1D(mu.t, m, q, geo.xmin, geo.xmax, geo.ymin, geo.ymax, RC, z)[
            : len(mu.t)
        ]
    )
    ADC = np.max(
        mu.Signal1D(mu.t, m, q, geo.xmin, geo.xmax, geo.ymin, geo.ymax, RC, z)[
            : len(mu.t)
        ]
    )

    # Determine intersection points of the line with the central pad
    x = []
    y = []
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
