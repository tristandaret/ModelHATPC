"""Headers package for the ModelHATPC project.

This package contains shared utilities and helper modules used by the
visualization and modelling scripts. Typical modules include:

- ModelUtils.py : physical model and compute functions (Compute0D, Compute1D)
- GridUtils.py  : shared plotting defaults and UI elements
- GeometryUtils.py : geometric helpers for tracks and grid layout
- ClustersUtils.py : interactive cluster/track visualization widgets

Importing this package makes it clearer that these helper modules live in a
single namespace; modules are still imported directly by name in scripts for
backwards compatibility.
"""

__all__ = [
    "ModelUtils",
    "GridUtils",
    "GeometryUtils",
    "ClustersUtils",
]

# Ensure the names listed in __all__ are present in the package namespace.
# Importing the submodules makes them available as attributes of the package.
from . import ClustersUtils, GeometryUtils, GridUtils, ModelUtils
