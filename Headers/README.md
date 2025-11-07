# Headers

Shared helper modules used by the ModelHATPC visualizations and utilities.

Important modules
- `ModelUtils.py` : Physics and model evaluation routines (Compute0D, Compute1D).
- `GridUtils.py`  : Shared plotting defaults (figure creation, time axis, colors).
- `GeometryUtils.py` : Geometric helpers for tracks and grid coordinates.
- `ClustersUtils.py` : Utilities to setup cluster/track interactive widgets.

These modules are lightweight and are intended to be imported directly by
scripts in `Model/` and `ClustersEvaluation/`.
