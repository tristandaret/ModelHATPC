# ModelHATPC

ModelHATPC is a Python project that implements a lineic charge deposit model for
High Angle Time Projection Chambers (HATPCs), derived from the ND280 detector
of the T2K experiment. The work extends previous point-like solutions by
solving the telegrapher's equation for lineic (track-like) deposits and
provides interactive visualizations to inspect per-cell time responses.

Key highlights
- Physics-focused: telegrapher's-equation based modelling for charge propagation.
- Interactive visualization: matplotlib-based GUIs to explore parameter space.
- Educational: demonstrates analytical and numerical techniques for detector modeling.

Quick start
-----------
1. Create a Python environment (Python 3.8+ recommended) and install dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt  # if you create one; otherwise install numpy, matplotlib
```

2. Run one of the interactive scripts (from the repository root):

```bash
python Model/GridPoint.py
python Model/GridTrack.py
```

What to read first
------------------
- `Headers/ModelUtils.py` — contains the physical model functions (Compute0D,
  Compute1D) and parameter definitions.
- `Model/` — top-level interactive scripts that demonstrate the model.
- `Headers/` — shared plotting and geometry utilities.

Repository layout
-----------------
- `Model/` - interactive visualizations (GridPoint, GridTrack)
- `Headers/` - helper modules and shared plotting defaults
- `ClustersEvaluation/`, `LUT/`, `ThesisPlots/` - analysis, LUT generation and
  plotting utilities used during development

LaTeX requirement for plots
-----------------
Some figures are rendered using `matplotlib` with `text.usetex=True`.
To reproduce the plots, install a LaTeX distribution (e.g. TeX Live) and ensure `latex`, `dvipng`, and `ghostscript` are in the PATH.
