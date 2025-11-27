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

Code hygiene & docstrings
-------------------------
This repository follows strict code-hygiene practices to keep the codebase
readable, consistent and recruiter-friendly:

- **Pre-commit hooks:** A `.pre-commit-config.yaml` is provided and configures
  formatting and checks (Black, isort, flake8 + bugbear, mypy and common
  pre-commit hooks) so commits are linted and formatted automatically.
- **Developer deps:** Development tools are listed in `requirements-dev.txt`
  and the project `pyproject.toml` includes an optional `dev` extras section
  to make onboarding and CI straightforward.
- **Runtime deps:** Runtime dependencies (e.g. `numpy`, `scipy`, `matplotlib`)
  are declared in `requirements.txt` and `pyproject.toml` metadata.
- **Docstrings:** Code in `Headers/`, `Model/` and other modules follows a
  consistent docstring convention (PEP 257 / Google-style) so functions,
  classes and modules are documented and easy to understand.

How to run the checks locally
-----------------------------
1. Install development tools: `pip install -r requirements-dev.txt`.
2. Install and activate pre-commit hooks: `pre-commit install`.
3. Run the linters/formatters manually if desired: `pre-commit run --all-files`.

These measures ensure the repository is well-formatted, statically checked,
and documented — a signal of good engineering practices for reviewers and
recruiters.
