# ModelHATPC — Track-based charge propagation model for particle detectors

ModelHATPC is an engineering-grade Python codebase that models lineic (track) charge deposits and their time response in High-Angle Time Projection Chambers (HATPCs), gaseous detectors installed in the near detector of the T2K experiment, studying the physics of neutrinos. The project combines analytical and numerical methods, interactive visualization, and reproducible tooling to demonstrate physics modeling and software engineering best practices.
**Highlights & Capabilities**
- **Applied physics & numerical methods:** Telegrapher's-equation–based modelling for distributed line charges; numerical integration and signal processing.
- **Scientific Python stack:** Practical, production-oriented use of `numpy`, `scipy`, and `matplotlib` for simulation and visualization.
- **Engineering practices:** Reproducible environments (Docker + Makefile), pre-commit hooks, static checking (`mypy`, `flake8`), and clear project layout.
- **Visualization & UX:** Interactive scripts that help explore model parameters and inspect per-cell time responses.

**Architecture**
- **Model:** `Model/GridTrack.py`, `Model/GridPoint.py` — interactive entry points for track and point-like deposit visualizations. Parameters can be adjusted using sliders from `matplotlib.widgets`.
- **Physics & utilities:** `Headers/ModelUtils.py`, `Headers/GeometryUtils.py`, `Headers/GridUtils.py` — core model functions and helpers.
- **LUT & analysis:** `LUT/` and `ClustersEvaluation/` — lookup table tools and analysis scripts used during development.

**Skills demonstrated**
- **Languages & tooling:** Python 3.8+, zsh-friendly scripts, Docker, Make.
- **Numerical analysis:** partial and ordinary differential equation modelling, numerical integration, signal processing and stability considerations.
- **Libraries:** `numpy`, `scipy`, `matplotlib` (interactive GUIs), standard testing and CI tooling.
- **Software engineering:** Modular code layout, docstrings, typed code, pre-commit, static analysis, and containerized reproducibility.
- **DevOps & CI:** Multi-stage Dockerfile, Docker Compose, Make targets, and CI workflow examples for container builds and publishing.

Run the interactive `GridTrack` demo (recommended workflow)

1) Create and activate a virtual environment (macOS / zsh):

```bash
python3 -m venv .venv
source .venv/bin/activate
```

2) Install runtime dependencies:

```bash
pip install -r requirements.txt
```

3) From the repository root, run the GridTrack demo:

```bash
python Model/GridTrack.py
```

Notes:
- The script `Model/GridTrack.py` launches an interactive matplotlib-based visualization. Make sure your environment has access to a display (local desktop or X forwarding) for interactive plots.
- If plots rely on LaTeX rendering (`text.usetex=True`), install a TeX distribution (e.g., TeX Live) or disable LaTeX rendering in script settings.

Quick Docker option (reproducible run)

```bash
# Build the runtime image
make build

# Run a quick smoke test inside the container
docker run --rm modelhatpc:latest python -c "import sys; print('Python', sys.version)"
```
