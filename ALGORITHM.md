# ALGORITHM SUMMARY — ModelHATPC

**Purpose**
- Model the time response observed on a grid of readout cells produced by a
  lineic (track-like) charge deposit in a Time Projection Chamber (HATPC).
- Provide a small, well-structured codebase that: (1) implements the physical
  model, (2) exposes interactive visualizations for parameter exploration,
  and (3) supports LUT/analysis tooling used during development.

**Core model**
- Governing equation: telegrapher's equation (1D transmission-line-like PDE)
  is used to model signal propagation and dispersion along resistive/capacitive
  detector structures when a continuous (line) deposit is projected onto a
  readout grid.
- The code implements both analytic 0D solutions (point-like approximations)
  and numerical 1D solutions for lineic deposits. Parameters (diffusion,
  resistive/capacitive constants, shaping times) are factored into
  `Headers/ModelUtils.py` and exposed to visualization scripts.

**Numerical approach**
- Spatial discretization: the physical track is projected onto a regular grid
  (see `Headers/GridUtils.py`) which computes per-cell collected charge and
  timing offsets.
- Temporal response: analytical kernels derived from the telegrapher's
  solution are convolved with the deposited line-charge distribution; when an
  analytic convolution is not available the implementation uses stable numeric
  integration (quadrature via `numpy`/`scipy`) with careful sampling to avoid
  aliasing.
- Stability & performance: numerical parameters (time step, integration window
  and sampling density) are chosen adaptively in the demo scripts to trade
  visual interactivity for accuracy. Long-running numerical runs are kept
  separate from the interactive path to preserve responsiveness.

**Software structure**
- `Headers/ModelUtils.py`: core model functions (Compute0D / Compute1D),
  parameter defaults, and kernels. This is the single place to inspect the
  physics implementation.
- `Headers/GridUtils.py` / `Headers/GeometryUtils.py`: mapping from continuous
  tracks to discrete grid cells, geometric transforms, and helper sampling
  routines.
- `Model/GridTrack.py`: interactive demo that composes the model and grid
  utilities, provides UI controls for parameters, and displays per-cell time
  responses using `matplotlib`.
- `LUT/` and `ClustersEvaluation/`: tools to precompute look-up tables and run
  evaluation/analysis pipelines used for tuning and validation.

**How to reason about correctness**
- Sanity checks: compare 0D (point-like) analytical results against numeric
  integrals over small-length line deposits; this verifies kernel normalization
  and time-centering.
- Convergence tests: refine temporal sampling and measure L2 differences in the
  convolved signals; when differences fall below a chosen tolerance the
  numeric configuration is considered stable for visualization.
- Physical validation: where available, compare model outputs to reference
  waveforms or LUT entries produced during development (see `LUT/`).

**Complexity & performance**
- Per-evaluation cost: dominated by convolution/integration cost which is
  O(Nt * Nc) where Nt is temporal samples and Nc is number of contributing
  cells/sources. Optimizations include analytic kernels (when available),
  adaptive sampling, and caching LUTs for repeated parameter combinations.
- Practical limits: the interactive demos are tuned for small grids (tens to
  low hundreds of cells) and interactive frame rates. Large batch runs should
  use dedicated scripts under `LUT/` which can be parallelized.
