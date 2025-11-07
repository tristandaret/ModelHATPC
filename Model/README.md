# Model

Contains interactive visualization scripts that demonstrate the HATPC model.

Files
- `GridPoint.py` : Visualize the response to a point-like deposit. Uses
  `Compute0D` from `Headers/ModelUtils.py`.
- `GridTrack.py` : Visualize the response to a lineic / track deposit. Uses
  `Compute1D` from `Headers/ModelUtils.py`.

How to run

From repository root, run:

```bash
python Model/GridPoint.py
python Model/GridTrack.py
```

Both scripts open an interactive matplotlib window with sliders for model
parameters (RC, drift, angle/impact, etc.) and a small map showing the drop or
track position.
