"""Grid plotting defaults and shared UI elements.

This module centralizes plotting defaults and shared matplotlib objects
used across the visualization scripts. It defines the time axis, plotting
colors and units depending on `vartype` (Signal, Charge or Current), creates
the main `fig, axs` grid of subplots and provides a `Reset` button used by
interactive figures.

Exports
-------
- fig, axs : matplotlib.figure.Figure and Axes array
- t, timescale, timeunit : time axis and units
- varcolor, varforeground, varxminplot, varxmaxplot, varyminplot, varymaxplot
- scalefactor, dim, unit, ylabel : plotting labels and scaling
- button : matplotlib Button instance to reset sliders

The plotting scripts import this module and then build the rest of the UI
around these shared objects.
"""

from sys import path

path.append("Headers/")
import warnings
from typing import Any, List

import GeometryUtils as geo
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Button

# ========================= Variable Setup =========================
# Variable types: Signal, Charge, Current
vartype = "Signal"

if vartype == "Signal":
    timescale = 1
    timeunit = "ns"
    t = np.linspace(1, 3000, 500)
    varcolor = "C1"
    varforeground = "black"
    varxminplot = 0
    varxmaxplot = 3000
    varyminplot = -200
    varymaxplot = 1000
    scalefactor = 1
    dim = "ADC$_{max}$"
    unit = ""
    ylabel = "ADC count"
elif vartype == "Charge":
    tmax = int(5e3)  # Maximum time in ns
    timescale = 1000  # Timescale factor to get in µs
    timeunit = "µs"
    ntsteps = int(5e3)
    t = np.linspace(1, tmax, ntsteps)
    varcolor = "C3"
    varforeground = "white"
    varxminplot = 0
    varxmaxplot = tmax
    varyminplot = -3
    varymaxplot = 30
    scalefactor = 1
    dim = "Q$_{max}$"
    unit = "fC"
    ylabel = "Charge (fC)"
elif vartype == "Current":
    timescale = 1
    timeunit = "ns"
    tmax = int(3000)
    ntsteps = int(3000)
    t = np.linspace(1, tmax, ntsteps)
    varcolor = "C0"
    varforeground = "white"
    varxminplot = 0
    varxmaxplot = tmax
    varyminplot = -20
    varymaxplot = 20
    scalefactor = int(1e3)  # Scale factor to convert from µA to nA
    dim = "I$_{max}$"
    unit = "nA"
    ylabel = "Current (nA)"
else:
    # Fallback: ensure all plotting variables are defined and notify the user.
    warnings.warn(
        f"Unknown vartype '{vartype}', using default 'Signal' settings.",
        UserWarning,
        stacklevel=2,
    )
    timescale = 1
    timeunit = "ns"
    t = np.linspace(1, 3000, 500)
    varcolor = "C1"
    varforeground = "black"
    varxminplot = 0
    varxmaxplot = 3000
    varyminplot = -200
    varymaxplot = 1000
    scalefactor = 1
    dim = "ADC$_{max}$"
    unit = ""
    ylabel = "ADC count"
    ylabel = "Current (nA)"

# ========================= Plot Setup =========================
fig, axs = plt.subplots(geo.nY, geo.nX, sharex=True, sharey=True, figsize=(12, 9))
fig.subplots_adjust(left=0.2, right=0.98, bottom=0.07, top=0.98)
lines: List[Any] = []
texts: List[Any] = []
# Number of ticks on x-axis
nxticks = int(t.max() // 1000) + 1  # Number of ticks on x-axis
for ax in axs.flat:
    ax.set_xticks(np.linspace(varxminplot, varxmaxplot / timescale, int(nxticks)))


# ========================= Reset Button Setup =========================
resetax = fig.add_axes((0.035, 0.05, 0.085, 0.06))
button = Button(resetax, "Reset", hovercolor="0.975")
button.label.set_fontsize(20)
