import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# Use LaTeX font rendering
plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
    }
)

xmin = -0.5
ymax = 3.5

pmns_matrix = np.array(
    [[0.824, 0.546, 0.149], [0.371, 0.598, 0.700], [0.395, 0.573, 0.692]]
)

ckm_matrix = np.array(
    [[0.974, 0.225, 0.0037], [0.225, 0.973, 0.0418], [0.0086, 0.041, 0.999]]
)


# Normalize and get colormap
cmapPMNS = plt.get_cmap("Blues")
cmapCKM = plt.get_cmap("Reds")

fig, ax = plt.subplots()
ax.set_aspect("equal")
ax.set_xlim(xmin, 3)
ax.set_ylim(0, ymax)

# Draw the PMNS matrix
for i in range(3):
    for j in range(3):
        val = np.sqrt(pmns_matrix[i, j])
        color_val = cmapPMNS(val)
        valedge = min(val + 0.2, 0.9999)
        if valedge == 0.9999:
            color_edge = "black"
        else:
            color_edge = cmapPMNS(valedge)
        center_x = j - xmin
        center_y = 2.5 - i
        rounding_size = val / 10
        ax.add_patch(
            FancyBboxPatch(
                (center_x - val / 2, center_y - val / 2),
                val,
                val,
                boxstyle=f"round,pad=0,rounding_size={rounding_size}",
                facecolor=color_val,
                edgecolor=color_edge,
            )
        )

# Row labels (left side)
row_labels = [r"$\nu_e$", r"$\nu_\mu$", r"$\nu_\tau$"]
for i, label in enumerate(row_labels):
    ax.text(-0.15, 3 - i - 0.5, label, va="center", ha="right", fontsize=30)

# Column labels (top)
col_labels = [r"$\nu_1$", r"$\nu_2$", r"$\nu_3$"]
for j, label in enumerate(col_labels):
    ax.text(j + 0.5, 3.15, label, ha="center", va="bottom", fontsize=30)

ax.axis("off")
plt.savefig("PMNS_matrix.pdf", format="pdf", bbox_inches="tight")

# Draw the CKM matrix
fig, ax = plt.subplots()
ax.set_aspect("equal")
ax.set_xlim(xmin, 3)
ax.set_ylim(0, ymax)
# Draw the CKM matrix
for i in range(3):
    for j in range(3):
        val = np.sqrt(ckm_matrix[i, j])
        color_val = cmapCKM(val)
        valedge = min(val + 0.2, 0.9999)
        if valedge == 0.9999:
            color_edge = "black"
        else:
            color_edge = cmapCKM(valedge)
        center_x = j - xmin
        center_y = 2.5 - i
        rounding_size = val / 10
        ax.add_patch(
            FancyBboxPatch(
                (center_x - val / 2, center_y - val / 2),
                val,
                val,
                boxstyle=f"round,pad=0,rounding_size={rounding_size}",
                facecolor=color_val,
                edgecolor=color_edge,
            )
        )
# Row labels (left side)
row_labels = [r"$u$", r"$c$", r"$t$"]
for i, label in enumerate(row_labels):
    ax.text(-0.15, 3 - i - 0.5, label, va="center", ha="right", fontsize=30)
# Column labels (top)
col_labels = [r"$d$", r"$s$", r"$b$"]
for j, label in enumerate(col_labels):
    ax.text(j + 0.5, 3.15, label, ha="center", va="bottom", fontsize=30)
ax.axis("off")
plt.savefig("CKM_matrix.pdf", format="pdf", bbox_inches="tight")
