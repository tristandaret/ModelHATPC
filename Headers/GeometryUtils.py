import numpy as np

##### GEOMETRY SETTINGS ######


### GEOMETRY PARAMETERS ###
nX = 3  # number of colums
xwidth = 11.28  # mm ; width of xmin pad
xc = xwidth / 2  # mm ; horizontal center of pad
xleft = (
    -(nX // 2) * xwidth
)  # Most left position of the superpad considered (wrt left of LP)
xright = (
    nX // 2 + 1
) * xwidth  # Most right position of the superpad considered (wrt left of LP)
xmin = 0  # mm ; bottom border of leading pad
xmax = xwidth  # mm ; top border of leading pad

nY = 3  # number of rows
ywidth = 10.19  # mm ; height of xmin pad
yc = ywidth / 2  # mm ; vertical center of pad
ylow = (
    -(nY // 2) * ywidth
)  # Lowest position of the superpad considered (wrt bottom of LP)
yhigh = (
    nY // 2 + 1
) * ywidth  # Highest position of the superpad considered (wrt bottom of LP)
ymin = 0  # mm ; left border of leading pad
ymax = ywidth  # mm ; right border of leading pad

diag = np.sqrt(xwidth**2 + ywidth**2)


#### GEOMETRY FUNCTIONS ####


def X(phi_rad, d, y):
    if phi_rad == 0:
        return -np.inf
    return (
        y - (d - np.sin(phi_rad) * xc + np.cos(phi_rad) * yc) / np.cos(phi_rad)
    ) / np.tan(phi_rad)


def Y(phi_rad, d, x):
    if phi_rad == 90:
        return -np.inf
    return np.tan(phi_rad) * x + (
        d - np.sin(phi_rad) * xc + np.cos(phi_rad) * yc
    ) / np.cos(phi_rad)


# Compute parameters of track in cluster
def compute_line_params(phi_deg, d):
    phi_rad = phi_deg / 180 * np.pi
    m = np.tan(phi_rad)
    q = (np.cos(phi_rad) * yc - np.sin(phi_rad) * xc + d) / np.cos(phi_rad)
    return m, q, phi_rad


# Cluster lengths
def ClusterLengths(phi_rad, d):
    out_vert = 0
    out_diag = 0
    out_cros = 0
    r_diag = 0
    r_vert = 0
    r_cros = 0
    L = 0

    y_xmin = Y(phi_rad, d, xmin)
    y_xmax = Y(phi_rad, d, xmax)
    y_xleft = Y(phi_rad, d, xleft)
    y_xright = Y(phi_rad, d, xright)
    x_ymin = X(phi_rad, d, ymin)
    x_ymax = X(phi_rad, d, ymax)
    x_ylow = X(phi_rad, d, ylow)
    x_yhigh = X(phi_rad, d, yhigh)

    # Whole map
    # In
    if ylow > y_xleft:
        (x0_eram, y0_eram) = (x_ylow, ylow)
    else:
        (x0_eram, y0_eram) = (xleft, y_xleft)
    # Out
    if yhigh > y_xright:
        (x1_eram, y1_eram) = (xright, y_xright)
    else:
        (x1_eram, y1_eram) = (x_yhigh, yhigh)

    # Leading Pad & Diagonal
    # In
    if ymax >= y_xmin > ymin:
        (x0_diag, y0_diag) = (xmin, y_xmin)
    elif y_xmax >= ymin > y_xmin:
        (x0_diag, y0_diag) = (x_ymin, ymin)
    else:
        out_diag = 1
    # Out
    if ymax > y_xmax:
        (x1_diag, y1_diag) = (xmax, y_xmax)
    else:
        (x1_diag, y1_diag) = (x_ymax, ymax)

    # Vertical
    # left-right
    if ylow < y_xmin < yhigh and ylow < y_xmax < yhigh:
        (x0_vert, y0_vert) = (xmin, y_xmin)
        (x1_vert, y1_vert) = (xmax, y_xmax)
    # left-top
    elif ylow < y_xmin < yhigh and yhigh < y_xmax:
        (x0_vert, y0_vert) = (xmin, y_xmin)
        (x1_vert, y1_vert) = (x_yhigh, yhigh)
    # left-bottom
    elif ylow < y_xmin < yhigh and y_xmax < ylow:
        (x0_vert, y0_vert) = (xmin, y_xmin)
        (x1_vert, y1_vert) = (x_ylow, ylow)
    # bottom-top
    elif y_xmin < ylow and yhigh < y_xmax:
        (x0_vert, y0_vert) = (x_ylow, ylow)
        (x1_vert, y1_vert) = (x_yhigh, yhigh)
    # bottom-right
    elif y_xmin < ylow and y_xmax < yhigh:
        (x0_vert, y0_vert) = (x_ylow, ylow)
        (x1_vert, y1_vert) = (xright, y_xright)
    else:
        out_vert = 1

    # Half-cross
    # In
    if y_xmin > yhigh:
        out_cros = 1  # pass above
    elif yhigh > y_xmin >= ymin:
        (x0_cros, y0_cros) = (xmin, y_xmin)  # enter from the left
    elif y_xright > ymin >= y_xmin:
        (x0_cros, y0_cros) = (x_ymin, ymin)  # enter from below
    elif ymin > y_xright:
        out_cros = 1  # pass below
    # Out
    if y_xmax > yhigh:
        (x1_cros, y1_cros) = (x_yhigh, yhigh)  # exit on the top
    elif yhigh > y_xmax >= ymax:
        (x1_cros, y1_cros) = (xmax, y_xmax)  # exit on the top right part
    elif y_xright >= ymax > y_xmax:
        (x1_cros, y1_cros) = (x_ymax, ymax)  # exit on the right top part
    elif ymax > y_xright:
        (x1_cros, y1_cros) = (xright, y_xright)  # exit on the right

    # Compute lengths
    if out_vert == 0:
        r_vert = np.sqrt((y1_vert - y0_vert) ** 2 + (x1_vert - x0_vert) ** 2)
    if out_diag == 0:
        r_diag = np.sqrt((y1_diag - y0_diag) ** 2 + (x1_diag - x0_diag) ** 2)
    if out_cros == 0:
        r_cros = np.sqrt((y1_cros - y0_cros) ** 2 + (x1_cros - x0_cros) ** 2)

    # Compute the length of the whole map
    L = np.sqrt((y1_eram - y0_eram) ** 2 + (x1_eram - x0_eram) ** 2)

    return (r_diag, r_vert, r_cros, L)
