from sys import path

path.append("Headers/")
from ModelUtils import *
import time
from ROOT import TFile, TTree
from array import array

start_time = time.time()
# Heatmaps --------------------------------------------------------------------------------------------
t = np.linspace(1, 1000, 250)  # ns | start at 1 to avoid sigma = 0
# Computations
ETF = lambdaG * ETF(t)
nphi = 250
nd = 250
nZ = 101

arr_r = np.full((nd, nphi), np.nan)
v_Dt = [310, 350]
v_RC = [112, 158]
v_d = np.linspace(0, diag / 2, nd)
v_phi = np.linspace(0, 90, nphi)
v_Z = np.linspace(0, 1000, nZ)

out_file = TFile("LUT/dEdx_XP_LUT_test.root", "RECREATE")
out_tree = TTree("outTree", "LUT")
Dt_array = array("f", [0])
RC_array = array("f", [0])
phi_array = array("f", [0])
d_array = array("f", [0])
z_array = array("f", [0])
weight_array = array("f", [0])
out_tree.Branch("transDiff", Dt_array, "transDiff/F")
out_tree.Branch("RC", RC_array, "RC/F")
out_tree.Branch("angle", phi_array, "angle/F")
out_tree.Branch("impact_param", d_array, "impact_param/F")
out_tree.Branch("drift_dist", z_array, "drift_dist/F")
out_tree.Branch("weight", weight_array, "weight/F")


# Make Length map
phi_index = 0
for phi in v_phi:
    phi_rad = phi / 180 * np.pi
    d_index = 0
    for d in v_d:

        # Determine the length of the track across the central pad
        x = []
        y = []

        y_xmin = Y(phi_rad, d, xmin)
        y_xmax = Y(phi_rad, d, xmax)
        x_ymin = X(phi_rad, d, ymin)
        x_ymax = X(phi_rad, d, ymax)

        if ymin <= y_xmin < ymax:
            x.append(xmin)
            y.append(y_xmin)

        if ymin <= y_xmax < ymax:
            x.append(xmax)
            y.append(y_xmax)

        if xmin <= x_ymin < xmax:
            x.append(x_ymin)
            y.append(ymin)

        if xmin <= x_ymax < xmax:
            x.append(x_ymax)
            y.append(ymax)

        L = 0
        if len(x) == 2:
            L = np.sqrt((y[1] - y[0]) ** 2 + (x[1] - x[0]) ** 2)
        arr_r[d_index, phi_index] = L
        d_index += 1
    phi_index += 1
print(f"Length map done in {time.time()-start_time:.1f} seconds")


# Make LUT for all (RC, z, phi, d)
LUT_time = time.time()
for Dt in v_Dt:
    print(f"Dt = {Dt:.2e} sqrt(mm)/cm")
    Dt_array[0] = Dt
    Dt_time = time.time()

    for RC in v_RC:
        print(f"RC = {RC:.0f} ns/mm²")
        RC_array[0] = RC
        RC_time = time.time()

        for z in v_Z:
            print(f"z = {z:.0f} mm")
            phi_index = 0
            z_array[0] = z

            for phi in v_phi:
                print(f"phi = {phi:.1f}°")
                d_index = 0
                phi_array[0] = phi
                if phi == 0:
                    phi = 1e-5
                if phi == 90:
                    phi = 90 - 1e-5
                phi_rad = phi / 180 * np.pi

                for d in v_d:
                    # print(f"d = {d:.2f} mm")
                    d_array[0] = d
                    m = np.tan(phi_rad)  # slope
                    q = (np.cos(phi_rad) * yc - np.sin(phi_rad) * xc + d) / np.cos(
                        phi_rad
                    )  # intercept

                    ETFr = arr_r[d_index, phi_index] * np.max(ETF)

                    ADC = np.max(
                        Signal1D(
                            t,
                            m,
                            q,
                            xmin,
                            xmax,
                            ymin,
                            ymax,
                            RC,
                            z,
                            Dt / np.power(10, 7 / 2),
                        )[: len(t)]
                    )
                    weight_array[0] = ETFr / ADC
                    if not (np.isnan(weight_array[0])) and weight_array[0] > 0:
                        out_tree.Fill()

                    d_index += 1
            phi_index += 1
        print(f"RC done in {time.time()-RC_time:.1f} seconds")
    print(f"Dt done in {time.time()-Dt_time:.1f} seconds")
print(f"LUT done in {time.time()-LUT_time:.1f} seconds")

print(f"Total time: {time.time()-start_time:.1f} seconds")

out_tree.Write()
out_file.Close()
