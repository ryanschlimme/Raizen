# 08 August 2023
# Ryan Schlimme

# Generate averaged PSD for all systems: microphone, single photodiode, balanced photodetection, and Sagnac interferometer.

#Sagnac_name = r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\Sagnac\iter_0.tdms"
#BPD_name = r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\SplitBeam\iter_0.tdms"
#PD_name = r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\Telescope\iter_0.tdms"

I = [145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160]

Sagnac_name_index = [r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\Sagnac\MinDetect\phi" + str(i) + ".tdms" for i in I]
BPD_name_index = [r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\SplitBeam\MinDetect\phi" + str(i) + ".tdms" for i in I]
PD_name_index = [r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\Telescope\MinDetect\phi" + str(i) + ".tdms" for i in I]

Sagnac_name = Sagnac_name_index[0]
BPD_name = BPD_name_index[0]
PD_name = PD_name_index[0]

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfinv
cmap = plt.colormaps('tab20')

sys.path.append(r"C:\Users\ryans\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS, PSD
from acoustic_entrainment import mic_response	


def local_detrend(col, tmin = None, tmax = None, inplace = False) -> None:
	for c in col.collection:
		t, x = c.time_gate(tmin = tmin, tmax = tmax)
		x_bar = np.mean(x)
		m, b = np.polyfit(t, x, 1)
		if inplace:
			c.x = c.x - (m * c.t) - b
	return None


fig, ax = plt.subplots(1,1, sharex = True)
cutoff = 3e6

# Sagnac
L = CollectionTDMS(Sagnac_name)
L.set_collection("X")
M = CollectionTDMS(Sagnac_name)
M.set_collection("Y")
local_detrend(L, tmin = 170e-6, tmax = 400e-6, inplace = True)
local_detrend(M, tmin = 170e-6, tmax = 400e-6, inplace = True)
#L.apply("detrend", mode = "linear", inplace = True)
#M.apply("detrend", mode = "linear", inplace = True)
#L.apply("calibrate", cal = -1/0.00044, inplace = True)
#M.apply("shift", tau = -19.5e-6, inplace = True)
#M.apply("lowpass", cutoff = 1e6, inplace = True)
#L.apply("lowpass", cutoff = 1e6, inplace = True)
Npts = int(L.r/ (2 * cutoff))			# Nyquist Criterion given cutoff frequency
L.apply("bin_average", Npts = Npts, inplace = True)
M.apply("bin_average", Npts = Npts, inplace = True)
# M.apply("correct", response = mic_response, recollect = True) 
# M.apply("calibrate", cal = 1/0.00068, inplace = True)

Sagnacf, SagnacPSD = L.average("PSD", tmin = 437e-6, tmax = 457e-6, detrend = False)
Micf, MicPSD = M.average("PSD", tmin = 414e-6, tmax = 448e-6, detrend = False)


# BPD
L = CollectionTDMS(BPD_name)
L.set_collection("X")
local_detrend(L, tmin = 170e-6, tmax = 400e-6, inplace = True)
#L.apply("detrend", mode = "linear", inplace = True)
#L.apply("calibrate", cal = 1/0.0002, inplace = True)
#L.apply("lowpass", cutoff = 1e6, inplace = True)
L.apply("bin_average", Npts = Npts, inplace = True)

BPDf, BPDPSD = L.average("PSD", tmin = 437e-6, tmax = 457e-6, detrend = False)


# PD
L = CollectionTDMS(PD_name)
L.set_collection("X")
local_detrend(L, tmin = 170e-6, tmax = 400e-6, inplace = True)
#L.apply("detrend", mode = "linear", inplace = True)
#L.apply("calibrate", cal = -1/0.000001, inplace = True)
#L.apply("lowpass", cutoff = 1e6, inplace = True)
L.apply("bin_average", Npts = Npts, inplace = True)

PDf, PDPSD = L.average("PSD", tmin = 437e-6, tmax = 457e-6, detrend = False)

# Plotting
ax.loglog(Micf, MicPSD, color = cmap[0])
ax.loglog(Micf, MicPSDnoise, color = cmap[1]) 	# odd numbers are faded
ax.loglog(Sagnacf, SagnacPSD)
ax.loglog(BPDf, BPDPSD)
ax.loglog(PDf, PDPSD)
ax.legend(["Mic", "Sagnac", "BPD", "PD"])
plt.xlabel("Frequency (Hz)")
plt.title("PSD of Four Methods")

plt.show()