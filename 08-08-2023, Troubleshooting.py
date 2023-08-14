# 08 August 2023
# Ryan Schlimme

# Troubleshooting PSD

Sagnac_name = r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\Sagnac\iter_0.tdms"
BPD_name = r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\SplitBeam\iter_0.tdms"
PD_name = r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\Telescope\iter_0.tdms"

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfinv

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
Npts = [16, 64, 256]

for N in Npts:
	L = CollectionTDMS(Sagnac_name)
	L.set_collection("X")
	local_detrend(L, tmin = 170e-6, tmax = 400e-6, inplace = True)
	L.apply("calibrate", cal = -1/0.00044, inplace = True)
	shot = L.collection[1]
	#Npts = L.r/ (2 * fc)			# Nyquist Criterion given cutoff frequency
	before = shot.x.size
	shot.bin_average(Npts = N, inplace = True)
	after = shot.x.size
	print(before/after)
	print(N)
	shot.plot(tmin = 400e-6, tmax = 469.999e-6, ax = ax)
plt.legend([str(i) for i in cutoffs])
#plt.show()