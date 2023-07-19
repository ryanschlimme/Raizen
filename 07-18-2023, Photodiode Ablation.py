# 18 July 2023
# Ryan Schlimme

# Analyzing photodiode ablation data at 19 J. 500 shots (index 1 to 501). Aggregating pulse data to try to recover pulse. Distance between ablation and active arm: 100 mm. May try reducing to 3-5 mm.

f_name = r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230718\Photodiode_fire\iter_0.tdms" # create a variable pointing to file (change Ryan Schlimme to ryans)

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(r"C:\Users\ryans\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS
from acoustic_entrainment import mic_response						

N = [0]
fig, ax = plt.subplots(1, 1)

M_cutoff = 200e3
L_cutoff = 1e6


for n in N:
# Using low pass filtering AND bin averaging w/ mic correction
	L = CollectionTDMS(f_name)
	L.set_collection("X")
	M = CollectionTDMS(f_name)
	M.set_collection("Y")
	L.apply("detrend", mode = "linear", inplace = True)
	M.apply("detrend", mode = "linear", inplace = True)
	L.apply("calibrate", cal = -1/0.0000004, inplace = True)		# No longer need to invert signal
	M.apply("shift", tau = -92e-6, inplace = True)
	M.apply("lowpass", cutoff = M_cutoff, inplace = True)
	L.apply("lowpass", cutoff = L_cutoff, inplace = True)
	M_Npts = L.r/ (2 * M_cutoff)			# Neiquist Criterion given cutoff frequency
	L_Npts = L.r/ (2 * L_cutoff)
	L.apply("bin_average", Npts = L_Npts, inplace = True)
	M.apply("bin_average", Npts = M_Npts, inplace = True)
	M.apply("correct", response = mic_response, recollect = True) 	# Calibration of microphone
	# M.apply("calibrate", cal = 1/0.00068, inplace = True)
	L.aggrigate(collection_slice = slice(1, 501, 1))
	L.agg.plot(tmin=450e-6, tmax = 500e-6, ax = ax, c = "r", tunit = "us")
	M.aggrigate(collection_slice = slice(1, 501, 1))
	M.agg.plot(tmin=450e-6, tmax = 500e-6, ax = ax, c = "b", tunit = "us")
	for i in list(range(2, 3)):
		Li = L.collection[i]
		Mi = M.collection[i]
		Li.plot(tmin = 440e-6, tmax = 500e-6, ax = ax, c = "crimson", tunit = "us")
		Mi.plot(tmin = 440e-6, tmax = 500e-6, ax = ax, c = "C0", tunit = "us")
	string = "\n" + str(19) + " J"
	ax.set_title(string, fontsize = 9)
plt.suptitle("Acoustic Detection by Photodiode")
plt.show()

 