# 10 July 2023
# Ryan Schlimme

# Comparing time domain laser pulses to calibrated microphone reading at energies scaling from 12 to 19 J in 1 J increments. Applying known calibration to microphone to transfer to pressure signal. All in Sagnac interferometer with telescope and balanced photodetection.

f_name_index = [r"C:\Users\ryans\OneDrive\Desktop\Research\Data\Sagnac_tele_ene_scan_501shots\iter_" + str(i) + ".tdms" for i in range(8)] # create a variable pointing to file (change Ryan Schlimme to ryans)

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(r"C:\Users\ryans\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS
from acoustic_entrainment import mic_response						

N = list(range(8))
fig, axes = plt.subplots(4,2, sharey = True, sharex = True)

M_cutoff = 200e3
L_cutoff = 1e6


for n, ax in zip(N, axes.flatten()):
# Using low pass filtering AND bin averaging w/ mic correction
	f_name = f_name_index[n]
	L = CollectionTDMS(f_name)
	L.set_collection("X")
	M = CollectionTDMS(f_name)
	M.set_collection("Y")
	L.apply("detrend", mode = "linear", inplace = True)
	M.apply("detrend", mode = "linear", inplace = True)
	L.apply("calibrate", cal = -1/0.001, inplace = True)		# No longer need to invert signal
	M.apply("shift", tau = -90e-6, inplace = True)
	M.apply("lowpass", cutoff = M_cutoff, inplace = True)
	L.apply("lowpass", cutoff = L_cutoff, inplace = True)
	M_Npts = L.r/ (2 * M_cutoff)			# Neiquist Criterion given cutoff frequency
	L_Npts = L.r/ (2 * L_cutoff)
	L.apply("bin_average", Npts = L_Npts, inplace = True)
	M.apply("bin_average", Npts = M_Npts, inplace = True)
	M.apply("correct", response = mic_response, recollect = True) 	# Calibration of microphone
	# M.apply("calibrate", cal = 1/0.00068, inplace = True)
	L.aggrigate(collection_slice = slice(2, 500, 1))
	L.agg.plot(tmin=450e-6, tmax = 490e-6, ax = ax, c = "r", tunit = "us")
	M.aggrigate(collection_slice = slice(2, 500, 1))
	M.agg.plot(tmin=450e-6, tmax = 490e-6, ax = ax, c = "b", tunit = "us")
	for i in list(range(2, 3)):
		Li = L.collection[i]
		Mi = M.collection[i]
		Li.plot(tmin = 450e-6, tmax = 490e-6, ax = ax, c = "crimson", tunit = "us")
		Mi.plot(tmin = 450e-6, tmax = 490e-6, ax = ax, c = "C0", tunit = "us")
	string = "\n" + str(n + 12) + " J"
	ax.set_title(string, fontsize = 9)
plt.suptitle("Acoustic Detection by Sagnac Interferometer")
plt.show()

