# 24 July 2023
# Ryan Schlimme

# Analyzing inidvidual shots of Sagnac with fire rate of 1 Hz and 10 Hz overlayed to qualitatively determine any dependence on amplitude.

f_name_1 = r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\20230721\Sagnac_19J_1hz\iter_0.tdms" # create a variable pointing to file (change Ryan Schlimme to ryans)

f_name_10 = r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\20230721\Sagnac_19J_10hz\iter_0.tdms" # create a variable pointing to file (change Ryan Schlimme to ryans)

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS
from acoustic_entrainment import mic_response						

N = [0]
fig, ax = plt.subplots(1, 1)

M_cutoff = 200e3
L_cutoff = 1e6


for n in N:
# Using low pass filtering AND bin averaging w/ mic correction
	L = CollectionTDMS(f_name_1)
	L.set_collection("X")
	M = CollectionTDMS(f_name_1)
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
	#L.aggrigate(collection_slice = slice(1, 501, 1))
	#L.agg.plot(tmin=450e-6, tmax = 500e-6, ax = ax, c = "r", tunit = "us")
	#M.aggrigate(collection_slice = slice(1, 501, 1))
	#M.agg.plot(tmin=450e-6, tmax = 500e-6, ax = ax, c = "b", tunit = "us")
	for i in list(range(3, 4)):
		Li = L.collection[i]
		Mi = M.collection[i]
		Li.plot(tmin = 440e-6, tmax = 500e-6, ax = ax, c = "r", tunit = "us")
		#Mi.plot(tmin = 440e-6, tmax = 500e-6, ax = ax, c = "b", tunit = "us")
	#string = "\n" + str(19) + " J"
	#ax.set_title(string, fontsize = 9)
	L = CollectionTDMS(f_name_10)
	L.set_collection("X")
	M = CollectionTDMS(f_name_10)
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
	#L.aggrigate(collection_slice = slice(1, 501, 1))
	#L.agg.plot(tmin=450e-6, tmax = 500e-6, ax = ax, c = "r", tunit = "us")
	#M.aggrigate(collection_slice = slice(1, 501, 1))
	#M.agg.plot(tmin=450e-6, tmax = 500e-6, ax = ax, c = "b", tunit = "us")
	for i in list(range(3, 4)):
		Li = L.collection[i]
		Mi = M.collection[i]
		Li.plot(tmin = 440e-6, tmax = 500e-6, ax = ax, c = "b", tunit = "us")
		#Mi.plot(tmin = 440e-6, tmax = 500e-6, ax = ax, c = "b", tunit = "us")
	#string = "\n" + str(19) + " J"
	#ax.set_title(string, fontsize = 9)
plt.suptitle("Acoustic Detection by Sagnac Interferometer")
plt.title("Comparing Fire Rates", fontsize = 9)
plt.legend(["1 Hz", "10 Hz"])
plt.show()

 