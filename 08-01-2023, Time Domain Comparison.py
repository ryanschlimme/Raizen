# 01 August 2023
# Ryan Schlimme

# Comparing time domain laser pulses to calibrated microphone reading at 19 J. Applying known calibration to microphone to transfer to pressure signal. All in Sagnac interferometer with telescope and balanced photodetection. Attempt 2 to correct for potential misfire

Sagnac_name = r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\20230801\Sagnac\iter_0.tdms" # create a variable pointing to file (change Ryan Schlimme to ryans)
BFD_name = r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\20230801\SplitBeam\iter_0.tdms"
PD_name = r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\20230801\Telescope\iter_0.tdms"

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS
from acoustic_entrainment import mic_response						

N = [0]
fig, ax = plt.subplots(1,1, sharey = True, sharex = True)

M_cutoff = 200e3
L_cutoff = 1e6


for n in N:
# Using low pass filtering AND bin averaging w/ mic correction
	L = CollectionTDMS(Sagnac_name)
	L.set_collection("X")
	M = CollectionTDMS(Sagnac_name)
	M.set_collection("Y")
	L.apply("detrend", mode = "linear", inplace = True)
	M.apply("detrend", mode = "linear", inplace = True)
	L.apply("calibrate", cal = -1/0.00044, inplace = True)
	M.apply("shift", tau = -19.5e-6, inplace = True)
	M.apply("lowpass", cutoff = M_cutoff, inplace = True)
	L.apply("lowpass", cutoff = L_cutoff, inplace = True)
	M_Npts = L.r/ (2 * M_cutoff)			# Neiquist Criterion given cutoff frequency
	L_Npts = L.r/ (2 * L_cutoff)
	L.apply("bin_average", Npts = L_Npts, inplace = True)
	M.apply("bin_average", Npts = M_Npts, inplace = True)
	M.apply("correct", response = mic_response, recollect = True) 	# Calibration of microphone
	# M.apply("calibrate", cal = 1/0.00068, inplace = True)
	L.aggrigate(collection_slice = slice(2, 500, 1))
	L.agg.plot(tmin=420e-6, tmax = 470e-6, ax = ax, tunit = "us")
	M.aggrigate(collection_slice = slice(2, 500, 1))
	M.agg.plot(tmin=420e-6, tmax = 470e-6, ax = ax, c = "b", tunit = "us")
	for i in list(range(2, 3)):
		Li = L.collection[i]
		Mi = M.collection[i]
		#Li.plot(tmin = 420e-6, tmax = 470e-6, ax = ax, c = "crimson", tunit = "us")
		#Mi.plot(tmin = 420e-6, tmax = 470e-6, ax = ax, c = "C0", tunit = "us")

for n in N:
# Using low pass filtering AND bin averaging w/ mic correction
	L = CollectionTDMS(BFD_name)
	L.set_collection("X")
	L.apply("detrend", mode = "linear", inplace = True)
	L.apply("calibrate", cal = 1/0.0002, inplace = True)
	L.apply("lowpass", cutoff = L_cutoff, inplace = True)
	L_Npts = L.r/ (2 * L_cutoff)
	L.apply("bin_average", Npts = L_Npts, inplace = True)
	L.aggrigate(collection_slice = slice(2, 500, 1))
	L.agg.plot(tmin=420e-6, tmax = 470e-6, ax = ax, tunit = "us")
	for i in list(range(2, 3)):
		Li = L.collection[i]
		Mi = M.collection[i]
		#Li.plot(tmin = 420e-6, tmax = 470e-6, ax = ax, c = "crimson", tunit = "us")
		#Mi.plot(tmin = 420e-6, tmax = 470e-6, ax = ax, c = "C0", tunit = "us")

for n in N:
# Using low pass filtering AND bin averaging w/ mic correction
	L = CollectionTDMS(PD_name)
	L.set_collection("X")
	L.apply("detrend", mode = "linear", inplace = True)
	L.apply("calibrate", cal = -1/0.0000002, inplace = True)	
	L.apply("lowpass", cutoff = L_cutoff, inplace = True)
	L_Npts = L.r/ (2 * L_cutoff)
	L.apply("bin_average", Npts = L_Npts, inplace = True)
	L.aggrigate(collection_slice = slice(2, 500, 1))
	L.agg.plot(tmin=420e-6, tmax = 470e-6, ax = ax, tunit = "us")
	for i in list(range(2, 3)):
		Li = L.collection[i]
		Mi = M.collection[i]
		#Li.plot(tmin = 420e-6, tmax = 470e-6, ax = ax, c = "crimson", tunit = "us")
		#Mi.plot(tmin = 420e-6, tmax = 470e-6, ax = ax, c = "C0", tunit = "us")


string = "\n" + str(19) + " J"
plt.title(string, fontsize = 9)
plt.suptitle("Acoustic Detection by Sagnac Interferometer")
plt.legend(["Sagnac", "BFD", "PD", "Tele"])
plt.show()

