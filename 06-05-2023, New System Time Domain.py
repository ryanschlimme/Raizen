# 05 June 2023
# Ryan Schlimme

# Comparing time domain laser pulses to calibrated microphone reading at energies scaling from 12 to 19 J in 1 J increments. Applying known calibration to microphone to transfer to pressure signal. All in newly built setup of larger detection length.

f_name_index = [r"C:\Users\ryans\OneDrive\Desktop\Research\Data\laserX_microphoneY_energy12-19_N8\iter_" + str(i) + ".tdms" for i in range(8)] # create a variable pointing to file (change Ryan Schlimme to ryans)

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(r"C:\Users\ryans\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS
from acoustic_entrainment import mic_response						

N = list(range(8))
fig, axes = plt.subplots(8,1)

for n, ax in zip(N, axes.flatten()):
# Using low pass filtering AND bin averaging w/ mic correction
	f_name = f_name_index[n]
	L = CollectionTDMS(f_name)
	L.set_collection("X")
	M = CollectionTDMS(f_name)
	M.set_collection("Y")
	L.apply("detrend", mode = "linear", inplace = True)
	M.apply("detrend", mode = "linear", inplace = True)
	L.apply("calibrate", cal = -1, inplace = True)
	M.apply("shift", tau = 35e-6, inplace = True)
	M.apply("lowpass", cutoff = 2.5e6, inplace = True)
	L.apply("lowpass", cutoff = 2.5e6, inplace = True)
	Npts = L.r/ (2 * 2.5e6)			# Neiquist Criterion given cutoff frequency
	L.apply("bin_average", Npts = Npts, inplace = True)
	M.apply("bin_average", Npts = Npts, inplace = True)
	L.aggrigate(collection_slice = slice(2, 500, 1))
	L.agg.plot(tmin=3.6e-4, tmax = 4e-4, ax = ax, c = "r", tunit = "us")
	M.aggrigate(collection_slice = slice(2, 500, 1))
	M.agg.plot(tmin=3.6e-4, tmax = 4e-4, ax = ax, c = "b", tunit = "us")
plt.show()