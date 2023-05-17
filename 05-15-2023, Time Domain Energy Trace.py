# 15 May 2023
# Ryan Schlimme

# Comparing time domain laser pulses to calibrated microphone reading at energies scaling from 14 to 19 J in 1 J increments. Applying known 0.00068 V/Pa offset to microphone to transfer to pressure signal. Still searching for laser transfer function.
# Noticed we are matching low frequecies well but exhibiting high frequency oscillations. Quantum noise? Or actual detection?

f_name_index = [r"C:\Users\ryans\OneDrive\Desktop\Research\Data\ene_scan\laser-X_microphone-Y0.00068V-per-Pa\iter_" + str(i) + ".tdms" for i in range(6)] # create a variable pointing to file (change Ryan Schlimme to ryans)

import sys
import matplotlib.pyplot as plt

sys.path.append(r"C:\Users\ryans\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS								# pull function from time_series module


N = list(range(6))

for n in N:
	f_name = f_name_index[n]
	L = CollectionTDMS(f_name)
	L.set_collection("X")
	M = CollectionTDMS(f_name)
	M.set_collection("Y")
##### Time Analysis to Compare Pulses ######
	L.apply("calibrate", cal = 30, inplace = True)			# Vertically stretch, detrend, and horizontally shift signals to align
	L.apply("detrend", mode = "linear", inplace = True)
	M.apply("calibrate", cal = 0.185, inplace = True)
	M.apply("detrend", mode = "linear", inplace = True)
	M.apply("shift", tau = 0.0000235, inplace = True)
	Npts = L.r / (2*500000)							# Bin_average to max frequency of 500 kHz
	L.apply("bin_average", Npts = Npts, inplace = True)
	M.apply("bin_average", Npts = Npts, inplace = True)
##### Aggregate Plots #####
	fig, ax = plt.subplots(1,1)
	L.aggrigate(collection_slice = slice(2, 100, 1))
	L.agg.plot(tmin=2.6e-4, tmax = 3e-4, ax = ax, c = "r", tunit = "us")
	M.aggrigate(collection_slice = slice(2, 100, 1))
	M.agg.plot(tmin=2.6e-4, tmax = 3e-4, ax = ax, c = "b", tunit = "us")
	plt.suptitle("Aggregate Pulses of Laser Ablation \n as Recorded by Microphone and Laser Deflection", fontsize = 11)
	string = "\n" + str(n + 14) + " J"
	plt.title(string, fontsize = 9)
	plt.xlabel("Time (us)")
	plt.ylabel("Voltage (V)")
	plt.legend(["Laser", "Microphone"], loc = "best")
plt.show()