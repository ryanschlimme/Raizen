# 01 August 2023
# Ryan Schlimme

# Comparing time domain laser pulses to calibrated microphone reading at 19 J. Applying known calibration to microphone to transfer to pressure signal. All in Sagnac interferometer with telescope and balanced photodetection. Attempt 2 to correct for potential misfire

Sagnac_name = r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\20230801\Sagnac\iter_0.tdms"
BPD_name = r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\20230801\SplitBeam\iter_0.tdms"
PD_name = r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\20230801\Telescope\iter_0.tdms"

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS
from acoustic_entrainment import mic_response						


def Phi(p):
    return np.sqrt(2) * erfinv(2*p-1)

def expected_max(N):
    mun = Phi(1-1/N)
    sigman = Phi(1-1/(N*np.e)) - mun
    return (mun + sigman*0.577)

def std_max(N):
    mun = Phi(1-1/N) + mu
    sigman = Phi(1-1/(N*np.e)) - mun
    return sigman*np.pi*np.sqrt(1/6)


def local_detrend(col, tmin = None, tmax = None, inplace = False) -> None:
	for c in col.collection:
		t, x = c.time_gate(tmin = tmin, tmax = tmax)
		x_bar = np.mean(x)
		m, b = np.polyfit(t, x, 1)
		if inplace:
			c.x = c.x - (m * c.t) - b
	return None


N = [0]
fig, ax = plt.subplots(1,1, sharey = True, sharex = True)

M_cutoff = 1e6
L_cutoff = 1e6

size = (3.375, 2.5)

for n in N:
# Using low pass filtering AND bin averaging w/ mic correction
	L = CollectionTDMS(Sagnac_name)
	L.set_collection("X")
	M = CollectionTDMS(Sagnac_name)
	M.set_collection("Y")
	local_detrend(L, tmin = 170e-6, tmax = 400e-6, inplace = True)
	local_detrend(M, tmin = 170e-6, tmax = 400e-6, inplace = True)
	L.apply("calibrate", cal = -1/0.00044, inplace = True)
	M.apply("shift", tau = -19.5e-6, inplace = True)
	M.apply("lowpass", cutoff = M_cutoff, inplace = True)
	L.apply("lowpass", cutoff = L_cutoff, inplace = True)
	M_Npts = int(L.r/ (2 * M_cutoff))			# Neiquist Criterion given cutoff frequency
	L_Npts = int(L.r/ (2 * L_cutoff))
	L.apply("bin_average", Npts = L_Npts, inplace = True)
	M.apply("bin_average", Npts = M_Npts, inplace = True)
	M.apply("correct", response = mic_response, recollect = True) 	# Calibration of microphone
	# M.apply("calibrate", cal = 1/0.00068, inplace = True)
	L.aggrigate(collection_slice = slice(2, 500, 1))
	L.agg.plot(tmin=400e-6, tmax = 470e-6, ax = ax, figsize = size, tunit = "us")
	M.aggrigate(collection_slice = slice(2, 500, 1))
	M.agg.plot(tmin=400e-6, tmax = 470e-6, ax = ax, figsize = size, c = "b", tunit = "us")
	for i in list(range(2, 3)):
		Li = L.collection[i]
		Mi = M.collection[i]
		#Li.plot(tmin = 400e-6, tmax = 470e-6, ax = ax, tunit = "us")
		#Mi.plot(tmin = 420e-6, tmax = 470e-6, ax = ax, c = "C0", tunit = "us")

for n in N:
# Using low pass filtering AND bin averaging w/ mic correction
	L = CollectionTDMS(BPD_name)
	L.set_collection("X")
	local_detrend(L, tmin = 170e-6, tmax = 400e-6, inplace = True)
	L.apply("calibrate", cal = 1/0.0002, inplace = True)
	L.apply("lowpass", cutoff = L_cutoff, inplace = True)
	L_Npts = int(L.r/ (2 * L_cutoff))
	L.apply("bin_average", Npts = L_Npts, inplace = True)
	L.aggrigate(collection_slice = slice(2, 500, 1))
	L.agg.plot(tmin=400e-6, tmax = 470e-6, ax = ax, figsize = size, tunit = "us")
	for i in list(range(2, 3)):
		Li = L.collection[i]
		Mi = M.collection[i]
		#Li.plot(tmin = 400e-6, tmax = 470e-6, ax = ax, c = "red", tunit = "us")
		#Mi.plot(tmin = 420e-6, tmax = 470e-6, ax = ax, c = "C0", tunit = "us")

for n in N:
# Using low pass filtering AND bin averaging w/ mic correction
	L = CollectionTDMS(PD_name)
	L.set_collection("X")
	local_detrend(L, tmin = 170e-6, tmax = 400e-6, inplace = True)
	L.apply("calibrate", cal = -1/0.0000002, inplace = True)	
	L.apply("lowpass", cutoff = L_cutoff, inplace = True)
	L_Npts = int(L.r/ (2 * L_cutoff))
	L.apply("bin_average", Npts = L_Npts, inplace = True)
	L.aggrigate(collection_slice = slice(2, 500, 1))
	L.agg.plot(tmin=400e-6, tmax = 470e-6, ax = ax, figsize = size, tunit = "us")
	for i in list(range(2, 3)):
		Li = L.collection[i]
		Mi = M.collection[i]
		#Li.plot(tmin = 400e-6, tmax = 470e-6, ax = ax, color = "r", tunit = "us")
		#Mi.plot(tmin = 420e-6, tmax = 470e-6, ax = ax, c = "C0", tunit = "us")


string = "\n" + str(19) + " J"
plt.title(string, fontsize = 9)
plt.suptitle("Acoustic Detection by Four Methods")
plt.legend(["Sagnac", "Mic", "BPD", "PD"])
plt.show()

#plt.savefig("Time-Domain Pulses.png")