# 01 August 2023
# Ryan Schlimme

# Investigating signal to noise ratio of four acoustic measurement systems.

Sagnac_name = r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\20230801\Sagnac\iter_0.tdms" # create a variable pointing to file (change Ryan Schlimme to ryans)
BPD_name = r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\20230801\SplitBeam\iter_0.tdms"
PD_name = r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\20230801\Telescope\iter_0.tdms"

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS
from acoustic_entrainment import mic_response

fig, ax = plt.subplots(1,1)

fc_list = list(np.geomspace(10000, 3000000, 10))

M_SNR = []	# Mic
Sagnac_SNR = []	# Sagnac
BPD_SNR = []	# balanced photodetection
PD_SNR = []	# photodiode

for fc in fc_list:
	L = CollectionTDMS(Sagnac_name)
	L.set_collection("X")
	M = CollectionTDMS(Sagnac_name)
	M.set_collection("Y")
	L.apply("detrend", mode = "linear", inplace = True)
	M.apply("detrend", mode = "linear", inplace = True)
	L.apply("calibrate", cal = -1/0.00044, inplace = True)
	M.apply("shift", tau = -19.5e-6, inplace = True)
	M.apply("lowpass", cutoff = fc, inplace = True)
	L.apply("lowpass", cutoff = fc, inplace = True)
	Npts = L.r/ (2 * fc)			# Neiquist Criterion given cutoff frequency
	L.apply("bin_average", Npts = Npts, inplace = True)
	M.apply("bin_average", Npts = Npts, inplace = True)
	M.apply("correct", response = mic_response, recollect = True) 	# Calibration of microphone
	for i in list(range(1, 11)):
		Li = L.collection[i]
		Mi = M.collection[i]
		L_maxV = []
		L_RMS = []
		M_maxV = []
		M_RMS = []
		L_maxV.append(max(Li.time_gate(tmin = 420e-6, tmax = 470e-6)[1]))
		L_RMS.append(np.sqrt(np.mean(Li.time_gate(tmin = 370e-6, tmax = 420e-6)[1] **  2)))
		M_maxV.append(max(Mi.time_gate(tmin = 420e-6, tmax = 470e-6)[1]))
		M_RMS.append(np.sqrt(np.mean(Mi.time_gate(tmin = 370e-6, tmax = 420e-6)[1] **  2)))
	Sagnac_SNR.append(np.mean(np.array(L_maxV) / np.array(L_RMS)))
	M_SNR.append(np.mean(np.array(M_maxV) / np.array(M_RMS)))


for fc in fc_list:
	L = CollectionTDMS(BPD_name)
	L.set_collection("X")
	L.apply("detrend", mode = "linear", inplace = True)
	L.apply("calibrate", cal = -1/0.00044, inplace = True)
	L.apply("lowpass", cutoff = fc, inplace = True)
	Npts = L.r/ (2 * fc)			# Neiquist Criterion given cutoff frequency
	L.apply("bin_average", Npts = Npts, inplace = True)
	for i in list(range(1, 11)):
		Li = L.collection[i]
		L_maxV = []
		L_RMS = []
		L_maxV.append(max(Li.time_gate(tmin = 420e-6, tmax = 470e-6)[1]))
		L_RMS.append(np.sqrt(np.mean(Li.time_gate(tmin = 370e-6, tmax = 420e-6)[1] **  2)))
	BPD_SNR.append(np.mean(np.array(L_maxV) / np.array(L_RMS)))

for fc in fc_list:
	L = CollectionTDMS(PD_name)
	L.set_collection("X")
	L.apply("detrend", mode = "linear", inplace = True)
	L.apply("calibrate", cal = -1/0.00000024, inplace = True)
	L.apply("lowpass", cutoff = fc, inplace = True)
	Npts = L.r/ (2 * fc)			# Neiquist Criterion given cutoff frequency
	L.apply("bin_average", Npts = Npts, inplace = True)
	for i in list(range(1, 11)):
		Li = L.collection[i]
		L_maxV = []
		L_RMS = []
		L_maxV.append(max(Li.time_gate(tmin = 420e-6, tmax = 470e-6)[1]))
		L_RMS.append(np.sqrt(np.mean(Li.time_gate(tmin = 370e-6, tmax = 420e-6)[1] **  2)))
	PD_SNR.append(np.mean(np.array(L_maxV) / np.array(L_RMS)))

ax.loglog(fc_list, Sagnac_SNR)
ax.loglog(fc_list, BPD_SNR)
ax.loglog(fc_list, PD_SNR)
ax.loglog(fc_list, M_SNR, color = "b")
plt.xlabel("Frequency Cutoff (fc)")
plt.ylabel("SNR (Max V / RMS Noise)")
plt.axhline(y=1, color = "black")
plt.legend(["Sagnac", "BPD", "Tele", "Mic", "SNR = 1"])
plt.show()

	