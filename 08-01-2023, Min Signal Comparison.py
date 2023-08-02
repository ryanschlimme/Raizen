# 01 August 2023
# Ryan Schlimme

# Investigating minimum detectable signal measured via Sagnac interferometer. Pulse laser firing at 12 J, creating thermoelastic acoustic event. Reducing phi on half wave plate to reduce power of pulse laser until signal to noise ratio reaches 1.

I = [145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160]

Sagnac_name_index = [r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\20230801\Sagnac\MinDetect\phi" + str(i) + ".tdms" for i in I]
BPD_name_index = [r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\20230801\SplitBeam\MinDetect\phi" + str(i) + ".tdms" for i in I]
PD_name_index = [r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\20230801\Telescope\MinDetect\phi" + str(i) + ".tdms" for i in I]

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS
from acoustic_entrainment import mic_response

fig, ax = plt.subplots(1,1)

M_cutoff = 200e3
Sagnac_cutoff = 1e6
BPD_cutoff = 1e6
PD_cutoff = 1e6

M_SNR = []	# Mic
Sagnac_SNR = []	# Sagnac
BPD_SNR = []	# balanced photodetection
PD_SNR = []	# photodiode

for n in Sagnac_name_index:
	L = CollectionTDMS(n)
	L.set_collection("X")
	M = CollectionTDMS(n)
	M.set_collection("Y")
	L.apply("detrend", mode = "linear", inplace = True)
	M.apply("detrend", mode = "linear", inplace = True)
	L.apply("calibrate", cal = -1/0.00044, inplace = True)		
	M.apply("shift", tau = -19.5e-6, inplace = True)
	M.apply("lowpass", cutoff = M_cutoff, inplace = True)
	L.apply("lowpass", cutoff = Sagnac_cutoff, inplace = True)
	M_Npts = L.r/ (2 * M_cutoff)			# Neiquist Criterion given cutoff frequency
	L_Npts = L.r/ (2 * Sagnac_cutoff)
	L.apply("bin_average", Npts = L_Npts, inplace = True)
	M.apply("bin_average", Npts = M_Npts, inplace = True)
	M.apply("correct", response = mic_response, recollect = True) 	# Calibration of microphone
	# M.apply("calibrate", cal = 1/0.00068, inplace = True)
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

for n in BPD_name_index:
	L = CollectionTDMS(n)
	L.set_collection("X")
	L.apply("detrend", mode = "linear", inplace = True)
	L.apply("calibrate", cal = -1/0.00044, inplace = True)		
	L.apply("lowpass", cutoff = BPD_cutoff, inplace = True)
	L_Npts = L.r/ (2 * BPD_cutoff)
	L.apply("bin_average", Npts = L_Npts, inplace = True)
	for i in list(range(1, 11)):
		Li = L.collection[i]
		L_maxV = []
		L_RMS = []
		L_maxV.append(max(Li.time_gate(tmin = 420e-6, tmax = 470e-6)[1]))
		L_RMS.append(np.sqrt(np.mean(Li.time_gate(tmin = 370e-6, tmax = 420e-6)[1] **  2)))
	BPD_SNR.append(np.mean(np.array(L_maxV) / np.array(L_RMS)))

for n in PD_name_index:
	L = CollectionTDMS(n)
	L.set_collection("X")
	L.apply("detrend", mode = "linear", inplace = True)
	L.apply("calibrate", cal = -1/0.00000024, inplace = True)		
	L.apply("lowpass", cutoff = PD_cutoff, inplace = True)
	L_Npts = L.r/ (2 * PD_cutoff)
	L.apply("bin_average", Npts = L_Npts, inplace = True)
	for i in list(range(1, 11)):
		Li = L.collection[i]
		L_maxV = []
		L_RMS = []
		L_maxV.append(max(Li.time_gate(tmin = 420e-6, tmax = 470e-6)[1]))
		L_RMS.append(np.sqrt(np.mean(Li.time_gate(tmin = 370e-6, tmax = 420e-6)[1] **  2)))
	PD_SNR.append(np.mean(np.array(L_maxV) / np.array(L_RMS)))

ax.plot(I, Sagnac_SNR)
ax.plot(I, BPD_SNR)
ax.plot(I, PD_SNR)
ax.plot(I, M_SNR, color = "b")
plt.xlabel("Phi")
plt.ylabel("SNR (Max V / RMS Noise)")
plt.axhline(y=1, color = "black")
plt.legend(["Sagnac", "BPD", "Tele", "Mic", "SNR = 1"])
plt.show()

	