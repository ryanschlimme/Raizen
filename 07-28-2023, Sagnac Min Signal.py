# 28 July 2023
# Ryan Schlimme

# Investigating minimum detectable signal measured via Sagnac interferometer. Pulse laser firing at 12 J, creating thermoelastic acoustic event. Reducing phi on half wave plate to reduce power of pulse laser until signal to noise ratio reaches 1.

I = [145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160]

f_name_index = [r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230728\Sagnac\MinDetect\phi" + str(i) + ".tdms" for i in I]

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(r"C:\Users\ryans\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS
from acoustic_entrainment import mic_response

fig, (ax0, ax1) = plt.subplots(1,2)

M_cutoff = 1e6
L_cutoff = 1e6

L_SNR = []
M_SNR = []

for n in f_name_index:
	L = CollectionTDMS(n)
	L.set_collection("X")
	M = CollectionTDMS(n)
	M.set_collection("Y")
	L.apply("detrend", mode = "linear", inplace = True)
	M.apply("detrend", mode = "linear", inplace = True)
	L.apply("calibrate", cal = -1/0.00044, inplace = True)		# No longer need to invert signal
	M.apply("shift", tau = -11e-6, inplace = True)
	M.apply("lowpass", cutoff = M_cutoff, inplace = True)
	L.apply("lowpass", cutoff = L_cutoff, inplace = True)
	M_Npts = L.r/ (2 * M_cutoff)			# Neiquist Criterion given cutoff frequency
	L_Npts = L.r/ (2 * L_cutoff)
	L.apply("bin_average", Npts = L_Npts, inplace = True)
	M.apply("bin_average", Npts = M_Npts, inplace = True)
	M.apply("correct", response = mic_response, recollect = True) 	# Calibration of microphone
	# M.apply("calibrate", cal = 1/0.00068, inplace = True)
	for i in list(range(1, 11)):
		Li = L.collection[i]
		Mi = M.collection[i]
		L.aggrigate(collection_slice = slice(2, 500, 1))
		L.agg.plot(tmin=370e-6, tmax = 470e-6, ax = ax0, c = "r", tunit = "us")
		M.aggrigate(collection_slice = slice(2, 500, 1))
		M.agg.plot(tmin=370e-6, tmax = 470e-6, ax = ax0, c = "b", tunit = "us")
		L_maxV = []
		L_RMS = []
		M_maxV = []
		M_RMS = []
		L_maxV.append(max(Li.time_gate(tmin = 420e-6, tmax = 470e-6)[1]))
		L_RMS.append(np.sqrt(np.mean(Li.time_gate(tmin = 370e-6, tmax = 420e-6)[1] **  2)))
		M_maxV.append(max(Mi.time_gate(tmin = 420e-6, tmax = 470e-6)[1]))
		M_RMS.append(np.sqrt(np.mean(Mi.time_gate(tmin = 370e-6, tmax = 420e-6)[1] **  2)))
	L_SNR.append(np.mean(np.array(L_maxV) / np.array(L_RMS)))
	M_SNR.append(np.mean(np.array(M_maxV) / np.array(M_RMS)))
ax1.plot(I, L_SNR, color = "r")
ax1.plot(I, M_SNR, color = "b")
ax1.set_xlabel("Phi")
ax1.set_ylabel("SNR (Max V / RMS Noise)")
ax1.axhline(y=1, color = "black")
ax1.legend(["Laser", "Mic", "SNR = 1"])
plt.show()

	