# 01 August 2023
# Ryan Schlimme

# Investigating signal to noise ratio of four acoustic measurement systems.

Sagnac_name = r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\Sagnac\iter_0.tdms"
BPD_name = r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\SplitBeam\iter_0.tdms"
PD_name = r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\Telescope\iter_0.tdms"

#Sagnac_name = r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\Sagnac\MinDetect\phi145.tdms"
#BPD_name = r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\SplitBeam\MinDetect\phi145.tdms"
#PD_name = r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\Telescope\MinDetect\phi145.tdms"

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfinv

sys.path.append(r"C:\Users\ryans\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS
from acoustic_entrainment import mic_response

fig, ax = plt.subplots(1,1)

fc_list = list(np.linspace(20000, 2500000, 25))


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

def stats(M_SNR):
	print("Maximum SNR: ", max(M_SNR))
	print("Corresponding fc: ", fc_list[M_SNR.index(max(M_SNR))] / 100000, " kHz")
	print("Last Value: ", M_SNR[-1]) 
	return None


Sagnac_SNR = []
M_SNR = []
BPD_SNR = []
PD_SNR = []

for fc in fc_list:
	L = CollectionTDMS(Sagnac_name)
	L.set_collection("X")
	local_detrend(L, tmin = 170e-6, tmax = 400e-6, inplace = True)
	#L.apply("detrend", mode = "linear", inplace = True)
	L.apply("calibrate", cal = -1/0.00044, inplace = True)
	L.apply("lowpass", cutoff = fc, inplace = True)
	Npts = L.r/ (2 * fc)			# Neiquist Criterion given cutoff frequency
	L.apply("bin_average", Npts = Npts, inplace = True)
	L_maxV = []
	L_RMS = []
	for i in list(range(1, 51)):
		Li = L.collection[i]
		L_maxV.append(max(Li.time_gate(tmin = 400e-6, tmax = 470e-6)[1]))
		L_RMS.append(np.std(Li.time_gate(tmin = 170e-6, tmax = 400e-6)[1]))
	Sagnac_SNR.append(np.mean(np.array(L_maxV) / np.array(L_RMS) / expected_max(len(Li.time_gate(tmin = 170e-6, tmax = 400e-6)[1])))

for fc in fc_list:
	M = CollectionTDMS(Sagnac_name)
	M.set_collection("Y")
	local_detrend(M, tmin = 170e-6, tmax = 400e-6, inplace = True)
	#M.apply("detrend", mode = "linear", inplace = True)
	M.apply("shift", tau = -19.5e-6, inplace = True)
	M.apply("lowpass", cutoff = fc, inplace = True)
	Npts = M.r/ (2 * fc)			# Neiquist Criterion given cutoff frequency
	M.apply("bin_average", Npts = Npts, inplace = True)
	M.apply("correct", response = mic_response, recollect = True) 	# Calibration of microphone
	M_maxV = []
	M_RMS = []
	for i in list(range(1, 51)):
		Mi = M.collection[i]
		M_maxV.append(max(Mi.time_gate(tmin = 400e-6, tmax = 470e-6)[1]))
		M_RMS.append(np.std(Mi.time_gate(tmin = 170e-6, tmax = 400e-6)[1]))
	M_SNR.append(np.mean(np.array(M_maxV) / np.array(M_RMS) / expected_max(len(Mi.time_gate(tmin = 170e-6, tmax = 400e-6)[1]))))

for fc in fc_list:
	L = CollectionTDMS(BPD_name)
	L.set_collection("X")
	local_detrend(L, tmin = 170e-6, tmax = 400e-6, inplace = True)
	#L.apply("detrend", mode = "linear", inplace = True)
	L.apply("calibrate", cal = 1/0.0002, inplace = True)
	L.apply("lowpass", cutoff = fc, inplace = True)
	Npts = L.r/ (2 * fc)			# Neiquist Criterion given cutoff frequency
	L.apply("bin_average", Npts = Npts, inplace = True)
	L_maxV = []
	L_RMS = []
	for i in list(range(1, 51)):
		Li = L.collection[i]
		L_maxV.append(max(Li.time_gate(tmin = 400e-6, tmax = 470e-6)[1]))
		L_RMS.append(np.std(Li.time_gate(tmin = 170e-6, tmax = 400e-6)[1]))
	BPD_SNR.append(np.mean(np.array(L_maxV) / np.array(L_RMS) / expected_max(len(Li.time_gate(tmin = 170e-6, tmax = 400e-6)[1]))))

for fc in fc_list:
	L = CollectionTDMS(PD_name)
	L.set_collection("X")
	local_detrend(L, tmin = 170e-6, tmax = 400e-6, inplace = True)
	#L.apply("detrend", mode = "linear", inplace = True)
	L.apply("calibrate", cal = -1/0.00000024, inplace = True)
	L.apply("lowpass", cutoff = fc, inplace = True)
	Npts = L.r/ (2 * fc)			# Neiquist Criterion given cutoff frequency
	L.apply("bin_average", Npts = Npts, inplace = True)
	L_maxV = []
	L_RMS = []
	for i in list(range(1, 51)):
		Li = L.collection[i]
		L_maxV.append(max(Li.time_gate(tmin = 400e-6, tmax = 470e-6)[1]))
		L_RMS.append(np.std(Li.time_gate(tmin = 170e-6, tmax = 400e-6)[1]))
	PD_SNR.append(np.mean(np.array(L_maxV) / np.array(L_RMS) / expected_max(len(Li.time_gate(tmin = 170e-6, tmax = 400e-6)[1]))))

ax.plot(fc_list, Sagnac_SNR)
ax.plot(fc_list, BPD_SNR)
ax.plot(fc_list, PD_SNR)
ax.plot(fc_list, M_SNR)
plt.axhline(y=1, color = "black")

plt.xlabel("Frequency Cutoff (fc)")
plt.ylabel("SNR (Max V / St Dev Noise)")

plt.legend(["Sagnac", "BPD", "Tele", "Mic", "SNR = 1"])

plt.suptitle("SNR Comparison of Four Methods")
plt.title("Flash Lamp Energy: 19 J", fontsize = 9)

plt.show()
	