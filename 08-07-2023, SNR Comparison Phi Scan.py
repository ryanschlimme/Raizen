# 07 August 2023
# Ryan Schlimme

# Investigating signal to noise ratio of four acoustic measurement systems with various signal strengths using parallel processing.

I = [145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160]

Sagnac_name_index = [r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\Sagnac\MinDetect\phi" + str(i) + ".tdms" for i in I]
BPD_name_index = [r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\SplitBeam\MinDetect\phi" + str(i) + ".tdms" for i in I]
PD_name_index = [r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\Telescope\MinDetect\phi" + str(i) + ".tdms" for i in I]

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfinv
from joblib import Parallel, delayed

sys.path.append(r"C:\Users\ryans\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
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


fig, axes = plt.subplots(4, 4, sharex = True, sharey = "row")
fc_list = list(np.linspace(20000, 2500000, 25))

Iteration = [c for c in range(len(I))]
	
def Sagnac_loop(z):
	Sagnac_name = Sagnac_name_index[z]
	Sagnac_SNR_inter = []
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
		Sagnac_SNR_inter.append(np.mean(np.array(L_maxV) / np.array(L_RMS) / expected_max(len(Li.time_gate(tmin = 170e-6, tmax = 400e-6)[1]))))
	return Sagnac_SNR_inter

def M_loop(z):
	Sagnac_name = Sagnac_name_index[z]
	M_SNR_inter = []
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
		M_SNR_inter.append(np.mean(np.array(M_maxV) / np.array(M_RMS) / expected_max(len(Mi.time_gate(tmin = 170e-6, tmax = 400e-6)[1]))))
	return M_SNR_inter

def BPD_loop(z):
	BPD_name = BPD_name_index[z]
	BPD_SNR_inter = []
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
		BPD_SNR_inter.append(np.mean(np.array(L_maxV) / np.array(L_RMS) / expected_max(len(Li.time_gate(tmin = 170e-6, tmax = 400e-6)[1]))))
	return BPD_SNR_inter


def PD_loop(z):
	PD_name = PD_name_index[z]
	PD_SNR_inter = []
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
		PD_SNR_inter.append(np.mean(np.array(L_maxV) / np.array(L_RMS) / expected_max(len(Li.time_gate(tmin = 170e-6, tmax = 400e-6)[1]))))
	return PD_SNR_inter

Sagnac_SNR = Parallel(n_jobs = -1)(delayed(Sagnac_loop)(i) for i in Iteration)
M_SNR = Parallel(n_jobs = -1)(delayed(M_loop)(i) for i in Iteration)
BPD_SNR = Parallel(n_jobs = -1)(delayed(BPD_loop)(i) for i in Iteration)
PD_SNR = Parallel(n_jobs = -1)(delayed(PD_loop)(i) for i in Iteration)

for z, ax, i in zip(Iteration, axes.flatten(), I):
	ax.plot(fc_list, Sagnac_SNR[z])
	ax.plot(fc_list, BPD_SNR[z])
	ax.plot(fc_list, PD_SNR[z])
	ax.plot(fc_list, M_SNR[z])
	ax.axhline(y=1, color = "black")
	ax.set_title("Phi: " + str(i))

#plt.xlabel("Frequency Cutoff (fc)")
#plt.ylabel("SNR (Max V / St Dev Noise)")

plt.legend(["Sagnac", "BPD", "Tele", "Mic", "SNR = 1"])

plt.suptitle("SNR Comparison of Four Methods")
#plt.title("Flash Lamp Energy: 19 J", fontsize = 9)

plt.show()
