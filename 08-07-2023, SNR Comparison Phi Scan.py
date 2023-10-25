# 07 August 2023
# Ryan Schlimme

# Investigating signal to noise ratio of four acoustic measurement systems with various signal strengths using parallel processing.

I = [152]
# I = [145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160]

Sagnac_name_index = [r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\20230801\Sagnac\MinDetect\phi" + str(i) + ".tdms" for i in I]
BPD_name_index = [r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\20230801\SplitBeam\MinDetect\phi" + str(i) + ".tdms" for i in I]
PD_name_index = [r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\20230801\Telescope\MinDetect\phi" + str(i) + ".tdms" for i in I]

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfinv
from joblib import Parallel, delayed

sys.path.append(r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\brownian\src") 
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

mu = 0

fig, axes = plt.subplots(1, 2, sharex = False)
fc_list = list(np.linspace(20000, 2500000, 50))

Iteration = [c for c in range(len(I))]
	
def Sagnac2(z, fc, tmin = 170e-6, tmax = 400e-6, cal = -1/0.00044):
	Sagnac_name = r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\20230801\Sagnac\MinDetect\phi" + str(z) + ".tdms"
	L = CollectionTDMS(Sagnac_name)
	L.set_collection("X")
	local_detrend(L, tmin = tmin, tmax = tmax, inplace = True)
	#L.apply("detrend", mode = "linear", inplace = True)
	L.apply("calibrate", cal = cal, inplace = True)
	L.apply("lowpass", cutoff = fc, inplace = True)
	Npts = int(L.r/ (2 * fc))			# Nyquist Criterion given cutoff frequency
	L.apply("bin_average", Npts = Npts, inplace = True)
	L_maxV = []
	L_RMS = []
	for i in list(range(1, 51)):
		Li = L.collection[i]
		L_maxV.append(max(Li.time_gate(tmin = 400e-6, tmax = 470e-6)[1]))
		L_RMS.append(np.std(Li.time_gate(tmin = 170e-6, tmax = 350e-6)[1]))
	SNR = np.mean(np.array(L_maxV) / np.array(L_RMS) / expected_max(len(Li.time_gate(tmin = 170e-6, tmax = 350e-6)[1])))
	STD = std_max(len(Li.time_gate(tmin = 170e-6, tmax = 350e-6)[1])) / expected_max(len(Li.time_gate(tmin = 170e-6, tmax = 350e-6)[1]))
	return SNR, STD

args = []

for i in I:
	for fc in fc_list:
		args.append((i, fc))

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
		Npts = int(M.r/ (2 * fc))			# Nyquist Criterion given cutoff frequency
		M.apply("bin_average", Npts = Npts, inplace = True)
		M.apply("correct", response = mic_response, recollect = True) 	# Calibration of microphone
		M_maxV = []
		M_RMS = []
		for i in list(range(1, 51)):
			Mi = M.collection[i]
			M_maxV.append(max(Mi.time_gate(tmin = 400e-6, tmax = 470e-6)[1]))
			M_RMS.append(np.std(Mi.time_gate(tmin = 170e-6, tmax = 350e-6)[1]))
		M_SNR_inter.append(np.mean(np.array(M_maxV) / np.array(M_RMS) / expected_max(len(Mi.time_gate(tmin = 170e-6, tmax = 350e-6)[1]))))
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
		Npts = int(L.r/ (2 * fc))			# Neiquist Criterion given cutoff frequency
		L.apply("bin_average", Npts = Npts, inplace = True)
		L_maxV = []
		L_RMS = []
		for i in list(range(1, 51)):
			Li = L.collection[i]
			L_maxV.append(max(Li.time_gate(tmin = 400e-6, tmax = 470e-6)[1]))
			L_RMS.append(np.std(Li.time_gate(tmin = 170e-6, tmax = 350e-6)[1]))
		BPD_SNR_inter.append(np.mean(np.array(L_maxV) / np.array(L_RMS) / expected_max(len(Li.time_gate(tmin = 170e-6, tmax = 350e-6)[1]))))
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
		Npts = int(L.r/ (2 * fc))			# Neiquist Criterion given cutoff frequency
		L.apply("bin_average", Npts = Npts, inplace = True)
		L_maxV = []
		L_RMS = []
		for i in list(range(1, 51)):
			Li = L.collection[i]
			L_maxV.append(max(Li.time_gate(tmin = 400e-6, tmax = 470e-6)[1]))
			L_RMS.append(np.std(Li.time_gate(tmin = 170e-6, tmax = 350e-6)[1]))
		PD_SNR_inter.append(np.mean(np.array(L_maxV) / np.array(L_RMS) / expected_max(len(Li.time_gate(tmin = 170e-6, tmax = 350e-6)[1]))))
	return PD_SNR_inter

SagnacSNRSTDlist = np.array(Parallel(n_jobs = -1)(delayed(Sagnac2)(*arg) for arg in args))
SagnacSNRSTDlist = SagnacSNRSTDlist.reshape(len(I), len(fc_list), 2)
M_SNR = Parallel(n_jobs = -1)(delayed(M_loop)(i) for i in Iteration)
BPD_SNR = Parallel(n_jobs = -1)(delayed(BPD_loop)(i) for i in Iteration)
PD_SNR = Parallel(n_jobs = -1)(delayed(PD_loop)(i) for i in Iteration)

print("Max Sagnac", max(SagnacSNRSTDlist[0,:,0]))
print("Max BPD", max(BPD_SNR[0]))
print("Max Mic", max(M_SNR[0]))
print("Max PD", max(PD_SNR[0]))

for i, z in enumerate(I):
	ax = axes.flatten()[i]
	ax.plot(fc_list, SagnacSNRSTDlist[i, :, 0], label = "Sagnac")
	d1 = SagnacSNRSTDlist[i, :, 1]
	ones = np.ones_like(d1)
	#ax.errorbar(fc_list, ones, yerr = d1)
	ax.fill_between(fc_list, ones - d1, ones + d1, color = "k", alpha = 0.4)

for z, ax, i in zip(Iteration, axes.flatten(), I):
	ax.plot(fc_list, BPD_SNR[z], label = "BPD")
	ax.plot(fc_list, PD_SNR[z], label = "PD")
	ax.plot(fc_list, M_SNR[z], label = "Mic")
	ax.set_title("Phi: " + str(i))

#plt.xlabel("Frequency Cutoff (fc)")
#plt.ylabel("SNR (Max V / St Dev Noise)")

ax.legend()

plt.suptitle("SNR Comparison of Four Methods")
#plt.title("Flash Lamp Energy: 19 J", fontsize = 9)

plt.show()