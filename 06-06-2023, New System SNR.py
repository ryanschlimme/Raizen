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
fig, axes = plt.subplots(1, 8, sharey = True)
fc_list = list(np.geomspace(10000, 3000000, 10))

for n, ax in zip(N, axes.flatten()):
# Using low pass filtering OR bin averaging w/ mic correction
	f_name = f_name_index[n]
# Initialize dummy arrays for SNR results as a function of fc
	L_SNR_array = np.array([])
	M_SNR_array = np.array([])
	Lmax_V_array = np.array([])
	Mmax_V_array = np.array([])
	L_RMS_array = np.array([])
	M_RMS_array = np.array([])
	L_SNR_std_array = np.array([])
	M_SNR_std_array = np.array([])
	for fc in fc_list:
		L = CollectionTDMS(f_name)
		L.set_collection("X")
		M = CollectionTDMS(f_name)
		M.set_collection("Y")
		L.apply("detrend", mode = "linear", inplace = True)
		M.apply("detrend", mode = "linear", inplace = True)
		L.apply("calibrate", cal = -1, inplace = True)
		M.apply("shift", tau = 35e-6, inplace = True)
		M.apply("lowpass", cutoff = fc, inplace = True)
		L.apply("lowpass", cutoff = fc, inplace = True)
		Npts = L.r/ (2 * fc)			# Neiquist Criterion given cutoff frequency
		L.apply("bin_average", Npts = Npts, inplace = True)
		M.apply("bin_average", Npts = Npts, inplace = True)
		M.apply("correct", response = mic_response, recollect = True) 	# Calibration of microphone
# Initialize dummy lists for Lmax_V, Mmax_V, L_RMS, M_RMS
		Lmax_V = []
		Mmax_V = []
		L_RMS = []
		M_RMS = []
# Calculating RMS for each 
		for shot_laser, shot_micro in zip(L.collection[2:], M.collection[2:]):
			Lmax_V.append(max(shot_laser.time_gate(tmin = 360e-6, tmax = 420e-6)[1]))
			L_RMS.append(np.sqrt(np.mean(shot_laser.time_gate(tmin = 300e-6, tmax = 360e-6)[1] ** 2)))
			Mmax_V.append(max(shot_micro.time_gate(tmin = 360e-6, tmax = 420e-6)[1]))
			M_RMS.append(np.sqrt(np.mean(shot_micro.time_gate(tmin = 300e-6, tmax = 360e-6)[1] ** 2)))
# Calculating SNRmean and st dev = Lmax_V / L_RMS_array, standard deviation
		L_SNR = np.mean(np.array(Lmax_V) / np.array(L_RMS))
		M_SNR = np.mean(np.array(Mmax_V) / np.array(M_RMS))
		L_SNR_array = np.append(L_SNR_array, L_SNR)
		M_SNR_array = np.append(M_SNR_array, M_SNR)
		Lmax_V_array = np.append(Lmax_V_array, np.mean(Lmax_V))
		Mmax_V_array = np.append(Mmax_V_array, np.mean(Mmax_V))
		L_RMS_array = np.append(L_RMS_array, np.mean(L_RMS))
		M_RMS_array = np.append(M_RMS_array, np.mean(M_RMS))
	# fig, (ax0, ax1) = plt.subplots(1,2)
	# ax0.plot(fc_list, Lmax_V_array)
	#ax0.plot(fc_list, Mmax_V_array)
	#ax1.plot(fc_list, L_RMS_array)
	#ax1.plot(fc_list, M_RMS_array)
	#ax0.legend(["Lmax", "Mmax"], loc = "best")
	#ax1.legend(["Lrms", "Mrms"], loc = "best")
	L_SNR_std_array = np.append(L_SNR_std_array, np.std(np.array(Lmax_V) / np.array(L_RMS)))
	M_SNR_std_array = np.append(M_SNR_std_array, np.std(np.array(Mmax_V) / np.array(M_RMS)))
	
# Generate SNR(fc) plots
	# np.linspace(10000, 
	# np.geomspace(start, stop, npts) #Loglog scale 
	ax.loglog(fc_list, L_SNR_array, color = "r")
	ax.loglog(fc_list, M_SNR_array, color = "b")
	plt.suptitle("Laser/Mic SNR Comparison")
	string = "\n" + str(n + 12) + " J"
	ax.set_title(string, fontsize = 9)
	fig.supxlabel("Cutoff Frequency ($fc$)")
	fig.supylabel("Avg of (Max Signal / RMS of Noise)")
	plt.legend(["Laser SNR", "Mic SNR"], loc = "best")
plt.show()