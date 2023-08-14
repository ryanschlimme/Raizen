# 18 May 2023
# Ryan Schlimme

f_name_index = [r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\ene_scan_laserX_microphoneY\iter_" + str(i) + ".tdms" for i in range(6)] # create a variable pointing to file (change Ryan Schlimme to ryans)

import sys
import matplotlib.pyplot as plt
import numpy as np

sys.path.append(r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS								# pull function from time_series module
from acoustic_entrainment import mic_response
from brownian import logbin_func

N = list(range(6))
fc_list = list(range(10000, 1000000, 25))

fig, axes = plt.subplots(1, 6, sharey = True)

for n, ax in zip(N, axes.flatten()):
# Using low pass filtering OR bin averaging w/ mic correction
	f_name = f_name_index[n]	
# Initialize dummy arrays for SNR results as a function of fc
	L_SNR_array = np.array([])
	M_SNR_array = np.array([])
	L_SNR_std_array = np.array([])
	M_SNR_std_array = np.array([])
	for fc in fc_list:
		L = CollectionTDMS(f_name)
		L.set_collection("X")
		M = CollectionTDMS(f_name)
		M.set_collection("Y")
		L.apply("detrend", mode = "linear", inplace = True)
		M.apply("detrend", mode = "linear", inplace = True)
		Npts = L.r/ (2 * fc)			# Neiquist Criterion given cutoff frequency
		# L.apply("bin_average", Npts = Npts, inplace = True)
		# M.apply("bin_average", Npts = Npts, inplace = True)
		L.apply("lowpass", cutoff = fc, inplace = True)
		M.apply("lowpass", cutoff = fc, inplace = True)
		M.apply("correct", response = mic_response, recollect = True) 	# Calibration of microphone
# Initialize dummy lists for Lmax_V, Mmax_V, L_RMS, M_RMS
		Lmax_V = []
		Mmax_V = []
		L_RMS = []
		M_RMS = []

# Calculating RMS for each 
		for shot_laser, shot_micro in zip(L.collection[2:72], M.collection[2:72]):
			Lmax_V.append(max(shot_laser.time_gate(tmin = 270e-6, tmax = 300e-6)[1]))
			L_RMS.append(np.sqrt(np.mean(shot_laser.time_gate(tmin = 0, tmax = 270e-6)[1] ** 2)))
			Mmax_V.append(max(shot_micro.time_gate(tmin = 270e-6, tmax = 330e-6)[1]))
			M_RMS.append(np.sqrt(np.mean(shot_micro.time_gate(tmin = 0, tmax = 270e-6)[1] ** 2)))
		
		for shot_laser, shot_micro in zip(L.collection, M.collection):
			Lmax_V.append(max(shot_laser.time_gate(tmin = 300e-6, tmax = None)[1]))
			L_RMS.append(np.sqrt(np.mean(shot_laser.time_gate(tmin = 0, tmax = 270e-6)[1] ** 2)))
			Mmax_V.append(max(shot_micro.time_gate(tmin = 330e-6, tmax = None)[1]))
			M_RMS.append(np.sqrt(np.mean(shot_micro.time_gate(tmin = 0, tmax = 270e-6)[1] ** 2)))
# Calculating SNRmean and st dev = Lmax_V / L_RMS_array, standard deviation
		L_SNR = np.mean(np.array(Lmax_V) / np.array(L_RMS))
		M_SNR = np.mean(np.array(Mmax_V) / np.array(M_RMS))
		L_SNR_array = np.append(L_SNR_array, L_SNR)
		M_SNR_array = np.append(M_SNR_array, M_SNR)
	L_SNR_std_array = np.append(L_SNR_std_array, np.std(np.array(Lmax_V) / np.array(L_RMS)))
	M_SNR_std_array = np.append(M_SNR_std_array, np.std(np.array(Mmax_V) / np.array(M_RMS)))

# Generate SNR(fc) plots
	# np.linspace(10000, 
	# np.geomspace(start, stop, npts) #Loglog scale 
	ax.plot(fc_list, L_SNR_array, color = "r")
	ax.plot(fc_list, M_SNR_array, color = "b")
	plt.suptitle("Laser/Mic SNR Comparison")
	string = "\n" + str(n + 14) + " J"
	ax.set_title(string, fontsize = 9)
	fig.supxlabel("Cutoff Frequency ($fc$)")
	fig.supylabel("Avg of (Max Signal / RMS of Noise)")
	plt.legend(["Laser SNR", "Mic SNR"], loc = "best")
plt.show()
	
# Historgram Exploration
	#fig, (ax0, ax1, ax2, ax3) = plt.subplots(1, 4)
	#ax0.hist(L_RMS, bins = 15)
	#ax1.hist(M_RMS, bins = 15)
	#ax2.hist(Lmax_V, bins = 15)
	#ax3.hist(Mmax_V, bins = 15)
	#ax0.set_title("Laser RMS")
	#ax1.set_title("Microphone RMS")
	#ax2.set_title("Laser Max Signal")
	#ax3.set_title("Microphone Max Signal")
	#plt.show()