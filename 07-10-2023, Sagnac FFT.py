# 10 July 2023
# Ryan Schlimme

# Generating Fourier spectra and signal to noise ratios at energies scaling from 12 to 19 J in 1 J increments. Applying known calibration to microphone to transfer to pressure signal. Applying qualitative calibration to laser to match troughs. All in Sagnac interferometer.

f_name_index = [r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\Sagnac_500shot_ene_scan\iter_" + str(i) + ".tdms" for i in range(8)] # create a variable pointing to file (change Ryan Schlimme to ryans)

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS
from acoustic_entrainment import mic_response						

N = list(range(8))

for n in N:
# Using low pass filtering OR bin averaging w/ mic correction
	f_name = f_name_index[n]
	L = CollectionTDMS(f_name)
	L.set_collection("X")
	M = CollectionTDMS(f_name)
	M.set_collection("Y")
	L.apply("detrend", mode = "linear", inplace = True)
	M.apply("detrend", mode = "linear", inplace = True)
	L.apply("calibrate", cal = 1/0.002, inplace = True)
	M.apply("shift", tau = -85e-6, inplace = True)
	#M.apply("lowpass", cutoff = 2.5e6, inplace = True)
	#L.apply("lowpass", cutoff = 2.5e6, inplace = True)
	Npts = L.r/ (2 * 2.5e6)			# Neiquist Criterion given cutoff frequency
	#L.apply("bin_average", Npts = Npts, inplace = True)
	#M.apply("bin_average", Npts = Npts, inplace = True)
	M.apply("correct", response = mic_response, recollect = True) 
	# M.apply("calibrate", cal = 1/0.00068, inplace = True)
	Ldata = np.array([shot.time_gate(tmin = 450e-6, tmax = 490e-6)[1] for shot in L.collection])
	Mdata = np.array([shot.time_gate(tmin = 450e-6, tmax = 490e-6)[1] for shot in M.collection])
	Ldata_noise = np.array([shot.time_gate(tmin = 410e-6, tmax = 450e-6)[1] for shot in L.collection])
	Mdata_noise = np.array([shot.time_gate(tmin = 410e-6, tmax = 450e-6)[1] for shot in M.collection]) 
	numcols = len(Ldata[0])
	avgMfft = 0
	for i in Mdata:
		avgMfft += np.fft.rfft(i) / numcols
	avgLfft = 0
	for i in Ldata:
		avgLfft += np.fft.rfft(i) / numcols
	avgMfft_noise = 0
	for i in Mdata_noise:
		avgMfft_noise += np.fft.rfft(i) / numcols
	avgLfft_noise = 0
	for i in Ldata_noise:
		avgLfft_noise += np.fft.rfft(i) / numcols
	Lfreq = np.fft.rfftfreq(numcols, 1/L.r)
	LSNR = np.abs(avgLfft) ** 2 / np.abs(avgLfft_noise) ** 2
	MSNR = np.abs(avgMfft) ** 2 / np.abs(avgMfft_noise) ** 2
	fig, (ax0, ax1) = plt.subplots(1,2, sharex = True, figsize = [8, 4])
	string = "Flash Lamp Energy: " + str(n + 12) + " J"
	plt.suptitle(string)
	ax0.set_title("FFT Comparison", fontsize = 10.5)
	ax1.set_title("SNR Comparison", fontsize = 10.5)
	ax0.loglog(Lfreq[1:], np.abs(avgLfft[1:]) ** 2, color = "r")
	ax0.loglog(Lfreq[1:], np.abs(avgMfft[1:]) ** 2, color = "b")
	ax0.loglog(Lfreq[1:], np.abs(avgLfft_noise[1:]) ** 2, ls = "--", color = "r")
	ax0.loglog(Lfreq[1:], np.abs(avgMfft_noise[1:]) ** 2, ls = "--", color = "b")
	ax0.legend(["Laser Signal", "Mic Signal", "Laser Noise", "Mic Noise"], loc = "lower left")
	ax1.loglog(Lfreq[1:], LSNR[1:], color = "r")
	ax1.loglog(Lfreq[1:], MSNR[1:], color = "b")
	ax1.legend(["Laser SNR", "Mic SNR"], loc = "upper right")
	fig.supxlabel("Frequency (Hz)")
	plt.xlim([Lfreq[1], Lfreq[-1]])
	plt.subplots_adjust(left = 0.074, bottom = 0.11, right = 0.947, top = 0.872, wspace = 0.2, hspace = 0.2)
	plt.show()