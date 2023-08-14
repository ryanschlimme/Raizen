# 07 August 2023
# Ryan Schlimme

# Generating Fourier spectra and signal to noise ratios. Applying known calibration to microphone to transfer to pressure signal. Applying qualitative calibration to laser to match troughs. Comparing all acoustoptical systems to microphone.

Sagnac_name = r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\Sagnac\iter_0.tdms"
BPD_name = r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\SplitBeam\iter_0.tdms"
PD_name = r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\Telescope\iter_0.tdms"

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfinv

sys.path.append(r"C:\Users\ryans\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS
from acoustic_entrainment import mic_response						


def local_detrend(col, tmin = None, tmax = None, inplace = False) -> None:
	for c in col.collection:
		t, x = c.time_gate(tmin = tmin, tmax = tmax)
		x_bar = np.mean(x)
		m, b = np.polyfit(t, x, 1)
		if inplace:
			c.x = c.x - (m * c.t) - b
	return None


fig, axes = plt.subplots(3,2, sharex = True, sharey = True, figsize = [8, 4])
cutoff = 2.1e6


# Sagnac and Mic
L = CollectionTDMS(Sagnac_name)
L.set_collection("X")
M = CollectionTDMS(Sagnac_name)
M.set_collection("Y")
local_detrend(L, tmin = 170e-6, tmax = 400e-6, inplace = True)
local_detrend(M, tmin = 170e-6, tmax = 400e-6, inplace = True)
#L.apply("detrend", mode = "linear", inplace = True)
#M.apply("detrend", mode = "linear", inplace = True)
L.apply("calibrate", cal = -1/0.00044, inplace = True)
M.apply("shift", tau = -19.5e-6, inplace = True)
#M.apply("lowpass", cutoff = 1e6, inplace = True)
#L.apply("lowpass", cutoff = 1e6, inplace = True)
Npts = L.r/ (2 * cutoff)			# Nyquist Criterion given cutoff frequency
L.apply("bin_average", Npts = Npts, inplace = True)
M.apply("bin_average", Npts = Npts, inplace = True)
M.apply("correct", response = mic_response, recollect = True) 
# M.apply("calibrate", cal = 1/0.00068, inplace = True)
Ldata = np.array([shot.time_gate(tmin = 400e-6, tmax = 470e-6)[1] for shot in L.collection])
Mdata = np.array([shot.time_gate(tmin = 400e-6, tmax = 470e-6)[1] for shot in M.collection])
Ldata_noise = np.array([shot.time_gate(tmin = 330e-6, tmax = 400e-6)[1] for shot in L.collection])
Mdata_noise = np.array([shot.time_gate(tmin = 330e-6, tmax = 400e-6)[1] for shot in M.collection]) 
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
axes[0,0].set_title("Sagnac FFT", fontsize = 10.5)
axes[0,1].set_title("Sagnac SNR", fontsize = 10.5)
axes[0,0].loglog(Lfreq[1:], np.abs(avgLfft[1:]) ** 2, color = "r")
axes[0,0].loglog(Lfreq[1:], np.abs(avgMfft[1:]) ** 2, color = "b")
axes[0,0].loglog(Lfreq[1:], np.abs(avgLfft_noise[1:]) ** 2, ls = "--", color = "crimson")
axes[0,0].loglog(Lfreq[1:], np.abs(avgMfft_noise[1:]) ** 2, ls = "--", color = "C0")
axes[0,0].legend(["Laser Signal", "Mic Signal", "Laser Noise", "Mic Noise"], loc = "lower left")
axes[0,1].loglog(Lfreq[1:], LSNR[1:], color = "r")
axes[0,1].loglog(Lfreq[1:], MSNR[1:], color = "b")
axes[0,1].legend(["Laser SNR", "Mic SNR"], loc = "upper right")
fig.supxlabel("Frequency (Hz)")
plt.xlim([Lfreq[1], Lfreq[-1]])


# BPD
L = CollectionTDMS(BPD_name)
L.set_collection("X")
local_detrend(L, tmin = 170e-6, tmax = 400e-6, inplace = True)
#L.apply("detrend", mode = "linear", inplace = True)
L.apply("calibrate", cal = 1/0.0002, inplace = True)
#L.apply("lowpass", cutoff = 1e6, inplace = True)
Npts = L.r/ (2 * cutoff)			# Nyquist Criterion given cutoff frequency
L.apply("bin_average", Npts = Npts, inplace = True)
Ldata = np.array([shot.time_gate(tmin = 400e-6, tmax = 470e-6)[1] for shot in L.collection])
Ldata_noise = np.array([shot.time_gate(tmin = 330e-6, tmax = 400e-6)[1] for shot in L.collection])
numcols = len(Ldata[0])
avgLfft = 0
for i in Ldata:
	avgLfft += np.fft.rfft(i) / numcols
avgLfft_noise = 0
for i in Ldata_noise:
	avgLfft_noise += np.fft.rfft(i) / numcols
Lfreq = np.fft.rfftfreq(numcols, 1/L.r)
LSNR = np.abs(avgLfft) ** 2 / np.abs(avgLfft_noise) ** 2
axes[1,0].set_title("BPD FFT", fontsize = 10.5)
axes[1,1].set_title("BPD SNR", fontsize = 10.5)
axes[1,0].loglog(Lfreq[1:], np.abs(avgLfft[1:]) ** 2, color = "r")
axes[1,0].loglog(Lfreq[1:], np.abs(avgMfft[1:]) ** 2, color = "b")
axes[1,0].loglog(Lfreq[1:], np.abs(avgLfft_noise[1:]) ** 2, ls = "--", color = "crimson")
axes[1,0].loglog(Lfreq[1:], np.abs(avgMfft_noise[1:]) ** 2, ls = "--", color = "C0")
axes[1,1].loglog(Lfreq[1:], LSNR[1:], color = "r")
axes[1,1].loglog(Lfreq[1:], MSNR[1:], color = "b")


# PD
L = CollectionTDMS(PD_name)
L.set_collection("X")
local_detrend(L, tmin = 170e-6, tmax = 400e-6, inplace = True)
#L.apply("detrend", mode = "linear", inplace = True)
L.apply("calibrate", cal = -1/0.00000024, inplace = True)
#L.apply("lowpass", cutoff = 1e6, inplace = True)
Npts = L.r/ (2 * cutoff)			# Nyquist Criterion given cutoff frequency
L.apply("bin_average", Npts = Npts, inplace = True)
Ldata = np.array([shot.time_gate(tmin = 400e-6, tmax = 470e-6)[1] for shot in L.collection])
Ldata_noise = np.array([shot.time_gate(tmin = 330e-6, tmax = 400e-6)[1] for shot in L.collection])
numcols = len(Ldata[0])
avgLfft = 0
for i in Ldata:
	avgLfft += np.fft.rfft(i) / numcols
avgLfft_noise = 0
for i in Ldata_noise:
	avgLfft_noise += np.fft.rfft(i) / numcols
Lfreq = np.fft.rfftfreq(numcols, 1/L.r)
LSNR = np.abs(avgLfft) ** 2 / np.abs(avgLfft_noise) ** 2
axes[2,0].set_title("PD FFT", fontsize = 10.5)
axes[2,1].set_title("PD SNR", fontsize = 10.5)
axes[2,0].loglog(Lfreq[1:], np.abs(avgLfft[1:]) ** 2, color = "r")
axes[2,0].loglog(Lfreq[1:], np.abs(avgMfft[1:]) ** 2, color = "b")
axes[2,0].loglog(Lfreq[1:], np.abs(avgLfft_noise[1:]) ** 2, ls = "--", color = "crimson")
axes[2,0].loglog(Lfreq[1:], np.abs(avgMfft_noise[1:]) ** 2, ls = "--", color = "C0")
axes[2,1].loglog(Lfreq[1:], LSNR[1:], color = "r")
axes[2,1].loglog(Lfreq[1:], MSNR[1:], color = "b")

plt.show()