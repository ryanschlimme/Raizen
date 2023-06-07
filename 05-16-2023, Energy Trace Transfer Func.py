# 16 May 2023
# Ryan Schlimme

# Attempt to identify transfer function using the following methodology

# Correct mic signal 0.00068 V/Pa
# Aggregate in time domain using some subset of the data n < nshot????
# Fourier transform microphone M(f) and laser L(f)
# Ensemble average the two FFTs 
# Estimate transfer function Hbar ~ Ensemble avg of (M(f)/L(f))
# Fit transfer function H(f) = A(f)e^i(phi(f)) = Re(f) + iIm(f) (fit amplitude and phase separately)
#		np.abs(Hbar) * e ** (i * np.arg(Hbar)), need to decide what np function to use for angle (arctan2 or angle maybe)
# Using fit result H(f), correct laser voltage signal


f_name_index = [r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\ene_scan_laserX_microphoneY\iter_" + str(i) + ".tdms" for i in range(6)] # create a variable pointing to file (change Ryan Schlimme to ryans)

import sys
import matplotlib.pyplot as plt
from scipy.fft import rfft, rfftfreq, irfft
from scipy.optimize import curve_fit
import numpy as np

sys.path.append(r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS								# pull function from time_series module
from acoustic_entrainment import mic_response
from brownian import logbin_func
import scipy

N = list(range(6))

def A(x, e, f, g, h):
	A = e * x ** 3 + f * x ** 2 + g * x + h
	return A

def NewResponse(s):		# From Diaci Paper
	return -2/1.003*(0.06/341)*s*scipy.special.kn(1, 0.06/341*s)*np.exp(0.06/341*s)

def Phi(f, a, b, c):
	Phi = a * f ** 2 + b * f + c
	return Phi

fig, (ax0, ax1) = plt.subplots(1, 2)

for n in N:
	f_name = f_name_index[n]
	L = CollectionTDMS(f_name)
	L.set_collection("X")
	M = CollectionTDMS(f_name)
	M.set_collection("Y")
##### Time Analysis to Compare Pulses ######
	L.apply("detrend", mode = "linear", inplace = True)
	M.apply("detrend", mode = "linear", inplace = True)
	L.apply("time_gate", tmin = 2.8e-4, tmax = 3.0e-4, inplace = True)
	M.apply("time_gate", tmin = 2.8e-4, tmax = 3.0e-4, inplace = True)
	# M.apply("shift", tau = 0.0000235, inplace = True)
	Npts = L.r / (2*200000)							# Bin_average to max frequency of 500 kHz
	L.apply("bin_average", Npts = Npts, inplace = True)
	M.apply("bin_average", Npts = Npts, inplace = True)
	M.apply("correct", response = mic_response, recollect = True)
	Mdata = np.array([shot.x for shot in M.collection]) # 2d array, each row is a different shot
	Ldata = np.array([shot.x for shot in L.collection])
	numcols = len(Mdata[0])
	numrows = len(Mdata)
	# print(divisors(numrows))

	k = 10 	# number of rows in each average	

# Initialize Dummy Paramters
	kcopy = k
	numavg = numrows / k
	index = 0
	Mavg = np.array([np.empty(numcols)])
	Lavg = np.array([np.empty(numcols)])
	k_new = k

# Loop Through Array Taking Averages
	while numavg > 0:
		avgd_row = np.array([np.mean(Mdata[index:k_new, :], axis = 0)])
		Mavg = np.append(Mavg, avgd_row, axis = 0)
		index = index + k
		k_new = k_new + k
		numavg = numavg - 1
# Remove Dummy 1st Row
	Mavg = np.delete(Mavg, 0, axis = 0)

# Initialize Dummy Paramters
	k = kcopy
	numavg = numrows / k
	index = 0
	Mavg = np.array([np.empty(numcols)])
	Lavg = np.array([np.empty(numcols)])
	k_new = k
# Loop Through Array Taking Averages
	while numavg > 0:
		avgd_row = np.array([np.mean(Ldata[index:k_new, :], axis = 0)])
		Lavg = np.append(Lavg, avgd_row, axis = 0)
		index = index + k
		k_new = k_new + k
		numavg = numavg - 1
# Remove Dummy 1st Row
	Lavg = np.delete(Lavg, 0, axis = 0)

##### FFT #####
	avgMfft = 0
	for n in Mavg:
		avgMfft += rfft(n) / numcols
	avgLfft = 0
	for n in Lavg:
		avgLfft += rfft(n) / numcols
	Lfreq = rfftfreq(len(Ldata[0]), 1/L.r)
##### Transfer Function #####
	H = avgLfft / avgMfft
# Log Bin Averaging Freq/Amplitude/Phase
	LogFreq = logbin_func(Lfreq, Npts = 25) #10-50
	LogAmp = logbin_func(np.abs(H), Npts = 25) #10-50
	LogPhi = logbin_func(np.angle(H), Npts = 25)
	# LogAmpUncertainty = logbin_func(np.abs(H), Npts = 20, func = np.std)
	ax0.loglog(LogFreq[1:], LogAmp[1:])
	# popt, pcov = curve_fit(NewResponse, LogFreq[1:], LogAmp[1:])
	# popt, pcov = curve_fit(A, LogFreq[1:], LogAmp[1:], sigma = LogAmpUncertainty[1:], absolute_sigma = True)
	# print(popt)
	# print(np.sqrt(np.diag(pcov)))
	xvals = np.linspace(1e3, 2e5, 100000)
	ax0.loglog(xvals, NewResponse(xvals))
	ax0.set_title("Amplitude Fit")
	ax1.semilogx(LogFreq[1:], LogPhi[1:])
	popt, pcov = curve_fit(Phi, LogFreq[1:], LogPhi[1:])
	ax1.semilogx(xvals, Phi(xvals, popt[0], popt[1], popt[2]))
	ax1.set_title("Phase Fit")
plt.show()