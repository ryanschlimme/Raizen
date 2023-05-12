# 11 May 2023
# Ryan Schlimme

# Characterizing sound profile of laser ablation as measured by 0.0068 V/Pa microphone and laser deflection detection system. 


##### Initialize TDMS File #####
f_name = r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\ablation pulse.tdms" # create a variable pointing to file (change Ryan Schlimme to ryans)
import sys													# allows access to path.append
sys.path.append(r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS 								# pull function from time_series module
L = CollectionTDMS(f_name)
L.set_collection("X")
M = CollectionTDMS(f_name)
M.set_collection("Y")

import matplotlib.pyplot as plt

##### Analyze Power Spectral Density of Laser and Microphone Outputs #####
Npts = L.r / (2*200000)
L.apply("lowpass", cutoff = 200000, inplace = True) 		# To set a frequency cutoff, we adjust Npts via the Neiquist Criterion fc = fs'/2 = 200, fs' = fs/Npts, fc = fs/2Npts, fs can be found by C.r, this downsamples our data to reduce the correlation of nearby points
# Non downsampling our data: take a signal and FFT to frequency space. Multiply by transfer function (AKA: filter, admittance, impedance). Inverse FFT to yield filtered new x(t) signal
L.apply("bin_average", Npts = Npts, inplace = True)
L.average("PSD", taumax = 20e-3)
M.apply("lowpass", cutoff = 200000, inplace = True)
M.apply("bin_average", Npts = Npts, inplace = True)
M.average("PSD", taumax = 20e-3)

plt.loglog(L.freq, L.psd, color = "red")
plt.loglog(M.freq, M.psd, color = "blue")
plt.title("PSD of Microphone and Laser Measurements")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude (g^2/Hz)")
plt.legend(["Laser", "Microphone"], loc = "best")
plt.axvline(x = 2e5, color = "black")
plt.text(2.2e5, 1e-9, "200 kHz")
plt.show()