# 11 May 2023
# Ryan Schlimme

# Characterizing sound profile of laser ablation as measured by 0.006 V/Pa microphone and laser deflection detection system. 


##### Initialize TDMS File #####
f_name = r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\ablation pulse.tdms" # create a variable pointing to file (change Ryan Schlimme to ryans)
import sys													# allows access to path.append
sys.path.append(r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS 								# pull function from time_series module
L = CollectionTDMS(f_name)
L.set_collection("X")
M = CollectionTDMS(f_name)
M.set_collection("Y")

##### Analyze Power Spectral Density of Laser and Microphone Outputs #####
L.apply("bin_average", Npts = 50, inplace = True)
L.average("PSD", taumax = 1e-3)
M.apply("bin_average", Npts = 50, inplace = True)
M.average("PSD", taumax = 1e-3)

import matplotlib.pyplot as plt
plt.loglog(L.freq, L.psd, color = "red")
plt.loglog(M.freq, M.psd, color = "blue")
plt.title("PSD of Microphone and Laser Measurements")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude (g^2/Hz)")
plt.legend(["Laser", "Microphone"], loc = "best")
plt.axvline(x = 20000, color = "black")
plt.text(23000, 1e-8, "20 kHz")
plt.ylim([1e-12,1])
plt.show()