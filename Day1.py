# 10 May 2023
# Ryan Schlimme
# This test script takes a downloaded tdms file from LABVIEW of random laser noise and analyzes its frequency spectrum

f_name = r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\iter_0.tdms"

import sys
import matplotlib.pyplot as plt
sys.path.append(r"\Users\Ryan Schlimme\OneDrive\Desktop\Research\brownian\src")
from time_series import CollectionTDMS
C = CollectionTDMS(f_name) 		# class to apply various functions to time series data


C.set_collection("Y") # select channel Y
C.apply("bin_average", Npts = 100, inplace = True) # reduces data size by factor of 100 by averaging over Npts
C.average("PSD", taumax = 20e-3) #averages and sets equal to C.freq and C.psd and C.Navg_psd (how many averages), taumax splits time domain to generate more averages

print(C.Navg_psd)
plt.loglog(C.freq, C.psd)
plt.show()
