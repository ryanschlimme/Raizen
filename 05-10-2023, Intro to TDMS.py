# 10 May 2023
# Ryan Schlimme
# This test script takes a downloaded tdms file (stored in the Data folder) from LABVIEW of random laser noise and analyzes its frequency spectrum
# It provides necessary code to pull from Dr Hillberry's brownian folder of functions/modules

##### Setup Code for Brownian Functions #####
f_name = r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\iter_0.tdms" # create a variable pointing to file (change Ryan Schlimme to ryans depending on PC vs laptop)

import sys
sys.path.append(r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\brownian\src") # append path to brownian src folder containing necessary modules and functions (change Ryan Schlimme to ryans depending on PC vs laptop)

from time_series import CollectionTDMS # pull function from time_series module
C = CollectionTDMS(f_name) 		# class to apply various functions to time series data
C.set_collection("Y") # select channel Y


##### Analysis Code #####
C.apply("bin_average", Npts = 100, inplace = True) # reduces data size by factor of 100 by averaging over Npts
C.average("PSD", taumax = 20e-3) #averages and sets equal to C.freq and C.psd and C.Navg_psd (how many averages), taumax splits time domain to generate more averages

print(C.Navg_psd)

import matplotlib.pyplot as plt
plt.loglog(C.freq, C.psd)
plt.show()
