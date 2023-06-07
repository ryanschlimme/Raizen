# 22 May 2023
# Ryan Schlimme

# Calibration Data

import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

##### Phi Fit #####
phi_data = pd.read_csv(r"C:\Users\ryans\OneDrive\Desktop\Research\Calibration Data\05-19-Phi.csv", sep = ",")

phi_raw = [float(i) for i in phi_data["phi"].values]
E_raw = [float(i) for i in phi_data["E"].values]

def phicos(phi, A, C, D):
	return A * np.abs(np.cos(np.radians(phi) + C)) + D

popt, pcov = curve_fit(phicos, phi_raw, E_raw)
print("Fitted Parameters: ", popt)
print("Parameter Standard Errors: ", np.sqrt(np.diag(pcov)))

xvals = [float(i) for i in range(200)]
plt.plot(phi_raw, E_raw)
plt.plot(xvals, phicos(xvals, popt[0], popt[1], popt[2]))
plt.show()


##### Freq Max Fit #####
freq_max_data = pd.read_csv(r"C:\Users\ryans\OneDrive\Desktop\Research\Calibration Data\05-19-FreqMax.csv", sep = ",")

freq_raw = [float(i) for i in freq_max_data["f"].values[:-2]]
E_raw = [float(i) for i in freq_max_data["f"].values[:-2]]

def freq(f, a, b, c):
	return a * np.square(f) + b * f + c

popt, pcov = curve_fit(freq, freq_raw, E_raw)
xvals = [float(i) for i in range(60)]
plt.plot(freq_raw, E_raw)
plt.plot(xvals, freq(xvals, popt[0], popt[1], popt[2]))
plt.show()
	