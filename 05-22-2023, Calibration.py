# 22 May 2023
# Ryan Schlimme

# Calibration Data

import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

##### Phi Fit #####
phi_data = pd.read_csv(r"C:\Users\ryans\OneDrive\Desktop\Calibration Data\05-19-Phi.csv", sep = ",")

phi_raw = [float(i) for i in phi_data["phi"].values]
P_raw = [float(i) for i in phi_data["P"].values]

def phicos(phi, A, C, D):
	return A * np.abs(np.cos(np.radians(phi) + C)) + D

bounds = [(30, 0.5, 30), (60, 3, 200)] 
popt, pcov = curve_fit(phicos, phi_raw, P_raw, [30, 1.5, 70], bounds = bounds)
print("Fitted Parameters: ", popt)
print("Parameter Standard Errors: ", np.sqrt(np.diag(pcov)))

xvals = [float(i) for i in range(200)]
plt.plot(phi_raw, P_raw)
plt.plot(xvals, phicos(xvals, popt[0], popt[1], popt[2]))
plt.show()