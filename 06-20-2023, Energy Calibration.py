# 20 June 2023
# Ryan Schlimme

# Calibration Data

import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

##### Phi Fit 16J #####
phi_data = pd.read_csv(r"C:\Users\ryans\OneDrive\Desktop\Research\Calibration Data\phi_scan_16J.csv", sep = ",")

phi_raw = [float(np.radians(i)) for i in phi_data["Phi (degrees)"].values]
E_raw = [float(i) for i in phi_data["Energy (mJ)"].values]

def phicos(phi, A, C):
	return A * (np.cos(phi + C))**2

guess = [20, 3]
bounds = ((20, 0), (30, np.pi))
popt, pcov = curve_fit(phicos, phi_raw, E_raw, p0 = guess, bounds = bounds)
print("Phi Scan for 16J")
print("Fitted Parameters: ", popt)
print("Parameter Standard Errors: ", np.sqrt(np.diag(pcov)))

xvals = [float(i) for i in np.linspace(0.3,2.8,200)]
plt.plot(phi_raw, E_raw)
plt.plot(xvals, phicos(xvals, popt[0], popt[1]))
plt.legend(["Raw Data", "Fit"])
plt.title("Phi Fit 16J Flash Lamp")
plt.xlabel("Phi (in radians)")
plt.ylabel("Measured Energy (in mJ)")


##### Phi Fit 18J #####
phi_data = pd.read_csv(r"C:\Users\ryans\OneDrive\Desktop\Research\Calibration Data\phi_scan_18J.csv", sep = ",")

phi_raw = [float(np.radians(i)) for i in phi_data["Phi (degrees)"].values]
E_raw = [float(i) for i in phi_data["Energy (mJ)"].values]

guess = [35, 3]
bounds = ((30, 0), (40, np.pi))
popt, pcov = curve_fit(phicos, phi_raw, E_raw, p0 = guess, bounds = bounds)
print("Phi Scan for 18J")
print("Fitted Parameters: ", popt)
print("Parameter Standard Errors: ", np.sqrt(np.diag(pcov)))

xvals = [float(i) for i in np.linspace(0.3,2.8,200)]
plt.plot(phi_raw, E_raw)
plt.plot(xvals, phicos(xvals, popt[0], popt[1]))
plt.legend(["Raw Data 16J", "Fit 16J", "Raw Data 18", "Fit 18J"])
plt.title("Phi Fit")
plt.xlabel("Phi (in radians)")
plt.ylabel("Measured Energy (in mJ)")
plt.show()


##### Energy Fit  #####
E_data = pd.read_csv(r"C:\Users\ryans\OneDrive\Desktop\Research\Calibration Data\energy_phi91.csv", sep = ",")

E_in_raw = [float(i) for i in E_data["EnergyIn (J)"].values]
E_out_raw = [float(i) for i in E_data["EnergyOut (mJ)"].values]
T_rate_raw = [float(i) for i in E_data["TransmissionRate"].values]

def Energy(E_in, A, B):
	return A * np.array(E_in) + B

popt, pcov = curve_fit(Energy, E_in_raw, E_out_raw)
mean = np.mean(T_rate_raw)
sd = np.std(T_rate_raw)

print("Energy Scan")
print("Fitted Parameters: ", popt)
print("Parameter Standard Errors: ", np.sqrt(np.diag(pcov)))
print()
print("Transmission Rate Statistics")
print("Mean: ", mean, " Standard Deviation: ", sd)

xvals = [float(i) for i in range(12, 20)]
plt.plot(E_in_raw, E_out_raw)
plt.plot(xvals, Energy(xvals, popt[0], popt[1]))
plt.legend(["Raw Data", "Fit"])
plt.title("Energy Scan")
plt.xlabel("Flash Lamp Energy (in J)")
plt.ylabel("Measured Energy (in mJ)")
plt.show()	