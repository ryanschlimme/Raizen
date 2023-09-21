# 20 September 2023
# Ryan Schlimme

# Generates figures of PSD comparison using PSD.csv file. Format of csv file is by columns, with each column representing 
# [Micf, MicPSD, Sagnacf, SagnacPSD, BPDf, BPDPSD, PDf, PDPSD]

import pandas as pd
from matplotlib import pyplot as plt

data = pd.read_csv("PSD.csv", sep = ',')

Micf = [float(i) for i in data["Micf"].values]
MicPSD = [float(i) for i in data["MicPSD"].values]
Sagnacf = [float(i) for i in data["Sagnacf"].values]
SagnacPSD = [float(i) for i in data["SagnacPSD"].values]
BPDf = [float(i) for i in data["BPDf"].values]
BPDPSD = [float(i) for i in data["BPDPSD"].values]

#plt.figure(figsize = (3.375, 2.5))

# Plotting
plt.plot(Micf, MicPSD)
#ax.loglog(Micf, MicPSDnoise) 	# odd numbers are faded
plt.plot(Sagnacf, SagnacPSD)
plt.plot(BPDf, BPDPSD)
plt.yscale("log")
plt.legend(["Mic", "Sagnac", "BPD", "PD"])
plt.xlabel("Frequency (Hz)")
plt.title("PSD of Four Methods")

#plt.show()
plt.savefig("PSD Comparison.png")