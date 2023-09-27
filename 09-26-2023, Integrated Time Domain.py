# 26 September 2023
# Ryan Schlimme

# Integrated voltage signal

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(r"C:\Users\ryans\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS
from acoustic_entrainment import mic_response
from scipy.integrate import cumtrapz

Sagnac_name = r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\Sagnac\iter_0.tdms"

def Phi(p):
    return np.sqrt(2) * erfinv(2*p-1)

def expected_max(N):
    mun = Phi(1-1/N)
    sigman = Phi(1-1/(N*np.e)) - mun
    return (mun + sigman*0.577)

def std_max(N):
    mun = Phi(1-1/N) + mu
    sigman = Phi(1-1/(N*np.e)) - mun
    return sigman*np.pi*np.sqrt(1/6)


def local_detrend(col, tmin = None, tmax = None, inplace = False) -> None:
	for c in col.collection:
		t, x = c.time_gate(tmin = tmin, tmax = tmax)
		x_bar = np.mean(x)
		m, b = np.polyfit(t, x, 1)
		if inplace:
			c.x = c.x - (m * c.t) - b
	return None


N = [0]
fig, ax = plt.subplots(1,1, sharey = True, sharex = True)

M_cutoff = 1e6
L_cutoff = 1e6

size = (3.375, 2.5)

# Using low pass filtering AND bin averaging w/ mic correction
L = CollectionTDMS(Sagnac_name)
L.set_collection("X")

intArray = np.zeros((len(L.collection[0].x), 1))

for i in range(1, len(L.collection)):
    Li = L.collection[i]
    dx = 1/Li.r
    intSignal = cumtrapz(Li.x, x= None, dx=dx)
    intArray = np.append(intArray, intSignal, axis = 1)

print(intArray)