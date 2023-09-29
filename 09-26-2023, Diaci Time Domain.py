# 26 September 2023
# Ryan Schlimme

# Integrated voltage signal

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(r"C:\Users\ryans\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS
from acoustic_entrainment import mic_response
from brownian import get_sound_speed
from scipy.special import k1e

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


def Diaci(s, dist = 0.10, temp = 20, n0 = 1.00029):
    speed = get_sound_speed(T = temp, RH = 0.5, p = 99e3)
    return -2/n0*dist/speed*s*k1e(s*dist/speed)


N = [0]
fig, ax = plt.subplots(1,1, sharey = True, sharex = True)

M_cutoff = 200e3
L_cutoff = 1e6

size = (3.375, 2.5)

# Using low pass filtering AND bin averaging w/ mic correction
L = CollectionTDMS(Sagnac_name)
M = CollectionTDMS(Sagnac_name)
L.set_collection("X")
M.set_collection("Y")
M.apply("detrend", mode = "linear", inplace = True)
L.apply("detrend", mode = "linear", inplace = True)
L.apply("calibrate", cal = -1/0.00033, inplace = True)
L.apply("lowpass", cutoff = L_cutoff, inplace = True)
M.apply("lowpass", cutoff = M_cutoff, inplace = True)
L_Npts = int(L.r/ (2 * L_cutoff))
M_Npts = int(M.r/ (2 * M_cutoff))
L.apply("bin_average", Npts = L_Npts, inplace = True)
M.apply("bin_average", Npts = M_Npts, inplace = True)
M.apply("correct", response = mic_response, recollect = True)
M.apply("shift", tau = -20e-6, inplace = True)
L.apply("correct", response = Diaci)
L.aggrigate(collection_slice = slice(1,500))
L.agg.plot(tmin = 400e-6, tmax = 470e-6, ax = ax)
M.aggrigate(collection_slice = slice(1,500))
M.agg.plot(tmin = 400e-6, tmax = 470e-6, ax = ax)

#intArray = np.zeros((len(L.agg.x), 1))

Ldata = [shot.time_gate(tmin = 400e-6, tmax = 470e-6)[1] for shot in L.collection]
Ltime = [shot.time_gate(tmin = 400e-6, tmax = 470e-6)[0] for shot in L.collection]



#DiaciSignal = np.fft.irfft(np.fft.rfft(Ldata[1])[1:]/Diaci(np.fft.rfftfreq(len(Ldata[1]))[1:]))
#DiaciSignal = [i - np.mean(DiaciSignal) for i in DiaciSignal]

# for i in range(1, len(L.agg.x)):
#     Li = L.agg
#     dx = 1/Li.r
#     intSignal = cumtrapz(Li.x, x= None, dx=dx)
#     #intArray = np.append(intArray, intSignal, axis = 1)

#intSignal = [(i) for i in DiaciSignal]

#ax.plot(intSignal)
#M.collection[1].plot(tmin = 400e-6, tmax = 470e-6, ax = ax)
plt.show()

#print(intArray)