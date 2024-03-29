# 18 July 2023
# Ryan Schlimme

# Generate averaged PSD of noise for all systems: microphone, photodiode, telescope, and Sagnac interferometer.

mic_name_index = [r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\Sagnac\Noise\iter_" + str(i) + ".tdms" for i in range(10)] # create a variable pointing to file (change Ryan Schlimme to ryans)

photo_name_index = [r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\SplitBeam\Noise\iter_" + str(i) + ".tdms" for i in range(10)] # create a variable pointing to file (change Ryan Schlimme to ryans)

tele_name_index = [r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\Telescope\Noise\iter_" + str(i) + ".tdms" for i in range(10)] # create a variable pointing to file (change Ryan Schlimme to ryans)

Sagnac_name_index = [r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\Sagnac\Noise\iter_" + str(i) + ".tdms" for i in range(10)] # create a variable pointing to file (change Ryan Schlimme to ryans)

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(r"C:\Users\ryans\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS, PSD
from acoustic_entrainment import mic_response	

N = list(range(10))

mic_PSD = 0
mic_Navg = 0
photo_PSD = 0
photo_Navg = 0
tele_PSD = 0
tele_Navg = 0
Sagnac_PSD = 0
Sagnac_Navg = 0

for n in N:
	M = CollectionTDMS(mic_name_index[n])
	M.set_collection("Y")
	mic_f, PSD = M.average("PSD", taumax = 20e-3)
	mic_Navg += M.Navg_psd
	mic_PSD += PSD / len(N)
	L = CollectionTDMS(photo_name_index[n])
	L.set_collection("X")
	photo_f, PSD = L.average("PSD", taumax = 20e-3)
	photo_PSD += PSD / len(N) 
	photo_Navg += L.Navg_psd
	L = CollectionTDMS(tele_name_index[n])
	L.set_collection("X")
	tele_f, PSD = L.average("PSD", taumax = 20e-3)
	tele_PSD += PSD / len(N) 
	tele_Navg += L.Navg_psd
	L = CollectionTDMS(Sagnac_name_index[n])
	L.set_collection("X")
	Sagnac_f, PSD = L.average("PSD", taumax = 20e-3)
	Sagnac_Navg += L.Navg_psd
	Sagnac_PSD += PSD / len(N) 

plt.loglog(mic_f, mic_PSD)
plt.loglog(photo_f, photo_PSD)
plt.loglog(tele_f, tele_PSD)
plt.loglog(Sagnac_f, Sagnac_PSD)
plt.title("PSD of Four Systems")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude (g^2/Hz)")
plt.legend(["Mic", "BPD", "PD", "Sagnac"], loc = "best")
plt.show()
	
		
