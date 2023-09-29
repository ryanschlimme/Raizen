from scipy.special import k1e
import numpy as np
from matplotlib import pyplot as plt

import sys
sys.path.append(r"C:\Users\ryans\OneDrive\Desktop\Research\brownian\src") 
from brownian import get_sound_speed
from time_series import CollectionTDMS

def Diaci(s, dist = 0.10, temp = 20, n0 = 1.00029):
    speed = get_sound_speed(T = temp, RH = 0.5, p = 99e3)
    return 2/n0*dist/speed*s*k1e(s*dist/speed)


Sagnac_name = r"C:\Users\ryans\OneDrive\Desktop\Research\Data\20230801\Sagnac\iter_0.tdms"
L = CollectionTDMS(Sagnac_name)
L.set_collection("X")
avgTemp = np.mean(L.T)

speed = get_sound_speed(T = avgTemp, RH = 0.5, p = 99e3)
print(speed)

dist = np.linspace(0.1, 1, 10)
s = np.linspace(1, 1e6, 10000)

fig, ax = plt.subplots(1,1)

for i in dist:
    ax.plot(s, Diaci(s, i))

plt.show()