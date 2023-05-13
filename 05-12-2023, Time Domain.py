# 12 May 2023
# Ryan Schlimme

# Comparing aggregate time domain pulse data between microphone and laser deflection methods. Searching for proper scaling and offsets. Noticed inability of microphone to record reflected sound waves whereas laser includes sensitivity to any disturbance. Laser readout looks similar to microsphere trap voltage readout. Will try to relate to voltage readouts of acoustically-deflected laser to generate transfer function.


f_name = r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\Data\ablation pulse.tdms" # create a variable pointing to file (change Ryan Schlimme to ryans)
import sys													# allows access to path.append
sys.path.append(r"C:\Users\Ryan Schlimme\OneDrive\Desktop\Research\brownian\src") 	# append path to brownian src folder (change Ryan Schlimme to ryans)
from time_series import CollectionTDMS								# pull function from time_series module
L = CollectionTDMS(f_name)
L.set_collection("X")
M = CollectionTDMS(f_name)
M.set_collection("Y")


import matplotlib.pyplot as plt

##### Time Analysis to Compare Pulses ######
L.apply("calibrate", cal = 2.2, inplace = True)			# Vertically stretch, detrend, and horizontally shift signals to align
L.apply("detrend", mode = "linear", inplace = True)
M.apply("calibrate", cal = 0.185, inplace = True)
M.apply("detrend", mode = "linear", inplace = True)
M.apply("shift", tau = 0.00002, inplace = True)
Npts = L.r / (2*500000)							# Bin_average to max frequency of 500 kHz
L.apply("bin_average", Npts = Npts, inplace = True)
M.apply("bin_average", Npts = Npts, inplace = True)


##### Individual Pulse Plots #####
I = list(range(2,11))
for i in I:
	fig, ax = plt.subplots(1)
	Li = L.collection[i]
	Mi = M.collection[i]
	Li.plot(tmin=2.8e-4, tmax = 3e-4, ax = ax, c = "r")
	Mi.plot(tmin=2.8e-4, tmax = 3e-4, ax = ax, c = "b")
	plt.suptitle("Acoustic Signal as Recorded by Microphone and Laser Deflection")
	pulse_num = "Pulse #" + str(i)
	plt.title(pulse_num)
	x_ticks = [0.0002800, 0.0002825, 0.0002850, 0.0002875, 0.0002900, 0.0002925, 0.0002950, 0.0002975, 0.0003000]
	x_labels = [280.0, 282.5, 285.0, 287.5, 290.0, 292.5, 295.0, 297.5, 300.0]
	plt.xticks(ticks = x_ticks, labels = x_labels)
	plt.xlabel("Time (us)")
	plt.ylabel("Voltage (V)")
	plt.legend(["Laser", "Microphone"], loc = "best")
plt.show()

##### Aggregate Plots #####
fig, ax = plt.subplots(1,1)
L.aggrigate(collection_slice = slice(2, 11, 1))
L.agg.plot(tmin=2.8e-4, tmax = 3e-4, ax = ax, c = "r")
M.aggrigate(collection_slice = slice(2, 11, 1))
M.agg.plot(tmin=2.8e-4, tmax = 3e-4, ax = ax, c = "b")
plt.title("Aggregate Pulses of Laser Ablation \n as Recorded by Microphone and Laser Deflection")
x_ticks = [0.0002800, 0.0002825, 0.0002850, 0.0002875, 0.0002900, 0.0002925, 0.0002950, 0.0002975, 0.0003000]
x_labels = [280.0, 282.5, 285.0, 287.5, 290.0, 292.5, 295.0, 297.5, 300.0]
plt.xticks(ticks = x_ticks, labels = x_labels)
plt.xlabel("Time (us)")
plt.ylabel("Voltage (V)")
plt.legend(["Laser", "Microphone"], loc = "best")
plt.show()