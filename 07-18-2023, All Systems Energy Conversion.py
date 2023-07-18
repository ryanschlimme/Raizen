# 18 July 2023
# Ryan Schlimme

# Analyzing three systems. Converting FLE/Phi to Pulse Energy.

import numpy as np


# Inputs flash lamp energy (FLE) in Joules and phi in degrees. Outputs the pulse energy (pulse_e) in mJ.
def pulse_e(FLE, phi):
	return 5.82 * (FLE - 70.01e-3) * np.cos(np.radians(phi) + 1.74) ** 2


##### Common Measurements #####

19J_Ablation_E = round(pulse_e(19, 82), 3)		# 110.076 mJ

min_mic_phi = [105, 118, 120, 122, 124, 126, 128, 130]
min_mic = []

for i in min_mic_phi:
	min_mic.append(round(pulse_e(12, i), 3))

print(min_mic)


##### Sagnac Interferometer #####

Sagnac_min_laser_phi = [125, 126, 128, 130, 132, 134, 136, 138, 140]
Sagnac_min_laser = []

for i in Sagnac_min_laser_phi:
	Sagnac_min_laser.append(round(pulse_e(12, i), 3))

print(Sagnac_min_laser)


##### Telescope #####

Tele_min_laser_phi = [105, 107, 109, 112, 115]
Tele_min_laser = []

for i in Tele_min_laser_phi:
	Tele_min_laser.append(round(pulse_e(12.5, i), 3))	# FLE increased to 12.5 J

print(Tele_min_laser)


##### Photodiode #####
