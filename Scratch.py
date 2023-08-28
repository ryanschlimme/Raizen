import numpy as np

def Energy(f, phi):
    return (5.82*f-70.01)*(np.cos(np.radians(phi)+1.74))**2

print(Energy(19, 145))