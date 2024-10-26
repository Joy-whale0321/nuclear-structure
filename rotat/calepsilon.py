# calculate_Qn.py

import numpy as np

def calculate_Qn(phi_values, n):
    Qn_real = sum(np.cos(n * phi) for phi in phi_values)
    Qn_imag = sum(np.sin(n * phi) for phi in phi_values)
    
    Qn = Qn_real + 1j * Qn_imag
    return Qn
