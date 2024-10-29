# calculate_Qn&corr
import numpy as np

def calculate_Qn(phi_values, n):
    Qn_real = sum(np.cos(n * phi) for phi in phi_values)
    Qn_imag = sum(np.sin(n * phi) for phi in phi_values)
    
    Qn = Qn_real + 1j * Qn_imag
    return Qn

def calculate_two_particle_cumulant(phi_values, n):
    M = len(phi_values)

    Qn = calculate_Qn(phi_values, n)

    two_particle_cumulant = (abs(Qn)**2 - M) / (M * (M - 1))
    return two_particle_cumulant

def calculate_four_particle_cumulant(phi_values, n):
    M = len(phi_values)

    Qn = calculate_Qn(phi_values, n)
    Q2n = calculate_Qn(phi_values, 2 * n)

    term1 = abs(Qn)**4
    term2 = abs(Q2n)**2
    term3 = 2 * np.real(Q2n * Qn.conjugate() * Qn.conjugate())
    term4 = 2 * ( 2 * (M - 2) * (abs(Qn)**2) - M * (M - 3) )

    four_particle_cumulant = (term1 + term2 - term3 - term4) / (M * (M - 1) * (M - 2) * (M - 3))
    return four_particle_cumulant

def calculate_two_particle_correlation(phi_values, n):
    two_cumulant = calculate_two_particle_cumulant(phi_values, n)
    
    two_correlation = two_cumulant
    return two_correlation

def calculate_four_particle_correlation(phi_values, n):
    two_cumulant = calculate_two_particle_cumulant(phi_values, n)
    four_cumulant = calculate_four_particle_cumulant(phi_values, n)
    
    four_correlation = four_cumulant - 2 * ((two_cumulant)**2)
    return four_correlation