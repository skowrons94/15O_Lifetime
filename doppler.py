import numpy as np

def calculate_doppler_beta(Eg0_val, energy_offset=0.0, angle_offset=0.0, fixed_beta=0.07):
    """
    Calculate Doppler-shifted gamma energies and corresponding beta values.
    """
    # Create theta range
    theta_vals = np.linspace(0, np.pi, 500)
    
    # Direct calculation of Doppler-shifted energy with fixed beta
    # Eg = Eg0 * sqrt(1 - beta^2) / (1 - beta * cos(theta))
    gamma_factor = np.sqrt(1 - fixed_beta**2)
    eg_vals_fixed_beta = Eg0_val * gamma_factor / (1 - fixed_beta * np.cos(theta_vals))
    
    # Now solve for beta using the measured Eg values and offset Eg0
    new_eg0_val = Eg0_val + energy_offset
    angle_offset_rad = np.deg2rad(angle_offset)
    
    # Calculate beta from: Eg = Eg0 * sqrt(1-beta^2) / (1 - beta*cos(theta))
    # Rearranging: r = Eg/Eg0, solve for beta
    r = eg_vals_fixed_beta / new_eg0_val
    cos_theta_shifted = np.cos(theta_vals + angle_offset_rad)
    
    # From the equation: r^2 * (1 - beta*cos(theta))^2 = 1 - beta^2
    # This gives: beta^2 * (r^2*cos^2(theta) - 1) - 2*r^2*beta*cos(theta) + (r^2 - 1) = 0
    # But your formula uses: beta = (r^2*cos(theta) ± sqrt(...)) / (1 + r^2*cos^2(theta))
    
    r2 = r * r
    r2_cos_theta = r2 * cos_theta_shifted
    r2_cos2_theta = r2 * cos_theta_shifted * cos_theta_shifted
    
    discriminant = np.sqrt(1 + r2_cos2_theta - r2)
    denom = 1 + r2_cos2_theta
    
    beta_vals1 = (r2_cos_theta + discriminant) / denom
    beta_vals2 = (r2_cos_theta - discriminant) / denom
    
    return beta_vals1, beta_vals2, theta_vals

def calculate_doppler_energy(Eg0_val, beta, theta_deg, ft=1):
    """
    Calculate Doppler-shifted gamma energy given beta, theta, and original energy.
    
    Uses the relativistic Doppler formula:
        Eg = Eg0 * sqrt(1 - beta^2) / (1 - beta * cos(theta))
    
    Parameters
    ----------
    Eg0_val : float or array-like
        Original (rest frame) gamma energy in keV (or any unit).
    beta : float or array-like
        Velocity of the recoiling nucleus as a fraction of c.
    theta_deg : float or array-like
        Detection angle in degrees (0 = forward, 180 = backward).
    
    Returns
    -------
    Eg_shifted : float or np.ndarray
        Doppler-shifted gamma energy in the same units as Eg0_val.
    """
    theta_rad = np.deg2rad(theta_deg)
    gamma_factor = np.sqrt(1 - ft**2 * beta**2)
    Eg_shifted = Eg0_val * gamma_factor / (1 - ft * beta * np.cos(theta_rad))
    return Eg_shifted