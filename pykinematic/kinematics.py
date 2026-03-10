"""
pykinematic - Relativistic Kinematics Calculator
Calculates relativistic kinematics for light ejectiles in nuclear physics reactions
"""

import math
import os

import numpy as np


def get_mass(particle, masstable_path=None):
    """
    Get mass of a particle from the mass table.

    Parameters:
    -----------
    particle : str
        Particle identifier (e.g., "12C", "1H", "23Na")
    masstable_path : str, optional
        Path to masstable.txt file. If None, uses default location.

    Returns:
    --------
    float
        Mass in atomic mass units (amu)
    """
    if masstable_path is None:
        # Default to masstable.txt in the same directory as this module
        module_dir = os.path.dirname(os.path.abspath(__file__))
        masstable_path = os.path.join(module_dir, "masstable.txt")

    with open(masstable_path, "r") as mass_table:
        line_increment = 0
        lines = mass_table.readlines()

        for text in lines:
            if line_increment > 38:
                mass_excess_flag = False
                content_index_str = 0
                content_index_flo = 0
                content_index = 0

                for content in text.split():
                    try:
                        num = float(content)
                        content_index_flo += 1
                        content_index += 1

                        if mass_excess_flag:
                            mass_excess = num
                            mass_excess_flag = False

                    except ValueError:
                        if content_index_str == 0:
                            mass = int(prev_float)
                            element = content
                            mass_excess_flag = True

                        content_index_str += 1
                        content_index += 1

                    prev_float = num

                if str(mass) + element == particle:
                    return ((mass_excess / 1000) / 931.5) + mass

            line_increment += 1

    raise ValueError(f"Particle {particle} not found in mass table")


def qvalue(beam, target, ejectile, recoil, ex=0.0, masstable_path=None):
    """
    Calculate Q-value for a nuclear reaction.

    Parameters:
    -----------
    beam : str
        Beam particle identifier
    target : str
        Target particle identifier
    ejectile : str
        Ejectile particle identifier
    recoil : str
        Recoil particle identifier
    ex : float, optional
        Excitation energy in MeV (default: 0.0)
    masstable_path : str, optional
        Path to masstable.txt file

    Returns:
    --------
    float
        Q-value in MeV
    """
    if masstable_path is None:
        module_dir = os.path.dirname(os.path.abspath(__file__))
        masstable_path = os.path.join(module_dir, "masstable.txt")

    with open(masstable_path, "r") as mass_table:
        line_increment = 0
        lines = mass_table.readlines()

        for text in lines:
            if line_increment > 38:
                mass_excess_flag = False
                content_index_str = 0
                content_index_flo = 0
                content_index = 0

                for content in text.split():
                    try:
                        num = float(content)
                        content_index_flo += 1
                        content_index += 1

                        if mass_excess_flag:
                            mass_excess = num
                            mass_excess_flag = False

                    except ValueError:
                        if content_index_str == 0:
                            mass = int(prev_float)
                            element = content
                            mass_excess_flag = True

                        content_index_str += 1
                        content_index += 1

                    prev_float = num

                if str(mass) + element == beam:
                    mass_excess_b = mass_excess

                if str(mass) + element == target:
                    mass_excess_t = mass_excess

                if str(mass) + element == ejectile:
                    mass_excess_e = mass_excess

                if str(mass) + element == recoil:
                    mass_excess_r = mass_excess

            line_increment += 1

    q_val = ((mass_excess_b + mass_excess_t - mass_excess_r - mass_excess_e) / 1000.0)
    return q_val - float(ex)


def condition(angle, m1, m2, m3, m4, q, beam_energy):
    """
    Check if kinematics is possible for given angle.

    Parameters:
    -----------
    angle : float
        Laboratory angle in degrees
    m1 : float
        Beam mass in MeV/c^2
    m2 : float
        Target mass in MeV/c^2
    m3 : float
        Ejectile mass in MeV/c^2
    m4 : float
        Recoil mass in MeV/c^2
    q : float
        Q-value in MeV
    beam_energy : float
        Beam kinetic energy in MeV

    Returns:
    --------
    float
        Discriminant value (positive if kinematics is possible)
    """
    psi = float(angle) * (math.pi / 180.0)
    t1 = float(beam_energy)

    et = t1 + m1 + m2
    p1 = math.sqrt((t1**2) + (2 * m1 * t1))
    a = (2 * m2 * t1) + (2 * m1 * m3) + (2 * m2 * m3) + (2 * q * (m1 + m2 - m3)) - (q**2)
    b = (et**2) - ((p1**2) * (math.cos(psi)**2))

    return ((a**2) - (4 * (m3**2) * b))


def energy_calc(angle, m1, m2, m3, m4, q, beam_energy, branch='+'):
    """
    Calculate ejectile energy at a given angle.

    Parameters:
    -----------
    angle : float
        Laboratory angle in degrees
    m1 : float
        Beam mass in MeV/c^2
    m2 : float
        Target mass in MeV/c^2
    m3 : float
        Ejectile mass in MeV/c^2
    m4 : float
        Recoil mass in MeV/c^2
    q : float
        Q-value in MeV
    beam_energy : float
        Beam kinetic energy in MeV
    branch : str, optional
        '+' for forward branch, '-' for backward branch (default: '+')

    Returns:
    --------
    float
        Ejectile kinetic energy in MeV
    """
    psi = float(angle) * (math.pi / 180.0)
    t1 = float(beam_energy)

    et = t1 + m1 + m2
    p1 = math.sqrt((t1**2) + (2 * m1 * t1))
    a = (2 * m2 * t1) + (2 * m1 * m3) + (2 * m2 * m3) + (2 * q * (m1 + m2 - m3)) - (q**2)
    b = (et**2) - ((p1**2) * (math.cos(psi)**2))

    discriminant = (a**2) - (4 * (m3**2) * b)
    sqrt_term = math.sqrt(discriminant)

    if branch == '+':
        t3 = (1 / (2 * b)) * ((et * a) + (p1 * math.cos(psi) * sqrt_term)) - m3
    else:  # branch == '-'
        t3 = (1 / (2 * b)) * ((et * a) - (p1 * math.cos(psi) * sqrt_term)) - m3

    return t3


def energy_to_beta(kinetic_energy, mass):
    """
    Convert kinetic energy to beta (v/c).

    Parameters:
    -----------
    kinetic_energy : float
        Kinetic energy in MeV
    mass : float
        Rest mass in MeV/c^2.
        For a particle in an excited state, use (ground_state_mass + Ex).

    Returns:
    --------
    float
        Beta = v/c (between 0 and 1)
    """
    if kinetic_energy <= 0:
        return 0.0
    total_energy = kinetic_energy + mass
    momentum = math.sqrt(kinetic_energy**2 + 2 * mass * kinetic_energy)
    return momentum / total_energy


def excitation_energy(beam, target, ejectile, recoil, beam_energy, ejectile_energy, angle, masstable_path=None):
    """
    Calculate the excitation energy of the recoil from measured ejectile energy and angle.

    Uses energy-momentum conservation to reconstruct the invariant mass of the
    recoil nucleus, then extracts the excitation energy as Ex = M4_reconstructed - m4.

    Parameters:
    -----------
    beam : str
        Beam particle identifier (e.g., "12C")
    target : str
        Target particle identifier (e.g., "12C")
    ejectile : str
        Ejectile particle identifier (e.g., "4He")
    recoil : str
        Recoil particle identifier (e.g., "20Ne")
    beam_energy : float
        Beam kinetic energy in MeV
    ejectile_energy : float or numpy array
        Measured ejectile kinetic energy in MeV
    angle : float or numpy array
        Laboratory angle of the ejectile in degrees
    masstable_path : str, optional
        Path to masstable.txt file

    Returns:
    --------
    float or numpy array
        Excitation energy in MeV. Returns NaN for unphysical points.
    """
    # Get masses in MeV/c^2
    m1 = get_mass(beam, masstable_path) * 931.5
    m2 = get_mass(target, masstable_path) * 931.5
    m3 = get_mass(ejectile, masstable_path) * 931.5
    m4 = get_mass(recoil, masstable_path) * 931.5

    t3 = np.asarray(ejectile_energy, dtype=float)
    psi = np.asarray(angle, dtype=float) * (np.pi / 180.0)

    # Total energy and beam momentum in lab frame
    t1 = float(beam_energy)
    e_total = t1 + m1 + m2
    p1 = math.sqrt(t1**2 + 2 * m1 * t1)

    # Ejectile total energy and momentum
    e3 = t3 + m3
    p3 = np.sqrt(t3**2 + 2 * m3 * t3)

    # Recoil 4-momentum from conservation
    e4 = e_total - e3
    p4x = p1 - p3 * np.cos(psi)
    p4y = p3 * np.sin(psi)
    p4_sq = p4x**2 + p4y**2

    # Invariant mass of the recoil
    m4_reconstructed_sq = e4**2 - p4_sq
    # Use NaN for unphysical points instead of raising
    m4_reconstructed = np.where(m4_reconstructed_sq >= 0,
                                np.sqrt(np.maximum(m4_reconstructed_sq, 0)),
                                np.nan)

    return m4_reconstructed - m4


def recoil_kinematics(beam, target, ejectile, recoil, beam_energy, ejectile_energy, angle, ex=0.0, masstable_path=None):
    """
    Reconstruct the recoil kinetic energy and lab angle from measured ejectile energy and angle.

    Uses momentum conservation to reconstruct the recoil 3-momentum, then computes
    the kinetic energy assuming the recoil has the given excitation energy (which sets
    its effective invariant mass as m4 + ex).

    Parameters:
    -----------
    beam : str
        Beam particle identifier (e.g., "12C")
    target : str
        Target particle identifier (e.g., "12C")
    ejectile : str
        Ejectile particle identifier (e.g., "4He")
    recoil : str
        Recoil particle identifier (e.g., "20Ne")
    beam_energy : float
        Beam kinetic energy in MeV
    ejectile_energy : float or numpy array
        Measured ejectile kinetic energy in MeV
    angle : float or numpy array
        Ejectile laboratory angle in degrees
    ex : float, optional
        Recoil excitation energy in MeV (default: 0.0)
    masstable_path : str, optional
        Path to masstable.txt file

    Returns:
    --------
    tuple of (float or numpy array, float or numpy array)
        (recoil_energy, recoil_angle) where recoil_energy is the recoil kinetic
        energy in MeV and recoil_angle is the recoil lab angle in degrees.
        Returns NaN for unphysical points (negative kinetic energy).
    """
    m1 = get_mass(beam, masstable_path) * 931.5
    m2 = get_mass(target, masstable_path) * 931.5
    m3 = get_mass(ejectile, masstable_path) * 931.5
    m4 = get_mass(recoil, masstable_path) * 931.5

    t3 = np.asarray(ejectile_energy, dtype=float)
    psi = np.asarray(angle, dtype=float) * (np.pi / 180.0)

    # Beam momentum in lab frame
    t1 = float(beam_energy)
    p1 = math.sqrt(t1**2 + 2 * m1 * t1)

    # Ejectile momentum
    p3 = np.sqrt(np.maximum(t3**2 + 2 * m3 * t3, 0.0))

    # Recoil 3-momentum from momentum conservation
    p4x = p1 - p3 * np.cos(psi)
    p4y = p3 * np.sin(psi)
    p4_sq = p4x**2 + p4y**2

    # Effective recoil mass includes excitation energy
    m4_eff = m4 + float(ex)

    # Recoil kinetic energy: T4 = sqrt(|p4|^2 + m4_eff^2) - m4_eff
    t4 = np.sqrt(p4_sq + m4_eff**2) - m4_eff

    # Recoil lab angle in degrees
    theta4 = np.degrees(np.arctan2(p4y, p4x))

    # Mark unphysical points (negative kinetic energy) as NaN
    unphysical = t4 < 0
    t4 = np.where(unphysical, np.nan, t4)
    theta4 = np.where(unphysical, np.nan, theta4)

    # Return scalars if scalar input was given
    if t4.ndim == 0:
        return float(t4), float(theta4)
    return t4, theta4


def calculate_kinematics(beam, target, ejectile, recoil, beam_energy, ex=0.0, angle_step=1.0, masstable_path=None):
    """
    Calculate complete kinematics for all angles, including both forward and backward branches.

    Parameters:
    -----------
    beam : str
        Beam particle identifier (e.g., "12C")
    target : str
        Target particle identifier (e.g., "12C")
    ejectile : str
        Ejectile particle identifier (e.g., "1H")
    recoil : str
        Recoil particle identifier (e.g., "23Na")
    beam_energy : float
        Beam kinetic energy in MeV
    ex : float, optional
        Excitation energy in MeV (default: 0.0)
    angle_step : float, optional
        Angular step size in degrees (default: 1.0)
    masstable_path : str, optional
        Path to masstable.txt file

    Returns:
    --------
    dict with keys:
        'theta_ejectile' : list of ejectile lab angles in degrees
        'energy_ejectile' : list of ejectile kinetic energies in MeV
        'theta_recoil' : list of recoil lab angles in degrees
        'energy_recoil' : list of recoil kinetic energies in MeV
    """
    # Get masses in MeV/c^2
    m1 = get_mass(beam, masstable_path) * 931.5
    m2 = get_mass(target, masstable_path) * 931.5
    m3 = get_mass(ejectile, masstable_path) * 931.5
    m4 = get_mass(recoil, masstable_path) * 931.5

    # Calculate Q-value (already includes excitation energy)
    q = qvalue(beam, target, ejectile, recoil, ex, masstable_path)

    # Beam momentum
    p1 = math.sqrt(beam_energy**2 + 2 * m1 * beam_energy)

    theta_ej = []
    energy_ej = []
    theta_rec = []
    energy_rec = []

    def _append_point(angle_deg, t3):
        """Compute recoil quantities and append a full kinematic point."""
        # Recoil kinetic energy from energy conservation
        t4 = beam_energy + q - t3
        if t4 < 0:
            return

        # Ejectile momentum
        p3 = math.sqrt(t3**2 + 2 * m3 * t3)

        # Recoil angle from momentum conservation
        psi = angle_deg * (math.pi / 180.0)
        p4x = p1 - p3 * math.cos(psi)
        p4y = p3 * math.sin(psi)

        theta4_rad = math.atan2(p4y, p4x)
        theta4_deg = theta4_rad * (180.0 / math.pi)

        theta_ej.append(angle_deg)
        energy_ej.append(t3)
        theta_rec.append(theta4_deg)
        energy_rec.append(t4)

    # First pass: calculate forward branch (+ solution) from 0 to 180 degrees
    max_angle_forward = 0
    angle = 0.0
    curve_closes = False

    while angle < 180.0:
        if condition(angle, m1, m2, m3, m4, q, beam_energy) > 0:
            energy = energy_calc(angle, m1, m2, m3, m4, q, beam_energy, branch='+')
            # Only keep physical solutions (positive energy)
            if energy > 0:
                _append_point(angle, energy)
                max_angle_forward = angle
        else:
            # Discriminant became negative - curve closes
            curve_closes = True
            break
        angle += angle_step

    # Only calculate backward branch if the curve closes (stopped before 180 degrees)
    # If max_angle_forward is close to 180, it's an open curve
    if curve_closes and max_angle_forward < 170.0:
        # Second pass: calculate backward branch (- solution) from max angle down to 0
        # We go backwards to create a closed curve
        angle = max_angle_forward

        while angle >= 0:
            if condition(angle, m1, m2, m3, m4, q, beam_energy) > 0:
                energy = energy_calc(angle, m1, m2, m3, m4, q, beam_energy, branch='-')
                # Only keep physical solutions (positive energy and different from forward branch)
                if energy > 0:
                    # Check if this is a distinct solution from the forward branch
                    forward_energy_at_angle = None
                    for j, ang in enumerate(theta_ej):
                        if abs(ang - angle) < angle_step * 0.5:
                            forward_energy_at_angle = energy_ej[j]
                            break

                    # If significantly different from forward solution, include it
                    if forward_energy_at_angle is None or abs(energy - forward_energy_at_angle) > 0.01:
                        _append_point(angle, energy)
            angle -= angle_step

    return {
        'theta_ejectile': theta_ej,
        'energy_ejectile': energy_ej,
        'theta_recoil': theta_rec,
        'energy_recoil': energy_rec,
    }
