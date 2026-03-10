"""
pykinematic - Relativistic Kinematics Calculator

A Python package for calculating relativistic kinematics for light ejectiles
in nuclear physics reactions.

Example usage:
    from pykinematic import calculate_kinematics, energy_to_beta, get_mass

    kin = calculate_kinematics(
        beam="12C",
        target="12C",
        ejectile="1H",
        recoil="23Na",
        beam_energy=4.0,
        ex=0.0
    )
    # kin['theta_ejectile'], kin['energy_ejectile']
    # kin['theta_recoil'], kin['energy_recoil']

    # Convert ejectile energies to betas
    m3 = get_mass("1H") * 931.5  # mass in MeV/c^2
    betas = [energy_to_beta(e, m3) for e in kin['energy_ejectile']]
"""

from .kinematics import (
    calculate_kinematics,
    recoil_kinematics,
    excitation_energy,
    energy_to_beta,
    energy_calc,
    condition,
    qvalue,
    get_mass
)

__version__ = "1.0.0"
__all__ = [
    "calculate_kinematics",
    "excitation_energy",
    "energy_to_beta",
    "energy_calc",
    "condition",
    "qvalue",
    "get_mass"
]
