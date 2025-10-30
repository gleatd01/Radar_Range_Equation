"""Radar Range Equation Package"""

from .main import vars, \
                    equations, \
                    solve, \
                    convert, \
                    redefine_variable, \
                    v_phi, \
                    calculate_Theta, \
                    v_phi_full, \
                    estimate_phi_hat, \
                    sigma_phi_amplitude, \
                    sigma_phi_phase, \
                    sigma_phi_time, \
                    db_to_linear

__all__ = ["vars",
           "equations",
           "solve",
           "convert",
           "redefine_variable",
           "v_phi",
           "calculate_Theta",
           "v_phi_full",
           "estimate_phi_hat",
           "sigma_phi_amplitude",
           "sigma_phi_phase",
           "sigma_phi_time",
           "db_to_linear"
           ]