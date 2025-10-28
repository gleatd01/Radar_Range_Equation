"""Radar Range Equation Package"""

from .main import vars, \
                    equations, \
                    solve, \
                    convert_to_degrees, \
                    convert_to_db, \
                    convert_m_to_mi, \
                    convert_w_to_kw, \
                    convert_ft_to_m, \
                    redefine_variable

__all__ = ["vars",
           "equations",
           "solve",
           "convert_to_degrees",
           "convert_to_db",
           "convert_m_to_mi",
           "convert_w_to_kw",
           "convert_ft_to_m",
           "redefine_variable"
           ]