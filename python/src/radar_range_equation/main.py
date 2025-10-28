import sympy
import scipy
import numpy as np
import sympy.physics.units as units
from sympy import symbols, Symbol, pprint
from sympy.physics.units import convert_to
from sympy.physics.units import speed_of_light

class vars:
    # Default physical constants and other useful defaults
    c = Symbol('c')                     # speed of light (symbolic)
    c = scipy.constants.c               # speed of light (m/s)
    k = Symbol('k')                     # Boltzmann constant (symbolic)
    k = scipy.constants.Boltzmann       # Boltzmann constant (J/K)
    pi = Symbol('pi')                   # pi (symbolic)
    pi = scipy.constants.pi             # pi (numeric)
    pi4 = Symbol('pi4')                 # 4*pi (symbolic)
    pi4 = 4 * scipy.constants.pi        # 4*pi (numeric)
    x = Symbol('x')                     # generic variable for conversions (symbolic)
    f = Symbol('f')                     # frequency (symbolic)
    T_0 = Symbol('T_0')                 # reference temperature (symbolic)
    T_0 = 290                           # reference temperature (K) — example
    A_e = Symbol('A_e')                 # effective aperture (symbolic)
    A = Symbol('A')                     # antenna area (symbolic)
    D_h = Symbol('D_h')                 # horizontal antenna dimension (symbolic)
    D_v = Symbol('D_v')                 # vertical antenna dimension (symbolic)
    D = Symbol('D')                     # antenna diameter dimension (symbolic)
    eta = Symbol('eta')                 # antenna efficiency (symbolic)
    G_t = Symbol('G_t')                 # transmit antenna gain (symbolic)
    G_t_dB = Symbol('G_t_dB')           # transmit antenna gain in dB (symbolic)
    G_r = Symbol('G_r')                 # receive antenna gain (symbolic)
    G_r_dB = Symbol('G_r_dB')           # receive antenna gain in dB (symbolic)
    # figure out a way to make all wavelength variable names point to the same value
    wavelength = Symbol('lambda')       # wavelength (symbolic)
    R = Symbol('R')                     # range to the fourth power (symbolic)
    P_t = Symbol('P_t')                 # transmit power (symbolic)
    S_min = Symbol('S_min')             # minimum detectable signal (symbolic)
    sigma = Symbol('sigma')             # radar cross section (symbolic)
    theta_B = Symbol('theta_B')         # beamwidth (symbolic)

class equations:
    # Create purely symbolic names (no vars.*) to keep equations symbolic for later substitution
    c_sym = Symbol('c')
    f_sym = Symbol('f')
    eta_sym = Symbol('eta')
    D_h_sym = Symbol('D_h')
    D_v_sym = Symbol('D_v')
    x_sym = Symbol('x')
    R_sym = Symbol('R')
    A_e_sym = Symbol('A_e')
    wavelength_sym = Symbol('lambda')
    G_t_sym = Symbol('G_t')
    P_t_sym = Symbol('P_t')
    sigma_sym = Symbol('sigma')
    S_min_sym = Symbol('S_min')
    R_max_sym = Symbol('R_max')
    pi_sym = Symbol('pi')
    theta_B_sym = Symbol('theta_B')

    # Symbolic equations without vars prefix
    A_e = sympy.Eq(A_e_sym, eta_sym * D_h_sym * D_v_sym)
    wavelength = sympy.Eq(wavelength_sym, c_sym / f_sym)
    G_t = sympy.Eq(G_t_sym, 4 * sympy.pi * A_e_sym / (wavelength_sym ** 2))
    Linear_to_dB = sympy.Eq(x_sym, 10 * sympy.log(x_sym) / sympy.log(10))
    R4 = sympy.Eq(R_sym**4, (P_t_sym * G_t_sym**2 * wavelength_sym**2 * sigma_sym) / ((4 * sympy.pi)**3 * S_min_sym))
    P_t = sympy.Eq(P_t_sym, ((4 * pi_sym)**3 * S_min_sym * R_sym**4) / (G_t_sym**2 * wavelength_sym**2 * sigma_sym))
    R_max = sympy.Eq(R_max_sym, sympy.Pow((P_t_sym * (G_t_sym**2) * (wavelength_sym**2) * sigma_sym) / ((4 * pi_sym)**3 * S_min_sym), sympy.Rational(1, 4), evaluate=False), evaluate=False)
    theta_B = sympy.Eq(theta_B_sym, 65 * sympy.pi / 180 * (wavelength_sym / D_h_sym))
class solve:
    def __init__():
        pass

    def A_sphere():
        return (vars.sigma / np.pi) ** (1/2)
    
    def sigma_sphere():
        return np.pi * vars.A ** 2

    def theta_B():
        return 65 * vars.pi / 180 * (vars.wavelength / vars.D_h)

    def A_e_rect():
        return vars.eta * vars.D_h * vars.D_v

    def A_e_circ():
        return vars.eta * vars.pi * (vars.D / 2) ** 2

    def wavelength():
        return vars.c / vars.f

    def G_t():
        return 4 * vars.pi * vars.A_e / (vars.wavelength ** 2)
    
    def R4():
        value = (vars.P_t * vars.G_t ** 2 * vars.G_r * vars.wavelength ** 2 * vars.sigma) / ( (vars.pi4) ** 3 * vars.S_min )
        return value
    
    def P_t():
        sym_expr = sympy.solve(equations.P_t, equations.P_t.lhs)[0]
        def _s(v):
            return v if isinstance(v, sympy.Basic) else sympy.Float(v)
        subs_map = {
            equations.S_min_sym: _s(vars.S_min),
            equations.R_sym: _s(vars.R),
            equations.G_t_sym: _s(vars.G_t),
            equations.wavelength_sym: _s(vars.wavelength),
            equations.sigma_sym: _s(vars.sigma),
            equations.pi_sym: _s(vars.pi)
        }
        value_sym = sym_expr.subs(subs_map)
        value_simpl = sympy.simplify(value_sym)
        # Return as a native Python float
        return float(value_simpl.evalf())
    
    def R_max():
        # Use the symbolic R_max equation, substitute current variable values
        # and then evaluate numerically. This preserves the symbolic derivation
        # while returning a clean numeric result (no leftover 2**(1/4) factors).
        sym_expr = sympy.solve(equations.R_max, equations.R_max.lhs)[0]

        # Helper to convert Python numbers to sympy Floats but leave sympy types alone
        def _s(v):
            return v if isinstance(v, sympy.Basic) else sympy.Float(v)

        subs_map = {
            equations.P_t_sym: _s(vars.P_t),
            equations.G_t_sym: _s(vars.G_t),
            equations.wavelength_sym: _s(vars.wavelength),
            equations.sigma_sym: _s(vars.sigma),
            equations.pi_sym: _s(vars.pi),
            equations.S_min_sym: _s(vars.S_min)
        }

        value_sym = sym_expr.subs(subs_map)
        value_simpl = sympy.simplify(value_sym)
        # Return as a native Python float
        return float(value_simpl.evalf())

def convert_to_degrees(value_radians):
    """
    Converts radians to degrees.
    Args:
        value_radians (float): The value in radians.
    Returns:
        float: The value in degrees.
    """
    return value_radians * (180 / np.pi)

def convert_to_db(value_linear):
    """
    Converts a linear scale value to decibels (dB).
    Args:
        value_linear (float): The value in linear scale.
    Returns:
        float: The value in decibels (dB).
    """
    return np.log(value_linear)/np.log(10)*10

def convert_m_to_mi(value_meters):
    """
    Converts meters to miles.
    Args:
        value_meters (float): The value in meters.
    Returns:
        float: The value in miles.
    """
    return value_meters / 1609.34

def convert_w_to_kw(value_watts):
    """
    Converts watts to kilowatts.
    Args:
        value_watts (float): The value in watts.
    Returns:
        float: The value in kilowatts.
    """
    return value_watts / 1000.0

def convert_ft_to_m(value_feet):
    """
    Converts feet to meters.
    Args:
        value_feet (float): The value in feet.
    Returns:
        float: The value in meters.
    """
    return value_feet * 0.3048

def redefine_variable(var_name, new_value):
    """
    Redefines a global variable within the 'vars' namespace.
    Args:
        var_name (str): The name of the variable to redefine (e.g., "lambda").
        new_value: The new value to assign to the variable.
    """
    setattr(vars, var_name, new_value)
# Demonstration of usage from a separate script or module:

if __name__ == '__main__':  # Only runs when the script is executed directly
    vars.f = 430*10**6  # Frequency in Hz
    vars.eta = .6
    vars.D = convert_ft_to_m(60)  # Antenna diameter in meters
    vars.R = 3.844*10**8  # Range in meters (example: distance to the Moon)
    vars.sigma = 6.64*10**11
    vars.A_e = solve.A_e_circ()
    vars.S_min = 1.5*10**-16
    vars.wavelength = solve.wavelength()
    pprint(vars.A_e)
    vars.G_t = solve.G_t()
    pprint(vars.G_t)
    pprint(equations.P_t)
    vars.P_t = solve.P_t()
    pprint(vars.P_t)

    vars.A=1.738*10**6  # Antenna radius in meters
    vars.sigma = solve.sigma_sphere()
    pprint(vars.sigma)

    vars.sigma = 6.64*10**11
    vars.A = solve.A_sphere()
    pprint(vars.A)