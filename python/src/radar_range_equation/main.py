import sympy
import scipy
import numpy as np
import sympy.physics.units as units
from sympy import symbols, Symbol, pprint
from sympy.physics.units import convert_to
from sympy.physics.units import speed_of_light

class vars:
    # Default physical constants and other useful defaults
    c = Symbol('c')  # speed of light (symbolic)
    c = scipy.constants.c  # speed of light (m/s)
    k = Symbol('k')  # Boltzmann constant (symbolic)
    k = scipy.constants.Boltzmann  # Boltzmann constant (J/K)
    pi = Symbol('pi')  # pi (symbolic)
    pi = scipy.constants.pi  # pi (numeric)
    pi4 = Symbol('pi4')  # 4*pi (symbolic)
    pi4 = 4 * scipy.constants.pi  # 4*pi (numeric)
    x = Symbol('x')  # generic variable for conversions (symbolic)
    f = Symbol('f')  # frequency (symbolic)
    T_0 = Symbol('T_0')  # reference temperature (symbolic)
    T_0 = 290                      # reference temperature (K) — example
    A_e = Symbol('A_e')  # effective aperture (symbolic)
    A = Symbol('A')      # antenna area (symbolic)
    D_h = Symbol('D_h')  # horizontal antenna dimension (symbolic)
    D_v = Symbol('D_v')  # vertical antenna dimension (symbolic)
    eta = Symbol('eta')  # antenna efficiency (symbolic)
    # between 0 and 1, unitless
    G_t = Symbol('G_t')  # transmit antenna gain (symbolic)
    G_t_dB = Symbol('G_t_dB')  # transmit antenna gain in dB (symbolic)
    G_r = Symbol('G_r')  # receive antenna gain (symbolic)
    G_r_dB = Symbol('G_r_dB')  # receive antenna gain in dB (symbolic)
    # figure out a way to make all wavelength variable names point to the same value
    wavelength = Symbol('lambda')  # wavelength (symbolic)
    R = Symbol('R')  # range to the fourth power (symbolic)
    P_t = Symbol('P_t')  # transmit power (symbolic)
    S_min = Symbol('S_min')  # minimum detectable signal (symbolic)
    sigma = Symbol('sigma')  # radar cross section (symbolic)

class equations:
    A_e = sympy.Eq(vars.A_e, vars.eta * vars.D_h * vars.D_v)
    wavelength = sympy.Eq(vars.wavelength, vars.c / vars.f)
    G_t = sympy.Eq(vars.G_t, 4 * sympy.pi * vars.A_e / (vars.wavelength ** 2))
    Linear_to_dB = sympy.Eq(vars.x, 10 * sympy.log(vars.x) / sympy.log(10))
    R4 = sympy.Eq(vars.R**4, (vars.P_t * vars.G_t ** 2 * vars.wavelength ** 2 * vars.sigma) / ( (4 * sympy.pi) ** 3 * vars.S_min ))
    P_t = sympy.Eq(vars.P_t, ((4 * sympy.pi) ** 3 * vars.S_min * vars.R ** 4) / (vars.G_t ** 2 * vars.wavelength ** 2 * vars.sigma))

class solve:
    def __init__():
        pass

    def A_e():
        return vars.eta * vars.D_h * vars.D_v

    def wavelength():
        return vars.c / vars.f

    def G_t():
        return 4 * vars.pi * vars.A_e / (vars.wavelength ** 2)
    
    def R4():
        value = (vars.P_t * vars.G_t ** 2 * vars.G_r * vars.wavelength ** 2 * vars.sigma) / ( (vars.pi4) ** 3 * vars.S_min )
        return value
    
    def P_t():
        value = ((vars.pi4) ** 3 * vars.S_min * vars.R ** 4) / (vars.G_t ** 2 * vars.wavelength ** 2 * vars.sigma)
        return value

def convert_to_db(value_linear):
    """
    Converts a linear scale value to decibels (dB).
    Args:
        value_linear (float): The value in linear scale.
    Returns:
        float: The value in decibels (dB).
    """
    return np.log(value_linear)/np.log(10)*10

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
    # Operate on local module objects instead of importing the package to avoid circular/package import issues
    # Use setattr to set 'lambda' since it's a reserved keyword
    vars.f = 1300e6
    pprint(f"Frequency (RRE.vars.f): {vars.f} Hz")
    pprint(f"Speed of Light (RRE.vars.c): {vars.c} m/s")
    vars.wavelength = solve.wavelength()
    vars.eta = 0.65
    vars.D_h = 12
    vars.D_v = 4
    vars.A_e = solve.A_e()
    pprint(equations.A_e)
    pprint(f"Effective Aperture (RRE.vars.A_e): {vars.A_e} m²")
    pprint(equations.wavelength)
    pprint(f"Wavelength (RRE.vars.wavelength): {vars.wavelength} m")

    pprint({equations.G_t})
    vars.G_t = solve.G_t()


    vars.G_t_dB = convert_to_db(vars.G_t)

    pprint(f"Transmit Antenna Gain (RRE.vars.G_t): {vars.G_t}")
    pprint(equations.Linear_to_dB)
    pprint(f"Transmit Antenna Gain in dB (RRE.vars.G_t_dB): {vars.G_t_dB} dB")
    pprint(equations.P_t)
    vars.S_min = 10e-13
    pprint(f"Minimum Detectable Signal (RRE.vars.S_min): {vars.S_min} W")
    vars.sigma = 1
    pprint(f"Radar Cross Section (RRE.vars.sigma): {vars.sigma} m²")
    vars.R = 200
    pprint(f"Range (RRE.vars.R): {vars.R} m")
    pprint(f"Range to the Fourth Power (RRE.vars.R**4): {vars.R**4} m⁴")
    vars.P_t = solve.P_t()
    pprint(f"Transmit Power (RRE.vars.P_t): {vars.P_t} W")
