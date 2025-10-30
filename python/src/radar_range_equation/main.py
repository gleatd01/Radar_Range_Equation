import sympy
import scipy
import numpy as np
import math
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
    R_max = Symbol('R_max')             # maximum range (symbolic)
    R_un = Symbol('R_un')                   # unnormalized range (symbolic)
    f_m = Symbol('f_m')                   # modulation frequency (symbolic)
    f_bu = Symbol('f_bu')                 # upper band frequency (symbolic)
    f_bd = Symbol('f_bd')                 # lower band frequency (symbolic)
    f_r = Symbol('f_r')                   # radar operating frequency (symbolic)
    f_d = Symbol('f_d')                   # frequency deviation (symbolic)
    deltaf = Symbol('Delta f')            # frequency difference (symbolic, with space)
    v = Symbol('v')                     # velocity (symbolic)
    #R(hat)_max
    latex = False  # Set to True for LaTeX-style variable names
    if latex == True:
        R_hat_max = Symbol(r"\hat{R}_{max}")  # normalized maximum range (symbolic/latex)
    else:
        R_hat_max = Symbol("R\u0302_max")     # normalized maximum range (symbolic)

v = vars()  # Create a global instance of vars for easy access

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
    R_un_sym = Symbol('R_un')
    f_m_sym = Symbol('f_m')
    f_bu_sym = Symbol('f_bu')
    f_bd_sym = Symbol('f_bd')
    f_r_sym = Symbol('f_r')
    f_d_sym = Symbol('f_d')
    f_0_sym = Symbol('f_0')
    deltaf_sym = Symbol('Delta f')
    v_sym = Symbol('v')
    # Angle estimation symbols
    phi_sym = Symbol('phi')
    phi_s_sym = Symbol('phi_s')
    Theta_sym = Symbol('Theta')
    Delta_sym = Symbol('Delta')
    Sigma_sym = Symbol('Sigma')
    S_N_sym = Symbol('S_N')
    d_sym = Symbol('d')
    B_sym = Symbol('B')

    # Symbolic equations without vars prefix
    A_e = sympy.Eq(A_e_sym, eta_sym * D_h_sym * D_v_sym)
    wavelength = sympy.Eq(wavelength_sym, c_sym / f_sym)
    G_t = sympy.Eq(G_t_sym, 4 * sympy.pi * A_e_sym / (wavelength_sym ** 2))
    Linear_to_dB = sympy.Eq(x_sym, 10 * sympy.log(x_sym) / sympy.log(10))
    R4 = sympy.Eq(R_sym**4, (P_t_sym * G_t_sym**2 * wavelength_sym**2 * sigma_sym) / ((4 * sympy.pi)**3 * S_min_sym))
    P_t = sympy.Eq(P_t_sym, ((4 * pi_sym)**3 * S_min_sym * R_sym**4) / (G_t_sym**2 * wavelength_sym**2 * sigma_sym))
    R_max = sympy.Eq(R_max_sym, sympy.Pow((P_t_sym * (G_t_sym**2) * (wavelength_sym**2) * sigma_sym) / ((4 * pi_sym)**3 * S_min_sym), sympy.Rational(1, 4), evaluate=False), evaluate=False)
    theta_B = sympy.Eq(theta_B_sym, 65 * sympy.pi / 180 * (wavelength_sym / D_h_sym))
    R = sympy.Eq(R_sym, (c_sym*f_r_sym)/(4*f_m_sym*deltaf_sym))
    v = sympy.Eq(v_sym, -(c_sym/f_sym)*(f_d_sym/2))
    f_m = sympy.Eq(f_m_sym, c_sym/(2*R_un_sym))
    f_0 = sympy.Eq(f_0_sym, 2*f_m_sym*deltaf_sym)
    f_r = sympy.Eq(f_r_sym, .5*(f_bu_sym+f_bd_sym))
    f_d = sympy.Eq(f_d_sym, .5*(f_bu_sym-f_bd_sym))
    # Angle estimation equations
    Theta = sympy.Eq(Theta_sym, (4 * sympy.log(2)) / (theta_B_sym**2))
    v_phi = sympy.Eq(v_sym, sympy.exp(-Theta_sym * (phi_sym - phi_s_sym)**2))
    v_phi_full = sympy.Eq(v_sym, sympy.exp((4 * sympy.log(2) * (phi_sym - phi_s_sym)**2) / (theta_B_sym**2)))
    phi_hat = sympy.Eq(phi_sym, (Delta_sym / Sigma_sym) * (theta_B_sym**2 / (8 * sympy.log(2) * phi_s_sym)))
    sigma_phi_amp = sympy.Eq(sigma_sym, (theta_B_sym**2 * sympy.sqrt(1 / S_N_sym)) / (8 * sympy.sqrt(2) * phi_s_sym * sympy.log(2)))
    sigma_phi_phase = sympy.Eq(sigma_sym, (wavelength_sym / (2 * sympy.pi * d_sym)) * sympy.sqrt(1 / S_N_sym))
    sigma_phi_time = sympy.Eq(sigma_sym, c_sym / (d_sym * B_sym))
    dB_to_linear = sympy.Eq(x_sym, 10**(x_sym / 10))

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
    
    def R():
        sym_expr = sympy.solve(equations.R, equations.R.lhs)[0]

        # Helper to convert Python numbers to sympy Floats but leave sympy types alone
        def _s(v):
            return v if isinstance(v, sympy.Basic) else sympy.Float(v)

        subs_map = {
            equations.c_sym: _s(vars.c),
            equations.f_r_sym: _s(vars.f_r),
            equations.f_m_sym: _s(vars.f_m),
            equations.deltaf_sym: _s(vars.deltaf)
        }

        value_sym = sym_expr.subs(subs_map)
        value_simpl = sympy.simplify(value_sym)
        # Return as a native Python float
        return float(value_simpl.evalf())
    
    def v():
        sym_expr = sympy.solve(equations.v, equations.v.lhs)[0]

        # Helper to convert Python numbers to sympy Floats but leave sympy types alone
        def _s(v):
            return v if isinstance(v, sympy.Basic) else sympy.Float(v)

        subs_map = {
            equations.c_sym: _s(vars.c),
            equations.f_sym: _s(vars.f),
            equations.f_d_sym: _s(vars.f_d)
        }

        value_sym = sym_expr.subs(subs_map)
        value_simpl = sympy.simplify(value_sym)
        # Return as a native Python float
        return float(value_simpl.evalf())
    
    def f_m():
        sym_expr = sympy.solve(equations.f_m, equations.f_m.lhs)[0]

        # Helper to convert Python numbers to sympy Floats but leave sympy types alone
        def _s(v):
            return v if isinstance(v, sympy.Basic) else sympy.Float(v)

        subs_map = {
            equations.c_sym: _s(vars.c),
            equations.R_un_sym: _s(vars.R_un)
        }

        value_sym = sym_expr.subs(subs_map)
        value_simpl = sympy.simplify(value_sym)
        # Return as a native Python float
        return float(value_simpl.evalf())
    
    def f_r():
        sym_expr = sympy.solve(equations.f_r, equations.f_r.lhs)[0]

        # Helper to convert Python numbers to sympy Floats but leave sympy types alone
        def _s(v):
            return v if isinstance(v, sympy.Basic) else sympy.Float(v)

        subs_map = {
            equations.f_bu_sym: _s(vars.f_bu),
            equations.f_bd_sym: _s(vars.f_bd)
        }

        value_sym = sym_expr.subs(subs_map)
        value_simpl = sympy.simplify(value_sym)
        # Return as a native Python float
        return float(value_simpl.evalf())
    
    def f_d():
        sym_expr = sympy.solve(equations.f_d, equations.f_d.lhs)[0]

        # Helper to convert Python numbers to sympy Floats but leave sympy types alone
        def _s(v):
            return v if isinstance(v, sympy.Basic) else sympy.Float(v)

        subs_map = {
            equations.f_bu_sym: _s(vars.f_bu),
            equations.f_bd_sym: _s(vars.f_bd)
        }

        value_sym = sym_expr.subs(subs_map)
        value_simpl = sympy.simplify(value_sym)
        # Return as a native Python float
        return float(value_simpl.evalf())
    
    def f_0():
        sym_expr = sympy.solve(equations.f_0, equations.f_0.lhs)[0]

        # Helper to convert Python numbers to sympy Floats but leave sympy types alone
        def _s(v):
            return v if isinstance(v, sympy.Basic) else sympy.Float(v)

        subs_map = {
            equations.f_m_sym: _s(vars.f_m),
            equations.deltaf_sym: _s(vars.deltaf)
        }

        value_sym = sym_expr.subs(subs_map)
        value_simpl = sympy.simplify(value_sym)
        # Return as a native Python float
        return float(value_simpl.evalf())
    
    # Angle estimation solve methods
    def calculate_Theta(theta_B):
        """Calculates the Theta parameter from the 3-dB beamwidth.
        
        Args:
            theta_B: 3-dB beamwidth in degrees
            
        Returns:
            Theta parameter
        """
        return (4 * math.log(2)) / (theta_B**2)
    
    def v_phi(phi, phi_s, Theta):
        """Calculates the Gaussian beam approximation.
        
        Args:
            phi: Angle in degrees
            phi_s: Center angle in degrees
            Theta: Beam parameter
            
        Returns:
            Gaussian beam approximation value
        """
        return math.exp(-Theta * (phi - phi_s)**2)
    
    def v_phi_full(phi, phi_s, theta_B):
        """Calculates v(phi) using theta_B directly.
        
        Args:
            phi: Angle in degrees
            phi_s: Center angle in degrees
            theta_B: 3-dB beamwidth in degrees
            
        Returns:
            Beam approximation value
        """
        # Note: -4 * log(0.5) = 4 * log(2) (mathematically equivalent)
        # Using 4 * log(2) for better performance
        return math.exp((4 * math.log(2) * (phi - phi_s)**2) / (theta_B**2))
    
    def estimate_phi_hat(Delta, Sigma, theta_B, phi_s):
        """Calculates the linear processor angle estimate.
        
        Args:
            Delta: Delta signal
            Sigma: Sigma signal
            theta_B: 3-dB beamwidth in degrees
            phi_s: Center angle in degrees
            
        Returns:
            Angle estimate in degrees
        """
        return (Delta / Sigma) * (theta_B**2 / (8 * math.log(2) * phi_s))
    
    def sigma_phi_amplitude(theta_B, S_N, phi_s):
        """Calculates the angle standard deviation for amplitude comparison.
        
        Args:
            theta_B: 3-dB beamwidth in degrees
            S_N: Signal-to-Noise ratio (linear, not dB)
            phi_s: Center angle in degrees
            
        Returns:
            Angle standard deviation in degrees
        """
        return (theta_B**2 * math.sqrt(1 / S_N)) / (8 * math.sqrt(2) * phi_s * math.log(2))
    
    def sigma_phi_phase(lambda_, d, S_N):
        """Calculates the angle standard deviation for phase comparison.
        
        Args:
            lambda_: Wavelength in meters
            d: Antenna separation in meters
            S_N: Signal-to-Noise ratio (linear, not dB)
            
        Returns:
            Angle standard deviation in radians
        """
        return (lambda_ / (2 * math.pi * d)) * math.sqrt(1 / S_N)
    
    def sigma_phi_time(c, d, B):
        """Calculates the angle standard deviation for time comparison.
        
        Args:
            c: Speed of light in m/s
            d: Antenna separation in meters
            B: Bandwidth in Hz
            
        Returns:
            Angle standard deviation in radians
        """
        return c / (d * B)
    
    def db_to_linear(value_db):
        """Converts SNR from dB to linear.
        
        Args:
            value_db: Value in dB
            
        Returns:
            Linear value
        """
        return 10**(value_db / 10)

class convert: # add alias con for convenience
    pass
    
    def rad_to_deg(value_radians):
        return value_radians * (180 / np.pi)
    
    def deg_to_rad(value_degrees):
        return value_degrees * (np.pi / 180)
    
    def lin_to_db(value_linear):
        return np.log(value_linear)/np.log(10)*10
    
    def m_to_mi(value_meters):
        return value_meters / 1609.34
    
    def w_to(value_w, target_unit):
        conversion_factors = {
            'kw': 1 / 1000.0,
            'mw': 1 * 10**6
        }
        return value_w * conversion_factors.get(target_unit.lower(), 1)

    def w_from(value, source_unit):
        conversion_factors = {
            'kw': 1 / 1000.0,
            'mw': 1 * 10**6
        }           
        return value / conversion_factors.get(source_unit.lower(), 1)
    
    def ft_to_m(value_feet):
        return value_feet * 0.3048

    def hz_to(value_hz, target_unit):
        conversion_factors = {
            'khz': 1 * 10**3,
            'mhz': 1 * 10**6,
            'ghz': 1 * 10**9
        }
        return value_hz / conversion_factors.get(target_unit.lower(), 1)
    
    def hz_from(value, source_unit):
        conversion_factors = {
            'khz': 1 * 10**3,
            'mhz': 1 * 10**6,
            'ghz': 1 * 10**9
        }
        return value * conversion_factors.get(source_unit.lower(), 1)
    
    def m_to(value_meters, target_unit):
        conversion_factors = {
            'km': 1 / 1000.0
        }
        return value_meters * conversion_factors.get(target_unit.lower(), 1)

    def m_from(value, source_unit):
        conversion_factors = {
            'km': 1 / 1000.0
        }
        return value / conversion_factors.get(source_unit.lower(), 1)
    
con = convert()  # alias for convenience

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
    v = vars()
    v.f_bu = Symbol('f_bu')
    pprint(v.f_bu)
    v.f_bu = convert.hz_from(1510, 'mhz')
    pprint(convert.hz_to(v.f_bu, 'mhz'))
    v.f_bd = Symbol('f_bd')
    pprint(v.f_bd)
    v.f_bd = convert.hz_from(1490, 'mhz')
    pprint(convert.hz_to(v.f_bd, 'mhz'))
    v.f_r = Symbol('f_r')
    pprint(v.f_r)
    v.f_r=.5*(v.f_bu+v.f_bd)
    pprint(convert.hz_to(v.f_r, 'mhz'))
    v.f_d = Symbol('f_d')
    pprint(v.f_d)
    v.f_d = .5*(v.f_bd-v.f_bu)
    pprint(convert.hz_to(v.f_d, 'mhz'))
    # plain word with space
    v.deltaf = Symbol('delta f')
    pprint(v.deltaf)