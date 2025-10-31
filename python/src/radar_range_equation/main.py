import sympy
import scipy
import numpy as np
import math
import sympy.physics.units as units
from sympy import symbols, Symbol, pprint, exp, sqrt, log
from sympy.physics.units import convert_to
from sympy.physics.units import speed_of_light

class vars:
    # =========================================================================
    # BASE/COMMON PHYSICAL CONSTANTS AND VARIABLES
    # =========================================================================
    
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
    wavelength = Symbol('lambda')       # wavelength (symbolic)
    v = Symbol('v')                     # velocity (symbolic)
    
    # Antenna parameters
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
    theta_B = Symbol('theta_B')         # beamwidth (symbolic)
    
    # Radar equation parameters
    R = Symbol('R')                     # range (symbolic)
    R_max = Symbol('R_max')             # maximum range (symbolic)
    P_t = Symbol('P_t')                 # transmit power (symbolic)
    S_min = Symbol('S_min')             # minimum detectable signal (symbolic)
    sigma = Symbol('sigma')             # radar cross section (symbolic)

    # =========================================================================
    # TOPIC 07: DOPPLER CW RADAR VARS
    # =========================================================================
    
    f_doppler = Symbol('f_doppler')     # Doppler frequency shift (symbolic)
    f_if = Symbol('f_if')               # Intermediate frequency (symbolic)
    f_obs = Symbol('f_obs')             # Observed frequency at receiver (symbolic)
    T_cpi = Symbol('T_cpi')             # Coherent processing time (symbolic)
    delta_v = Symbol('delta_v')         # Velocity resolution (symbolic)

    # =========================================================================
    # TOPIC 08: CWFM RADAR VARS
    # =========================================================================
    
    R_un = Symbol('R_un')               # unnormalized range (symbolic)
    f_m = Symbol('f_m')                 # modulation frequency (symbolic)
    f_bu = Symbol('f_bu')               # upper band frequency (symbolic)
    f_bd = Symbol('f_bd')               # lower band frequency (symbolic)
    f_r = Symbol('f_r')                 # radar operating frequency (symbolic)
    f_d = Symbol('f_d')                 # frequency deviation (symbolic)
    deltaf = Symbol('Delta f')          # frequency difference (symbolic, with space)

    # =========================================================================
    # TOPIC 09: PULSED RADAR VARS
    # =========================================================================
    
    f_p = Symbol('f_p')                 # Pulse Repetition Frequency (PRF) (symbolic)
    T_p = Symbol('T_p')                 # Pulse Repetition Interval (PRI) (symbolic)
    t_delay = Symbol('t_delay')         # Round-trip time delay (symbolic)
    duty_cycle = Symbol('duty_cycle')   # Duty cycle (symbolic)
    tau = Symbol('tau')                 # Pulse width (symbolic)
    n_p = Symbol('n_p')                 # Number of pulses integrated (symbolic)
    t_scan = Symbol('t_scan')           # Antenna scan time (symbolic)
    S_N_1 = Symbol('S_N_1')             # Single pulse SNR (linear) (symbolic)
    S_N_1_dB = Symbol('S_N_1_dB')       # Single pulse SNR (dB) (symbolic)
    S_N_n = Symbol('S_N_n')             # Integrated SNR (linear) (symbolic)
    E_i = Symbol('E_i')                 # Integration efficiency (symbolic)

    # =========================================================================
    # TOPIC 10: DIRECTION FINDING VARS
    # =========================================================================
    
    phi = Symbol('phi')                 # Angle (symbolic)
    phi_s = Symbol('phi_s')             # Squint angle (symbolic)
    Theta = Symbol('Theta')             # Gaussian beam parameter (symbolic)
    v_phi = Symbol('v_phi')             # Gaussian beam voltage (symbolic)
    Delta = Symbol('Delta')             # Difference signal (symbolic)
    Sigma = Symbol('Sigma')             # Sum signal (symbolic)
    phi_hat = Symbol('phi_hat')         # Angle estimate (symbolic)
    sigma_phi = Symbol('sigma_phi')     # Angle standard deviation (symbolic)
    S_N = Symbol('S_N')                 # Signal-to-Noise ratio (linear) (symbolic)
    S_N_dB = Symbol('S_N_dB')           # Signal-to-Noise ratio (dB) (symbolic)
    d = Symbol('d')                     # Element separation (symbolic)
    B = Symbol('B')                     # Bandwidth (symbolic)

    # =========================================================================
    # TOPIC 11: PULSE COMPRESSION VARS
    # =========================================================================
    
    delta_r = Symbol('delta_r')         # Range resolution (symbolic)
    gamma = Symbol('gamma')             # Chirp rate (symbolic)
    PCR = Symbol('PCR')                 # Pulse Compression Ratio (symbolic)
    R_offset = Symbol('R_offset')       # Range offset from reference (symbolic)
    f_range_tone = Symbol('f_range_tone') # IF frequency from dechirp (symbolic)

    # =========================================================================
    # SPECIAL VARIABLES
    # =========================================================================
    
    latex = False  # Set to True for LaTeX-style variable names
    if latex == True:
        R_hat_max = Symbol(r"\hat{R}_{max}")  # normalized maximum range (symbolic/latex)
    else:
        R_hat_max = Symbol("R\u0302_max")    # normalized maximum range (symbolic)

v = vars()  # Create a global instance of vars for easy access

class equations:
    # --- Base/Common Symbols ---
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
    v_sym = Symbol('v')

    # --- Topic 07: Doppler CW Radar Symbols ---
    f_doppler_sym = Symbol('f_doppler')
    f_if_sym = Symbol('f_if')
    f_obs_sym = Symbol('f_obs')
    T_cpi_sym = Symbol('T_cpi')
    delta_v_sym = Symbol('delta_v')

    # --- Topic 08: CWFM Radar Symbols ---
    R_un_sym = Symbol('R_un')
    f_m_sym = Symbol('f_m')
    f_bu_sym = Symbol('f_bu')
    f_bd_sym = Symbol('f_bd')
    f_r_sym = Symbol('f_r')
    f_d_sym = Symbol('f_d')
    f_0_sym = Symbol('f_0')
    deltaf_sym = Symbol('Delta f')
    
    # --- Topic 09: Pulsed Radar Symbols ---
    f_p_sym = Symbol('f_p')
    T_p_sym = Symbol('T_p')
    t_delay_sym = Symbol('t_delay')
    duty_cycle_sym = Symbol('duty_cycle')
    tau_sym = Symbol('tau')
    n_p_sym = Symbol('n_p')
    t_scan_sym = Symbol('t_scan')
    S_N_1_sym = Symbol('S_N_1')
    S_N_n_sym = Symbol('S_N_n')
    E_i_sym = Symbol('E_i')

    # --- Topic 10: Direction Finding Symbols ---
    phi_sym = Symbol('phi')
    phi_s_sym = Symbol('phi_s')
    Theta_sym = Symbol('Theta')
    v_phi_sym = Symbol('v_phi')
    Delta_sym = Symbol('Delta')
    Sigma_sym = Symbol('Sigma')
    phi_hat_sym = Symbol('phi_hat')
    sigma_phi_sym = Symbol('sigma_phi')
    S_N_sym = Symbol('S_N')
    S_N_dB_sym = Symbol('S_N_dB')
    d_sym = Symbol('d')
    B_sym = Symbol('B')

    # --- Topic 11: Pulse Compression Symbols ---
    delta_r_sym = Symbol('delta_r')
    gamma_sym = Symbol('gamma')
    PCR_sym = Symbol('PCR')
    R_offset_sym = Symbol('R_offset')
    f_range_tone_sym = Symbol('f_range_tone')

    # =========================================================================
    # EQUATIONS
    # =========================================================================

    # --- Base/Common Equations ---
    A_e = sympy.Eq(A_e_sym, eta_sym * D_h_sym * D_v_sym)
    wavelength = sympy.Eq(wavelength_sym, c_sym / f_sym)
    G_t = sympy.Eq(G_t_sym, 4 * sympy.pi * A_e_sym / (wavelength_sym ** 2))
    Linear_to_dB = sympy.Eq(x_sym, 10 * sympy.log(x_sym) / sympy.log(10))
    R4 = sympy.Eq(R_sym**4, (P_t_sym * G_t_sym**2 * wavelength_sym**2 * sigma_sym) / ((4 * sympy.pi)**3 * S_min_sym))
    P_t = sympy.Eq(P_t_sym, ((4 * pi_sym)**3 * S_min_sym * R_sym**4) / (G_t_sym**2 * wavelength_sym**2 * sigma_sym))
    R_max = sympy.Eq(R_max_sym, sympy.Pow((P_t_sym * (G_t_sym**2) * (wavelength_sym**2) * sigma_sym) / ((4 * pi_sym)**3 * S_min_sym), sympy.Rational(1, 4), evaluate=False), evaluate=False)
    theta_B = sympy.Eq(theta_B_sym, 65 * sympy.pi / 180 * (wavelength_sym / D_h_sym)) # Gaussian approx

    # --- Topic 07: Doppler CW Radar Equations ---
    eq_f_doppler = sympy.Eq(f_doppler_sym, -2 * v_sym / wavelength_sym)
    eq_v_from_doppler = sympy.Eq(v_sym, -wavelength_sym * f_doppler_sym / 2)
    eq_f_obs_if = sympy.Eq(f_obs_sym, f_if_sym + f_doppler_sym)
    eq_delta_v = sympy.Eq(delta_v_sym, wavelength_sym / (2 * T_cpi_sym))
    
    # --- Topic 08: CWFM Radar Equations ---
    R_cwfm = sympy.Eq(R_sym, (c_sym*f_r_sym)/(4*f_m_sym*deltaf_sym))
    v_cwfm = sympy.Eq(v_sym, -(c_sym/f_sym)*(f_d_sym/2))
    f_m_cwfm = sympy.Eq(f_m_sym, c_sym/(2*R_un_sym))
    f_0_cwfm = sympy.Eq(f_0_sym, 2*f_m_sym*deltaf_sym)
    f_r_cwfm = sympy.Eq(f_r_sym, .5*(f_bu_sym+f_bd_sym))
    f_d_cwfm = sympy.Eq(f_d_sym, .5*(f_bu_sym-f_bd_sym))

    # --- Topic 09: Pulsed Radar & Range Ambiguity Equations ---
    eq_R_un_from_fp = sympy.Eq(R_un_sym, c_sym / (2 * f_p_sym))
    eq_fp_from_R_un = sympy.Eq(f_p_sym, c_sym / (2 * R_un_sym))
    eq_T_p = sympy.Eq(T_p_sym, 1 / f_p_sym)
    eq_R_from_time = sympy.Eq(R_sym, c_sym * t_delay_sym / 2)
    eq_tau_from_duty = sympy.Eq(tau_sym, T_p_sym * duty_cycle_sym)
    eq_n_p = sympy.Eq(n_p_sym, f_p_sym * t_scan_sym * (theta_B_sym / (2 * pi_sym)))
    eq_S_N_n_coherent = sympy.Eq(S_N_n_sym, n_p_sym * S_N_1_sym)
    eq_S_N_n_noncoherent_Ei = sympy.Eq(S_N_n_sym, E_i_sym * n_p_sym * S_N_1_sym)
    eq_E_i_sqrt_n = sympy.Eq(E_i_sym, 1 / sqrt(n_p_sym))
    eq_S_N_n_noncoherent = eq_S_N_n_noncoherent_Ei.subs(E_i_sym, eq_E_i_sqrt_n.rhs)
    eq_theta_B_skolnik = sympy.Eq(theta_B_sym, 1.2 * wavelength_sym / D_h_sym)

    # --- Topic 10: Direction Finding Equations ---
    S_N_from_dB = sympy.Eq(S_N_sym, 10**(S_N_dB_sym / 10))
    v_phi = sympy.Eq(v_phi_sym, exp(-Theta_sym * (phi_sym - phi_s_sym)**2))
    Theta = sympy.Eq(Theta_sym, (4 * log(2)) / theta_B_sym**2)
    phi_hat_amp = sympy.Eq(phi_hat_sym, (Delta_sym / Sigma_sym) * (theta_B_sym**2 / (8 * log(2) * phi_s_sym)))
    sigma_phi_amp = sympy.Eq(sigma_phi_sym, (theta_B_sym**2 * sqrt(1 / S_N_sym)) / (8 * sqrt(2) * phi_s_sym * log(2)))
    sigma_phi_phase = sympy.Eq(sigma_phi_sym, (wavelength_sym / (2 * pi_sym * d_sym)) * sqrt(1 / S_N_sym))
    sigma_phi_time = sympy.Eq(sigma_phi_sym, c_sym / (d_sym * B_sym))

    # --- Topic 11: Pulse Compression Equations ---
    eq_delta_r_uncompressed = sympy.Eq(delta_r_sym, c_sym * tau_sym / 2)
    eq_B_chirp = sympy.Eq(B_sym, gamma_sym * tau_sym)
    eq_delta_r_compressed = sympy.Eq(delta_r_sym, c_sym / (2 * B_sym))
    eq_PCR_1 = sympy.Eq(PCR_sym, tau_sym * B_sym)
    eq_PCR_2 = sympy.Eq(PCR_sym, (tau_sym**2) * gamma_sym)
    eq_f_range_tone = sympy.Eq(f_range_tone_sym, -gamma_sym * (2 * R_offset_sym / c_sym))

class solve:
    def __init__():
        pass
    
    # =========================================================================
    # HELPER FUNCTION
    # =========================================================================
    
    def _solver(equation, solve_for_sym):
        """Internal helper to solve a SymPy equation and return a callable function."""
        try:
            # Solve the equation for the desired symbol
            sym_expr = sympy.solve(equation, solve_for_sym)[0]
        except (IndexError, Exception) as e:
            error_msg = f"Error: Could not solve equation for {solve_for_sym}: {e}"
            print(error_msg)
            def error_function():
                raise ValueError(error_msg)
            return error_function
            
        # Helper to convert Python numbers to sympy Floats but leave sympy types alone
        def _s(val):
            return val if isinstance(val, sympy.Basic) else sympy.Float(val)

        # Get all symbols used in the expression
        free_symbols = sym_expr.free_symbols
        
        def calculate():
            """Dynamically created solver function."""
            # Access the global v instance directly
            global v
            
            subs_map = {}
            for s in free_symbols:
                if hasattr(v, s.name):
                    subs_map[s] = _s(getattr(v, s.name))
                # Handle special cases like 'lambda'
                elif s.name == 'lambda' and hasattr(v, 'wavelength'):
                    subs_map[s] = _s(getattr(v, 'wavelength'))
                # Handle 'Delta f'
                elif str(s) == 'Delta f' and hasattr(v, 'deltaf'):
                    subs_map[s] = _s(getattr(v, 'deltaf'))
                else:
                    # This will catch sympy constants like pi, log(2), sqrt(2)
                    # which don't need substitution.
                    pass 
            
            value_sym = sym_expr.subs(subs_map)
            value_simpl = sympy.simplify(value_sym)
            return float(value_simpl.evalf())
        
        return calculate

    # =========================================================================
    # BASE/COMMON SOLVERS
    # =========================================================================
    
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
    
    # SymPy-based solvers
    P_t = _solver(equations.P_t, equations.P_t.lhs)
    R_max = _solver(equations.R_max, equations.R_max.lhs)

    # =========================================================================
    # TOPIC 07: DOPPLER CW RADAR SOLVERS
    # =========================================================================
    
    f_doppler = _solver(equations.eq_f_doppler, equations.eq_f_doppler.lhs)
    v_from_doppler = _solver(equations.eq_v_from_doppler, equations.eq_v_from_doppler.lhs)
    f_obs_if = _solver(equations.eq_f_obs_if, equations.eq_f_obs_if.lhs)
    delta_v = _solver(equations.eq_delta_v, equations.eq_delta_v.lhs)

    # =========================================================================
    # TOPIC 08: CWFM RADAR SOLVERS
    # =========================================================================
    
    R_cwfm = _solver(equations.R_cwfm, equations.R_cwfm.lhs)
    v_cwfm = _solver(equations.v_cwfm, equations.v_cwfm.lhs)
    f_m_cwfm = _solver(equations.f_m_cwfm, equations.f_m_cwfm.lhs)
    f_r_cwfm = _solver(equations.f_r_cwfm, equations.f_r_cwfm.lhs)
    f_d_cwfm = _solver(equations.f_d_cwfm, equations.f_d_cwfm.lhs)
    f_0_cwfm = _solver(equations.f_0_cwfm, equations.f_0_cwfm.lhs)

    # =========================================================================
    # TOPIC 09: PULSED RADAR & RANGE AMBIGUITY SOLVERS
    # =========================================================================
    
    R_un_from_fp = _solver(equations.eq_R_un_from_fp, equations.eq_R_un_from_fp.lhs)
    fp_from_R_un = _solver(equations.eq_fp_from_R_un, equations.eq_fp_from_R_un.lhs)
    R_from_time = _solver(equations.eq_R_from_time, equations.eq_R_from_time.lhs)
    tau_from_duty = _solver(equations.eq_tau_from_duty, equations.eq_tau_from_duty.lhs)
    
    def S_N_n_coherent_dB():
        """Calculates integrated SNR (coherent) in dB."""
        S_N_1_lin = convert.db_to_lin(vars.S_N_1_dB)
        S_N_n_lin = vars.n_p * S_N_1_lin
        return convert.lin_to_db(S_N_n_lin)

    def S_N_n_noncoherent_dB():
        """Calculates integrated SNR (non-coherent, E_i=1/sqrt(n)) in dB."""
        S_N_1_lin = convert.db_to_lin(vars.S_N_1_dB)
        S_N_n_lin = np.sqrt(vars.n_p) * S_N_1_lin
        return convert.lin_to_db(S_N_n_lin)

    # =========================================================================
    # TOPIC 10: DIRECTION FINDING SOLVERS
    # =========================================================================
    
    # SymPy-based solvers
    S_N_from_dB = _solver(equations.S_N_from_dB, equations.S_N_from_dB.lhs)
    sigma_phi_amplitude = _solver(equations.sigma_phi_amp, equations.sigma_phi_amp.lhs)
    phi_s_from_sigma_amp = _solver(equations.sigma_phi_amp, equations.phi_s_sym)
    sigma_phi_phase = _solver(equations.sigma_phi_phase, equations.sigma_phi_phase.lhs)
    d_from_sigma_phase = _solver(equations.sigma_phi_phase, equations.d_sym)
    sigma_phi_time = _solver(equations.sigma_phi_time, equations.sigma_phi_time.lhs)
    B_from_sigma_time = _solver(equations.sigma_phi_time, equations.B_sym)

    # Legacy math-based methods (kept for backward compatibility)
    def calculate_Theta():
        """Calculates the Theta parameter from the 3-dB beamwidth.
        
        Uses vars.theta_B
        
        Returns:
            Theta parameter
        """
        return (4 * math.log(2)) / (vars.theta_B**2)
    
    def v_phi():
        """Calculates the Gaussian beam approximation.
        
        Uses vars.phi, vars.phi_s, vars.Theta
        
        Returns:
            Gaussian beam approximation value
        """
        return math.exp(-vars.Theta * (vars.phi - vars.phi_s)**2)
    
    def v_phi_full():
        """Calculates v(phi) using theta_B directly.
        
        Uses vars.phi, vars.phi_s, vars.theta_B
        
        Returns:
            Beam approximation value
        """
        return math.exp((4 * math.log(2) * (vars.phi - vars.phi_s)**2) / (vars.theta_B**2))
    
    def estimate_phi_hat():
        """Calculates the linear processor angle estimate.
        
        Uses vars.Delta, vars.Sigma, vars.theta_B, vars.phi_s
        
        Returns:
            Angle estimate in degrees
        """
        return (vars.Delta / vars.Sigma) * (vars.theta_B**2 / (8 * math.log(2) * vars.phi_s))
    
    def sigma_phi_amplitude():
        """Calculates the angle standard deviation for amplitude comparison.
        
        Uses vars.theta_B, vars.S_N, vars.phi_s
        
        Returns:
            Angle standard deviation in degrees
        """
        return (vars.theta_B**2 * math.sqrt(1 / vars.S_N)) / (8 * math.sqrt(2) * vars.phi_s * math.log(2))
    
    def sigma_phi_phase():
        """Calculates the angle standard deviation for phase comparison.
        
        Uses vars.wavelength, vars.d, vars.S_N
        
        Returns:
            Angle standard deviation in radians
        """
        return (vars.wavelength / (2 * math.pi * vars.d)) * math.sqrt(1 / vars.S_N)
    
    def sigma_phi_time():
        """Calculates the angle standard deviation for time comparison.
        
        Uses vars.c, vars.d, vars.B
        
        Returns:
            Angle standard deviation in radians
        """
        return vars.c / (vars.d * vars.B)
    
    def db_to_linear():
        """Converts SNR from dB to linear.
        
        Uses vars.x (as the dB value)
        
        Returns:
            Linear value
        """
        return 10**(vars.x / 10)

    # =========================================================================
    # TOPIC 11: PULSE COMPRESSION SOLVERS
    # =========================================================================
    
    delta_r_uncompressed = _solver(equations.eq_delta_r_uncompressed, equations.eq_delta_r_uncompressed.lhs)
    B_chirp = _solver(equations.eq_B_chirp, equations.eq_B_chirp.lhs)
    delta_r_compressed = _solver(equations.eq_delta_r_compressed, equations.eq_delta_r_compressed.lhs)
    PCR_from_B = _solver(equations.eq_PCR_1, equations.eq_PCR_1.lhs)
    PCR_from_gamma = _solver(equations.eq_PCR_2, equations.eq_PCR_2.lhs)
    f_range_tone = _solver(equations.eq_f_range_tone, equations.eq_f_range_tone.lhs)
    R_offset_from_tone = _solver(equations.eq_f_range_tone, equations.R_offset_sym)

class convert:  # add alias con for convenience
    pass
    
    def rad_to_deg(value_radians):
        return value_radians * (180 / np.pi)
    
    def deg_to_rad(value_degrees):
        return value_degrees * (np.pi / 180)
    
    def lin_to_db(value_linear):
        # Handle non-positive inputs to avoid math errors
        if value_linear <= 0:
            return -np.inf
        return np.log(value_linear)/np.log(10)*10
    
    def db_to_lin(value_db):
        return 10**(value_db / 10.0)
    
    def m_to_mi(value_meters):
        return value_meters / 1609.34
    
    def mi_to_m(value_miles):
        return value_miles * 1609.34
        
    def nmi_to_m(value_nmi):
        return value_nmi * 1852.0

    def m_to_nmi(value_m):
        return value_m / 1852.0

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
            'km': 1000.0,
            'nmi': 1852.0,
            'mi': 1609.34
        }
        return value * conversion_factors.get(source_unit.lower(), 1)
    
con = convert()  # alias for convenience

def redefine_variable(var_name, new_value):
    """
    Redefines a global variable within the 'vars' namespace.
    Args:
        var_name (str): The name of the variable to redefine (e.g., "lambda").
        new_value: The new value to assign to the variable.
    """
    global v
    setattr(vars, var_name, new_value)
    # Also set on the global v instance for consistency
    setattr(v, var_name, new_value)
# Demonstration of usage from a separate script or module:

if __name__ == '__main__':  # Only runs when the script is executed directly
    
    # --- Original Demo Code ---
    print("--- Original Demo (CWFM) ---")
    # Use the global v instance directly
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
    v.f_r = solve.f_r_cwfm() # Use solver
    pprint(convert.hz_to(v.f_r, 'mhz'))
    v.f_d = Symbol('f_d')
    pprint(v.f_d)
    v.f_d = solve.f_d_cwfm() # Use solver
    pprint(convert.hz_to(v.f_d, 'mhz'))
    # plain word with space
    v.deltaf = Symbol('delta f')
    pprint(v.deltaf)

    # --- Topic 10: Direction Finding Demo ---
    print("\n" + "="*30)
    print("Topic 10: Direction Finding Demo")
    print("="*30 + "\n")

    # --- Problem 2: Amplitude Comparison ---
    print("--- Problem 2: Amplitude Comparison ---")
    redefine_variable('S_N_dB', 10)  # 10 dB
    redefine_variable('S_N', solve.S_N_from_dB()) # Convert to linear
    redefine_variable('phi_s', 5.0)   # 5 degrees
    redefine_variable('theta_B', 10.0) # 10 degrees
    print(f"Given: S/N = {v.S_N_dB} dB ({v.S_N:.2f} linear), Squint Angle = {v.phi_s} deg, Beamwidth = {v.theta_B} deg")
    redefine_variable('sigma_phi', solve.sigma_phi_amplitude())
    print(f"Calculated DOA Accuracy (sigma_phi): {v.sigma_phi:.2f} degrees")
    redefine_variable('sigma_phi', 0.5) # New target accuracy
    print(f"\nTarget DOA Accuracy: {v.sigma_phi} degrees")
    redefine_variable('phi_s', solve.phi_s_from_sigma_amp())
    print(f"Calculated required Squint Angle: {v.phi_s:.2f} degrees")
    
    # --- Problem 3a: Phase Comparison ---
    print("\n--- Problem 3a: Phase Comparison ---")
    redefine_variable('d', 2.0)     # 2 meters
    redefine_variable('f', convert.hz_from(10, 'ghz')) # 10 GHz
    redefine_variable('wavelength', solve.wavelength())
    redefine_variable('S_N_dB', 8)    # 8 dB
    redefine_variable('S_N', solve.S_N_from_dB())
    print(f"Given: d = {v.d} m, f = {convert.hz_to(v.f, 'ghz')} GHz (lambda = {v.wavelength:.3f} m), S/N = {v.S_N_dB} dB ({v.S_N:.2f} linear)")
    sigma_phi_rad = solve.sigma_phi_phase()
    sigma_phi_deg = convert.rad_to_deg(sigma_phi_rad)
    print(f"Calculated DOA Accuracy (sigma_phi): {sigma_phi_rad:.2e} rad ({sigma_phi_deg:.3f} deg)")
    target_sigma_phi_deg = 0.5
    target_sigma_phi_rad = convert.deg_to_rad(target_sigma_phi_deg)
    redefine_variable('sigma_phi', target_sigma_phi_rad)
    print(f"\nTarget DOA Accuracy: {target_sigma_phi_deg} degrees ({v.sigma_phi:.4f} rad)")
    redefine_variable('d', solve.d_from_sigma_phase())
    print(f"Calculated required Separation 'd': {v.d:.2f} m")

    # --- Problem 3b: Time Comparison ---
    print("\n--- Problem 3b: Time Comparison ---")
    redefine_variable('d', 2.0)     # 2 meters
    redefine_variable('B', convert.hz_from(200, 'mhz')) # 200 MHz
    print(f"Given: d = {v.d} m, B = {convert.hz_to(v.B, 'mhz')} MHz")
    sigma_phi_rad_time = solve.sigma_phi_time()
    sigma_phi_deg_time = convert.rad_to_deg(sigma_phi_rad_time)
    print(f"Calculated DOA Accuracy (sigma_phi): {sigma_phi_rad_time:.2f} rad ({sigma_phi_deg_time:.2f} deg)")
    target_sigma_phi_deg_2 = 0.5
    target_sigma_phi_rad_2 = convert.deg_to_rad(target_sigma_phi_deg_2)
    redefine_variable('sigma_phi', target_sigma_phi_rad_2)
    print(f"\nTarget DOA Accuracy: {target_sigma_phi_deg_2} degrees ({v.sigma_phi:.4f} rad)")
    redefine_variable('B', solve.B_from_sigma_time())
    print(f"Calculated required Bandwidth 'B': {convert.hz_to(v.B, 'ghz'):.2f} GHz")

    # --- Topic 07: Doppler CW Radar Demo ---
    print("\n" + "="*30)
    print("Topic 07: Doppler CW Radar Demo (Prob 2)")
    print("="*30 + "\n")
    redefine_variable('f', convert.hz_from(10, 'ghz')) # 10 GHz
    redefine_variable('wavelength', solve.wavelength())
    redefine_variable('f_if', convert.hz_from(100, 'mhz')) # 100 MHz IF
    redefine_variable('v', -100.0) # 100 m/s closing
    redefine_variable('T_cpi', 10e-3) # 10 ms
    
    redefine_variable('f_doppler', solve.f_doppler())
    print(f"Given: f = {convert.hz_to(v.f, 'ghz')} GHz (lambda = {v.wavelength:.3f} m), v = {v.v} m/s (closing)")
    print(f"Calculated Doppler Shift (f_doppler): {convert.hz_to(v.f_doppler, 'khz'):.2f} kHz")
    redefine_variable('f_obs', solve.f_obs_if())
    print(f"Given IF = {convert.hz_to(v.f_if, 'mhz')} MHz, Observed Freq = {convert.hz_to(v.f_obs, 'mhz'):.4f} MHz")
    
    redefine_variable('f_doppler', 10e3) # 10 kHz shift (Prob 2c)
    redefine_variable('v', solve.v_from_doppler())
    print(f"\nGiven Doppler Shift = {convert.hz_to(v.f_doppler, 'khz')} kHz, Calculated Velocity = {v.v:.2f} m/s (separating)")
    
    redefine_variable('delta_v', solve.delta_v())
    print(f"\nGiven CPI Time = {v.T_cpi * 1000} ms, Velocity Resolution = {v.delta_v:.2f} m/s")

    # --- Topic 08: CWFM Radar Demo ---
    print("\n" + "="*30)
    print("Topic 08: CWFM Radar Demo (Prob 2)")
    print("="*30 + "\n")
    redefine_variable('f', convert.hz_from(35, 'ghz')) # 35 GHz
    redefine_variable('wavelength', solve.wavelength())
    redefine_variable('f_m', 100.0) # 100 Hz
    redefine_variable('deltaf', convert.hz_from(30, 'mhz')) # 30 MHz
    redefine_variable('f_bu', convert.hz_from(85, 'khz')) # 85 kHz
    redefine_variable('f_bd', convert.hz_from(75, 'khz')) # 75 kHz
    
    redefine_variable('f_r', solve.f_r_cwfm())
    redefine_variable('f_d', solve.f_d_cwfm())
    print(f"Given: f_bu = {convert.hz_to(v.f_bu, 'khz')} kHz, f_bd = {convert.hz_to(v.f_bd, 'khz')} kHz")
    print(f"Calculated: f_r = {convert.hz_to(v.f_r, 'khz')} kHz, f_d = {convert.hz_to(v.f_d, 'khz')} kHz")
    
    redefine_variable('R', solve.R_cwfm())
    redefine_variable('v', solve.v_cwfm()) # Uses f_d, f
    print(f"Calculated Range = {v.R:.0f} m, Calculated Velocity = {v.v:.2f} m/s") # Matches 2000m, 21.5 m/s

    # --- Topic 09/10: Pulsed Radar & Range Ambiguity Demo ---
    print("\n" + "="*30)
    print("Topic 09/10: Pulsed Radar & Ambiguity Demo")
    print("="*30 + "\n")
    
    print("--- Problem 3 (Topic 10) ---")
    redefine_variable('R_un', convert.nmi_to_m(60)) # 60 nmi
    redefine_variable('f_p', solve.fp_from_R_un())
    print(f"Given R_un = 60 nmi, Calculated PRF (f_p) = {convert.hz_to(v.f_p, 'khz'):.2f} kHz")
    
    redefine_variable('R_un', convert.m_from(10, 'km')) # 10 km
    redefine_variable('T_p', 2 * v.R_un / v.c)
    redefine_variable('duty_cycle', 0.10) # 10%
    redefine_variable('tau', solve.tau_from_duty())
    print(f"\nGiven R_un = 10 km (T_p = {v.T_p * 1e6:.2f} us), Duty = 10%")
    print(f"Calculated Pulse Width (tau) = {v.tau * 1e6:.2f} us")
    
    print("\n--- Problem 7 (Topic 09) ---")
    redefine_variable('S_N_1_dB', 15.0) # 15 dB
    redefine_variable('n_p', 25) # 25 pulses
    print(f"Given: (S/N)_1 = {v.S_N_1_dB} dB, n_p = {v.n_p}")
    S_N_n_coh = solve.S_N_n_coherent_dB()
    S_N_n_noncoh = solve.S_N_n_noncoherent_dB()
    print(f"Calculated (S/N)_n (Coherent) = {S_N_n_coh:.2f} dB") # 15 + 10*log10(25) = 28.98 dB
    print(f"Calculated (S/N)_n (Non-Coherent, 1/sqrt(n)) = {S_N_n_noncoh:.2f} dB") # 15 + 5*log10(25) = 21.99 dB

    print("\n--- Problem 4 (Topic 10) - Range Ambiguity ---")
    prfs = [1000, 1250, 1500] # Hz
    detections = {
        1000: [20, 100], # km
        1250: [10, 80],  # km
        1500: [20, 50]   # km
    }
    max_range_km = 300 # Search up to 300 km
    unambiguous_ranges_km = {f: convert.m_to(v.c / (2 * f), 'km') for f in prfs}
    print(f"PRFs (Hz): {prfs}")
    print(f"Unambiguous Ranges (km): {[round(r,1) for r in unambiguous_ranges_km.values()]}")
    
    possible_ranges = {}
    for prf, ru in unambiguous_ranges_km.items():
        possible_ranges[prf] = set()
        for R_measured in detections[prf]:
            for n in range(int(max_range_km / ru) + 1):
                R_true = R_measured + n * ru
                if R_true <= max_range_km:
                    possible_ranges[prf].add(round(R_true))
                    
    print(f"Possible Ranges (km) for PRF 1000: {sorted(list(possible_ranges[1000]))}")
    print(f"Possible Ranges (km) for PRF 1250: {sorted(list(possible_ranges[1250]))}")
    print(f"Possible Ranges (km) for PRF 1500: {sorted(list(possible_ranges[1500]))}")
    
    # Find intersection
    true_targets = possible_ranges[1000].intersection(possible_ranges[1250]).intersection(possible_ranges[1500])
    print(f"True Target Ranges (km): {true_targets}") # Matches {170, 200} km (rounding)
    
    # --- Topic 11: Pulse Compression Demo ---
    print("\n" + "="*30)
    print("Topic 11: Pulse Compression Demo (Prob 1 & 4)")
    print("="*30 + "\n")
    
    redefine_variable('tau', 20e-6) # 20 us
    redefine_variable('delta_r', solve.delta_r_uncompressed())
    print(f"Given: Uncompressed Pulse Width (tau) = {v.tau * 1e6:.0f} us")
    print(f"Calculated Uncompressed Range Res = {v.delta_r:.0f} m") # 3000m
    
    redefine_variable('gamma', 10e6 / 1e-6) # 10 MHz/us = 10e12 Hz/s
    redefine_variable('B', solve.B_chirp())
    redefine_variable('delta_r', solve.delta_r_compressed())
    print(f"\nGiven: Chirp Rate (gamma) = 10 MHz/us, tau = {v.tau * 1e6:.0f} us")
    print(f"Calculated Bandwidth (B) = {convert.hz_to(v.B, 'mhz')} MHz")
    print(f"Calculated Compressed Range Res = {v.delta_r:.2f} m") # 0.75m
    
    redefine_variable('PCR', solve.PCR_from_B())
    print(f"Calculated PCR (from B*tau) = {v.PCR:.0f}") # 4000
    
    redefine_variable('R_offset', 100) # 100m
    redefine_variable('f_range_tone', solve.f_range_tone())
    print(f"\nGiven: Range Offset = {v.R_offset} m")
    print(f"Calculated Dechirp Range Tone = {convert.hz_to(v.f_range_tone, 'mhz'):.2f} MHz") # -6.67 MHz