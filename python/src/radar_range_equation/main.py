"""Radar Range Equation Package - Core Module.

This module provides a comprehensive toolkit for radar range equation calculations,
including symbolic and numeric solutions for various radar types (CW, CWFM, pulsed),
direction finding, and pulse compression techniques.

The module is organized into five main classes:
    - vars: Container for physical constants and radar parameters
    - equations: Symbolic SymPy equations for radar calculations
    - solve: Numeric solver functions for radar problems
    - convert: Unit conversion utilities
    - analysis: Analysis helpers for pulse parsing, integration, jammers

Example:
    >>> import radar_range_equation as RRE
    >>> RRE.vars.f = 10e9  # 10 GHz
    >>> RRE.vars.wavelength = RRE.solve.wavelength()
    >>> print(RRE.vars.wavelength)
"""

import sympy
import scipy
import numpy as np
import math
import sympy.physics.units as units
from sympy import symbols, Symbol, pprint, exp, sqrt, log
from sympy.physics.units import convert_to
from sympy.physics.units import speed_of_light

class vars:
    """Container for radar system variables and physical constants.
    
    This class provides a centralized namespace for all radar-related variables,
    including physical constants, antenna parameters, radar equation parameters,
    and topic-specific variables for different radar types.
    
    Attributes:
        c (float): Speed of light in m/s (from scipy.constants.c)
        k (float): Boltzmann constant in J/K (from scipy.constants.Boltzmann)
        pi (float): Mathematical constant pi
        pi4 (float): 4*pi constant
        f (Symbol): Frequency (symbolic)
        wavelength (Symbol): Wavelength (symbolic, use getattr/setattr for 'lambda')
        
        Antenna Parameters:
            A_e (Symbol): Effective aperture
            A (Symbol): Antenna area
            D_h (Symbol): Horizontal antenna dimension
            D_v (Symbol): Vertical antenna dimension
            D (Symbol): Antenna diameter
            eta (Symbol): Antenna efficiency
            G_t (Symbol): Transmit antenna gain
            G_r (Symbol): Receive antenna gain
            theta_B (Symbol): Beamwidth
            
        Radar Equation Parameters:
            R (Symbol): Range
            R_max (Symbol): Maximum range
            P_t (Symbol): Transmit power
            S_min (Symbol): Minimum detectable signal
            sigma (Symbol): Radar cross section
            
        Doppler CW Radar (Topic 07):
            f_doppler (Symbol): Doppler frequency shift
            f_if (Symbol): Intermediate frequency
            f_obs (Symbol): Observed frequency at receiver
            T_cpi (Symbol): Coherent processing time
            delta_v (Symbol): Velocity resolution
            
        CWFM Radar (Topic 08):
            R_un (Symbol): Unnormalized range
            f_m (Symbol): Modulation frequency
            f_bu (Symbol): Upper band frequency
            f_bd (Symbol): Lower band frequency
            f_r (Symbol): Radar operating frequency
            f_d (Symbol): Frequency deviation
            
        Pulsed Radar (Topic 09):
            f_p (Symbol): Pulse Repetition Frequency (PRF)
            T_p (Symbol): Pulse Repetition Interval (PRI)
            tau (Symbol): Pulse width
            n_p (Symbol): Number of pulses integrated
            S_N_1 (Symbol): Single pulse SNR
            
        Direction Finding (Topic 10):
            phi (Symbol): Angle
            phi_s (Symbol): Squint angle
            d (Symbol): Element separation
            S_N (Symbol): Signal-to-Noise ratio
            B (Symbol): Bandwidth
            
        Pulse Compression (Topic 11):
            delta_r (Symbol): Range resolution
            gamma (Symbol): Chirp rate
            PCR (Symbol): Pulse Compression Ratio
    """
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
    g = 9.80665                         # Gravitational acceleration (m/s^2)
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
    # TOPIC 12: CHAFF VARS
    # =========================================================================
    L_fiber = Symbol('L_fiber')         # Chaff Fiber Length (m)
    D_fiber = Symbol('D_fiber')         # Chaff Fiber Diameter (m)
    V_ch = Symbol('V_ch')               # Volume of a single chaff fiber (m^3)
    V_box = Symbol('V_box')             # Volume of chaff cartridge (m^3)
    Fill_ratio = Symbol('Fill_ratio')   # Chaff cartridge fill ratio (dimensionless)
    N_fiber = Symbol('N_fiber')         # Number of chaff fibers (dimensionless)
    sigma_ch = Symbol('sigma_ch')       # Average RCS of chaff (m^2)
    zeta_ch = Symbol('zeta_ch')         # Chaff dispersion constant (s)

    # =========================================================================
    # TOPIC 13: NOISE JAMMING VARS
    # =========================================================================
    Pj = Symbol('Pj')                   # Jammer Transmit Power (W)
    Gj = Symbol('Gj')                   # Jammer Antenna Gain (linear)
    Bj = Symbol('Bj')                   # Jammer Bandwidth (Hz)
    Lossj = Symbol('Lossj')             # Jammer Loss (linear)
    S_J_ratio = Symbol('S/J')           # Signal-to-Jammer Ratio (linear)
    R_bt = Symbol('R_bt')               # Burnthrough Range (m)

    # =========================================================================
    # TOPIC 14: GATED NOISE VARS
    # =========================================================================
    R_tgt = Symbol('R_tgt')             # Target Range (m)
    R_gn_start_offset = Symbol('R_gn_start_offset') # Range before target to start noise (m)
    t_tgt_2way = Symbol('t_tgt_2way')   # Two-way time to target (s)
    t_gn_start_release = Symbol('t_gn_start_release') # Gated Noise start time (s)
    Delta_R_mask = Symbol('Delta_R_mask') # Total masking range (m)

    # =========================================================================
    # TOPIC 15: FALSE TARGET GENERATION VARS
    # =========================================================================
    v_tgt = Symbol('v_tgt')             # Target range rate (m/s)
    R_ft = Symbol('R_ft')               # False Target Range (m)
    v_ft = Symbol('v_ft')               # False Target range rate (m/s)
    f_D_tgt = Symbol('f_D_tgt')         # Doppler Frequency Target (Hz)
    f_D_ft = Symbol('f_D_ft')           # Doppler Frequency False Target (Hz)
    Delta_t_ft = Symbol('Delta_t_ft')   # Time Delay for False Target (s)
    Delta_f_ft = Symbol('Delta_f_ft')   # Frequency Shift for False Target (Hz)

    # =========================================================================
    # TOPIC 16: RADAR TRACKING / FALSE TRACKS VARS
    # =========================================================================
    P_density = Symbol('P_density')     # Power Density (W/m^2)
    Pj_emulated = Symbol('Pj_emulated') # Jammer Power to emulate target RCS (W)

    # =========================================================================
    # TOPIC 17: GATE STEALING VARS
    # =========================================================================

    alpha = Symbol('alpha')             # Target acceleration (effective for range) (symbolic)
    T_time = Symbol('T')                # Time duration (symbolic)
    Delta_r_max = Symbol('Delta_r_max') # Maximum required range offset (symbolic)
    Delta_r_t = Symbol('Delta_r(t)')    # Range offset profile (symbolic)
    rho_v = Symbol('rho_v')             # Velocity resolution (Doppler bin size) (symbolic)
    n_gate_r = Symbol('n_gate_r')       # Range gate size (resolution cells) (symbolic)
    n_gate_v = Symbol('n_gate_v')       # Velocity gate size (resolution cells) (symbolic)
    Delta_v_max = Symbol('Delta_v_max') # Maximum required velocity offset (symbolic)
    a_accel = Symbol('a')               # Target acceleration for velocity (symbolic)
    Delta_v_t = Symbol('Delta_v(t)')    # Velocity offset profile (symbolic)

    # =========================================================================
    # TOPIC 18: CROSS-EYE VARS
    # =========================================================================

    L_cross = Symbol('L')               # Cross-eye aperture separation (symbolic)
    a_gain_ratio = Symbol('a_ratio')    # Jammer gain ratio (J1/J2) (symbolic)
    J_1 = Symbol('J_1')                 # Jammer channel 1 gain (symbolic)
    J_2 = Symbol('J_2')                 # Jammer channel 2 gain (symbolic)
    phi_hat_ce = Symbol('phi_hat_ce')   # Cross-eye angle error (symbolic)
    S_phi_bar = Symbol('S_phi_bar')     # Normalized Monopulse Slope (symbolic)

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
    """Symbolic SymPy equations for radar calculations.
    
    This class contains symbolic representations of radar equations using SymPy.
    Each equation is stored as a SymPy Eq object that can be manipulated symbolically
    or solved numerically using the solve class.
    
    The equations are organized by topic:
        - Base/Common: Fundamental radar equations (wavelength, gain, range)
        - Topic 07: Doppler CW radar equations
        - Topic 08: CWFM radar equations
        - Topic 09: Pulsed radar equations
        - Topic 10: Direction finding equations
        - Topic 11: Pulse compression equations
    
    Attributes:
        A_e (Eq): Effective aperture equation
        wavelength (Eq): Wavelength equation (lambda = c/f)
        G_t (Eq): Transmit antenna gain equation
        R_max (Eq): Maximum radar range equation
        Linear_to_dB (Eq): Linear to dB conversion equation
        
        Doppler CW Radar:
            eq_f_doppler (Eq): Doppler frequency shift equation
            eq_v_from_doppler (Eq): Velocity from Doppler shift
            eq_delta_v (Eq): Velocity resolution equation
            
        CWFM Radar:
            R_cwfm (Eq): Range equation for CWFM radar
            v_cwfm (Eq): Velocity equation for CWFM radar
            
        Pulsed Radar:
            eq_R_un_from_fp (Eq): Unambiguous range from PRF
            eq_T_p (Eq): Pulse repetition interval equation
            eq_S_N_n_coherent (Eq): Coherent integration SNR
            
        Direction Finding:
            phi_hat_amp (Eq): Angle estimate for amplitude comparison
            sigma_phi_amp (Eq): Angle accuracy for amplitude comparison
            sigma_phi_phase (Eq): Angle accuracy for phase comparison
            sigma_phi_time (Eq): Angle accuracy for time comparison
            
        Pulse Compression:
            eq_delta_r_uncompressed (Eq): Uncompressed range resolution
            eq_delta_r_compressed (Eq): Compressed range resolution
            eq_PCR_1 (Eq): Pulse compression ratio equation
    
    Example:
        >>> from radar_range_equation import equations
        >>> print(equations.wavelength)
        Eq(lambda, c/f)
    """
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

    # --- Topic 12: Chaff Symbols ---
    L_fiber_sym = Symbol('L_fiber')
    D_fiber_sym = Symbol('D_fiber')
    V_ch_sym = Symbol('V_ch')
    V_box_sym = Symbol('V_box')
    Fill_ratio_sym = Symbol('Fill_ratio')
    N_fiber_sym = Symbol('N_fiber')
    sigma_ch_sym = Symbol('sigma_ch')
    zeta_ch_sym = Symbol('zeta_ch')
    
    # --- Topic 13: Noise Jamming Symbols ---
    Pj_sym = Symbol('Pj')
    Gj_sym = Symbol('Gj')
    Bj_sym = Symbol('Bj')
    Lossj_sym = Symbol('Lossj')
    S_J_ratio_sym = Symbol('S/J')
    R_bt_sym = Symbol('R_bt')
    
    # --- Topic 14: Gated Noise Symbols ---
    R_tgt_sym = Symbol('R_tgt')
    R_gn_start_offset_sym = Symbol('R_gn_start_offset')
    t_tgt_2way_sym = Symbol('t_tgt_2way')
    t_gn_start_release_sym = Symbol('t_gn_start_release')
    Delta_R_mask_sym = Symbol('Delta_R_mask')
    
    # --- Topic 15: False Target Generation Symbols ---
    v_tgt_sym = Symbol('v_tgt')
    R_ft_sym = Symbol('R_ft')
    v_ft_sym = Symbol('v_ft')
    f_D_tgt_sym = Symbol('f_D_tgt')
    f_D_ft_sym = Symbol('f_D_ft')
    Delta_t_ft_sym = Symbol('Delta_t_ft')
    Delta_f_ft_sym = Symbol('Delta_f_ft')
    
    # --- Topic 16: Radar Tracking / False Tracks Symbols ---
    P_density_sym = Symbol('P_density')
    Pj_emulated_sym = Symbol('Pj_emulated')

    # --- Topic 17: Gate Stealing Symbols ---
    alpha_sym = Symbol('alpha')
    g_sym = Symbol('g')
    T_time_sym = Symbol('T')
    Delta_r_max_sym = Symbol('Delta_r_max')
    Delta_r_t_sym = Symbol('Delta_r(t)')
    rho_v_sym = Symbol('rho_v')
    n_gate_r_sym = Symbol('n_gate_r')
    n_gate_v_sym = Symbol('n_gate_v')
    Delta_v_max_sym = Symbol('Delta_v_max')
    a_accel_sym = Symbol('a')
    Delta_v_t_sym = Symbol('Delta_v(t)')
    t_delay_sym = Symbol('t')

    # --- Topic 18: Cross-Eye Symbols ---
    L_cross_sym = Symbol('L')
    a_gain_ratio_sym = Symbol('a_ratio')
    J_1_sym = Symbol('J_1')
    J_2_sym = Symbol('J_2')
    phi_hat_ce_sym = Symbol('phi_hat_ce')
    S_phi_bar_sym = Symbol('S_phi_bar')

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

    # --- Topic 12: Chaff Equations ---
    # Chaff Fiber Length
    eq_L_fiber = sympy.Eq(L_fiber_sym, wavelength_sym / 2)
    # Volume of a Single Chaff Fiber (approximated as cylinder)
    eq_V_ch = sympy.Eq(V_ch_sym, (sympy.pi * L_fiber_sym * D_fiber_sym**2) / 4)
    # Number of Fibers in Box
    eq_N_fiber = sympy.Eq(N_fiber_sym, (V_box_sym * Fill_ratio_sym) / V_ch_sym)
    # RCS of Chaff Cloud (Time-dependent)
    eq_sigma_ch_t = sympy.Eq(sigma_ch_sym, 0.15 * N_fiber_sym * wavelength_sym**2 * (1 - exp(-t_delay_sym / zeta_ch_sym)))
    # Max RCS of Chaff Cloud (t -> inf)
    eq_sigma_ch_max = sympy.Eq(sigma_ch_sym, 0.15 * N_fiber_sym * wavelength_sym**2)
    
    # --- Topic 13: Noise Jamming Equations ---
    # Signal-to-Jammer (S/J) Ratio for Barrage Noise
    eq_S_J_ratio = sympy.Eq(S_J_ratio_sym, (P_t_sym * G_t_sym * sigma_sym * n_p_sym * Bj_sym) / (4 * sympy.pi * R_sym**2 * Pj_sym * Gj_sym * Lossj_sym * B_sym))
    # Burnthrough Range
    eq_R_bt = sympy.Eq(R_bt_sym**2, (P_t_sym * G_t_sym * sigma_sym * n_p_sym * Bj_sym) / (4 * sympy.pi * Pj_sym * Gj_sym * Lossj_sym * B_sym * S_min_sym))
    
    # --- Topic 14: Gated Noise Equations ---
    # Two-way time to target
    eq_t_tgt_2way = sympy.Eq(t_tgt_2way_sym, 2 * R_tgt_sym / c_sym)
    # Gated Noise Start Release Time
    eq_t_gn_start_release = sympy.Eq(t_gn_start_release_sym, (2 * (R_tgt_sym - R_gn_start_offset_sym) / c_sym) - (tau_sym / 2))
    
    # --- Topic 15: False Target Generation Equations ---
    # Target Doppler Frequency
    eq_f_D_tgt = sympy.Eq(f_D_tgt_sym, -2 * v_tgt_sym / wavelength_sym)
    # False Target Doppler Frequency
    eq_f_D_ft = sympy.Eq(f_D_ft_sym, -2 * v_ft_sym / wavelength_sym)
    # Time Delay to apply
    eq_Delta_t_ft = sympy.Eq(Delta_t_ft_sym, (2 * R_ft_sym / c_sym) - (2 * R_sym / c_sym))
    # Frequency Shift to apply
    eq_Delta_f_ft = sympy.Eq(Delta_f_ft_sym, f_D_ft_sym - f_D_tgt_sym)

    # --- Topic 16: Radar Tracking / False Tracks Equations ---
    # Power Density at Jammer/Target
    eq_P_density = sympy.Eq(P_density_sym, (P_t_sym * G_t_sym) / (4 * sympy.pi * R_sym**2))
    # Required Jammer Power to Emulate RCS
    eq_Pj_emulated = sympy.Eq(Pj_emulated_sym, (P_density_sym * sigma_sym) / Gj_sym)

    # --- Topic 17: Gate Stealing Equations ---
    # Velocity Resolution (Doppler bin size)
    eq_rho_v = sympy.Eq(rho_v_sym, wavelength_sym / (2 * T_cpi_sym))
    # Required Max Range Offset (Gate Size * Resolution)
    eq_Delta_r_max_gate = sympy.Eq(Delta_r_max_sym, n_gate_r_sym * delta_r_sym)
    # Time (T) from Max Range Offset (Delta_r_max = 0.5 * alpha * T^2)
    eq_T_from_Delta_r = sympy.Eq(Delta_r_max_sym, sympy.Rational(1, 2) * alpha_sym * T_time_sym**2)
    # Range Offset Profile (Delta_r(t) = 0.5 * alpha * t^2)
    eq_Delta_r_t = sympy.Eq(Delta_r_t_sym, sympy.Rational(1, 2) * alpha_sym * t_delay_sym**2)
    # Required Max Velocity Offset (Gate Size * Resolution)
    eq_Delta_v_max_gate = sympy.Eq(Delta_v_max_sym, n_gate_v_sym * rho_v_sym)
    # Time (T) from Max Velocity Offset (Delta_v_max = a * T)
    eq_T_from_Delta_v = sympy.Eq(Delta_v_max_sym, a_accel_sym * T_time_sym)
    # Velocity Offset Profile (Delta_v(t) = a * t)
    eq_Delta_v_t = sympy.Eq(Delta_v_t_sym, a_accel_sym * t_delay_sym)

    # --- Topic 18: Cross-Eye Equations ---
    # Jammer Gain Ratio (a)
    eq_J_ratio = sympy.Eq(a_gain_ratio_sym, J_1_sym / J_2_sym)
    # Apparent Cross-Eye Angle Error (Amplitude Monopulse Approximation)
    eq_phi_hat_ce_amp = sympy.Eq(phi_hat_ce_sym, sympy.Rational(1, 2) * (L_cross_sym / R_sym) * (1 + a_gain_ratio_sym) / (1 - a_gain_ratio_sym))

class solve:
    """Numeric solver functions for radar calculations.
    
    This class provides methods to solve radar equations numerically using values
    from the vars class. Functions either compute values directly using Python/NumPy
    or use the _solver helper to create callable functions from symbolic equations.
    
    The solver functions are organized by topic:
        - Base/Common: Fundamental calculations (wavelength, gain, range)
        - Topic 07: Doppler CW radar solvers
        - Topic 08: CWFM radar solvers
        - Topic 09: Pulsed radar solvers
        - Topic 10: Direction finding solvers
        - Topic 11: Pulse compression solvers
    
    Example:
        >>> import radar_range_equation as RRE
        >>> RRE.vars.c = 3e8
        >>> RRE.vars.f = 10e9
        >>> wavelength = RRE.solve.wavelength()
        >>> print(f"Wavelength: {wavelength} m")
    """
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
        def _s(v):
            return v if isinstance(v, sympy.Basic) else sympy.Float(v)

        # Get all symbols used in the expression
        free_symbols = sym_expr.free_symbols
        
        def calculate():
            """Dynamically created solver function."""
            subs_map = {}
            for s in free_symbols:
                if hasattr(vars, s.name):
                    subs_map[s] = _s(getattr(vars, s.name))
                # Handle special cases like 'lambda'
                elif s.name == 'lambda' and hasattr(vars, 'wavelength'):
                    subs_map[s] = _s(getattr(vars, 'wavelength'))
                # Handle 'Delta f'
                elif str(s) == 'Delta f' and hasattr(vars, 'deltaf'):
                    subs_map[s] = _s(getattr(vars, 'deltaf'))
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
        """Calculate the radius of a sphere from its radar cross section.
        
        Uses vars.sigma (radar cross section in m²).
        
        Returns:
            float: Radius of the sphere in meters.
        
        Example:
            >>> import radar_range_equation as RRE
            >>> RRE.vars.sigma = 3.14159
            >>> radius = RRE.solve.A_sphere()
        """
        return (vars.sigma / np.pi) ** (1/2)
    
    def sigma_sphere():
        """Calculate the radar cross section of a sphere from its area.
        
        Uses vars.A (antenna area in m²).
        
        Returns:
            float: Radar cross section in m².
        
        Example:
            >>> import radar_range_equation as RRE
            >>> RRE.vars.A = 1.0
            >>> rcs = RRE.solve.sigma_sphere()
        """
        return np.pi * vars.A ** 2

    def theta_B():
        """Calculate the 3-dB beamwidth using Gaussian approximation.
        
        Uses vars.wavelength (wavelength in m) and vars.D_h (horizontal antenna 
        dimension in m). Uses the approximation: theta_B = 65° * pi/180 * lambda/D_h.
        
        Returns:
            float: Beamwidth in radians.
        
        Example:
            >>> import radar_range_equation as RRE
            >>> RRE.vars.wavelength = 0.03
            >>> RRE.vars.D_h = 1.0
            >>> beamwidth = RRE.solve.theta_B()
        """
        return 65 * vars.pi / 180 * (vars.wavelength / vars.D_h)

    def A_e_rect():
        """Calculate effective aperture for a rectangular antenna.
        
        Uses vars.eta (antenna efficiency), vars.D_h (horizontal dimension in m),
        and vars.D_v (vertical dimension in m).
        
        Returns:
            float: Effective aperture in m².
        
        Example:
            >>> import radar_range_equation as RRE
            >>> RRE.vars.eta = 0.6
            >>> RRE.vars.D_h = 2.0
            >>> RRE.vars.D_v = 1.5
            >>> aperture = RRE.solve.A_e_rect()
        """
        return vars.eta * vars.D_h * vars.D_v

    def A_e_circ():
        """Calculate effective aperture for a circular antenna.
        
        Uses vars.eta (antenna efficiency) and vars.D (antenna diameter in m).
        
        Returns:
            float: Effective aperture in m².
        
        Example:
            >>> import radar_range_equation as RRE
            >>> RRE.vars.eta = 0.6
            >>> RRE.vars.D = 2.0
            >>> aperture = RRE.solve.A_e_circ()
        """
        return vars.eta * vars.pi * (vars.D / 2) ** 2

    def wavelength():
        """Calculate wavelength from frequency.
        
        Uses vars.c (speed of light in m/s) and vars.f (frequency in Hz).
        
        Returns:
            float: Wavelength in meters.
        
        Example:
            >>> import radar_range_equation as RRE
            >>> RRE.vars.c = 3e8
            >>> RRE.vars.f = 10e9
            >>> wl = RRE.solve.wavelength()
        """
        return vars.c / vars.f

    def G_t():
        """Calculate transmit antenna gain.
        
        Uses vars.pi (pi constant), vars.A_e (effective aperture in m²),
        and vars.wavelength (wavelength in m).
        
        Returns:
            float: Antenna gain (dimensionless).
        
        Example:
            >>> import radar_range_equation as RRE
            >>> RRE.vars.A_e = 1.0
            >>> RRE.vars.wavelength = 0.03
            >>> gain = RRE.solve.G_t()
        """
        return 4 * vars.pi * vars.A_e / (vars.wavelength ** 2)
    
    def R4():
        """Calculate R^4 from the radar range equation.
        
        Uses vars.P_t (transmit power), vars.G_t (transmit gain), vars.G_r (receive gain),
        vars.wavelength (wavelength in m), vars.sigma (radar cross section in m²),
        vars.pi4 (4*pi), and vars.S_min (minimum detectable signal).
        
        Returns:
            float: R^4 value. Take the fourth root to get range in meters.
        
        Example:
            >>> import radar_range_equation as RRE
            >>> RRE.vars.P_t = 1000
            >>> RRE.vars.G_t = 1000
            >>> RRE.vars.G_r = 1000
            >>> r4 = RRE.solve.R4()
            >>> r = r4 ** 0.25  # Get actual range
        """
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

    # =========================================================================
    # TOPIC 12: CHAFF SOLVERS
    # =========================================================================
    @staticmethod
    def L_fiber():
        """Calculate Chaff Fiber Length (lambda/2).
        Uses vars.wavelength (m).
        Returns: float: Fiber length (m)."""
        return vars.wavelength / 2

    @staticmethod
    def V_ch():
        """Calculate Volume of a single chaff fiber (cylinder).
        Uses vars.L_fiber (m) and vars.D_fiber (m).
        Returns: float: Fiber volume (m^3)."""
        # V_ch = pi * L * D^2 / 4
        return (vars.pi * vars.L_fiber * vars.D_fiber**2) / 4

    @staticmethod
    def N_fiber():
        """Calculate Number of fibers in a cartridge.
        Uses vars.V_box (m^3), vars.Fill_ratio (dim), and vars.V_ch (m^3).
        Returns: float: Number of fibers (dimensionless)."""
        # N_fiber = (V_box * Fill_ratio) / V_ch
        if vars.V_ch == 0:
            raise ValueError("Chaff fiber volume (V_ch) cannot be zero.")
        return (vars.V_box * vars.Fill_ratio) / vars.V_ch

    @staticmethod
    def sigma_ch_t(t_s):
        """Calculate RCS of Chaff Cloud at time t.
        Uses vars.N_fiber (dim), vars.wavelength (m), vars.zeta_ch (s), and t_s (s).
        Returns: float: RCS (m^2)."""
        # sigma_ch(t) = 0.15 * N * lambda^2 * (1 - e^(-t/zeta))
        return 0.15 * vars.N_fiber * vars.wavelength**2 * (1 - math.exp(-t_s / vars.zeta_ch))
        
    # =========================================================================
    # TOPIC 13: NOISE JAMMING SOLVERS
    # =========================================================================
    @staticmethod
    def S_J_ratio():
        """Calculate Signal-to-Jammer (S/J) Ratio for Barrage Noise.
        Uses Pt, Gt, sigma, N_p, Bj, R, Pj, Gj, Lossj, B.
        Returns: float: S/J Ratio (linear)."""
        # S/J = (Pt * Gt * sigma * N_p * Bj) / (4 * pi * R^2 * Pj * Gj * Lossj * B)
        numerator = vars.Pt * vars.Gt * vars.sigma * vars.N_p * vars.Bj
        denominator = 4 * vars.pi * vars.R**2 * vars.Pj * vars.Gj * vars.Lossj * vars.B
        if denominator == 0:
            raise ValueError("Denominator is zero. Check Pj, R, Gj, Lossj, or B.")
        return numerator / denominator

    @staticmethod
    def R_burnthrough():
        """Calculate Burnthrough Range (R_bt).
        Uses Pt, Gt, sigma, N_p, Bj, Pj, Gj, Lossj, B, S_min.
        Returns: float: Burnthrough Range (m)."""
        # R_bt = sqrt( (Pt * Gt * sigma * N_p * Bj) / (4 * pi * Pj * Gj * Lossj * B * S_min) )
        numerator = vars.Pt * vars.Gt * vars.sigma * vars.N_p * vars.Bj
        denominator = 4 * vars.pi * vars.Pj * vars.Gj * vars.Lossj * vars.B * vars.S_min
        if denominator <= 0:
            return np.inf
        return math.sqrt(numerator / denominator)

    # =========================================================================
    # TOPIC 14: GATED NOISE SOLVERS
    # =========================================================================
    @staticmethod
    def t_tgt_2way():
        """Calculate Two-way Time of Flight to Target.
        Uses vars.R_tgt (m) and vars.c (m/s).
        Returns: float: Time (s)."""
        return 2 * vars.R_tgt / vars.c

    @staticmethod
    def t_gn_start_release():
        """Calculate Gated Noise Start Release Time (relative to T=0 radar pulse).
        Uses vars.R_tgt, vars.R_gn_start_offset, vars.c, vars.tau.
        Returns: float: Time (s)."""
        # t_gn_start_release = 2 * (R_tgt - R_gn_start_offset) / c - tau / 2
        R_start = vars.R_tgt - vars.R_gn_start_offset
        return (2 * R_start / vars.c) - (vars.tau / 2)

    # =========================================================================
    # TOPIC 15: FALSE TARGET GENERATION SOLVERS
    # =========================================================================
    @staticmethod
    def f_D_tgt():
        """Calculate Target Doppler Frequency.
        Uses vars.v_tgt (m/s) and vars.wavelength (m).
        Returns: float: Doppler frequency (Hz)."""
        return -2 * vars.v_tgt / vars.wavelength

    @staticmethod
    def f_D_ft():
        """Calculate False Target Doppler Frequency.
        Uses vars.v_ft (m/s) and vars.wavelength (m).
        Returns: float: Doppler frequency (Hz)."""
        return -2 * vars.v_ft / vars.wavelength

    @staticmethod
    def Delta_t_ft():
        """Calculate Time Delay to apply for False Target Generation.
        Uses vars.R_ft (m), vars.R (m), and vars.c (m/s).
        Returns: float: Time delay (s)."""
        # Delta_t_ft = 2 * (R_ft - R) / c
        return 2 * (vars.R_ft - vars.R) / vars.c

    @staticmethod
    def Delta_f_ft():
        """Calculate Frequency Shift to apply for False Target Generation.
        Uses vars.f_D_ft (Hz) and vars.f_D_tgt (Hz).
        Returns: float: Frequency shift (Hz)."""
        # Delta_f_ft = f_D_ft - f_D_tgt
        return vars.f_D_ft - vars.f_D_tgt

    # =========================================================================
    # TOPIC 16: RADAR TRACKING / FALSE TRACKS SOLVERS
    # =========================================================================
    @staticmethod
    def P_density():
        """Calculate Power Density at Target.
        Uses vars.Pt (W), vars.Gt (dim), and vars.R (m).
        Returns: float: Power density (W/m^2)."""
        # P_density = (Pt * Gt) / (4 * pi * R^2)
        if vars.R == 0:
            raise ValueError("Range (R) cannot be zero.")
        return (vars.Pt * vars.Gt) / (4 * vars.pi * vars.R**2)

    @staticmethod
    def Pj_emulated():
        """Calculate Jammer Power required to emulate a target RCS (sigma).
        Uses vars.P_density (W/m^2), vars.sigma (m^2), and vars.Gj (dim).
        Returns: float: Jammer transmit power (W)."""
        # Pj_emulated = (P_density * sigma) / Gj
        if vars.Gj == 0:
            raise ValueError("Jammer gain (Gj) cannot be zero.")
        return (vars.P_density * vars.sigma) / vars.Gj

    # =========================================================================
    # TOPIC 17: GATE STEALING SOLVERS
    # =========================================================================

    @staticmethod
    def rho_v():
        """Calculate Velocity Resolution (Doppler bin size).
        Uses vars.wavelength (m) and vars.T_cpi (s).
        Returns: float: Velocity resolution (m/s)."""
        # rho_v = lambda / (2 * T_cpi)
        if not hasattr(vars, 'T_cpi') or not vars.T_cpi:
            raise ValueError("T_cpi (Coherent Processing Interval) is not set or zero.")
        return (vars.wavelength / (2 * vars.T_cpi))

    @staticmethod
    def Delta_r_max_gate():
        """Calculate the maximum required range offset to exit the gate.
        Uses vars.n_gate_r (cells) and vars.delta_r (m).
        Returns: float: Maximum range offset (m)."""
        # Delta_r_max = n_gate_r * delta_r (Requires delta_r to be set)
        return vars.n_gate_r * vars.delta_r

    @staticmethod
    def T_from_Delta_r():
        """Calculate the time required (T) to achieve range offset (Delta_r_max) at constant acceleration (alpha).
        Uses vars.Delta_r_max (m) and vars.alpha (m/s^2).
        Returns: float: Time duration (s)."""
        # T = sqrt((2 * Delta_r_max) / alpha)
        if vars.alpha <= 0:
            raise ValueError("Acceleration (alpha) must be positive.")
        return math.sqrt((2 * vars.Delta_r_max) / vars.alpha)

    @staticmethod
    def Delta_v_max_gate():
        """Calculate the maximum required velocity offset to exit the gate.
        Uses vars.n_gate_v (cells) and vars.rho_v (m/s).
        Returns: float: Maximum velocity offset (m/s)."""
        # Delta_v_max = n_gate_v * rho_v (Requires rho_v to be set)
        return vars.n_gate_v * vars.rho_v

    @staticmethod
    def T_from_Delta_v():
        """Calculate the time required (T) to achieve velocity offset (Delta_v_max) at constant acceleration (a_accel).
        Uses vars.Delta_v_max (m/s) and vars.a_accel (m/s^2).
        Returns: float: Time duration (s)."""
        # T = Delta_v_max / a_accel
        if vars.a_accel == 0:
            raise ValueError("Acceleration (a_accel) cannot be zero.")
        return vars.Delta_v_max / vars.a_accel

    # =========================================================================
    # TOPIC 18: CROSS-EYE SOLVERS
    # =========================================================================

    @staticmethod
    def phi_hat_cross_eye_amp():
        """Calculate the apparent cross-eye angle error (amplitude monopulse approximation).
        Uses vars.L_cross (m), vars.R (m), and vars.a_gain_ratio (J1/J2, dimensionless).
        Returns: float: Angle error in radians."""
        # phi_hat = (L/(2R)) * (1+a)/(1-a)
        if vars.R == 0:
             raise ValueError("Range (R) cannot be zero.")
        if abs(vars.a_gain_ratio - 1.0) < 1e-9: # Check for a very close to 1
            # The denominator (1-a) approaches zero, causing the angle error to approach infinity (saturating the tracker)
            return np.inf 
        return (vars.L_cross / (2 * vars.R)) * ((1 + vars.a_gain_ratio) / (1 - vars.a_gain_ratio))

    @staticmethod
    def L_cross_from_phi_hat():
        """Calculate the required cross-eye aperture separation (L) to produce a specific angle error (phi_hat_ce).
        Uses vars.phi_hat_ce (rad), vars.R (m), and vars.a_gain_ratio (J1/J2).
        Returns: float: Aperture separation (m)."""
        # Rearrange: L = 2 * R * phi_hat_ce * (1-a)/(1+a)
        if abs(vars.a_gain_ratio + 1.0) < 1e-9: # Check for a very close to -1
            return np.inf # Denominator approaches zero
        return 2 * vars.R * vars.phi_hat_ce * ((1 - vars.a_gain_ratio) / (1 + vars.a_gain_ratio))

class convert:  # add alias con for convenience
    """Unit conversion utilities for radar calculations.
    
    This class provides a comprehensive set of unit conversion functions for
    angles, power, frequency, distance, and signal levels commonly used in
    radar engineering.
    
    Available conversions:
        - Angles: rad_to_deg, deg_to_rad
        - Power/Signal: lin_to_db, db_to_lin
        - Distance: m_to_mi, mi_to_m, m_to_nmi, nmi_to_m, ft_to_m, m_to, m_from
        - Power: w_to, w_from
        - Frequency: hz_to, hz_from
    
    Example:
        >>> from radar_range_equation import convert
        >>> deg = convert.rad_to_deg(3.14159)
        >>> print(f"{deg} degrees")
    """
    pass
    
    def rad_to_deg(value_radians):
        """Convert radians to degrees.
        
        Args:
            value_radians (float): Angle in radians.
        
        Returns:
            float: Angle in degrees.
        
        Example:
            >>> from radar_range_equation import convert
            >>> deg = convert.rad_to_deg(3.14159)
        """
        return value_radians * (180 / np.pi)
    
    def deg_to_rad(value_degrees):
        """Convert degrees to radians.
        
        Args:
            value_degrees (float): Angle in degrees.
        
        Returns:
            float: Angle in radians.
        
        Example:
            >>> from radar_range_equation import convert
            >>> rad = convert.deg_to_rad(180)
        """
        return value_degrees * (np.pi / 180)
    
    def lin_to_db(value_linear):
        """Convert linear value to decibels (dB).
        
        Args:
            value_linear (float): Linear value (must be positive).
        
        Returns:
            float: Value in dB. Returns -inf for non-positive inputs.
        
        Example:
            >>> from radar_range_equation import convert
            >>> db = convert.lin_to_db(10)  # Returns 10 dB
        """
        # Handle non-positive inputs to avoid math errors
        if value_linear <= 0:
            return -np.inf
        return np.log(value_linear)/np.log(10)*10
    
    def db_to_lin(value_db):
        """Convert decibels (dB) to linear value.
        
        Args:
            value_db (float): Value in dB.
        
        Returns:
            float: Linear value.
        
        Example:
            >>> from radar_range_equation import convert
            >>> lin = convert.db_to_lin(10)  # Returns 10.0
        """
        return 10**(value_db / 10.0)
    
    def m_to_mi(value_meters):
        """Convert meters to miles.
        
        Args:
            value_meters (float): Distance in meters.
        
        Returns:
            float: Distance in miles.
        
        Example:
            >>> from radar_range_equation import convert
            >>> miles = convert.m_to_mi(1609.34)  # Returns ~1 mile
        """
        return value_meters / 1609.34
    
    def mi_to_m(value_miles):
        """Convert miles to meters.
        
        Args:
            value_miles (float): Distance in miles.
        
        Returns:
            float: Distance in meters.
        
        Example:
            >>> from radar_range_equation import convert
            >>> meters = convert.mi_to_m(1)  # Returns 1609.34 meters
        """
        return value_miles * 1609.34
        
    def nmi_to_m(value_nmi):
        """Convert nautical miles to meters.
        
        Args:
            value_nmi (float): Distance in nautical miles.
        
        Returns:
            float: Distance in meters.
        
        Example:
            >>> from radar_range_equation import convert
            >>> meters = convert.nmi_to_m(1)  # Returns 1852 meters
        """
        return value_nmi * 1852.0

    def m_to_nmi(value_m):
        """Convert meters to nautical miles.
        
        Args:
            value_m (float): Distance in meters.
        
        Returns:
            float: Distance in nautical miles.
        
        Example:
            >>> from radar_range_equation import convert
            >>> nmi = convert.m_to_nmi(1852)  # Returns ~1 nautical mile
        """
        return value_m / 1852.0

    def w_to(value_w, target_unit):
        """Convert watts to other power units.
        
        Args:
            value_w (float): Power in watts.
            target_unit (str): Target unit ('kw' for kilowatts, 'mw' for milliwatts).
        
        Returns:
            float: Power in the target unit.
        
        Example:
            >>> from radar_range_equation import convert
            >>> kw = convert.w_to(1000, 'kw')  # Returns 1.0 kW
        """
        conversion_factors = {
            'kw': 1 / 1000.0,
            'mw': 1 * 10**6
        }
        return value_w * conversion_factors.get(target_unit.lower(), 1)

    def w_from(value, source_unit):
        """Convert from other power units to watts.
        
        Args:
            value (float): Power in the source unit.
            source_unit (str): Source unit ('kw' for kilowatts, 'mw' for milliwatts).
        
        Returns:
            float: Power in watts.
        
        Example:
            >>> from radar_range_equation import convert
            >>> watts = convert.w_from(1, 'kw')  # Returns 1000.0 W
        """
        conversion_factors = {
            'kw': 1 / 1000.0,
            'mw': 1 * 10**6
        }
        return value / conversion_factors.get(source_unit.lower(), 1)
    
    def ft_to_m(value_feet):
        """Convert feet to meters.
        
        Args:
            value_feet (float): Distance in feet.
        
        Returns:
            float: Distance in meters.
        
        Example:
            >>> from radar_range_equation import convert
            >>> meters = convert.ft_to_m(1)  # Returns 0.3048 meters
        """
        return value_feet * 0.3048

    def hz_to(value_hz, target_unit):
        """Convert hertz to other frequency units.
        
        Args:
            value_hz (float): Frequency in hertz.
            target_unit (str): Target unit ('khz', 'mhz', or 'ghz').
        
        Returns:
            float: Frequency in the target unit.
        
        Example:
            >>> from radar_range_equation import convert
            >>> ghz = convert.hz_to(10e9, 'ghz')  # Returns 10.0 GHz
        """
        conversion_factors = {
            'khz': 1 * 10**3,
            'mhz': 1 * 10**6,
            'ghz': 1 * 10**9
        }
        return value_hz / conversion_factors.get(target_unit.lower(), 1)
    
    def hz_from(value, source_unit):
        """Convert from other frequency units to hertz.
        
        Args:
            value (float): Frequency in the source unit.
            source_unit (str): Source unit ('khz', 'mhz', or 'ghz').
        
        Returns:
            float: Frequency in hertz.
        
        Example:
            >>> from radar_range_equation import convert
            >>> hz = convert.hz_from(10, 'ghz')  # Returns 10e9 Hz
        """
        conversion_factors = {
            'khz': 1 * 10**3,
            'mhz': 1 * 10**6,
            'ghz': 1 * 10**9
        }
        return value * conversion_factors.get(source_unit.lower(), 1)
    
    def m_to(value_meters, target_unit):
        """Convert meters to other distance units.
        
        Args:
            value_meters (float): Distance in meters.
            target_unit (str): Target unit ('km' for kilometers).
        
        Returns:
            float: Distance in the target unit.
        
        Example:
            >>> from radar_range_equation import convert
            >>> km = convert.m_to(1000, 'km')  # Returns 1.0 km
        """
        conversion_factors = {
            'km': 1 / 1000.0
        }
        return value_meters * conversion_factors.get(target_unit.lower(), 1)

    def m_from(value, source_unit):
        """Convert from other distance units to meters.
        
        Args:
            value (float): Distance in the source unit.
            source_unit (str): Source unit ('km', 'nmi', or 'mi').
        
        Returns:
            float: Distance in meters.
        
        Example:
            >>> from radar_range_equation import convert
            >>> meters = convert.m_from(1, 'km')  # Returns 1000.0 meters
        """
        conversion_factors = {
            'km': 1000.0,
            'nmi': 1852.0,
            'mi': 1609.34
        }
        return value * conversion_factors.get(source_unit.lower(), 1)
    
con = convert()  # alias for convenience


class analysis:
    """Analysis helpers: pulse parsing, integration gains, fluctuation loss,
    jammer and burnthrough calculations.

    These convenience functions wrap common calculation patterns used in the
    quiz solutions and examples. They use `vars`, `convert` and numpy utilities
    already defined in this module.
    """

    @staticmethod
    def pulse_timing(pulse_intervals_us):
        """Compute pulse width, PRI, PRF and duty cycle.

        Args:
            pulse_intervals_us (list of (start_us, end_us)): list of ON intervals in microseconds

        Returns:
            dict: {pulse_width_s, PRI_s, PRF_hz, duty_cycle, n_p}
        """
        if len(pulse_intervals_us) == 0:
            return {'pulse_width_s': 0.0, 'PRI_s': 0.0, 'PRF_hz': float('inf'), 'duty_cycle': 0.0, 'n_p': 0}
        widths_us = [end - start for (start, end) in pulse_intervals_us]
        pulse_width_us = widths_us[0]
        starts = [s for (s, e) in pulse_intervals_us]
        if len(starts) > 1:
            PRI_us = starts[1] - starts[0]
        else:
            PRI_us = 0.0
        pulse_width_s = pulse_width_us * 1e-6
        PRI_s = PRI_us * 1e-6 if PRI_us > 0 else 0.0
        PRF_hz = 1.0 / PRI_s if PRI_s > 0 else float('inf')
        duty_cycle = pulse_width_s / PRI_s if PRI_s > 0 else 0.0
        n_p = len(pulse_intervals_us)
        return {'pulse_width_s': pulse_width_s, 'PRI_s': PRI_s, 'PRF_hz': PRF_hz, 'duty_cycle': duty_cycle, 'n_p': n_p}

    @staticmethod
    def power_and_integration(amplitude, pulse_intervals_us, use='coherent'):
        """Compute peak power, average power and integration SNR gain.

        Args:
            amplitude (float): amplitude A of the transmit pulse (arbitrary units)
            pulse_intervals_us (list): list of (start,end) intervals in us
            use (str): 'coherent' or 'noncoherent' for integration gain model

        Returns:
            dict: {P_peak, P_avg, n_p, gain_linear, gain_db}
        """
        P_peak = (amplitude ** 2) / 2.0
        t = analysis.pulse_timing(pulse_intervals_us)
        P_avg = P_peak * t['duty_cycle']
        n = t['n_p']
        if n <= 0:
            gain_linear = 1.0
        else:
            if use == 'coherent':
                gain_linear = float(n)
            else:
                # simple non-coherent approximation (amplitude averaging)
                gain_linear = float(np.sqrt(n))
        gain_db = convert.lin_to_db(gain_linear)
        return {'P_peak': P_peak, 'P_avg': P_avg, 'n_p': n, 'gain_linear': gain_linear, 'gain_db': gain_db}

    @staticmethod
    def effective_independent_looks(pulse_start_times_us, correlation_time_us):
        """Estimate number of effective independent looks given correlation time.

        Groups pulses whose start times are within `correlation_time_us` of any
        pulse already in a group. Returns the clusters and n_e (number of 
        independent looks).
        
        Args:
            pulse_start_times_us (list): List of pulse start times in microseconds
            correlation_time_us (float): Correlation time in microseconds
            
        Returns:
            dict: {'clusters': list of clusters, 'n_e': number of independent looks}
        """
        clusters = []
        for t in pulse_start_times_us:
            placed = False
            for cl in clusters:
                # Check if this pulse is within correlation time of ANY pulse in cluster
                if any(abs(t - existing_t) < correlation_time_us for existing_t in cl):
                    cl.append(t)
                    placed = True
                    break
            if not placed:
                clusters.append([t])
        return {'clusters': clusters, 'n_e': len(clusters)}

    @staticmethod
    def fluctuation_loss_total(L_single_db, n_e):
        """Combine single-look fluctuation loss across n_e independent looks.

        Simple method (default): L_total_db = L_single_db - 10*log10(n_e)
        (This is an approximate power-averaging rule; actual behavior depends on Swerling model.)
        """
        if n_e <= 0:
            return L_single_db
        return L_single_db - 10.0 * np.log10(n_e)

    @staticmethod
    def monostatic_received_power(Pt, Gt_linear, Gr_linear, wavelength_m, sigma_m2, R_m):
        """Monostatic single-pulse received power (Friis-based radar equation).

        Pr = Pt * Gt * Gr * lambda^2 * sigma / ((4*pi)^3 * R^4)
        """
        return (Pt * Gt_linear * Gr_linear * (wavelength_m ** 2) * sigma_m2) / ((4.0 * np.pi) ** 3 * (R_m ** 4))

    @staticmethod
    def processed_signal_power(Pr_single, n_p, PCR_db):
        """Apply coherent integration (n_p) and pulse compression (PCR_db) to single-pulse power."""
        PCR_lin = convert.db_to_lin(PCR_db)
        return Pr_single * float(n_p) * PCR_lin

    @staticmethod
    def jammer_received_power(Pj, G_j_tx_db, G_radar_on_jam_db, wavelength_m, R_jr_m, band_ratio=1.0):
        """Compute jammer power received in the radar band.

        Assumes free-space propagation and that only 1/band_ratio of jammer power
        falls into the radar bandwidth.
        """
        Gjt = convert.db_to_lin(G_j_tx_db)
        Grj = convert.db_to_lin(G_radar_on_jam_db)
        return Pj * Gjt * Grj * (wavelength_m ** 2) / ((4.0 * np.pi * R_jr_m) ** 2) * (1.0 / band_ratio)

    @staticmethod
    def signal_to_jammer(S_processed, P_jr):
        """Return S/J linear and dB given processed signal power and jammer received power."""
        if P_jr <= 0:
            return {'SJR_linear': np.inf, 'SJR_db': np.inf}
        s = S_processed / P_jr
        return {'SJR_linear': s, 'SJR_db': convert.lin_to_db(s)}

    @staticmethod
    def burnthrough_range(Pj, G_j_tx_db, G_radar_on_jam_db, wavelength_m, band_ratio, S_processed):
        """Solve for jammer range R_jr where processed signal equals jammer power in band.

        Returns range in meters.
        """
        Gjt = convert.db_to_lin(G_j_tx_db)
        Grj = convert.db_to_lin(G_radar_on_jam_db)
        numer = Pj * Gjt * Grj * (wavelength_m ** 2)
        den = (4.0 * np.pi) ** 2 * (1.0 / band_ratio) * S_processed
        if den <= 0:
            return np.inf
        return float(np.sqrt(numer / den))


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
    
    # --- Original Demo Code ---
    print("--- Original Demo (CWFM) ---")
    v = vars()  # Create a reference to the global vars instance
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

    # =========================================================================
    # Topic 12: Chaff Demo
    # =========================================================================
    print("\n" + "="*30)
    print("Topic 12: Chaff Demo")
    print("="*30)

    # Example 1: 10 GHz radar
    redefine_variable('f', 10e9)
    redefine_variable('wavelength', solve.wavelength())
    redefine_variable('L_fiber', solve.L_fiber())
    # Fiber diameter is 1/1000th of length
    redefine_variable('D_fiber', v.L_fiber / 1000) 
    redefine_variable('V_ch', solve.V_ch())

    # Box size: 5cm x 5cm x 2cm, 60% fill ratio
    redefine_variable('V_box', 0.05 * 0.05 * 0.02)
    redefine_variable('Fill_ratio', 0.60)
    redefine_variable('N_fiber', solve.N_fiber())

    print(f"Fiber Length (L_fiber): {v.L_fiber*1e3:.2f} mm")
    print(f"Single Fiber Volume (V_ch): {v.V_ch:.2e} m^3")
    print(f"Total Fibers (N_fiber): {v.N_fiber:.2e} fibers")

    # RCS at T=infinity
    redefine_variable('sigma_ch', 0.15 * v.N_fiber * v.wavelength**2)
    print(f"Max Chaff RCS (sigma_ch_max): {v.sigma_ch:.0f} m^2")
    print(f"Max Chaff RCS (dBsm): {convert.lin_to_db(v.sigma_ch):.2f} dBsm")

    # =========================================================================
    # Topic 13: Noise Jamming Demo
    # =========================================================================
    print("\n" + "="*30)
    print("Topic 13: Noise Jamming Demo (Burnthrough)")
    print("="*30)

    # Based on Topic 13 Example 1
    redefine_variable('f', 12e9)
    redefine_variable('wavelength', solve.wavelength())
    redefine_variable('R', 10e3)              # Range 10 km
    redefine_variable('Pt', 1500.0)           # P_t = 1500 W
    redefine_variable('Gt', convert.db_to_lin(33)) # G_t = 33 dB
    redefine_variable('sigma', convert.db_to_lin(13)) # sigma = 13 dBsm
    redefine_variable('N_p', 1.0)             # N_p = 1
    redefine_variable('B', 67e6)              # B_r = 67 MHz (from 15 ns pulse width)
    redefine_variable('S_min', convert.db_to_lin(10)) # Assume 10dB S/N required
    redefine_variable('Pj', 500.0)            # P_j = 500 W
    redefine_variable('Gj', convert.db_to_lin(3)) # G_j = 3 dB
    redefine_variable('Bj', 6e9)              # B_j = 6 GHz
    redefine_variable('Lossj', convert.db_to_lin(2)) # Loss_j = 2 dB

    redefine_variable('S_J_ratio', solve.S_J_ratio())
    print(f"S/J Ratio: {v.S_J_ratio:.3f} (Linear)")
    print(f"S/J Ratio: {convert.lin_to_db(v.S_J_ratio):.2f} dB")

    redefine_variable('R_bt', solve.R_burnthrough())
    print(f"Burnthrough Range (R_bt): {v.R_bt / 1e3:.2f} km")

    # =========================================================================
    # Topic 14: Gated Noise Demo
    # =========================================================================
    print("\n" + "="*30)
    print("Topic 14: Gated Noise Demo (Timing)")
    print("="*30)

    # Based on Topic 14 Example 1
    redefine_variable('R_tgt', 50e3)          # Range 50 km
    redefine_variable('tau', 30e-6)           # Pulse Width 30 us
    redefine_variable('R_gn_start_offset', 600) # Start 600m before target
    
    redefine_variable('t_tgt_2way', solve.t_tgt_2way())
    redefine_variable('t_gn_start_release', solve.t_gn_start_release())

    print(f"Target 2-way time (t_tgt_2way): {v.t_tgt_2way*1e6:.2f} µs")
    print(f"Noise Start Release Time (t_gn_start_release): {v.t_gn_start_release*1e6:.2f} µs")

    # =========================================================================
    # Topic 15: False Target Generation Demo
    # =========================================================================
    print("\n" + "="*30)
    print("Topic 15: False Target Generation Demo")
    print("="*30)
    
    # Based on Topic 15 Example 1
    redefine_variable('f', 10e9)
    redefine_variable('wavelength', solve.wavelength())
    redefine_variable('R', 50e3)              # Target Range 50 km
    redefine_variable('v_tgt', -120.0)        # Target range rate -120 m/s
    redefine_variable('R_ft', 60e3)           # False Target Range 60 km
    redefine_variable('v_ft', 60.0)           # False Target range rate 60 m/s

    redefine_variable('f_D_tgt', solve.f_D_tgt())
    redefine_variable('f_D_ft', solve.f_D_ft())

    redefine_variable('Delta_t_ft', solve.Delta_t_ft())
    redefine_variable('Delta_f_ft', solve.Delta_f_ft())

    print(f"Target Doppler (f_D_tgt): {v.f_D_tgt:.0f} Hz")
    print(f"False Target Doppler (f_D_ft): {v.f_D_ft:.0f} Hz")
    print(f"Required Time Delay (Delta_t_ft): {v.Delta_t_ft*1e6:.2f} µs")
    print(f"Required Freq Shift (Delta_f_ft): {v.Delta_f_ft/1e3:.1f} kHz")

    # =========================================================================
    # Topic 16: Radar Tracking / False Tracks Demo
    # =========================================================================
    print("\n" + "="*30)
    print("Topic 16: False Tracks Demo (Jammer Power)")
    print("="*30)

    # Based on Topic 16 Example 1
    redefine_variable('f', 10e9)
    redefine_variable('R', 50e3)              # Range 50 km
    redefine_variable('Pt', 3000.0)           # P_t = 3 kW
    redefine_variable('Gt', convert.db_to_lin(25)) # G_t = 25 dB
    redefine_variable('sigma', convert.db_to_lin(10)) # Target RCS 10 dBsm
    redefine_variable('Gj', convert.db_to_lin(3)) # Jammer Gain 3 dB

    redefine_variable('P_density', solve.P_density())
    print(f"Power Density at Target: {v.P_density:.2e} W/m^2")

    redefine_variable('Pj_emulated', solve.Pj_emulated())
    print(f"Jammer Power to Emulate RCS: {v.Pj_emulated*1e3:.2f} mW")

    # =========================================================================
    # Topic 17: Gate Stealing Demo
    # =========================================================================
    print("\n" + "="*30)
    print("Topic 17: Gate Stealing Demo (Range and Velocity)")
    print("="*30)
    
    # 1. Range Stealing Example (Range Gate Pull-Off)
    print("\n--- Range Gate Pull-Off (RGPO) ---")
    redefine_variable('f', 3e9)                         # 3 GHz
    redefine_variable('wavelength', solve.wavelength())
    redefine_variable('alpha', 2 * v.g)                  # Target Accel = 2g
    redefine_variable('tau', 20e-6)                      # 20 us
    redefine_variable('gamma', 0.5e12)                   # 0.5 MHz/us = 0.5e12 Hz/s
    redefine_variable('B', solve.B_chirp())
    redefine_variable('delta_r', solve.delta_r_compressed()) # Compressed Range Res

    print(f"Accel (alpha): {v.alpha:.2f} m/s^2")
    print(f"Compressed Range Res (delta_r): {v.delta_r:.2f} m")

    # Required Max Range Offset of at least 200m (Delta_r_max = 200m)
    redefine_variable('Delta_r_max', 200) # Directly set the requirement
    
    redefine_variable('T_time', solve.T_from_Delta_r())
    print(f"Time to achieve 200m offset (T): {v.T_time:.2f} s")
    print(f"Range Offset Profile: Delta_r(t) = 0.5 * {v.alpha:.2f} * t^2")
    
    # 2. Velocity Stealing Example (Velocity Gate Pull-Off)
    print("\n--- Velocity Gate Pull-Off (VGPO) ---")
    redefine_variable('f', 12e9)
    redefine_variable('wavelength', solve.wavelength())
    redefine_variable('T_cpi', 60e-3)                    # 60 ms
    redefine_variable('a_accel', 1.5 * v.g)              # Target Accel = 1.5g

    redefine_variable('rho_v', solve.rho_v())
    print(f"Doppler Bin Size (rho_v): {v.rho_v:.3f} m/s")

    # Max velocity offset is two gate widths (2 * 9 cells * rho_v)
    redefine_variable('n_gate_v', 9) 
    redefine_variable('Delta_v_max', 2 * v.n_gate_v * v.rho_v) # Max offset is two gate widths
    print(f"Max Velocity Offset (Delta_v_max): {v.Delta_v_max:.2f} m/s")
    
    redefine_variable('T_time', solve.T_from_Delta_v())
    print(f"Time to achieve max offset (T): {v.T_time:.2f} s")
    print(f"Velocity Offset Profile: Delta_v(t) = {v.a_accel:.2f} * t")

    # =========================================================================
    # Topic 18: Cross-Eye Demo
    # =========================================================================
    print("\n" + "="*30)
    print("Topic 18: Cross-Eye Demo")
    print("="*30)
    
    # Cross-Eye Example (based on Topic 18, Prob 1b)
    # Required parameters: R, phi_hat_ce, a_gain_ratio
    redefine_variable('R', 30e3)                         # Target Range (30 km)
    redefine_variable('a_gain_ratio', -0.9)              # J1/J2 = -0.9
    
    # Assume the required angle error is 0.0035 rad (5 times the example's sigma_phi)
    redefine_variable('phi_hat_ce', 0.0035) 

    redefine_variable('L_cross', solve.L_cross_from_phi_hat())
    print(f"Target Range (R): {v.R / 1e3:.0f} km")
    print(f"Gain Ratio (a_ratio): {v.a_gain_ratio}")
    print(f"Required Aperture Separation (L_cross): {v.L_cross:.2f} m") 
    
    # Reverse check: Calculate the angle error for a large aperture (L=138m from solution)
    redefine_variable('L_cross', 138.0) 
    redefine_variable('phi_hat_ce', solve.phi_hat_cross_eye_amp())
    print(f"\nGiven L_cross=138m, Angle Error (phi_hat_ce): {v.phi_hat_ce * 1e3:.2f} mrad ({convert.rad_to_deg(v.phi_hat_ce):.4f} deg)")