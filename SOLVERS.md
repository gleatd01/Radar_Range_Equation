## Solver Reference
This section contains all solver functions and Solvable instances available in the `radar_range_equation.solve` module.

**Total Solvers**: 76

### Base/Common
#### `A_e_circ` (Function)
Calculate effective aperture for a circular antenna.

Uses vars.eta (antenna efficiency) and vars.D (antenna diameter in m).

Returns:
    float: Effective aperture in m².

Example:
    >>> import radar_range_equation as RRE
    >>> RRE.vars.eta = 0.6
    >>> RRE.vars.D = 2.0
    >>> aperture = RRE.solve.A_e_circ()

#### `A_e_rect` (Function)
Calculate effective aperture for a rectangular antenna.

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

#### `A_sphere` (Function)
Calculate the radius of a sphere from its radar cross section.

Uses vars.sigma (radar cross section in m²).

Returns:
    float: Radius of the sphere in meters.

Example:
    >>> import radar_range_equation as RRE
    >>> RRE.vars.sigma = 3.14159
    >>> radius = RRE.solve.A_sphere()

#### `G_t` (Solvable)
A dynamic, introspective solver for a single radar parameter.

    This class encapsulates the logic to solve for a single parameter using
    one or more symbolic equations. It can report its status (what's needed
    to solve) and perform the calculation when all required inputs are available.

    The Solvable class interacts with the vars module to check which variables
    are defined and uses sympy to perform symbolic and numeric calculations.

    Attributes:
        target_symbol (sympy.Symbol): The variable to be solved for
        equation_list (list): List of sympy.Eq equations that can solve for the target

    Example:
        >>> import radar_range_equation as RRE
        >>> # In interactive console: RRE.solve.G_t shows status
        >>> # RRE.solve.G_t() executes calculation


#### `G_t_func` (Function)
Calculate transmit antenna gain (legacy direct calculation).

Uses vars.pi (pi constant), vars.A_e (effective aperture in m²),
and vars.wavelength (wavelength in m).

Returns:
    float: Antenna gain (dimensionless).

Example:
    >>> import radar_range_equation as RRE
    >>> RRE.vars.A_e = 1.0
    >>> RRE.vars.wavelength = 0.03
    >>> gain = RRE.solve.G_t_func()

#### `P_t` (Solvable)
A dynamic, introspective solver for a single radar parameter.

    This class encapsulates the logic to solve for a single parameter using
    one or more symbolic equations. It can report its status (what's needed
    to solve) and perform the calculation when all required inputs are available.

    The Solvable class interacts with the vars module to check which variables
    are defined and uses sympy to perform symbolic and numeric calculations.

    Attributes:
        target_symbol (sympy.Symbol): The variable to be solved for
        equation_list (list): List of sympy.Eq equations that can solve for the target

    Example:
        >>> import radar_range_equation as RRE
        >>> # In interactive console: RRE.solve.G_t shows status
        >>> # RRE.solve.G_t() executes calculation


#### `R4` (Function)
Calculate R^4 from the radar range equation.

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

#### `R_max` (Solvable)
A dynamic, introspective solver for a single radar parameter.

    This class encapsulates the logic to solve for a single parameter using
    one or more symbolic equations. It can report its status (what's needed
    to solve) and perform the calculation when all required inputs are available.

    The Solvable class interacts with the vars module to check which variables
    are defined and uses sympy to perform symbolic and numeric calculations.

    Attributes:
        target_symbol (sympy.Symbol): The variable to be solved for
        equation_list (list): List of sympy.Eq equations that can solve for the target

    Example:
        >>> import radar_range_equation as RRE
        >>> # In interactive console: RRE.solve.G_t shows status
        >>> # RRE.solve.G_t() executes calculation


#### `db_to_linear` (Function)
Converts SNR from dB to linear.

Uses vars.x (as the dB value)

Returns:
    Linear value

#### `sigma_sphere` (Function)
Calculate the radar cross section of a sphere from its area.

Uses vars.A (antenna area in m²).

Returns:
    float: Radar cross section in m².

Example:
    >>> import radar_range_equation as RRE
    >>> RRE.vars.A = 1.0
    >>> rcs = RRE.solve.sigma_sphere()

#### `theta_B` (Function)
Calculate the 3-dB beamwidth using Gaussian approximation.

Uses vars.wavelength (wavelength in m) and vars.D_h (horizontal antenna 
dimension in m). Uses the approximation: theta_B = 65° * pi/180 * lambda/D_h.

Returns:
    float: Beamwidth in radians.

Example:
    >>> import radar_range_equation as RRE
    >>> RRE.vars.wavelength = 0.03
    >>> RRE.vars.D_h = 1.0
    >>> beamwidth = RRE.solve.theta_B()

#### `wavelength` (Solvable)
A dynamic, introspective solver for a single radar parameter.

    This class encapsulates the logic to solve for a single parameter using
    one or more symbolic equations. It can report its status (what's needed
    to solve) and perform the calculation when all required inputs are available.

    The Solvable class interacts with the vars module to check which variables
    are defined and uses sympy to perform symbolic and numeric calculations.

    Attributes:
        target_symbol (sympy.Symbol): The variable to be solved for
        equation_list (list): List of sympy.Eq equations that can solve for the target

    Example:
        >>> import radar_range_equation as RRE
        >>> # In interactive console: RRE.solve.G_t shows status
        >>> # RRE.solve.G_t() executes calculation


#### `wavelength_func` (Function)
Calculate wavelength from frequency (legacy direct calculation).

Uses vars.c (speed of light in m/s) and vars.f (frequency in Hz).

Returns:
    float: Wavelength in meters.

Example:
    >>> import radar_range_equation as RRE
    >>> RRE.vars.c = 3e8
    >>> RRE.vars.f = 10e9
    >>> wl = RRE.solve.wavelength_func()


### Doppler CW Radar
#### `delta_v` (Solvable)
A dynamic, introspective solver for a single radar parameter.

    This class encapsulates the logic to solve for a single parameter using
    one or more symbolic equations. It can report its status (what's needed
    to solve) and perform the calculation when all required inputs are available.

    The Solvable class interacts with the vars module to check which variables
    are defined and uses sympy to perform symbolic and numeric calculations.

    Attributes:
        target_symbol (sympy.Symbol): The variable to be solved for
        equation_list (list): List of sympy.Eq equations that can solve for the target

    Example:
        >>> import radar_range_equation as RRE
        >>> # In interactive console: RRE.solve.G_t shows status
        >>> # RRE.solve.G_t() executes calculation


#### `f_doppler` (Solvable)
A dynamic, introspective solver for a single radar parameter.

    This class encapsulates the logic to solve for a single parameter using
    one or more symbolic equations. It can report its status (what's needed
    to solve) and perform the calculation when all required inputs are available.

    The Solvable class interacts with the vars module to check which variables
    are defined and uses sympy to perform symbolic and numeric calculations.

    Attributes:
        target_symbol (sympy.Symbol): The variable to be solved for
        equation_list (list): List of sympy.Eq equations that can solve for the target

    Example:
        >>> import radar_range_equation as RRE
        >>> # In interactive console: RRE.solve.G_t shows status
        >>> # RRE.solve.G_t() executes calculation


#### `f_obs_if` (Solvable)
A dynamic, introspective solver for a single radar parameter.

    This class encapsulates the logic to solve for a single parameter using
    one or more symbolic equations. It can report its status (what's needed
    to solve) and perform the calculation when all required inputs are available.

    The Solvable class interacts with the vars module to check which variables
    are defined and uses sympy to perform symbolic and numeric calculations.

    Attributes:
        target_symbol (sympy.Symbol): The variable to be solved for
        equation_list (list): List of sympy.Eq equations that can solve for the target

    Example:
        >>> import radar_range_equation as RRE
        >>> # In interactive console: RRE.solve.G_t shows status
        >>> # RRE.solve.G_t() executes calculation


#### `v_from_doppler` (Solvable)
A dynamic, introspective solver for a single radar parameter.

    This class encapsulates the logic to solve for a single parameter using
    one or more symbolic equations. It can report its status (what's needed
    to solve) and perform the calculation when all required inputs are available.

    The Solvable class interacts with the vars module to check which variables
    are defined and uses sympy to perform symbolic and numeric calculations.

    Attributes:
        target_symbol (sympy.Symbol): The variable to be solved for
        equation_list (list): List of sympy.Eq equations that can solve for the target

    Example:
        >>> import radar_range_equation as RRE
        >>> # In interactive console: RRE.solve.G_t shows status
        >>> # RRE.solve.G_t() executes calculation



### CWFM Radar
#### `R_cwfm` (Solvable)
A dynamic, introspective solver for a single radar parameter.

    This class encapsulates the logic to solve for a single parameter using
    one or more symbolic equations. It can report its status (what's needed
    to solve) and perform the calculation when all required inputs are available.

    The Solvable class interacts with the vars module to check which variables
    are defined and uses sympy to perform symbolic and numeric calculations.

    Attributes:
        target_symbol (sympy.Symbol): The variable to be solved for
        equation_list (list): List of sympy.Eq equations that can solve for the target

    Example:
        >>> import radar_range_equation as RRE
        >>> # In interactive console: RRE.solve.G_t shows status
        >>> # RRE.solve.G_t() executes calculation


#### `f_0_cwfm` (Solvable)
A dynamic, introspective solver for a single radar parameter.

    This class encapsulates the logic to solve for a single parameter using
    one or more symbolic equations. It can report its status (what's needed
    to solve) and perform the calculation when all required inputs are available.

    The Solvable class interacts with the vars module to check which variables
    are defined and uses sympy to perform symbolic and numeric calculations.

    Attributes:
        target_symbol (sympy.Symbol): The variable to be solved for
        equation_list (list): List of sympy.Eq equations that can solve for the target

    Example:
        >>> import radar_range_equation as RRE
        >>> # In interactive console: RRE.solve.G_t shows status
        >>> # RRE.solve.G_t() executes calculation


#### `f_d_cwfm` (Solvable)
A dynamic, introspective solver for a single radar parameter.

    This class encapsulates the logic to solve for a single parameter using
    one or more symbolic equations. It can report its status (what's needed
    to solve) and perform the calculation when all required inputs are available.

    The Solvable class interacts with the vars module to check which variables
    are defined and uses sympy to perform symbolic and numeric calculations.

    Attributes:
        target_symbol (sympy.Symbol): The variable to be solved for
        equation_list (list): List of sympy.Eq equations that can solve for the target

    Example:
        >>> import radar_range_equation as RRE
        >>> # In interactive console: RRE.solve.G_t shows status
        >>> # RRE.solve.G_t() executes calculation


#### `f_m_cwfm` (Solvable)
A dynamic, introspective solver for a single radar parameter.

    This class encapsulates the logic to solve for a single parameter using
    one or more symbolic equations. It can report its status (what's needed
    to solve) and perform the calculation when all required inputs are available.

    The Solvable class interacts with the vars module to check which variables
    are defined and uses sympy to perform symbolic and numeric calculations.

    Attributes:
        target_symbol (sympy.Symbol): The variable to be solved for
        equation_list (list): List of sympy.Eq equations that can solve for the target

    Example:
        >>> import radar_range_equation as RRE
        >>> # In interactive console: RRE.solve.G_t shows status
        >>> # RRE.solve.G_t() executes calculation


#### `f_r_cwfm` (Solvable)
A dynamic, introspective solver for a single radar parameter.

    This class encapsulates the logic to solve for a single parameter using
    one or more symbolic equations. It can report its status (what's needed
    to solve) and perform the calculation when all required inputs are available.

    The Solvable class interacts with the vars module to check which variables
    are defined and uses sympy to perform symbolic and numeric calculations.

    Attributes:
        target_symbol (sympy.Symbol): The variable to be solved for
        equation_list (list): List of sympy.Eq equations that can solve for the target

    Example:
        >>> import radar_range_equation as RRE
        >>> # In interactive console: RRE.solve.G_t shows status
        >>> # RRE.solve.G_t() executes calculation


#### `v_cwfm` (Solvable)
A dynamic, introspective solver for a single radar parameter.

    This class encapsulates the logic to solve for a single parameter using
    one or more symbolic equations. It can report its status (what's needed
    to solve) and perform the calculation when all required inputs are available.

    The Solvable class interacts with the vars module to check which variables
    are defined and uses sympy to perform symbolic and numeric calculations.

    Attributes:
        target_symbol (sympy.Symbol): The variable to be solved for
        equation_list (list): List of sympy.Eq equations that can solve for the target

    Example:
        >>> import radar_range_equation as RRE
        >>> # In interactive console: RRE.solve.G_t shows status
        >>> # RRE.solve.G_t() executes calculation



### Pulsed Radar
#### `R_from_time` (Function)
Dynamically created solver function.

#### `R_un_from_fp` (Function)
Dynamically created solver function.

#### `S_N_n_coherent_dB` (Function)
Calculates integrated SNR (coherent) in dB.

#### `S_N_n_noncoherent_dB` (Function)
Calculates integrated SNR (non-coherent, E_i=1/sqrt(n)) in dB.

#### `calculate_Theta` (Function)
Calculates the Theta parameter from the 3-dB beamwidth.

Uses vars.theta_B

Returns:
    Theta parameter

#### `fp_from_R_un` (Function)
Dynamically created solver function.

#### `tau_from_duty` (Function)
Dynamically created solver function.


### Direction Finding
#### `B_from_sigma_time` (Function)
Dynamically created solver function.

#### `L_cross_from_phi_hat` (Function)
Calculate the required cross-eye aperture separation (L) to produce a specific angle error (phi_hat_ce).
Uses vars.phi_hat_ce (rad), vars.R (m), and vars.a_gain_ratio (J1/J2).
Returns: float: Aperture separation (m).

#### `S_N_from_dB` (Function)
Dynamically created solver function.

#### `d_from_sigma_phase` (Function)
Dynamically created solver function.

#### `estimate_phi_hat` (Function)
Calculates the linear processor angle estimate.

Uses vars.Delta, vars.Sigma, vars.theta_B, vars.phi_s

Returns:
    Angle estimate in degrees

#### `phi_hat_cross_eye_amp` (Function)
Calculate the apparent cross-eye angle error (amplitude monopulse approximation).
Uses vars.L_cross (m), vars.R (m), and vars.a_gain_ratio (J1/J2, dimensionless).
Returns: float: Angle error in radians.

#### `phi_s_from_sigma_amp` (Function)
Dynamically created solver function.

#### `sigma_phi_amplitude` (Function)
Calculates the angle standard deviation for amplitude comparison.

Uses vars.theta_B, vars.S_N, vars.phi_s

Returns:
    Angle standard deviation in degrees

#### `sigma_phi_phase` (Function)
Calculates the angle standard deviation for phase comparison.

Uses vars.wavelength, vars.d, vars.S_N

Returns:
    Angle standard deviation in radians

#### `sigma_phi_time` (Function)
Calculates the angle standard deviation for time comparison.

Uses vars.c, vars.d, vars.B

Returns:
    Angle standard deviation in radians

#### `v_phi` (Function)
Calculates the Gaussian beam approximation.

Uses vars.phi, vars.phi_s, vars.Theta

Returns:
    Gaussian beam approximation value

#### `v_phi_full` (Function)
Calculates v(phi) using theta_B directly.

Uses vars.phi, vars.phi_s, vars.theta_B

Returns:
    Beam approximation value


### Pulse Compression
#### `B_chirp` (Function)
Dynamically created solver function.

#### `PCR_from_B` (Function)
Dynamically created solver function.

#### `PCR_from_gamma` (Function)
Dynamically created solver function.

#### `R_offset_from_tone` (Function)
Dynamically created solver function.

#### `delta_r_compressed` (Function)
Dynamically created solver function.

#### `delta_r_uncompressed` (Function)
Dynamically created solver function.

#### `f_range_tone` (Function)
Dynamically created solver function.


### Chaff
#### `L_fiber` (Function)
Calculate Chaff Fiber Length (lambda/2).
Uses vars.wavelength (m).
Returns: float: Fiber length (m).

#### `N_fiber` (Function)
Calculate Number of fibers in a cartridge.
Uses vars.V_box (m^3), vars.Fill_ratio (dim), and vars.V_ch (m^3).
Returns: float: Number of fibers (dimensionless).

#### `V_ch` (Function)
Calculate Volume of a single chaff fiber (cylinder).
Uses vars.L_fiber (m) and vars.D_fiber (m).
Returns: float: Fiber volume (m^3).

#### `sigma_ch_t` (Function)
Calculate RCS of Chaff Cloud at time t.
Uses vars.N_fiber (dim), vars.wavelength (m), vars.zeta_ch (s), and t_s (s).
Returns: float: RCS (m^2).


### Noise Jamming
#### `R_burnthrough` (Function)
Calculate Burnthrough Range (R_bt).
Uses Pt, Gt, sigma, n_p, Bj, Pj, Gj, Lossj, B, S_min.
Returns: float: Burnthrough Range (m).

#### `S_J_ratio` (Function)
Calculate Signal-to-Jammer (S/J) Ratio for Barrage Noise.
Uses Pt, Gt, sigma, n_p, Bj, R, Pj, Gj, Lossj, B.
Returns: float: S/J Ratio (linear).


### Gated Noise
#### `t_gn_start_release` (Function)
Calculate Gated Noise Start Release Time (relative to T=0 radar pulse).
Uses vars.R_tgt, vars.R_gn_start_offset, vars.c, vars.tau.
Returns: float: Time (s).

#### `t_tgt_2way` (Function)
Calculate Two-way Time of Flight to Target.
Uses vars.R_tgt (m) and vars.c (m/s).
Returns: float: Time (s).


### False Target Generation
#### `Delta_f_ft` (Function)
Calculate Frequency Shift to apply for False Target Generation.
Uses vars.f_D_ft (Hz) and vars.f_D_tgt (Hz).
Returns: float: Frequency shift (Hz).

#### `Delta_t_ft` (Function)
Calculate Time Delay to apply for False Target Generation.
Uses vars.R_ft (m), vars.R (m), and vars.c (m/s).
Returns: float: Time delay (s).

#### `f_D_ft` (Function)
Calculate False Target Doppler Frequency.
Uses vars.v_ft (m/s) and vars.wavelength (m).
Returns: float: Doppler frequency (Hz).

#### `f_D_tgt` (Function)
Calculate Target Doppler Frequency.
Uses vars.v_tgt (m/s) and vars.wavelength (m).
Returns: float: Doppler frequency (Hz).


### Radar Tracking
#### `P_density` (Function)
Calculate Power Density at Target.
Uses vars.P_t (W), vars.G_t (dim), and vars.R (m).
Returns: float: Power density (W/m^2).

#### `Pj_emulated` (Function)
Calculate Jammer Power required to emulate a target RCS (sigma).
Uses vars.P_density (W/m^2), vars.sigma (m^2), and vars.Gj (dim).
Returns: float: Jammer transmit power (W).


### Gate Stealing
#### `Delta_r_max_gate` (Function)
Calculate the maximum required range offset to exit the gate.
Uses vars.n_gate_r (cells) and vars.delta_r (m).
Returns: float: Maximum range offset (m).

#### `Delta_v_max_gate` (Function)
Calculate the maximum required velocity offset to exit the gate.
Uses vars.n_gate_v (cells) and vars.rho_v (m/s).
Returns: float: Maximum velocity offset (m/s).

#### `T_from_Delta_r` (Function)
Calculate the time required (T) to achieve range offset (Delta_r_max) at constant acceleration (alpha).
Uses vars.Delta_r_max (m) and vars.alpha (m/s^2).
Returns: float: Time duration (s).

#### `T_from_Delta_v` (Function)
Calculate the time required (T) to achieve velocity offset (Delta_v_max) at constant acceleration (a_accel).
Uses vars.Delta_v_max (m/s) and vars.a_accel (m/s^2).
Returns: float: Time duration (s).

#### `rho_v` (Function)
Calculate Velocity Resolution (Doppler bin size).
Uses vars.wavelength (m) and vars.T_cpi (s).
Returns: float: Velocity resolution (m/s).


### Legacy RGPO
#### `delta_R_pull_from_time` (Function)
Dynamically created solver function.

#### `delta_t_pull_from_range` (Function)
Dynamically created solver function.

#### `gate_bias_error` (Function)
Dynamically created solver function.

#### `n_pulses_to_capture` (Function)
Dynamically created solver function.

#### `rgpo_capture_analysis` (Function)
Analyze when/if RGPO captures the tracking gate.

Args:
    R_true: True target range (m)
    R_jammer_start: Initial jammer range (usually equals R_true) (m)
    delta_R_pull: Range pull increment per pulse (m/pulse)
    gate_width: Range gate width (m)
    n_pulses_max: Maximum pulses to simulate

Returns:
    dict with analysis results

#### `rgpo_delay_profile` (Function)
Generate a Range Gate Pull-Off delay profile.

Args:
    R_initial: Initial range (m) where jammer matches true target
    R_final: Final range (m) to pull the gate to
    delta_R_per_pulse: Range increment per pulse (m/pulse)
    n_pulses: Number of pulses (optional, calculated if not provided)

Returns:
    dict with keys:
        - 'pulse_number': List of pulse indices
        - 'range_m': List of false target ranges (m)
        - 'delay_us': List of time delays (μs)
        - 'n_pulses': Total number of pulses
        - 'total_time_s': Total time if PRF known

#### `rgpo_max_pull_rate` (Function)
Calculate maximum safe RGPO pull rate.

The pull rate must be slow enough that the gate doesn't lose lock.
A common rule of thumb is to move no more than 10-20% of the gate
width per pulse.

Args:
    delta_R_gate: Range gate width (m)
    factor: Fraction of gate width to move per pulse (default 0.1)

Returns:
    Maximum delta_R per pulse (m/pulse)

#### `tracking_bandwidth` (Function)
Dynamically created solver function.

