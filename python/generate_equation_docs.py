#!/usr/bin/env python3
"""Generate LaTeX documentation for all equations.

This script extracts all equations from the equations module and formats them
as LaTeX for documentation in README.md.
"""

import sys
import sympy


def generate_equation_docs():
    """Generate LaTeX documentation for all equations."""
    try:
        from radar_range_equation import equations
        
        print("=" * 70)
        print("Radar Range Equation - Equation Reference")
        print("=" * 70)
        print()
        
        # Get all equation attributes from equations class
        equation_names = [name for name in dir(equations) 
                         if not name.startswith('_') and not name.endswith('_sym')]
        
        # Filter to only actual equations (sympy.Eq objects)
        equation_attrs = []
        for name in equation_names:
            attr = getattr(equations, name)
            if isinstance(attr, sympy.Eq):
                equation_attrs.append((name, attr))
        
        # Sort by name
        equation_attrs.sort(key=lambda x: x[0])
        
        # Group equations by topic
        topics = {
            "Base/Common": [],
            "Doppler CW Radar": [],
            "CWFM Radar": [],
            "Pulsed Radar": [],
            "Direction Finding": [],
            "Pulse Compression": [],
            "Chaff": [],
            "Noise Jamming": [],
            "Gated Noise": [],
            "False Target Generation": [],
            "Radar Tracking": [],
            "Gate Stealing": [],
            "Cross-Eye": [],
            "Legacy RGPO": []
        }
        
        for name, eq in equation_attrs:
            # Categorize equations based on naming
            if name.startswith('eq_f_doppler') or name.startswith('eq_v_from_doppler') or \
               name.startswith('eq_delta_v') or name.startswith('eq_f_obs'):
                topics["Doppler CW Radar"].append((name, eq))
            elif '_cwfm' in name or name in ['R_cwfm', 'v_cwfm', 'f_0_cwfm', 'f_d_cwfm', 'f_m_cwfm', 'f_r_cwfm']:
                topics["CWFM Radar"].append((name, eq))
            elif name.startswith('eq_R_un') or name.startswith('eq_fp') or name.startswith('eq_T_p') or \
                 name.startswith('eq_S_N_n') or name.startswith('eq_n_p') or name.startswith('eq_tau') or \
                 name.startswith('eq_R_from_time') or name.startswith('eq_E_i') or name == 'eq_theta_B_skolnik':
                topics["Pulsed Radar"].append((name, eq))
            elif 'phi' in name or name.startswith('sigma_phi') or name == 'S_N_from_dB' or \
                 name == 'v_phi' or name == 'Theta':
                topics["Direction Finding"].append((name, eq))
            elif name.startswith('eq_delta_r') or name.startswith('eq_PCR') or \
                 name.startswith('eq_B_chirp') or name.startswith('eq_f_range_tone'):
                topics["Pulse Compression"].append((name, eq))
            elif name.startswith('eq_L_fiber') or name.startswith('eq_V_ch') or \
                 name.startswith('eq_N_fiber') or name.startswith('eq_sigma_ch'):
                topics["Chaff"].append((name, eq))
            elif name.startswith('eq_S_J') or name.startswith('eq_R_bt'):
                topics["Noise Jamming"].append((name, eq))
            elif name.startswith('eq_t_gn') or name.startswith('eq_t_tgt'):
                topics["Gated Noise"].append((name, eq))
            elif name.startswith('eq_f_D') or name.startswith('eq_Delta_t_ft') or name.startswith('eq_Delta_f_ft'):
                topics["False Target Generation"].append((name, eq))
            elif name.startswith('eq_P_density') or name.startswith('eq_Pj_emulated'):
                topics["Radar Tracking"].append((name, eq))
            elif name.startswith('eq_rho_v') or name.startswith('eq_Delta_r_max') or \
                 name.startswith('eq_Delta_v_max') or name.startswith('eq_T_from') or \
                 name.startswith('eq_Delta_r_t') or name.startswith('eq_Delta_v_t'):
                topics["Gate Stealing"].append((name, eq))
            elif name.startswith('eq_J_ratio') or name.startswith('eq_phi_hat_ce'):
                topics["Cross-Eye"].append((name, eq))
            elif name.startswith('eq_delta_t_pull') or name.startswith('eq_delta_R') or \
                 name.startswith('eq_pull_rate') or name.startswith('eq_n_pulses_rgpo') or \
                 name.startswith('eq_gate_bias') or name.startswith('eq_BW_track'):
                topics["Legacy RGPO"].append((name, eq))
            else:
                topics["Base/Common"].append((name, eq))
        
        # Generate markdown output
        output = []
        output.append("## Equations Reference\n")
        output.append("This section contains all equations available in the `radar_range_equation` package.\n")
        
        for topic, equations_list in topics.items():
            if not equations_list:
                continue
                
            output.append(f"\n### {topic}\n")
            
            for name, eq in equations_list:
                # Convert to LaTeX
                latex = sympy.latex(eq)
                
                # Format the equation in markdown with LaTeX
                output.append(f"**{name}**:\n")
                output.append(f"```math\n{latex}\n```\n")
        
        # Write to file
        with open('EQUATIONS.md', 'w') as f:
            f.write('\n'.join(output))
        
        print(f"✓ Generated documentation for {len(equation_attrs)} equations")
        print(f"✓ Documentation written to EQUATIONS.md")
        
        # Also print to stdout for README integration
        print("\n" + "=" * 70)
        print("Markdown output (for README.md):")
        print("=" * 70)
        print('\n'.join(output))
        
        return 0
        
    except ImportError as e:
        print(f"✗ ImportError: {e}")
        return 1
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(generate_equation_docs())
