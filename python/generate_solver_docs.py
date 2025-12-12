#!/usr/bin/env python3
"""Generate documentation for all solvers.

This script extracts all solvers from the solve module and formats them
for documentation in SOLVERS.md.
"""

import sys
import inspect


def generate_solver_docs():
    """Generate documentation for all solvers."""
    try:
        from radar_range_equation import solve, Solvable
        
        print("=" * 70)
        print("Radar Range Equation - Solver Reference")
        print("=" * 70)
        print()
        
        # Get all solver attributes from solve class
        solver_names = [name for name in dir(solve) 
                       if not name.startswith('_')]
        
        # Categorize solvers
        solver_info = []
        for name in solver_names:
            attr = getattr(solve, name)
            solver_type = None
            docstring = None
            
            if isinstance(attr, Solvable):
                solver_type = "Solvable"
                # Get docstring from the Solvable if available
                docstring = getattr(attr, '__doc__', None)
            elif callable(attr):
                solver_type = "Function"
                docstring = inspect.getdoc(attr)
            
            if solver_type:
                solver_info.append({
                    'name': name,
                    'type': solver_type,
                    'doc': docstring or "No documentation available."
                })
        
        # Sort by name
        solver_info.sort(key=lambda x: x['name'])
        
        # Group solvers by topic
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
        
        for solver in solver_info:
            name = solver['name']
            
            # Categorize solvers based on naming
            if 'doppler' in name or 'delta_v' == name or 'v_from_doppler' in name or 'f_obs' in name:
                topics["Doppler CW Radar"].append(solver)
            elif '_cwfm' in name or name in ['R_cwfm', 'v_cwfm', 'f_0_cwfm', 'f_d_cwfm', 'f_m_cwfm', 'f_r_cwfm']:
                topics["CWFM Radar"].append(solver)
            elif 'R_un' in name or 'fp_from' in name or 'T_p' in name or \
                 'S_N_n' in name or 'tau_from' in name or 'R_from_time' in name or 'calculate_Theta' in name:
                topics["Pulsed Radar"].append(solver)
            elif 'phi' in name or 'sigma_phi' in name or name == 'S_N_from_dB' or \
                 name == 'v_phi' or 'd_from_sigma' in name or 'B_from_sigma' in name or \
                 'estimate_phi' in name or 'phi_s_from' in name:
                topics["Direction Finding"].append(solver)
            elif 'delta_r' in name or 'PCR' in name or 'B_chirp' in name or \
                 'f_range_tone' in name or 'R_offset' in name:
                topics["Pulse Compression"].append(solver)
            elif 'L_fiber' in name or 'V_ch' in name or 'N_fiber' in name or 'sigma_ch' in name:
                topics["Chaff"].append(solver)
            elif 'S_J' in name or 'R_burnthrough' in name or name == 'R_burnthrough':
                topics["Noise Jamming"].append(solver)
            elif 't_gn' in name or 't_tgt' in name:
                topics["Gated Noise"].append(solver)
            elif 'f_D' in name or 'Delta_t_ft' in name or 'Delta_f_ft' in name:
                topics["False Target Generation"].append(solver)
            elif 'P_density' in name or 'Pj_emulated' in name:
                topics["Radar Tracking"].append(solver)
            elif 'rho_v' in name or 'Delta_r_max' in name or 'Delta_v_max' in name or \
                 'T_from_Delta' in name:
                topics["Gate Stealing"].append(solver)
            elif 'cross_eye' in name or 'L_cross' in name:
                topics["Cross-Eye"].append(solver)
            elif 'rgpo' in name or 'delta_t_pull' in name or 'delta_R_pull' in name or \
                 'pull_rate' in name or 'n_pulses_to_capture' in name or \
                 'gate_bias' in name or 'tracking_bandwidth' in name:
                topics["Legacy RGPO"].append(solver)
            else:
                topics["Base/Common"].append(solver)
        
        # Generate markdown output
        output = []
        output.append("## Solver Reference\n")
        output.append("This section contains all solver functions and Solvable instances available in the `radar_range_equation.solve` module.\n")
        output.append(f"\n**Total Solvers**: {len(solver_info)}\n")
        
        for topic, solvers in topics.items():
            if not solvers:
                continue
                
            output.append(f"\n### {topic}\n")
            
            for solver in solvers:
                output.append(f"#### `{solver['name']}` ({solver['type']})\n")
                
                # Format docstring
                doc_lines = solver['doc'].split('\n')
                for line in doc_lines:
                    if line.strip():
                        output.append(f"{line}\n")
                    else:
                        output.append("\n")
                
                output.append("\n")
        
        # Write to file
        with open('SOLVERS.md', 'w') as f:
            f.write(''.join(output))
        
        print(f"✓ Generated documentation for {len(solver_info)} solvers")
        print(f"✓ Documentation written to SOLVERS.md")
        
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
    sys.exit(generate_solver_docs())
