#!/usr/bin/env python3
"""Test suite for Topics 12-18 (Electronic Attack techniques).

This script tests the newly integrated functionality for:
- Topic 12: Chaff
- Topic 13: Noise Jamming
- Topic 14: Gated Noise
- Topic 15: False Target Generation
- Topic 16: Radar Tracking / False Tracks
- Topic 17: Gate Stealing
- Topic 18: Cross-Eye
"""

import sys
from math import isclose

def test_topics_12_18():
    """Run tests for Topics 12-18.
    
    Returns 0 on success, 1 on failure.
    """
    try:
        import radar_range_equation as RRE

        print("Testing Topics 12-18: Electronic Attack Techniques")
        print("=" * 60)

        # =====================================================================
        # Topic 12: Chaff Tests
        # =====================================================================
        print("\n✓ Testing Topic 12: Chaff")
        RRE.vars.wavelength = 0.03  # 10 GHz
        RRE.vars.L_fiber = RRE.solve.L_fiber()
        assert isclose(RRE.vars.L_fiber, 0.015, rel_tol=1e-6), f"L_fiber calculation failed: {RRE.vars.L_fiber}"
        
        RRE.vars.D_fiber = RRE.vars.L_fiber / 1000
        RRE.vars.V_ch = RRE.solve.V_ch()
        assert RRE.vars.V_ch > 0, "V_ch must be positive"
        
        RRE.vars.V_box = 0.05 * 0.05 * 0.02
        RRE.vars.Fill_ratio = 0.60
        RRE.vars.N_fiber = RRE.solve.N_fiber()
        assert RRE.vars.N_fiber > 1e6, f"N_fiber seems low: {RRE.vars.N_fiber}"
        
        RRE.vars.zeta_ch = 0.5
        sigma_t = RRE.solve.sigma_ch_t(1.0)  # At t=1s
        assert sigma_t > 0, "Chaff RCS must be positive"
        print(f"  L_fiber={RRE.vars.L_fiber*1e3:.2f}mm, N_fiber={RRE.vars.N_fiber:.2e}")

        # =====================================================================
        # Topic 13: Noise Jamming Tests
        # =====================================================================
        print("\n✓ Testing Topic 13: Noise Jamming")
        RRE.vars.Pt = 1500
        RRE.vars.Gt = RRE.convert.db_to_lin(33)
        RRE.vars.sigma = RRE.convert.db_to_lin(13)
        RRE.vars.N_p = 1
        RRE.vars.Bj = 6e9
        RRE.vars.R = 10e3
        RRE.vars.Pj = 500
        RRE.vars.Gj = RRE.convert.db_to_lin(3)
        RRE.vars.Lossj = RRE.convert.db_to_lin(2)
        RRE.vars.B = 67e6
        
        s_j = RRE.solve.S_J_ratio()
        assert s_j > 0, "S/J ratio must be positive"
        assert s_j < 1, "S/J ratio should be less than 1 for this jamming scenario"
        
        RRE.vars.S_min = RRE.convert.db_to_lin(10)
        r_bt = RRE.solve.R_burnthrough()
        assert r_bt > 0, "Burnthrough range must be positive"
        assert r_bt < RRE.vars.R, "Burnthrough range should be less than current range"
        print(f"  S/J={s_j:.3f}, R_bt={r_bt/1e3:.2f}km")

        # =====================================================================
        # Topic 14: Gated Noise Tests
        # =====================================================================
        print("\n✓ Testing Topic 14: Gated Noise")
        RRE.vars.c = 3e8
        RRE.vars.R_tgt = 50e3
        RRE.vars.tau = 30e-6
        RRE.vars.R_gn_start_offset = 600
        
        t_2way = RRE.solve.t_tgt_2way()
        assert isclose(t_2way, 2 * 50e3 / 3e8, rel_tol=1e-6), "Two-way time calculation failed"
        
        t_gn = RRE.solve.t_gn_start_release()
        assert t_gn < t_2way, "Gated noise should start before target echo"
        print(f"  t_tgt_2way={t_2way*1e6:.2f}µs, t_gn_start={t_gn*1e6:.2f}µs")

        # =====================================================================
        # Topic 15: False Target Generation Tests
        # =====================================================================
        print("\n✓ Testing Topic 15: False Target Generation")
        RRE.vars.wavelength = 0.03
        RRE.vars.v_tgt = -120  # Closing
        RRE.vars.R = 50e3
        RRE.vars.R_ft = 60e3
        RRE.vars.v_ft = 60  # Opening
        
        RRE.vars.f_D_tgt = RRE.solve.f_D_tgt()
        assert RRE.vars.f_D_tgt > 0, "Target Doppler should be positive for closing target"
        
        RRE.vars.f_D_ft = RRE.solve.f_D_ft()
        assert RRE.vars.f_D_ft < 0, "False target Doppler should be negative for opening target"
        
        dt = RRE.solve.Delta_t_ft()
        assert dt > 0, "Time delay for false target should be positive"
        
        df = RRE.solve.Delta_f_ft()
        assert df < 0, "Frequency shift should be negative"
        print(f"  f_D_tgt={RRE.vars.f_D_tgt:.0f}Hz, f_D_ft={RRE.vars.f_D_ft:.0f}Hz")

        # =====================================================================
        # Topic 16: Radar Tracking / False Tracks Tests
        # =====================================================================
        print("\n✓ Testing Topic 16: Radar Tracking / False Tracks")
        RRE.vars.Pt = 3000
        RRE.vars.Gt = RRE.convert.db_to_lin(25)
        RRE.vars.R = 50e3
        RRE.vars.pi = 3.14159
        
        p_d = RRE.solve.P_density()
        assert p_d > 0, "Power density must be positive"
        
        RRE.vars.P_density = p_d
        RRE.vars.sigma = RRE.convert.db_to_lin(10)
        RRE.vars.Gj = RRE.convert.db_to_lin(3)
        pj_em = RRE.solve.Pj_emulated()
        assert pj_em > 0, "Jammer power must be positive"
        print(f"  P_density={p_d:.2e}W/m², Pj_emulated={pj_em*1e3:.2f}mW")

        # =====================================================================
        # Topic 17: Gate Stealing Tests
        # =====================================================================
        print("\n✓ Testing Topic 17: Gate Stealing")
        RRE.vars.wavelength = 0.025  # 12 GHz
        RRE.vars.T_cpi = 60e-3
        rho_v = RRE.solve.rho_v()
        assert rho_v > 0, "Velocity resolution must be positive"
        
        RRE.vars.n_gate_r = 10
        RRE.vars.delta_r = 15
        dr_max = RRE.solve.Delta_r_max_gate()
        assert dr_max == 150, f"Delta_r_max_gate calculation failed: {dr_max}"
        
        RRE.vars.Delta_r_max = 200
        RRE.vars.g = 9.80665
        RRE.vars.alpha = 2 * RRE.vars.g
        t_dr = RRE.solve.T_from_Delta_r()
        assert t_dr > 0, "Time to achieve offset must be positive"
        
        RRE.vars.rho_v = rho_v
        RRE.vars.n_gate_v = 9
        dv_max = RRE.solve.Delta_v_max_gate()
        assert dv_max > 0, "Max velocity offset must be positive"
        
        RRE.vars.Delta_v_max = dv_max * 2
        RRE.vars.a_accel = 1.5 * RRE.vars.g
        t_dv = RRE.solve.T_from_Delta_v()
        assert t_dv > 0, "Time to achieve velocity offset must be positive"
        print(f"  rho_v={rho_v:.3f}m/s, T_from_Delta_r={t_dr:.2f}s")

        # =====================================================================
        # Topic 18: Cross-Eye Tests
        # =====================================================================
        print("\n✓ Testing Topic 18: Cross-Eye")
        RRE.vars.L_cross = 138
        RRE.vars.R = 30e3
        RRE.vars.a_gain_ratio = -0.9
        
        phi_hat = RRE.solve.phi_hat_cross_eye_amp()
        assert phi_hat > 0, "Cross-eye angle error must be positive"
        
        RRE.vars.phi_hat_ce = 0.0035
        l_cross = RRE.solve.L_cross_from_phi_hat()
        assert l_cross > 0, "Aperture separation must be positive"
        
        # Test saturation case (a_gain_ratio = 1)
        RRE.vars.a_gain_ratio = 1.0
        phi_hat_sat = RRE.solve.phi_hat_cross_eye_amp()
        assert phi_hat_sat == float('inf'), "Should saturate to infinity when a=1"
        print(f"  phi_hat={phi_hat*1e3:.2f}mrad, L_cross={l_cross:.2f}m")

        print("\n" + "=" * 60)
        print("✓ All Topics 12-18 tests passed successfully!")
        print("=" * 60)
        return 0

    except AssertionError as e:
        print(f"\n✗ Test failed: {e}")
        return 1
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(test_topics_12_18())
