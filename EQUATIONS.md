## Equations Reference

This section contains all equations available in the `radar_range_equation` package.


### Base/Common

**A_e**:

```math
A_{e} = D_{h} D_{v} \eta
```

**G_t**:

```math
G_{t} = \frac{4 \pi A_{e}}{\lambda^{2}}
```

**Linear_to_dB**:

```math
x = \frac{10 \log{\left(x \right)}}{\log{\left(10 \right)}}
```

**P_t**:

```math
P_{t} = \frac{64 R^{4} S_{min} \pi^{3}}{G_{t}^{2} \lambda^{2} \sigma}
```

**R4**:

```math
R^{4} = \frac{G_{t}^{2} P_{t} \lambda^{2} \sigma}{64 \pi^{3} S_{min}}
```

**R_max**:

```math
R_{max} = \sqrt[4]{\frac{G_{t}^{2} P_{t} \lambda^{2} \sigma}{64 S_{min} \pi^{3}}}
```

**theta_B**:

```math
\theta_{B} = \frac{13 \pi \lambda}{36 D_{h}}
```

**wavelength**:

```math
\lambda = \frac{c}{f}
```


### Doppler CW Radar

**eq_delta_v**:

```math
\delta_{v} = \frac{\lambda}{2 T_{cpi}}
```

**eq_f_doppler**:

```math
f_{doppler} = - \frac{2 v}{\lambda}
```

**eq_f_obs_if**:

```math
f_{obs} = f_{doppler} + f_{if}
```

**eq_v_from_doppler**:

```math
v = - \frac{f_{doppler} \lambda}{2}
```


### CWFM Radar

**R_cwfm**:

```math
R = \frac{c f_{r}}{4 Delta f f_{m}}
```

**f_0_cwfm**:

```math
f_{0} = 2 Delta f f_{m}
```

**f_d_cwfm**:

```math
f_{d} = - 0.5 f_{bd} + 0.5 f_{bu}
```

**f_m_cwfm**:

```math
f_{m} = \frac{c}{2 R_{un}}
```

**f_r_cwfm**:

```math
f_{r} = 0.5 f_{bd} + 0.5 f_{bu}
```

**v_cwfm**:

```math
v = - \frac{c f_{d}}{2 f}
```


### Pulsed Radar

**eq_E_i_sqrt_n**:

```math
E_{i} = \frac{1}{\sqrt{n_{p}}}
```

**eq_R_from_time**:

```math
R = \frac{c t_{delay}}{2}
```

**eq_R_un_from_fp**:

```math
R_{un} = \frac{c}{2 f_{p}}
```

**eq_S_N_n_coherent**:

```math
S_{N n} = S_{N 1} n_{p}
```

**eq_S_N_n_noncoherent**:

```math
S_{N n} = S_{N 1} \sqrt{n_{p}}
```

**eq_S_N_n_noncoherent_Ei**:

```math
S_{N n} = E_{i} S_{N 1} n_{p}
```

**eq_T_p**:

```math
T_{p} = \frac{1}{f_{p}}
```

**eq_fp_from_R_un**:

```math
f_{p} = \frac{c}{2 R_{un}}
```

**eq_n_p**:

```math
n_{p} = \frac{f_{p} t_{scan} \theta_{B}}{2 \pi}
```

**eq_n_pulses_rgpo**:

```math
n_{pulses capture} = \frac{R_{false} - R_{true}}{\delta_{R pull}}
```

**eq_tau_from_duty**:

```math
\tau = T_{p} duty_{cycle}
```

**eq_theta_B_skolnik**:

```math
\theta_{B} = \frac{1.2 \lambda}{D_{h}}
```


### Direction Finding

**S_N_from_dB**:

```math
S_{N} = 10^{\frac{S_{N dB}}{10}}
```

**Theta**:

```math
\Theta = \frac{4 \log{\left(2 \right)}}{\theta_{B}^{2}}
```

**eq_phi_hat_ce_amp**:

```math
\phi_{hat ce} = \frac{L \left(a_{ratio} + 1\right)}{2 R \left(1 - a_{ratio}\right)}
```

**phi_hat_amp**:

```math
\phi_{hat} = \frac{\Delta \theta_{B}^{2}}{8 \Sigma \phi_{s} \log{\left(2 \right)}}
```

**sigma_phi_amp**:

```math
\sigma_{\phi} = \frac{\sqrt{2} \theta_{B}^{2} \sqrt{\frac{1}{S_{N}}}}{16 \phi_{s} \log{\left(2 \right)}}
```

**sigma_phi_phase**:

```math
\sigma_{\phi} = \frac{\lambda \sqrt{\frac{1}{S_{N}}}}{2 d \pi}
```

**sigma_phi_time**:

```math
\sigma_{\phi} = \frac{c}{B d}
```

**v_phi**:

```math
v_{\phi} = e^{- \Theta \left(\phi - \phi_{s}\right)^{2}}
```


### Pulse Compression

**eq_B_chirp**:

```math
B = \gamma \tau
```

**eq_PCR_1**:

```math
PCR = B \tau
```

**eq_PCR_2**:

```math
PCR = \gamma \tau^{2}
```

**eq_delta_r_compressed**:

```math
\delta_{r} = \frac{c}{2 B}
```

**eq_delta_r_uncompressed**:

```math
\delta_{r} = \frac{c \tau}{2}
```

**eq_f_range_tone**:

```math
f_{range tone} = - \frac{2 R_{offset} \gamma}{c}
```


### Chaff

**eq_L_fiber**:

```math
L_{fiber} = \frac{\lambda}{2}
```

**eq_N_fiber**:

```math
N_{fiber} = \frac{Fill_{ratio} V_{box}}{V_{ch}}
```

**eq_V_ch**:

```math
V_{ch} = \frac{\pi D_{fiber}^{2} L_{fiber}}{4}
```

**eq_sigma_ch_max**:

```math
\sigma_{ch} = 0.15 N_{fiber} \lambda^{2}
```

**eq_sigma_ch_t**:

```math
\sigma_{ch} = 0.15 N_{fiber} \lambda^{2} \left(1 - e^{- \frac{t_{delay}}{\zeta_{ch}}}\right)
```


### Noise Jamming

**eq_R_bt**:

```math
R_{bt}^{2} = \frac{Bj G_{t} P_{t} n_{p} \sigma}{4 \pi B Gj Lossj Pj S_{min}}
```

**eq_S_J_ratio**:

```math
S/J = \frac{Bj G_{t} P_{t} n_{p} \sigma}{4 \pi B Gj Lossj Pj R^{2}}
```


### Gated Noise

**eq_t_gn_start_release**:

```math
t_{gn start release} = - \frac{\tau}{2} + \frac{- 2 R_{gn start offset} + 2 R_{tgt}}{c}
```

**eq_t_tgt_2way**:

```math
t_{tgt 2way} = \frac{2 R_{tgt}}{c}
```


### False Target Generation

**eq_Delta_f_ft**:

```math
\Delta_{f ft} = f_{D ft} - f_{D tgt}
```

**eq_Delta_t_ft**:

```math
\Delta_{t ft} = - \frac{2 R}{c} + \frac{2 R_{ft}}{c}
```

**eq_f_D_ft**:

```math
f_{D ft} = - \frac{2 v_{ft}}{\lambda}
```

**eq_f_D_tgt**:

```math
f_{D tgt} = - \frac{2 v_{tgt}}{\lambda}
```


### Radar Tracking

**eq_P_density**:

```math
P_{density} = \frac{G_{t} P_{t}}{4 \pi R^{2}}
```

**eq_Pj_emulated**:

```math
Pj_{emulated} = \frac{P_{density} \sigma}{Gj}
```


### Gate Stealing

**eq_Delta_r_max_gate**:

```math
\Delta_{r max} = \delta_{r} n_{gate r}
```

**eq_Delta_r_t**:

```math
\Delta_{r(t)} = \frac{\alpha t_{delay}^{2}}{2}
```

**eq_Delta_v_max_gate**:

```math
\Delta_{v max} = n_{gate v} \rho_{v}
```

**eq_Delta_v_t**:

```math
\Delta_{v(t)} = a t_{delay}
```

**eq_T_from_Delta_r**:

```math
\Delta_{r max} = \frac{T^{2} \alpha}{2}
```

**eq_T_from_Delta_v**:

```math
\Delta_{v max} = T a
```

**eq_rho_v**:

```math
\rho_{v} = \frac{\lambda}{2 T_{cpi}}
```


### Cross-Eye

**eq_J_ratio**:

```math
a_{ratio} = \frac{J_{1}}{J_{2}}
```


### Legacy RGPO

**eq_BW_track**:

```math
BW_{track} = \frac{\alpha_{track} f_{p}}{2 \pi}
```

**eq_delta_R_from_time**:

```math
\delta_{R pull} = \frac{c \delta_{t pull}}{2}
```

**eq_delta_t_pull**:

```math
\delta_{t pull} = \frac{2 \delta_{R pull}}{c}
```

**eq_gate_bias**:

```math
gate_{bias} = R_{gate} - R_{true}
```

**eq_pull_rate**:

```math
pull_{rate} = \delta_{R pull}
```
