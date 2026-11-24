# P error table, errBudget_P07_100m

100M events, P_true = 0.7, 4 blocks, 25 GeV, LCRC. 16295 fit tracks.

**P = 0.470 ± 0.086 (stat) ± 0.036 (syst)**, pull −2.5 (total error 0.093).

**Total dilution:** D = D_contam = 0.780

**Total stat:** σ_D = 0.0032

**Total syst:** δ_D = D × (1 − f) × 0.415 = 0.071

**Corrected P, D_contam only = 0.603 ± 0.110 (stat) ± 0.072 (syst)**, pull −0.7 (total error 0.131).

**Corrected P, D_φ and D_contam = 0.634 ± 0.116 (stat) ± 0.075 (syst)**, pull −0.5 (total error 0.138).



## Errors

| Type | Group | Item | On P | Gaussian? | In the code's ±? |
|---|---|---|---|---|---|
| Stat | σ_P | errFit (likelihood) | 0.0858 | Yes | Yes |
| Stat | σ_AP | AP spread over tracks | 0.0023 | Yes | Yes |
| Syst | δ_AP | model coefficients A, k | 0.0202 | No (envelope, no covariance) | Yes |
| Syst | δ_AP | beam momentum spread (conservative) | 0.0167 | Yes | Yes |
| Stat | σ_θ | random statistical error from the finite number of tracks at each est neutron θ | 0.0009 | Yes | Yes (macro edited, not rerun) |
| Syst | δ_θ | est neutron θ range | 0.0015 | No (half-width of a flat range) | Yes |
| Stat | σ_φ | random statistical error from the finite number of tracks at each est neutron φ | 0.0010 | Yes | Yes (macro edited, not rerun) |
| Syst | δ_φ | est neutron φ range (the full size of the φ dilution correction, P × (1/D − 1), D = 0.950) | 0.0245 | No (envelope, full size of the correction) | Yes (macro edited, not rerun) |
| Stat | | **σ_stat** | **0.0859** | Yes | Yes |
| Syst | | **δ_syst** | **0.036** | No | Partly |

## Dilution (sim only: each needs the true θ, φ or CE label)
systemic

| Type | Group | Item | On P | Applied? |
|---|---|---|---|---|
| Dilution | D_contam | non-CE tracks flat in φ: CE 12706, background 3589, of 16295 fit tracks | CE fraction 0.780 | No |
| Stat | σ_D,contam | binomial error on the CE fraction, ±0.0032; on the counts ±52 CE and ±52 background (total fixed at 16295) | 0.0026 | No |
| Syst | δ_D,contam | uncertainty on the CE fraction and so on the CE count and background count: CE cross section from the errors on A (1828.5 ± 27.7), α, g0, g1, c0, c1, no correlations, 41% at 25 GeV. P/f × (1 − f) × 0.415 | 0.058 | No |
| Dilution | D_θ | est n frame block θ − true n frame θ (smearBias −0.00405 ± 0.00023) (the θ version of D_φ, but an additive shift on A_eff, not a factor) | +0.0151 | No (the code applies it, taken out here) |
| Stat | σ_D,θ | error on the mean of est n frame θ − true n frame θ (smearBias), the θ version of σ_D,φ | 0.0009 | No |
| Syst | δ_D,θ | model dependence of smearBias (sim physics vs real data), at most the size of the correction | not computed, at most 0.0151 | No |
| Dilution | D_φ | est n frame φ − true n frame φ, D = 0.950 ± 0.002 (a factor on P, ÷ D) | +0.0245 | No |
| Stat | σ_D,φ | error on the D mean, the φ version of σ_D,θ | 0.0010 | No |
| Syst | δ_D,φ | model dependence of D (sim physics vs real data), at most the size of the correction | not computed, at most 0.0245 | No |
| Dilution | D_AP | t mismatch, A_model / A_analysis = 0.9966 | +0.0016 | No |
| Stat | σ_D,AP | error on the ratio from the per-track A_model − A_analysis (100M run errBudget_P0_100m: 0.00037 on the ratio, times P 0.470) | 0.00017 | No |
| Syst | δ_D,AP | model dependence of the t mismatch, at most the size of the correction | not computed, at most 0.0016 | No |
| Dilution | D_t | t = (p_beam sin θ)², shift on A_eff from the p_beam and θ errors, already in the errors table above | see δ_AP beam momentum spread 0.0167 and δ_θ 0.0015 | In the code's ± |
| Stat | σ_D,t | error from the p_beam error and the θ error propagated through t, already in the errors table above | 0.0167 (p_beam), 0.0015 (θ) | In the code's ± |
| Syst | δ_D,t | model dependence of the t shift (p_beam and θ errors) | not computed | No |
| Dilution | D_track | number of recon tracks whose hits all share the same trackID / number of recon tracks | not computed | No |
| Stat | σ_D,track | random statistical error from the finite number of tracks | not computed | No |
| Syst | δ_D,track | model dependence of the track reconstruction (sim physics vs real data) | not computed | No |

## Fits

| Fit | P | Stat error | (P − 0.7) / σ_stat |
|---|---|---|---|
| Total, real data only (the P above) | 0.470 | 0.086 | −2.7 |
| Total, also ÷ D_φ 0.950 and ÷ CE fraction 0.780 (sim only) | 0.634 | 0.116 | −0.6 |
| CE only tracks, est angle (12706, sim only) | 0.721 | 0.125 | +0.2 |
| CE only tracks, true angle (12706, sim only) | 0.718 | 0.125 | +0.1 |

The CE only rows still include the θ smear correction. The total reads low because 22% of its tracks are non-CE and flat in φ.

## Closure

Closure is the check that the corrected P matches the P put into the sim: P − P_true, against the total error.

| Sample | P − P_true | Total error | Pull |
|---|---|---|---|
| Total, real data only | −0.230 | 0.093 | −2.5 |
| Total, also ÷ D_φ and ÷ CE fraction (sim only) | −0.066 | 0.126 | −0.5 |

It closes only after the purity and φ dilution corrections, so the −2.5σ miss is that dilution, not a fit problem.

## How each value is calculated

### Errors (the table above)

Here A_eff = −0.1298 is the mean of AP(θ_i) over the 16295 fit tracks, with AP(θ) = −A (1 − e^(−k p)) √t, t = (p sin θ)², p the nominal beam momentum and θ the est n frame block θ.

- **errFit (σ_P):** unbinned likelihood in φ, L = Π (1 + P A_eff cos φ_i) / norm. P comes from Newton's method started at 0 and limited to [−1, 1]. errFit = 1 / √(−d²lnL/dP²) at the maximum, about √(2/N) / |A_eff| = 0.0854 (fitted 0.0886).
- **σ_AP:** StatErr(A_eff) = std(AP_i) / √N = 0.00064, the spread of AP over the fit tracks. On P: |P| × StatErr / |A_eff| = 0.0023.
- **σ_θ, σ_φ:** the statistics of the smearBias and D means (σ_D,θ and σ_D,φ below), 0.0009 and 0.0010 on P. They are added to σ_stat in `PErrorBudget.cpp` (edited, not rerun).
- **δ_model:** A_eff is recomputed at (A + σ_A, k + σ_k) and at (A − σ_A, k − σ_k), with A = 0.7079 ± 0.0303 and k = 0.3860 ± 0.0523. The relative error is half the difference over A_eff (kModelAFrac = 0.043), and on P it is |P| × 0.043 = 0.0202.
- **δ_beam:** A_eff is recomputed at p + σ_p and p − σ_p. The relative error is half the difference over A_eff (kBeamAFrac = 0.0355), and on P it is |P| × 0.0355 = 0.0167.
- **δ_θ:** |P| × half the A_eff difference between the two ends of each track's est neutron θ range, 0.0015.
- **δ_φ:** the full size of the φ dilution correction, P × (1/D − 1) = 0.0245 (D = 0.950). The D difference between the ends of the φ range alone gave about 0, so the correction itself is taken as the uncertainty. Added in `PErrorBudget.cpp` as delta_phi_range (edited, not rerun).
- **How they join in the code's ±:** rel = √(rel_spread² + rel_beam² + rel_model² + rel_θrange²), the A part is |P| × rel / (1 − rel²), and the quoted error is √(errFit² + A part²). The code's θ range term is the est neutron θ range (0.0015 on P).
- **In this table:** σ_stat = √(errFit² + σ_AP² + σ_θ² + σ_φ²) = 0.0859, δ_syst = √(δ_model² + δ_beam² + δ_θ,est² + δ_φ²) = 0.036, and the total is √(σ_stat² + δ_syst²) = 0.093. Each range term is |P| × half the A_eff (or D) difference between the two ends of the track's range.

### Dilution (the table above)

All of these use the sim truth, so none exist on real data.

- **D_θ:** smearBias = mean AP(est n frame block θ) − mean AP(true n frame θ) over the CE recon tracks, −0.00405. The correction is P = P_raw × A_eff / (A_eff − smearBias), a shift of +0.0151.
- **σ_D,θ:** std of the per-track AP difference / √N = 0.00023 on A_eff. On P: |P| × 0.00023 / |A_eff| = 0.0009.
- **D_φ:** D = mean cos(φ_est − φ_true) over the CE recon tracks in the true neutron frame, 0.950 ± std/√N = 0.002. The shift on P is P × (1/D − 1) = +0.0245.
- **σ_D,φ:** |P| × 0.002 / 0.950 = 0.0010.
- **D_contam:** f = CE only tracks / all fit tracks = 12706 / 16295 = 0.780. If the non-CE tracks are flat in φ, P_true ≈ P / f.
- **σ_D,contam:** √(f (1 − f) / N) = 0.0032, so ±52 on the CE and background counts (N × 0.0032, total fixed at 16295). On the purity-corrected P (0.634): 0.634 × 0.0032 / 0.780 = 0.0026.
- **δ_D,contam:** the CE cross section of `G4ChargeExchangeNP.cc`, σ ∝ A (s/s0)^(2α−2) (1 + G) / C with G = 1 + g0 + g1 L, C = 1 + c0 + c1 L, L = ln(s/s0), s0 = 10 GeV², s = 50.4 GeV² at 25 GeV. The relative error is the log-derivative times each parameter error (A 1828.5 ± 27.7, α, g0, g1, c0, c1), added in quadrature with no correlations, 0.415. A higher cross section raises f, so on P it is P / f × (1 − f) × 0.415 = 0.058 (0.055 for P / D_contam alone).
- **D_AP (t mismatch):** A_analysis = mean AP at the true θ with t = (p sin θ)², A_model = mean LastAP, the AP the model used for each CE. A_model / A_analysis = 0.9966, and the size on P is P × (1 − ratio) = 0.0016.
- **σ_D,AP:** std of the per-track A_model − A_analysis / √N, over the mean A_analysis. The 100M P = 0 run (errBudget_P0_100m) gives 0.00037 on the ratio, so on P 0.470 × 0.00037 = 0.00017.
- **D_t:** t = (p_beam sin θ)². Its p_beam and θ errors are the δ_beam (0.0167) and δ_θ (0.0015) rows above, so σ_D,t is those values.
- **D_track:** recon tracks whose hits all share one trackID / all recon tracks, counted in `F8_TrackerFits.cpp` (pure and mixed tracks, D_track.txt). σ_D,track = √(D (1 − D) / N). Job 7871605 on errBudget_P0_100m is running, so not computed yet.
- **δ_D,θ, δ_D,φ, δ_D,AP:** not computed, at most the size of each correction (0.0151, 0.0245, 0.0016). A real number needs sims with other physics.
- **δ_D,t, δ_D,track:** not computed, they need sims with other physics.
- **Total dilution:** D = D_contam = 0.780, σ_D = 0.0032, δ_D = D × (1 − f) × 0.415 = 0.071. The corrected P is P / D, with its stat and syst the real data ones divided by D plus δ_D,contam in quadrature on the syst.


## Notes

- The θ and φ range rows are the real-data errors on the angles. The measured est n frame − true n frame differences (D_θ, D_φ) hold the same effect in the sim, but need truth, so they are not added on top.
- The one-hit-per-plane and other cut rows are left out. The raw run summary counted them, so its δ_syst (0.479) is not this table's.
- Not computed: efficiency, false asymmetry, seed road, plane weights, kTubeSigmas, H/C target mix.
