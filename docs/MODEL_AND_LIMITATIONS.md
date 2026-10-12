# Mathematical model, implementation choices and limitations

## What is taken from the underlying review

The manuscript *Personalized Neoantigen Vaccine Combined with CAR-T Therapy for Solid Tumor Immunotherapy*, section 4.3, describes 3 state variables and the conceptual interactions:

- `T(t)` = normalized tumor burden
- `C(t)` = CAR-T population
- `E(t)` = vaccine-driven endogenous effector T-cell population

The implemented equations match the structure stated in the paper:

\[
\dot T = rT(1-T/K)-k_C CT-k_E ET
\]
\[
\dot C = \alpha_C(1+\gamma E) \frac{T}{h_C+T} C-\delta_C C+u_C(t)
\]
\[
\dot E = \alpha_V V(t)+\beta CT-\delta_E E
\]

`beta` is a hypothesis relating CAR-T-mediated tumor killing/antigen release to additional effector activity; `gamma` assumes vaccine-mediated immunity increases the CAR-T growth term. **Both encode assumed beneficial interactions**, so a lower tumor burden in the simulated combined arm does not independently establish real synergy.

## What this code adds (post-paper computational implementation)

- Explicit normalized demonstration parameter values in `include/TumorModel.h` (not fitted)
- Discontinuous, piecewise-constant input functions approximated by rate pulses: `V(t)` of duration 2 model-time units, `u_C(t)` of duration 1
- RK4 and Euler solvers with integration steps ending exactly at each event time, so an RK4 step never crosses a discontinuity in input rates
- A fixed-horizon trapezoidal integral of tumor burden (AUC)
- A deliberately arbitrary illustrative ranking objective `AUC + 0.5 * final_tumor`
- Discrete brute-force candidate scheduling grid with a bounded Top-K max heap
- Sensitivity examples for beta and gamma

There is **no delay-differential equation** for delayed vaccine response. An input's modeled effect starts during the input pulse via `alphaV V(t)`; biological immune priming delay is not explicitly represented. Therefore, this model should *not* claim to identify biologically meaningful ideal injection intervals.

## Numerical reproducibility

- Baseline parameters and initial conditions are written in code, not estimated from data.
- The main run uses RK4, time step 0.05, 40 model-time units, sampled every 0.5 units.
- `ctest` checks event boundaries, positivity (for tested scenarios), ranking, and the solver against the analytic logistic-growth solution when both treatment populations are zero.
- `--method euler` provides a simple numerical-method comparison.

## Limitations (do not omit in reports or applications)

1. The 3-variable equations simplify complex immune and tumor processes. They omit spatial gradients, mutation/antigen heterogeneity, detailed pharmacokinetics, toxicity, host variability and explicit time-to-prime immune responses.
2. Parameter values are normalized, hypothetical, dimensionally illustrative and **not biologically calibrated**. "Day" labels are just model time units, not clinically meaningful days.
3. The hypothesized positive feedback terms may bias the model toward combined-treatment benefit; sensitivity runs explore dependence on those assumptions but do not validate them.
4. Top-K optimization is with respect to an arbitrary mathematical score only and **does not optimize actual patient outcomes**.
5. No experiment, data fitting, patient study, evidence of clinical safety, or clinical prediction was performed.
6. The paper's manually constructed Figure 2 is not a numerical solution; calculated trajectories in this project may differ, and differences should be discussed rather than concealed.

## How to write about the project accurately

Good: "I used C++ to numerically implement a three-state conceptual immunotherapy model and explored how schedule choices and feedback assumptions changed its trajectories."

Not good: "I developed an AI tool that predicts optimal CAR-T administration schedules for cancer patients."

Only say "I implemented" if you personally did the implementation or substantial independently verifiable work. With external assistance, describe the assistance accurately.
