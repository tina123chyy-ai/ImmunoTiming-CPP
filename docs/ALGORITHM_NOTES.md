# C++ and data structures notes / 算法说明

## 1. ODE integration

The state vector is `x=[T,C,E]`, and the derivative function computes `f(x,u)`.

- Euler: `x_next = x + h*f(x,u)` (global error O(h))
- RK4: four slopes `k1,k2,k3,k4`; `x_next = x + h*(k1+2k2+2k3+k4)/6` (global error O(h^4) for sufficiently smooth segments)

Because vaccine/CAR-T inputs use step functions, the solver integrates **up to** each start/end event and then updates the rate. It does not apply a single RK4 step across an input discontinuity.

## 2. Event scheduler

Each dose becomes two events: `(start,+rate)` and `(end,-rate)`. A `std::priority_queue` ordered by smallest event time manages events. Adding 2m events for m doses costs O(m log m); processing is O(m log m). The dose list represents synthetic model inputs, NOT real injection regimens.

## 3. Top-K scenarios

For each pair `(vaccine_day, cart_day)` in the grid, run the model and compute

`score = trapezoidal_AUC(T) + 0.5 * T_final`.

Maintain a max heap of at most K candidates; when a better candidate is observed, evict the worst. For N candidates, selection bookkeeping is O(N log K). **Numerical simulation is additional** and usually dominates cost.

## 4. Try these exercises yourself

1. Set `beta=0` and explain the meaning of the removed feedback.
2. Set `gamma=0` and compare to the previous experiment.
3. Add a second vaccine pulse and verify that the scheduler handles overlapping dose periods.
4. Reduce the step size from 0.1 to 0.05 to 0.025 and measure convergence.
5. Make a `--mode parameter_sweep` that stores all candidates, not only Top-K.
6. Add a unit test validating two simultaneous events at the same time.
7. Implement a small sensitivity comparison for different carrying capacities `K`.

## 5. Limitations of the ranking objective

Optimizing only tumor AUC creates a particular preference for early control. A real decision would require far more information (toxicity, safety, pharmacology, uncertainty and patient-specific clinical constraints). This project is not appropriate for real treatment planning.
