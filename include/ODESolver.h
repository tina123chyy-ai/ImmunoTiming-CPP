#pragma once
#include "TreatmentScheduler.h"
#include <vector>

namespace immuno {

enum class Method { Euler, RK4 };
struct Point { double time; State state; };
struct SimulationOptions {
    double end_time{40.0};
    double step{0.05};
    double output_interval{0.5};
    Method method{Method::RK4};
};
struct SimulationResult { std::vector<Point> points; };

SimulationResult simulate(const Parameters& params, const State& initial,
                          const std::vector<Dose>& doses, const SimulationOptions& options);
} // namespace immuno
