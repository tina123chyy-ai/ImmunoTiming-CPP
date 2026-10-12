#pragma once
#include "ODESolver.h"
#include <vector>

namespace immuno {
struct Metrics {
    double tumor_auc{0.0};
    double final_tumor{0.0};
    double peak_cart{0.0};
    double peak_effector{0.0};
    double objective{0.0};
};
struct Candidate {
    double vaccine_day{0.0};
    double cart_day{0.0};
    Metrics metrics;
};
Metrics compute_metrics(const SimulationResult& result);
std::vector<Candidate> top_k_schedules(const Parameters& p, const State& init,
    const SimulationOptions& opt, const std::vector<double>& vaccine_days,
    const std::vector<double>& cart_days, int k);
} // namespace immuno
