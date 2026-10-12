#include "ScheduleRanker.h"
#include <algorithm>
#include <queue>
#include <stdexcept>

namespace immuno {
Metrics compute_metrics(const SimulationResult& res) {
    if (res.points.empty()) throw std::invalid_argument("Empty trajectory");
    Metrics m;
    for (size_t i=0;i<res.points.size();++i) {
        const auto& curr=res.points[i];
        m.peak_cart=std::max(m.peak_cart,curr.state.cart);
        m.peak_effector=std::max(m.peak_effector,curr.state.effector);
        if (i>0) {
            const auto& prev=res.points[i-1];
            m.tumor_auc += (curr.time-prev.time)*(prev.state.tumor+curr.state.tumor)*0.5;
        }
    }
    m.final_tumor=res.points.back().state.tumor;
    // Example, arbitrarily selected objective for model exploration ONLY.
    // On a 40-time-unit horizon, AUC dominates final burden; this is NOT a clinical objective.
    m.objective = m.tumor_auc + 0.5*m.final_tumor;
    return m;
}

std::vector<Candidate> top_k_schedules(const Parameters& p, const State& init,
    const SimulationOptions& opt, const std::vector<double>& vaccine_days,
    const std::vector<double>& cart_days, int k) {
    if (k<=0) throw std::invalid_argument("k must be positive");
    struct WorseFirst {
        bool operator()(const Candidate& a,const Candidate& b) const {
            return a.metrics.objective < b.metrics.objective;
        }
    };
    std::priority_queue<Candidate,std::vector<Candidate>,WorseFirst> best;
    for (double vd : vaccine_days) {
        for (double cd : cart_days) {
            const auto doses=treatment_plan(true,true,vd,cd);
            Candidate x{vd,cd,compute_metrics(simulate(p,init,doses,opt))};
            if (static_cast<int>(best.size())<k) best.push(x);
            else if (x.metrics.objective < best.top().metrics.objective) {
                best.pop(); best.push(x);
            }
        }
    }
    std::vector<Candidate> out;
    while (!best.empty()) {out.push_back(best.top());best.pop();}
    std::sort(out.begin(),out.end(),[](const Candidate& a,const Candidate& b){
        return a.metrics.objective < b.metrics.objective;
    });
    return out;
}
} // namespace immuno
