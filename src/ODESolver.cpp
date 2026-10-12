#include "ODESolver.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace immuno {
namespace {
constexpr double eps = 1e-9;
State plus(const State& a, const State& b, double scale) {
    return {a.tumor + scale*b.tumor, a.cart + scale*b.cart,
            a.effector + scale*b.effector};
}
State checked(const State& x) {
    if (!std::isfinite(x.tumor) || !std::isfinite(x.cart) || !std::isfinite(x.effector))
        throw std::runtime_error("Numerical overflow: reduce time step or parameters");
    // Small numerical undershoots are clipped. Strongly negative values flag an unstable step.
    if (x.tumor < -1e-5 || x.cart < -1e-5 || x.effector < -1e-5)
        throw std::runtime_error("Negative state: numerical step may be too large");
    return {std::max(0.0,x.tumor), std::max(0.0,x.cart), std::max(0.0,x.effector)};
}
State one_step(const State& x, const Parameters& p, const TreatmentRates& u,
               double h, Method method) {
    const State k1 = derivatives(x,p,u);
    if (method == Method::Euler) return checked(plus(x,k1,h));
    const State k2 = derivatives(plus(x,k1,0.5*h),p,u);
    const State k3 = derivatives(plus(x,k2,0.5*h),p,u);
    const State k4 = derivatives(plus(x,k3,h),p,u);
    return checked({x.tumor+h*(k1.tumor+2*k2.tumor+2*k3.tumor+k4.tumor)/6.0,
                    x.cart+h*(k1.cart+2*k2.cart+2*k3.cart+k4.cart)/6.0,
                    x.effector+h*(k1.effector+2*k2.effector+2*k3.effector+k4.effector)/6.0});
}
} // namespace

SimulationResult simulate(const Parameters& params, const State& initial,
                          const std::vector<Dose>& doses, const SimulationOptions& opt) {
    validate(params);
    if (!(opt.end_time>0.0 && opt.step>0.0 && opt.output_interval>0.0) ||
        !std::isfinite(opt.end_time) || !std::isfinite(opt.step) || !std::isfinite(opt.output_interval))
        throw std::invalid_argument("end_time, step and output_interval must be positive finite numbers");
    if (initial.tumor<0.0 || initial.cart<0.0 || initial.effector<0.0 ||
        !std::isfinite(initial.tumor) || !std::isfinite(initial.cart) || !std::isfinite(initial.effector))
        throw std::invalid_argument("initial states must be finite and nonnegative");

    TreatmentScheduler schedule(doses);
    schedule.advance_to(0.0);
    double t=0.0;
    double next_output=opt.output_interval;
    State x=initial;
    SimulationResult result;
    result.points.push_back({t,x});
    // Step exactly to output times or rate-event boundaries to avoid RK4 crossing discontinuities.
    while (t < opt.end_time-eps) {
        double boundary = std::min({opt.end_time, t+opt.step, next_output, schedule.next_event_time()});
        if (boundary <= t+eps) {
            schedule.advance_to(t);
            if (next_output <= t+eps) {
                if (std::abs(result.points.back().time-t)>eps) result.points.push_back({t,x});
                next_output+=opt.output_interval;
            }
            continue;
        }
        x=one_step(x,params,schedule.rates(),boundary-t,opt.method);
        t=boundary;
        schedule.advance_to(t);
        if (t>=next_output-eps || t>=opt.end_time-eps) {
            if (std::abs(result.points.back().time-t)>eps) result.points.push_back({t,x});
            while (next_output <= t+eps) next_output+=opt.output_interval;
        }
    }
    return result;
}
} // namespace immuno
