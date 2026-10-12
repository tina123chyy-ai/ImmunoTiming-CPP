#include "TreatmentScheduler.h"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace immuno {

TreatmentScheduler::TreatmentScheduler(const std::vector<Dose>& doses) {
    for (const Dose& d : doses) {
        if (!std::isfinite(d.start) || !std::isfinite(d.duration) ||
            !std::isfinite(d.rate) || d.start < 0.0 || d.duration <= 0.0 || d.rate < 0.0) {
            throw std::invalid_argument("Dose start/rate must be >=0 and duration >0; all must be finite");
        }
        events_.push(Event{d.start, d.type, d.rate});
        events_.push(Event{d.start+d.duration, d.type, -d.rate});
    }
}

void TreatmentScheduler::advance_to(double time) {
    constexpr double eps = 1e-10;
    while (!events_.empty() && events_.top().time <= time + eps) {
        const Event ev = events_.top(); events_.pop();
        if (ev.type == TreatmentType::Vaccine) active_.vaccine += ev.rate_delta;
        else active_.cart += ev.rate_delta;
    }
    if (std::abs(active_.vaccine) < eps) active_.vaccine = 0.0;
    if (std::abs(active_.cart) < eps) active_.cart = 0.0;
}

double TreatmentScheduler::next_event_time() const {
    return events_.empty() ? std::numeric_limits<double>::infinity() : events_.top().time;
}

std::vector<Dose> treatment_plan(bool vaccine, bool cart, double vaccine_day,
                                 double cart_day, double vaccine_rate, double cart_rate) {
    std::vector<Dose> out;
    // Box-shaped rate pulses over finite intervals, not instantaneous clinical doses.
    if (vaccine) out.push_back({TreatmentType::Vaccine, vaccine_day, 2.0, vaccine_rate});
    if (cart) out.push_back({TreatmentType::Cart, cart_day, 1.0, cart_rate});
    return out;
}
} // namespace immuno
