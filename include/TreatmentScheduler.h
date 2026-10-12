#pragma once
#include "TumorModel.h"
#include <queue>
#include <vector>

namespace immuno {

enum class TreatmentType { Vaccine, Cart };

struct Dose {
    TreatmentType type{TreatmentType::Vaccine};
    double start{0.0};
    double duration{1.0};
    double rate{0.0};
};

struct Event {
    double time{0.0};
    TreatmentType type{TreatmentType::Vaccine};
    double rate_delta{0.0};
};

struct EarlierFirst {
    bool operator()(const Event& a, const Event& b) const {
        return a.time > b.time;
    }
};

class TreatmentScheduler {
public:
    explicit TreatmentScheduler(const std::vector<Dose>& doses);
    void advance_to(double time);   // apply every event up through time
    double next_event_time() const; // infinity if no event remains
    TreatmentRates rates() const { return active_; }
    bool empty() const { return events_.empty(); }
private:
    std::priority_queue<Event, std::vector<Event>, EarlierFirst> events_;
    TreatmentRates active_{};
};

std::vector<Dose> treatment_plan(bool use_vaccine, bool use_cart,
                                 double vaccine_day, double cart_day,
                                 double vaccine_rate=0.20, double cart_rate=0.18);
} // namespace immuno
