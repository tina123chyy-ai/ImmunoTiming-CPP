#include "TumorModel.h"
#include <cmath>
#include <stdexcept>

namespace immuno {

State derivatives(const State& x, const Parameters& p, const TreatmentRates& u) {
    const double t = x.tumor, c = x.cart, e = x.effector;
    const double tumor_change = p.growth_rate*t*(1.0-t/p.carrying_capacity)
        - p.cart_killing*c*t - p.effector_killing*e*t;
    const double cart_change = p.cart_expansion*(1.0+p.synergy_gamma*e)
        * (t/(p.antigen_half_sat+t))*c - p.cart_decay*c + u.cart;
    const double effector_change = p.vaccine_activation*u.vaccine
        + p.antigen_feedback_beta*c*t - p.effector_decay*e;
    return {tumor_change, cart_change, effector_change};
}

void validate(const Parameters& p) {
    const double vals[] = {p.growth_rate,p.carrying_capacity,p.cart_killing,
        p.effector_killing,p.cart_expansion,p.antigen_half_sat,p.cart_decay,
        p.synergy_gamma,p.vaccine_activation,p.antigen_feedback_beta,p.effector_decay};
    for (double v : vals) {
        if (!std::isfinite(v) || v < 0.0) throw std::invalid_argument("Model parameters must be finite and nonnegative");
    }
    if (p.carrying_capacity <= 0.0 || p.antigen_half_sat <= 0.0) {
        throw std::invalid_argument("carrying_capacity and antigen_half_sat must be positive");
    }
}

std::string parameter_description() {
    return "Illustrative normalized parameters only; NOT fitted to biological/clinical data.";
}
} // namespace immuno
