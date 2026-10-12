#pragma once
#include <string>

namespace immuno {

struct State {
    double tumor{0.25};
    double cart{0.0};
    double effector{0.0};
};

struct Parameters {
    double growth_rate{0.15};       // r: logistic tumor growth
    double carrying_capacity{1.0}; // K: tumor burden scale
    double cart_killing{0.45};      // kC
    double effector_killing{0.25};  // kE
    double cart_expansion{0.90};    // alphaC
    double antigen_half_sat{0.30}; // hC
    double cart_decay{0.34};       // deltaC
    double synergy_gamma{0.38};    // gamma: hypothesized E->C enhancement
    double vaccine_activation{0.90}; // alphaV
    double antigen_feedback_beta{0.16}; // beta: hypothesized C*T->E feedback
    double effector_decay{0.20};   // deltaE
};

struct TreatmentRates {
    double vaccine{0.0};  // V(t): normalized vaccine administration rate
    double cart{0.0};     // uC(t): normalized CAR-T infusion rate
};

State derivatives(const State& x, const Parameters& p, const TreatmentRates& u);
void validate(const Parameters& p);
std::string parameter_description();

} // namespace immuno
