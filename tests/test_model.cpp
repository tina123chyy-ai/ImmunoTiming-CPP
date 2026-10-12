#include "ScheduleRanker.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace immuno;
namespace {
void require(bool ok,const char* why) {if(!ok) throw std::runtime_error(why);}
bool near(double x,double y,double tol=1e-6){return std::abs(x-y)<=tol;}
}
int main() {
    try {
        Parameters p;
        State init{0.2,0.0,0.0};
        SimulationOptions opt;
        opt.end_time=10; opt.step=0.01; opt.output_interval=0.25;
        const auto control=simulate(p,init,{},opt);
        require(control.points.back().state.tumor>init.tumor,"No-treatment tumor must grow from below K");
        require(near(control.points.back().state.cart,0),"No CAR-T input should leave C=0");
        require(near(control.points.back().state.effector,0),"No immune input should leave E=0");
        const double exp_rt=std::exp(p.growth_rate*opt.end_time);
        const double logistic=p.carrying_capacity/(1.0+(p.carrying_capacity/init.tumor-1.0)/exp_rt);
        require(near(control.points.back().state.tumor,logistic,1e-7),"RK4 must solve no-treatment logistic growth");

        TreatmentScheduler sched({Dose{TreatmentType::Vaccine,2.0,1.0,0.2},
                                  Dose{TreatmentType::Cart,2.5,0.5,0.4}});
        sched.advance_to(1.0);require(near(sched.rates().vaccine,0),"No premature vaccine");
        sched.advance_to(2.0);require(near(sched.rates().vaccine,0.2),"Vaccine event start");
        sched.advance_to(2.5);require(near(sched.rates().cart,0.4),"CAR-T event start");
        sched.advance_to(3.0);require(near(sched.rates().vaccine,0) && near(sched.rates().cart,0),"Event end");

        const auto treated=simulate(p,init,treatment_plan(true,true,2,4),opt);
        require(treated.points.back().state.cart>=0 && treated.points.back().state.effector>=0,"States nonnegative");
        require(treated.points.size()==control.points.size(),"Same output grid");
        require(treated.points.back().state.tumor<control.points.back().state.tumor,"Expected behavior for toy parameters");
        const auto m=compute_metrics(treated);
        require(m.tumor_auc>0 && m.objective>0,"Metrics computed");

        const auto ranked=top_k_schedules(p,init,opt,{0,2,4},{1,3,5},3);
        require(ranked.size()==3,"Top K count");
        for(size_t i=1;i<ranked.size();++i)
            require(ranked[i].metrics.objective>=ranked[i-1].metrics.objective,"Top K sorted");
        bool thrown=false;
        try {auto bad=simulate(p,init,{},SimulationOptions{10,-1,1,Method::RK4});(void)bad;}
        catch(const std::invalid_argument&){thrown=true;}
        require(thrown,"Invalid steps rejected");
        // Sensitivity spot check: beta is only a *hypothetical* feedback mechanism.
        Parameters zero=p, positive=p;
        zero.antigen_feedback_beta=0;positive.antigen_feedback_beta=0.3;
        const auto a=simulate(zero,init,treatment_plan(true,true,2,4),opt);
        const auto b=simulate(positive,init,treatment_plan(true,true,2,4),opt);
        require(b.points.back().state.effector>=a.points.back().state.effector,"Feedback should affect effector levels in example");
        std::cout << "PASS: solver, logistic benchmark, scheduler, state constraints, metrics, Top-K, validation, feedback\n";
        return 0;
    } catch(const std::exception& e){
        std::cerr << "FAIL: " << e.what() << '\n';return 1;
    }
}
