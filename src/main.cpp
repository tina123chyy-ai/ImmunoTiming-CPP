#include "Output.h"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace immuno;
namespace {
void usage() {
    std::cout << R"(ImmunoTiming-CPP: conceptual cancer-immunotherapy model (NOT clinical advice)

Usage:
  immunotiming_cli [--mode baseline|single|rank|sensitivity]
                   [--out results] [--vaccine-day 2] [--cart-day 10]
                   [--days 40] [--step 0.05] [--method rk4|euler] [--k 5]

Examples:
  ./immunotiming_cli --mode baseline --out results/baseline
  ./immunotiming_cli --mode single --vaccine-day 2 --cart-day 10
  ./immunotiming_cli --mode rank --k 5 --out results/rank
  ./immunotiming_cli --mode sensitivity --out results/sensitivity

All parameters and doses are invented/normalized. Results only illustrate model behavior.
)";
}
double parse_positive(const std::string& s, const std::string& name, bool allow_zero=false) {
    size_t read=0;
    double v=std::stod(s,&read);
    if (read!=s.size() || !std::isfinite(v) || (allow_zero ? v<0.0 : v<=0.0))
        throw std::invalid_argument("Invalid value for "+name+": "+s);
    return v;
}
}
int main(int argc,char* argv[]) {
    try {
        std::string mode="baseline", out="results", method="rk4";
        double vaccine_day=2.0, cart_day=10.0;
        int k=5;
        SimulationOptions opt;
        for(int i=1;i<argc;++i) {
            std::string arg=argv[i];
            if(arg=="--help"||arg=="-h") {usage();return 0;}
            if (i+1>=argc) throw std::invalid_argument("Missing value after "+arg);
            std::string value=argv[++i];
            if(arg=="--mode") mode=value;
            else if(arg=="--out") out=value;
            else if(arg=="--vaccine-day") vaccine_day=parse_positive(value,arg,true);
            else if(arg=="--cart-day") cart_day=parse_positive(value,arg,true);
            else if(arg=="--days") opt.end_time=parse_positive(value,arg);
            else if(arg=="--step") opt.step=parse_positive(value,arg);
            else if(arg=="--method") method=value;
            else if(arg=="--k") {
                double v=parse_positive(value,arg);
                if(v>1000 || std::floor(v)!=v) throw std::invalid_argument("k must be a positive integer <=1000");
                k=static_cast<int>(v);
            } else throw std::invalid_argument("Unknown argument: "+arg);
        }
        if(method=="rk4") opt.method=Method::RK4;
        else if(method=="euler") opt.method=Method::Euler;
        else throw std::invalid_argument("--method must be rk4 or euler");
        if(out.empty()) throw std::invalid_argument("--out must not be empty");

        const Parameters p;
        const State init{0.25,0.0,0.0};
        const std::filesystem::path folder(out);
        std::filesystem::create_directories(folder);
        if (mode=="baseline") {
            std::vector<std::pair<std::string,Metrics>> summary;
            const struct {const char* name; bool vaccine; bool cart;} settings[] = {
                {"none",false,false},{"vaccine_only",true,false},
                {"cart_only",false,true},{"combined",true,true}};
            for (const auto& a:settings) {
                auto res=simulate(p,init,treatment_plan(a.vaccine,a.cart,vaccine_day,cart_day),opt);
                write_trajectory(folder/(std::string(a.name)+".csv"),res);
                summary.emplace_back(a.name,compute_metrics(res));
            }
            write_metrics(folder/"summary.csv",summary);
            std::cout << "Created four trajectories and summary in " << folder << '\n';
        } else if(mode=="single") {
            const auto res=simulate(p,init,treatment_plan(true,true,vaccine_day,cart_day),opt);
            write_trajectory(folder/"single.csv",res);
            write_metrics(folder/"single_metrics.csv",{{"combined",compute_metrics(res)}});
            std::cout << "Created single.csv and single_metrics.csv in " << folder << '\n';
        } else if(mode=="rank") {
            const std::vector<double> vdays{0,4,8,12};
            const std::vector<double> cdays{0,2,4,6,8,10,12,14,16,18,20};
            auto ranked=top_k_schedules(p,init,opt,vdays,cdays,k);
            write_candidates(folder/"top_k.csv",ranked);
            if(!ranked.empty()) {
                const auto& top=ranked.front();
                const auto res=simulate(p,init,treatment_plan(true,true,top.vaccine_day,top.cart_day),opt);
                write_trajectory(folder/"top_candidate.csv",res);
                std::cout << "Best *model-only* candidate: vaccine=" << top.vaccine_day
                          << ", CAR-T=" << top.cart_day << ", objective=" << top.metrics.objective << '\n';
            }
            std::cout << "Created ranked candidates in " << folder << '\n';
        } else if(mode=="sensitivity") {
            std::vector<std::pair<std::string,Metrics>> summary;
            for(double beta : {0.0,0.08,0.16,0.32}) {
                for(double gamma : {0.0,0.19,0.38,0.76}) {
                    Parameters q=p;
                    q.antigen_feedback_beta=beta;
                    q.synergy_gamma=gamma;
                    const auto res=simulate(q,init,treatment_plan(true,true,vaccine_day,cart_day),opt);
                    std::string label="beta_"+std::to_string(beta)+"_gamma_"+std::to_string(gamma);
                    summary.emplace_back(label,compute_metrics(res));
                }
            }
            write_metrics(folder/"sensitivity.csv",summary);
            std::cout << "Created 16-model sensitivity.csv in " << folder << '\n';
        } else throw std::invalid_argument("--mode must be baseline, single, rank or sensitivity");
        std::cout << "Reminder: hypotheses and rates are illustrative, not validated clinical predictions.\n";
    } catch(const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\nUse --help for instructions.\n";
        return 1;
    }
}
