#include "Output.h"
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace immuno {
namespace {
std::ofstream make_file(const std::filesystem::path& path) {
    if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Could not open output file: "+path.string());
    out << std::setprecision(12);
    return out;
}
}
void write_trajectory(const std::filesystem::path& path, const SimulationResult& data) {
    auto out=make_file(path);
    out << "time,tumor,cart,effector\n";
    for(const auto& p:data.points) out << p.time << ',' << p.state.tumor << ','
                                        << p.state.cart << ',' << p.state.effector << '\n';
}
void write_metrics(const std::filesystem::path& path,
                   const std::vector<std::pair<std::string,Metrics>>& rows) {
    auto out=make_file(path);
    out << "scenario,tumor_auc,final_tumor,peak_cart,peak_effector,objective\n";
    for(const auto& row: rows) {
        const auto& m=row.second;
        out << row.first << ',' << m.tumor_auc << ',' << m.final_tumor << ','
            << m.peak_cart << ',' << m.peak_effector << ',' << m.objective << '\n';
    }
}
void write_candidates(const std::filesystem::path& path,const std::vector<Candidate>& rows) {
    auto out=make_file(path);
    out << "rank,vaccine_day,cart_day,tumor_auc,final_tumor,objective\n";
    for(size_t i=0;i<rows.size();++i) {
        const auto& c=rows[i];
        out << i+1 << ',' << c.vaccine_day << ',' << c.cart_day << ','
            << c.metrics.tumor_auc << ',' << c.metrics.final_tumor << ','
            << c.metrics.objective << '\n';
    }
}
} // namespace immuno
