#pragma once
#include "ScheduleRanker.h"
#include <filesystem>
#include <string>
#include <vector>

namespace immuno {
void write_trajectory(const std::filesystem::path& path, const SimulationResult& data);
void write_metrics(const std::filesystem::path& path,
                   const std::vector<std::pair<std::string,Metrics>>& rows);
void write_candidates(const std::filesystem::path& path, const std::vector<Candidate>& rows);
} // namespace immuno
