// Algorithms B.4.2 getWalks and B.4.3 isPlanar.
#pragma once

#include <vector>

#include "core/graph.hpp"

namespace apex {

std::vector<std::vector<int>> getWalks(const Graph& g);
bool isPlanar(const Graph& g);

}  // namespace apex
