// Algorithms B.4.4-B.4.7: separating-cycle detection.
#pragma once

#include "core/graph.hpp"

namespace apex {

// Algorithm B.4.4.  Requires isPlanar(z) == true.
bool hasSeparatingCycle(const Graph& z);

}  // namespace apex
