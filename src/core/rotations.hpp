// Algorithm B.2.10 fromVRotations (the digon-aware generalisation of Algorithm A.5.1).
#pragma once

#include <utility>
#include <vector>

#include "core/graph.hpp"

namespace apex {

// rotations[i] is the clockwise rotation of neighbours around v_i; -1 marks the boundary.
// digons lists the endpoint pairs of the digons.  Throws std::runtime_error when a vertex is
// adjacent to the same vertex twice outside a digon, or when two rotations disagree.
// Degree ranges are left at [1, INF]; callers overwrite them.
Graph fromVRotations(int N, const std::vector<std::vector<int>>& rotations,
                     const std::vector<std::pair<int, int>>& digons);

}  // namespace apex
