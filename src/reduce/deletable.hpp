// Algorithms B.1.3-B.1.5: enumerating the deletable edge sets F (Definition 3.8) of a
// multi-boundary island.
#pragma once

#include <vector>

#include "island/island.hpp"

namespace apex {

// Algorithm B.1.3. Each returned vector has length equal to the total edge count (ring +
// pendant + other) and is true exactly at the indices of a deletable F (always a subset of the
// "other" -- non-ring, non-pendant -- edges, since ring and pendant edges are always retained).
std::vector<std::vector<char>> getDeletableEdgeSet(const MultiBoundaryIsland& island);

}  // namespace apex
