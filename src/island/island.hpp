// Algorithms B.4.11-B.4.16: turning an outer extension into a multi-boundary island.
#pragma once

#include <array>
#include <vector>

#include "core/graph.hpp"

namespace apex {

// Algorithm B.4.11.  K-hat must be connected with every boundary vertex of degree one
// (Claim 7.2).
Graph freeCompletionFromOuterExtension(const Graph& khat);

// A multi-boundary island, per the "Multi Boundary Island Format" of FORMAT.md.
struct MultiBoundaryIsland {
    int n = 0;                                  // number of degree-2/3 vertices
    std::vector<int> ringSizes;                  // R_1, ..., R_k = |F_R(I)| entries
    int numPendant = 0;                          // M: dummy edges, one per degree-2 vertex
    std::vector<std::array<int, 3>> incidentEdges;  // per vertex, its 3 incident edge indices
};

// Algorithm B.4.12 (uses B.4.13-B.4.16 internally).  Every closed walk of S has length 2 (a
// digon face, giving a degree-2 island vertex) or 3 (a triangle face, degree-3 vertex); every
// open (boundary) walk of length > 1 is a ring, and every open walk of length 1 is a degree-1
// pendant marker.
MultiBoundaryIsland islandFromFreeCompletion(const Graph& s);

}  // namespace apex
