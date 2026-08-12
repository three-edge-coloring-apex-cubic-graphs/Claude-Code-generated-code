// Reconstructs the actual planar graph I (with genuine degree-one leaves for both the ring
// edges of Definition 3.1 and the auxiliary E_pendant edges of B.1.3-B.1.5) from the compact
// "N vertices x 3 incident edge indices" MultiBoundaryIsland encoding of FORMAT.md. Needed to
// trace faces F(I) for the deletable-edge-set algorithms; the edge-colouring search itself
// (Definition 3.6/3.9) only needs the compact incidentEdges array directly.
#pragma once

#include <vector>

#include "core/graph.hpp"
#include "island/island.hpp"

namespace apex {

struct IslandGraph {
    Graph g;
    std::vector<int> dartToEdge;  // per dart, the FORMAT.md edge index it corresponds to
    int totalEdges = 0;
    int ringEdgeCount = 0;     // sum of ringSizes: edge indices [0, ringEdgeCount)
    int pendantEdgeCount = 0;  // edge indices [ringEdgeCount, ringEdgeCount + pendantEdgeCount)
};

IslandGraph buildIslandGraph(const MultiBoundaryIsland& island);

// Faces of I other than the ring faces F_R(I) (Algorithm B.1.3/B.1.4's "F(I) \ F_R(I)"), each
// given as the list of FORMAT.md edge indices on its boundary walk (edges bounding the face
// twice, e.g. pendant/bridge edges, appear twice).
std::vector<std::vector<int>> nonRingFaces(const IslandGraph& ig);

}  // namespace apex
