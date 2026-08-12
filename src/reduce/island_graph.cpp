#include "reduce/island_graph.hpp"

#include <map>
#include <stdexcept>

#include "core/rotations.hpp"
#include "core/walks.hpp"

namespace apex {

IslandGraph buildIslandGraph(const MultiBoundaryIsland& island) {
    IslandGraph ig;
    int ringEdgeCount = 0;
    for (int r : island.ringSizes) ringEdgeCount += r;
    ig.ringEdgeCount = ringEdgeCount;
    ig.pendantEdgeCount = island.numPendant;

    int totalEdges = 0;
    for (const auto& t : island.incidentEdges)
        for (int e : t) totalEdges = std::max(totalEdges, e + 1);
    ig.totalEdges = totalEdges;

    // occurrences[e] = list of real vertex indices whose incidentEdges triple mentions e.
    std::vector<std::vector<int>> occurrences(totalEdges);
    for (int v = 0; v < island.n; ++v)
        for (int e : island.incidentEdges[v]) occurrences[e].push_back(v);

    std::vector<int> leafId(totalEdges, NIL);
    int nextLeaf = island.n;
    for (int e = 0; e < totalEdges; ++e) {
        if (occurrences[e].size() == 1) leafId[e] = nextLeaf++;
    }
    const int totalVertices = nextLeaf;

    // Every closed 2-walk (digon face) of the free completion that produced this island was
    // already turned into a single degree-2 island vertex with a padding pendant edge by
    // islandFromFreeCompletion (B.4.12/B.4.16) -- so I itself, as reconstructed here, should
    // never have two of its own real vertices joined by a genuine double edge. fromVRotations
    // throws "repeated neighbour outside a digon" if that assumption is ever wrong, so no
    // digon list is passed here.
    std::vector<std::pair<int, int>> digons;

    std::vector<std::vector<int>> rot(totalVertices);
    // Track, for the (a, b) neighbour slot at real vertex a, which edge index produced it, so
    // darts can be mapped back to edge indices after fromVRotations builds the Graph.
    std::map<std::pair<int, int>, int> slotEdge;
    for (int v = 0; v < island.n; ++v) {
        std::vector<int> r;
        r.reserve(3);
        for (int e : island.incidentEdges[v]) {
            int nb;
            if (occurrences[e].size() == 1) {
                nb = leafId[e];
            } else {
                nb = occurrences[e][0] == v ? occurrences[e][1] : occurrences[e][0];
            }
            r.push_back(nb);
            slotEdge[{v, nb}] = e;
        }
        rot[v] = std::move(r);
    }
    for (int e = 0; e < totalEdges; ++e) {
        if (occurrences[e].size() != 1) continue;
        const int leaf = leafId[e];
        const int v = occurrences[e][0];
        rot[leaf] = {v};
        slotEdge[{leaf, v}] = e;
    }

    ig.g = fromVRotations(totalVertices, rot, digons);
    ig.dartToEdge.assign(ig.g.nd(), NIL);
    for (int d = 0; d < ig.g.nd(); ++d) {
        const int a = ig.g.head[d];
        const int b = ig.g.tail(d);
        auto it = slotEdge.find({a, b});
        if (it == slotEdge.end())
            throw std::runtime_error("buildIslandGraph: missing slot-to-edge mapping");
        ig.dartToEdge[d] = it->second;
    }
    return ig;
}

std::vector<std::vector<int>> nonRingFaces(const IslandGraph& ig) {
    const auto walks = getWalks(ig.g);
    std::vector<std::vector<int>> faces;
    for (const auto& w : walks) {
        std::vector<int> edges;
        edges.reserve(w.size());
        bool isRing = false;
        for (int d : w) {
            const int e = ig.dartToEdge[d];
            edges.push_back(e);
            if (e < ig.ringEdgeCount) isRing = true;
        }
        if (!isRing) faces.push_back(std::move(edges));
    }
    return faces;
}

}  // namespace apex
