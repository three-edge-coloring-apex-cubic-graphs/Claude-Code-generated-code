#include "island/island.hpp"

#include <map>

#include "core/walks.hpp"
#include "hom/hom.hpp"

namespace apex {

// Algorithm B.4.11.
Graph freeCompletionFromOuterExtension(const Graph& khat) {
    const auto walks = getWalks(khat);
    Graph z = khat;
    Map phi = identityMap(z);
    for (const auto& w : walks) {
        for (int e : w) {
            const int e0 = phi.d[e];
            if (z.pred[z.rev[e0]] != NIL) continue;
            const int e1 = z.rev[z.succ[e0]];
            if (z.succ[e1] == NIL) {
                addBoundaryDartsDirectly(z, z.rev[e1], e0);
                continue;
            }
            const int e2 = z.rev[z.succ[e1]];
            if (z.succ[e2] == NIL) {
                Map phi1 = linkIncidenceListEnds(z, z.rev[e0], e2);
                phi = compose(phi1, phi);
                continue;
            }
            // Unreachable: freeCompletionFromOuterExtension is only called on a valid outer
            // extension, whose only boundary vertices are degree-one outer endpoints.
        }
    }
    return z;
}

namespace {

// Algorithm B.4.13.
std::vector<int> indexBoundaryEdges(const Graph& s, const std::vector<std::vector<int>>& fs,
                                    std::vector<int>& dartsToEdge, int& m) {
    std::vector<int> ringSizes;
    for (const auto& w : fs) {
        if (s.succ[w[0]] == NIL && w.size() > 1) {
            for (int e : w) {
                dartsToEdge[e] = m;
                dartsToEdge[s.rev[e]] = m;
                ++m;
            }
            ringSizes.push_back(static_cast<int>(w.size()));
        }
    }
    return ringSizes;
}

// Algorithm B.4.14.
std::pair<int, std::map<int, int>> indexPendantEdges(const Graph& s,
                                                      const std::vector<std::vector<int>>& fs,
                                                      std::vector<int>& dartsToEdge, int& m) {
    int nPendant = 0;
    std::map<int, int> digonToPendant;  // keyed by the walk's index in fs
    for (int wi = 0; wi < static_cast<int>(fs.size()); ++wi) {
        const auto& w = fs[wi];
        const int e0 = w[0];
        if (s.succ[e0] == NIL && w.size() == 1) {
            dartsToEdge[e0] = m;
            dartsToEdge[s.rev[e0]] = m;
            ++m;
            ++nPendant;
        }
        if (s.succ[e0] != NIL && w.size() == 2) {
            digonToPendant[wi] = m;
            ++m;
            ++nPendant;
        }
    }
    return {nPendant, digonToPendant};
}

// Algorithm B.4.15.
void indexOtherEdges(const Graph& s, const std::vector<std::vector<int>>& fs,
                     std::vector<int>& dartsToEdge, int& m) {
    for (const auto& w : fs) {
        if (s.succ[w[0]] == NIL) continue;
        for (int e : w) {
            if (dartsToEdge[e] == NIL && dartsToEdge[s.rev[e]] == NIL) {
                dartsToEdge[e] = m;
                dartsToEdge[s.rev[e]] = m;
                ++m;
            }
        }
    }
}

// Algorithm B.4.16.
MultiBoundaryIsland constructIsland(const Graph& s, const std::vector<std::vector<int>>& fs,
                                    const std::vector<int>& dartsToEdge,
                                    std::vector<int> ringSizes, int nPendantEdge,
                                    const std::map<int, int>& digonToPendantEdge) {
    MultiBoundaryIsland island;
    island.ringSizes = std::move(ringSizes);
    island.numPendant = nPendantEdge;
    for (int wi = 0; wi < static_cast<int>(fs.size()); ++wi) {
        const auto& w = fs[wi];
        if (s.succ[w[0]] == NIL) continue;
        if (w.size() == 2) {
            island.incidentEdges.push_back(
                {dartsToEdge[w[0]], dartsToEdge[w[1]], digonToPendantEdge.at(wi)});
        } else if (w.size() == 3) {
            island.incidentEdges.push_back(
                {dartsToEdge[w[0]], dartsToEdge[w[1]], dartsToEdge[w[2]]});
        }
    }
    island.n = static_cast<int>(island.incidentEdges.size());
    return island;
}

}  // namespace

// Algorithm B.4.12.
MultiBoundaryIsland islandFromFreeCompletion(const Graph& s) {
    const auto fs = getWalks(s);
    int m = 0;
    std::vector<int> dartsToEdge(s.nd(), NIL);
    std::vector<int> ringSizes = indexBoundaryEdges(s, fs, dartsToEdge, m);
    auto [nPendantEdge, digonToPendantEdge] = indexPendantEdges(s, fs, dartsToEdge, m);
    indexOtherEdges(s, fs, dartsToEdge, m);
    return constructIsland(s, fs, dartsToEdge, std::move(ringSizes), nPendantEdge,
                           digonToPendantEdge);
}

}  // namespace apex
