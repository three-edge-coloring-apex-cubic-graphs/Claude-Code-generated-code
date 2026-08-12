#include "reduce/coloring.hpp"

#include <functional>
#include <set>

namespace apex {
namespace {

struct VertexConstraint {
    std::vector<int> edges;  // active edge indices at this vertex, size >= 2
    bool requireEqual;       // true: all must share one colour; false: pairwise distinct
};

// Generic backtracking three-edge-colouring search over `domain` (the edges that need a
// colour, in ascending order), subject to `constraints` (each one vertex's active-edge group).
// Returns the set of distinct colourings restricted to `ringEdgeCount` leading ring positions.
//
// computeCI/computeCIModuloF only need the *set* of boundary colourings that some completion
// achieves, not every completion achieving each one -- and both callers put the ring edges
// first in `domain` (indices 0 until ringEdgeCount), followed by the "other" edges. So once the
// recursion is past the ring prefix, it only needs to know a completion of the remaining "other"
// edges *exists* for the ring pattern fixed so far; it can stop at the first one instead of
// enumerating every internal completion that happens to project to the same boundary. For
// islands with a lot of internal colouring freedom this is the dominant cost (confirmed live: a
// worker thread stuck on one of the largest islands in I was sampled 22-39 recursive frames deep
// inside this search, well past the point where the boundary was already fixed) -- skipping
// same-boundary alternatives is a large, correctness-preserving speedup, not an approximation:
// the ring-restricted *set* returned is identical either way, only the redundant enumeration of
// how the interior achieves it is dropped.
std::vector<std::vector<int>> search(int totalEdges, const std::vector<int>& domain,
                                     const std::vector<VertexConstraint>& constraints,
                                     int ringEdgeCount) {
    std::vector<std::vector<int>> groupsOf(totalEdges);
    for (size_t g = 0; g < constraints.size(); ++g)
        for (int e : constraints[g].edges) groupsOf[e].push_back(static_cast<int>(g));

    std::vector<int> colour(totalEdges, -1);
    std::set<std::vector<int>> results;

    // Returns true if some completion of domain[di..] was found (and, when di > ringEdgeCount's
    // position in domain, that's all the caller needs -- further alternatives for the same ring
    // prefix are skipped).
    std::function<bool(size_t)> rec = [&](size_t di) -> bool {
        if (di == domain.size()) {
            results.emplace(colour.begin(), colour.begin() + ringEdgeCount);
            return true;
        }
        const int e = domain[di];
        const bool pastRingPrefix = e >= ringEdgeCount;
        for (int c = 0; c < 3; ++c) {
            bool ok = true;
            for (int gi : groupsOf[e]) {
                const auto& grp = constraints[gi];
                for (int f : grp.edges) {
                    if (f == e || colour[f] < 0) continue;
                    const bool same = (colour[f] == c);
                    if (grp.requireEqual ? !same : same) {
                        ok = false;
                        break;
                    }
                }
                if (!ok) break;
            }
            if (!ok) continue;
            colour[e] = c;
            const bool found = rec(di + 1);
            colour[e] = -1;
            if (found && pastRingPrefix) return true;
        }
        return false;
    };
    rec(0);
    return {results.begin(), results.end()};
}

int ringEdgeCountOf(const MultiBoundaryIsland& island) {
    int r = 0;
    for (int x : island.ringSizes) r += x;
    return r;
}

int totalEdgesOf(const MultiBoundaryIsland& island) {
    int m = 0;
    for (const auto& t : island.incidentEdges)
        for (int e : t) m = std::max(m, e + 1);
    return m;
}

}  // namespace

std::vector<std::vector<int>> computeCI(const MultiBoundaryIsland& island) {
    const int ringEdgeCount = ringEdgeCountOf(island);
    const int totalEdges = totalEdgesOf(island);

    std::vector<VertexConstraint> constraints;
    constraints.reserve(island.n);
    for (int v = 0; v < island.n; ++v) {
        constraints.push_back({{island.incidentEdges[v][0], island.incidentEdges[v][1],
                                island.incidentEdges[v][2]},
                               false});
    }
    std::vector<int> domain(totalEdges);
    for (int e = 0; e < totalEdges; ++e) domain[e] = e;
    return search(totalEdges, domain, constraints, ringEdgeCount);
}

std::vector<std::vector<int>> computeCIModuloF(const MultiBoundaryIsland& island,
                                               const std::vector<char>& deletedOther) {
    const int ringEdgeCount = ringEdgeCountOf(island);
    const int pendantEnd = ringEdgeCount + island.numPendant;
    const int totalEdges = totalEdgesOf(island);

    std::vector<VertexConstraint> constraints;
    constraints.reserve(island.n);
    std::vector<char> inDomain(totalEdges, 0);
    for (int e = 0; e < ringEdgeCount; ++e) inDomain[e] = 1;

    for (int v = 0; v < island.n; ++v) {
        std::vector<int> active;
        int kv = 0;
        for (int e : island.incidentEdges[v]) {
            if (e < pendantEnd) {
                if (e < ringEdgeCount) active.push_back(e);
                continue;  // pendant: not part of E(I), excluded from the search.
            }
            if (deletedOther[e]) {
                ++kv;
                continue;  // e in F: excluded from the search.
            }
            active.push_back(e);
            inDomain[e] = 1;
        }
        if (active.size() >= 2) constraints.push_back({std::move(active), kv == 1});
    }

    std::vector<int> domain;
    for (int e = 0; e < totalEdges; ++e)
        if (inDomain[e]) domain.push_back(e);
    return search(totalEdges, domain, constraints, ringEdgeCount);
}

}  // namespace apex
