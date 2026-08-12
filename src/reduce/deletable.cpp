#include "reduce/deletable.hpp"

#include <array>
#include <set>

#include "reduce/island_graph.hpp"

namespace apex {
namespace {

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

// Algorithm B.1.5.
bool isDeletable(const std::vector<char>& boolExists, const std::vector<std::vector<int>>& faces,
                 const std::vector<std::set<int>>& etoF, int totalEIandPendant) {
    for (const auto& f : faces) {
        int deletion = 0;
        std::set<int> counted;
        for (int e : f)
            if (!boolExists[e] && counted.insert(e).second) ++deletion;
        if (deletion >= 3) return true;
    }
    for (int e = 0; e < totalEIandPendant; ++e) {
        std::set<int> deletion;
        for (int f : etoF[e]) {
            for (int g : faces[f])
                if (!boolExists[g]) deletion.insert(g);
        }
        if (static_cast<int>(deletion.size()) == 4) return true;
    }
    return false;
}

// Algorithm B.1.4.
struct DeletableSearch {
    int totalEIandPendant;  // |E(I) u E_pendant|: edges 0..this-1 are ring/pendant/other-in-I
    const std::vector<std::vector<std::array<int, 2>>>& pairsAt;
    const std::vector<std::vector<int>>& faces;
    const std::vector<std::set<int>>& etoF;
    std::vector<std::vector<char>>& existsList;
    static constexpr int kMaxDeleteCount = 4;

    bool place(std::vector<int>& exists, int e0) {
        std::vector<int> q = {e0};
        size_t qi = 0;
        while (qi < q.size()) {
            const int e = q[qi++];
            for (const auto& fg : pairsAt[e]) {
                const int f = fg[0], g = fg[1];
                const std::array<int, 3> trio = {e, f, g};
                int unsetCount = 0, unsetH = -1, count1 = 0;
                for (int h : trio) {
                    if (exists[h] == -1) {
                        ++unsetCount;
                        unsetH = h;
                    } else if (exists[h] == 1) {
                        ++count1;
                    }
                }
                if (unsetCount == 0 && count1 == 1) return false;
                if (unsetCount == 1 && count1 <= 1) {
                    exists[unsetH] = count1;
                    q.push_back(unsetH);
                }
            }
        }
        return true;
    }

    void recurse(std::vector<int> exists, int e) {
        int deleteCount = 0;
        for (int x : exists)
            if (x == 0) ++deleteCount;
        if (deleteCount > kMaxDeleteCount) return;

        while (e < totalEIandPendant && exists[e] >= 0) ++e;
        if (e == totalEIandPendant) {
            std::vector<char> boolExists(totalEIandPendant);
            for (int h = 0; h < totalEIandPendant; ++h) boolExists[h] = exists[h] == 1;
            bool keep = false;
            if (deleteCount > 0 && deleteCount < kMaxDeleteCount) {
                keep = true;
            } else if (deleteCount == kMaxDeleteCount &&
                      isDeletable(boolExists, faces, etoF, totalEIandPendant)) {
                keep = true;
            }
            if (keep) {
                // Public convention (see deletable.hpp): true = deleted (in F), the complement
                // of boolExists's "true = retained" (Place/isDeletable's pseudocode convention).
                std::vector<char> inF(totalEIandPendant);
                for (int h = 0; h < totalEIandPendant; ++h) inF[h] = !boolExists[h];
                existsList.push_back(std::move(inF));
            }
            return;
        }
        for (int a = 0; a <= 1; ++a) {
            std::vector<int> exists2 = exists;
            exists2[e] = a;
            if (place(exists2, e)) recurse(exists2, e + 1);
        }
    }
};

}  // namespace

// Algorithm B.1.3.
std::vector<std::vector<char>> getDeletableEdgeSet(const MultiBoundaryIsland& island) {
    const int ringEdgeCount = ringEdgeCountOf(island);
    const int retainedCount = ringEdgeCount + island.numPendant;
    const int totalEdges = totalEdgesOf(island);
    // E(I) u E_pendant excludes only the "other" edges NOT already counted; totalEdges already
    // covers exactly E(I) u E_pendant (there is nothing beyond "other" edges in this encoding).
    const int totalEIandPendant = totalEdges;

    std::vector<std::vector<std::array<int, 2>>> pairsAt(totalEdges);
    for (int v = 0; v < island.n; ++v) {
        const auto& t = island.incidentEdges[v];
        pairsAt[t[0]].push_back({t[1], t[2]});
        pairsAt[t[1]].push_back({t[0], t[2]});
        pairsAt[t[2]].push_back({t[0], t[1]});
    }

    IslandGraph ig = buildIslandGraph(island);
    std::vector<std::vector<int>> faces = nonRingFaces(ig);
    std::vector<std::set<int>> etoF(totalEdges);
    for (size_t f = 0; f < faces.size(); ++f)
        for (int e : faces[f]) etoF[e].insert(static_cast<int>(f));

    std::vector<int> exists(totalEdges, -1);
    for (int i = 0; i < retainedCount; ++i) exists[i] = 1;

    std::vector<std::vector<char>> existsList;
    DeletableSearch search{totalEIandPendant, pairsAt, faces, etoF, existsList};
    search.recurse(std::move(exists), 0);
    return existsList;
}

}  // namespace apex
