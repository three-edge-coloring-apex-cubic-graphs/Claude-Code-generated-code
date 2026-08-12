#include "island/outer_extension.hpp"

#include "core/rotations.hpp"

namespace apex {

Graph outerExtension(const Configuration& c) {
    const int numInternal = c.N - c.R;
    auto toNew = [&](int orig) { return orig - c.R; };  // valid for orig >= R

    int totalVertices = numInternal;
    for (int u = c.R; u < c.N; ++u)
        for (int nb : c.rotations[u])
            if (nb < c.R) ++totalVertices;

    std::vector<std::vector<int>> rot(totalVertices);
    int nextPendant = numInternal;
    for (int u = c.R; u < c.N; ++u) {
        std::vector<int> r;
        for (int nb : c.rotations[u]) {
            if (nb >= c.R) {
                r.push_back(toNew(nb));
            } else {
                const int p = nextPendant++;
                r.push_back(p);
                rot[p] = {toNew(u), NIL};  // the outer endpoint: one dart back to u, open
            }
        }
        rot[toNew(u)] = std::move(r);
    }

    std::vector<std::pair<int, int>> digons;
    for (auto [a, b] : c.digons) digons.emplace_back(toNew(a), toNew(b));

    Graph g = fromVRotations(totalVertices, rot, digons);
    for (int v = 0; v < g.nv; ++v) {
        if (v < numInternal) {
            g.dlo[v] = c.dlo[c.R + v];
            g.dhi[v] = c.dhi[c.R + v];
        } else {
            // Outer endpoint: the genuine open boundary vertex, degree range [5, infinity).
            g.dlo[v] = 5;
            g.dhi[v] = INF;
        }
    }
    return g;
}

}  // namespace apex
