#include "kbar/kbar.hpp"

#include <algorithm>
#include <stdexcept>
#include <tuple>

#include "core/rotations.hpp"

namespace apex {
namespace {

// Algorithm B.2.8: the cut-vertices together with their two ring neighbours.
std::vector<std::tuple<int, int, int>> findCutTuples(const Configuration& c) {
    std::vector<std::tuple<int, int, int>> p;
    for (int i = c.R; i < c.N; ++i) {
        const auto& rot = c.rotations[i];
        const int d = static_cast<int>(rot.size());
        std::vector<int> ur;
        int t = 0;
        for (int j = 0; j < d; ++j) {
            const int k1 = rot[j];
            if (k1 != NIL && k1 < c.R && std::find(ur.begin(), ur.end(), k1) == ur.end())
                ur.push_back(k1);
            const int k2 = rot[(j + 1) % d];
            if (k1 != NIL && k1 < c.R && (k2 == NIL || k2 >= c.R)) ++t;
        }
        if (t >= 2 && ur.size() != 2)
            throw std::runtime_error(c.name + ": not a normal configuration");
        if (t == 2 && ur.size() == 2) p.emplace_back(i, ur[0], ur[1]);
    }
    return p;
}

// Algorithm A.6.4.
int maximumDegreeDart(const Graph& g) {
    int f = NIL;
    std::pair<int, int> best{0, 0};
    for (int e = 0; e < g.nd(); ++e) {
        const int y = g.head[e], x = g.tail(e);
        if (g.dlo[y] != g.dhi[y] || g.dlo[x] != g.dhi[x]) continue;
        const std::pair<int, int> de{g.dlo[y], g.dlo[x]};
        if (de > best) {
            best = de;
            f = e;
        }
    }
    return f;
}

// Algorithm B.2.9.
Graph removeRing(const Configuration& c, const std::vector<int>& adjacentCutvertex) {
    std::vector<int> old2new(c.N, NIL);
    int newId = 0;
    for (int i = 0; i < c.N; ++i) {
        if (i < c.R && adjacentCutvertex[i] == NIL) continue;
        old2new[i] = newId++;
    }
    const int n1 = newId;

    std::vector<std::vector<int>> newRot(n1);
    for (int i = 0; i < c.N; ++i) {
        if (i < c.R && adjacentCutvertex[i] == NIL) continue;
        for (const int j : c.rotations[i]) {
            const bool drop = (j == NIL) || (i < c.R && j != adjacentCutvertex[i]) ||
                              (j < c.R && i != adjacentCutvertex[j]);
            newRot[old2new[i]].push_back(drop ? NIL : old2new[j]);
        }
    }

    std::vector<std::pair<int, int>> newDigons;
    for (auto [a, b] : c.digons) {
        if (old2new[a] == NIL || old2new[b] == NIL) continue;
        newDigons.emplace_back(old2new[a], old2new[b]);
    }

    Graph g = fromVRotations(n1, newRot, newDigons);
    for (int i = 0; i < c.R; ++i) {
        if (adjacentCutvertex[i] == NIL) continue;
        const int k = old2new[i];
        int d = 0;
        for (const int x : newRot[k])
            if (x != NIL) ++d;
        g.dlo[k] = d + 1;
        g.dhi[k] = INF;
    }
    for (int i = c.R; i < c.N; ++i) {
        const int k = old2new[i];
        g.dlo[k] = c.dlo[i];
        g.dhi[k] = c.dhi[i];
    }
    return g;
}

}  // namespace

std::vector<ConfEntry> extendFromCutVertices(const Configuration& c) {
    const auto p = findCutTuples(c);
    std::vector<ConfEntry> out;
    for (int s = 0; s < (1 << p.size()); ++s) {
        std::vector<int> adjacentCutvertex(c.R, NIL);
        for (size_t i = 0; i < p.size(); ++i) {
            const auto [v, a, b] = p[i];
            if (s >> i & 1)
                adjacentCutvertex[a] = v;
            else
                adjacentCutvertex[b] = v;
        }
        Graph g = removeRing(c, adjacentCutvertex);
        out.push_back({g, maximumDegreeDart(g)});
    }
    return out;
}

ConfSet buildKBar(const std::vector<Configuration>& confs) {
    std::vector<ConfEntry> all;
    for (const auto& c : confs) {
        for (auto& e : extendFromCutVertices(c)) {
            Graph m = mirror(e.g);
            all.push_back(e);
            all.push_back({m, e.dart});
        }
    }
    return makeConfSet(std::move(all));
}

}  // namespace apex
