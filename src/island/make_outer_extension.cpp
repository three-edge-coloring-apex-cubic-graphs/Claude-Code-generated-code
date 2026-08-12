#include "island/make_outer_extension.hpp"

#include <array>
#include <deque>
#include <optional>

namespace apex {
namespace {

// Algorithm B.4.9.
std::optional<std::array<int, 4>> findFourDarts(const Graph& z) {
    for (int e0 = 0; e0 < z.nd(); ++e0) {
        if (z.succ[e0] == NIL) continue;
        std::array<int, 4> e{e0, NIL, NIL, NIL};
        bool foundBoundaryVertex = false;
        for (int i = 0; i < 3; ++i) {
            if ((i == 1 || i == 2) && z.succ[e[i]] == NIL) {
                foundBoundaryVertex = true;
                break;
            }
            e[i + 1] = z.rev[z.succ[e[i]]];
        }
        if (!foundBoundaryVertex && e0 != e[2] && e0 != e[3]) return e;
    }
    return std::nullopt;
}

// Algorithm B.4.10.
std::vector<HomImage> ensureOuterExtension(const Graph& z, const std::array<int, 4>& e) {
    auto z3 = freeHomomorphismAndEnforceSingleDigonIncidence(z, {{e[0], e[3]}},
                                                             Boundary::PseudoEmbedding);
    auto z2 = freeHomomorphismAndEnforceSingleDigonIncidence(z, {{e[0], e[2]}},
                                                             Boundary::PseudoEmbedding);
    std::vector<HomImage> out = std::move(z3);
    for (auto& r : z2) out.push_back(std::move(r));
    return out;
}

}  // namespace

// Algorithm B.4.8.
std::vector<HomImage> makeOuterExtension(const Graph& z) {
    std::vector<HomImage> s;
    std::deque<HomImage> q;
    q.push_back({z, identityMap(z)});
    while (!q.empty()) {
        HomImage cur = std::move(q.front());
        q.pop_front();
        auto a = findFourDarts(cur.g);
        if (a) {
            for (auto& r : ensureOuterExtension(cur.g, *a))
                q.push_back({std::move(r.g), compose(r.phi, cur.phi)});
        } else {
            s.push_back(std::move(cur));
        }
    }
    return s;
}

}  // namespace apex
