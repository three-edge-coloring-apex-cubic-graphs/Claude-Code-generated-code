#include "core/walks.hpp"

namespace apex {

std::vector<std::vector<int>> getWalks(const Graph& g) {
    std::vector<char> visited(g.nd(), 0);
    std::vector<std::vector<int>> walks;
    for (int e = 0; e < g.nd(); ++e) {
        if (visited[e]) continue;
        std::vector<int> w;
        int cur = e;
        do {
            w.push_back(cur);
            visited[cur] = 1;
            if (g.succ[cur] == NIL) {
                cur = g.rev[firstDart(g, g.head[cur])];
            } else {
                cur = g.rev[g.succ[cur]];
            }
        } while (cur != e);
        walks.push_back(std::move(w));
    }
    return walks;
}

bool isPlanar(const Graph& g) {
    const auto walks = getWalks(g);
    // Only vertices that actually carry a dart participate in the embedding.
    std::vector<char> used(g.nv, 0);
    for (int e = 0; e < g.nd(); ++e) used[g.head[e]] = 1;
    int nv = 0;
    for (int v = 0; v < g.nv; ++v) nv += used[v];
    return nv - g.nd() / 2 + static_cast<int>(walks.size()) == 2;
}

}  // namespace apex
