#include "core/rotations.hpp"

#include <stdexcept>

namespace apex {

Graph fromVRotations(int N, const std::vector<std::vector<int>>& rotations,
                     const std::vector<std::pair<int, int>>& digons) {
    std::vector<char> isDigon(static_cast<size_t>(N) * N, 0);
    for (auto [a, b] : digons) {
        isDigon[static_cast<size_t>(a) * N + b] = 1;
        isDigon[static_cast<size_t>(b) * N + a] = 1;
    }

    Graph g;
    for (int i = 0; i < N; ++i) g.addVertex();

    std::vector<int> darts(static_cast<size_t>(N) * N, NIL);
    std::vector<int> digonDarts(static_cast<size_t>(N) * N, NIL);
    auto at = [N](int a, int b) { return static_cast<size_t>(a) * N + b; };

    for (int a = 0; a < N; ++a) {
        for (int b : rotations[a]) {
            if (b == NIL) continue;
            if (darts[at(a, b)] != NIL)
                throw std::runtime_error("fromVRotations: repeated neighbour outside a digon");
            darts[at(a, b)] = g.addDart();
            if (isDigon[at(a, b)]) digonDarts[at(a, b)] = g.addDart();
        }
    }

    for (int a = 0; a < N; ++a) {
        const auto& rot = rotations[a];
        const int size = static_cast<int>(rot.size());
        for (int i = 0; i < size; ++i) {
            const int b = rot[i];
            if (b == NIL) continue;
            if (darts[at(b, a)] == NIL)
                throw std::runtime_error("fromVRotations: rotations disagree");

            const int s = rot[i + 1 < size ? i + 1 : 0];
            const int succDart = (s == NIL) ? NIL : darts[at(a, s)];
            const int p = rot[i > 0 ? i - 1 : size - 1];
            int predDart = NIL;
            if (p != NIL) predDart = isDigon[at(a, p)] ? digonDarts[at(a, p)] : darts[at(a, p)];

            const int e = darts[at(a, b)];
            if (isDigon[at(a, b)]) {
                const int e1 = digonDarts[at(a, b)];
                g.head[e] = a;
                g.rev[e] = digonDarts[at(b, a)];
                g.succ[e] = e1;
                g.pred[e] = predDart;
                g.head[e1] = a;
                g.rev[e1] = darts[at(b, a)];
                g.succ[e1] = succDart;
                g.pred[e1] = e;
            } else {
                g.head[e] = a;
                g.rev[e] = darts[at(b, a)];
                g.succ[e] = succDart;
                g.pred[e] = predDart;
            }
        }
    }
    return g;
}

}  // namespace apex
