#include "island/planarity.hpp"

#include <functional>

namespace apex {
namespace {

// Algorithm B.4.5.
std::vector<std::vector<int>> enumCycles(const Graph& z, int l) {
    std::vector<char> visitedV(z.nv, 0), visitedD(z.nd(), 0);
    std::vector<std::vector<int>> cycles;
    std::vector<int> path;

    std::function<void()> dfs = [&]() {
        if (static_cast<int>(path.size()) > l) return;
        const int e = path.back();
        const int h = z.head[e];
        for (int e1 = 0; e1 < z.nd(); ++e1) {
            if (z.head[e1] != h || e1 == e) continue;
            const int f = z.rev[e1];
            if (f == path.front()) {
                cycles.push_back(path);
                continue;
            }
            if (visitedD[f] || visitedV[z.head[f]]) continue;
            path.push_back(f);
            visitedV[z.head[f]] = 1;
            dfs();
            visitedV[z.head[f]] = 0;
            path.pop_back();
        }
    };

    for (int e = 0; e < z.nd(); ++e) {
        path = {e};
        visitedV[z.head[e]] = 1;
        dfs();
        visitedV[z.head[e]] = 0;
        visitedD[e] = 1;
    }
    return cycles;
}

enum class Label { None, L, R };

// Algorithm B.4.6.
std::vector<Label> labelDarts(const Graph& z, const std::vector<int>& c) {
    std::vector<Label> ld(z.nd(), Label::None);
    const int l = static_cast<int>(c.size());
    for (int d : c) {
        ld[d] = Label::L;
        ld[z.rev[d]] = Label::R;
    }

    std::function<void(int, Label)> propagate = [&](int e, Label delta) {
        ld[e] = delta;
        for (int f : {z.succ[e], z.pred[e], z.rev[e]}) {
            if (f != NIL && ld[f] == Label::None) propagate(f, delta);
        }
    };

    for (int i = 0; i < l; ++i) {
        const int di = c[i];
        const int diNext = c[(i + 1) % l];
        if (z.succ[di] != NIL && z.succ[di] != z.rev[diNext]) propagate(z.succ[di], Label::L);
        const int revNext = z.rev[diNext];
        if (z.succ[revNext] != NIL && z.succ[revNext] != di) propagate(z.succ[revNext], Label::R);
    }
    return ld;
}

// Algorithm B.4.7.
std::pair<int, int> numSeparatedVertices(const Graph& z, const std::vector<int>& c,
                                         const std::vector<Label>& ld) {
    std::vector<char> inC(z.nv, 0);
    for (int d : c) inC[z.head[d]] = 1;

    std::vector<char> inVL(z.nv, 0), inVR(z.nv, 0);
    for (int d = 0; d < z.nd(); ++d) {
        const int h = z.head[d];
        if (inC[h]) continue;
        if (ld[d] == Label::L) inVL[h] = 1;
        else if (ld[d] == Label::R) inVR[h] = 1;
    }

    int nLInner = 0, nLBoundary = 0, nRInner = 0, nRBoundary = 0;
    for (int v = 0; v < z.nv; ++v) {
        const bool boundary = isBoundaryVertex(z, v);
        if (inVL[v]) (boundary ? nLBoundary : nLInner)++;
        if (inVR[v]) (boundary ? nRBoundary : nRInner)++;
    }
    const int nL = nLBoundary > 0 ? nLInner + 1 : nLInner;
    const int nR = nRBoundary > 0 ? nRInner + 1 : nRInner;
    return {nL, nR};
}

}  // namespace

// Algorithm B.4.4.
bool hasSeparatingCycle(const Graph& z) {
    for (const auto& c : enumCycles(z, 4)) {
        const auto ld = labelDarts(z, c);
        const auto [nL, nR] = numSeparatedVertices(z, c, ld);
        const int len = static_cast<int>(c.size());
        if (len <= 3 && nL > 0 && nR > 0) return true;
        if (len == 4 && nL > 2 && nR > 2) return true;
    }
    return false;
}

}  // namespace apex
