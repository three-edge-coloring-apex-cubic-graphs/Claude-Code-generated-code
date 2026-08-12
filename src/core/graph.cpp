#include "core/graph.hpp"

#include <algorithm>
#include <deque>
#include <stdexcept>

namespace apex {

Map identityMap(const Graph& g) {
    Map m;
    m.v.resize(g.nv);
    m.d.resize(g.nd());
    for (int i = 0; i < g.nv; ++i) m.v[i] = i;
    for (int i = 0; i < g.nd(); ++i) m.d[i] = i;
    return m;
}

Map compose(const Map& second, const Map& first) {
    Map m;
    m.v.resize(first.v.size());
    m.d.resize(first.d.size());
    for (size_t i = 0; i < first.v.size(); ++i) m.v[i] = second.v[first.v[i]];
    for (size_t i = 0; i < first.d.size(); ++i) m.d[i] = second.d[first.d[i]];
    return m;
}

std::vector<std::vector<int>> incidence(const Graph& g) {
    std::vector<std::vector<int>> out(g.nv);
    std::vector<char> seen(g.nd(), 0);
    // Walk each incidence list from its first dart so the lists come out in rotation order.
    for (int e = 0; e < g.nd(); ++e) {
        if (seen[e] || g.pred[e] != NIL) continue;
        for (int f = e; f != NIL; f = g.succ[f]) {
            seen[f] = 1;
            out[g.head[f]].push_back(f);
        }
    }
    for (int e = 0; e < g.nd(); ++e) {
        if (seen[e]) continue;
        int f = e;
        do {
            seen[f] = 1;
            out[g.head[f]].push_back(f);
            f = g.succ[f];
        } while (f != NIL && f != e);
    }
    return out;
}

int dartDegree(const Graph& g, int v) {
    int c = 0;
    for (int e = 0; e < g.nd(); ++e)
        if (g.head[e] == v) ++c;
    return c;
}

bool isBoundaryVertex(const Graph& g, int v) {
    for (int e = 0; e < g.nd(); ++e)
        if (g.head[e] == v && g.succ[e] == NIL) return true;
    return false;
}

int firstDart(const Graph& g, int v) {
    for (int e = 0; e < g.nd(); ++e)
        if (g.head[e] == v && g.pred[e] == NIL) return e;
    return NIL;
}

int lastDart(const Graph& g, int v) {
    for (int e = 0; e < g.nd(); ++e)
        if (g.head[e] == v && g.succ[e] == NIL) return e;
    return NIL;
}

bool hasLoop(const Graph& g) {
    for (int e = 0; e < g.nd(); ++e)
        if (g.head[e] == g.head[g.rev[e]]) return true;
    return false;
}

Graph disjointUnion(const Graph& a, const Graph& b) {
    Graph g = a;
    const int vo = a.nv, dof = a.nd();
    g.nv = a.nv + b.nv;
    g.dlo.insert(g.dlo.end(), b.dlo.begin(), b.dlo.end());
    g.dhi.insert(g.dhi.end(), b.dhi.begin(), b.dhi.end());
    auto shift = [&](int x, int off) { return x == NIL ? NIL : x + off; };
    for (int e = 0; e < b.nd(); ++e) {
        g.head.push_back(shift(b.head[e], vo));
        g.rev.push_back(shift(b.rev[e], dof));
        g.succ.push_back(shift(b.succ[e], dof));
        g.pred.push_back(shift(b.pred[e], dof));
    }
    return g;
}

Graph mirror(const Graph& g) {
    Graph m = g;
    m.succ = g.pred;
    m.pred = g.succ;
    return m;
}

std::string canonicalKey(const Graph& g, int root) {
    std::vector<int> dnew(g.nd(), NIL), vnew(g.nv, NIL);
    std::vector<int> dorder, vorder;
    std::deque<int> q;

    auto visitDart = [&](int e) {
        if (e == NIL || dnew[e] != NIL) return;
        dnew[e] = static_cast<int>(dorder.size());
        dorder.push_back(e);
        const int h = g.head[e];
        if (vnew[h] == NIL) {
            vnew[h] = static_cast<int>(vorder.size());
            vorder.push_back(h);
        }
        q.push_back(e);
    };

    visitDart(root);
    while (!q.empty()) {
        const int e = q.front();
        q.pop_front();
        visitDart(g.rev[e]);
        visitDart(g.succ[e]);
        visitDart(g.pred[e]);
    }

    std::string s;
    s.reserve(dorder.size() * 12);
    auto put = [&](int x) {
        s += std::to_string(x);
        s += ',';
    };
    put(static_cast<int>(vorder.size()));
    put(static_cast<int>(dorder.size()));
    for (int e : dorder) {
        put(vnew[g.head[e]]);
        put(g.rev[e] == NIL ? NIL : dnew[g.rev[e]]);
        put(g.succ[e] == NIL ? NIL : dnew[g.succ[e]]);
        put(g.pred[e] == NIL ? NIL : dnew[g.pred[e]]);
    }
    s += '|';
    for (int v : vorder) {
        put(g.dlo[v]);
        put(g.dhi[v]);
    }
    return s;
}

}  // namespace apex
