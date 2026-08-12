#include "hom/hom.hpp"

#include <algorithm>
#include <deque>
#include <set>
#include <string>

#include "core/union_find.hpp"

namespace apex {

bool gIntersection(int lo, int hi, int lo2, int hi2) {
    return std::max(lo, lo2) <= std::min(hi, hi2);
}

bool gInclude(int lo, int hi, int lo2, int hi2) { return lo <= lo2 && hi2 <= hi; }

// ---------------------------------------------------------------- Algorithm A.2.1

std::optional<Map> homomorphism(const Graph& z, int e, const Graph& zs, int es,
                                bool (*g)(int, int, int, int)) {
    Map phi;
    phi.v.assign(z.nv, NIL);
    phi.d.assign(z.nd(), NIL);
    std::deque<std::pair<int, int>> q;
    q.emplace_back(e, es);
    while (!q.empty()) {
        auto [f, fs] = q.front();
        q.pop_front();
        if (phi.d[f] != NIL) {
            if (phi.d[f] != fs) return std::nullopt;
            continue;
        }
        phi.d[f] = fs;
        const int h = z.head[f], hs = zs.head[fs];
        if (phi.v[h] != NIL && phi.v[h] != hs) return std::nullopt;
        phi.v[h] = hs;
        if (!g(z.dlo[h], z.dhi[h], zs.dlo[hs], zs.dhi[hs])) return std::nullopt;
        q.emplace_back(z.rev[f], zs.rev[fs]);
        if (z.succ[f] != NIL && zs.succ[fs] == NIL) return std::nullopt;
        if (z.succ[f] != NIL) q.emplace_back(z.succ[f], zs.succ[fs]);
        if (z.pred[f] != NIL && zs.pred[fs] == NIL) return std::nullopt;
        if (z.pred[f] != NIL) q.emplace_back(z.pred[f], zs.pred[fs]);
    }
    return phi;
}

bool homomorphismExists(const Graph& z, int e, const Graph& zs, int es,
                        bool (*g)(int, int, int, int)) {
    // Cheap rejection on the two endpoints of the root pair before touching the scratch space.
    {
        const int h = z.head[e], hs = zs.head[es];
        if (!g(z.dlo[h], z.dhi[h], zs.dlo[hs], zs.dhi[hs])) return false;
        const int t = z.tail(e), ts = zs.tail(es);
        if (!g(z.dlo[t], z.dhi[t], zs.dlo[ts], zs.dhi[ts])) return false;
    }

    // Stamped scratch buffers: a generation counter avoids re-clearing between calls.
    thread_local std::vector<int> imgV, imgD, stampV, stampD;
    thread_local std::vector<std::pair<int, int>> queue;
    thread_local int generation = 0;
    ++generation;
    if (static_cast<int>(imgV.size()) < z.nv) imgV.resize(z.nv);
    if (static_cast<int>(imgD.size()) < z.nd()) imgD.resize(z.nd());
    if (static_cast<int>(stampV.size()) < z.nv) stampV.resize(z.nv, 0);
    if (static_cast<int>(stampD.size()) < z.nd()) stampD.resize(z.nd(), 0);

    queue.clear();
    queue.emplace_back(e, es);
    for (size_t qi = 0; qi < queue.size(); ++qi) {
        const auto [f, fs] = queue[qi];
        if (stampD[f] == generation) {
            if (imgD[f] != fs) return false;
            continue;
        }
        stampD[f] = generation;
        imgD[f] = fs;
        const int h = z.head[f], hs = zs.head[fs];
        if (stampV[h] == generation) {
            if (imgV[h] != hs) return false;
        } else {
            stampV[h] = generation;
            imgV[h] = hs;
            if (!g(z.dlo[h], z.dhi[h], zs.dlo[hs], zs.dhi[hs])) return false;
        }
        queue.emplace_back(z.rev[f], zs.rev[fs]);
        if (z.succ[f] != NIL) {
            if (zs.succ[fs] == NIL) return false;
            queue.emplace_back(z.succ[f], zs.succ[fs]);
        }
        if (z.pred[f] != NIL) {
            if (zs.pred[fs] == NIL) return false;
            queue.emplace_back(z.pred[f], zs.pred[fs]);
        }
    }
    return true;
}

// ---------------------------------------------------------------- Algorithm A.3.1

namespace {

// A UnionFind view over caller-owned storage, so the parent arrays can be reused.  Same
// orientation as core/union_find.hpp: unite(x, y) makes root(y) the representative.
struct ScratchUnionFind {
    std::vector<int>& parent;
    int root(int x) {
        while (parent[x] != x) x = parent[x] = parent[parent[x]];
        return x;
    }
    void unite(int x, int y) { parent[root(x)] = root(y); }
    bool same(int x, int y) { return root(x) == root(y); }
};

std::pair<Graph, Map> freeHomTriangulation(const Graph& z,
                                           const std::vector<std::pair<int, int>>& requests) {
    // Scratch buffers reused across calls: this routine runs tens of millions of times inside
    // the cartwheel enumeration and the vectors dominate its cost otherwise.
    thread_local std::vector<int> head, rev, succ, pred, vnew, dnew, parentV, parentD;
    thread_local std::vector<std::pair<int, int>> q;
    head = z.head;
    rev = z.rev;
    succ = z.succ;
    pred = z.pred;
    parentV.resize(z.nv);
    parentD.resize(z.nd());
    for (int i = 0; i < z.nv; ++i) parentV[i] = i;
    for (int i = 0; i < z.nd(); ++i) parentD[i] = i;
    ScratchUnionFind ufV{parentV}, ufD{parentD};

    q.assign(requests.begin(), requests.end());
    for (size_t qi = 0; qi < q.size(); ++qi) {
        const auto [e, f] = q[qi];
        if (ufD.same(e, f)) continue;
        if (!ufV.same(head[e], head[f])) ufV.unite(head[e], head[f]);
        const int es = ufD.root(e), fs = ufD.root(f);
        ufD.unite(es, fs);  // fs is now the representative
        q.emplace_back(rev[es], rev[fs]);
        if (succ[es] != NIL && succ[fs] != NIL) q.emplace_back(succ[es], succ[fs]);
        if (pred[es] != NIL && pred[fs] != NIL) q.emplace_back(pred[es], pred[fs]);
        if (succ[es] != NIL && succ[fs] == NIL) succ[fs] = succ[es];
        if (pred[es] != NIL && pred[fs] == NIL) pred[fs] = pred[es];
    }

    vnew.assign(z.nv, NIL);
    dnew.assign(z.nd(), NIL);
    int nv = 0, nd = 0;
    for (int v = 0; v < z.nv; ++v)
        if (ufV.root(v) == v) vnew[v] = nv++;
    for (int d = 0; d < z.nd(); ++d)
        if (ufD.root(d) == d) dnew[d] = nd++;

    Graph g;
    g.nv = nv;
    g.dlo.assign(nv, 1);
    g.dhi.assign(nv, INF);
    g.head.resize(nd);
    g.rev.resize(nd);
    g.succ.resize(nd);
    g.pred.resize(nd);
    auto mapDart = [&](int d) { return d == NIL ? NIL : dnew[ufD.root(d)]; };
    for (int d = 0; d < z.nd(); ++d) {
        if (ufD.root(d) != d) continue;
        const int nd2 = dnew[d];
        g.head[nd2] = vnew[ufV.root(head[d])];
        g.rev[nd2] = mapDart(rev[d]);
        g.succ[nd2] = mapDart(succ[d]);
        g.pred[nd2] = mapDart(pred[d]);
    }

    Map phi;
    phi.v.resize(z.nv);
    phi.d.resize(z.nd());
    for (int v = 0; v < z.nv; ++v) phi.v[v] = vnew[ufV.root(v)];
    for (int d = 0; d < z.nd(); ++d) phi.d[d] = dnew[ufD.root(d)];
    return {g, phi};
}

}  // namespace

// ---------------------------------------------------------------- Algorithm A.4.1

std::optional<HomImage> dartIdentification(const Graph& z,
                                           const std::vector<std::pair<int, int>>& requests) {
    auto [g, phi] = freeHomTriangulation(z, requests);
    if (hasLoop(g)) return std::nullopt;  // a loop error
    for (int v = 0; v < z.nv; ++v) {
        const int vs = phi.v[v];
        const int lo = std::max(g.dlo[vs], z.dlo[v]);
        const int hi = std::min(g.dhi[vs], z.dhi[v]);
        if (lo > hi) return std::nullopt;  // a degree-mismatch error
        g.dlo[vs] = lo;
        g.dhi[vs] = hi;
    }
    return HomImage{g, phi};
}

// ------------------------------------------- Algorithms B.2.4, B.2.6 (shared building blocks)

void addBoundaryDartsDirectly(Graph& z, int efirst, int elast) {
    const int u = z.tail(efirst), w = z.tail(elast);
    const int revFirst = z.rev[efirst], revLast = z.rev[elast];
    const int f = z.addDart();
    const int gd = z.addDart();
    z.head[f] = u;
    z.rev[f] = gd;
    z.succ[f] = NIL;
    z.pred[f] = revFirst;
    z.head[gd] = w;
    z.rev[gd] = f;
    z.succ[gd] = revLast;
    z.pred[gd] = NIL;
    z.succ[revFirst] = f;
    z.pred[revLast] = gd;
}

Map linkIncidenceListEnds(Graph& z, int eufirst, int ewlast) {
    const int u = z.head[eufirst];
    const int w = z.head[ewlast];
    z.pred[eufirst] = ewlast;
    z.succ[ewlast] = eufirst;
    if (u == w) return identityMap(z);

    for (int e = 0; e < z.nd(); ++e)
        if (z.head[e] == u) z.head[e] = w;
    z.dlo[w] = std::max(z.dlo[w], z.dlo[u]);
    z.dhi[w] = std::min(z.dhi[w], z.dhi[u]);

    // Remove u from the vertex set, renumbering the survivors.
    std::vector<int> vnew(z.nv, NIL);
    int id = 0;
    std::vector<int> lo, hi;
    for (int v = 0; v < z.nv; ++v) {
        if (v == u) continue;
        vnew[v] = id++;
        lo.push_back(z.dlo[v]);
        hi.push_back(z.dhi[v]);
    }
    vnew[u] = vnew[w];
    for (int e = 0; e < z.nd(); ++e) z.head[e] = vnew[z.head[e]];
    z.nv = id;
    z.dlo = lo;
    z.dhi = hi;

    Map phi;
    phi.v = vnew;
    phi.d.resize(z.nd());
    for (int e = 0; e < z.nd(); ++e) phi.d[e] = e;
    return phi;
}

// --------------------------------------- Algorithms B.2.1, B.2.2, B.2.3, B.2.5

namespace {

// Algorithm B.2.3 (= Algorithm A.4.8).
std::optional<Graph> addBoundaryDarts(const Graph& z, int v) {
    const int ef = firstDart(z, v), el = lastDart(z, v);
    if (z.tail(ef) == z.tail(el)) return std::nullopt;
    Graph g = z;
    g.pred[ef] = el;
    g.succ[el] = ef;
    addBoundaryDartsDirectly(g, ef, el);
    return g;
}

// Algorithm B.2.5.
std::optional<HomImage> identifyNeighbors(const Graph& z, int v) {
    const int ef = firstDart(z, v), el = lastDart(z, v);
    const int tf = z.tail(ef), tl = z.tail(el);
    if (!gIntersection(z.dlo[tf], z.dhi[tf], z.dlo[tl], z.dhi[tl])) return std::nullopt;
    Graph g = z;
    g.pred[ef] = el;
    g.succ[el] = ef;
    Map phi = linkIncidenceListEnds(g, g.rev[el], g.rev[ef]);
    if (hasLoop(g)) return std::nullopt;
    return HomImage{g, phi};
}

std::vector<HomImage> boundaryCompletions(const Graph& z, int v, Boundary policy) {
    std::vector<HomImage> out;
    switch (policy) {
        case Boundary::PseudoEmbedding: {  // Algorithm B.2.1
            Graph g = z;
            const int ef = firstDart(z, v), el = lastDart(z, v);
            g.pred[ef] = el;
            g.succ[el] = ef;
            out.push_back({g, identityMap(z)});
            break;
        }
        case Boundary::PseudoTriangulation: {  // Algorithm B.2.3
            if (auto g = addBoundaryDarts(z, v)) out.push_back({*g, identityMap(z)});
            break;
        }
        case Boundary::PseudoTriangulationWithDigons: {  // Algorithm B.2.2
            if (auto g = addBoundaryDarts(z, v)) out.push_back({*g, identityMap(z)});
            if (auto r = identifyNeighbors(z, v)) out.push_back(*r);
            break;
        }
    }
    return out;
}

// Algorithm A.4.5.
bool innerSubdegreeError(const Graph& z) {
    std::vector<int> deg(z.nv, 0);
    std::vector<char> boundary(z.nv, 0);
    for (int e = 0; e < z.nd(); ++e) {
        ++deg[z.head[e]];
        if (z.succ[e] == NIL) boundary[z.head[e]] = 1;
    }
    for (int v = 0; v < z.nv; ++v)
        if (!boundary[v] && deg[v] < z.dlo[v]) return true;
    return false;
}

// Algorithm A.4.6.
int vertexSingleDegreeIssue(const Graph& z) {
    std::vector<int> deg(z.nv, 0);
    std::vector<char> boundary(z.nv, 0);
    for (int e = 0; e < z.nd(); ++e) {
        ++deg[z.head[e]];
        if (z.succ[e] == NIL) boundary[z.head[e]] = 1;
    }
    for (int v = 0; v < z.nv; ++v) {
        if (z.dlo[v] != z.dhi[v]) continue;
        if (z.dlo[v] < deg[v]) return v;
        if (boundary[v] && deg[v] == z.dlo[v]) return v;
    }
    return NIL;
}

// Algorithm A.4.7, generalised to return a set (Algorithm B.2.2 may offer two completions).
std::vector<HomImage> fixSingleDegreeIssue(const Graph& z, int v, Boundary policy) {
    const int deg = dartDegree(z, v);
    if (z.dlo[v] < deg) {
        int e = isBoundaryVertex(z, v) ? firstDart(z, v) : NIL;
        if (e == NIL)
            for (int d = 0; d < z.nd(); ++d)
                if (z.head[d] == v) {
                    e = d;
                    break;
                }
        int f = e;
        for (int i = 0; i < z.dlo[v]; ++i) f = z.succ[f];
        if (auto a = dartIdentification(z, {{e, f}})) return {*a};
        return {};
    }
    if (isBoundaryVertex(z, v) && z.dlo[v] == deg) return boundaryCompletions(z, v, policy);
    return {};
}

// Algorithm A.4.9.
std::optional<std::pair<Graph, Graph>> singleOutLowerDegree(const Graph& z) {
    for (int v = 0; v < z.nv; ++v) {
        if (z.dlo[v] < z.dhi[v] && z.dlo[v] <= dartDegree(z, v)) {
            Graph a = z, b = z;
            a.dhi[v] = z.dlo[v];
            b.dlo[v] = z.dlo[v] + 1;
            return std::make_pair(a, b);
        }
    }
    return std::nullopt;
}

std::string stateKey(const HomImage& h) {
    std::string s;
    auto put = [&](int x) {
        s += std::to_string(x);
        s += ',';
    };
    put(h.g.nv);
    put(h.g.nd());
    for (int e = 0; e < h.g.nd(); ++e) {
        put(h.g.head[e]);
        put(h.g.rev[e]);
        put(h.g.succ[e]);
        put(h.g.pred[e]);
    }
    for (int v = 0; v < h.g.nv; ++v) {
        put(h.g.dlo[v]);
        put(h.g.dhi[v]);
    }
    s += '|';
    for (int x : h.phi.v) put(x);
    s += '|';
    for (int x : h.phi.d) put(x);
    return s;
}

// Algorithm A.4.4.
std::vector<HomImage> resolveDegreeIssues(const Graph& z, Boundary policy) {
    std::vector<HomImage> out;
    std::set<std::string> seen;
    std::deque<HomImage> q;
    q.push_back({z, identityMap(z)});
    while (!q.empty()) {
        HomImage cur = std::move(q.front());
        q.pop_front();
        if (innerSubdegreeError(cur.g)) continue;
        const int v = vertexSingleDegreeIssue(cur.g);
        if (v != NIL) {
            for (auto& a : fixSingleDegreeIssue(cur.g, v, policy))
                q.push_back({a.g, compose(a.phi, cur.phi)});
            continue;
        }
        if (auto b = singleOutLowerDegree(cur.g)) {
            q.push_back({b->first, cur.phi});
            q.push_back({b->second, cur.phi});
            continue;
        }
        if (seen.insert(stateKey(cur)).second) out.push_back(std::move(cur));
    }
    return out;
}

}  // namespace

// ---------------------------------------------------------------- Algorithm A.4.3

std::vector<HomImage> freeHomomorphism(const Graph& z,
                                       const std::vector<std::pair<int, int>>& requests,
                                       Boundary policy) {
    auto a = dartIdentification(z, requests);
    if (!a) return {};
    std::vector<HomImage> out;
    for (auto& r : resolveDegreeIssues(a->g, policy))
        out.push_back({std::move(r.g), compose(r.phi, a->phi)});
    return out;
}

// ---------------------------------------------------------------- Algorithms B.2.12, B.3.6

std::optional<std::pair<int, int>> twoDigonsIncidentWithSameVertex(const Graph& z) {
    const auto inc = incidence(z);
    for (int v = 0; v < z.nv; ++v) {
        int digonEdge = NIL;
        for (int e0 : inc[v]) {
            if (z.succ[e0] == NIL) continue;
            const int e1 = z.rev[z.succ[e0]];
            if (z.succ[e1] == NIL) continue;
            const int e2 = z.rev[z.succ[e1]];
            if (e0 != e2) continue;
            if (digonEdge == NIL)
                digonEdge = e0;
            else
                return std::make_pair(digonEdge, e0);
        }
    }
    return std::nullopt;
}

std::vector<std::pair<int, int>> enumDigons(const Graph& z) {
    std::set<std::pair<int, int>> d;
    for (int e0 = 0; e0 < z.nd(); ++e0) {
        if (z.succ[e0] == NIL) continue;
        const int e1 = z.rev[z.succ[e0]];
        if (z.succ[e1] == NIL) continue;
        const int e2 = z.rev[z.succ[e1]];
        if (e2 != e0) continue;
        int u1 = z.head[e0], u2 = z.head[z.rev[e0]];
        if (u1 > u2) std::swap(u1, u2);
        d.insert({u1, u2});
    }
    return {d.begin(), d.end()};
}

// ---------------------------------------------------------------- Algorithms B.2.11, B.2.13

std::vector<HomImage> enforceSingleDigonIncidence(const Graph& z, Boundary policy) {
    std::vector<HomImage> out;
    std::deque<HomImage> q;
    q.push_back({z, identityMap(z)});
    while (!q.empty()) {
        HomImage cur = std::move(q.front());
        q.pop_front();
        auto c = twoDigonsIncidentWithSameVertex(cur.g);
        if (c) {
            for (auto& r : freeHomomorphism(cur.g, {{c->first, c->second}}, policy))
                q.push_back({std::move(r.g), compose(r.phi, cur.phi)});
            continue;
        }
        out.push_back(std::move(cur));
    }
    return out;
}

bool freeHomomorphismAndEnforceSingleDigonIncidenceAny(
    const Graph& z, const std::vector<std::pair<int, int>>& requests, Boundary policy) {
    // Algorithms A.4.3, A.4.4 and B.2.11 fused into one worklist, carrying graphs only and
    // returning at the first image.  Whenever a degree-issue-free graph still has two digons at
    // a common vertex, Algorithm B.2.11 asks for the free homomorphic images respecting that
    // pair -- which is dartIdentification followed by resolveDegreeIssues again, i.e. exactly
    // what pushing the identified graph back onto this worklist does.
    auto a = dartIdentification(z, requests);
    if (!a) return false;

    std::vector<Graph> stack;
    stack.push_back(std::move(a->g));
    while (!stack.empty()) {
        Graph cur = std::move(stack.back());
        stack.pop_back();
        if (innerSubdegreeError(cur)) continue;
        const int v = vertexSingleDegreeIssue(cur);
        if (v != NIL) {
            for (auto& r : fixSingleDegreeIssue(cur, v, policy)) stack.push_back(std::move(r.g));
            continue;
        }
        if (auto b = singleOutLowerDegree(cur)) {
            stack.push_back(std::move(b->first));
            stack.push_back(std::move(b->second));
            continue;
        }
        const auto c = twoDigonsIncidentWithSameVertex(cur);
        if (!c) return true;
        if (auto d = dartIdentification(cur, {{c->first, c->second}}))
            stack.push_back(std::move(d->g));
    }
    return false;
}

std::vector<HomImage> freeHomomorphismAndEnforceSingleDigonIncidence(
    const Graph& z, const std::vector<std::pair<int, int>>& requests, Boundary policy) {
    std::vector<HomImage> out;
    for (auto& a : freeHomomorphism(z, requests, policy))
        for (auto& b : enforceSingleDigonIncidence(a.g, policy))
            out.push_back({std::move(b.g), compose(b.phi, a.phi)});
    return out;
}

}  // namespace apex
