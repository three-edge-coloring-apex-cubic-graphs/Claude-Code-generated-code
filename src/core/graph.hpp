// Dart representation of pseudo-embeddings / pseudo-triangulations (with digons).
//
// A dart e carries head(e), rev(e), succ(e), pred(e).  Following the convention of
// Algorithm A.5.1 / B.2.10, the dart darts[a][b] has head v_a and tail v_b, and succ
// walks the rotation around the head.  succ(e) == NIL / pred(e) == NIL marks the
// boundary of the incidence list of head(e).
#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace apex {

constexpr int NIL = -1;
// The paper writes the degree upper bound "infinity"; the file formats encode it as 0.
constexpr int INF = 1000000;

struct Graph {
    int nv = 0;
    std::vector<int> head, rev, succ, pred;  // one entry per dart
    std::vector<int> dlo, dhi;               // one entry per vertex: [delta^-, delta^+]

    int nd() const { return static_cast<int>(head.size()); }
    int tail(int e) const { return head[rev[e]]; }

    int addVertex(int lo = 1, int hi = INF) {
        dlo.push_back(lo);
        dhi.push_back(hi);
        return nv++;
    }
    int addDart() {
        head.push_back(NIL);
        rev.push_back(NIL);
        succ.push_back(NIL);
        pred.push_back(NIL);
        return nd() - 1;
    }
};

// A homomorphism V(Z) u D(Z) -> V(Z*) u D(Z*), stored as two index tables.
struct Map {
    std::vector<int> v, d;
};

Map identityMap(const Graph& g);
// Returns second o first.
Map compose(const Map& second, const Map& first);

// darts grouped by head vertex, each list in rotation order starting at the first dart
std::vector<std::vector<int>> incidence(const Graph& g);

int dartDegree(const Graph& g, int v);         // d_Z(v)
bool isBoundaryVertex(const Graph& g, int v);  // some dart with head v has succ == NIL
int firstDart(const Graph& g, int v);          // pred == NIL, or NIL if v is inner
int lastDart(const Graph& g, int v);           // succ == NIL, or NIL if v is inner
bool hasLoop(const Graph& g);

// Disjoint union; vertices/darts of b are shifted by a.nv / a.nd().
Graph disjointUnion(const Graph& a, const Graph& b);

// Mirror image (Algorithm A.6.5): swap pred and succ of every dart.
Graph mirror(const Graph& g);

// Canonical serialisation of the connected component of `root`, rooted at dart `root`.
// Two rooted graphs have equal keys iff they are isomorphic by an isomorphism sending
// root to root.  Used to deduplicate sets of combined rules / cartwheels.
std::string canonicalKey(const Graph& g, int root);

}  // namespace apex
