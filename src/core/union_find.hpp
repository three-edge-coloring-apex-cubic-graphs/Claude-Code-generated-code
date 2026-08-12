// Union-find as specified in Section A.3 of [IKM+26]: unite(x, y) makes root(y) the new
// root.  The orientation matters -- Algorithm A.3.1 relies on the second argument becoming
// the representative -- so this deliberately does not use union by rank.
#pragma once

#include <numeric>
#include <vector>

namespace apex {

struct UnionFind {
    std::vector<int> parent;

    explicit UnionFind(int n) : parent(n) { std::iota(parent.begin(), parent.end(), 0); }

    int root(int x) {
        while (parent[x] != x) x = parent[x] = parent[parent[x]];
        return x;
    }
    void unite(int x, int y) { parent[root(x)] = root(y); }
    bool same(int x, int y) { return root(x) == root(y); }
};

}  // namespace apex
