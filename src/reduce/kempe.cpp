#include "reduce/kempe.hpp"

#include <algorithm>
#include <map>
#include <set>

namespace apex {
namespace {

Matching normalized(Matching m) {
    for (auto& [a, b] : m)
        if (a > b) std::swap(a, b);
    std::sort(m.begin(), m.end());
    return m;
}

Matching shifted(const Matching& m, int c) {
    Matching r = m;
    for (auto& [a, b] : r) {
        a += c;
        b += c;
    }
    return r;
}

}  // namespace

// Algorithm B.1.1.
const std::vector<Matching>& getValidParens(int q) {
    static std::map<int, std::vector<Matching>> memo;
    if (auto it = memo.find(q); it != memo.end()) return it->second;

    std::vector<Matching>& p = memo[q];
    if (q == 1) {
        p.push_back({{0, 1}});
        return p;
    }
    // The paper writes every union here as a set union (P <- P u {...}); the same matching can
    // be produced by more than one split i (e.g. {(0,1),(2,3),(4,5)} respects both i=1 and i=2),
    // so this recursion needs explicit deduplication, not a plain append.
    std::set<Matching> seen;
    for (int i = 1; i <= q - 1; ++i) {
        const auto& ls = getValidParens(i);
        const auto& rs = getValidParens(q - i);
        for (const Matching& l : ls) {
            for (const Matching& r : rs) {
                Matching combined = l;
                const Matching rShifted = shifted(r, 2 * i);
                combined.insert(combined.end(), rShifted.begin(), rShifted.end());
                if (seen.insert(normalized(combined)).second) p.push_back(std::move(combined));
            }
        }
    }
    for (const Matching& m : getValidParens(q - 1)) {
        // Paper: (M + 1) union {{1, 2q}}, 1-indexed; 0-indexed that is (M+1) union {(0, 2q-1)}.
        // pdftotext rendered "2q" as a literal "2" here (the same class of subscript-dropping
        // error as elsewhere in this document); {1,2} would double-cover position 2 with M+1's
        // shifted range and leave position 2q uncovered, so it cannot be read literally.
        Matching combined = shifted(m, 1);
        combined.emplace_back(0, 2 * q - 1);
        if (seen.insert(normalized(combined)).second) p.push_back(std::move(combined));
    }
    return p;
}

// Algorithm B.1.2.
const std::vector<Matching>& getPlanarHalfKempes(int n) {
    static std::map<int, std::vector<Matching>> memo;
    if (auto it = memo.find(n); it != memo.end()) return it->second;

    std::vector<Matching>& r = memo[n];
    if (n == 1) {
        r.push_back({});
        return r;
    }

    std::set<Matching> a;  // all generated partial matchings, including redundant ones.
    std::vector<int> bIdx(n);
    for (int i = 0; i < n; ++i) bIdx[i] = i;

    for (int q = n / 2; q >= 1; --q) {
        std::vector<int> comb(2 * q);
        for (int i = 0; i < 2 * q; ++i) comb[i] = i;
        std::vector<char> inB(n, 0);
        while (true) {
            for (int i = 0; i < n; ++i) inB[i] = 0;
            for (int i : comb) inB[i] = 1;

            bool skip = false;
            for (int i = 0; i < n && !skip; ++i) {
                if (!inB[i] && !inB[(i + 1) % n]) skip = true;
            }
            if (!skip) {
                const std::vector<int> b = comb;  // already ascending: b_1 < ... < b_2q
                for (const Matching& m0 : getValidParens(q)) {
                    Matching m;
                    m.reserve(q);
                    for (const auto& [x, y] : m0) m.emplace_back(b[x], b[y]);
                    m = normalized(m);
                    a.insert(m);

                    bool redundant = false;
                    std::vector<char> matched(n, 0);
                    for (const auto& [x, y] : m) {
                        matched[x] = 1;
                        matched[y] = 1;
                    }
                    std::vector<int> un;
                    for (int i = 0; i < n; ++i)
                        if (!matched[i]) un.push_back(i);
                    for (size_t i = 0; i < un.size() && !redundant; ++i) {
                        for (size_t j = i + 1; j < un.size() && !redundant; ++j) {
                            Matching m1 = m;
                            m1.emplace_back(un[i], un[j]);
                            m1 = normalized(m1);
                            if (a.count(m1)) redundant = true;
                        }
                    }
                    if (!redundant) r.push_back(m);
                }
            }

            // Advance comb to the next 2q-subset of [0, n) in ascending combinatorial order.
            int k = 2 * q - 1;
            while (k >= 0 && comb[k] == n - (2 * q - k)) --k;
            if (k < 0) break;
            ++comb[k];
            for (int j = k + 1; j < 2 * q; ++j) comb[j] = comb[j - 1] + 1;
        }
    }
    return r;
}

void warmKempeCaches(int maxN) {
    for (int n = 1; n <= maxN; ++n) getPlanarHalfKempes(n);
}

}  // namespace apex
