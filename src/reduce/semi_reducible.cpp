#include "reduce/semi_reducible.hpp"

#include <algorithm>
#include <cstdint>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include "reduce/coloring.hpp"
#include "reduce/deletable.hpp"
#include "reduce/kempe.hpp"

namespace apex {
namespace {

// 3^24 ~= 2.8e11 would already be far too large; this cap keeps the dense array within a few
// hundred MB. Revisit once real ring-edge-total statistics are in (see HANDOFF.md).
constexpr uint64_t kMaxDenseSize = 60'000'000;

struct RingLayout {
    std::vector<int> offset;  // per boundary component, global start position in phi
    std::vector<int> size;    // per boundary component, r_i
    int total = 0;            // |R|
};

RingLayout makeLayout(const MultiBoundaryIsland& island) {
    RingLayout rl;
    int off = 0;
    for (int r : island.ringSizes) {
        rl.offset.push_back(off);
        rl.size.push_back(r);
        off += r;
    }
    rl.total = off;
    return rl;
}

std::vector<uint64_t> pow3Table(int m) {
    std::vector<uint64_t> p(m + 1);
    p[0] = 1;
    for (int i = 1; i <= m; ++i) p[i] = p[i - 1] * 3;
    return p;
}

uint64_t codeOf(const std::vector<int>& phi, const std::vector<uint64_t>& pow3) {
    uint64_t c = 0;
    for (size_t i = 0; i < phi.size(); ++i) c += static_cast<uint64_t>(phi[i]) * pow3[i];
    return c;
}

std::vector<int> decode(uint64_t code, int m) {
    std::vector<int> phi(m);
    for (int i = 0; i < m; ++i) {
        phi[i] = static_cast<int>(code % 3);
        code /= 3;
    }
    return phi;
}

// A unit of a semi-matching M: either a match (2 global positions, flipped together when the
// unit is chosen) or a singleton (1 global position, flipped alone). Definition 3.3 makes
// singletons full members of M, not inert padding: Definition 3.4's M' ranges over subsets of
// ALL of M's elements (matches and singletons alike), and a singleton {r} included in M' still
// has r in supp(M'), so its own colour flips too. A semi-matching with t matches and s
// singletons therefore has 2^(t+s) reachable colourings, not 2^t.
using Unit = std::vector<int>;

// Candidate semi-matchings (Definition 3.3) of the positions of one boundary component coloured
// {x, y} under phi, each given as its full list of units (matches + singletons), in global
// position indices. Per B.1.1/B.1.2's own remark, semi-matchings for a multi-boundary island are
// the concatenation of per-component semi-matchings, and the per-component match candidates are
// exactly GetPlanarHalfKempes applied to those positions, taken in their cyclic order within the
// ring (a chord's crossing status depends only on relative cyclic order, so relabelling the
// x/y-coloured subset 0..k-1 and mapping back is exact); the positions GetPlanarHalfKempes
// leaves unmatched become this candidate's singleton units.
//
// This result depends only on (offset, size, which local positions are coloured {x, y}) -- not
// on which of x or y each one is, nor on which island/phi produced the pattern. Across the huge
// number of (phi, x, y) triples checked in the fixpoint (millions, for the larger ring sizes in
// I), the same local pattern recurs constantly -- for a ring of size r there are only 2^r
// possible patterns, orders of magnitude fewer than 3^r * 3 calls -- so this is cached globally,
// keyed on (offset, local bitmask), guarded by a mutex since checkReducibility runs concurrently
// across islands (see verify_b3_reduce.cpp). Measured live on one of the largest islands in I
// (ring 13, ~1.5M active colourings in the first fixpoint round): rebuilding and sorting this
// list from scratch on every call was the dominant cost of a multi-hour-per-island runtime.
std::vector<std::shared_ptr<const std::vector<Unit>>> componentCandidatesCached(int offset,
                                                                                int size,
                                                                                uint32_t bitmask);

std::vector<std::vector<Unit>> componentCandidatesUncached(int offset,
                                                            const std::vector<int>& local) {
    std::vector<std::vector<Unit>> out;
    if (local.empty()) {
        out.push_back({});
        return out;
    }
    if (local.size() == 1) {
        out.push_back({Unit{offset + local[0]}});
        return out;
    }
    const auto& base = getPlanarHalfKempes(static_cast<int>(local.size()));
    out.reserve(base.size());
    for (const auto& m0 : base) {
        std::vector<char> matched(local.size(), 0);
        std::vector<Unit> units;
        units.reserve(local.size());
        for (auto& [a, b] : m0) {
            units.push_back({offset + local[a], offset + local[b]});
            matched[a] = matched[b] = 1;
        }
        for (size_t i = 0; i < local.size(); ++i)
            if (!matched[i]) units.push_back({offset + local[static_cast<int>(i)]});
        out.push_back(std::move(units));
    }
    // Fewer units means fewer switch combinations to check (2^|units|) and, since more matches
    // (vs. singletons) among the same local positions means fewer units, tends to be more likely
    // to survive too — trying those candidates first lets both outcomes (a witness found, or
    // every candidate exhausted) surface with less wasted work on the worst offenders.
    std::sort(out.begin(), out.end(), [](const std::vector<Unit>& a, const std::vector<Unit>& b) {
        return a.size() < b.size();
    });
    return out;
}

std::vector<std::shared_ptr<const std::vector<Unit>>> componentCandidatesCached(int offset,
                                                                                int size,
                                                                                uint32_t bitmask) {
    static std::mutex mu;
    static std::unordered_map<uint64_t,
                              std::shared_ptr<const std::vector<std::vector<Unit>>>>
        cache;
    const uint64_t key = (static_cast<uint64_t>(offset) << 40) ^
                         (static_cast<uint64_t>(size) << 32) ^ bitmask;

    std::shared_ptr<const std::vector<std::vector<Unit>>> full;
    {
        std::lock_guard<std::mutex> lock(mu);
        auto it = cache.find(key);
        if (it != cache.end()) full = it->second;
    }
    if (!full) {
        std::vector<int> local;
        for (int p = 0; p < size; ++p)
            if ((bitmask >> p) & 1u) local.push_back(p);
        auto computed = std::make_shared<std::vector<std::vector<Unit>>>(
            componentCandidatesUncached(offset, local));
        std::lock_guard<std::mutex> lock(mu);
        full = cache.emplace(key, std::move(computed)).first->second;
    }

    std::vector<std::shared_ptr<const std::vector<Unit>>> out;
    out.reserve(full->size());
    for (const auto& units : *full)
        out.push_back(std::shared_ptr<const std::vector<Unit>>(full, &units));
    return out;
}

std::vector<std::shared_ptr<const std::vector<Unit>>> componentCandidates(
    const std::vector<int>& phi, int offset, int size, int x, int y) {
    uint32_t bitmask = 0;
    for (int p = 0; p < size; ++p) {
        const int c = phi[offset + p];
        if (c == x || c == y) bitmask |= (1u << p);
    }
    return componentCandidatesCached(offset, size, bitmask);
}

// Whether every Kempe-switch (Definition 3.4: every subset M' of `units`, each unit flipping the
// colours x/y at all of its global positions together) of phi (code `baseCode`) stays inside
// `inS`. Visits the 2^|units| subsets in Gray-code order (each step toggles exactly one unit) so
// `code` is updated incrementally in O(1) amortised instead of rescanning all units per subset —
// this loop dominates checkReducibility's cost, so the difference matters at the ring sizes (up
// to 15) actually occurring in I.
bool allSwitchesSurvive(const std::vector<int>& phi, uint64_t baseCode,
                        const std::vector<uint8_t>& inS, const std::vector<uint64_t>& pow3,
                        const std::vector<Unit>& units, int x, int y) {
    const int t = static_cast<int>(units.size());
    std::vector<int64_t> delta(t);
    for (int i = 0; i < t; ++i) {
        int64_t d = 0;
        for (int p : units[i]) {
            const int cp = phi[p];
            d += static_cast<int64_t>((x + y - cp) - cp) * static_cast<int64_t>(pow3[p]);
        }
        delta[i] = d;
    }

    int64_t code = static_cast<int64_t>(baseCode);
    if (!inS[static_cast<uint64_t>(code)]) return false;
    const uint32_t total = 1u << t;
    for (uint32_t i = 1; i < total; ++i) {
        const int bit = __builtin_ctz(i);
        const uint32_t gray = i ^ (i >> 1);
        if ((gray >> bit) & 1u)
            code += delta[bit];
        else
            code -= delta[bit];
        if (!inS[static_cast<uint64_t>(code)]) return false;
    }
    return true;
}

// Whether phi has, for every colour pair {x, y}, a witnessing combined semi-matching (Cartesian
// product of per-component candidates) all of whose Kempe-switches stay inside the current S.
bool survives(const std::vector<int>& phi, uint64_t code, const std::vector<uint8_t>& inS,
             const std::vector<uint64_t>& pow3, const RingLayout& rl) {
    const int k = static_cast<int>(rl.size.size());
    for (int x = 0; x < 3; ++x) {
        for (int y = x + 1; y < 3; ++y) {
            std::vector<std::vector<std::shared_ptr<const std::vector<Unit>>>> perComponent(k);
            for (int i = 0; i < k; ++i)
                perComponent[i] = componentCandidates(phi, rl.offset[i], rl.size[i], x, y);

            bool foundWitness = false;
            std::vector<int> idx(k, 0);
            while (true) {
                std::vector<Unit> combined;
                for (int i = 0; i < k; ++i)
                    for (auto& u : *perComponent[i][idx[i]]) combined.push_back(u);
                if (allSwitchesSurvive(phi, code, inS, pow3, combined, x, y)) {
                    foundWitness = true;
                    break;
                }
                int i = k - 1;
                while (i >= 0) {
                    if (++idx[i] < static_cast<int>(perComponent[i].size())) break;
                    idx[i] = 0;
                    --i;
                }
                if (i < 0) break;  // exhausted every combination, k == 0 included
            }
            if (!foundWitness) return false;
        }
    }
    return true;
}

}  // namespace

ReducibilityResult checkReducibility(const MultiBoundaryIsland& island) {
    ReducibilityResult result;
    const RingLayout rl = makeLayout(island);
    const int m = rl.total;
    result.ringEdgeTotal = m;
    const auto pow3 = pow3Table(m);
    const uint64_t total = pow3[m];
    if (total > kMaxDenseSize)
        throw std::runtime_error("checkReducibility: ring space too large to enumerate densely");

    std::vector<uint8_t> inS(total, 1);
    for (const auto& phi : computeCI(island)) inS[codeOf(phi, pow3)] = 0;

    std::vector<uint64_t> active;
    active.reserve(total);
    for (uint64_t c = 0; c < total; ++c)
        if (inS[c]) active.push_back(c);

    bool changed = true;
    while (changed) {
        changed = false;
        std::vector<uint64_t> stillActive;
        stillActive.reserve(active.size());
        for (uint64_t code : active) {
            const auto phi = decode(code, m);
            if (survives(phi, code, inS, pow3, rl)) {
                stillActive.push_back(code);
            } else {
                inS[code] = 0;
                changed = true;
            }
        }
        active.swap(stillActive);
    }

    result.fixpointSize = active.size();
    result.semiDReducible = active.empty();
    if (!result.semiDReducible) {
        for (const auto& f : getDeletableEdgeSet(island)) {
            bool disjoint = true;
            for (const auto& phi : computeCIModuloF(island, f)) {
                if (inS[codeOf(phi, pow3)]) {
                    disjoint = false;
                    break;
                }
            }
            if (disjoint) {
                result.semiCReducible = true;
                result.witnessF = f;
                break;
            }
        }
    }
    return result;
}

}  // namespace apex
