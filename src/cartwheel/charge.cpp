#include "cartwheel/charge.hpp"

#include <algorithm>
#include <array>
#include <deque>
#include <mutex>
#include <numeric>

#include "hom/hom.hpp"

namespace apex {
namespace {

using Range = std::pair<int, int>;

// The tables of Algorithm B.3.5.  Every "[5, 8]"-looking entry here is really [5, infinity]:
// pdftotext renders the infinity sign as the digit 8, and the typeset page shows [5, oo].
// Only [7, 10] and [6, 6] / [6, 7] / [5, 5] are genuine finite bounds.
const std::array<Range, 5> kDistPairs{Range{0, 1}, Range{1, 0}, Range{1, 1}, Range{1, 2},
                                      Range{2, 1}};
const std::vector<std::vector<std::pair<Range, Range>>> kDegreePairs{
    {{{7, 10}, {5, INF}}},
    {{{5, INF}, {7, 10}}},
    {{{5, 5}, {5, 5}}, {{5, 5}, {6, 6}}, {{6, 6}, {5, 5}}, {{5, INF}, {5, INF}}},
    {{{5, 5}, {5, 5}}, {{5, 5}, {6, 7}}, {{6, 7}, {5, 5}}},
    {{{5, 5}, {5, 5}}, {{5, 5}, {6, 7}}, {{6, 7}, {5, 5}}},
};
const std::vector<std::vector<int>> kCharges{{4}, {4}, {5, 3, 3, 2}, {4, 2, 2}, {4, 2, 2}};

bool included(int lo, int hi, const Range& d) { return d.first <= lo && hi <= d.second; }

std::vector<int> distancesFrom(const Graph& g, int src) {
    std::vector<std::vector<int>> adj(g.nv);
    for (int e = 0; e < g.nd(); ++e) adj[g.tail(e)].push_back(g.head[e]);
    std::vector<int> dist(g.nv, -1);
    std::deque<int> q{src};
    dist[src] = 0;
    while (!q.empty()) {
        const int u = q.front();
        q.pop_front();
        for (const int w : adj[u])
            if (dist[w] < 0) {
                dist[w] = dist[u] + 1;
                q.push_back(w);
            }
    }
    return dist;
}

thread_local Counters tlCounters;
std::mutex counterMutex;
Counters globalCounters;
struct CounterFlusher {
    ~CounterFlusher() {
        std::lock_guard<std::mutex> lock(counterMutex);
        globalCounters.prunes += tlCounters.prunes;
        globalCounters.alwaysApply += tlCounters.alwaysApply;
        globalCounters.neverApply += tlCounters.neverApply;
        globalCounters.neverApplyDeep += tlCounters.neverApplyDeep;
        globalCounters.blocked += tlCounters.blocked;
        globalCounters.digonCharge += tlCounters.digonCharge;
        globalCounters.prunedByRule += tlCounters.prunedByRule;
        globalCounters.prunedByChargeEarly += tlCounters.prunedByChargeEarly;
        globalCounters.prunedByBlocked += tlCounters.prunedByBlocked;
        globalCounters.prunedByCharge += tlCounters.prunedByCharge;
    }
};
thread_local CounterFlusher tlFlusher;
// The flusher is only instantiated if it is referenced, so touch it on every prune.
void touchFlusher() { (void)&tlFlusher; }

}  // namespace

Counters takeCounters() {
    std::lock_guard<std::mutex> lock(counterMutex);
    Counters c = globalCounters;
    globalCounters = Counters{};
    return c;
}

void Context::finalise() {
    byChargeDesc.resize(combined.size());
    std::iota(byChargeDesc.begin(), byChargeDesc.end(), 0);
    std::sort(byChargeDesc.begin(), byChargeDesc.end(),
              [&](int a, int b) { return combined[a].r > combined[b].r; });
    maxCharge = 0;
    for (const auto& c : combined) maxCharge = std::max(maxCharge, c.r);
}

bool alwaysApply(const Graph& z, int e, const Graph& rg, int rdart) {
    ++tlCounters.alwaysApply;
    return homomorphismExists(rg, rdart, z, e, gInclude);
}

bool neverApply(const Graph& z, int e, const Graph& rg, int rdart) {
    ++tlCounters.neverApply;
    // Identifying e with rdart forces head(e) onto head(rdart) and tail(e) onto tail(rdart), so
    // disjoint degree ranges at either end make dartIdentification -- and hence the whole set of
    // Algorithm B.3.4 -- empty.  Checking that first skips most of the work.
    {
        const int h = z.head[e], rh = rg.head[rdart];
        if (!gIntersection(z.dlo[h], z.dhi[h], rg.dlo[rh], rg.dhi[rh])) return true;
        const int t = z.tail(e), rt = rg.tail(rdart);
        if (!gIntersection(z.dlo[t], z.dhi[t], rg.dlo[rt], rg.dhi[rt])) return true;
        ++tlCounters.neverApplyDeep;
    }
    const Graph joined = disjointUnion(z, rg);
    return !freeHomomorphismAndEnforceSingleDigonIncidenceAny(
        joined, {{e, rdart + z.nd()}}, Boundary::PseudoTriangulationWithDigons);
}

int amountOfChargeSend(const Graph& z, int e, const std::vector<Rule>& rules) {
    int a = 0;
    for (const auto& r : rules)
        if (alwaysApply(z, e, r.g, r.dart)) a += r.r;
    return a;
}

int amountOfPossibleChargeSend(const Graph& z, int e, const Context& ctx) {
    return amountOfPossibleChargeSendAbove(z, e, ctx, 0);
}

int amountOfPossibleChargeSendAbove(const Graph& z, int e, const Context& ctx, int threshold) {
    // Algorithm A.9.4 scans every combined rule and keeps the maximum charge among those that
    // do not never-apply.  Scanning in order of decreasing charge lets the first hit be the
    // answer; and when the caller only needs to distinguish "> threshold" from "<= threshold",
    // rules of charge at most `threshold` cannot change that verdict and are skipped.  With
    // threshold = 0 this returns exactly the value of Algorithm A.9.4, since the running
    // maximum there starts at 0.
    for (const int i : ctx.byChargeDesc) {
        const auto& c = ctx.combined[i];
        if (c.r <= threshold) break;
        if (!neverApply(z, e, c.g, c.dart)) return c.r;
    }
    return threshold < 0 ? threshold : 0;
}

int lowerBoundOfDigonCharge(const Cartwheel& c) {
    ++tlCounters.digonCharge;
    const auto dist = distancesFrom(c.g, c.center);
    int charge = 0;
    for (const auto& [u1, u2] : enumDigons(c.g)) {
        for (size_t i = 0; i < kDistPairs.size(); ++i) {
            if (dist[u1] != kDistPairs[i].first || dist[u2] != kDistPairs[i].second) continue;
            for (size_t j = 0; j < kDegreePairs[i].size(); ++j) {
                const auto& [d1, d2] = kDegreePairs[i][j];
                if (included(c.g.dlo[u1], c.g.dhi[u1], d1) &&
                    included(c.g.dlo[u2], c.g.dhi[u2], d2)) {
                    charge += kCharges[i][j];
                    break;
                }
            }
            break;
        }
    }
    return charge;
}

int upperBoundOfCharge(const Cartwheel& c, const std::vector<int>& applied, const Context& ctx) {
    const int d = static_cast<int>(c.spokes.size());
    int total = 10 * (6 - d);
    for (size_t j = 0; j < applied.size(); ++j) total += ctx.combined[applied[j]].r;
    for (size_t j = applied.size(); j < c.spokes.size(); ++j)
        total += amountOfPossibleChargeSend(c.g, c.spokes[j], ctx);
    for (const int e : c.spokes) total -= amountOfChargeSend(c.g, c.g.rev[e], ctx.rules);
    return total;
}

bool pruneByNonAssociatedRule(const Cartwheel& c, const std::vector<int>& applied,
                              const Context& ctx) {
    for (size_t j = 0; j < applied.size(); ++j) {
        const auto& flag = ctx.combined[applied[j]].flag;
        for (size_t k = 0; k < ctx.rules.size(); ++k) {
            if (flag[k]) continue;
            if (alwaysApply(c.g, c.spokes[j], ctx.rules[k].g, ctx.rules[k].dart)) return true;
        }
    }
    return false;
}

bool prune(const Cartwheel& c, const std::vector<int>& applied, const Context& ctx) {
    // Algorithm B.3.10 is the disjunction of three independent predicates, so evaluating them
    // in a cheaper order cannot change the answer.  The charge test is additionally computed
    // incrementally: every term still to be added is at most ctx.maxCharge and every term still
    // to be subtracted is at least 0, so a running bound that already falls to or below
    // lowerBoundOfDigonCharge proves the exact bound does too.  This is what makes the ~18M
    // wheels of degrees 7..11 tractable; it never accepts a cartwheel the literal test rejects.
    ++tlCounters.prunes;
    touchFlusher();
    if (pruneByNonAssociatedRule(c, applied, ctx)) { ++tlCounters.prunedByRule; return true; }

    const int d = static_cast<int>(c.spokes.size());
    const int digonCharge = lowerBoundOfDigonCharge(c);
    int total = 10 * (6 - d);
    for (const int j : applied) total += ctx.combined[j].r;
    int remaining = d - static_cast<int>(applied.size());

    for (const int e : c.spokes) {
        total -= amountOfChargeSend(c.g, c.g.rev[e], ctx.rules);
        if (total + ctx.maxCharge * remaining <= digonCharge) { ++tlCounters.prunedByChargeEarly; return true; }
    }
    ++tlCounters.blocked;
    if (blockedByReducibleConfiguration(c.g, c.center, ctx.kbar)) { ++tlCounters.prunedByBlocked; return true; }
    for (size_t j = applied.size(); j < c.spokes.size(); ++j) {
        --remaining;  // spokes still unaccounted for after this one
        const int threshold = digonCharge - total - ctx.maxCharge * remaining;
        const int rj = amountOfPossibleChargeSendAbove(c.g, c.spokes[j], ctx, threshold);
        if (rj <= threshold) { ++tlCounters.prunedByCharge; return true; }
        total += rj;
    }
    return false;
}

}  // namespace apex
