#include "cartwheel/verify.hpp"

#include <algorithm>
#include <atomic>
#include <deque>
#include <thread>

#include "hom/hom.hpp"

namespace apex {
namespace {

// Algorithm B.3.9.
std::vector<Cartwheel> concreteDegreeExceptTail(const Cartwheel& c) {
    std::vector<Cartwheel> out{c};
    for (int u = 0; u < c.g.nv; ++u) {
        if (c.g.dlo[u] == c.g.dhi[u] || c.g.dhi[u] == INF) continue;
        std::vector<Cartwheel> next;
        for (int d = c.g.dlo[u]; d <= c.g.dhi[u]; ++d)
            for (const auto& z : out) {
                Cartwheel z2 = z;
                z2.g.dlo[u] = d;
                z2.g.dhi[u] = d;
                next.push_back(std::move(z2));
            }
        out = std::move(next);
    }
    return out;
}

// Algorithm B.3.12.
bool shouldRefine(const Cartwheel& c, int i, const AuxiliaryRule& aux) {
    const int e = c.g.rev[c.spokes[i]];
    if (!alwaysApply(c.g, e, aux.base.g, aux.base.dart)) return false;
    for (const auto& sub : aux.cover)
        if (alwaysApply(c.g, e, sub.g, sub.dart)) return false;
    return true;
}

}  // namespace

std::vector<Cartwheel> enumPossibleBadWheels(int d, const Context& ctx,
                                             const std::function<void(long long, long long)>& progress) {
    // Algorithm B.3.14, streamed: the two wheel families are generated one member at a time and
    // only the survivors are kept, since for d = 11 the families together have ~14M members.
    const long long nWheel = wheelCandidateCount(d);
    const long long nDigon = digonIncidentWheelCount(d);
    const long long total = nWheel + nDigon;

    const unsigned nThreads = std::max(1u, std::thread::hardware_concurrency());
    constexpr long long kChunk = 2048;
    const long long nChunks = (total + kChunk - 1) / kChunk;

    std::vector<std::vector<Cartwheel>> perChunk(static_cast<size_t>(nChunks));
    std::atomic<long long> nextChunk{0}, doneChunks{0};
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < nThreads; ++t) {
        pool.emplace_back([&] {
            for (;;) {
                const long long ch = nextChunk.fetch_add(1);
                if (ch >= nChunks) return;
                auto& out = perChunk[static_cast<size_t>(ch)];
                const long long lo = ch * kChunk, hi = std::min(total, lo + kChunk);
                Cartwheel w;
                for (long long i = lo; i < hi; ++i) {
                    if (i < nWheel) {
                        if (!wheelCandidateAt(d, i, w)) continue;
                    } else {
                        w = digonIncidentWheelAt(d, i - nWheel);
                    }
                    if (prune(w, {}, ctx)) continue;
                    out.push_back(w);
                }
                const long long k = doneChunks.fetch_add(1) + 1;
                if (progress && k % 256 == 0) progress(k * kChunk, total);
            }
        });
    }
    for (auto& th : pool) th.join();

    std::vector<Cartwheel> out;
    for (auto& v : perChunk)
        for (auto& w : v) out.push_back(std::move(w));
    return out;
}

std::vector<Cartwheel> updateDegreeByRule(const Cartwheel& c, int dart, const Graph& rg,
                                          int rdart) {
    const Graph joined = disjointUnion(c.g, rg);
    std::vector<Cartwheel> out;
    for (auto& img : freeHomomorphismAndEnforceSingleDigonIncidence(
             joined, {{dart, rdart + c.g.nd()}}, Boundary::PseudoTriangulationWithDigons)) {
        Cartwheel star;
        star.g = std::move(img.g);
        star.center = img.phi.v[c.center];
        for (const int e : c.spokes) star.spokes.push_back(img.phi.d[e]);
        for (auto& z : concreteDegreeExceptTail(star)) out.push_back(std::move(z));
    }
    return out;
}

std::vector<FixedCartwheel> fixInRules(const Cartwheel& c, const Context& ctx) {
    std::vector<FixedCartwheel> cur{{c, {}}};
    for (size_t i = 0; i < c.spokes.size(); ++i) {
        std::vector<FixedCartwheel> next;
        for (const auto& f : cur) {
            for (size_t ri = 0; ri < ctx.combined.size(); ++ri) {
                const auto& rule = ctx.combined[ri];
                for (auto& z : updateDegreeByRule(f.c, f.c.spokes[i], rule.g, rule.dart)) {
                    std::vector<int> applied = f.applied;
                    applied.push_back(static_cast<int>(ri));
                    if (prune(z, applied, ctx)) continue;
                    next.push_back({std::move(z), std::move(applied)});
                }
            }
        }
        cur = std::move(next);
    }
    return cur;
}

std::vector<FixedCartwheel> fixOutRules(const std::vector<FixedCartwheel>& cd,
                                        const Context& ctx) {
    std::deque<FixedCartwheel> q(cd.begin(), cd.end());
    std::vector<FixedCartwheel> out;
    while (!q.empty()) {
        FixedCartwheel cur = std::move(q.front());
        q.pop_front();
        bool refined = false;
        for (size_t i = 0; i < cur.c.spokes.size() && !refined; ++i) {
            for (const auto& aux : ctx.auxiliary) {
                if (!shouldRefine(cur.c, static_cast<int>(i), aux)) continue;
                refined = true;
                // Algorithm B.3.13.
                for (const auto& sub : aux.cover) {
                    for (auto& z : updateDegreeByRule(cur.c, cur.c.g.rev[cur.c.spokes[i]], sub.g,
                                                      sub.dart)) {
                        if (prune(z, cur.applied, ctx)) continue;
                        q.push_back({std::move(z), cur.applied});
                    }
                }
                break;
            }
        }
        if (!refined) out.push_back(std::move(cur));
    }
    return out;
}

std::vector<FixedCartwheel> verifyNoBadCartwheels(const Cartwheel& c, const Context& ctx) {
    return fixOutRules(fixInRules(c, ctx), ctx);
}

}  // namespace apex
