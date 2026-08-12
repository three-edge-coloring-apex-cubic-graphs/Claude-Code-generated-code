#include "rules/combine.hpp"

#include <atomic>
#include <cstdio>
#include <mutex>
#include <thread>
#include <unordered_set>

#include "hom/hom.hpp"

namespace apex {
namespace {

// The trivial combined rule of Algorithm A.8.2 line 2: a single edge s0 t0, both endpoints of
// degree range [1, infinity), with the distinguished dart pointing at t0.
CombinedRule trivialRule(size_t nRules) {
    CombinedRule c;
    c.g.addVertex(1, INF);  // s0
    c.g.addVertex(1, INF);  // t0
    const int d0 = c.g.addDart();
    const int d1 = c.g.addDart();
    c.g.head[d0] = 1;
    c.g.rev[d0] = d1;
    c.g.head[d1] = 0;
    c.g.rev[d1] = d0;
    c.dart = d0;
    c.r = 0;
    c.flag.assign(nRules, 0);
    return c;
}

std::string key(const CombinedRule& c) {
    std::string s = canonicalKey(c.g, c.dart);
    s += '|';
    s += std::to_string(c.r);
    s += '|';
    for (char f : c.flag) s += static_cast<char>('0' + f);
    return s;
}

// Algorithm A.8.1.
std::vector<CombinedRule> addRuleToCombination(const CombinedRule& star, const Rule& rule,
                                               int ruleIndex, const ConfSet& k) {
    const Graph joined = disjointUnion(star.g, rule.g);
    const int e0 = star.dart;
    const int e1 = rule.dart + star.g.nd();

    std::vector<CombinedRule> out;
    for (auto& img : freeHomomorphismAndEnforceSingleDigonIncidence(
             joined, {{e0, e1}}, Boundary::PseudoTriangulationWithDigons)) {
        CombinedRule c;
        c.g = std::move(img.g);
        c.dart = img.phi.d[e0];
        c.r = star.r + rule.r;
        c.flag = star.flag;
        c.flag[ruleIndex] = 1;
        if (!k.confs.empty() &&
            blockedByReducibleConfiguration(c.g, c.g.head[c.dart], k))
            continue;
        out.push_back(std::move(c));
    }
    return out;
}

}  // namespace

std::vector<CombinedRule> combineRules(const std::vector<Rule>& rules, const ConfSet& k,
                                       bool verbose) {
    std::vector<CombinedRule> star{trivialRule(rules.size())};
    std::unordered_set<std::string> seen{key(star[0])};

    const unsigned nThreads = std::max(1u, std::thread::hardware_concurrency());

    for (size_t i = 0; i < rules.size(); ++i) {
        // Snapshot: workers read `prev` while results are appended to `star`.
        const std::vector<CombinedRule> prev = star;
        std::mutex mu;
        std::atomic<size_t> next{0};
        std::vector<std::thread> pool;
        for (unsigned t = 0; t < nThreads; ++t) {
            pool.emplace_back([&] {
                for (;;) {
                    const size_t j = next.fetch_add(1);
                    if (j >= prev.size()) return;
                    auto produced = addRuleToCombination(prev[j], rules[i], i, k);
                    std::lock_guard<std::mutex> lock(mu);
                    for (auto& c : produced) {
                        if (!seen.insert(key(c)).second) continue;
                        star.push_back(std::move(c));
                    }
                }
            });
        }
        for (auto& th : pool) th.join();
        if (verbose)
            std::printf("  after %s: |R*| = %zu\n", rules[i].name.c_str(), star.size());
    }
    return star;
}

}  // namespace apex
