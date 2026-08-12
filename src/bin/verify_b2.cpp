// Lemma B.2: run Algorithm B.3.16 and check that the assertion inside every call to
// Algorithm B.3.15 holds.
//
// Metrics.md expects enumPossibleBadWheels to return 4438 / 4939 / 2409 / 567 / 38 cartwheels
// for center degrees 7 / 8 / 9 / 10 / 11.
#include <atomic>
#include <chrono>
#include <cstdio>
#include <mutex>
#include <string>
#include <thread>

#include "cartwheel/verify.hpp"
#include "kbar/kbar.hpp"

using namespace apex;

namespace {

const int kExpected[5] = {4438, 4939, 2409, 567, 38};

Context loadContext() {
    Context ctx;
    std::vector<Configuration> confs;
    for (const auto& p : listFiles("configurations/K", ".conf"))
        confs.push_back(parseConfiguration(p));
    ctx.kbar = buildKBar(confs);
    for (const auto& p : listFiles("discharging-rules/R", ".rule"))
        ctx.rules.push_back(parseRule(p));
    for (const auto& p : listFiles("discharging-rules/R_auxiliary", ".rule_auxiliary"))
        ctx.auxiliary.push_back(parseAuxiliaryRule(p));
    ctx.combined = combineRules(ctx.rules, ctx.kbar);
    ctx.finalise();
    std::printf("|K-bar| = %zu, |R| = %zu, |R*^-K| = %zu, |R_auxiliary| = %zu\n",
                ctx.kbar.confs.size(), ctx.rules.size(), ctx.combined.size(),
                ctx.auxiliary.size());
    return ctx;
}

// Parallel map over a work list.
template <typename F>
void parallelFor(size_t n, F&& f) {
    const unsigned nThreads = std::max(1u, std::thread::hardware_concurrency());
    std::atomic<size_t> next{0};
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < nThreads; ++t)
        pool.emplace_back([&] {
            for (;;) {
                const size_t i = next.fetch_add(1);
                if (i >= n) return;
                f(i);
            }
        });
    for (auto& th : pool) th.join();
}

}  // namespace

int main(int argc, char** argv) {
    bool countsOnly = false;
    int only = 0;  // restrict to a single center degree, for quick iteration
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--counts-only") countsOnly = true;
        if (a == "--degree" && i + 1 < argc) only = std::stoi(argv[++i]);
    }

    const auto t0 = std::chrono::steady_clock::now();
    const Context ctx = loadContext();
    std::fflush(stdout);

    bool ok = true;
    for (int d = 7; d <= 11; ++d) {
        if (only && d != only) continue;
        const auto tA = std::chrono::steady_clock::now();
        const auto wheels = enumPossibleBadWheels(d, ctx, [&](long long k, long long n) {
            std::printf("      d=%2d scanned %lld/%lld wheels\n", d, k, n);
            std::fflush(stdout);
        });
        const auto tB = std::chrono::steady_clock::now();
        const int expected = kExpected[d - 7];
        std::printf("d=%2d: enumPossibleBadWheels -> %6zu  (expected %5d) %s   [%.1f s]\n", d,
                    wheels.size(), expected, wheels.size() == static_cast<size_t>(expected) ? "OK" : "MISMATCH",
                    std::chrono::duration<double>(tB - tA).count());
        std::fflush(stdout);
        {
            const auto ct = takeCounters();
            std::printf("      calls: prune=%lld alwaysApply=%lld neverApply=%lld (deep %lld) blocked=%lld\n"
                        "      pruned by: rule=%lld chargeEarly=%lld blocked=%lld charge=%lld\n",
                        ct.prunes, ct.alwaysApply, ct.neverApply, ct.neverApplyDeep, ct.blocked, ct.prunedByRule,
                        ct.prunedByChargeEarly, ct.prunedByBlocked, ct.prunedByCharge);
            std::fflush(stdout);
        }
        if (wheels.size() != static_cast<size_t>(expected)) ok = false;
        if (countsOnly) continue;

        std::atomic<size_t> bad{0}, done{0};
        std::mutex mu;
        parallelFor(wheels.size(), [&](size_t i) {
            const auto c = verifyNoBadCartwheels(wheels[i], ctx);
            if (!c.empty()) bad.fetch_add(c.size());
            const size_t k = done.fetch_add(1) + 1;
            if (k % 200 == 0) {
                std::lock_guard<std::mutex> lock(mu);
                std::printf("      verified %zu/%zu\n", k, wheels.size());
                std::fflush(stdout);
            }
        });
        const auto tC = std::chrono::steady_clock::now();
        std::printf("d=%2d: verifyNoBadCartwheels -> C is %s   [%.1f s]\n", d,
                    bad.load() == 0 ? "empty (assertion holds)" : "NON-EMPTY (ASSERTION FAILED)",
                    std::chrono::duration<double>(tC - tB).count());
        std::fflush(stdout);
        if (bad.load() != 0) ok = false;
    }

    std::printf("\ntotal elapsed: %.1f s\n%s\n",
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(),
                ok ? "MATCH" : "MISMATCH");
    return ok ? 0 : 1;
}
