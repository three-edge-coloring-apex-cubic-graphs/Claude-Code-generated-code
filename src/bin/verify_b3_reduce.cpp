// Lemma B.3, full claim: every multi-boundary island in I = allHomImages(K1, ...) U ... is
// semi-D-reducible or semi-C-reducible (Definitions 3.6/3.9), computed via the maximal
// semi-consistent subset of D_R \ C_I (Definition 3.5), the "standard iterative procedure" that
// CLAUDE.md notes is not published as pseudocode anywhere -- see HANDOFF.md and
// src/reduce/semi_reducible.cpp for the derivation.
//
// Islands are checked as they are produced by the same B.4.1 sweep verify_b3.cpp uses, so this
// binary reproduces both the island-construction counts (254/88393/20836/18/0 by |F_R(I)|) and
// the reducibility verdicts in a single pass, rather than needing to regenerate I twice.
//
// checkReducibility has no dependency on configuration order or Ksmaller (unlike allHomImages,
// which must see configurations strictly in filename order -- see the ground rules in
// HANDOFF.md), so it is embarrassingly parallel across islands -- but per-island cost varies
// enormously with ring size (confirmed: single-threaded, most islands take well under a second,
// but a handful take 40+ minutes each), and the expensive ones are NOT evenly spread across
// configurations. Parallelising only within each configuration's own island batch (as an earlier
// version of this file did) leaves most of 256 cores idle whenever a batch's few expensive
// islands haven't finished but its cheap ones have, even though later configurations' (already
// constructible) islands could otherwise keep every core busy. So construction (sequential, one
// configuration at a time, respecting Ksmaller ordering) and reducibility-checking (fully
// parallel) are decoupled here via a producer/consumer queue: the main thread constructs and
// enqueues islands as fast as it can, while a fixed pool of worker threads drains the queue
// continuously, drawing work from however far construction has gotten, not just the current
// configuration.
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <deque>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>

#include "hom/blocked.hpp"
#include "island/hom_images.hpp"
#include "island/outer_extension.hpp"
#include "io/parse.hpp"
#include "reduce/kempe.hpp"
#include "reduce/semi_reducible.hpp"

using namespace apex;

namespace {

const long long kExpected[5] = {254, 88393, 20836, 18, 0};

int pickExactDart(const Graph& g) {
    for (int e = 0; e < g.nd(); ++e) {
        const int h = g.head[e], t = g.tail(e);
        if (g.dlo[h] == g.dhi[h] && g.dlo[t] == g.dhi[t]) return e;
    }
    throw std::runtime_error("no dart with exact-degree endpoints");
}

struct QueueItem {
    MultiBoundaryIsland island;
    std::string confName;
    long long globalIndex = 0;
};

// A bounded blocking queue: producer blocks when full (caps memory; construction easily outruns
// checking otherwise) and consumers block when empty, until `close()` (no more items coming) and
// the queue drains.
class IslandQueue {
public:
    explicit IslandQueue(size_t cap) : cap_(cap) {}

    void push(QueueItem item) {
        std::unique_lock<std::mutex> lock(mu_);
        notFull_.wait(lock, [&] { return q_.size() < cap_; });
        q_.push_back(std::move(item));
        notEmpty_.notify_one();
    }

    bool pop(QueueItem& out) {
        std::unique_lock<std::mutex> lock(mu_);
        notEmpty_.wait(lock, [&] { return !q_.empty() || closed_; });
        if (q_.empty()) return false;
        out = std::move(q_.front());
        q_.pop_front();
        notFull_.notify_one();
        return true;
    }

    void close() {
        std::lock_guard<std::mutex> lock(mu_);
        closed_ = true;
        notEmpty_.notify_all();
    }

private:
    std::mutex mu_;
    std::condition_variable notEmpty_, notFull_;
    std::deque<QueueItem> q_;
    size_t cap_;
    bool closed_ = false;
};

}  // namespace

int main(int argc, char** argv) {
    int limit = 0;  // 0 = all
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--limit" && i + 1 < argc) limit = std::stoi(argv[++i]);
    }

    const unsigned nThreads = std::max(1u, std::thread::hardware_concurrency());
    const auto files = listFiles("configurations/K", ".conf");
    std::printf("|K| = %zu configurations, %u worker threads\n", files.size(), nThreads);

    // getValidParens/getPlanarHalfKempes memoize into non-thread-safe static maps; warm them
    // single-threaded before any concurrent checkReducibility calls (see kempe.hpp). 20 is a
    // safe margin over the largest ring size actually observed in I (15, confirmed by survey).
    warmKempeCaches(20);

    std::atomic<long long> countSemiD{0}, countSemiC{0}, countFailed{0}, countTooLarge{0};
    std::atomic<long long> countChecked{0};
    std::mutex printMu;

    IslandQueue queue(4 * nThreads);
    std::vector<std::thread> workers;
    for (unsigned t = 0; t < nThreads; ++t) {
        workers.emplace_back([&] {
            QueueItem item;
            while (queue.pop(item)) {
                try {
                    const auto r = checkReducibility(item.island);
                    if (r.semiDReducible) {
                        ++countSemiD;
                    } else if (r.semiCReducible) {
                        ++countSemiC;
                    } else {
                        ++countFailed;
                        std::lock_guard<std::mutex> lock(printMu);
                        std::printf(
                            "  ** NOT semi-D or semi-C reducible: %s island #%lld n=%d ring=%d "
                            "fixpoint=%llu **\n",
                            item.confName.c_str(), item.globalIndex, item.island.n,
                            r.ringEdgeTotal, static_cast<unsigned long long>(r.fixpointSize));
                        std::fflush(stdout);
                    }
                } catch (const std::exception& e) {
                    ++countTooLarge;
                    std::lock_guard<std::mutex> lock(printMu);
                    std::printf("  ! %s island #%lld n=%d rings=%zu: %s\n",
                               item.confName.c_str(), item.globalIndex, item.island.n,
                               item.island.ringSizes.size(), e.what());
                    std::fflush(stdout);
                }
                const long long done = ++countChecked;
                if (done % 2000 == 0) {
                    std::lock_guard<std::mutex> lock(printMu);
                    std::printf(
                        "    ... checked %lld islands so far (semiD=%lld semiC=%lld failed=%lld "
                        "tooLarge=%lld)\n",
                        done, countSemiD.load(), countSemiC.load(), countFailed.load(),
                        countTooLarge.load());
                    std::fflush(stdout);
                }
            }
        });
    }

    std::vector<ConfEntry> smallerEntries;
    ConfSet ksmaller = makeConfSet({});
    std::map<int, long long> byRing;
    long long total = 0;

    const auto t0 = std::chrono::steady_clock::now();
    const size_t n = limit > 0 ? std::min<size_t>(limit, files.size()) : files.size();
    for (size_t i = 0; i < n; ++i) {
        const auto tA = std::chrono::steady_clock::now();
        Configuration c = parseConfiguration(files[i]);
        Graph khat = outerExtension(c);
        auto islands = allHomImages(khat, ksmaller);

        for (auto& isl : islands) {
            byRing[static_cast<int>(isl.ringSizes.size())]++;
            ++total;
            queue.push(QueueItem{std::move(isl), c.name, total});
        }

        const int dart = pickExactDart(khat);
        smallerEntries.push_back({khat, dart});
        const Graph mir = mirror(khat);
        smallerEntries.push_back({mir, dart});
        ksmaller = makeConfSet(smallerEntries);

        const auto tB = std::chrono::steady_clock::now();
        {
            std::lock_guard<std::mutex> lock(printMu);
            std::printf(
                "[%4zu/%4zu] %-12s islands=%6zu constructed=%8lld checked=%8lld  "
                "(construction: %.2fs)\n",
                i + 1, n, c.name.c_str(), islands.size(), total, countChecked.load(),
                std::chrono::duration<double>(tB - tA).count());
            std::fflush(stdout);
        }
    }
    queue.close();
    for (auto& th : workers) th.join();

    std::printf("\n|F_R(I)|  count  expected\n");
    bool ok = true;
    for (int k = 0; k <= 3; ++k) {
        const long long got = byRing.count(k) ? byRing[k] : 0;
        std::printf("%8d  %5lld  %8lld  %s\n", k, got, kExpected[k],
                    got == kExpected[k] ? "OK" : "MISMATCH");
        if (limit == 0 && got != kExpected[k]) ok = false;
    }
    long long ge4 = 0;
    for (auto& [k, v] : byRing)
        if (k >= 4) ge4 += v;
    std::printf("     >=4  %5lld  %8lld  %s\n", ge4, kExpected[4],
                ge4 == kExpected[4] ? "OK" : "MISMATCH");
    if (limit == 0 && ge4 != kExpected[4]) ok = false;

    std::printf(
        "\ntotal islands = %lld\nsemi-D-reducible = %lld\nsemi-C-reducible = %lld\n"
        "NOT reducible = %lld\ntoo large to check = %lld\ntotal elapsed: %.1f s\n%s\n",
        total, countSemiD.load(), countSemiC.load(), countFailed.load(), countTooLarge.load(),
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(),
        limit == 0 ? (ok ? "ISLAND COUNTS MATCH" : "ISLAND COUNTS MISMATCH")
                   : "(partial run, --limit set)");
    return (limit == 0 && !ok) || countFailed.load() > 0 ? 1 : 0;
}
