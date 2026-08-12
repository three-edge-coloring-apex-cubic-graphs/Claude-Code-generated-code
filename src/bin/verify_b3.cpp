// Lemma B.3: construct I = union over K1, ..., K|K| (in filename order) of
// allHomImages(Ki, {K1, ..., K_{i-1}} and mirrors), and report the number of multi-boundary
// islands grouped by |F_R(I)|.
//
// Metrics.md expects 254 / 88393 / 20836 / 18 / 0 for |F_R(I)| = 0 / 1 / 2 / 3 / >= 4.
#include <chrono>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <string>

#include "hom/blocked.hpp"
#include "island/hom_images.hpp"
#include "island/outer_extension.hpp"
#include "io/parse.hpp"

using namespace apex;

namespace {

const long long kExpected[5] = {254, 88393, 20836, 18, 0};

// A dart with exact (dlo == dhi) degree at both endpoints, for ConfSet's bucketing; ring
// vertices carry the range [5, INF), so this always lands on an internal-internal dart.
int pickExactDart(const Graph& g) {
    for (int e = 0; e < g.nd(); ++e) {
        const int h = g.head[e], t = g.tail(e);
        if (g.dlo[h] == g.dhi[h] && g.dlo[t] == g.dhi[t]) return e;
    }
    throw std::runtime_error("no dart with exact-degree endpoints");
}

}  // namespace

int main(int argc, char** argv) {
    int limit = 0;  // 0 = all
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--limit" && i + 1 < argc) limit = std::stoi(argv[++i]);
    }

    const auto files = listFiles("configurations/K", ".conf");
    std::printf("|K| = %zu configurations\n", files.size());

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
        for (auto& isl : islands) byRing[static_cast<int>(isl.ringSizes.size())]++;
        total += islands.size();

        const int dart = pickExactDart(khat);
        smallerEntries.push_back({khat, dart});
        const Graph mir = mirror(khat);
        smallerEntries.push_back({mir, dart});
        ksmaller = makeConfSet(smallerEntries);

        const auto tB = std::chrono::steady_clock::now();
        std::printf("[%4zu/%4zu] %-12s islands=%6zu total=%8lld  (%.2fs)\n", i + 1, n,
                    c.name.c_str(), islands.size(), total,
                    std::chrono::duration<double>(tB - tA).count());
        std::fflush(stdout);
    }

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

    std::printf("\ntotal islands = %lld\ntotal elapsed: %.1f s\n%s\n", total,
                std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(),
                limit == 0 ? (ok ? "MATCH" : "MISMATCH") : "(partial run, --limit set)");
    return (limit == 0 && !ok) ? 1 : 0;
}
