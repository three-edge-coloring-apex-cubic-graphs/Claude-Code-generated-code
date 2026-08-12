// Lemma B.1: build R*^{-K} from the rule set R and the configuration set K-bar, and check
// that every combined rule R* satisfies r(R*) <= 5.
//
// Metrics.md expects |R*^{-K}| = 747 and a maximum charge of 5.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>

#include "hom/hom.hpp"
#include "kbar/kbar.hpp"
#include "rules/combine.hpp"

using namespace apex;

int main(int argc, char** argv) {
    bool verbose = false, noBlock = false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "-v") verbose = true;
        // Diagnostic only: skip the "blocked by a configuration in K-bar" filter, which must
        // change the answer -- otherwise the filter is not being exercised at all.
        if (a == "--no-block") noBlock = true;
    }
    const auto t0 = std::chrono::steady_clock::now();

    std::vector<Configuration> confs;
    for (const auto& p : listFiles("configurations/K", ".conf")) confs.push_back(parseConfiguration(p));
    std::printf("configurations: %zu\n", confs.size());

    ConfSet kbar = buildKBar(confs);
    std::printf("|K-bar| = %zu  (CONF_DEG_MAX = %d)\n", kbar.confs.size(), kbar.confDegMax);
    if (noBlock) kbar = ConfSet{};

    std::vector<Rule> rules;
    for (const auto& p : listFiles("discharging-rules/R", ".rule")) rules.push_back(parseRule(p));
    std::printf("rules: %zu\n", rules.size());

    const auto star = combineRules(rules, kbar, verbose);

    int maxCharge = 0, withDigon = 0;
    std::map<int, int> byCharge;
    for (const auto& c : star) {
        maxCharge = std::max(maxCharge, c.r);
        byCharge[c.r]++;
        if (!enumDigons(c.g).empty()) ++withDigon;
    }

    const auto t1 = std::chrono::steady_clock::now();
    std::printf("\n=== Lemma B.1 ===\n");
    std::printf("|R*^-K|                                  = %zu   (expected 747)\n", star.size());
    std::printf("max charge of a combined rule in R*^-K   = %d     (expected 5)\n", maxCharge);
    std::printf("combined rules containing a digon         = %d\n", withDigon);
    std::printf("charge distribution:");
    for (auto [r, n] : byCharge) std::printf(" %d:%d", r, n);
    std::printf("\nelapsed: %.1f s\n",
                std::chrono::duration<double>(t1 - t0).count());

    const bool ok = star.size() == 747 && maxCharge == 5;
    std::printf("%s\n", ok ? "MATCH" : "MISMATCH");
    return ok ? 0 : 1;
}
