// Parser smoke test: every configuration, rule and auxiliary rule must parse, and every free
// completion must be a planar near-triangulation whose only non-triangular face is the ring.
#include <cstdio>
#include <map>
#include <stdexcept>

#include "core/walks.hpp"
#include "io/parse.hpp"

using namespace apex;

int main() {
    int failures = 0;
    std::map<int, int> ringSizes;
    int maxDelta = 0;

    for (const auto& p : listFiles("configurations/K", ".conf")) {
        try {
            Configuration c = parseConfiguration(p);
            ringSizes[c.R]++;
            for (int v = c.R; v < c.N; ++v) maxDelta = std::max(maxDelta, c.dlo[v]);

            Graph g = freeCompletion(c);
            if (!isPlanar(g)) {
                std::printf("NOT PLANAR: %s\n", c.name.c_str());
                ++failures;
                continue;
            }
            // Faces: every face a triangle or a digon, except one outer face of size R.
            int outer = 0;
            for (const auto& w : getWalks(g)) {
                const int len = static_cast<int>(w.size());
                if (len == static_cast<int>(c.R) && c.R != 3) {
                    ++outer;
                } else if (len != 3 && len != 2) {
                    std::printf("BAD FACE %d in %s\n", len, c.name.c_str());
                    ++failures;
                }
            }
            if (c.R != 3 && outer != 1) {
                std::printf("outer-face count %d in %s (R=%d)\n", outer, c.name.c_str(), c.R);
                ++failures;
            }
        } catch (const std::exception& e) {
            std::printf("EXCEPTION %s: %s\n", p.c_str(), e.what());
            ++failures;
        }
    }

    std::printf("configurations parsed; ring sizes:");
    for (auto [r, n] : ringSizes) std::printf(" %d:%d", r, n);
    std::printf("\nmax delta_K = %d\n", maxDelta);

    int nRules = 0;
    for (const auto& p : listFiles("discharging-rules/R", ".rule")) {
        try {
            Rule r = parseRule(p);
            if (!isPlanar(r.g)) {
                std::printf("rule not planar: %s\n", r.name.c_str());
                ++failures;
            }
            ++nRules;
        } catch (const std::exception& e) {
            std::printf("EXCEPTION %s: %s\n", p.c_str(), e.what());
            ++failures;
        }
    }
    std::printf("rules parsed: %d\n", nRules);

    int nAux = 0, nCover = 0;
    for (const auto& p : listFiles("discharging-rules/R_auxiliary", ".rule_auxiliary")) {
        try {
            AuxiliaryRule a = parseAuxiliaryRule(p);
            nCover += static_cast<int>(a.cover.size());
            ++nAux;
        } catch (const std::exception& e) {
            std::printf("EXCEPTION %s: %s\n", p.c_str(), e.what());
            ++failures;
        }
    }
    std::printf("auxiliary rules parsed: %d (cover rules %d)\n", nAux, nCover);

    std::printf(failures ? "FAILURES: %d\n" : "all parser checks passed (%d failures)\n", failures);
    return failures ? 1 : 0;
}
