// Charges along an edge, bounds on the final charge, and pruning.
//
//   Algorithm A.9.1  alwaysApply
//   Algorithm B.3.4  neverApply       (the apex version: a free-homomorphism emptiness test)
//   Algorithm A.9.3  amountOfChargeSend
//   Algorithm A.9.4  amountOfPossibleChargeSend
//   Algorithm A.9.13 upperBoundOfCharge
//   Algorithm A.9.12 pruneByNonAssociatedRule
//   Algorithm B.3.5  lowerBoundOfDigonCharge
//   Algorithm B.3.10 prune
#pragma once

#include <vector>

#include "cartwheel/cartwheel.hpp"
#include "hom/blocked.hpp"
#include "io/parse.hpp"
#include "rules/combine.hpp"

namespace apex {

// Everything the cartwheel algorithms need besides the cartwheel itself.
struct Context {
    std::vector<Rule> rules;             // R
    std::vector<CombinedRule> combined;  // R*^{-K}
    std::vector<int> byChargeDesc;       // indices into `combined`, decreasing charge
    ConfSet kbar;                        // K-bar
    std::vector<AuxiliaryRule> auxiliary;  // R_auxiliary
    int maxCharge = 0;                     // max r(R*) over R*^{-K}

    void finalise();
};

bool alwaysApply(const Graph& z, int e, const Graph& rg, int rdart);
bool neverApply(const Graph& z, int e, const Graph& rg, int rdart);

int amountOfChargeSend(const Graph& z, int e, const std::vector<Rule>& rules);
int amountOfPossibleChargeSend(const Graph& z, int e, const Context& ctx);
// Exact when the result exceeds `threshold`; otherwise returns some value <= threshold.
int amountOfPossibleChargeSendAbove(const Graph& z, int e, const Context& ctx, int threshold);

int lowerBoundOfDigonCharge(const Cartwheel& c);
int upperBoundOfCharge(const Cartwheel& c, const std::vector<int>& applied, const Context& ctx);
bool pruneByNonAssociatedRule(const Cartwheel& c, const std::vector<int>& applied,
                              const Context& ctx);
bool prune(const Cartwheel& c, const std::vector<int>& applied, const Context& ctx);

// Call counters, for locating the cost of the enumeration.
struct Counters {
    long long prunes = 0, alwaysApply = 0, neverApply = 0, neverApplyDeep = 0, blocked = 0,
        digonCharge = 0;
    long long prunedByRule = 0, prunedByChargeEarly = 0, prunedByBlocked = 0, prunedByCharge = 0;
};
Counters takeCounters();

}  // namespace apex
