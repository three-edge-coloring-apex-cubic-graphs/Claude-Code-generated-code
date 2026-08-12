// Algorithms B.3.7-B.3.9 and B.3.11-B.3.16: enumerating and eliminating bad cartwheels.
#pragma once

#include <functional>
#include <vector>

#include "cartwheel/charge.hpp"

namespace apex {

// A cartwheel together with the combined rules already fixed on its spokes.
struct FixedCartwheel {
    Cartwheel c;
    std::vector<int> applied;  // indices into Context::combined
};

// Algorithm B.3.14.  `progress`, when set, is called occasionally with (done, total).
std::vector<Cartwheel> enumPossibleBadWheels(
    int d, const Context& ctx,
    const std::function<void(long long, long long)>& progress = {});

// Algorithm B.3.8.
std::vector<Cartwheel> updateDegreeByRule(const Cartwheel& c, int dart, const Graph& rg,
                                          int rdart);

// Algorithm B.3.7.
std::vector<FixedCartwheel> fixInRules(const Cartwheel& c, const Context& ctx);

// Algorithm B.3.11.
std::vector<FixedCartwheel> fixOutRules(const std::vector<FixedCartwheel>& cd, const Context& ctx);

// Algorithm B.3.15.  Returns the surviving set C, which the assertion requires to be empty.
std::vector<FixedCartwheel> verifyNoBadCartwheels(const Cartwheel& c, const Context& ctx);

}  // namespace apex
