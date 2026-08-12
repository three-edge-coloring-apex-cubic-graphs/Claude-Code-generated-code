// Definitions 3.5/3.6/3.9: the maximal semi-consistent subset of D_R \ C_I, computed by the
// "standard iterative procedure" CLAUDE.md asks for (not given as pseudocode anywhere): start
// from S = D_R \ C_I and repeatedly delete any phi for which some colour pair {x,y} has no
// semi-matching M (drawn from GetPlanarHalfKempes per boundary component, concatenated per
// B.1.1/B.1.2's own remark) whose every Kempe-switch M' subseteq M stays inside the current S;
// iterate to a fixpoint. I is semi-D-reducible iff the fixpoint is empty; semi-C-reducible by a
// deletable F iff the fixpoint is disjoint from C_{I-F}.
#pragma once

#include <cstdint>
#include <vector>

#include "island/island.hpp"

namespace apex {

struct ReducibilityResult {
    bool semiDReducible = false;
    bool semiCReducible = false;   // only meaningful if !semiDReducible
    std::vector<char> witnessF;    // the deletable F that worked, if semiCReducible
    uint64_t fixpointSize = 0;     // |maximal semi-consistent subset of D_R \ C_I|, for logging
    int ringEdgeTotal = 0;         // |R|, for logging / cost accounting
};

// Throws std::runtime_error if 3^|R| is too large to enumerate densely (see semi_reducible.cpp
// for the cap); callers should catch this to flag the island for separate handling.
ReducibilityResult checkReducibility(const MultiBoundaryIsland& island);

}  // namespace apex
