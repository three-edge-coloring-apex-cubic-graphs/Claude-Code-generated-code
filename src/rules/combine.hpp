// Algorithms A.8.1 and A.8.2: free combinations of discharging rules.
//
// Per Lemma B.1, the routine that computes free homomorphic images of pseudo-triangulations is
// replaced by Algorithm B.2.13 for pseudo-triangulations with digons enforcing single digon
// incidence.
#pragma once

#include <string>
#include <vector>

#include "hom/blocked.hpp"
#include "io/parse.hpp"

namespace apex {

struct CombinedRule {
    Graph g;
    int dart = NIL;  // the distinguished dart, head == the vertex receiving the charge
    int r = 0;       // amount of charge
    std::vector<char> flag;  // which rules of R were combined
};

std::vector<CombinedRule> combineRules(const std::vector<Rule>& rules, const ConfSet& k,
                                       bool verbose = false);

}  // namespace apex
