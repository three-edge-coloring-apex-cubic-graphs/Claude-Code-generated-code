// Parsers for the file formats described in FORMAT.md.
#pragma once

#include <string>
#include <utility>
#include <vector>

#include "core/graph.hpp"

namespace apex {

// A configuration file: only the internal vertices carry a line, so the rotations of the ring
// vertices are reconstructed from the (near-triangulated) free completion.
struct Configuration {
    std::string name;
    int N = 0, R = 0;
    std::vector<std::vector<int>> rotations;  // 0-indexed, length N, -1 marks the boundary
    std::vector<int> dlo, dhi;                // delta_K, defined on the internal vertices
    std::vector<std::pair<int, int>> digons;  // 0-indexed endpoints
};

struct Rule {
    std::string name;
    Graph g;
    int s = 0, t = 0;  // 0-indexed sender / receiver of the charge
    int r = 0;         // amount of charge
    int dart = NIL;    // the distinguished dart: head == t, tail == s
};

struct AuxiliaryRule {
    Rule base;
    std::vector<Rule> cover;  // the homomorphic cover {R_1, ..., R_k}
};

Configuration parseConfiguration(const std::string& path);
Rule parseRule(const std::string& path);
AuxiliaryRule parseAuxiliaryRule(const std::string& path);

// The free completion of a configuration, as a closed pseudo-embedding (ring included).
// Used to validate the reconstructed ring rotations.
Graph freeCompletion(const Configuration& c);

std::vector<std::string> listFiles(const std::string& dir, const std::string& suffix);

}  // namespace apex
