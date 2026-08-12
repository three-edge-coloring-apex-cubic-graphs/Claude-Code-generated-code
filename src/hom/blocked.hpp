// Algorithms A.6.6-A.6.8 (containConf) and A.7.1-A.7.2 (blockedByReducibleConfiguration).
#pragma once

#include <vector>

#include "core/graph.hpp"

namespace apex {

// A configuration of K-bar together with the special dart chosen by Algorithm A.6.4.
struct ConfEntry {
    Graph g;
    int dart = NIL;
};

struct ConfSet {
    std::vector<ConfEntry> confs;
    int confDegMax = 0;  // CONF_DEG_MAX: the largest fixed degree occurring in the set
    // confs indexed by the endpoint degrees (dy, dx) of their special dart, so that
    // containConf can iterate over the darts of the target instead of over all configurations.
    std::vector<std::vector<int>> byDeg;
    const std::vector<int>& bucket(int dy, int dx) const {
        return byDeg[dy * (confDegMax + 1) + dx];
    }
};

ConfSet makeConfSet(std::vector<ConfEntry> confs);

// Algorithm A.7.1.  `center` is the distinguished vertex c*; pass NIL for the uncentered
// variant used by Algorithm B.4.1 (see the preamble of Section B.4 of the apex paper).
bool blockedByReducibleConfiguration(const Graph& z, int center, const ConfSet& k);

}  // namespace apex
