// Algorithms B.2.7-B.2.9 (the apex variants of A.6.1-A.6.3) plus A.6.4 and A.6.5:
// construction of the configuration set K-bar.
#pragma once

#include <vector>

#include "hom/blocked.hpp"
#include "io/parse.hpp"

namespace apex {

// Algorithm B.2.7.
std::vector<ConfEntry> extendFromCutVertices(const Configuration& c);

// K-bar: every extension of every configuration, together with all mirror images.
ConfSet buildKBar(const std::vector<Configuration>& confs);

}  // namespace apex
