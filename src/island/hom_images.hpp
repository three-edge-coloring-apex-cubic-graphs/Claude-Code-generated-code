// Algorithm B.4.1.
#pragma once

#include "hom/blocked.hpp"
#include "island/island.hpp"

namespace apex {

std::vector<MultiBoundaryIsland> allHomImages(const Graph& khat, const ConfSet& ksmaller);

}  // namespace apex
