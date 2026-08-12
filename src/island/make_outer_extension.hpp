// Algorithms B.4.8-B.4.10.
#pragma once

#include "hom/hom.hpp"

namespace apex {

// Algorithm B.4.8.  Repeatedly resolves any dart e0 whose next three darts around its head
// (e1, e2, e3) neither close a triangle/digon (e0 == e2 or e0 == e3) nor hit the boundary, until
// every pseudo-embedding in the returned set is a genuine outer extension (Claim 7.7).
std::vector<HomImage> makeOuterExtension(const Graph& z);

}  // namespace apex
