// Enumerates the boundary colourings C_I (Definition 3.6) and C_{I-F} (Definition 3.9) of a
// multi-boundary island by directly searching for three-edge-colourings (resp. three-edge-
// colourings modulo a deletable edge set F) of the underlying graph I.
#pragma once

#include <vector>

#include "island/island.hpp"

namespace apex {

// Definition 3.6's C_I: boundary colourings phi: R -> {0,1,2} that extend to a genuine
// three-edge-colouring eta of I. Each returned vector has length equal to the number of ring
// edges (sum of ringSizes), in FORMAT.md's edge-index order (ring components concatenated in
// order, each in clockwise order). The auxiliary pendant edges are not part of I; padding a
// degree-two vertex's two real edges with a pendant third slot and requiring all three pairwise
// distinct is equivalent to requiring the two real edges distinct (the pendant's colour is then
// forced, not a free choice), so computeCI can use the same "all three slots pairwise distinct"
// rule uniformly at every vertex.
std::vector<std::vector<int>> computeCI(const MultiBoundaryIsland& island);

// Definition 3.9's C_{I-F}: boundary colourings extendable to a three-edge-colouring of I modulo
// the deletable edge set F, whose indices (FORMAT.md numbering, always among the "other" -- non-
// ring, non-pendant -- edges) are the true entries of deletedOther (size = total edge count;
// only entries at "other"-edge positions are consulted). Per Definition 3.8, two edges e, f of
// E(I) sharing a vertex v get equal colours iff v is incident to exactly one edge of F; pendant
// edges are not part of E(I) and are excluded from the search entirely.
std::vector<std::vector<int>> computeCIModuloF(const MultiBoundaryIsland& island,
                                               const std::vector<char>& deletedOther);

}  // namespace apex
