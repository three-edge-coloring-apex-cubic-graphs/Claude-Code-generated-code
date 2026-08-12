// The outer extension of a configuration (Section 4, paragraph defining K-hat; used as the input
// to Algorithm B.4.1 allHomImages).
//
// Construction: take the free completion of the configuration, then delete every ring vertex.
// Each internal vertex that was adjacent to a (now-deleted) ring vertex loses that edge; restore
// its degree to its specified delta_K by adding one pendant edge per lost connection, to a fresh
// vertex of its own. Every other internal-internal edge (and every digon) is unchanged. Internal
// vertices keep their exact delta_K (already given by the file); each new pendant vertex is the
// genuine open/boundary vertex, carrying the range [5, infinity) (Lemma 2.2's global minimum
// degree in G*), that Algorithm B.4.1's recursive free-homomorphism machinery grows.
#pragma once

#include "core/graph.hpp"
#include "io/parse.hpp"

namespace apex {

Graph outerExtension(const Configuration& c);

}  // namespace apex
