// Algorithms B.1.1-B.1.2: generating semi-matchings (Definition 3.3) for a single boundary.
#pragma once

#include <utility>
#include <vector>

namespace apex {

// A partial matching on a 0-indexed cyclically ordered point set. Each pair (a, b) has a < b.
using Matching = std::vector<std::pair<int, int>>;

// Algorithm B.1.1. The Catalan family of noncrossing perfect matchings on 2q points, 0-indexed
// as {0, ..., 2q-1} (the paper's [2q] = {1, ..., 2q}). Memoized; the returned reference is stable
// for the lifetime of the program.
const std::vector<Matching>& getValidParens(int q);

// Algorithm B.1.2. The set of nonredundant planar partial matchings on the cyclically ordered,
// 0-indexed set {0, ..., n-1} (the paper's [n] = {1, ..., n}). Memoized. n must be >= 1; callers
// with an empty position set (n == 0) should treat it as the single trivial candidate {} directly
// rather than calling this function.
const std::vector<Matching>& getPlanarHalfKempes(int n);

// Populates the getValidParens/getPlanarHalfKempes memo tables for every n in [1, maxN], single-
// threaded. Both functions memoize into function-local `static` maps that are NOT safe under
// concurrent first-population (a data race: e.g. two threads both missing the same cache entry
// can corrupt the map while inserting into it). Callers that query these functions from multiple
// threads (checkReducibility does, via componentCandidates) MUST call this once, single-threaded,
// for the largest n they will ever query, before starting concurrent work; every call thereafter
// is a read-only std::map lookup, which is safe.
void warmKempeCaches(int maxN);

}  // namespace apex
