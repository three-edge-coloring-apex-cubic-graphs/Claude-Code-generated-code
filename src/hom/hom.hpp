// Homomorphisms and free homomorphisms.
//
//   Algorithm A.2.1  homomorphism
//   Algorithm A.3.1  freeHomomorphismTriangulation
//   Algorithm A.4.1  dartIdentification
//   Algorithm A.4.3  freeHomomorphismConfiguration        (entry point: freeHomomorphism)
//   Algorithm A.4.4  resolveDegreeIssues
//   Algorithm A.4.5-A.4.9 the degree-issue subroutines
//   Algorithm B.2.1-B.2.6 the three boundary-completion policies
//   Algorithm B.2.11-B.2.13 enforcing single digon incidence
#pragma once

#include <optional>
#include <utility>
#include <vector>

#include "core/graph.hpp"

namespace apex {

struct HomImage {
    Graph g;
    Map phi;
};

// Which boundary-completion step the free-homomorphism algorithm uses.  Section 6.3 of the
// apex paper: the substitution happens at exactly one call site, inside fixSingleDegreeIssue.
enum class Boundary {
    PseudoEmbedding,               // Algorithm B.2.1
    PseudoTriangulation,           // Algorithm B.2.3 == Algorithm A.4.8
    PseudoTriangulationWithDigons  // Algorithm B.2.2
};

// Degree-range predicates used by Algorithm A.2.1.
bool gIntersection(int lo, int hi, int lo2, int hi2);
bool gInclude(int lo, int hi, int lo2, int hi2);

// Algorithm A.2.1.  Returns the homomorphism Z -> Zs with e |-> es, or nullopt.
std::optional<Map> homomorphism(const Graph& z, int e, const Graph& zs, int es,
                                bool (*g)(int, int, int, int));

// Algorithm A.2.1 when only existence is needed.  Same result as homomorphism(...).has_value(),
// but it reuses thread-local scratch space instead of allocating; this is the innermost loop of
// the cartwheel enumeration.
bool homomorphismExists(const Graph& z, int e, const Graph& zs, int es,
                        bool (*g)(int, int, int, int));

// Algorithm A.4.1.
std::optional<HomImage> dartIdentification(const Graph& z,
                                           const std::vector<std::pair<int, int>>& requests);

// Algorithm A.4.3 (called "FreeHomomorphism" throughout the apex paper).
std::vector<HomImage> freeHomomorphism(const Graph& z,
                                       const std::vector<std::pair<int, int>>& requests,
                                       Boundary policy);

// Algorithm B.2.11.
std::vector<HomImage> enforceSingleDigonIncidence(const Graph& z, Boundary policy);

// Algorithm B.2.13.
std::vector<HomImage> freeHomomorphismAndEnforceSingleDigonIncidence(
    const Graph& z, const std::vector<std::pair<int, int>>& requests, Boundary policy);

// True iff Algorithm B.2.13 would return a non-empty set.  Same answer, but it stops at the
// first image and tracks no homomorphisms, which is all Algorithm B.3.4 (neverApply) needs.
bool freeHomomorphismAndEnforceSingleDigonIncidenceAny(
    const Graph& z, const std::vector<std::pair<int, int>>& requests, Boundary policy);

// Algorithm B.2.12.  Returns a pair of darts lying in two digons at a common vertex.
std::optional<std::pair<int, int>> twoDigonsIncidentWithSameVertex(const Graph& z);

// Algorithm B.3.6.  The digons of Z, as endpoint pairs.
std::vector<std::pair<int, int>> enumDigons(const Graph& z);

// Algorithm B.2.4, exposed because Algorithm B.4.11 calls it directly.
void addBoundaryDartsDirectly(Graph& z, int efirst, int elast);
// Algorithm B.2.6, likewise.
Map linkIncidenceListEnds(Graph& z, int eufirst, int ewlast);

}  // namespace apex
