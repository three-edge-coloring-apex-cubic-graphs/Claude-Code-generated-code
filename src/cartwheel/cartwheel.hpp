// Algorithms B.3.1, B.3.2, B.3.3: wheels and cartwheels.
#pragma once

#include <utility>
#include <vector>

#include "core/graph.hpp"

namespace apex {

// CARTWHEEL_DEGREES from Section B.3.1.  The paper prints [(5,5),(6,6),(7,7),(8,8),(9,inf)];
// note that the fourth entry is a literal 8 and only the last is infinity.
inline const std::vector<std::pair<int, int>>& cartwheelDegrees() {
    static const std::vector<std::pair<int, int>> v{{5, 5}, {6, 6}, {7, 7}, {8, 8}, {9, INF}};
    return v;
}

// A cartwheel: a pseudo-triangulation with digons with degree ranges, a center vertex, and the
// images e_1, ..., e_d of the center darts in clockwise order.
struct Cartwheel {
    Graph g;
    int center = NIL;
    std::vector<int> spokes;
};

// Algorithm B.3.3.
Cartwheel generateCartwheel(int d, const std::vector<std::pair<int, int>>& degrees,
                            bool incidentDigon);

// Algorithm B.3.2.
std::vector<Cartwheel> enumWheels(int d);

// Algorithm B.3.1.
std::vector<Cartwheel> enumDigonIncidentWheels(int d);

// Index-addressable forms of the two enumerations above.  They produce exactly the same
// families, but one wheel at a time: for d = 11 the digon-incident family alone has 5^10
// members, far too many to hold in memory at once.
long long digonIncidentWheelCount(int d);  // 5^(d-1)
Cartwheel digonIncidentWheelAt(int d, long long idx);

// The candidate space explored by Algorithm B.3.2's i_lowest loops, of size sum_j (5-j)^(d-1).
// Returns false when `degrees` is not the lexicographically minimal rotation, i.e. exactly when
// the check on line 5 of Algorithm B.3.2 rejects it.
long long wheelCandidateCount(int d);
bool wheelCandidateAt(int d, long long idx, Cartwheel& out);

}  // namespace apex
