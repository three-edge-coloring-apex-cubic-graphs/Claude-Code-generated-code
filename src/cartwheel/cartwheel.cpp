#include "cartwheel/cartwheel.hpp"

#include "core/rotations.hpp"

namespace apex {

Cartwheel generateCartwheel(int d, const std::vector<std::pair<int, int>>& degrees,
                            bool incidentDigon) {
    const int n = incidentDigon ? d - 1 : d;
    std::vector<std::vector<int>> rotations(n + 1);
    for (int i = 1; i <= n; ++i) rotations[0].push_back(i);
    for (int i = 1; i <= n; ++i) {
        const int i1 = i < n ? i + 1 : 1;
        const int i2 = i > 1 ? i - 1 : n;
        rotations[i] = {i1, 0, i2, NIL};
    }
    std::vector<std::pair<int, int>> digons;
    if (incidentDigon) digons.emplace_back(0, 1);

    Cartwheel c;
    c.g = fromVRotations(n + 1, rotations, digons);
    c.g.dlo[0] = c.g.dhi[0] = d;
    for (int i = 1; i <= n; ++i) {
        c.g.dlo[i] = degrees[i - 1].first;
        c.g.dhi[i] = degrees[i - 1].second;
    }
    c.center = 0;
    c.spokes = incidence(c.g)[0];
    return c;
}

std::vector<Cartwheel> enumDigonIncidentWheels(int d) {
    const auto& cd = cartwheelDegrees();
    std::vector<Cartwheel> out;
    std::vector<std::pair<int, int>> degrees(d - 1);

    auto enumDegree = [&](auto&& self, int i) -> void {
        if (i == d - 1) {
            out.push_back(generateCartwheel(d, degrees, true));
            return;
        }
        for (size_t j = 0; j < cd.size(); ++j) {
            degrees[i] = cd[j];
            self(self, i + 1);
        }
    };
    enumDegree(enumDegree, 0);
    return out;
}

namespace {

long long ipow(long long base, int exp) {
    long long r = 1;
    for (int i = 0; i < exp; ++i) r *= base;
    return r;
}

// Line 5 of Algorithm B.3.2: reject unless `degrees` is its own lexicographically minimal
// rotation.  CARTWHEEL_DEGREES is strictly increasing, so comparing the pairs is the same as
// comparing their indices.
bool isMinimalRotation(const std::vector<std::pair<int, int>>& degrees) {
    const int d = static_cast<int>(degrees.size());
    for (int s = 1; s < d; ++s)
        for (int k = 0; k < d; ++k) {
            if (degrees[(k + s) % d] < degrees[k]) return false;
            if (degrees[k] < degrees[(k + s) % d]) break;
        }
    return true;
}

}  // namespace

long long digonIncidentWheelCount(int d) { return ipow(5, d - 1); }

Cartwheel digonIncidentWheelAt(int d, long long idx) {
    const auto& cd = cartwheelDegrees();
    std::vector<std::pair<int, int>> degrees(d - 1);
    for (int i = d - 2; i >= 0; --i) {
        degrees[i] = cd[idx % 5];
        idx /= 5;
    }
    return generateCartwheel(d, degrees, true);
}

long long wheelCandidateCount(int d) {
    long long total = 0;
    for (int j = 0; j < 5; ++j) total += ipow(5 - j, d - 1);
    return total;
}

bool wheelCandidateAt(int d, long long idx, Cartwheel& out) {
    const auto& cd = cartwheelDegrees();
    int j = 0;
    for (; j < 5; ++j) {
        const long long n = ipow(5 - j, d - 1);
        if (idx < n) break;
        idx -= n;
    }
    const int base = 5 - j;
    std::vector<std::pair<int, int>> degrees(d);
    degrees[0] = cd[j];
    for (int i = d - 1; i >= 1; --i) {
        degrees[i] = cd[j + idx % base];
        idx /= base;
    }
    if (!isMinimalRotation(degrees)) return false;
    out = generateCartwheel(d, degrees, false);
    return true;
}

std::vector<Cartwheel> enumWheels(int d) {
    const auto& cd = cartwheelDegrees();
    std::vector<Cartwheel> out;
    std::vector<std::pair<int, int>> degrees(d);

    auto enumDegree = [&](auto&& self, int i, int iLowest) -> void {
        if (i == d) {
            // Keep only the lexicographically minimal rotation of `degrees`.
            for (int s = 1; s < d; ++s) {
                for (int k = 0; k < d; ++k) {
                    const auto& a = degrees[(k + s) % d];
                    const auto& b = degrees[k];
                    if (a < b) return;
                    if (b < a) break;
                }
            }
            out.push_back(generateCartwheel(d, degrees, false));
            return;
        }
        for (size_t j = iLowest; j < cd.size(); ++j) {
            degrees[i] = cd[j];
            self(self, i + 1, iLowest);
        }
    };

    for (size_t j = 0; j < cd.size(); ++j) {
        degrees[0] = cd[j];
        enumDegree(enumDegree, 1, static_cast<int>(j));
    }
    return out;
}

}  // namespace apex
