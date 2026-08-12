#include "hom/blocked.hpp"

#include <algorithm>

#include "hom/hom.hpp"

namespace apex {
namespace {

// Algorithms A.6.6 + A.6.7 fused.  Algorithm A.6.7 buckets the darts of the target by their
// endpoint degrees so that each configuration only has to try the darts that match its special
// dart; here the same pairing is done from the other side, since the configuration set is fixed
// across millions of targets while each target is seen once.  `z` must carry a single concrete
// degree per vertex (dlo == dhi).
bool containConf(const Graph& z, int center, const ConfSet& k) {
    for (int fs = 0; fs < z.nd(); ++fs) {
        const int dy = z.dlo[z.head[fs]];
        const int dx = z.dlo[z.tail(fs)];
        if (dy > k.confDegMax || dx > k.confDegMax) continue;
        // A vertex of a configuration with delta > 8 may only map to the center.
        const bool centerOnly = (center != NIL && z.head[fs] != center);
        for (const int ci : k.bucket(dy, dx)) {
            if (centerOnly && dy > 8) continue;
            const auto& conf = k.confs[ci];
            if (homomorphismExists(conf.g, conf.dart, z, fs, gInclude)) return true;  // A.6.8
        }
    }
    return false;
}

// Algorithm A.7.2 fused with Algorithm A.7.1 so that a non-blocked representative aborts the
// enumeration immediately instead of materialising the whole product.
bool allRepresentativesBlocked(Graph& z, const std::vector<int>& lo, const std::vector<int>& hi,
                               int v, int center, const ConfSet& k) {
    if (v == z.nv) return containConf(z, center, k);

    std::vector<int> choices;
    const bool isCenter = (v == center);
    if (isCenter && hi[v] > k.confDegMax) {
        choices = {hi[v]};
    } else if (!isCenter && hi[v] > 8) {
        choices = {hi[v]};
    } else {
        for (int d = lo[v]; d <= hi[v]; ++d) choices.push_back(d);
    }
    for (const int d : choices) {
        z.dlo[v] = d;
        z.dhi[v] = d;
        if (!allRepresentativesBlocked(z, lo, hi, v + 1, center, k)) return false;
    }
    z.dlo[v] = lo[v];
    z.dhi[v] = hi[v];
    return true;
}

}  // namespace

ConfSet makeConfSet(std::vector<ConfEntry> confs) {
    ConfSet s;
    s.confs = std::move(confs);
    for (const auto& c : s.confs)
        for (int v = 0; v < c.g.nv; ++v)
            if (c.g.dlo[v] == c.g.dhi[v]) s.confDegMax = std::max(s.confDegMax, c.g.dlo[v]);

    const int n = s.confDegMax + 1;
    s.byDeg.assign(static_cast<size_t>(n) * n, {});
    for (size_t i = 0; i < s.confs.size(); ++i) {
        const auto& c = s.confs[i];
        const int dy = c.g.dlo[c.g.head[c.dart]];
        const int dx = c.g.dlo[c.g.tail(c.dart)];
        if (dy > s.confDegMax || dx > s.confDegMax) continue;
        s.byDeg[dy * n + dx].push_back(static_cast<int>(i));
    }
    return s;
}

bool blockedByReducibleConfiguration(const Graph& z, int center, const ConfSet& k) {
    if (k.confs.empty()) return false;
    Graph work = z;
    return allRepresentativesBlocked(work, z.dlo, z.dhi, 0, center, k);
}

}  // namespace apex
