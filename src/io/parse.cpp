#include "io/parse.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

#include "core/rotations.hpp"

namespace apex {
namespace {

struct Cursor {
    std::vector<std::vector<std::string>> lines;  // tokenised, blank lines kept as empty
    size_t i = 0;
    std::string path;

    bool atEnd() {
        while (i < lines.size() && lines[i].empty()) ++i;
        return i >= lines.size();
    }
    const std::vector<std::string>& next() {
        if (atEnd()) throw std::runtime_error("unexpected end of file: " + path);
        return lines[i++];
    }
    [[noreturn]] void fail(const std::string& what) const {
        throw std::runtime_error(path + ": " + what);
    }
};

Cursor readFile(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    Cursor c;
    c.path = path;
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream ss(line);
        std::vector<std::string> toks;
        std::string t;
        while (ss >> t) toks.push_back(t);
        c.lines.push_back(std::move(toks));
    }
    return c;
}

int toInt(const std::string& s) { return std::stoi(s); }

// A degree upper bound of 0 in a file means infinity.
int degHi(int x) { return x == 0 ? INF : x; }

// Reconstructs the rotation of each ring vertex of the free completion.  Every internal face of
// the free completion is a triangle, so for adjacent x, y the identity next_x(next_y(x)) = y
// holds; applying it in both directions around a ring vertex pins down its whole rotation, and
// the single remaining link (next(v_{i-1}) = v_{i+1}) is the corner at the outer face.
void reconstructRingRotations(Configuration& c) {
    std::vector<std::map<int, int>> nxt(c.R);
    std::vector<std::vector<int>> ringNbrs(c.R);

    for (int u = c.R; u < c.N; ++u) {
        const auto& rot = c.rotations[u];
        const int d = static_cast<int>(rot.size());
        for (int j = 0; j < d; ++j) {
            const int vi = rot[j];
            if (vi < 0 || vi >= c.R) continue;
            const int nextOfVi = rot[(j + 1) % d];
            const int prevOfVi = rot[(j - 1 + d) % d];
            ringNbrs[vi].push_back(u);
            nxt[vi][u] = prevOfVi;
            nxt[vi][nextOfVi] = u;
        }
    }

    for (int i = 0; i < c.R; ++i) {
        const int prevRing = (i - 1 + c.R) % c.R, nextRing = (i + 1) % c.R;
        nxt[i][prevRing] = nextRing;
        std::vector<int> rot;
        int cur = nextRing;
        for (size_t step = 0; step <= nxt[i].size(); ++step) {
            rot.push_back(cur);
            auto it = nxt[i].find(cur);
            if (it == nxt[i].end())
                throw std::runtime_error(c.name + ": incomplete ring rotation at vertex " +
                                         std::to_string(i));
            cur = it->second;
            if (cur == nextRing) break;
        }
        if (rot.size() != nxt[i].size())
            throw std::runtime_error(c.name + ": ring rotation is not a single cycle at vertex " +
                                     std::to_string(i));
        c.rotations[i] = rot;
    }
}

Rule parseRuleAt(Cursor& c, bool digonCountMandatory) {
    Rule r;
    r.name = c.path;
    const auto& hdr = c.next();
    if (hdr.size() != 4) c.fail("rule header must be 'N s t r'");
    const int N = toInt(hdr[0]);
    r.s = toInt(hdr[1]) - 1;
    r.t = toInt(hdr[2]) - 1;
    r.r = toInt(hdr[3]);

    std::vector<std::vector<int>> rotations(N);
    std::vector<int> lo(N), hi(N);
    for (int k = 0; k < N; ++k) {
        const auto& ln = c.next();
        if (ln.size() < 3) c.fail("short vertex line");
        const int idx = toInt(ln[0]) - 1;
        if (idx < 0 || idx >= N) c.fail("vertex index out of range");
        lo[idx] = toInt(ln[1]);
        hi[idx] = degHi(toInt(ln[2]));
        for (size_t j = 3; j < ln.size(); ++j) {
            const int a = toInt(ln[j]);
            rotations[idx].push_back(a == -1 ? NIL : a - 1);
        }
    }

    std::vector<std::pair<int, int>> digons;
    if (digonCountMandatory || !c.atEnd()) {
        const auto& ml = c.next();
        if (ml.size() != 1) c.fail("expected the digon count M");
        const int M = toInt(ml[0]);
        for (int k = 0; k < M; ++k) {
            const auto& dl = c.next();
            if (dl.size() != 2) c.fail("expected a digon endpoint pair");
            digons.emplace_back(toInt(dl[0]) - 1, toInt(dl[1]) - 1);
        }
    }

    r.g = fromVRotations(N, rotations, digons);
    for (int v = 0; v < N; ++v) {
        r.g.dlo[v] = lo[v];
        r.g.dhi[v] = hi[v];
    }
    for (int e = 0; e < r.g.nd(); ++e) {
        if (r.g.head[e] == r.t && r.g.tail(e) == r.s) {
            r.dart = e;
            break;
        }
    }
    if (r.dart == NIL) c.fail("no dart from s to t");
    return r;
}

}  // namespace

Configuration parseConfiguration(const std::string& path) {
    Cursor c = readFile(path);
    Configuration conf;
    conf.name = std::filesystem::path(path).filename().string();
    c.path = conf.name;

    const auto& hdr = c.next();
    if (hdr.size() != 2) c.fail("configuration header must be 'N R'");
    conf.N = toInt(hdr[0]);
    conf.R = toInt(hdr[1]);
    conf.rotations.assign(conf.N, {});
    conf.dlo.assign(conf.N, 0);
    conf.dhi.assign(conf.N, 0);

    std::vector<int> listedDegree(conf.N, 0);
    for (int k = conf.R; k < conf.N; ++k) {
        const auto& ln = c.next();
        if (ln.size() < 2) c.fail("short vertex line");
        const int idx = toInt(ln[0]) - 1;
        const int d = toInt(ln[1]);
        if (idx < conf.R || idx >= conf.N) c.fail("internal vertex index out of range");
        if (static_cast<int>(ln.size()) != d + 2) c.fail("neighbour count does not match d(v)");
        for (int j = 0; j < d; ++j) {
            const int a = toInt(ln[2 + j]);
            conf.rotations[idx].push_back(a == -1 ? NIL : a - 1);
        }
        listedDegree[idx] = d;
    }

    if (!c.atEnd()) {
        const auto& ml = c.next();
        if (ml.size() != 1) c.fail("expected the digon count M");
        const int M = toInt(ml[0]);
        for (int k = 0; k < M; ++k) {
            const auto& dl = c.next();
            if (dl.size() != 2) c.fail("expected a digon endpoint pair");
            conf.digons.emplace_back(toInt(dl[0]) - 1, toInt(dl[1]) - 1);
        }
    }
    if (!c.atEnd()) c.fail("trailing content");

    // delta_K(v) = d(v) + (number of incident digons).
    std::vector<int> digonCount(conf.N, 0);
    for (auto [a, b] : conf.digons) {
        ++digonCount[a];
        ++digonCount[b];
        if (a < conf.R || b < conf.R) c.fail("a digon endpoint lies on the ring");
    }
    for (int v = conf.R; v < conf.N; ++v)
        conf.dlo[v] = conf.dhi[v] = listedDegree[v] + digonCount[v];

    reconstructRingRotations(conf);
    return conf;
}

Graph freeCompletion(const Configuration& c) {
    Graph g = fromVRotations(c.N, c.rotations, c.digons);
    for (int v = 0; v < c.N; ++v) {
        if (v < c.R) {
            g.dlo[v] = dartDegree(g, v);
            g.dhi[v] = dartDegree(g, v);
        } else {
            g.dlo[v] = c.dlo[v];
            g.dhi[v] = c.dhi[v];
        }
    }
    return g;
}

Rule parseRule(const std::string& path) {
    Cursor c = readFile(path);
    c.path = std::filesystem::path(path).filename().string();
    Rule r = parseRuleAt(c, /*digonCountMandatory=*/false);
    r.name = c.path;
    return r;
}

AuxiliaryRule parseAuxiliaryRule(const std::string& path) {
    Cursor c = readFile(path);
    c.path = std::filesystem::path(path).filename().string();
    AuxiliaryRule a;
    a.base = parseRuleAt(c, /*digonCountMandatory=*/true);
    a.base.name = c.path;
    const auto& kl = c.next();
    if (kl.size() != 1) c.fail("expected the cover size k");
    const int k = toInt(kl[0]);
    for (int i = 0; i < k; ++i) {
        Rule sub = parseRuleAt(c, /*digonCountMandatory=*/true);
        sub.name = c.path + "#" + std::to_string(i + 1);
        a.cover.push_back(std::move(sub));
    }
    return a;
}

std::vector<std::string> listFiles(const std::string& dir, const std::string& suffix) {
    std::vector<std::string> out;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        const std::string p = entry.path().string();
        if (p.size() >= suffix.size() && p.compare(p.size() - suffix.size(), suffix.size(), suffix) == 0)
            out.push_back(p);
    }
    std::sort(out.begin(), out.end());
    return out;
}

}  // namespace apex
