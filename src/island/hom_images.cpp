#include "island/hom_images.hpp"

#include "core/walks.hpp"
#include "island/make_outer_extension.hpp"
#include "island/planarity.hpp"

namespace apex {

std::vector<MultiBoundaryIsland> allHomImages(const Graph& khat, const ConfSet& ksmaller) {
    std::vector<MultiBoundaryIsland> l;
    if (blockedByReducibleConfiguration(khat, NIL, ksmaller)) return l;

    if (isPlanar(khat) && !hasSeparatingCycle(khat)) {
        const Graph s = freeCompletionFromOuterExtension(khat);
        MultiBoundaryIsland island = islandFromFreeCompletion(s);
        if (island.numPendant <= 3) l.push_back(std::move(island));
    }

    std::vector<char> checked(khat.nd(), 0);
    for (int e = 0; e < khat.nd(); ++e) {
        for (int f = 0; f < khat.nd(); ++f) {
            if (e == f || checked[e] || checked[f]) continue;
            for (auto& img : freeHomomorphismAndEnforceSingleDigonIncidence(
                     khat, {{e, f}}, Boundary::PseudoEmbedding)) {
                for (auto& shat : makeOuterExtension(img.g)) {
                    auto m = allHomImages(shat.g, ksmaller);
                    for (auto& isl : m) l.push_back(std::move(isl));
                }
            }
        }
        checked[e] = 1;
        checked[khat.rev[e]] = 1;
    }
    return l;
}

}  // namespace apex
