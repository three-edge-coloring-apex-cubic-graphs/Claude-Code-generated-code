# Instructions to write source code to check "Three-edge coloring apex cubic graphs"

In the paper "Three-edge coloring apex cubic graphs" (`apex.pdf`), several lemmas are proved with the aid of computer calculations.

To check our computational work, the goal is to write source code corresponding to the pseudocode in Appendix B.1 to B.4 in the paper, together with the algorithms in [IKM+](https://arxiv.org/pdf/2603.24880) that are cited and required by our paper and use the resulting program to verify Lemma B.1, B.2, and B.3 in our paper.
You can use any programming language.

We provide files descrbing the configurations, rules, and auxiliary rules to verify the above lemmas in the following directories.

* Configurations: `configurations/K`
* Rules and auxiliary rules: `discharging-rules/R`, `discharging-rules/R_auxiliary`

The file formats are described in `FORMAT.md`.
Functions to parse these files are required to read each of these file types.


There are several additional remarks:

* Regarding semi-reducible checks, we provide pseudocode only for the newly introduced parts that differ from the existing implementation.
To perform the complete semi-D- and semi-C-reducible checks, one must additionally implement the standard iterative procedure used in existing reducibility checkers to compute the maximal semi-consistent subset of the nonextendable colorings.
* When implementing the algorithm for computing free homomorphic images, the procedure may terminate immediately if the intersection of the degree ranges within any merged vertex class is empty, as specified in Algorithm A.4.2 of [IKM+](https://arxiv.org/pdf/2603.24880). This early-termination rule improves efficiency without affecting the output of the algorithm.
* When verifying Lemma B.3, configurations in $\mathcal{K}$ must be processed in lexicographic order of their filenames. The configuration files are named `K001.conf`, `K002.conf`, $\ldots$, `K915.conf`. The reported value of $|\mathcal{I}|$ in `README.md` was computed using this ordering.
* To reproduce the numbers of multi-boundary islands reported in `README.md` for Lemma B.3, implement the pseudocode in the paper literally. In particular, do not introduce optimizations that change the order of processing, eliminate duplicate outputs, or otherwise alter the enumeration performed by the pseudocode. Such optimizations may preserve the final mathematical conclusion while producing different intermediate counts from those in `README.md`.
