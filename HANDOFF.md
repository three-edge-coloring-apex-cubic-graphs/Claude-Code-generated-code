# Implementation notes

Lemmas B.1, B.2, and B.3 of `apex.pdf` are all fully reproduced (see `README.md` for the
top-level status table, build/run instructions, and source map). This file is the detailed
technical record of getting there: the two ambiguities in the pseudocode/file formats that
Lemma B.3's island construction needed resolving, the derivation of the semi-D/C-reducibility
fixpoint procedure (unpublished anywhere, per `.claude/CLAUDE.md`), and the bugs -- one
correctness, one performance -- found and fixed while verifying it end-to-end.

## Ground rules carried over from the task

* Implement the pseudocode of Appendix B of `apex.pdf` and the algorithms of Appendix A of
  [IKM+26](https://arxiv.org/abs/2603.24880) that it cites, then check Lemma B.1, B.2, B.3.
* **The original source code of [IKM+26] must not be consulted.** Do not fetch
  `github.com/edge-coloring/reducibility_checker` or `github.com/near-linear-4ct/computer-checks`,
  which the paper links. Everything is written from the published pseudocode only.
* Lemma B.3 must process `configurations/K/K001.conf` ... `K915.conf` in lexicographic filename
  order, and must not reorder, deduplicate, or otherwise optimise the enumeration itself —
  those choices change the intermediate counts in `Metrics.md` even when the mathematical
  conclusion is unaffected. (Caching *reducibility verdicts* by a canonical form of an island is
  fine; caching that skips enumeration is not.)

## Setting up on a new machine

`.gitignore` excludes the inputs, so a clone is not enough:

```sh
git clone git@github.com:kappybar/apex-by-ClaudeCode.git && cd apex-by-ClaudeCode
# then place, from wherever you keep them:
#   apex.pdf                  the paper
#   configurations/           from git@github.com:kappybar/apex_confs.git   (has K/)
#   discharging-rules/        from git@github.com:kappybar/apex_confs.git   (has R/, R_auxiliary/)
curl -L -o ikm26.pdf https://arxiv.org/pdf/2603.24880      # reference only, gitignored

cmake -B build -DCMAKE_BUILD_TYPE=Release -S . && cmake --build build -j$(nproc)
ctest --test-dir build --output-on-failure                 # ~1 s: parsers + Lemma B.1
```

Sanity check that the inputs are in place: `test_parse` must report 915 configurations,
70 rules, 4 auxiliary rules, and `verify_b1` must print `|K-bar| = 2150`.

## What is done

Lemmas B.1, B.2, and B.3 are all fully reproduced, including B.3's semi-D/semi-C-reducibility
verification (B.1.1-B.1.5 plus the unpublished fixpoint procedure) over every one of the 109 501
islands in `I` — see "Lemma B.3 semi-reducibility verification" below for the derivation, a
critical correctness bug found and fixed during testing, and a critical performance bug found
and fixed afterward (37x speedup, needed to make the full sweep tractable).

| Lemma | Metric | Expected | Reproduced |
| --- | --- | --- | --- |
| B.1 | \|R*^-K\| | 747 | **747** |
| B.1 | max charge of a combined rule | 5 | **5** |
| B.2 | `enumPossibleBadWheels`, d = 7 | 4438 | **4438** |
| B.2 | `enumPossibleBadWheels`, d = 8 | 4939 | **4939** |
| B.2 | `enumPossibleBadWheels`, d = 9 | 2409 | **2409** |
| B.2 | `enumPossibleBadWheels`, d = 10 | 567 | **567** |
| B.2 | `enumPossibleBadWheels`, d = 11 | 38 | **38** |
| B.2 | `verifyNoBadCartwheels`: C = ∅, d = 7..11 | — | **holds at every degree** |
| B.3 | multi-boundary islands by \|F_R(I)\| | 254/88393/20836/18/0 | **254/88393/20836/18/0** |
| B.3 | every island semi-D- or semi-C-reducible | 109 501 / 109 501 | **109 501 / 109 501** (0 failures) |

`./build/verify_b2` (no flags) runs the full sweep in under two minutes on 256 cores. Matching
d = 10 and d = 11 also resolves the `kDegreePairs` ambiguity noted previously: the `[5,∞]`
reading of Algorithm B.3.5's third row (not the literal-8 `pdftotext` rendering) is confirmed
correct, since the literal-8 reading is now ruled out by two more data points, not just d = 7.

`verifyNoBadCartwheels`'s cost tracks the *surviving* wheel count at each degree
(4438/4939/2409/567/38), not the number of candidates `enumPossibleBadWheels` scans to find
them — so it stayed cheap (11.0/20.4/6.4/1.5/0.2 s) even at d = 11, where scanning alone takes
34 s. The memoisation lever mentioned in an earlier draft of this file was not needed.

## Lemma B.3 semi-reducibility: two full-scale runs, in agreement

`./build/verify_b3_reduce` (see below) implements island construction and the semi-D/C-reducibility
check in one pass, decoupled via a producer/consumer queue so construction and checking proceed
independently (see the performance notes below). Two full end-to-end runs over all 915
configurations have completed:

* An unoptimised run (before the `componentCandidates` caching fix below) took 450 649.9s
  (~5.2 days), run under `screen` for resilience across harness session restarts.
* A second, independent run with the caching fix applied was launched in parallel (deliberately
  kept running alongside the first rather than replacing it, to get an independent cross-check)
  and completed in 299 188.4s (~3.46 days).

Both runs agree exactly: **109 501/109 501 islands semi-D- or semi-C-reducible, 0 failures**
(94 333 semi-D, 15 168 semi-C in both), island-construction counts matching `Metrics.md` exactly
(254/88393/20836/18/0). This fully closes out Lemma B.3. The second run's aggregate speedup over
the first (~1.5x) is much smaller than the 37x measured on the single worst island in isolation
(see below) — most of both runs' wall time is dominated by the same handful of the largest
islands (the two ring-15 cases in particular), which even after the fix still take a long time in
absolute terms since they're single-threaded once the work queue narrows down to them; the fix's
main effect was making the run *tractable* (days instead of an unbounded multi-week tail), not
uniformly fast throughout.

## Lemma B.3 semi-reducibility verification — implemented, one critical bug found and fixed

`./build/verify_b3_reduce` builds `I` exactly as `verify_b3` does (so it also reproduces the
254/88393/20836/18/0 island-construction counts), and for each island runs `checkReducibility`
(`src/reduce/semi_reducible.cpp`), which needs B.1.1-B.1.5 (planar half-Kempe chains, deletable
edge sets — `src/reduce/kempe.*`, `deletable.*`, `island_graph.*`, `coloring.*`) plus one piece
that **is not published anywhere**: the iterative computation of the maximal semi-consistent
subset of the non-extendable colourings, derived from Definitions 3.5/3.6/3.9 of `apex.pdf`:

> Start with `S = D_R \ C_I`. Repeatedly delete any `φ ∈ S` for which some colour pair `{x,y}`
> admits **no** semi-matching `M` partitioning `{r : φ(r) ∈ {x,y}}` with every Kempe-switch of
> every `M' ⊆ M` landing back in `S`; iterate to a fixpoint. `I` is semi-D-reducible iff the
> fixpoint is empty, and semi-C-reducible by `F` iff the fixpoint is disjoint from `C_{I−F}`.

Candidate semi-matchings come from `GetPlanarHalfKempes(n)` (B.1.2), concatenated across
boundary components for multi-boundary islands — restricting the existential in Definition 3.5
to *only* `GetPlanarHalfKempes`'s "nonredundant" (geometrically maximal) family, rather than every
noncrossing partial matching, is necessary: searching over all noncrossing partial matchings
(including the trivial all-singleton one) was verified empirically to make the check vacuous —
nothing is ever removed from `S` — because a *correctly modelled* singleton is independently
flippable (see the bug below), and once every matching including the maximal ones is on the
table, the search space includes options for every φ that never fails. Ring sizes reach 15
(confirmed by a full survey: max single ring = max total ring-edge count across all boundary
components combined = 15; max ring count = 3; every 3-ring island has ring sizes exactly
`[2, 2, 2]`), so `3^r` colourings per island is the dominant cost; `checkReducibility` throws if
`3^(total ring edges)` exceeds a 60M-entry cap, which the survey confirms never triggers for `I`.

**Bug found and fixed: a semi-matching's singletons are independently flippable, not inert.**
Definition 3.3 makes a singleton `{r}` a full member of a semi-matching `M`, on equal footing
with its matches. Definition 3.4 then lets `M'` range over *all* subsets of `M` — matches *and*
singletons alike — and if a singleton `{r}` is included in `M'`, `r ∈ supp(M')`, so `r`'s own
colour flips too (there is no clause exempting singletons from the swap formula). The first
implementation treated singletons as inert padding (only matches were toggle-able), so a
candidate matching with `t` matches and `s` singletons was checked against `2^t` reachable
colourings instead of the correct `2^(t+s)` — silently skipping every combination that flips one
or more lone positions. This makes "survives" too easy to satisfy (fewer constraints checked than
Definition 3.5 requires), so the computed fixpoint comes out *too large* — in one confirmed case
(`K029.conf`, island #101: `n=22`, single ring of size 10) large enough that no candidate deletable
edge set among the 237 valid ones (verified complete and correct via an independent brute-force
re-check of Definition 3.8, and via cross-referencing against the user's own knowledge of this
example) could make it disjoint from `C_{I-F}`, even though a specific size-4 F
(`{16, 20, 29, 33}`, edges `(2,14) (4,17) (11,12) (15,16)`) is genuinely valid and was confirmed
by the user to work. The fix (`src/reduce/semi_reducible.cpp`'s `Unit`/`componentCandidates`
rework) treats every candidate matching's unmatched positions as additional singleton units,
each independently included/excluded in the `2^(t+s)`-way switch enumeration. Re-running the same
K029 case after the fix: 0 failures across the first 100 configurations, with K029's island #101
now correctly semi-C-reducible.

Finding this took an unusual amount of independent cross-validation, all of which came back
clean and pointed the remaining suspicion at the fixpoint's core semantics rather than any of its
supporting pieces — worth knowing if a similar discrepancy ever needs re-diagnosing:

* `getDeletableEdgeSet`: brute-forced Definition 3.8 from scratch (bypassing `Place`/`Recurse`
  entirely) on K029's island #101 — exact match on counts by size (28/332/2192/237).
* `computeCI`: independently re-verified (a second, differently-structured backtracking search)
  that a specific stuck boundary colouring has no direct completion — correctly absent from `C_I`.
* `computeCIModuloF`: independently re-verified that the same colouring *does* extend modulo
  `F={16,20,29,33}` — correctly present in `C_{I-F}`.
* `getPlanarHalfKempes`: cross-validated against an independently-coded brute-force "geometrically
  maximal noncrossing matching" enumerator for every `n` from 1 to 10 (the full range this island
  uses) — exact match every time.
* The fixpoint iteration itself: self-consistency checked (every element of the final computed
  fixpoint was re-verified to still have a witness against that same final set) — ruling out
  premature termination as the cause.

None of that found anything wrong, because the bug was in what "flippable" means, not in how any
individual piece was computed or combined.

Runtime is meaningfully higher after the fix (more switch combinations to check per candidate).
Per-island cost varies enormously with ring size and is *not* evenly spread across
configurations: most islands take well under a second, but occasional ones take 40+ minutes
single-threaded, and a configuration's few expensive islands can appear alongside many cheap
ones. An earlier version of `verify_b3_reduce` parallelised `checkReducibility` only within each
configuration's own island batch (`parallelFor`, matching `verify_b2.cpp`'s pattern) — this left
most of 256 cores idle whenever a batch's expensive islands hadn't finished but its cheap ones
had, since the next configuration's (already constructible) islands couldn't be picked up yet.
Decoupling construction from checking via a bounded producer/consumer queue (the main thread
constructs and enqueues islands as fast as it can; a fixed pool of worker threads drains the
queue continuously, drawing on however far construction has gotten rather than just the current
configuration's batch) measured a **3.5x wall-clock speedup** on the first 45 configurations
(9609.7s -> 2737.6s, identical results both times: 0 failures, semiD=119, semiC=236). Two
narrower micro-optimisations tried first (Gray-code incremental updates in the switch-mask loop
to avoid an O(t) rescan per mask; sorting each candidate's unit list so cheaper, more-likely-to-
survive candidates are tried first) measured as **no improvement** on their own -- the bottleneck
is core utilisation, not per-mask constant factors, so the queue is what actually matters; the
micro-optimisations were kept since they're free once correct, but don't expect much from
similar tuning without addressing utilisation first.

**`getValidParens`/`getPlanarHalfKempes` memoise into non-thread-safe function-local `static`
maps** (`src/reduce/kempe.cpp`) — call `warmKempeCaches(20)` once, single-threaded, before any
concurrent `checkReducibility` calls, or first-population races will corrupt results under
concurrency (this was mistaken for the K029 bug initially, before the singleton issue was found;
fixing it did not change K029's result, but it is still a real, necessary fix — verify_b3_reduce
already calls it).

### Second performance bug: componentCandidates rebuilt and re-sorted from scratch on every call

The queue redesign above fixed *core utilisation*; it did not fix the fact that a single large
island could still take hours on its own. Diagnosed by launching the real (unmodified)
`verify_b3_reduce` full sweep and, after it had been running for multiple days and had visibly
narrowed down to a handful of threads (most of 256 workers had already exited, meaning the queue
had drained to just the last few exceptionally expensive islands), attaching `gdb` to the live
process to sample worker call stacks — see the technique below, useful for any future "is it
actually stuck" question. That showed threads 22-39 frames deep in `computeCIModuloF`'s
backtracking colour search, which led to first suspecting (and fixing, harmlessly but with no
measured benefit) redundant enumeration in `src/reduce/coloring.cpp::search` once the ring
prefix is fixed. **That fix was real but not the bottleneck** -- a follow-up instrumented rebuild
that timed `computeCI`, each fixpoint round, and the semi-C fallback separately on the worst
island found live in the sweep (`K905.conf` island 0: `n=29`, ring 13, 6905.8s total) found:

```
computeCI:              67.0s
fixpoint round 1:     1765.3s   |active| 1,594,323 - 48,414 (C_I) -> 1,463,656
fixpoint round 2:     1608.0s   -> 1,255,430
fixpoint round 3:     1459.2s   -> 1,195,743
fixpoint round 4:     1408.9s   -> 1,195,743 (confirms convergence)
fixpoint total:      6241.3s   (90% of the wall time)
getDeletableEdgeSet:    0.1s   8367 candidates
semi-C search:         27.6s   found a working F on the *first* of 8367 tried
```

The fixpoint, not the semi-C fallback, is completely dominant. Within it,
`componentCandidates` (`src/reduce/semi_reducible.cpp`) rebuilds its whole candidate list --
including a full re-sort -- from scratch on every one of the millions of `(phi, x, y)` calls in
a round, even though the result depends only on `(offset, size, which local positions are
coloured {x, y})`, not on which of `x`/`y` each one is or which island produced the pattern. A
ring of size `r` has only `2^r` distinct local patterns, orders of magnitude fewer than the
`3^r * 3` colouring/colour-pair combinations that hit this function per round (for `r = 13`,
`2^13 = 8192` vs. ~4.6M calls in round 1 alone) -- so it is now cached globally (keyed on
`(offset, size, bitmask)`, mutex-guarded since `checkReducibility` runs concurrently across
islands), returning `shared_ptr<const vector<Unit>>` per candidate to keep cache hits allocation-
free. Re-measured on the same island: **6905.8s -> 186.6s, a 37x speedup**, with identical
correctness (same fixpoint size 1,195,743, same semi-C witness found). Also re-verified against
the full K001-K045 baseline: identical semiD=119/semiC=236/failed=0.

**Live-debugging technique, if a future run's progress needs checking without waiting for the
next log checkpoint:** `gdb -p <pid> -batch -ex "thread apply all bt 3"` samples every worker's
current call stack; comparing two samples taken with a real time gap between them (different
instruction-pointer addresses = genuine progress, not a hang) is a fast sanity check. Counting
distinct threads still present is also informative: as the bounded producer/consumer queue
drains near the end of a run, most of the 256 workers exit (queue empty + closed), so seeing only
a handful left is itself a sign the run has narrowed down to the last few hard islands, not that
something is wrong. Reading a specific live counter's exact value (e.g. `countChecked`) is
possible in principle via the binary's disassembly (locate the atomic instruction, resolve which
register/closure-offset it operates on, read it from a live thread's register/stack state) but
proved too fragile to rely on in a release build without debug symbols (register-unwind
information for deeply recursive frames was unreliable) -- rebuilding with `-g` would fix this,
at the cost of restarting whatever is currently running.

## Lemma B.3 island construction — done, matches exactly

```sh
./build/verify_b3              # ~2h50m single-threaded; run under `screen`/`tmux`, it is not
                                # parallelised and does not checkpoint
./build/verify_b3 --limit 30   # first N configurations only, for quick iteration
```

Result: 254 / 88393 / 20836 / 18 / 0 islands by `|F_R(I)|`, matching `Metrics.md` exactly
(109 501 total). Getting here needed two corrections beyond the pseudocode itself, both from
things the pseudocode/file-format genuinely don't state:

1. **The outer extension K-hat is not "the free completion with one corner left open."** Ring
   vertices have no `δ_K` anywhere in the `.conf` file format; the paper gives them the range
   `[5, ∞)` (Lemma 2.2's global minimum degree in G*), not a fixed value, so
   `δ_K^out(e) = δ_K(v) - d_{G(K)}(v)` can't be computed as a single number the way it can for
   internal vertices. The construction that reproduces the published counts: take the free
   completion, **delete every ring vertex**, then attach one fresh pendant vertex (range
   `[5, ∞)`, genuinely open) to each remaining internal vertex for every edge it lost to a
   deleted ring vertex — restoring it to its own already-exact `δ_K`. See
   `src/island/outer_extension.cpp` for the derivation in full; it does not depend on
   `freeCompletion()` or `reconstructRingRotations()` at all, since internal vertices' rotations
   are already given directly in the file.
2. **Algorithm B.4.6 `labelDarts`'s second branch was transcribed wrong from the PDF.** As
   captured, it reads `succ(pred(d_{i+1})) ... -> Propagate(succ(pred(d_{i+1})), R)`. `succ` and
   `pred` are inverses (M3 in Section 5), so `succ(pred(x)) = x` identically — the branch was
   propagating `R` starting from one of the cycle's own darts, and since `Propagate` sets
   `LD(e)` unconditionally on entry (no check for "already labelled"), it silently overwrote
   that dart's own label. Reading `pred` as an OCR mis-render of `rev` — symmetric with the
   first branch's `succ(d_i)` vs `rev(d_{i+1})` pattern — fixes it; see `src/island/planarity.cpp`.
   Symptom before the fix: `hasSeparatingCycle` fired on ordinary internal triangles/digons that
   aren't actually separating, which was catastrophic — it suppressed the vast majority of valid
   islands (observed totals under 4000 with the corrupted labelling, vs. 109 501 fixed).

Runtime is dominated by the (e, f) dart-pair search in `allHomImages` on the larger/later
configurations — some single configurations near the end of the K-set (e.g. K740-K750-ish) take
30-80s each. Not currently parallelised across configurations, because `Ksmaller` (and thus
`blockedByReducibleConfiguration`'s behaviour) must see exactly the configurations processed
strictly before the current one, in filename order (see the ground rules above) — parallelising
would need each config's call to snapshot `Ksmaller` at the right point rather than running
concurrently against a shared, growing one.

## Reusable pieces already in the tree

* `getWalks` / `isPlanar` (B.4.2, B.4.3) — `src/core/walks.cpp`
* `addBoundaryDartsDirectly` (B.2.4) and `linkIncidenceListEnds` (B.2.6) are exported from
  `src/hom/hom.hpp` precisely because Algorithm B.4.11 calls them directly
* `blockedByReducibleConfiguration` takes a **nullable** center, so the uncentered variant that
  B.4.1 requires is already available — pass `NIL`
* `freeHomomorphismAndEnforceSingleDigonIncidence` with `Boundary::PseudoEmbedding` is the
  variant B.4.1 and B.4.10 need (they operate on pseudo-embeddings, not triangulations)
* `canonicalKey(g, rootDart)` for deduplication / memoisation
* `fromVRotations` (B.2.10) and the `FORMAT.md` parsers
* `src/island/outer_extension.*`, `planarity.*`, `make_outer_extension.*`, `island.*`,
  `hom_images.*` — all of B.4.1-B.4.16, done and matching (see above)
* `src/reduce/kempe.*` (B.1.1-B.1.2), `deletable.*` (B.1.3-B.1.5), `island_graph.*` (reconstructs
  I's actual planar graph, with real degree-one leaves for ring/pendant edges, from the compact
  `MultiBoundaryIsland` encoding — needed for `deletable.cpp`'s face-boundary walks),
  `coloring.*` (`computeCI`/`computeCIModuloF`), `semi_reducible.*` (the fixpoint procedure and
  `checkReducibility`) — all of the semi-reducibility verification, done and correct (see above)

## Performance notes that are easy to undo by accident

Four optimisations in the B.2 path are load-bearing and all preserve the computed set exactly;
each is commented at its site. If a future change makes results drift, suspect these first:

1. `enumPossibleBadWheels` streams its two wheel families index-addressably rather than
   materialising `W ∪ W_digon` (d = 11 alone would need ~14 GB).
2. `prune` evaluates B.3.10's three predicates cheapest-first (they are a disjunction) and
   bounds the charge incrementally.
3. `neverApply` rejects on the root-pair degree intersection before building anything (82% of
   calls) and only asks whether *an* image exists, via
   `freeHomomorphismAndEnforceSingleDigonIncidenceAny`.
4. `homomorphismExists` and `freeHomTriangulation` reuse thread-local scratch buffers. These
   are **not reentrant** — do not call them recursively from within themselves.

`src/cartwheel/charge.hpp` exposes `takeCounters()`; `verify_b2` prints call counts and
prune reasons per degree, which is how the above were located. Keep it when profiling.
