// ============================================================
// distance.cpp -- linear-time-core weighted nearest-BFB-vector search
//
// Finds the highest-weight BFB count-vector "explaining" a noisy
// observed count-vector, under a probabilistic error model -- the
// noisy/distance variant of the BFB count-vector problem (Zakov,
// Kinsella & Bafna PNAS 2013; error models per Zakov & Bafna, JCB 2015,
// "Reconstructing breakage fusion bridge architectures using noisy
// copy numbers").
//
// Design note: the reference implementation's Signature.heaviestBFBVector
// additionally searches over variable-length sub-vectors (for finding
// the longest BFB-consistent *region* of a chromosome), using a
// parallel signature-update primitive (Signature.minIncrement) that
// checks palindrome-concatenation validity at *every* processed
// position -- needed to know, for each candidate sub-length, whether
// stopping there would already be valid.
//
// This port instead answers a narrower question matching
// bfbtools::bfb_distance's existing semantics: given a FIXED-length
// count-vector, find the closest (highest-weight) full-length BFB
// count-vector. For that question, validity only needs checking once,
// at the very end -- exactly the structure already implemented and
// validated in decision.cpp (0 mismatches across 20,000 random vectors
// against bfbtools, exact matches on all the paper's worked examples).
// Reusing that validated core here, rather than re-porting a second,
// structurally different implementation of the same combinatorics,
// substantially lowers transcription risk.
//
// At each position (processed centromere-to-telomere internally, same
// as decision.cpp), a *frontier* of candidate (signature, weight) pairs
// is maintained, explored bidirectionally around the observed count
// (always the maximum-weight candidate, by construction of the error
// models below). Dominated candidates -- lexicographically-higher
// signature AND lower-or-equal weight than some other frontier member
// -- are pruned, following the same Pareto-frontier argument as
// Signature.java's extend()/dominance pruning, which is what keeps the
// frontier size (and thus runtime) bounded in practice.
// ============================================================

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
using namespace Rcpp;

// ---- shared primitives (mirrors decision.cpp, with one change: the
// signature vector auto-grows as needed instead of relying on a
// pre-sized buffer, since the noisy search can explore candidate values
// well beyond the observed counts and a fixed size heuristic would risk
// silent out-of-bounds access). Kept file-local to avoid a
// cross-translation-unit header for this small amount of shared code.
// -----------------------------------------------------------------------
static inline int dig(int m) {
  int d = 0;
  while (m % 2 == 0) { m /= 2; ++d; }
  return d;
}

static inline void ensureSize(std::vector<int>& s, int idx) {
  if (idx >= (int)s.size()) s.resize(idx + 8, 0);  // grow with headroom
}

static int decisionFold(std::vector<int>& s, int prevCount, int r, int n) {
  if (prevCount != n) {
    int Delta_d = prevCount;
    int factor;
    int d;
    int digVal = dig(prevCount - n);

    ensureSize(s, std::max(r + 1, digVal) + 2);

    if (r < digVal) {
      factor = 1 << digVal;
      for (d = r + 1; d <= digVal; ++d) s[d] = 0;
      --d;
    } else {
      for (d = r + 1, factor = 1 << (r + 1); d > digVal; --d, factor >>= 1) {
        Delta_d -= (factor / 2) * std::abs(s[d - 1]);
      }
    }

    ensureSize(s, d + 2);

    if (prevCount < n) {
      ++s[d];
      s[d + 1] = (Delta_d + factor * std::abs(s[d]) - n) / (2 * factor);
      return d + 1;
    } else {
      for (; d > 0 && n < factor * std::max(s[d] + 1, 0) + Delta_d; --d, factor >>= 1) {
        Delta_d -= (factor >> 1) * std::abs(s[d - 1]);
      }
      if (d > 0 || n >= factor * std::max(s[d] + 1, 0) + Delta_d) {
        if (d == digVal) ++s[d];
        else s[d] += 2;

        int deltaDif = factor * std::abs(s[d]);
        ensureSize(s, d + 2);

        if (n >= Delta_d + deltaDif) {
          s[d + 1] = (Delta_d + deltaDif - n) / (2 * factor);
        } else {
          s[d] = (Delta_d - n) / factor;
          s[d + 1] = 0;
        }

        for (r = d + 1; s[r] == 0 && s[r - 1] <= 0; --r) ;
        return r;
      } else {
        return -1;
      }
    }
  } else {
    return r;
  }
}

static bool hasPalindromeConcatenation(const std::vector<int>& s, int n, int r) {
  if (n == 1 || s[0] == 0) return true;
  if (s[0] > 1) return false;
  int i = 1;
  for (; i <= r && s[i] == 0; ++i) ;
  return i < (int)s.size() && s[i] <= 0;
}

// ---- error models -----------------------------------------------------
// weight(trueCount, estimatedCount): the (unnormalized) plausibility of
// a candidate true count, given the observed estimate. Both models are
// maximized (weight == 1) exactly at trueCount == estimatedCount, which
// is what lets the search below start at the observed value and
// explore outward in both directions until weight decays away.

// Direct port of PoissonErrorModel: models observed counts as Poisson
// draws around a (trueCount + 0.5) mean (the +0.5 avoids a zero mean at
// trueCount == 0 and matches the reference implementation exactly),
// scored relative to the probability of observing the mode.
static const double HALF_LOG_PI = std::log(M_PI) / 2.0;
static const double factLans[] = {
  0, 0, 0.6931471805599453, 1.791759469228055, 3.1780538303479458,
  4.787491742782046, 6.579251212010101, 8.525161361065415,
  10.60460290274525, 12.801827480081469, 15.104412573075516, 17.502307845873887
};
static const int computedLans = 12;

static double factorialApprox(int n) {
  if (n < computedLans) return factLans[n];
  double a = n * std::log((double)n) - n;
  double b = std::log(n * (1.0 + 4.0 * n * (1.0 + 2.0 * n))) / 6.0;
  return a + b + HALF_LOG_PI;
}

static double poissonLogProb(double mean, int observed) {
  return observed * std::log(mean) - mean - factorialApprox(observed);
}

// Returns log-weight directly (avoids exp() then log() round-trips,
// and is what the search actually accumulates).
static double poissonLogWeight(int trueCount, int estimatedCount) {
  return poissonLogProb(trueCount + 0.5, estimatedCount) -
         poissonLogProb(trueCount + 0.5, trueCount);
}

// ---- frontier of Pareto-optimal (signature, weight) candidates --------
struct Candidate {
  std::vector<int> sig;
  int r;
  int prevCount;           // == cardinality of sig
  double logWeight;
  std::vector<int> counts; // chosen values so far, in processing order
};

static int sigCompare(const std::vector<int>& s1, int r1,
                       const std::vector<int>& s2, int r2) {
  int maxR = std::max(r1, r2);
  for (int d = 0; d <= maxR; ++d) {
    int v1 = (d < (int)s1.size() && d <= r1) ? s1[d] : 0;
    int v2 = (d < (int)s2.size() && d <= r2) ? s2[d] : 0;
    if (v1 != v2) return v1 - v2;
  }
  return 0;
}

// Inserts `cand` into the sorted-by-signature, Pareto-pruned frontier,
// mirroring Signature.java's extend()/dominance-pruning logic: a lower
// (or equal) signature with a higher (or equal) weight dominates.
// Returns false if the frontier's hard size cap was hit (caller should
// treat this as a budget-exceeded condition).
static bool insertPareto(std::vector<Candidate>& frontier, Candidate&& cand,
                          int maxFrontier) {
  int lo = 0, hi = (int)frontier.size();
  while (lo < hi) {
    int mid = (lo + hi) / 2;
    int c = sigCompare(cand.sig, cand.r, frontier[mid].sig, frontier[mid].r);
    if (c < 0) hi = mid; else if (c > 0) lo = mid + 1;
    else { lo = mid; break; }
  }
  int rank = lo;

  bool exactMatch = (rank < (int)frontier.size() &&
                      sigCompare(cand.sig, cand.r, frontier[rank].sig, frontier[rank].r) == 0);

  if (exactMatch) {
    if (frontier[rank].logWeight < cand.logWeight) frontier[rank] = std::move(cand);
    return true;
  }
  if (rank > 0 && frontier[rank - 1].logWeight >= cand.logWeight) {
    return true;  // dominated by the preceding (lower-signature) entry
  }
  if ((int)frontier.size() >= maxFrontier) {
    return false;  // budget exceeded
  }

  frontier.insert(frontier.begin() + rank, std::move(cand));

  int j = rank + 1;
  while (j < (int)frontier.size() && frontier[j].logWeight <= frontier[rank].logWeight) {
    frontier.erase(frontier.begin() + j);
  }
  return true;
}

// ---- top-level search ---------------------------------------------------
// Explores candidate true-counts m around each observed count,
// bidirectionally, stopping each direction once its per-step weight
// contribution decays below `minLogWeight` (cumulative, relative to the
// best-so-far candidate feeding into it). model 0 = NoErrorModel (exact
// match only -- reduces to decision_bfb_cpp's single deterministic
// path); model 1 = PoissonErrorModel.
//
// Returns a list with $weight, $vector -- or R_NilValue if nothing was
// found (can happen if minLogWeight is too strict), or a list with
// $budget_exceeded = TRUE if the frontier safety cap was hit before the
// search completed.
// [[Rcpp::export]]
List nearest_bfb_weighted_cpp(IntegerVector counts_r, int model, double minLogWeight,
                               int maxFrontier) {
  std::vector<int> counts = as<std::vector<int>>(counts_r);
  std::reverse(counts.begin(), counts.end());  // match internal processing order
  int k = (int)counts.size();

  std::vector<Candidate> frontier;
  Candidate init;
  init.sig.assign(4, 0);
  init.r = 0;
  init.prevCount = 0;
  init.logWeight = 0.0;
  frontier.push_back(std::move(init));

  for (int l = k - 1; l >= 0 && !frontier.empty(); --l) {
    int observed = counts[l];
    std::vector<Candidate> nextFrontier;

    for (const Candidate& c : frontier) {
      // Direction 0: m = observed, observed+1, observed+2, ...
      // Direction 1: m = observed-1, observed-2, ..., down to 1.
      for (int dir = 0; dir < 2; ++dir) {
        if (model == 0 && dir == 1) break;  // NoErrorModel: only m == observed

        int m = (dir == 0) ? observed : observed - 1;
        for (; m >= 1; m += (dir == 0 ? 1 : -1)) {
          double logW = (model == 1) ? poissonLogWeight(m, observed) : 0.0;
          double newLogWeight = c.logWeight + logW;
          if (newLogWeight < minLogWeight) break;  // weight decays monotonically outward

          std::vector<int> newSig = c.sig;
          int newR = decisionFold(newSig, c.prevCount, c.r, m);
          if (newR >= 0) {
            Candidate nc;
            nc.sig = std::move(newSig);
            nc.r = newR;
            nc.prevCount = m;
            nc.logWeight = newLogWeight;
            nc.counts = c.counts;
            nc.counts.push_back(m);
            if (!insertPareto(nextFrontier, std::move(nc), maxFrontier)) {
              return List::create(Named("budget_exceeded") = true);
            }
          }

          if (model == 0) break;  // single value only
        }
      }
    }

    frontier = std::move(nextFrontier);
  }

  double bestLogWeight = -std::numeric_limits<double>::infinity();
  const Candidate* best = nullptr;
  for (const Candidate& c : frontier) {
    if (hasPalindromeConcatenation(c.sig, c.prevCount, c.r) && c.logWeight > bestLogWeight) {
      bestLogWeight = c.logWeight;
      best = &c;
    }
  }

  if (best == nullptr) return R_NilValue;

  // best->counts was built by push_back() during processing (l = k-1
  // downto 0 over the internally-reversed array), which -- given the
  // single reversal applied on input -- already reconstructs the
  // original telomere-to-centromere order. (Verified directly: with
  // counts reversed once, counts[l] for l = k-1, k-2, ..., 0 steps
  // through the original array in its original order, c0, c1, ..., so
  // appending in processing order already yields [c0, c1, ..., ck-1].)
  return List::create(
    Named("weight") = std::exp(bestLogWeight),
    Named("vector") = wrap(best->counts)
  );
}
