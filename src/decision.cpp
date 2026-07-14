// ============================================================
// decision.cpp -- linear-time BFB-count-vector decision algorithm
//
// A direct port of BFBCalculator.decisionBFB / decisionFold /
// hasPalindromeConcatenation from the reference Java implementation:
//   Zakov, Kinsella, Bafna. "An algorithmic approach for
//   breakage-fusion-bridge detection in tumor genomes." PNAS 2013.
//   https://github.com/shay-zakov/BFBFinder (GPL-3.0)
//
// The algorithm tracks only an integer *signature* array -- it never
// constructs an actual BFB string or palindrome-collection object --
// which is what makes it O(N) instead of exponential: unlike
// BFB-Pivot/BFB-Tree (bfbtools package), there is no backtracking or
// branching at all; the signature update at each position is a fixed,
// constant-amortized-time computation.
//
// See Signature.java's class documentation (ported package) for the
// formal definition of a signature and its role in deciding whether a
// palindrome collection can be folded to a given size.
// ============================================================

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// Number of trailing zero bits of m (largest d such that 2^d | m).
static inline int dig(int m) {
  int d = 0;
  while (m % 2 == 0) {
    m /= 2;
    ++d;
  }
  return d;
}

// A generous upper bound on the signature length ever needed, given the
// largest count value present (mirrors Env.getMaxR).
static inline int getMaxR(const std::vector<int>& counts) {
  int maxValue = 1;
  for (int c : counts) if (c > maxValue) maxValue = c;
  return (int)std::ceil(std::log((double)(maxValue + 1)) / std::log(2.0)) + 3;
}

// Updates signature `s` (length r(B)+1, representing an l-block
// collection B of size prevCount) in place to become the signature of
// the (unique, lexicographically-minimal) collection obtained by
// folding B to size n -- or returns -1 if no such folding exists.
// Direct port of BFBCalculator.decisionFold; see that method's
// docstring (and the paper's SI, Algorithm SIGNATURE-FOLD) for the
// mathematical derivation -- this is a faithful line-by-line port, not
// a reimplementation from the proof, specifically to minimize risk of
// transcription bugs in genuinely intricate index arithmetic.
static int decisionFold(std::vector<int>& s, int prevCount, int r, int n) {
  if (prevCount != n) {
    int Delta_d = prevCount;
    int factor;
    int d;
    int digVal = dig(prevCount - n);

    if (r < digVal) {
      factor = 1 << digVal;
      for (d = r + 1; d <= digVal; ++d) s[d] = 0;
      --d;
    } else {
      for (d = r + 1, factor = 1 << (r + 1); d > digVal; --d, factor >>= 1) {
        Delta_d -= (factor / 2) * std::abs(s[d - 1]);
      }
    }

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

        if (n >= Delta_d + deltaDif) {
          s[d + 1] = (Delta_d + deltaDif - n) / (2 * factor);
        } else {
          s[d] = (Delta_d - n) / factor;
          s[d + 1] = 0;
        }

        for (r = d + 1; s[r] == 0 && s[r - 1] <= 0; --r) ;
        return r;
      } else {
        return -1;  // no folding to size n exists
      }
    }
  } else {
    return r;  // no change needed
  }
}

// Decides whether the (implicit) collection with the given signature
// and size n can be concatenated into a single palindrome (the final
// step of the decision algorithm, checked once after processing the
// whole count vector). Direct port of hasPalindromeConcatenation.
static bool hasPalindromeConcatenation(const std::vector<int>& s, int n, int r) {
  if (n == 1 || s[0] == 0) return true;
  if (s[0] > 1) return false;
  int i = 1;
  for (; i <= r && s[i] == 0; ++i) ;
  return s[i] <= 0;
}

//' Linear-time BFB count-vector decision
//'
//' Direct port of the linear-time algorithm of Zakov, Kinsella & Bafna
//' (PNAS 2013) for deciding whether a count-vector admits a BFB
//' schedule. Unlike bfbtools::bfb_tree()/admits_bfb_schedule() (which
//' implement the earlier, worst-case-exponential algorithms of
//' Kinsella & Bafna 2012), this algorithm is provably O(N) time, with
//' no backtracking or search budget ever needed.
//'
//' Note on array order: the reference Java implementation's internal
//' `counts[l]` array is processed centromere-to-telomere (it folds
//' counts[k-1] first, counts[0] last), the reverse of the
//' telomere-to-centromere convention used throughout bfbtools and in
//' Kinsella & Bafna (2012). This function reverses the input before
//' calling the ported algorithm, so callers can pass count-vectors in
//' the same (telomere-to-centromere) order used everywhere else.
//'
//' @param counts_r integer count-vector (all entries >= 1), telomere
//'   to centromere.
//' @return TRUE if the count-vector admits a BFB schedule, FALSE
//'   otherwise.
// [[Rcpp::export]]
bool decision_bfb_cpp(IntegerVector counts_r) {
  std::vector<int> counts = as<std::vector<int>>(counts_r);
  std::reverse(counts.begin(), counts.end());
  int k = (int)counts.size();
  int maxR = getMaxR(counts);
  std::vector<int> s(maxR + 4, 0);  // generous padding around Env.getMaxR's own margin

  int r = 0;
  int prevCount = 0;
  bool res = true;

  for (int l = k - 1; l >= 0; prevCount = counts[l], --l) {
    r = decisionFold(s, prevCount, r, counts[l]);
    if (r < 0) {
      res = false;
      break;
    }
  }

  if (res) {
    res = hasPalindromeConcatenation(s, counts[0], r);
  }
  return res;
}
