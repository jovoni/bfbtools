#' bfbtools: Linear-Time BFB Count-Vector Algorithms
#'
#' Provides linear-time algorithms (ported from Zakov, Kinsella & Bafna
#' 2013, and the reference implementation at
#' \url{https://github.com/shay-zakov/BFBFinder}) for deciding whether a
#' count-vector admits a BFB schedule, and for finding the closest
#' BFB-admitting count-vector to a noisy observed one. Sim-agnostic: no
#' dependency on any particular simulator (see the sibling packages
#' \code{bfbgill} for Gillespie-based simulation and \code{bridges} for
#' fitting/plotting/batch-detection glue code). All algorithmic work
#' happens in C++; R functions here are thin argument-checking wrappers.
#'
#' @docType package
#' @name bfbtools
#' @useDynLib bfbtools, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"

#' Validate a count-vector's basic shape
#' @noRd
.bfbtools_check_counts <- function(n_vec) {
  n_vec <- suppressWarnings(as.integer(n_vec))
  if (length(n_vec) < 1 || anyNA(n_vec) || any(n_vec < 1)) {
    stop("n_vec must be a vector of one or more positive integers", call. = FALSE)
  }
  n_vec
}

#' Does a count-vector admit a BFB schedule? (linear time)
#'
#' Decides whether a count-vector \code{n_vec} (telomere to centromere)
#' can be produced by the breakage-fusion-bridge mechanism, using the
#' O(N) algorithm of Zakov, Kinsella & Bafna (2013) -- in contrast to
#' the earlier, worst-case-exponential BFB-Pivot/BFB-Tree algorithms of
#' Kinsella & Bafna (2012) (an earlier implementation of this package
#' under those algorithms could need a search-step budget on hard
#' instances). This function never needs one: it always returns a
#' definitive TRUE/FALSE, in time proportional to the input size.
#'
#' @param n_vec integer count-vector (telomere to centromere), all
#'   entries >= 1, e.g. \code{c(6, 3, 5)}.
#' @return a single logical: TRUE if the count-vector admits a BFB
#'   schedule, FALSE otherwise.
#' @examples
#' admits_bfb(c(6, 3, 5))                      # TRUE
#' admits_bfb(c(14, 7, 18, 16, 9, 12))          # FALSE
#' admits_bfb(c(14, 7, 19, 17, 9, 13))          # TRUE (paper's Fig. 5 near-miss)
#' @export
admits_bfb <- function(n_vec) {
  n_vec <- .bfbtools_check_counts(n_vec)
  decision_bfb_cpp(n_vec)
}

#' Find the nearest BFB-admitting count-vector (noisy/weighted search)
#'
#' Finds the highest-weight BFB count-vector explaining an observed,
#' potentially noisy count-vector, using the linear-fold-core weighted
#' search (ported from Zakov, Kinsella & Bafna 2013's algorithm family;
#' see \code{src/distance.cpp} for the specific design). Unlike an
#' L1-ball neighbor search bounded by a fixed radius (which can miss the
#' true nearest vector or require tuning the radius -- the approach an
#' earlier, exponential-search-based implementation of this package
#' used), this performs an unbounded, Pareto-frontier-pruned search
#' under an explicit probabilistic error model, and finds the *exact*
#' highest-weight vector within the search's weight cutoff (not an
#' approximation).
#'
#' @param n_vec observed integer count-vector, telomere to centromere.
#' @param model error model: \code{"poisson"} (default; observed counts
#'   modeled as noisy Poisson draws around the true count) or
#'   \code{"none"} (no error tolerated -- equivalent to
#'   \code{\link{admits_bfb}}, provided for consistency/testing).
#' @param min_weight search cutoff: candidates whose cumulative weight
#'   would fall below this value are pruned (default 1e-6). Lower values
#'   search further from the observed vector at the cost of more work.
#' @param max_frontier safety cap on the number of live candidates
#'   tracked at once; returns \code{NA} (with a warning) rather than
#'   growing unboundedly if this is hit (default 200000).
#' @return a list with \code{weight} (in (0, 1], 1 = exact match) and
#'   \code{vector} (the nearest/corrected BFB count-vector) -- or
#'   \code{NULL} (with a warning) if nothing was found within
#'   \code{min_weight}, or \code{NA} if \code{max_frontier} was hit.
#' @examples
#' nearest_bfb(c(14, 7, 18, 16, 9, 12))  # paper's Fig. 5 near-miss example
#' @export
nearest_bfb <- function(n_vec, model = c("poisson", "none"),
                         min_weight = 1e-6, max_frontier = 200000) {
  n_vec <- .bfbtools_check_counts(n_vec)
  model <- match.arg(model)
  model_code <- if (model == "poisson") 1L else 0L

  res <- nearest_bfb_weighted_cpp(n_vec, model_code, log(min_weight), as.integer(max_frontier))

  if (!is.null(res$budget_exceeded)) {
    warning("nearest_bfb: frontier safety cap exceeded; try a larger max_frontier ",
            "or a stricter (larger) min_weight.", call. = FALSE)
    return(NA)
  }
  if (length(res) == 0) {
    warning("nearest_bfb: no BFB-admitting vector found within min_weight = ",
            min_weight, "; try lowering min_weight.", call. = FALSE)
    return(NULL)
  }
  list(weight = res$weight, vector = as.integer(res$vector))
}

#' Distance from a count-vector to the nearest BFB-admitting vector
#'
#' Convenience wrapper around \code{\link{nearest_bfb}} returning a
#' single numeric "distance" (\code{-log(weight)}: 0 for an exact match,
#' increasing as the observed vector requires more/larger corrections to
#' admit a BFB schedule). This is on the same scale as a per-position
#' negative log-likelihood under the chosen error model, not the
#' Canberra distance an earlier, exponential-search-based implementation
#' of this package used -- the two are not numerically comparable,
#' though they are both zero exactly when the input already admits BFB.
#'
#' @inheritParams nearest_bfb
#' @return a single numeric, or \code{NA} if nothing was found (see
#'   \code{\link{nearest_bfb}}).
#' @examples
#' bfb_distance_weighted(c(14, 7, 18, 16, 9, 12))
#' @export
bfb_distance_weighted <- function(n_vec, model = c("poisson", "none"),
                                   min_weight = 1e-6, max_frontier = 200000) {
  res <- nearest_bfb(n_vec, model, min_weight, max_frontier)
  if (is.null(res) || (length(res) == 1 && is.na(res))) return(NA_real_)
  -log(res$weight)
}
