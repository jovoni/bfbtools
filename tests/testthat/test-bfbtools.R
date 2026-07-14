test_that("admits_bfb reproduces the paper's own worked examples", {
  expect_true(admits_bfb(c(6, 3, 5)))
  expect_false(admits_bfb(c(14, 7, 18, 16, 9, 12)))
  expect_true(admits_bfb(c(14, 7, 19, 17, 9, 13)))
  expect_false(admits_bfb(c(7, 13, 6, 4, 19, 12)))
  expect_true(admits_bfb(c(8, 14, 6, 4, 20, 12)))
})

test_that("admits_bfb respects known necessary conditions (Kinsella & Bafna 2012)", {
  # Rule of One (Lemma 11): if n_i = 1 and n_j > 1 for some i < j, the
  # vector cannot admit a BFB schedule.
  expect_false(admits_bfb(c(1, 3)))
  expect_false(admits_bfb(c(5, 1, 2)))
  # Odd-Even Rule (Lemma 12): no odd count may precede an even count.
  expect_false(admits_bfb(c(3, 2)))
  expect_false(admits_bfb(c(5, 3, 4)))
})

test_that("admits_bfb is stable under the Lemma 7 doubling construction", {
  # [n1,...,nk] admits BFB iff [2n1,...,2nk,2] does
  v <- c(6, 3, 5)
  expect_true(admits_bfb(v))
  for (i in 1:5) {
    v <- c(2 * v, 2)
    expect_true(admits_bfb(v))
  }
})

test_that("admits_bfb handles large inputs in well under a second", {
  set.seed(3)
  for (k in c(100, 1000)) {
    nv <- sample(1:1000, k, replace = TRUE)
    t0 <- Sys.time()
    admits_bfb(nv)
    expect_lt(as.numeric(Sys.time() - t0), 1)
  }
})

test_that("nearest_bfb with model='none' reconstructs admits_bfb exactly", {
  expect_equal(nearest_bfb(c(6, 3, 5), model = "none")$vector, c(6, 3, 5))
  expect_equal(nearest_bfb(c(6, 3, 5), model = "none")$weight, 1)
  expect_null(suppressWarnings(nearest_bfb(c(14, 7, 18, 16, 9, 12), model = "none")))
})

test_that("nearest_bfb (Poisson) reproduces the paper's Figure 5 nearest vectors", {
  check_one <- function(v, expected_nn) {
    r <- nearest_bfb(v)
    expect_equal(r$vector, expected_nn)
  }
  check_one(c(14,7,18,16,9,12), c(14,7,19,17,9,13))
  check_one(c(7,13,6,4,19,12),  c(8,14,6,4,20,12))
  check_one(c(9,7,7,8,2,14),    c(8,8,8,8,2,14))
  check_one(c(16,17,14,18,1,19),c(16,18,14,18,2,19))
  check_one(c(20,1,7,3,10,6),   c(20,2,7,3,11,7))
})

test_that("nearest_bfb always returns a vector that genuinely admits BFB", {
  set.seed(1)
  for (i in 1:300) {
    k <- sample(2:8, 1)
    nv <- sample(1:20, k, replace = TRUE)
    r <- suppressWarnings(nearest_bfb(nv))
    if (is.null(r) || (length(r) == 1 && is.na(r))) next
    expect_true(admits_bfb(r$vector), info = paste(nv, collapse = ","))
  }
})

test_that("nearest_bfb returns the input itself (weight 1) when it already admits BFB", {
  r <- nearest_bfb(c(6, 3, 5))
  expect_equal(r$vector, c(6, 3, 5))
  expect_equal(r$weight, 1)
})

test_that("bfb_distance_weighted is 0 for admitting vectors and positive otherwise", {
  expect_equal(bfb_distance_weighted(c(6, 3, 5)), 0)
  expect_gt(bfb_distance_weighted(c(14, 7, 18, 16, 9, 12)), 0)
})

test_that("input validation rejects bad input", {
  expect_error(admits_bfb(c(0, 1, 2)))
  expect_error(admits_bfb(c(1, -1)))
  expect_error(admits_bfb(character(0)))
})
