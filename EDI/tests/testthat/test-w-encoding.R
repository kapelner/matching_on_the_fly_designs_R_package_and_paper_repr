expect_w_encoding = function(w, n = NULL) {
  expect_true(all(w %in% c(-1L, 1L)),
    info = paste("w contains values outside {-1, +1}:", paste(unique(w), collapse = ", ")))
  if (!is.null(n)) expect_length(w, n)
}

expect_ws_encoding = function(ws, n = NULL, r = NULL) {
  expect_true(all(ws %in% c(-1L, 1L)),
    info = paste("ws matrix contains values outside {-1, +1}:", paste(unique(as.vector(ws)), collapse = ", ")))
  if (!is.null(n)) expect_equal(nrow(ws), n)
  if (!is.null(r)) expect_equal(ncol(ws), r)
}

# ── helpers ─────────────────────────────────────────────────────────────────

run_seq_design = function(des, X) {
  for (i in seq_len(nrow(X))) {
    des$add_one_subject_to_experiment_and_assign(X[i, , drop = FALSE])
  }
  des
}

run_fixed_design = function(des, X) {
  des$add_all_subjects_to_experiment(X)
  des$assign_w_to_all_subjects()
  des
}

n   = 10L
n12 = 12L
X   = data.frame(x1 = rnorm(n),   x2 = rnorm(n))
X12 = data.frame(x1 = rnorm(n12), x2 = rnorm(n12))
# Xf: n=10 with strata and cluster columns for blocking/cluster designs
Xf = data.frame(
  x1      = rnorm(n),
  strata  = rep(c("A", "B"), each = n / 2),
  cluster = rep(1:5, each = 2)
)
# Xf12: n=12 with clusters nested within blocks (required by BlockedCluster)
Xf12 = data.frame(
  x1      = rnorm(n12),
  strata  = rep(c("A", "B"), each = n12 / 2),
  cluster = c(rep(1:3, each = 2), rep(4:6, each = 2))
)

# ── Sequential designs: get_w() ───────────────────────────────────────────────

test_that("DesignSeqOneByOneBernoulli get_w returns {-1,+1}", {
  set.seed(1)
  des = run_seq_design(DesignSeqOneByOneBernoulli$new(n = n, response_type = "continuous"), X)
  expect_w_encoding(des$get_w(), n)
})

test_that("DesignSeqOneByOneEfron get_w returns {-1,+1}", {
  set.seed(1)
  des = run_seq_design(DesignSeqOneByOneEfron$new(n = n, response_type = "continuous"), X)
  expect_w_encoding(des$get_w(), n)
})

test_that("DesignSeqOneByOneAtkinson get_w returns {-1,+1}", {
  set.seed(1)
  des = run_seq_design(DesignSeqOneByOneAtkinson$new(n = n, response_type = "continuous"), X)
  expect_w_encoding(des$get_w(), n)
})

test_that("DesignSeqOneByOneKK14 get_w returns {-1,+1}", {
  set.seed(1)
  des = run_seq_design(DesignSeqOneByOneKK14$new(n = n, response_type = "continuous"), X)
  expect_w_encoding(des$get_w(), n)
})

test_that("DesignSeqOneByOneSPBR get_w returns {-1,+1}", {
  set.seed(1)
  des = run_seq_design(
    DesignSeqOneByOneSPBR$new(strata_cols = "strata", n = n, response_type = "continuous"),
    Xf
  )
  expect_w_encoding(des$get_w(), n)
})

test_that("DesignSeqOneByOneiBCRD get_w returns {-1,+1}", {
  set.seed(1)
  des = run_seq_design(DesignSeqOneByOneiBCRD$new(n = n, response_type = "continuous"), X)
  expect_w_encoding(des$get_w(), n)
})

test_that("DesignSeqOneByOnePocockSimon get_w returns {-1,+1}", {
  set.seed(1)
  des = run_seq_design(
    DesignSeqOneByOnePocockSimon$new(strata_cols = "strata", n = n, response_type = "continuous"),
    Xf
  )
  expect_w_encoding(des$get_w(), n)
})

test_that("DesignSeqOneByOneUrn get_w returns {-1,+1}", {
  set.seed(1)
  des = run_seq_design(DesignSeqOneByOneUrn$new(n = n, response_type = "continuous"), X)
  expect_w_encoding(des$get_w(), n)
})

# ── Fixed designs: get_w() and draw_ws_according_to_design() ─────────────────

test_that("DesignFixedBernoulli get_w and draw_ws return {-1,+1}", {
  set.seed(1)
  des = run_fixed_design(DesignFixedBernoulli$new(n = n, response_type = "continuous"), X)
  expect_w_encoding(des$get_w(), n)
  expect_ws_encoding(des$draw_ws_according_to_design(5L), n, 5L)
})

test_that("DesignFixedRerandomization get_w and draw_ws return {-1,+1}", {
  set.seed(1)
  des = run_fixed_design(DesignFixedRerandomization$new(n = n, response_type = "continuous"), X)
  expect_w_encoding(des$get_w(), n)
  expect_ws_encoding(des$draw_ws_according_to_design(5L), n, 5L)
})

test_that("DesignFixedBlocking get_w and draw_ws return {-1,+1}", {
  set.seed(1)
  des = run_fixed_design(
    DesignFixedBlocking$new(n = n, response_type = "continuous",
                            strata_cols = "strata", equal_block_sizes = FALSE),
    Xf
  )
  expect_w_encoding(des$get_w(), n)
  expect_ws_encoding(des$draw_ws_according_to_design(5L), n, 5L)
})

test_that("DesignFixedBinaryMatch get_w and draw_ws return {-1,+1}", {
  set.seed(1)
  des = run_fixed_design(DesignFixedBinaryMatch$new(n = n, response_type = "continuous"), X)
  expect_w_encoding(des$get_w(), n)
  expect_ws_encoding(des$draw_ws_according_to_design(5L), n, 5L)
})

test_that("DesignFixedGreedy get_w and draw_ws return {-1,+1}", {
  set.seed(1)
  des = run_fixed_design(DesignFixedGreedy$new(n = n, response_type = "continuous"), X)
  expect_w_encoding(des$get_w(), n)
  expect_ws_encoding(des$draw_ws_according_to_design(5L), n, 5L)
})

test_that("DesignFixedOptimalBlocks get_w and draw_ws return {-1,+1}", {
  set.seed(1)
  des = run_fixed_design(DesignFixedOptimalBlocks$new(n = n12, response_type = "continuous"), X12)
  expect_w_encoding(des$get_w(), n12)
  expect_ws_encoding(des$draw_ws_according_to_design(5L), n12, 5L)
})

test_that("DesignFixedMatchingGreedyPairSwitching get_w and draw_ws return {-1,+1}", {
  set.seed(1)
  des = run_fixed_design(DesignFixedMatchingGreedyPairSwitching$new(n = n12, response_type = "continuous"), X12)
  expect_w_encoding(des$get_w(), n12)
  expect_ws_encoding(des$draw_ws_according_to_design(5L), n12, 5L)
})

test_that("DesignFixediBCRD get_w and draw_ws return {-1,+1}", {
  set.seed(1)
  des = run_fixed_design(DesignFixediBCRD$new(n = n, response_type = "continuous"), X)
  expect_w_encoding(des$get_w(), n)
  expect_ws_encoding(des$draw_ws_according_to_design(5L), n, 5L)
})

test_that("DesignFixedCluster get_w and draw_ws return {-1,+1}", {
  set.seed(1)
  des = run_fixed_design(
    DesignFixedCluster$new(n = n, response_type = "continuous", cluster_col = "cluster"),
    Xf
  )
  expect_w_encoding(des$get_w(), n)
  expect_ws_encoding(des$draw_ws_according_to_design(5L), n, 5L)
})

test_that("DesignFixedBlockedCluster get_w and draw_ws return {-1,+1}", {
  set.seed(1)
  des = run_fixed_design(
    DesignFixedBlockedCluster$new(n = n12, response_type = "continuous",
                                  strata_cols = "strata", cluster_col = "cluster"),
    Xf12
  )
  expect_w_encoding(des$get_w(), n12)
  expect_ws_encoding(des$draw_ws_according_to_design(5L), n12, 5L)
})
