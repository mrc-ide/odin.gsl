context("example")

test_that("confirm that dde example confirms to package", {
  skip_if_not_installed("dde")

  ## Here's the reference model
  path <- system.file("examples/sir.R", package = "odin.gsl", mustWork = TRUE)
  gen <- odin::odin_(path)
  mod <- gen()

  y <- mod$initial(0)
  steps <- seq(0, 100, by = 4)
  replicate <- 20

  ## Run with odin/dde
  seed <- 1
  set.seed(seed)
  cmp <- mod$run(steps, y, replicate = replicate)

  res <- dde_example(replicate, steps, y, threads = 4, user=list(seed=seed))
  out <- attr(res, "output")
  attr(res, "output") <- NULL

  ## First iteration is easy to visualise and reasona about:
  expect_equal(res[, , 1], t(unname(cmp[, 2:4, 1])))
  expect_equal(out[, , 1], t(unname(cmp[, 5:6, 1])))

  ## 10th only happens if everything is legit
  expect_equal(res[, , 10], t(unname(cmp[, 2:4, 10])))
  expect_equal(out[, , 10], t(unname(cmp[, 5:6, 10])))

  ## Whole object at once:
  expect_equal(res, unname(aperm(cmp[, 2:4, ], c(2, 1, 3))))
  expect_equal(out, unname(aperm(cmp[, 5:6, ], c(2, 1, 3))))
})
