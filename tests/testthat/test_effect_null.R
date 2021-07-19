context("effect_null - class of object")

test_that("manyglm object", {

  expect_error(
    effect_null(fit1.cord, term_cont)
  )

})

context("effect_null matrix")

test_that("term_cont", {

  for (fit in list(fit1.glm, fit2.glm)) {
    returned = effect_null(fit, term_cont)

    expect_equal(rownames(returned), rownames(fit$coefficients))
    expect_equal(colnames(returned), colnames(fit$coefficients))
    expect_equal(sum(returned[term_cont,]), 0)
    expect_true(sum(returned["(Intercept)",]) != 0)
  }

})

test_that("term_factors", {

  for (fit in list(fit_fac2.glm, fit_fac4.glm)) {

    term_fac = labels(terms(fit))
    returned = effect_null(fit, term_fac)

    expect_equal(rownames(returned), rownames(fit$coefficients))
    expect_equal(colnames(returned), colnames(fit$coefficients))
    expect_equal(sum(returned[2:nrow(returned),]), 0)
    expect_true(sum(returned["(Intercept)",]) != 0)
  }

})

test_that("coeffs == coeffs0, one covariate", {

  for (fit in list(fit1.glm, fit_fac2.glm)) {
    term = labels(terms(fit))
    returned = effect_null(fit, term)

    whichRow = which(startsWith(rownames(returned), term))
    returned0 = t(matrix(returned[-whichRow,]))
    colnames(returned0) = colnames(returned)
    expect_equal(returned0, fit0.glm$coefficients)
  }

})

test_that("coeffs == coeffs0, two covariates", {

  returned = effect_null(fit2.glm, term_cont)

  whichRow = match(term_cont, rownames(returned))
  returned = returned[-whichRow,]
  expect_equal(returned, fit2_null.glm$coefficients)

})

context("effect_null - term")

test_that("interaction term", {

  expect_error(
    effect_null(fit_int.glm, "bare.sand:Treatment2B")
  )

})

test_that("type of term", {

  expect_error(
    effect_null(fit1.glm, NULL)
  )

})

context("effect_null - object")

test_that("model without intercept", {

  expect_error(
    effect_null(fit00.glm, term_cont)
  )

})
