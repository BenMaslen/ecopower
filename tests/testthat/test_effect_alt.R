context("effect_alt - class of object")

test_that("manyglm object", {

  expect_error(
    effect_alt(fit1.cord, effect_size, increasers, decreasers, term_cont)
  )

})

context("effect_alt matrix")

test_that("term_cont", {

  for (fit in list(fit1.glm, fit2.glm)) {
    returned = effect_alt(fit, effect_size, increasers, decreasers, term_cont)

    expect_equal(rownames(returned), rownames(fit$coefficients))
    expect_equal(colnames(returned), colnames(fit$coefficients))
  }

})

test_that("term_factors", {

  for (fit in list(fit_fac2.glm, fit_fac4.glm)) {

    term_fac = labels(terms(fit))
    returned = effect_alt(fit, effect_size, increasers, decreasers, term_fac)

    expect_equal(rownames(returned), rownames(fit$coefficients))
    expect_equal(colnames(returned), colnames(fit$coefficients))
  }

})

context("effect_alt - effect size")

test_that("increasers, decreasers and no_effect - term_cont", {

  for (fit in list(fit1.glm, fit2.glm)) {
    returned = effect_alt(fit, effect_size, increasers, decreasers, term_cont)

    for (j in 1:length(increasers)) {
      expect_equal(unname(returned[term_cont,][increasers][j]), log(1.5))
    }

    for (j in 1:length(decreasers)) {
      expect_equal(unname(returned[term_cont,][decreasers][j]), -log(1.5))
    }

    i = nrow(returned)
    expect_equal(sum(returned[i,no_effect]), 0)
  }

})

test_that("increasers, decreasers and no_effect - nlevels = 2", {

  term = "Treatment2"
  for (fit in list(fit_fac2.glm, fit_mth.glm)) {
    returned = effect_alt(fit, effect_size, increasers, decreasers, term)
    i = nrow(returned)

    for (j in 1:length(increasers)) {
      expect_equal(unname(returned[,increasers][i,j]), log(1.5^1))
    }

    for (j in 1:length(decreasers)) {
      expect_equal(unname(returned[,decreasers][i,j]), -log(1.5^1))
    }

    expect_equal(sum(returned[i,no_effect]), 0)
  }

})

test_that("increasers, decreasers and no_effect - nlevels = 4", {

  term = "Treatment4"
  returned = effect_alt(fit_fac4.glm, effect_size, increasers, decreasers, term)
  rowNames = row.names(returned)[startsWith(row.names(returned), term)]

  for (i in 1:length(rowNames)) {
    for (j in 1:length(increasers)) {
      expect_equal(unname(returned[rowNames[i],][increasers][j]), log(1.5^i))
    }
  }

  for (i in 1:length(rowNames)) {
    for (j in 1:length(decreasers)) {
      expect_equal(unname(returned[rowNames[i],][decreasers][j]), -log(1.5^i))
    }
  }

  expect_equal(sum(returned[2:4,no_effect]), 0)

})

context("effect_alt - K")

test_that("K = term levels - 1", {

  for (K in c(2,5)) {

    expect_error(
      effect_alt(fit_fac4.glm, effect_size, increasers, decreasers, "Treatment4", K)
    )
  }

})

test_that("input new vector K", {

  K = c(2,3,1)
  term = "Treatment4"
  returned = effect_alt(fit_fac4.glm, effect_size, increasers, decreasers, term, K)
  rowNames = row.names(returned)[startsWith(row.names(returned), term)]

  for (i in 1:length(rowNames)) {
    for (j in 1:length(increasers)) {
      expect_equal(unname(returned[rowNames[i],][increasers][j]), log(1.5^K[i]))
    }
  }

  for (i in 1:length(rowNames)) {
    for (j in 1:length(decreasers)) {
      expect_equal(unname(returned[rowNames[i],][decreasers][j]), -log(1.5^K[i]))
    }
  }

  expect_equal(sum(returned[2:4,no_effect]), 0)

})
