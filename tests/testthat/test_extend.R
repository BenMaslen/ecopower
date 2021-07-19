context("extend - class of object")

test_that("cord object", {

  expect_error(
    extend(fit1.glm)
  )

})

context("compute_N")

test_that("n_replicate overwrites N", {

  expect_warning(
    extend(fit1.cord, n_replicate=2)
  )

})

context("check extended data")

test_that("N simulation", {

  for (fit in list(fit1.cord, fit2.cord, fit_fac4.cord, fit_mth.cord, fit_mix.cord)) {

    returned = extend(fit, N=nrow(spiddat)*2)$data
    expect_equal(nrow(returned), 56)
    expect_equal(
      colnames(returned)[1:ncol(spiddat)], colnames(spiddat)
    )
    expect_equal(
      colnames(returned)[(ncol(spiddat)+1):(ncol(spiddat)+ncol(X))],
      colnames(X)
    )
  }

})

test_that("n_replicate simulation", {

  for (fit in list(fit1.cord, fit2.cord, fit_fac4.cord, fit_mth.cord, fit_mix.cord)) {

    returned = extend(fit, N=NULL, n_replicate=3)$data
    expect_equal(nrow(returned), 84)
    expect_equal(
      colnames(returned)[1:ncol(spiddat)], colnames(spiddat)
    )
    expect_equal(
      colnames(returned)[(ncol(spiddat)+1):(ncol(spiddat)+ncol(X))],
      colnames(X)
    )
  }

})

context("intercept model")

test_that("do.fit=FALSE, no data specified", {

  returned = extend(fit0.cord, N=nrow(spiddat)*2)$data
  returned1 = extend(fit0.cord, N=NULL, n_replicate=2)$data
  expect_equal(nrow(returned), 56)
  expect_equal(nrow(returned), nrow(returned1))
  expect_equal(colnames(returned), colnames(spiddat))
  expect_equal(colnames(returned), colnames(returned1))

})

test_that("do.fit=FALSE, data specified", {

  returned = extend(fit0_X.cord, N=nrow(spiddat)*2)$data
  returned1 = extend(fit0_X.cord, N=NULL, n_replicate=2)$data
  expect_equal(nrow(returned), 56)
  expect_equal(nrow(returned), nrow(returned1))
  expect_equal(colnames(returned), colnames(spiddat))
  expect_equal(colnames(returned), colnames(returned1))

})

test_that("do.fit=TRUE, no data specified", {

  returned = extend(fit0.cord, N=nrow(spiddat)*2, do.fit=TRUE)$data
  returned1 = extend(fit0.cord, N=NULL, n_replicate=2, do.fit=TRUE)$data
  expect_equal(nrow(returned), 56)
  expect_equal(nrow(returned), nrow(returned1))
  expect_equal(colnames(returned), colnames(spiddat))
  expect_equal(colnames(returned), colnames(returned1))

})

test_that("do.fit=TRUE, data specified", {

  returned = extend(fit0_X.cord, N=nrow(spiddat)*2, do.fit=TRUE)$data
  returned1 = extend(fit0_X.cord, N=NULL, n_replicate=2, do.fit=TRUE)$data
  expect_equal(nrow(returned), 28)
  expect_equal(nrow(returned), nrow(returned1))
  expect_equal(colnames(returned), colnames(X))
  expect_equal(colnames(returned), colnames(returned1))

})

context("newdata - new factor")

test_that("N simulation - new factor", {

  for (fit in list(fit_fac2.cord, fit_mth.cord)) {

    returned = extend(fit, N=30, newdata=Xnew)$data
    expect_equal(nrow(returned), 30)
    expect_equal(
      colnames(returned)[1:ncol(spiddat)], colnames(spiddat)
    )
    expect_equal(
      colnames(returned)[(ncol(spiddat)+1):(ncol(spiddat)+ncol(X))],
      colnames(X)
    )
  }

})

test_that("n_replicate simulation - new factor", {

  for (fit in list(fit_fac2.cord, fit_mth.cord)) {

    returned = extend(fit, N=NULL, newdata=Xnew, n_replicate=2)$data
    expect_equal(nrow(returned), 56)
    expect_equal(
      colnames(returned)[1:ncol(spiddat)], colnames(spiddat)
    )
    expect_equal(
      colnames(returned)[(ncol(spiddat)+1):(ncol(spiddat)+ncol(X))],
      colnames(X)
    )
  }

})

context("newdata - new vector")

test_that("N simulation - new vector, imexplicit N", {

  returned = extend(fit_vec.cord, newdata=Xvec2)$data
  expect_equal(nrow(returned), 28)
    expect_equal(
    colnames(returned)[1:ncol(spiddat)], colnames(spiddat)
  )
  expect_equal(
    colnames(returned)[(ncol(spiddat)+1):(ncol(spiddat)+ncol(Xvec))],
    colnames(Xvec)
  )

})

test_that("N simulation and n_replicate - new vector", {

  returned = extend(fit_vec.cord, N=30, newdata=Xvec2)$data
  expect_equal(nrow(returned), 30)
  expect_false(any(is.null(returned)))
  expect_equal(
    colnames(returned)[1:ncol(spiddat)], colnames(spiddat)
  )
  expect_equal(
    colnames(returned)[(ncol(spiddat)+1):(ncol(spiddat)+ncol(Xvec))],
    colnames(Xvec)
  )

})

test_that("n_replicate simulation - new vector", {

  returned = extend(fit_vec.cord, N=NULL, n_replicate=15, newdata=Xvec2)$data
  expect_equal(nrow(returned), 30)
  expect_false(any(is.null(returned)))
  expect_equal(
    colnames(returned)[1:ncol(spiddat)], colnames(spiddat)
  )
  expect_equal(
    colnames(returned)[(ncol(spiddat)+1):(ncol(spiddat)+ncol(Xvec))],
    colnames(Xvec)
  )

})

context("simulate smaller size")

test_that("N simulation", {

  returned = extend(fit1.cord, N=2)$data
  expect_equal(nrow(returned), 2)
  expect_false(any(is.null(returned)))
  expect_equal(
    colnames(returned)[1:ncol(spiddat)], colnames(spiddat)
  )
  expect_equal(
    colnames(returned)[(ncol(spiddat)+1):(ncol(spiddat)+ncol(X))],
    colnames(X)
  )

})

test_that("n_replicate simulation", {

  returned = extend(fit_vec.cord, N=NULL, n_replicate=1)$data
  expect_equal(nrow(returned), 2)
  expect_false(any(is.null(returned)))
  expect_equal(
    colnames(returned)[1:ncol(spiddat)], colnames(spiddat)
  )
  expect_equal(
    colnames(returned)[(ncol(spiddat)+1):(ncol(spiddat)+ncol(Xvec))],
    colnames(Xvec)
  )

})

test_that("N simulation - newdata", {

  returned = extend(fit_mth.cord, N=2, newdata=Xnew)$data
  expect_equal(nrow(returned), 2)
  expect_false(any(is.null(returned)))
  expect_equal(
    colnames(returned)[1:ncol(spiddat)], colnames(spiddat)
  )
  expect_equal(
    colnames(returned)[(ncol(spiddat)+1):(ncol(spiddat)+ncol(X))],
    colnames(X)
  )

})

test_that("n_replicate simulation - newdata", {

  returned = extend(fit_vec.cord, N=NULL, n_replicate=4, newdata=Xvec2)$data
  expect_equal(nrow(returned), 8)
  expect_false(any(is.null(returned)))
  expect_equal(
    colnames(returned)[1:ncol(spiddat)], colnames(spiddat)
  )
  expect_equal(
    colnames(returned)[(ncol(spiddat)+1):(ncol(spiddat)+ncol(Xvec))],
    colnames(Xvec)
  )

})

test_that("intercept model - with or without data specified", {

  for (fit in list(fit0.cord, fit0_X.cord)) {
    returned = extend(fit, N=14)$data
    returned1 = extend(fit, N=NULL, n_replicate=0.5)$data
    expect_equal(nrow(returned), 14)
    expect_equal(nrow(returned), nrow(returned1))
    expect_equal(colnames(returned), colnames(spiddat))
    expect_equal(colnames(returned), colnames(returned1))
  }

})

context("get_new_fit")

test_that("do.fit=TRUE returns manyglm object", {

  returned = extend(fit1.cord, do.fit=TRUE)
  expect_is(returned, "manyglm")

})

test_that("class of returned object data", {

  returned = extend(fit1.cord, do.fit=TRUE)
  expect_is(returned$data, "data.frame")

})

test_that("do.fit=FALSE returns a data.frame", {

  returned = extend(fit1.cord)$data
  expect_is(returned, "data.frame")

})
