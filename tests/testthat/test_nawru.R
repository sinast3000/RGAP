library(RGAP)

data("gap")
tsList <- amecoData2input(gap[["France"]], alpha = 0.65)
exoType <- initializeExo(varNames = "ws", D = 2, L = 0)
model <- NAWRUmodel(tsl = tsList, exoType = exoType)

test_that("NAWRU model", {

  expect_s3_class(model, "NAWRUmodel")
  expect_equal(is.NAWRUmodel(model, return.logical = TRUE), TRUE)
  expect_snapshot_output(x = model, cran = FALSE)

})

test_that("NAWRU MLE fit", {

  skip_on_cran()
  f <-  fit(model)

  expect_s3_class(f, "NAWRUfit")
  expect_equal(RGAP:::is.NAWRUfit(f, return.logical = TRUE), TRUE)
  expect_snapshot_output(x = f, cran = FALSE)

})

test_that("NAWRU bayesian fit", {

  skip_on_cran()
  f <-  fit(model)
  set.seed(5)
  fBayes <- fit(model = model, method = "bayesian", R = 1000, thin = 2, MLEfit = f)

  expect_s3_class(fBayes, "NAWRUfit")
  expect_equal(RGAP:::is.NAWRUfit(fBayes, return.logical = TRUE), TRUE)
  expect_snapshot_output(x = fBayes, cran = FALSE)

})
