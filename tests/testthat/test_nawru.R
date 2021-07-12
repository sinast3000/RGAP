library(RGAP)

data("gap")
tsList <- amecoData2input(gap[["France"]], alpha = 0.65)
exoType <- initializeExo(varNames = "ws")
exoType[1, , "difference"] <- 2
exoType[1, , "lag"] <- 0
model <- NAWRUmodel(tsl = tsList, exoType = exoType)

test_that("NAWRU model", {

  expect_s3_class(model, "NAWRUmodel")
  expect_equal(is.NAWRUmodel(model, return.logical = TRUE), TRUE)
  expect_snapshot_output(x = model, cran = FALSE)

})

test_that("NAWRU MLE fit", {

  skip_on_cran()
  fit <-  fitNAWRU(model)

  expect_s3_class(fit, "NAWRUfit")
  expect_equal(RGAP:::is.NAWRUfit(fit, return.logical = TRUE), TRUE)
  expect_snapshot_output(x = fit, cran = FALSE)

})

test_that("NAWRU bayesian fit", {

  skip_on_cran()
  fit <-  fitNAWRU(model)
  set.seed(5)
  fitBayes <- fitNAWRU(model = model, method = "bayesian", R = 1000, thin = 2, MLEfit = fit)

  expect_s3_class(fitBayes, "NAWRUfit")
  expect_equal(RGAP:::is.NAWRUfit(fitBayes, return.logical = TRUE), TRUE)
  expect_snapshot_output(x = fitBayes, cran = FALSE)

})
