library(RGAP)

data("gap")
tsList <- amecoData2input(gap[["Italy"]], alpha = 0.65)
model <- TFPmodel(tsl = tsList, cycle = "RAR2")

test_that("TFP model", {

  expect_s3_class(model, "TFPmodel")
  expect_equal(is.TFPmodel(model, return.logical = TRUE), TRUE)
  expect_snapshot_output(x = model, cran = FALSE)

})

test_that("TFP MLE fit", {

  skip_on_cran()
  fit <- fitTFP(model = model)

  expect_s3_class(fit, "TFPfit")
  expect_equal(RGAP:::is.TFPfit(fit, return.logical = TRUE), TRUE)
  expect_snapshot_output(x = fit, cran = FALSE)

})

test_that("TFP bayesian fit", {

  skip_on_cran()
  fit <- fitTFP(model = model)
  set.seed(5)
  fitBayes <- fitTFP(model = model, method = "bayesian", R = 1000, thin = 2, MLEfit = fit)

  expect_s3_class(fitBayes, "TFPfit")
  expect_equal(RGAP:::is.TFPfit(fitBayes, return.logical = TRUE), TRUE)
  expect_snapshot_output(x = fitBayes, cran = FALSE)

})
