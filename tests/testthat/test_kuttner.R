library(RGAP)

data("gap")
tsList <- as.list(gap[["Netherlands"]][,c("cpih","gdp")])
tsList$infl <- diff(tsList$cpih)
model <- KuttnerModel(tsl = tsList)
parRestr <- initializeRestr(model, type = "hp")


test_that("Kuttner model", {

  expect_s3_class(model, "KuttnerModel")
  expect_equal(is.KuttnerModel(model, return.logical = TRUE), TRUE)
  expect_snapshot_output(x = model, cran = FALSE)

})

test_that("Kuttner MLE fit", {

  skip_on_cran()
  fit <- fitKuttner(model, parRestr = parRestr)

  expect_s3_class(fit, "KuttnerFit")
  expect_equal(RGAP:::is.KuttnerFit(fit, return.logical = TRUE), TRUE)
  expect_snapshot_output(x = fit, cran = FALSE)

})


