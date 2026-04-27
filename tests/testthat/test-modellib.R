test_that("modellib", {
  expect_output(
    modellib(),
    regexp = "PK_1cmt"
  )
  modelTest <- suppressMessages(modellib(name = "PK_1cmt", eta = c("ka", "vc"), reserr = "addSd"))

  # Expected ini and model blocks (no metadata preamble, so robust to future
  # additions of description/reference/units above ini())
  expected <- function() {
    ini({
      lka <- 0.45
      label("Absorption rate (Ka)")
      lcl <- 1
      label("Clearance (CL)")
      lvc <- 3.45
      label("Central volume of distribution (V)")
      CcAddSd <- c(0, 1)
      etaKa ~ 0.1
      etaVc ~ 0.1
    })
    model({
      ka <- exp(lka + etaKa)
      cl <- exp(lcl)
      vc <- exp(lvc + etaVc)
      Cc <- linCmt()
      Cc ~ add(CcAddSd)
    })
  }

  expect_equal(findBlock(as.function(modelTest), "ini"),   findBlock(expected, "ini"))
  expect_equal(findBlock(as.function(modelTest), "model"), findBlock(expected, "model"))
})
