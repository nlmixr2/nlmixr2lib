# An inhibition factor 1 - Imax*C/(C + IC50) stays non-negative at every
# concentration only if Imax <= 1. Above 1 it turns negative at high
# concentration and drives the inhibited production or loss rate negative.
# Every inhibition template's default Imax must therefore lie in (0, 1]. A
# template with no Imax parameter inhibits completely, which is Imax = 1. The
# templates are taken from the registry, so one added later is covered too.

test_that("every indirect inhibition template's default Imax is in (0, 1]", {
  skip_on_cran()
  templates <- grep("^indirect_.*_inhi_", nlmixr2lib::modeldb$name, value = TRUE)
  expect_gt(length(templates), 0L)
  for (nm in templates) {
    ui <- rxode2::as.rxUi(nlmixr2lib::readModelDb(nm))
    theta <- ui$theta
    if ("limax" %in% names(theta)) {
      imax <- exp(theta[["limax"]])
      expect_gt(imax, 0, label = paste(nm, "default Imax"))
      expect_lte(imax, 1, label = paste(nm, "default Imax"))
    } else {
      # Without an Imax parameter the model must not compute one either.
      expect_false("imax" %in% ui$mv0$lhs, label = paste(nm, "computes imax without estimating it"))
    }
  }
})
