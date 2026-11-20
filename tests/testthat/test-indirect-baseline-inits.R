# Every indirect-response template starts at its drug-free baseline. Solved
# without drug from its own default initial conditions, the response, and in
# the precursor templates the precursor pool, must not move, except in the
# circadian templates, whose response must follow the circadian baseline from
# t = 0. Without drug each is a linear system with one steady state, so a
# solution that does not move can only be that steady state. The templates are
# taken from the registry, so a template added later is covered too.

test_that("every indirect_* template starts at its drug-free baseline", {
  skip_on_cran()
  templates <- grep("^indirect_", nlmixr2lib::modeldb$name, value = TRUE)
  # A template file the registry missed would drop out of the loop unnoticed.
  onDisk <- list.files(system.file("modeldb", package = "nlmixr2lib"), pattern = "^indirect_.*[.]R$", recursive = TRUE)
  expect_setequal(templates, sub("[.]R$", "", basename(onDisk)))
  # Update this count when an indirect-response template is added or removed.
  expect_identical(length(templates), 24L)
  times <- seq(0, 48, by = 0.5)
  for (nm in templates) {
    ui <- rxode2::as.rxUi(nlmixr2lib::readModelDb(nm))
    sol <- rxode2::rxSolve(
      ui,
      events = rxode2::et(times),
      atol = 1e-12,
      rtol = 1e-10,
      returnType = "data.frame",
      addDosing = FALSE
    )
    r0 <- sol$effect[[1]]
    expect_gt(r0, 0, label = paste(nm, "response at t = 0"))
    if (startsWith(nm, "indirect_circ_")) {
      theta <- ui$theta
      rm <- exp(theta[["lrm"]])
      baseline <- rm + exp(theta[["lra"]]) * cos(0.2616 * (sol$time - exp(theta[["ltacro"]])))
      drift <- max(abs(sol$effect - baseline)) / rm
    } else {
      drift <- max(abs(sol$effect - r0)) / r0
    }
    expect_lt(drift, 1e-8, label = paste(nm, "relative drift of the response from its baseline"))
    if ("precursor1" %in% names(sol)) {
      p0 <- sol$precursor1[[1]]
      expect_lt(
        max(abs(sol$precursor1 - p0)) / p0,
        1e-8,
        label = paste(nm, "relative drift of the precursor pool")
      )
    }
  }
})
