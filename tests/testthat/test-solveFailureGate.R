# A failed subject, as rxode2 returns it: every state NA at every time, while
# the parameter column keeps its value.
gateSim <- function(failedId = 2L) {
  d <- data.frame(
    id = rep(1:3, each = 3),
    time = rep(c(0, 1, 2), 3),
    Cc = c(1, 0.5, 0.25, 2, 1, 0.5, 3, 1.5, 0.75),
    central = c(10, 5, 2.5, 20, 10, 5, 30, 15, 7.5),
    cl = rep(c(1, 2, 3), each = 3)
  )
  d[d$id %in% failedId, c("Cc", "central")] <- NA_real_
  d
}

test_that(".solveFailedSubjects() names the subjects whose states are all NA", {
  expect_identical(.solveFailedSubjects(gateSim(2L), states = "central"), "id=2")
  expect_identical(
    .solveFailedSubjects(gateSim(c(1L, 3L)), states = "central"),
    c("id=1", "id=3")
  )
  expect_identical(
    .solveFailedSubjects(gateSim(integer(0)), states = "central"),
    character(0)
  )
})

test_that(".solveFailedSubjects() does not count a subject with some finite rows", {
  d <- gateSim(integer(0))
  d$central[d$id == 2L & d$time == 2] <- NA_real_
  expect_identical(.solveFailedSubjects(d, states = "central"), character(0))
})

test_that(".solveFailedSubjects() needs every state column NA", {
  d <- gateSim(integer(0))
  d$central[d$id == 2L] <- NA_real_
  d$depot <- 1
  expect_identical(
    .solveFailedSubjects(d, states = c("depot", "central")),
    character(0)
  )
  d$depot[d$id == 2L] <- NA_real_
  expect_identical(.solveFailedSubjects(d, states = c("depot", "central")), "id=2")
})

test_that(".solveFailedSubjects() falls back to calculated variables, then gives up", {
  d <- gateSim(2L)
  d$central <- NULL
  expect_identical(
    .solveFailedSubjects(d, states = "central", lhs = "Cc"),
    "id=2"
  )
  expect_identical(
    .solveFailedSubjects(d, states = "central"),
    character(0)
  )
  # A model without states cannot fail a solve, whatever its outputs hold.
  expect_identical(
    .solveFailedSubjects(d, states = character(0), lhs = "Cc"),
    character(0)
  )
})

test_that(".solveFailedSubjects() keys on sim.id and id together", {
  d <- rbind(gateSim(integer(0)), gateSim(3L))
  d$sim.id <- rep(1:2, each = 9)
  expect_identical(
    .solveFailedSubjects(d, states = "central"),
    "sim.id=2,id=3"
  )
  # A single-subject solve has neither column.
  one <- gateSim(1L)[1:3, c("time", "central")]
  expect_identical(.solveFailedSubjects(one, states = "central"), "subject 1")
})

test_that(".solveFailedSubjects() accepts a matrix and ignores empty or non-tabular values", {
  m <- as.matrix(gateSim(2L))
  expect_identical(.solveFailedSubjects(m, states = "central"), "id=2")
  expect_identical(
    .solveFailedSubjects(gateSim(2L)[0, ], states = "central"),
    character(0)
  )
  expect_identical(.solveFailedSubjects(NULL, states = "central"), character(0))
  expect_identical(.solveFailedSubjects(1:3, states = "central"), character(0))
})

test_that(".solveFailureGateDecision() fires on a failed subject, not on a clean solve", {
  fired <- .solveFailureGateDecision(gateSim(2L), states = "central", expectFailed = FALSE)
  expect_identical(fired, list(failed = "id=2", nSubjects = 3L, fire = TRUE))
  clean <- .solveFailureGateDecision(gateSim(integer(0)), states = "central", expectFailed = FALSE)
  expect_identical(clean, list(failed = character(0), nSubjects = 0L, fire = FALSE))
})

test_that(".solveFailureGateDecision() does not fire with the opt-in", {
  optedIn <- .solveFailureGateDecision(gateSim(2L), states = "central", expectFailed = TRUE)
  expect_identical(optedIn, list(failed = "id=2", nSubjects = 3L, fire = FALSE))
  withr::with_options(list(nlmixr2lib.expectFailedSubjects = TRUE), {
    expect_false(.solveFailureGateDecision(gateSim(2L), states = "central")$fire)
  })
  withr::with_options(list(nlmixr2lib.expectFailedSubjects = NULL), {
    expect_true(.solveFailureGateDecision(gateSim(2L), states = "central")$fire)
  })
})

test_that(".solveFailureGateMessage() names the article, the count and the subjects", {
  d <- .solveFailureGateDecision(gateSim(2L), states = "central", expectFailed = FALSE)
  expect_identical(
    .solveFailureGateMessage("articles/Foo_2020_bar", d),
    paste0(
      "articles/Foo_2020_bar: rxSolve() could not solve 1 of 3 subject(s) ",
      "(id=2); rxode2 replaced their output with NA. Fix the model or the ",
      "article. If the failures are expected, count and assert them in the ",
      "article and set options(nlmixr2lib.expectFailedSubjects = TRUE) around ",
      "that solve."
    )
  )
  many <- list(failed = paste0("id=", 1:12), nSubjects = 20L, fire = TRUE)
  expect_match(
    .solveFailureGateMessage("a", many),
    "could not solve 12 of 20 subject(s) (id=1 id=2 id=3 id=4 id=5 id=6 id=7 id=8 id=9 id=10 and 2 more)",
    fixed = TRUE
  )
})

test_that("the installed gate stops a solve with a failed subject", {
  skip_on_cran()
  # x' = a x^2 from x(0) = 1 blows up at t = 1/a, so the subject with a = 2
  # cannot be solved past t = 0.5 and the others can.
  mod <- rxode2::rxode2({
    x(0) <- 1
    d/dt(x) <- a * x^2
    y <- x + 1
  })
  ev <- rxode2::et(seq(0, 2, by = 0.25)) |> rxode2::et(id = 1:3)
  bad <- data.frame(id = 1:3, a = c(0, 2, -1))
  good <- data.frame(id = 1:3, a = c(0, 0.1, -1))

  .installSolveFailureGate("gate-test")
  withr::defer(.removeSolveFailureGate())

  expect_error(
    suppressWarnings(rxode2::rxSolve(mod, bad, ev)),
    class = "nlmixr2libSolveFailure",
    regexp = "gate-test: rxSolve() could not solve 1 of 3 subject(s) (id=2)",
    fixed = TRUE
  )
  # The error is recorded, so a render that caught it still fails.
  expect_length(.solveFailureGateFired(), 1L)

  # Every return type is judged.
  expect_error(
    suppressWarnings(rxode2::rxSolve(mod, bad, ev, returnType = "data.frame")),
    class = "nlmixr2libSolveFailure"
  )

  sim <- rxode2::rxSolve(mod, good, ev)
  expect_false(anyNA(sim$x))

  withr::with_options(list(nlmixr2lib.expectFailedSubjects = TRUE), {
    sim <- suppressWarnings(rxode2::rxSolve(mod, bad, ev))
  })
  expect_identical(unique(as.character(sim$id[is.na(sim$x)])), "2")
  # The opted-in solve is recorded as accepted, not fired.
  expect_length(.solveFailureGateFired(), 2L)
  expect_length(.solveFailureGateState$records, 3L)
})
