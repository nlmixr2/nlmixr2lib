# Every compartment an article doses must be in its model's `dosing`.
#
# The articles are the only place most models are dosed, so a dose target
# they use but `dosing` leaves out is invisible to the solve probe and the
# dead-state check, and to anyone reading the registry. helper-articleDoses.R
# reads the dose records from each article's code.

# Doses that set or displace an endogenous state, to check a model's
# feedback or attractor, rather than give a drug, tracer or challenge. These
# are not dose targets, so they stay out of `dosing`.
articleStateDoses <- rbind(
  data.frame(
    model = "Decrane_2023_oxyfluorfen_rat",
    state = "t4_serum",
    reason = "a signed bolus that displaces serum T4 from its baseline"
  ),
  data.frame(
    model = "Ousey_2026_plozasiran",
    state = c("apoc3", "tg"),
    reason = "one baseline's worth added, to check the biomarker returns to baseline"
  ),
  data.frame(
    model = "Snelder_2013_cardiovascular_rat",
    state = "co",
    reason = "a signed bolus that displaces cardiac output from its baseline"
  ),
  data.frame(
    model = "Snelder_2014_cardiovascular_rat",
    state = c("hr", "svt", "tpr"),
    reason = "a bolus that displaces heart rate, stroke volume or total peripheral resistance from its baseline"
  ),
  data.frame(
    model = "Tetschke_2018_erythropoiesis",
    state = "thb",
    reason = "a negative bolus: the hemoglobin a blood donation removes"
  ),
  data.frame(
    model = "Thorsted_2019_endotoxin",
    state = c("il6", "tnf"),
    reason = "a bolus that doubles the cytokine at baseline"
  ),
  data.frame(
    model = "Wang_2023_zoledronicAcid_mbma",
    state = "BMD_LS",
    reason = "a signed bolus that displaces bone mineral density from its baseline"
  ),
  data.frame(
    model = "Yin_2024_soticlestat",
    state = "hc24",
    reason = "a signed bolus that displaces 24S-hydroxycholesterol from its baseline"
  ),
  data.frame(
    model = "Zhang_2025_epoetinAlfa_erythroferrone_cia_rat",
    state = "erythrocytes",
    reason = "a transfusion or phlebotomy analogue that displaces the red-cell pool"
  )
)

# Doses the reader credits to a registry model that the article solves with
# another model. Each record is named so the entry can be checked by hand.
articleMiscreditedDoses <- data.frame(
  model = "Aoki_2024_intratarget_microdosing_pbpk",
  state = "target",
  reason = paste(
    "`ev_auth` in the xcheck chunk doses state 9 (`cmt = 9L`, the free receptor",
    "X_FreeR) of the inline `authors` model, inside a loop that also solves the",
    "packaged model; the packaged model sets `target(0) <- rtot` instead"
  )
)

# Models whose `vignette` names an article that does not exist. No article
# solves either of them.
articleMissingArticles <- data.frame(
  model = c("Brekkan_2018_pegfilgrastim", "Hansson_2013_sunitinib_os"),
  reason = c(
    "the Ibrahim 2019 article lists it, nothing doses it",
    "the Hansson 2013 myelosuppression and Centanni 2025 articles name it, nothing doses it"
  )
)

# Dose records whose dose rows the reader cannot tell from their observation
# rows, so an observation's compartment reads as dosed. Each article doses
# only model states.
articleUnpairedRecords <- data.frame(
  article = c(rep("Knibbe_2009_morphine", 3), "Wu_2025_paracetamol"),
  cmt = c("Cc", "Cc_m3g", "Cc_m6g", "Cc"),
  reason = c(
    rep("`cmt = dplyr::case_when(MDV == 1 ~ \"central\", CMT == 1 ~ \"Cc\", ...)`: doses are the MDV == 1 rows", 3),
    "`cmt = ifelse(evid == 1L, \"central\", \"Cc\")`: doses are the evid == 1 rows"
  )
)

# Writes R code into a one-chunk article and reads its doses.
.dosesOfCode <- function(code, models = "M") {
  rmd <- tempfile(fileext = ".Rmd")
  on.exit(unlink(rmd))
  writeLines(c("```{r}", code, "```"), rmd)
  d <- articleDoseCmts(rmd, models)
  d[order(d$model, d$cmt), c("model", "cmt", "owner")]
}

.doseCmtsOfCode <- function(code) {
  sort(unique(.dosesOfCode(code)$cmt))
}

test_that("a dose record is read from its evid and amt", {
  expect_identical(.doseCmtsOfCode('ev <- rxode2::et(amt = 100, cmt = "depot")'), "depot")
  expect_identical(.doseCmtsOfCode('ev <- rxode2::et(dose = 100, cmt = "central")'), "central")
  # Observations, zero or NA amounts and evid 5 resets are not doses.
  expect_identical(.doseCmtsOfCode('ev <- rxode2::et(time = 1:3, cmt = "central")'), character(0))
  expect_identical(.doseCmtsOfCode('ev <- rxode2::et(amt = 0, evid = 1, cmt = "depot2")'), character(0))
  expect_identical(.doseCmtsOfCode('ev <- data.frame(time = 1, amt = NA_real_, cmt = "Cc")'), character(0))
  expect_identical(.doseCmtsOfCode('ev <- data.frame(time = 1, evid = 5L, amt = 0, cmt = "cyp")'), character(0))
  expect_identical(.doseCmtsOfCode('ev <- data.frame(time = 1, evid = 4L, amt = 10, cmt = "depot")'), "depot")
  expect_identical(.doseCmtsOfCode("ev <- rxode2::et(amt = 10, cmt = 2L)"), "#2")
})

test_that("only the dose rows of a vector record count", {
  expect_identical(
    .doseCmtsOfCode(
      c(
        "ev <- data.frame(time = c(0, t), evid = c(1L, rep(0L, n)),",
        '  amt = c(d, rep(NA_real_, n)), cmt = c("depot", rep("central", n)))'
      )
    ),
    "depot"
  )
  expect_identical(
    .doseCmtsOfCode(
      c(
        "ev <- mutate(g, amt = ifelse(time == 0, 120, NA_real_), evid = ifelse(time == 0, 1L, 0L),",
        '  cmt = ifelse(time == 0, "depot", "central"))'
      )
    ),
    "depot"
  )
  expect_identical(
    .doseCmtsOfCode(
      c(
        "ev <- data.frame(amt = c(160, 160, 160, rep(NA_real_, n)), evid = c(1L, 1L, 1L, rep(0L, n)),",
        '  cmt = c(rep("depot", 3), rep("central", n)))'
      )
    ),
    "depot"
  )
  expect_identical(
    .doseCmtsOfCode(
      c(
        "ev <- data.frame(amt = rep(c(a, rep(NA_real_, n)), k), evid = rep(c(1L, rep(0L, n)), k),",
        '  cmt = rep(c("central", rep("Cc", n)), k))'
      )
    ),
    "central"
  )
  # Dose rows first, observation rows split differently between columns.
  expect_identical(
    .doseCmtsOfCode(
      c(
        "ev <- data.frame(amt = c(5, rep(NA_real_, 2 * n)), evid = c(1L, rep(0L, 2 * n)),",
        '  cmt = c("depot", rep(a, n), rep(b, n)))'
      )
    ),
    "depot"
  )
  # An optional element present in every column.
  expect_identical(
    .doseCmtsOfCode(
      c(
        "ev <- tibble(evid = c(if (gj) 5L else integer(0), 1L, rep(0L, n)),",
        "  amt = c(if (gj) 0 else numeric(0), 50, rep(NA_real_, n)),",
        '  cmt = c(if (gj) "cyp3a4" else character(0), "depot", rep("central", n)))'
      )
    ),
    "depot"
  )
  expect_identical(
    .doseCmtsOfCode('ev <- tibble::tribble(~time, ~amt, ~evid, ~cmt, 0, 0.5, 1L, "depot", 12, NA, 0L, "central")'),
    "depot"
  )
})

test_that("a compartment is read through pipes, variables and helper functions", {
  expect_identical(
    .doseCmtsOfCode('ev <- tidyr::crossing(cmt = c("depot", "depot2")) |> dplyr::mutate(time = 0, amt = d, evid = 1L)'),
    c("depot", "depot2")
  )
  expect_identical(.doseCmtsOfCode('ev <- d %>% dplyr::mutate(amt = d, evid = 1L, cmt = "depot")'), "depot")
  # An event table's own `cmt` without `amt` is an observation, not the
  # compartment of an `amt` added downstream.
  expect_identical(
    .doseCmtsOfCode(
      c(
        'ev <- as.data.frame(rxode2::et(rxode2::et(amt = 1, cmt = "stomach"), 1:2, cmt = "venous")) |>',
        "  dplyr::mutate(amt = amt * d)"
      )
    ),
    "stomach"
  )
  expect_identical(
    .doseCmtsOfCode(c('route_cmt <- if (iv) "venous" else "stomach"', "ev <- rxode2::et(amt = 1, cmt = route_cmt)")),
    c("stomach", "venous")
  )
  expect_identical(
    .doseCmtsOfCode(c("slot <- ifelse(iv == 1, 6L, 1L)", "ev <- data.frame(amt = 1, evid = 1L, cmt = slot)")),
    c("#1", "#6")
  )
  expect_identical(
    .doseCmtsOfCode(c(
      'arms <- tibble::tribble(~arm, ~cmt, "iv", "central", "im", "depot")',
      "row <- arms[i, ]",
      "ev <- data.frame(amt = 1, evid = 1L, cmt = row$cmt)"
    )),
    # Every string of the table; those that are not states of the model are
    # dropped when the doses are checked.
    c("central", "depot", "im", "iv")
  )
  expect_identical(.doseCmtsOfCode('ev <- data.frame(time = 0, amt = 5); ev$cmt <- "depot2"'), "depot2")
  # A subscript names no compartment; the vector it subsets does.
  expect_identical(.doseCmtsOfCode('ev <- rxode2::et(amt = 1, cmt = d["cmt"])'), character(0))
  expect_identical(
    .doseCmtsOfCode('ev <- rxode2::et(amt = 1, cmt = c("depot", "central")[[k]])'),
    c("central", "depot")
  )
  expect_identical(
    .doseCmtsOfCode(c(
      "run <- function(mod, amt, cmt) rxode2::rxSolve(mod, rxode2::et(amt = amt, cmt = cmt))",
      'a <- run(m, 1, cmt = "depot_ip")',
      'b <- run(m, 1, "central")'
    )),
    c("central", "depot_ip")
  )
})

test_that("each dose is credited to the model it is solved with", {
  d <- .dosesOfCode(
    c(
      'm1 <- readModelDb("A")',
      'm2 <- readModelDb("B")',
      'ev <- rxode2::et(amt = 1, cmt = "depot")',
      "s <- rxode2::rxSolve(m2, ev)",
      's1 <- rxode2::rxSolve(m1, rxode2::et(amt = 1, cmt = "central"))',
      'ev2 <- rxode2::et(amt = 1, cmt = "depot2")'
    ),
    models = c("A", "B")
  )
  expect_identical(d$model, c("A", "A", "B"))
  expect_identical(d$cmt, c("central", "depot2", "depot"))
  expect_identical(d$owner, c("referenced", "recent", "used"))
  # Before any model is named, a dose is credited to every model.
  d <- .dosesOfCode('ev <- rxode2::et(amt = 1, cmt = "depot")', models = c("A", "B"))
  expect_identical(d$model, c("A", "B"))
  expect_identical(d$owner, c("all", "all"))
  # A dose solved with a model the article defines itself is credited to none,
  # even when its parameters come from a registry model.
  d <- .dosesOfCode(
    c(
      'm1 <- readModelDb("A")',
      "th <- m1$theta",
      "ref <- rxode2::rxode2({ ntr <- 1; d/dt(depot_brain) <- -depot_brain })",
      'ev <- rxode2::et(amt = 1, cmt = "depot_brain")',
      "r <- rxode2::rxSolve(ref, params = th, events = ev)"
    ),
    models = "A"
  )
  expect_identical(nrow(d), 0L)
})

test_that("a model whose article is missing is an error, not a skip", {
  articleDir <- file.path(tempfile(), "articles")
  dir.create(articleDir, recursive = TRUE)
  on.exit(unlink(dirname(articleDir), recursive = TRUE))
  writeLines(c("```{r}", "x <- 1", "```"), file.path(articleDir, "Here.Rmd"))
  db <- data.frame(name = c("A", "B", "C"), vignette = c("Here", "Gone", NA), dosing = NA_character_)
  noExempt <- data.frame(model = character(0), state = character(0))
  expect_error(articleDoseGaps(articleDir, db, noExempt), "names no article.*: B$")
  expect_identical(nrow(articleDoseGaps(articleDir, db, noExempt, noArticle = "B")$gaps), 0L)
  # A model listed as having no article must still lack one.
  expect_error(articleDoseGaps(articleDir, db, noExempt, noArticle = c("A", "B")), "whose article exists: A$")
  # An article in the vignettes directory above the articles counts.
  writeLines(c("```{r}", "x <- 1", "```"), file.path(dirname(articleDir), "Gone.Rmd"))
  expect_identical(nrow(articleDoseGaps(articleDir, db, noExempt)$gaps), 0L)
})

test_that("a dosed name that is a state of no registry model is reported", {
  skip_on_cran()
  articleDir <- file.path(tempfile(), "articles")
  dir.create(articleDir, recursive = TRUE)
  on.exit(unlink(dirname(articleDir), recursive = TRUE))
  writeLines(
    c(
      "```{r}",
      'm <- readModelDb("PK_1cmt")',
      's1 <- rxode2::rxSolve(m, rxode2::et(amt = 1, cmt = "Cc"))',
      's2 <- rxode2::rxSolve(m, rxode2::et(amt = 1, cmt = "centrl"))',
      's3 <- rxode2::rxSolve(m, rxode2::et(amt = 1, cmt = "central"))',
      "```"
    ),
    file.path(articleDir, "Typo.Rmd")
  )
  db <- data.frame(name = "PK_1cmt", vignette = "Typo", dosing = "depot,central")
  res <- articleDoseGaps(articleDir, db, data.frame(model = character(0), state = character(0)))
  expect_identical(nrow(res$gaps), 0L)
  expect_identical(sort(res$unknown$cmt), c("Cc", "centrl"))
})

test_that("every compartment an article doses is in its model's `dosing`", {
  skip_on_cran()
  articleDir <- testthat::test_path("..", "..", "vignettes", "articles")
  skip_if(!dir.exists(articleDir), "vignettes/articles not present (installed package)")
  exempt <- rbind(articleStateDoses, articleMiscreditedDoses)
  res <- articleDoseGaps(articleDir, nlmixr2lib::modeldb, exempt, noArticle = articleMissingArticles$model)
  expect_identical(
    nrow(res$gaps),
    0L,
    info = paste(
      sprintf(
        "%s: %s doses %s (%s)",
        res$gaps$article,
        res$gaps$model,
        res$gaps$state,
        substr(res$gaps$record, 1, 120)
      ),
      collapse = "\n"
    )
  )
  # An exemption the articles no longer need, or whose state `dosing` now
  # lists, is stale.
  stale <- exempt[!paste(exempt$model, exempt$state) %in% paste(res$exempted$model, res$exempted$state), ]
  expect_identical(nrow(stale), 0L, info = paste(stale$model, stale$state, collapse = "\n"))
  # A dosed name that is a state of no registry model is a typo or an
  # observable, unless the reader is known to misread that record.
  found <- paste(res$unknown$article, res$unknown$cmt)
  listed <- paste(articleUnpairedRecords$article, articleUnpairedRecords$cmt)
  unknown <- res$unknown[!found %in% listed, ]
  expect_identical(
    nrow(unknown),
    0L,
    info = paste(
      sprintf("%s doses %s (%s)", unknown$article, unknown$cmt, substr(unknown$record, 1, 120)),
      collapse = "\n"
    )
  )
  expect_identical(setdiff(listed, found), character(0))
})
