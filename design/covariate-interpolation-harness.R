# usage: Rscript design/covariate-interpolation-harness.R <article-name> <out-dir>
# Environment: NLMIXR2LIB_DIR, the repository root (default: the working
# directory); NLMIXR2LIB_LIB, an R library holding nlmixr2lib built from that
# root (optional, put first on the library path).
#
# Runs one article's R code with every rxSolve() call routed through .rxS().
# .rxS() returns the article's own result unchanged. When a covariate the model
# reads varies within a subject in that call's events, it also solves up to 20
# such subjects twice, with covsInterpolation "locf" and "nocb", using the
# subjects' realised etas on a zeroRe() model so both solves are deterministic,
# and records how far the two outputs differ.
args <- commandArgs(TRUE)
art <- args[1]
outdir <- normalizePath(args[2])
WT <- normalizePath(Sys.getenv("NLMIXR2LIB_DIR", getwd()))
lib <- Sys.getenv("NLMIXR2LIB_LIB")
if (nzchar(lib)) {
  .libPaths(c(lib, .libPaths()))
}
grDevices::pdf(NULL)
setwd(file.path(WT, "vignettes", "articles"))

.study <- new.env(parent = emptyenv())
.study$calls <- list()
.study$t0 <- Sys.time()

std_cols <- c(
  "id",
  "time",
  "evid",
  "amt",
  "cmt",
  "rate",
  "dur",
  "ii",
  "addl",
  "ss",
  "dv",
  "mdv",
  "cens",
  "limit",
  "dvid",
  "sim.id",
  "resetno",
  "tad",
  "nlmixrrowid"
)

is_ev <- function(x) {
  (is.data.frame(x) || inherits(x, "rxEt")) &&
    any(tolower(names(as.data.frame(x))) == "time") &&
    any(tolower(names(as.data.frame(x))) %in% c("evid", "amt"))
}

find_events <- function(a) {
  if (!is.null(a$events) && is_ev(a$events)) {
    return("events")
  }
  nm <- names(a)
  if (is.null(nm)) {
    nm <- rep("", length(a))
  }
  for (k in seq_along(a)) {
    if (nm[k] %in% c("", "params") && is_ev(a[[k]])) return(k)
  }
  NULL
}

model_covs <- function(object) {
  mv <- tryCatch(rxode2::rxModelVars(object), error = function(e) NULL)
  if (is.null(mv)) {
    mv <- tryCatch(rxode2::rxModelVars(rxode2::rxode2(object)), error = function(e) NULL)
  }
  if (is.null(mv)) {
    return(character(0))
  }
  mv$params
}

varying <- function(ev, covs) {
  idc <- names(ev)[tolower(names(ev)) == "id"][1]
  if (is.na(idc)) {
    idc <- NULL
  }
  cols <- setdiff(intersect(names(ev), covs), names(ev)[tolower(names(ev)) %in% std_cols])
  out <- character(0)
  for (cc in cols) {
    v <- ev[[cc]]
    if (!is.numeric(v) && !is.logical(v) && !is.factor(v) && !is.character(v)) {
      next
    }
    g <- if (is.null(idc)) rep(1, nrow(ev)) else ev[[idc]]
    nd <- tapply(v, g, function(z) length(unique(z[!is.na(z)])))
    if (any(nd > 1, na.rm = TRUE)) out <- c(out, cc)
  }
  list(cols = out, idc = idc)
}

# Time from each change of a varying covariate back to the record before it,
# per subject: how long the two conventions disagree for.
change_gaps <- function(ev, cols, idc, tcol) {
  g <- if (is.null(idc)) rep(1, nrow(ev)) else ev[[idc]]
  gaps <- c()
  for (s in unique(g)) {
    e <- ev[g == s, , drop = FALSE]
    e <- e[order(e[[tcol]]), , drop = FALSE]
    if (nrow(e) < 2) {
      next
    }
    for (cc in cols) {
      v <- e[[cc]]
      for (k in 2:nrow(e)) {
        if (!is.na(v[k]) && !is.na(v[k - 1]) && !identical(v[k], v[k - 1]) && e[[tcol]][k] > e[[tcol]][k - 1]) {
          gaps <- c(gaps, e[[tcol]][k] - e[[tcol]][k - 1])
        }
      }
    }
  }
  gaps
}

main_col <- function(d) {
  for (c in c("Cc", "Cp", "cp", "CP", "C", "conc", "Conc", "Ccentral")) {
    if (c %in% names(d)) return(c)
  }
  NA_character_
}

solve_with <- function(model, a, interp) {
  rxode2::rxSetSeed(1)
  as.data.frame(suppressWarnings(do.call(rxode2::rxSolve, c(list(model), a, list(covsInterpolation = interp)))))
}

subject_auc <- function(d, mc) {
  xid <- if ("id" %in% names(d)) d$id else rep(1, nrow(d))
  vapply(
    split(seq_len(nrow(d)), xid),
    function(i) {
      tt <- d$time[i]
      cc <- d[[mc]][i]
      if (length(tt) < 2) {
        return(NA_real_)
      }
      sum(diff(tt) * (utils::head(cc, -1) + utils::tail(cc, -1)) / 2)
    },
    0
  )
}

# x and y are the locf and nocb solves of the same subjects.
diff_summary <- function(x, y, ev, ids) {
  skip <- c(names(ev), std_cols, "sim", "ipredSim", "rxLambda", "rxYj", "rxLow", "rxHi")
  num <- setdiff(names(x)[vapply(x, is.numeric, TRUE)], skip)
  num <- num[!grepl("^(eta|rx|err)", num)]
  sc <- vapply(
    num,
    function(cn) {
      m <- max(abs(x[[cn]]), na.rm = TRUE)
      if (!is.finite(m) || m == 0) {
        return(0)
      }
      max(abs(x[[cn]] - y[[cn]]), na.rm = TRUE) / m
    },
    0
  )
  out <- list(
    status = "compared",
    n_ids = if (is.null(ids)) 1L else length(ids),
    worst_col = if (length(sc)) names(sc)[which.max(sc)] else NA,
    worst_scaled = if (length(sc)) max(sc) else NA
  )
  mc <- main_col(x)
  out$main_col <- mc
  if (!is.na(mc)) {
    out$main_scaled <- sc[[mc]]
    rr <- abs(subject_auc(y, mc) / subject_auc(x, mc) - 1)
    rr <- rr[is.finite(rr)]
    out$auc_rel_median <- if (length(rr)) stats::median(rr) else NA
    out$auc_rel_max <- if (length(rr)) max(rr) else NA
  }
  out
}

.study_compare <- function(object, res, ...) {
  a <- list(...)
  k <- find_events(a)
  if (is.null(k)) {
    return(list(status = "no events found"))
  }
  ev <- as.data.frame(a[[k]])
  covs <- model_covs(object)
  vr <- varying(ev, covs)
  rec <- list(
    status = "no within-subject variation",
    varying = vr$cols,
    article_interp = if (is.null(a$covsInterpolation)) "default (locf)" else a$covsInterpolation[1],
    n_rows = nrow(ev)
  )
  if (!length(vr$cols)) {
    return(rec)
  }
  idc <- vr$idc
  tcol <- names(ev)[tolower(names(ev)) == "time"][1]
  gaps <- change_gaps(ev, vr$cols, idc, tcol)
  rec$gap_median <- if (length(gaps)) stats::median(gaps) else NA
  rec$gap_max <- if (length(gaps)) max(gaps) else NA
  rec$n_changes <- length(gaps)
  # Up to 20 subjects in which a varying covariate changes.
  ids <- if (is.null(idc)) {
    NULL
  } else {
    chg <- unique(ev[[idc]][unlist(lapply(vr$cols, function(cc) {
      nd <- tapply(ev[[cc]], ev[[idc]], function(z) length(unique(z[!is.na(z)])))
      ev[[idc]] %in% names(nd)[nd > 1]
    }))])
    utils::head(chg, 20)
  }
  ev2 <- if (is.null(ids)) ev else ev[ev[[idc]] %in% ids, , drop = FALSE]
  a2 <- a
  a2[[k]] <- ev2
  a2$nSub <- NULL
  a2$nStud <- NULL
  a2$covsInterpolation <- NULL
  model2 <- object
  pr <- tryCatch(as.data.frame(res$params), error = function(e) NULL)
  eta_cols <- if (is.null(pr)) character(0) else grep("^eta", names(pr), value = TRUE)
  if (length(eta_cols) && !is.null(ids)) {
    model2 <- tryCatch(rxode2::zeroRe(object), error = function(e) NULL)
    if (is.null(model2)) {
      rec$status <- "random effects, zeroRe() failed"
      return(rec)
    }
    pidc <- names(pr)[tolower(names(pr)) == "id"][1]
    if (is.na(pidc) || "sim.id" %in% names(pr)) {
      rec$status <- "random effects, params not per id"
      return(rec)
    }
    pr2 <- pr[pr[[pidc]] %in% ids, , drop = FALSE]
    if (!is.null(a2$params) && is.data.frame(a2$params)) {
      p0 <- a2$params
      p0id <- names(p0)[tolower(names(p0)) == "id"][1]
      if (!is.na(p0id)) {
        p0 <- p0[p0[[p0id]] %in% ids, , drop = FALSE]
        pr2 <- merge(p0, pr2[, c(pidc, eta_cols), drop = FALSE], by.x = p0id, by.y = pidc)
      }
    }
    a2$params <- pr2
    a2$omega <- NULL
    a2$sigma <- NULL
    a2$thetaMat <- NULL
  } else if (!is.null(a2$params) && is.data.frame(a2$params) && !is.null(ids)) {
    p0id <- names(a2$params)[tolower(names(a2$params)) == "id"][1]
    if (!is.na(p0id)) a2$params <- a2$params[a2$params[[p0id]] %in% ids, , drop = FALSE]
  }
  t1 <- Sys.time()
  x <- tryCatch(solve_with(model2, a2, "locf"), error = function(e) conditionMessage(e))
  y <- tryCatch(solve_with(model2, a2, "nocb"), error = function(e) conditionMessage(e))
  rec$compare_secs <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  if (is.character(x) || is.character(y)) {
    rec$status <- "comparison solve failed"
    rec$error <- if (is.character(x)) x else y
    return(rec)
  }
  if (nrow(x) != nrow(y)) {
    rec$status <- "comparison rows differ"
    return(rec)
  }
  utils::modifyList(rec, diff_summary(x, y, ev, ids))
}

.rxS <- function(object, ...) {
  res <- rxode2::rxSolve(object, ...)
  rec <- tryCatch(.study_compare(object, res, ...), error = function(e) {
    list(status = "harness error", error = conditionMessage(e))
  })
  .study$calls[[length(.study$calls) + 1]] <- rec
  res
}

save_out <- function(status) {
  saveRDS(
    list(
      article = art,
      status = status,
      secs = as.numeric(difftime(Sys.time(), .study$t0, units = "secs")),
      calls = .study$calls
    ),
    file.path(outdir, paste0(art, ".rds"))
  )
}

rmd <- file.path(WT, "vignettes", "articles", paste0(art, ".Rmd"))
purled <- tempfile(fileext = ".R")
invisible(knitr::purl(rmd, output = purled, documentation = 0, quiet = TRUE))
code <- readLines(purled, warn = FALSE)
code <- gsub("rxode2::rxSolve\\b", ".rxS", code, perl = TRUE)
code <- gsub("(?<![[:alnum:]_.:])rxSolve(?=\\s*\\()", ".rxS", code, perl = TRUE)
writeLines(code, purled)
st <- tryCatch(
  {
    source(purled, local = globalenv(), echo = FALSE, max.deparse.length = 0)
    "ok"
  },
  error = function(e) paste("article error:", conditionMessage(e))
)
save_out(st)
cat(art, st, length(.study$calls), "calls\n")
