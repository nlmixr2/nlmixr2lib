#!/usr/bin/env Rscript
# Parallel pkgdown site build for nlmixr2lib's GitHub Actions workflow.
#
# Temporary workaround until pkgdown gains native worker support upstream.
# When that PR lands and is released, this script can be deleted and the
# pkgdown.yaml workflow reverted to a single build_site_github_pages() call.
#
# Strategy: render every vignette in parallel via a hand-rolled scheduler over
# parallel::mcparallel + mccollect (Linux fork), with FAIL-FAST cancellation
# on the first vignette failure (sends SIGTERM to all in-flight siblings).
# After the parallel pass, call pkgdown::build_site_github_pages(lazy = TRUE,
# clean = FALSE) so pkgdown handles the non-article sections (reference, news,
# search, ...) and skips already-rendered articles via its mtime lazy check.
#
# The PKGDOWN_WORKERS env var is the single managed knob. Blank => use
# parallel::detectCores(); a positive integer overrides (e.g. to cap memory).

env <- Sys.getenv("PKGDOWN_WORKERS", "")
ncores <- if (nzchar(env)) as.integer(env) else parallel::detectCores()
stopifnot(is.finite(ncores), ncores >= 1L)
Sys.setenv(PKGDOWN_WORKERS = ncores)

# Single-syscall logger; writes < PIPE_BUF (4 KB) are atomic on Linux, so
# concurrent workers will not interleave mid-line.
log_line <- function(...) {
  cat(sprintf("[%s] %s\n",
              format(Sys.time(), "%H:%M:%S"),
              paste0(...)),
      file = stderr())
}

log_line("PKGDOWN_WORKERS = ", ncores,
         "  (detectCores() = ", parallel::detectCores(), ")")

pkg <- pkgdown::as_pkgdown(".")

# init_site() copies assets/CSS into docs/. Required because the wrap-up
# build_site_github_pages(lazy = TRUE) call skips init_site when lazy.
pkgdown::init_site(pkg)

articles <- pkg$vignettes$name[pkg$vignettes$type == "rmd"]
log_line("Rendering ", length(articles), " articles in parallel...")
t0 <- Sys.time()

# Spawn a forked worker rendering one article. Returns the mcparallel job.
launch_one <- function(name) {
  parallel::mcparallel({
    t_start <- Sys.time()
    tryCatch({
      pkgdown::build_article(name = name, pkg = pkg,
                             new_process = FALSE, quiet = TRUE)
      list(ok = TRUE, name = name,
           secs = as.numeric(difftime(Sys.time(), t_start, units = "secs")))
    }, error = function(e) {
      out <- pkg$vignettes$file_out[match(name, pkg$vignettes$name)]
      try(unlink(file.path(pkg$dst_path, out)), silent = TRUE)
      list(ok = FALSE, name = name, err = conditionMessage(e))
    })
  }, name = name)
}

pending <- articles
# running: vignette name -> mcparallel job descriptor.
# Keyed by vignette name because mccollect() returns its result list keyed by
# the job's 'name' argument (which we set to the vignette name); keying running
# the same way lets us simply do `running[[result_name]] <- NULL` to deregister
# a finished worker.
running <- list()

# Kill all in-flight workers and stop the script with a clear error.
abort_with <- function(failed_name, err) {
  if (length(running) > 0) {
    log_line("Killing ", length(running), " in-flight worker(s) and aborting")
    for (other_name in names(running)) {
      try(tools::pskill(running[[other_name]]$pid, tools::SIGTERM),
          silent = TRUE)
    }
    try(parallel::mccollect(unname(running), wait = TRUE), silent = TRUE)
  }
  stop("Vignette build failed: ", failed_name, " -- ", err, call. = FALSE)
}

while (length(running) < ncores && length(pending) > 0) {
  nm <- pending[1]
  pending <- pending[-1]
  job <- launch_one(nm)
  log_line("START ", nm, "  (pid ", job$pid, ")")
  running[[nm]] <- job
}

while (length(running) > 0) {
  res <- parallel::mccollect(unname(running), wait = FALSE, timeout = 0.5)
  if (is.null(res)) next
  for (article_name in names(res)) {
    r <- res[[article_name]]
    running[[article_name]] <- NULL
    failed <- is.null(r) ||
              inherits(r, "try-error") ||
              (is.list(r) && isFALSE(r$ok))
    if (failed) {
      err <- if (inherits(r, "try-error")) attr(r, "condition")$message
             else if (is.list(r) && !is.null(r$err)) r$err
             else "(silent crash)"
      log_line("FAIL  ", article_name, "  -- ", err)
      abort_with(article_name, err)
    }
    log_line("DONE  ", article_name, "  (", sprintf("%.1fs", r$secs), ")")
  }
  while (length(running) < ncores && length(pending) > 0) {
    nm <- pending[1]
    pending <- pending[-1]
    job <- launch_one(nm)
    log_line("START ", nm, "  (pid ", job$pid, ")")
    running[[nm]] <- job
  }
}

log_line("Parallel article render finished in ", format(Sys.time() - t0))

# Build remaining sections (reference, home, news, sitemap, search, navbar,
# GitHub-pages metadata). lazy = TRUE makes build_articles() short-circuit on
# articles whose HTML is newer than their Rmd; clean = FALSE preserves the
# articles we just rendered.
log_line("Building non-article sections (reference, news, sitemap, ...)")
t1 <- Sys.time()
pkgdown::build_site_github_pages(
  new_process = FALSE,
  install     = FALSE,
  clean       = FALSE,
  lazy        = TRUE
)
log_line("Wrap-up build finished in ", format(Sys.time() - t1))
log_line("Total pkgdown build time: ", format(Sys.time() - t0))
