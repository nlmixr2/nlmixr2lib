# Locale-independent ordering for anything that feeds a committed artifact.
#
# R collates with LC_COLLATE, and `sort()`, `order()`, `list.files()` and
# `Sys.glob()` all follow it. C and en_US.UTF-8 disagree about case: C puts
# every upper-case initial ahead of every lower-case one ("ABT 102 (Othman
# 2013)" < "Abacavir (Archary 2019)"), en_US interleaves them. The generated
# artifacts are therefore a different byte stream depending on the locale of
# the machine that regenerated them -- most of the `modeldb` rows and most of
# the `_pkgdown.yml` navbar entries move between the two -- which shows up as
# a huge spurious diff that hides the real change.
#
# `method = "radix"` collates in the C locale whatever LC_COLLATE says. Every
# ordering in the generator path goes through these wrappers, so a rebuild is
# a zero diff on any machine. `list.files()` and `Sys.glob()` sort internally
# and have no such option, so their results are re-sorted.
#
# The gate is `tests/testthat/test-collate.R`: it re-runs the generators under
# a non-C collation and requires byte-identical output, and it rejects a bare
# `sort()`/`order()`/`list.files()`/`Sys.glob()` reaching the generator path.
.collateSort <- function(x) {
  sort(x, method = "radix")
}

.collateOrder <- function(...) {
  order(..., method = "radix")
}

.collateListFiles <- function(...) {
  .collateSort(list.files(...))
}

.collateGlob <- function(paths) {
  .collateSort(Sys.glob(paths))
}
