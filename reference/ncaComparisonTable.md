# Build a side-by-side NCA comparison table (simulated vs. reference)

Combines model-predicted NCA results with a reference set (e.g. paper
values, prior estimates, alternative model runs, observed data) into a
single tidy comparison frame ready for \`knitr::kable()\` or interactive
review. Each row is one NCA parameter \\\times\\ optional grouping
level; columns show the reference value, simulated value, and percent
difference. Rows whose discrepancy exceeds \`tolerance_pct\` are flagged
with a trailing asterisk and the function attaches a footnote string to
the result.

## Usage

``` r
ncaComparisonTable(
  simulated,
  reference,
  by = NULL,
  params = NULL,
  tolerance_pct = 20,
  units = NULL,
  label_first_column = "NCA parameter"
)
```

## Arguments

- simulated:

  Model-predicted NCA. One of: a \`PKNCAresults\` object, a
  \`data.frame\` with columns \`PPTESTCD\` and \`PPORRES\` plus any
  grouping columns named in \`by\`, or a wide \`data.frame\` (one row
  per group, one column per PKNCA code).

- reference:

  Reference NCA. Same shapes as \`simulated\`. A wide \`data.frame\` is
  the typical form when transcribing values from a paper table.

- by:

  Optional character vector of grouping column name(s) (e.g.
  \`"treatment"\`, \`c("cohort","weight_band")\`). When \`NULL\`, the
  inputs are treated as a single ungrouped comparison.

- params:

  Optional character vector of PKNCA codes restricting which parameters
  appear in the output. Defaults to the intersection of codes present in
  both inputs.

- tolerance_pct:

  Numeric scalar; rows with \`\| attaches a footnote string via the
  \`"footnote"\` attribute on the result. Pass \`Inf\` to disable
  flagging.

- units:

  Optional named character vector keyed by PKNCA code, passed through to
  \[ncaParamLabel()\].

- label_first_column:

  Header for the parameter column. Defaults to \`"NCA parameter"\`.

## Value

A \`data.frame\` whose first column is the friendly parameter label
(under \`label_first_column\`); when \`by\` is non-NULL, the grouping
column(s) follow; the final three columns are \`Reference\`,
\`Simulated\`, and \` parameter order then by group. Carries a
\`"footnote"\` attribute when any row exceeds the tolerance.

## Details

Designed for any nlmixr2 user comparing predicted NCA to a reference –
not solely a vignette utility. Accepts a \`PKNCAresults\` object, the
\`\$result\` data.frame, or a tidy long frame with columns \`PPTESTCD\`
and \`PPORRES\` for \`simulated\`; the same shapes plus a wide
\`data.frame\` (one row per group, one column per PKNCA parameter code)
for \`reference\`. Multiple rows per group in \`simulated\` (the typical
per-subject case from \[PKNCA::pk.nca()\]) are aggregated to the group
level using \[stats::median()\] before joining; pass a pre-aggregated
frame if a different summary statistic is required.

## Author

Bill Denney

## Examples

``` r
simulated <- data.frame(
  treatment = rep(c("50 mg", "100 mg"), each = 3),
  PPTESTCD  = rep(c("cmax", "tmax", "auclast"), 2),
  PPORRES   = c(15.2, 2.0, 96.0, 29.1, 2.1, 191.0)
)
reference <- data.frame(
  treatment = c("50 mg", "100 mg"),
  cmax      = c(14.8, 28.5),
  tmax      = c(2.0,  2.1),
  auclast   = c(95.0, 190.0)
)
tbl <- ncaComparisonTable(
  simulated, reference,
  by = "treatment",
  units = c(cmax = "ug/mL", auclast = "ug*h/mL", tmax = "h")
)
tbl
#>       NCA parameter treatment Reference Simulated % diff
#> 1      Cmax (ug/mL)     50 mg      14.8      15.2  +2.7%
#> 2      Cmax (ug/mL)    100 mg      28.5      29.1  +2.1%
#> 3          Tmax (h)     50 mg         2         2  +0.0%
#> 4          Tmax (h)    100 mg       2.1       2.1  +0.0%
#> 5 AUClast (ug*h/mL)     50 mg        95        96  +1.1%
#> 6 AUClast (ug*h/mL)    100 mg       190       191  +0.5%
attr(tbl, "footnote")
#> NULL
```
