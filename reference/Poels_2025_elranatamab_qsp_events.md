# Build an event table for the Poels 2025 elranatamab QSP model

\`Poels_2025_elranatamab_qsp\` reproduces the dose-to-dose attenuation
of cytokine release described by Poels et al. (2025) Supplementary Eq
37, which rescales the cumulative cytokine exposure state \`cauc\` at
the start of every dosing interval:

## Usage

``` r
Poels_2025_elranatamab_qsp_events(
  dose_time,
  dose_mg,
  obs_time,
  mw = 148500,
  id = 1L,
  cytokine_reset = TRUE,
  obs_cmt = "central"
)
```

## Arguments

- dose_time:

  Numeric vector of dose times, in hours.

- dose_mg:

  Numeric vector of subcutaneous doses, in mg. Either length 1
  (recycled) or the same length as \`dose_time\`.

- obs_time:

  Numeric vector of observation times, in hours.

- mw:

  Elranatamab molecular weight in g/mol. Defaults to 148500, the
  "approximately 148.5 kDa" of the ELREXFIO US prescribing information
  (Description section).

- id:

  Subject identifier written to the \`id\` column.

- cytokine_reset:

  Logical; when \`TRUE\` (the default) the Eq 37 reset records are
  included. Set to \`FALSE\` to simulate the model without the
  dose-to-dose attenuation of cytokine release.

- obs_cmt:

  Compartment used for observation records. Must be an ODE state of the
  model; every algebraic observable (\`Cc\`, \`mProtein\`,
  \`tumorBurden\`, ...) is returned as a column at those rows regardless
  of which state is named.

## Value

A \`data.frame\` with columns \`id\`, \`time\`, \`amt\`, \`evid\` and
\`cmt\`, ordered by time, suitable for \[rxode2::rxSolve()\].

## Details

\$\$C\_{auc,N}(0) = \frac{C\_{auc,N-1}(\tau)}{5\left(1 -
\frac{(N+1)^2}{1.3^2 + (N+1)^2}\right)}\$\$

where \\N\\ is the number of doses given. That is a discrete state
reset, which cannot be written inside an rxode2 ODE. This helper
therefore emits it as an \`evid = 6\` (multiply) record on \`cauc\` at
each dose time, carrying the reciprocal of the Eq 37 divisor as its
\`amt\`, alongside the subcutaneous dose records and the requested
observation records.

Doses are supplied in mg and converted to pmol using the elranatamab
molecular weight, because the model carries \`depot\` as a pmol amount
(Poels 2025 Supplementary Eq 1 divides the absorption flux by \`Vc\`).

## References

Poels KE, Elmeliegy M, Hibma J, Wang D, Musante CJ, Shtylla B.
Leveraging quantitative systems pharmacology modeling for elranatamab
regimen optimization in relapsed or refractory multiple myeloma. npj
Syst Biol Appl. 2025;11:102.
[doi:10.1038/s41540-025-00585-z](https://doi.org/10.1038/s41540-025-00585-z)

## Examples

``` r
# MagnetisMM-3 Cohort A: 12/32/76 mg on C1D1/C1D4/C1D8, then 76 mg weekly
dose_time <- c(0, 72, 168, seq(336, by = 168, length.out = 4))
dose_mg <- c(12, 32, rep(76, 5))
ev <- Poels_2025_elranatamab_qsp_events(
  dose_time = dose_time,
  dose_mg = dose_mg,
  obs_time = seq(0, 1344, by = 24)
)
head(ev)
#>   id time          amt evid     cmt
#> 1  1    0 8.080808e+04    1   depot
#> 2  1    0 6.733728e-01    6    cauc
#> 3  1    0           NA    0 central
#> 4  1   24           NA    0 central
#> 5  1   48           NA    0 central
#> 6  1   72 2.154882e+05    1   depot
```
