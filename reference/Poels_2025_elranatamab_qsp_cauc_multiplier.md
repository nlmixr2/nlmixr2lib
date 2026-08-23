# Eq 37 cumulative-cytokine-exposure multiplier

Reciprocal of the Poels 2025 Supplementary Eq 37 divisor, \`5 \* (1 -
(n + 1)^2 / (1.3^2 + (n + 1)^2))\`, evaluated at dose index \`n\`.
Exported so the reset applied by \[Poels_2025_elranatamab_qsp_events()\]
can be inspected and tested directly.

## Usage

``` r
Poels_2025_elranatamab_qsp_cauc_multiplier(n)
```

## Arguments

- n:

  Integer vector of dose indices (1 for the first dose).

## Value

Numeric vector of multipliers applied to the \`cauc\` state.

## Examples

``` r
Poels_2025_elranatamab_qsp_cauc_multiplier(1:5)
#> [1] 0.6733728 1.2650888 2.0934911 3.1585799 4.4603550
```
