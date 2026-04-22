# Convert coefficient of variation (percent) to between-subject variance

\`cvToBsv()\` converts a log-normal coefficient of variation (expressed
as a percent) to the variance on the log scale, \`omega^2 =
log((CV/100)^2 + 1)\`. \`bsvToCv()\` is the inverse transform.

## Usage

``` r
cvToBsv(cv)

bsvToCv(bsv)
```

## Arguments

- cv:

  The coefficient of variation, expressed as a percent.

- bsv:

  The between-subject variability on the variance (omega^2) scale.

## Value

The between-subject variability on the variance (omega^2) scale.

The coefficient of variation on the percent scale.

## Functions

- `bsvToCv()`: Convert log-scale between-subject variance back to a
  percent coefficient of variation.
