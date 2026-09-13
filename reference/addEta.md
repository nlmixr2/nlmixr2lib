# Add random effects to a model

Add random effects to a model

## Usage

``` r
addEta(
  ui,
  eta,
  priorName = getOption("nlmixr2lib.priorEta", TRUE),
  etaCombineType = c("default", "snake", "camel", "dot", "blank"),
  model
)
```

## Arguments

- ui:

  The model as a function

- eta:

  vector with the parameters to add random effects (sometimes referred
  to as inter-individual variability, IIV) on

- priorName:

  logical, if TRUE, the parameter name specified in \`eta\` will be used
  to add the eta value prior name is used instead of the left handed
  side of the equation.

- etaCombineType:

  the option for the how to combine the eta with the parameter name. Can
  be: "default", "snake", "camel", "dot", "blank"

- model:

  Deprecated alias for `ui`. Supplying `model` instead of `ui` still
  works but emits a deprecation warning.

## Value

The model with eta added to the requested parameters

## Author

Bill Denney, Richard Hooijmaijers & Matthew L. Fidler

## Examples

``` r
library(rxode2)
#> rxode2 5.1.6 using 2 threads (see ?getRxThreads)
#>   no cache: create with `rxCreateCache()`
readModelDb("PK_1cmt") |> addEta("ka")
#>  
#>  
#>  
#>  
#> → Adding eta to lka instead of ka due to mu-referencing
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> ℹ promote `etaKa` to between subject variability with initial estimate 0.1
#> Error in rbind(deparse.level, ...): numbers of columns of arguments do not match
```
