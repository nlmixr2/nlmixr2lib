# Fake blank Cc for creating PD only models

Fake blank Cc for creating PD only models

## Usage

``` r
fakeCc(fun, ..., cc = "Cc")
```

## Arguments

- fun:

  function to that requires Cc

- ...:

  arguments sent PD function

- cc:

  character name of the concentration in the central compartment that
  will be faked to allow models that require Cc to change to models with
  Cc as a covariate

## Value

Model where Cc is a covariate

## Author

Matthew L. Fidler

## Examples

``` r
fakeCc(addDirectLin) |> convertEmaxHill()
#>  
#>  
#>  ── rxode2-based Pred model ───────────────────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>  effectSd     lEmax     lEC50       lgg 
#>  0.100000  0.100000  0.100000 -2.302585 
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     ini({
#>         effectSd <- c(0, 0.1)
#>         label("additive error for effect")
#>         lEmax <- 0.1
#>         label("Maximum effect (Emax)")
#>         lEC50 <- 0.1
#>         label("Concentration of 50% Emax (Emax)")
#>         lgg <- -2.30258509299405
#>         label("logit-constrained Hill coefficient g")
#>     })
#>     model({
#>         Emax <- exp(lEmax)
#>         EC50 <- exp(lEC50)
#>         g <- expit(lgg, 0.1, 10)
#>         effect <- Emax * Cc^g/(Cc^g + EC50^g)
#>         effect ~ add(effectSd)
#>     })
#> }
```
