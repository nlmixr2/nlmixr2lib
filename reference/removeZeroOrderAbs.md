# Remove zero-order absorption from a model

Removes the modeled duration on the central compartment, leaving the
dose as an intravenous bolus input. The duration variable is taken from
the right hand side of the `dur(central)` line: when it is defined by an
assignment in the model (the `tk0 <- exp(ltk0)` added by
\[addZeroOrderAbs()\] or the `durCentral` added by \[addDur()\]) the
variable and its initial estimate are dropped too; a bare estimated
parameter that is not used elsewhere in the model is dropped from the
initial estimates. First-order absorption can be restored with
\[addDepot()\].

## Usage

``` r
removeZeroOrderAbs(ui, central = "central")
```

## Arguments

- ui:

  The model as a function (or something convertible to an rxUi object)

- central:

  central compartment name

## Value

a model where the zero-order absorption is removed

## See also

Other absorption:
[`addSecondAbsorption()`](https://nlmixr2.github.io/nlmixr2lib/reference/addSecondAbsorption.md),
[`addTransit()`](https://nlmixr2.github.io/nlmixr2lib/reference/addTransit.md),
[`addWeibullAbs()`](https://nlmixr2.github.io/nlmixr2lib/reference/addWeibullAbs.md),
[`addZeroOrderAbs()`](https://nlmixr2.github.io/nlmixr2lib/reference/addZeroOrderAbs.md),
[`convertAbsForceLongerDelay()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertAbsForceLongerDelay.md),
[`convertAbsSequential()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertAbsSequential.md),
[`removeTransit()`](https://nlmixr2.github.io/nlmixr2lib/reference/removeTransit.md)

## Author

Matthew L. Fidler

## Examples

``` r

readModelDb("PK_1cmt_des") |> addZeroOrderAbs() |> removeZeroOrderAbs()
#>  
#>  
#> Warning: 'depot' removed for zero-order absorption model
#>  ── rxode2-based free-form 1-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lcl    lvc propSd 
#>   1.00   3.45   0.50 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1          central
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     compartmentData <- list(depot = list(analyte = "drug", units = NA_character_, 
#>         specimen = "administration site", verified = FALSE), 
#>         central = list(analyte = "drug", units = NA_character_, 
#>             specimen = "plasma", verified = FALSE))
#>     dosing <- c("central", "depot")
#>     reference <- "nlmixr2lib template"
#>     units <- list(time = "time_unit", dosing = "dose_unit", concentration = "conc_unit/vol_unit")
#>     ini({
#>         lcl <- 1
#>         label("Clearance (CL)")
#>         lvc <- 3.45
#>         label("Central volume of distribution (V)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>     })
#>     model({
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         kel <- cl/vc
#>         d/dt(central) <- -kel * central
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }
```
