# Remove a peripheral compartment from a model

Removes \`peripheral2\` (or \`peripheral1\` when it is the only one),
rewiring the central compartment and dropping the associated parameters.
\`peripheral1\` cannot be removed while \`peripheral2\` is still
present.

## Usage

``` r
removePeriph(ui, n = NULL, central = "central", model)
```

## Arguments

- ui:

  The model as a function (or something convertible to an rxUi object)

- n:

  which peripheral to remove: 1 or 2; defaults to the highest-numbered
  one present

- central:

  central compartment name

- model:

  Deprecated alias for `ui`. Supplying `model` instead of `ui` still
  works but emits a deprecation warning.

## Value

a model with the peripheral compartment removed

## See also

Other distribution:
[`addPeriph()`](https://nlmixr2.github.io/nlmixr2lib/reference/addPeriph.md)

## Author

Matthew L. Fidler

## Examples

``` r

readModelDb("PK_2cmt_des") |> removePeriph()
#>  
#>  
#>  ── rxode2-based free-form 2-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka    lcl    lvc propSd 
#>   0.45   1.00   3.00   0.50 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     compartmentData <- list(depot = list(analyte = "drug", units = NA_character_, 
#>         specimen = "administration site", verified = FALSE), 
#>         central = list(analyte = "drug", units = NA_character_, 
#>             specimen = "plasma", verified = FALSE), peripheral1 = list(analyte = "drug", 
#>             units = NA_character_, specimen = "plasma", verified = FALSE))
#>     reference <- "nlmixr2lib template"
#>     units <- list(time = "time_unit", dosing = "dose_unit", concentration = "conc_unit/vol_unit")
#>     ini({
#>         lka <- 0.45
#>         label("Absorption rate (Ka)")
#>         lcl <- 1
#>         label("Clearance (CL)")
#>         lvc <- 3
#>         label("Central volume of distribution (V)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>     })
#>     model({
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         kel <- cl/vc
#>         d/dt(depot) <- -ka * depot
#>         d/dt(central) <- ka * depot - kel * central
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }
```
