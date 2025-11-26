# Converts first order absorption model to Weibull absorption model

Converts first order absorption model to Weibull absorption model

## Usage

``` r
addWeibullAbs(
  ui,
  ntransit,
  central = "central",
  depot = "depot",
  transit = "transit",
  wa = "wa",
  wb = "wb",
  ka = "ka",
  ktr = "ktr"
)
```

## Arguments

- ui:

  The model as a function

- ntransit:

  the number of transit compartments to be added

- central:

  central compartment name

- depot:

  depot compartment name

- transit:

  the transit compartment prefix

- wa:

  Weibull alpha parameter name

- wb:

  Weibull beta parameter name

- ka:

  absorption rate parameter name

- ktr:

  the parameter name for the transit compartment rate

## Value

model where first order absorption is changed to Weibull absorption
model

## See also

Other absorption:
[`addTransit()`](https://nlmixr2.github.io/nlmixr2lib/reference/addTransit.md),
[`removeTransit()`](https://nlmixr2.github.io/nlmixr2lib/reference/removeTransit.md)

## Author

Matthew L. Fidler

## Examples

``` r
readModelDb("PK_1cmt_des") |>
  addWeibullAbs()
#>  
#>  
#>  ── rxode2-based free-form 2-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lcl    lvc propSd    lwa    lwb 
#>   1.00   3.45   0.50   0.10   0.10 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     dosing <- c("central", "depot")
#>     ini({
#>         lcl <- 1
#>         label("Clearance (CL)")
#>         lvc <- 3.45
#>         label("Central volume of distribution (V)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lwa <- 0.1
#>         label("Weibull absorption alpha (wa)")
#>         lwb <- 0.1
#>         label("Weibull absorption beta (wa)")
#>     })
#>     model({
#>         wa <- exp(lwa)
#>         wb <- exp(lwb)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         kel <- cl/vc
#>         d/dt(depot) <- -(wb/wa) * (tad0(depot)/wa)^(wb - 1) * 
#>             depot
#>         d/dt(central) <- (wb/wa) * (tad0(depot)/wa)^(wb - 1) * 
#>             depot - kel * central
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }
```
