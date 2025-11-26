# To remove transit compartments from the model

To remove transit compartments from the model

## Usage

``` r
removeTransit(
  ui,
  ntransit,
  central = "central",
  depot = "depot",
  transit = "transit",
  ktr = "ktr",
  ka = "ka"
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

  The number of transit compartments to remove

- ktr:

  the parameter name for the transit compartment rate

- ka:

  absorption rate parameter name

## Value

rxode2 model with transit compartment removed

## See also

Other absorption:
[`addTransit()`](https://nlmixr2.github.io/nlmixr2lib/reference/addTransit.md),
[`addWeibullAbs()`](https://nlmixr2.github.io/nlmixr2lib/reference/addWeibullAbs.md)

## Examples

``` r
# In this example the transit is added and then a few are removed

readModelDb("PK_1cmt_des") |> addTransit(4) |> removeTransit(3)
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka    lcl    lvc propSd   lktr 
#>   0.45   1.00   3.45   0.50   0.10 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2         transit1
#> 3                  3          central
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     dosing <- c("central", "depot")
#>     ini({
#>         lka <- 0.45
#>         label("Absorption rate (Ka)")
#>         lcl <- 1
#>         label("Clearance (CL)")
#>         lvc <- 3.45
#>         label("Central volume of distribution (V)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lktr <- 0.1
#>         label("First order transition rate (ktr)")
#>     })
#>     model({
#>         ktr <- exp(lktr)
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         kel <- cl/vc
#>         d/dt(depot) <- -ktr * depot
#>         d/dt(transit1) <- ktr * depot - ktr * transit1
#>         d/dt(central) <- ka * transit1 - kel * central
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

readModelDb("PK_1cmt_des") |> addTransit(4) |> removeTransit()
#>  
#>  
#>  ── rxode2-based free-form 2-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka    lcl    lvc propSd 
#>   0.45   1.00   3.45   0.50 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     dosing <- c("central", "depot")
#>     ini({
#>         lka <- 0.45
#>         label("Absorption rate (Ka)")
#>         lcl <- 1
#>         label("Clearance (CL)")
#>         lvc <- 3.45
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
