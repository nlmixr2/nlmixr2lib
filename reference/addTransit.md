# To add transit compartments to the model

To add transit compartments to the model

## Usage

``` r
addTransit(
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

  the transit compartment prefix

- ktr:

  the parameter name for the transit compartment rate

- ka:

  absorption rate parameter name

## Value

a model with transit compartment added

This matches

\`dose-\>a0-\>a1-\>abs cmt-\>central\`

But \`a0\` is depot so dosing records labeled depot do not need to be
changed

The abs cmt becomes the last "transit" compartment

This is simply for convenience

## See also

Other absorption:
[`addWeibullAbs()`](https://nlmixr2.github.io/nlmixr2lib/reference/addWeibullAbs.md),
[`removeTransit()`](https://nlmixr2.github.io/nlmixr2lib/reference/removeTransit.md)

## Examples

``` r
readModelDb("PK_1cmt_des") |> addTransit(3)
#>  
#>  
#>  ── rxode2-based free-form 5-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka    lcl    lvc propSd   lktr 
#>   0.45   1.00   3.45   0.50   0.10 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2         transit1
#> 3                  3         transit2
#> 4                  4         transit3
#> 5                  5          central
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
#>         d/dt(transit2) <- ktr * transit1 - ktr * transit2
#>         d/dt(transit3) <- ktr * transit2 - ka * transit3
#>         d/dt(central) <- ka * transit3 - kel * central
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }
```
