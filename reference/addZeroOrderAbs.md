# Convert a model to zero-order absorption (Monolix `Tk0`)

Zero-order absorption releases the dose into the central compartment at
a constant rate over a modeled duration `Tk0`, matching Monolix's
`absorption(type=2, Tk0)` oral route. Following monolix2rx's
translation, this is expressed as a modeled duration on the central
compartment (`dur(central) <- tk0`), so a depot compartment (and its
`ka`) is removed first. Dosing records for the zero-order route must
request the modeled duration in the event table (`RATE = -2`; see
`et(rate=-2)`); ordinary bolus records bypass the modeled duration.

## Usage

``` r
addZeroOrderAbs(
  ui,
  central = "central",
  depot = "depot",
  transit = "transit",
  ktr = "ktr",
  ka = "ka",
  tk0 = "tk0"
)
```

## Arguments

- ui:

  The model as a function (or something convertible to an rxUi object)

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

- tk0:

  zero-order absorption duration parameter name (Monolix's `Tk0`)

## Value

a model with zero-order absorption

## Details

Lag time and bioavailability can be combined with \[addLag()\] and
\[addLogitBioavailability()\].

## See also

Other absorption:
[`addSecondAbsorption()`](https://nlmixr2.github.io/nlmixr2lib/reference/addSecondAbsorption.md),
[`addTransit()`](https://nlmixr2.github.io/nlmixr2lib/reference/addTransit.md),
[`addWeibullAbs()`](https://nlmixr2.github.io/nlmixr2lib/reference/addWeibullAbs.md),
[`convertAbsForceLongerDelay()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertAbsForceLongerDelay.md),
[`convertAbsSequential()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertAbsSequential.md),
[`removeTransit()`](https://nlmixr2.github.io/nlmixr2lib/reference/removeTransit.md),
[`removeZeroOrderAbs()`](https://nlmixr2.github.io/nlmixr2lib/reference/removeZeroOrderAbs.md)

## Author

Matthew L. Fidler

## Examples

``` r

readModelDb("PK_1cmt_des") |> addZeroOrderAbs()
#>  
#>  
#> Warning: 'depot' removed for zero-order absorption model
#>  ── rxode2-based free-form 1-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lcl    lvc propSd   ltk0 
#>   1.00   3.45   0.50   0.10 
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
#>         ltk0 <- 0.1
#>         label("Zero-order absorption duration (Tk0)")
#>     })
#>     model({
#>         tk0 <- exp(ltk0)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         kel <- cl/vc
#>         d/dt(central) <- -kel * central
#>         dur(central) <- tk0
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

# a model without a depot gets the same modeled duration input
readModelDb("PK_2cmt_no_depot") |> addZeroOrderAbs()
#>  
#>  
#>  ── rxode2-based free-form 2-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lcl    lvc    lvp     lq propSd   ltk0 
#>    1.0    3.0    5.0    0.1    0.5    0.1 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1          central
#> 2                  2      peripheral1
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     compartmentData <- list(central = list(analyte = "drug", 
#>         units = NA_character_, specimen = "plasma", verified = FALSE), 
#>         peripheral1 = list(analyte = "drug", units = NA_character_, 
#>             specimen = "plasma", verified = FALSE))
#>     reference <- "nlmixr2lib template"
#>     units <- list(time = "time_unit", dosing = "dose_unit", concentration = "conc_unit/vol_unit")
#>     ini({
#>         lcl <- 1
#>         label("Clearance (CL)")
#>         lvc <- 3
#>         label("Central volume of distribution (V)")
#>         lvp <- 5
#>         label("Peripheral volume of distribution (Vp)")
#>         lq <- 0.1
#>         label("Intercompartmental clearance (Q)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         ltk0 <- 0.1
#>         label("Zero-order absorption duration (Tk0)")
#>     })
#>     model({
#>         tk0 <- exp(ltk0)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         vp <- exp(lvp)
#>         q <- exp(lq)
#>         kel <- cl/vc
#>         k12 <- q/vc
#>         k21 <- q/vp
#>         d/dt(central) <- -kel * central - k12 * central + k21 * 
#>             peripheral1
#>         dur(central) <- tk0
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }
```
