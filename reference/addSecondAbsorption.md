# Add a second absorption path for double-absorption models

Adds a parallel absorption path (\`depot2\`) to a model that already has
first-order absorption through \`depot\`, for Monolix-style
double-absorption models: one dose record is split at translation time
into two parallel inputs with a logit-parameterized \`F1\`
apportionment. The second path can itself be zero-order (\`type =
"zero"\`, a modeled \`dur(depot2)\` input), first-order (\`type =
"first"\`, \`ka2\`), or a transit chain (\`delay = "transit"\`), with an
optional lag.

## Usage

``` r
addSecondAbsorption(
  ui,
  type = c("first", "zero"),
  delay = c("none", "lag", "transit"),
  n = NULL,
  central = "central",
  depot = "depot",
  depot2 = "depot2",
  ka = "ka",
  ka2 = "ka2",
  tk0 = "tk0",
  f1 = 0.7
)
```

## Arguments

- ui:

  The model as a function (or something convertible to an rxUi object)

- type:

  second-path input type: \`"zero"\` for a zero-order (modeled-duration)
  input or \`"first"\` for a first-order (\`ka2\`) input

- delay:

  second-path delay: \`"none"\`, \`"lag"\` (an \`alag\` on the second
  path), or \`"transit"\` (a transit chain feeding the second path)

- n:

  number of transit compartments when \`delay = "transit"\`

- central:

  central compartment name

- depot:

  depot compartment name

- depot2:

  name of the second depot compartment

- ka:

  absorption rate parameter name

- ka2:

  name of the second first-order absorption rate

- tk0:

  name of the zero-order duration on the second path

- f1:

  initial fraction of the dose entering the first path, in (0,1)

## Value

a model with two parallel absorption paths

## Details

Only one second absorption path is supported per model; applying this
twice raises an error.

## See also

Other absorption:
[`addTransit()`](https://nlmixr2.github.io/nlmixr2lib/reference/addTransit.md),
[`addWeibullAbs()`](https://nlmixr2.github.io/nlmixr2lib/reference/addWeibullAbs.md),
[`addZeroOrderAbs()`](https://nlmixr2.github.io/nlmixr2lib/reference/addZeroOrderAbs.md),
[`convertAbsForceLongerDelay()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertAbsForceLongerDelay.md),
[`convertAbsSequential()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertAbsSequential.md),
[`removeTransit()`](https://nlmixr2.github.io/nlmixr2lib/reference/removeTransit.md),
[`removeZeroOrderAbs()`](https://nlmixr2.github.io/nlmixr2lib/reference/removeZeroOrderAbs.md)

## Author

Matthew L. Fidler

## Examples

``` r

# simultaneous zero-order + first-order (Monolix double absorption)
readModelDb("PK_1cmt_des") |>
  addSecondAbsorption(type = "first", delay = "lag", f1 = 0.7)
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>         lka         lcl         lvc      propSd        lka2  llagDepot2 
#>   0.4500000   1.0000000   3.4500000   0.5000000   0.1000000   0.1000000 
#> logitfDepot 
#>   0.8472979 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2           depot2
#> 3                  3          central
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
#>         lka <- 0.45
#>         label("Absorption rate (Ka)")
#>         lcl <- 1
#>         label("Clearance (CL)")
#>         lvc <- 3.45
#>         label("Central volume of distribution (V)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lka2 <- 0.1
#>         label("First order absorption rate (ka2)")
#>         llagDepot2 <- 0.1
#>         logitfDepot <- 0.847297860387204
#>         label("Fraction of dose absorbed from depot (fDepot)")
#>     })
#>     model({
#>         fDepot <- expit(logitfDepot, 0, 1)
#>         lagDepot2 <- exp(llagDepot2)
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         kel <- cl/vc
#>         splitInfusionBolus(depot, depot, depot2)
#>         d/dt(depot) <- -ka * depot
#>         f(depot) <- fDepot
#>         ka2 <- exp(lka2)
#>         d/dt(depot2) <- -ka2 * depot2
#>         f(depot2) <- 1 - fDepot
#>         lag(depot2) <- lagDepot2
#>         d/dt(central) <- ka * depot - kel * central + ka2 * depot2
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }
```
