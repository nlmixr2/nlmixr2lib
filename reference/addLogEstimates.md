# Add log estimates to a model

Add log estimates to a model

## Usage

``` r
addLogEstimates(ui, vars, extraLines = NULL, beforeCmt = NULL)
```

## Arguments

- ui:

  rxode2 model

- vars:

  estimates to add they will be parameterized as:

  `var <- exp(lvar)`

  where `var` is the variable name in the model and `lvar` is the log
  transformed variable that will be estimated

- extraLines:

  this is a list of additional lines to add to the model just after the
  variables are defined. It must be `NULL` or a list of `language`
  objects.

- beforeCmt:

  if the model is compartmental you can specify the preferred names
  where the estimates and extra lines are added before

## Value

rxode2 model with log estimates added (and possibly extra lines)

## Author

Matthew L. Fidler

## Examples

``` r
# Change the transformation of the PK model from cl to k

readModelDb("PK_3cmt_des") |>
  removeLinesAndInis(c("kel", "k12", "k21", "k13", "k31", "vc")) |>
  addLogEstimates(c("kel", "k12", "k21", "k13", "k31", "vc"))
#>  
#>  
#>  ── rxode2-based free-form 4-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka propSd   lkel   lk12   lk21   lk13   lk31    lvc 
#>   0.45   0.50   0.10   0.10   0.10   0.10   0.10   0.10 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#> 3                  3      peripheral1
#> 4                  4      peripheral2
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     ini({
#>         lka <- 0.45
#>         label("Absorption rate (Ka)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lkel <- 0.1
#>         lk12 <- 0.1
#>         lk21 <- 0.1
#>         lk13 <- 0.1
#>         lk31 <- 0.1
#>         lvc <- 0.1
#>     })
#>     model({
#>         kel <- exp(lkel)
#>         k12 <- exp(lk12)
#>         k21 <- exp(lk21)
#>         k13 <- exp(lk13)
#>         k31 <- exp(lk31)
#>         vc <- exp(lvc)
#>         ka <- exp(lka)
#>         d/dt(depot) <- -ka * depot
#>         d/dt(central) <- ka * depot - kel * central - k12 * central + 
#>             k21 * peripheral1 - k13 * central + k31 * peripheral2
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         d/dt(peripheral2) <- k13 * central - k31 * peripheral2
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

# You can also label the parameters by using a named character
# vector with the names of the parameters representing the
# variables and the values representing the labels:

readModelDb("PK_3cmt_des") |>
  removeLinesAndInis(c("kel", "k12", "k21", "k13", "k31", "vc")) |>
  addLogEstimates(c(kel="elimination", k12="k12 constant",
                    k21="k21 constant",
                    k13="k13 constant",
                    k31="k31 constant",
                    vc="volume of central compartment"))
#>  
#>  
#>  ── rxode2-based free-form 4-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka propSd   lkel   lk12   lk21   lk13   lk31    lvc 
#>   0.45   0.50   0.10   0.10   0.10   0.10   0.10   0.10 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#> 3                  3      peripheral1
#> 4                  4      peripheral2
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     ini({
#>         lka <- 0.45
#>         label("Absorption rate (Ka)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lkel <- 0.1
#>         label("elimination")
#>         lk12 <- 0.1
#>         label("k12 constant")
#>         lk21 <- 0.1
#>         label("k21 constant")
#>         lk13 <- 0.1
#>         label("k13 constant")
#>         lk31 <- 0.1
#>         label("k31 constant")
#>         lvc <- 0.1
#>         label("volume of central compartment")
#>     })
#>     model({
#>         kel <- exp(lkel)
#>         k12 <- exp(lk12)
#>         k21 <- exp(lk21)
#>         k13 <- exp(lk13)
#>         k31 <- exp(lk31)
#>         vc <- exp(lvc)
#>         ka <- exp(lka)
#>         d/dt(depot) <- -ka * depot
#>         d/dt(central) <- ka * depot - kel * central - k12 * central + 
#>             k21 * peripheral1 - k13 * central + k31 * peripheral2
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         d/dt(peripheral2) <- k13 * central - k31 * peripheral2
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }
```
