# To convert from infusion/intravenous administration to first-order oral absorption

To convert from infusion/intravenous administration to first-order oral
absorption

## Usage

``` r
addDepot(ui, central = "central", depot = "depot", ka = "ka")
```

## Arguments

- ui:

  The model as a function (or something convertible to an rxUi object)

- central:

  central compartment name

- depot:

  depot compartment name

- ka:

  absorption rate parameter name

## Value

a model with the depot added

## Examples

``` r
# most of the examples in the model library already have a depot
# the PK_2cmt_no_depot is an exception
readModelDb("PK_2cmt_no_depot")  |> addDepot()
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lcl    lvc    lvp     lq propSd    lka 
#>    1.0    3.0    5.0    0.1    0.5    0.1 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#> 3                  3      peripheral1
#>  ── Model (Normalized Syntax): ── 
#> function() {
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
#>         lka <- 0.1
#>         label("First order absorption rate (ka)")
#>     })
#>     model({
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         vp <- exp(lvp)
#>         q <- exp(lq)
#>         kel <- cl/vc
#>         k12 <- q/vc
#>         k21 <- q/vp
#>         d/dt(depot) <- -ka * depot
#>         d/dt(central) <- -kel * central - k12 * central + k21 * 
#>             peripheral1 + ka * depot
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }
```
