# To convert from first order oral absorption to IV/Intravenous

To convert from first order oral absorption to IV/Intravenous

## Usage

``` r
removeDepot(ui, central = "central", depot = "depot", ka = "ka")
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

Returns a model with the depot from a first order absorption model
removed

## Examples

``` r
readModelDb("PK_1cmt_des") |> removeDepot()
#>  
#>  
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
#>     dosing <- c("central", "depot")
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
