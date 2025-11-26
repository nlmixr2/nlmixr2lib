# Removes lines and inis from a model

Removes lines and inis from a model

## Usage

``` r
removeLinesAndInis(ui, vars)
```

## Arguments

- ui:

  A rxode2 model

- vars:

  A character vector of variables to remove

## Value

model with rxode2 lines and any estimates associate with lines removed

## Author

Matthew L. Fidler

## Examples

``` r
readModelDb("PK_3cmt_des") |> removeLinesAndInis(c("kel", "k12", "k21"))
#>  
#>  
#>  ── rxode2-based free-form 4-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka   lvp2    lq2 propSd 
#>   0.45   8.00   0.50   0.50 
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
#>         lvp2 <- 8
#>         label("Second peripheral volume of distribution (Vp2)")
#>         lq2 <- 0.5
#>         label("Second intercompartmental clearance (Q2)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>     })
#>     model({
#>         ka <- exp(lka)
#>         vp2 <- exp(lvp2)
#>         q2 <- exp(lq2)
#>         k13 <- q2/vc
#>         k31 <- q2/vp2
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
