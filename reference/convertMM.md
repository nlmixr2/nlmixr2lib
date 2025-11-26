# Convert models from linear elimination to Michaelis-Menten elimination

Convert models from linear elimination to Michaelis-Menten elimination

## Usage

``` r
convertMM(
  ui,
  central = "central",
  elimination = "kel",
  vm = "vm",
  km = "km",
  vc = "vc"
)
```

## Arguments

- ui:

  model to convert

- central:

  the central compartment where the elimination is present

- elimination:

  variable for the elimination constant in the model

- vm:

  variable name for Vmax in the model

- km:

  variable name for Km in the model

- vc:

  variable name for Vc in the model

## Value

new model changing linear elimination to Michaelis-Menten elimination

## Author

Matthew L. Fidler

## Examples

``` r
readModelDb("PK_1cmt_des") |> convertMM()
#>  
#>  
#>  ── rxode2-based free-form 2-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka    lvc propSd    lvm    lkm 
#>   0.45   3.45   0.50   0.10   0.10 
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
#>         lvc <- 3.45
#>         label("Central volume of distribution (V)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lvm <- 0.1
#>         lkm <- 0.1
#>     })
#>     model({
#>         vm <- exp(lvm)
#>         km <- exp(lkm)
#>         ka <- exp(lka)
#>         vc <- exp(lvc)
#>         d/dt(depot) <- -ka * depot
#>         d/dt(central) <- ka * depot - (vm * central/vc)/(km + 
#>             central/vc)
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

readModelDb("PK_2cmt_des") |> convertMM()
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka    lvc    lvp     lq propSd    lvm    lkm 
#>   0.45   3.00   5.00   0.10   0.50   0.10   0.10 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1            depot
#> 2                  2          central
#> 3                  3      peripheral1
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     ini({
#>         lka <- 0.45
#>         label("Absorption rate (Ka)")
#>         lvc <- 3
#>         label("Central volume of distribution (V)")
#>         lvp <- 5
#>         label("Peripheral volume of distribution (Vp)")
#>         lq <- 0.1
#>         label("Intercompartmental clearance (Q)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lvm <- 0.1
#>         lkm <- 0.1
#>     })
#>     model({
#>         vm <- exp(lvm)
#>         km <- exp(lkm)
#>         ka <- exp(lka)
#>         vc <- exp(lvc)
#>         vp <- exp(lvp)
#>         q <- exp(lq)
#>         k12 <- q/vc
#>         k21 <- q/vp
#>         d/dt(depot) <- -ka * depot
#>         d/dt(central) <- ka * depot - (vm * central/vc)/(km + 
#>             central/vc) - k12 * central + k21 * peripheral1
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

readModelDb("PK_3cmt_des") |> convertMM()
#>  
#>  
#>  ── rxode2-based free-form 4-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka    lvc    lvp   lvp2     lq    lq2 propSd    lvm    lkm 
#>   0.45   3.00   5.00   8.00   0.10   0.50   0.50   0.10   0.10 
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
#>         lvc <- 3
#>         label("Central volume of distribution (V)")
#>         lvp <- 5
#>         label("Peripheral volume of distribution (Vp)")
#>         lvp2 <- 8
#>         label("Second peripheral volume of distribution (Vp2)")
#>         lq <- 0.1
#>         label("Intercompartmental clearance (Q)")
#>         lq2 <- 0.5
#>         label("Second intercompartmental clearance (Q2)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lvm <- 0.1
#>         lkm <- 0.1
#>     })
#>     model({
#>         vm <- exp(lvm)
#>         km <- exp(lkm)
#>         ka <- exp(lka)
#>         vc <- exp(lvc)
#>         vp <- exp(lvp)
#>         vp2 <- exp(lvp2)
#>         q <- exp(lq)
#>         q2 <- exp(lq2)
#>         k12 <- q/vc
#>         k21 <- q/vp
#>         k13 <- q2/vc
#>         k31 <- q2/vp2
#>         d/dt(depot) <- -ka * depot
#>         d/dt(central) <- ka * depot - (vm * central/vc)/(km + 
#>             central/vc) - k12 * central + k21 * peripheral1 - 
#>             k13 * central + k31 * peripheral2
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         d/dt(peripheral2) <- k13 * central - k31 * peripheral2
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }

readModelDb("PK_3cmt_des") |> removeDepot() |> convertMM()
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lvc    lvp   lvp2     lq    lq2 propSd    lvm    lkm 
#>    3.0    5.0    8.0    0.1    0.5    0.5    0.1    0.1 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1          central
#> 2                  2      peripheral1
#> 3                  3      peripheral2
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     ini({
#>         lvc <- 3
#>         label("Central volume of distribution (V)")
#>         lvp <- 5
#>         label("Peripheral volume of distribution (Vp)")
#>         lvp2 <- 8
#>         label("Second peripheral volume of distribution (Vp2)")
#>         lq <- 0.1
#>         label("Intercompartmental clearance (Q)")
#>         lq2 <- 0.5
#>         label("Second intercompartmental clearance (Q2)")
#>         propSd <- c(0, 0.5)
#>         label("Proportional residual error (fraction)")
#>         lvm <- 0.1
#>         lkm <- 0.1
#>     })
#>     model({
#>         vm <- exp(lvm)
#>         km <- exp(lkm)
#>         vc <- exp(lvc)
#>         vp <- exp(lvp)
#>         vp2 <- exp(lvp2)
#>         q <- exp(lq)
#>         q2 <- exp(lq2)
#>         k12 <- q/vc
#>         k21 <- q/vp
#>         k13 <- q2/vc
#>         k31 <- q2/vp2
#>         d/dt(central) <- -(vm * central/vc)/(km + central/vc) - 
#>             k12 * central + k21 * peripheral1 - k13 * central + 
#>             k31 * peripheral2
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         d/dt(peripheral2) <- k13 * central - k31 * peripheral2
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>     })
#> }
```
