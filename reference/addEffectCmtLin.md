# Add effect compartment

Add effect compartment

## Usage

``` r
addEffectCmtLin(
  ui,
  ke0 = "ke0",
  cc = "Cc",
  ce = "Ce",
  ek = "Ek",
  effect = "effect"
)
```

## Arguments

- ui:

  rxode2 model

- ke0:

  This is the effect compartment keo rate

- cc:

  the concentration value

- ce:

  This is the concentration in the effect compartment

- ek:

  simulation linear constant

- effect:

  the effect variable that will be modeled

## Value

a model with an effect compartment attached

## Author

Matthew L. Fidler

## Examples

``` r
readModelDb("PK_2cmt_no_depot") |>
   addEffectCmtLin()
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>      lcl      lvc      lvp       lq   propSd     lke0      uEk effectSd 
#>      1.0      3.0      5.0      0.1      0.5      0.1      0.1      0.1 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1          central
#> 2                  2      peripheral1
#> 3                  3               Ce
#>  ── Multiple Endpoint Model ($multipleEndpoint): ──  
#>     variable                   cmt                   dvid*
#> 1     Cc ~ …     cmt='Cc' or cmt=4     dvid='Cc' or dvid=1
#> 2 effect ~ … cmt='effect' or cmt=5 dvid='effect' or dvid=2
#>   * If dvids are outside this range, all dvids are re-numered sequentially, ie 1,7, 10 becomes 1,2,3 etc
#> 
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
#>         lke0 <- 0.1
#>         label("effect compartment rate (ke0)")
#>         uEk <- 0.1
#>         label("untransformed linear slope (Ek)")
#>         effectSd <- c(0, 0.1)
#>         label("additive error for effect")
#>     })
#>     model({
#>         ke0 <- exp(lke0)
#>         Ek <- uEk
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         vp <- exp(lvp)
#>         q <- exp(lq)
#>         kel <- cl/vc
#>         k12 <- q/vc
#>         k21 <- q/vp
#>         d/dt(central) <- -kel * central - k12 * central + k21 * 
#>             peripheral1
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>         d/dt(Ce) <- ke0 * (Cc - Ce)
#>         effect <- Ce * Ek
#>         effect ~ add(effectSd)
#>     })
#> }

# Can also be changed to the more typical Emax with constant (estimated) baselie
readModelDb("PK_2cmt_no_depot") |>
  addEffectCmtLin() |>
  convertEmaxHill() |>
  addBaselineConst()
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>       lcl       lvc       lvp        lq    propSd      lke0  effectSd     lEmax 
#>  1.000000  3.000000  5.000000  0.100000  0.500000  0.100000  0.100000  0.100000 
#>     lEC50       lgg       uEb 
#>  0.100000 -2.302585  0.100000 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1          central
#> 2                  2      peripheral1
#> 3                  3               Ce
#>  ── Multiple Endpoint Model ($multipleEndpoint): ──  
#>     variable                   cmt                   dvid*
#> 1     Cc ~ …     cmt='Cc' or cmt=4     dvid='Cc' or dvid=1
#> 2 effect ~ … cmt='effect' or cmt=5 dvid='effect' or dvid=2
#>   * If dvids are outside this range, all dvids are re-numered sequentially, ie 1,7, 10 becomes 1,2,3 etc
#> 
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
#>         lke0 <- 0.1
#>         label("effect compartment rate (ke0)")
#>         effectSd <- c(0, 0.1)
#>         label("additive error for effect")
#>         lEmax <- 0.1
#>         label("Maximum effect (Emax)")
#>         lEC50 <- 0.1
#>         label("Concentration of 50% Emax (Emax)")
#>         lgg <- -2.30258509299405
#>         label("logit-constrained Hill coefficient g")
#>         uEb <- 0.1
#>         label("untransformed constant baseline (Eb)")
#>     })
#>     model({
#>         Eb <- uEb
#>         Emax <- exp(lEmax)
#>         EC50 <- exp(lEC50)
#>         g <- expit(lgg, 0.1, 10)
#>         ke0 <- exp(lke0)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         vp <- exp(lvp)
#>         q <- exp(lq)
#>         kel <- cl/vc
#>         k12 <- q/vc
#>         k21 <- q/vp
#>         d/dt(central) <- -kel * central - k12 * central + k21 * 
#>             peripheral1
#>         d/dt(peripheral1) <- k12 * central - k21 * peripheral1
#>         Cc <- central/vc
#>         Cc ~ prop(propSd)
#>         d/dt(Ce) <- ke0 * (Cc - Ce)
#>         effect <- Ce * Ek + Eb
#>         effect ~ add(effectSd)
#>     })
#> }
```
