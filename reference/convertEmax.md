# Convert linear effect to Emax effect

Convert linear effect to Emax effect

## Usage

``` r
convertEmax(
  ui,
  emax = "Emax",
  ec50 = "EC50",
  imax = "Imax",
  ic50 = "IC50",
  ek = c("Ik", "Ek"),
  cc = c("Ec", "Cc")
)
```

## Arguments

- ui:

  rxode2 model

- emax:

  Emax parameter

- ec50:

  EC50 parameter

- imax:

  Imax parameter used when input model contains "Ik" instead of "Ek"

- ic50:

  IC50 parameter used when input model contains "Ik" instead of "Ek"

- ek:

  simulation linear constant

- cc:

  the concentration value

## Value

Model with the linear effect converted to an Emax effect

## See also

Other PD:
[`addBaseline1exp()`](https://nlmixr2.github.io/nlmixr2lib/reference/addBaseline1exp.md),
[`addBaselineConst()`](https://nlmixr2.github.io/nlmixr2lib/reference/addBaselineConst.md),
[`addBaselineExp()`](https://nlmixr2.github.io/nlmixr2lib/reference/addBaselineExp.md),
[`addBaselineLin()`](https://nlmixr2.github.io/nlmixr2lib/reference/addBaselineLin.md),
[`addDirectLin()`](https://nlmixr2.github.io/nlmixr2lib/reference/addDirectLin.md),
[`convertLogLin()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertLogLin.md),
[`convertQuad()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertQuad.md)

## Author

Matthew L. Fidler

## Examples

``` r
readModelDb("PK_2cmt_no_depot") |>
  addIndirectLin(stim="in") |>
  convertEmax()
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>      lcl      lvc      lvp       lq   propSd     lkin    lkout effectSd 
#>      1.0      3.0      5.0      0.1      0.5      0.1      0.1      0.1 
#>    lEmax    lEC50 
#>      0.1      0.1 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1          central
#> 2                  2      peripheral1
#> 3                  3                R
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
#>         lkin <- 0.1
#>         label("zero order response production(kin)")
#>         lkout <- 0.1
#>         label("first order rate response loss (kout)")
#>         effectSd <- c(0, 0.1)
#>         label("additive error for effect")
#>         lEmax <- 0.1
#>         label("Maximum effect (Emax)")
#>         lEC50 <- 0.1
#>         label("Concentration of 50% Emax (Emax)")
#>     })
#>     model({
#>         Emax <- exp(lEmax)
#>         EC50 <- exp(lEC50)
#>         kin <- exp(lkin)
#>         kout <- exp(lkout)
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
#>         R(0) <- kin/kout
#>         d/dt(R) <- kin * (1 + Emax * Cc/(Cc + EC50)) - kout * 
#>             R
#>         effect <- R
#>         effect ~ add(effectSd)
#>     })
#> }

# When emax=1
readModelDb("PK_2cmt_no_depot") |>
  addIndirectLin(stim="in") |>
  convertEmax(emax=1)
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>      lcl      lvc      lvp       lq   propSd     lkin    lkout effectSd 
#>      1.0      3.0      5.0      0.1      0.5      0.1      0.1      0.1 
#>    lEC50 
#>      0.1 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1          central
#> 2                  2      peripheral1
#> 3                  3                R
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
#>         lkin <- 0.1
#>         label("zero order response production(kin)")
#>         lkout <- 0.1
#>         label("first order rate response loss (kout)")
#>         effectSd <- c(0, 0.1)
#>         label("additive error for effect")
#>         lEC50 <- 0.1
#>         label("Concentration of 50% 1 (1)")
#>     })
#>     model({
#>         EC50 <- exp(lEC50)
#>         kin <- exp(lkin)
#>         kout <- exp(lkout)
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
#>         R(0) <- kin/kout
#>         d/dt(R) <- kin * (1 + Cc/(Cc + EC50)) - kout * R
#>         effect <- R
#>         effect ~ add(effectSd)
#>     })
#> }
```
