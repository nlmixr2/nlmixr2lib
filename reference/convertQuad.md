# Convert linear elimination to quadratic elimination

Convert linear elimination to quadratic elimination

## Usage

``` r
convertQuad(ui, ek = c("Ik", "Ek"), cc = c("Ce", "Cc"), ek2 = "Ek2")
```

## Arguments

- ui:

  rxode2 model

- ek:

  simulation linear constant

- cc:

  the concentration value

- ek2:

  quadratic coefficient

## Value

model with linear effect converted to quadratic effect

## See also

Other PD:
[`addBaseline1exp()`](https://nlmixr2.github.io/nlmixr2lib/reference/addBaseline1exp.md),
[`addBaselineConst()`](https://nlmixr2.github.io/nlmixr2lib/reference/addBaselineConst.md),
[`addBaselineExp()`](https://nlmixr2.github.io/nlmixr2lib/reference/addBaselineExp.md),
[`addBaselineLin()`](https://nlmixr2.github.io/nlmixr2lib/reference/addBaselineLin.md),
[`addDirectLin()`](https://nlmixr2.github.io/nlmixr2lib/reference/addDirectLin.md),
[`convertEmax()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertEmax.md),
[`convertLogLin()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertLogLin.md)

## Author

Matthew L. Fidler

## Examples

``` r
readModelDb("PK_2cmt_no_depot") |>
  addIndirectLin(stim="out") |>
  convertQuad()
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>      lcl      lvc      lvp       lq   propSd     lkin    lkout      lEk 
#>      1.0      3.0      5.0      0.1      0.5      0.1      0.1      0.1 
#> effectSd     uEk2 
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
#>         lEk <- 0.1
#>         label("linear effect constant (Ek)")
#>         effectSd <- c(0, 0.1)
#>         label("additive error for effect")
#>         uEk2 <- 0.1
#>         label("untransformed quadratic slope (Ek2)")
#>     })
#>     model({
#>         Ek2 <- uEk2
#>         kin <- exp(lkin)
#>         kout <- exp(lkout)
#>         Ek <- exp(lEk)
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
#>         d/dt(R) <- kin - kout * R * (1 + (Ek * Cc + Ek2 * Cc^2))
#>         effect <- R
#>         effect ~ add(effectSd)
#>     })
#> }

readModelDb("PK_2cmt_no_depot") |>
  addDirectLin() |>
  convertQuad()
#>  
#>  
#>  ── rxode2-based free-form 2-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>      lcl      lvc      lvp       lq   propSd      uEk effectSd     uEk2 
#>      1.0      3.0      5.0      0.1      0.5      0.1      0.1      0.1 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1          central
#> 2                  2      peripheral1
#>  ── Multiple Endpoint Model ($multipleEndpoint): ──  
#>     variable                   cmt                   dvid*
#> 1     Cc ~ …     cmt='Cc' or cmt=3     dvid='Cc' or dvid=1
#> 2 effect ~ … cmt='effect' or cmt=4 dvid='effect' or dvid=2
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
#>         uEk <- 0.1
#>         label("untransformed slope (Ek)")
#>         effectSd <- c(0, 0.1)
#>         label("additive error for effect")
#>         uEk2 <- 0.1
#>         label("untransformed quadratic slope (Ek2)")
#>     })
#>     model({
#>         Ek2 <- uEk2
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
#>         effect <- Ek * Cc + Ek2 * Cc^2
#>         effect ~ add(effectSd)
#>     })
#> }

readModelDb("PK_2cmt_no_depot") |>
 addEffectCmtLin() |>
 convertQuad()
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>      lcl      lvc      lvp       lq   propSd     lke0      uEk effectSd 
#>      1.0      3.0      5.0      0.1      0.5      0.1      0.1      0.1 
#>     uEk2 
#>      0.1 
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
#>         uEk2 <- 0.1
#>         label("untransformed quadratic slope (Ek2)")
#>     })
#>     model({
#>         Ek2 <- uEk2
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
#>         effect <- Ek * Ce + Ek2 * Ce^2
#>         effect ~ add(effectSd)
#>     })
#> }
```
