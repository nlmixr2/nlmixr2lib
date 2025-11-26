# Converts a linear effect to a log-linear effect

Converts a linear effect to a log-linear effect

## Usage

``` r
convertLogLin(ui, ek = c("Ik", "Ek"), cc = c("Ce", "Cc"))
```

## Arguments

- ui:

  rxode2 model

- ek:

  simulation linear constant

- cc:

  the concentration value

## Value

model converted from linear to log-linear effect

## See also

Other PD:
[`addBaseline1exp()`](https://nlmixr2.github.io/nlmixr2lib/reference/addBaseline1exp.md),
[`addBaselineConst()`](https://nlmixr2.github.io/nlmixr2lib/reference/addBaselineConst.md),
[`addBaselineExp()`](https://nlmixr2.github.io/nlmixr2lib/reference/addBaselineExp.md),
[`addBaselineLin()`](https://nlmixr2.github.io/nlmixr2lib/reference/addBaselineLin.md),
[`addDirectLin()`](https://nlmixr2.github.io/nlmixr2lib/reference/addDirectLin.md),
[`convertEmax()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertEmax.md),
[`convertQuad()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertQuad.md)

## Author

Matthew L. Fidler

## Examples

``` r
readModelDb("PK_2cmt_no_depot") |>
  addDirectLin() |>
  convertLogLin()
#>  
#>  
#>  ── rxode2-based free-form 2-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>      lcl      lvc      lvp       lq   propSd      uEk effectSd 
#>      1.0      3.0      5.0      0.1      0.5      0.1      0.1 
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
#>     })
#>     model({
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
#>         effect <- Ek * log(Cc)
#>         effect ~ add(effectSd)
#>     })
#> }

readModelDb("PK_2cmt_no_depot") |>
  addIndirectLin(stim="out") |>
  convertLogLin()
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>      lcl      lvc      lvp       lq   propSd     lkin    lkout      lEk 
#>      1.0      3.0      5.0      0.1      0.5      0.1      0.1      0.1 
#> effectSd 
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
#>         lEk <- 0.1
#>         label("linear effect constant (Ek)")
#>         effectSd <- c(0, 0.1)
#>         label("additive error for effect")
#>     })
#>     model({
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
#>         d/dt(R) <- kin - kout * R * (1 + Ek * log(Cc))
#>         effect <- R
#>         effect ~ add(effectSd)
#>     })
#> }
```
