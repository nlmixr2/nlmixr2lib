# Add an estimated baseline linear constant

Add an estimated baseline linear constant

## Usage

``` r
addBaselineLin(ui, effect = "effect", eb = "Eb", time = "time")
```

## Arguments

- ui:

  rxode2 model

- effect:

  the effect variable that will be modeled

- eb:

  baseline constant parameter

- time:

  the time or other variable used for baseline decay

## Value

model with baseline linear constant

## See also

Other PD:
[`addBaseline1exp()`](https://nlmixr2.github.io/nlmixr2lib/reference/addBaseline1exp.md),
[`addBaselineConst()`](https://nlmixr2.github.io/nlmixr2lib/reference/addBaselineConst.md),
[`addBaselineExp()`](https://nlmixr2.github.io/nlmixr2lib/reference/addBaselineExp.md),
[`addDirectLin()`](https://nlmixr2.github.io/nlmixr2lib/reference/addDirectLin.md),
[`convertEmax()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertEmax.md),
[`convertLogLin()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertLogLin.md),
[`convertQuad()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertQuad.md)

## Author

Matthew L. Fidler

## Examples

``` r
readModelDb("PK_2cmt_no_depot") |>
  addDirectLin() |>
  convertQuad() |>
  addBaselineLin()
#>  
#>  
#>  ── rxode2-based free-form 2-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>      lcl      lvc      lvp       lq   propSd      uEk effectSd     uEk2 
#>      1.0      3.0      5.0      0.1      0.5      0.1      0.1      0.1 
#>      uEb 
#>      0.1 
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
#>         uEb <- 0.1
#>         label("untransformed constant baseline (Eb)")
#>     })
#>     model({
#>         Eb <- uEb
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
#>         effect <- Ek * Cc + Ek2 * Cc^2 + Eb * time
#>         effect ~ add(effectSd)
#>     })
#> }
```
