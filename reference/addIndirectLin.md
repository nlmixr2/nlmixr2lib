# Add linear indirect response model

Add linear indirect response model

## Usage

``` r
addIndirectLin(
  ui,
  stim = c("in", "out"),
  inhib = c("in", "out"),
  ek = "Ek",
  ik = "Ik",
  kin = "kin",
  kout = "kout",
  cc = "Cc",
  R = "R",
  effect = "effect"
)
```

## Arguments

- ui:

  rxode2 model

- stim:

  what type of stimulation indirect response model:

  \- \`in\`: stimulation of input

  \- \`out\`: stimulation of output

- inhib:

  what type of inhibition indirect response model:

  \- \`in\`: inhibition of input

  \- \`out\`: inhibition of output

- ek:

  simulation linear constant

- ik:

  inhibition linear constant

- kin:

  this is the kin parameter name

- kout:

  this is the kout parameter name

- cc:

  the concentration value

- R:

  drug response compartment

- effect:

  the effect variable that will be modeled

## Value

model with linear indirect effect added

Note that while linear indirect effects are not common, it allows an
easier hook to produce other standard effect curves like Emax/Imax,
Hill, etc.

## See also

Other Indirect response:
[`addIndirect()`](https://nlmixr2.github.io/nlmixr2lib/reference/addIndirect.md),
[`convertKinR0()`](https://nlmixr2.github.io/nlmixr2lib/reference/convertKinR0.md)

## Author

Matthew L. Fidler

## Examples

``` r
readModelDb("PK_2cmt_no_depot") |> addIndirectLin(stim="in")
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
#>         d/dt(R) <- kin * (1 + Ek * Cc) - kout * R
#>         effect <- R
#>         effect ~ add(effectSd)
#>     })
#> }

readModelDb("PK_2cmt_no_depot") |> addIndirectLin(stim="out")
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
#>         d/dt(R) <- kin - kout * R * (1 + Ek * Cc)
#>         effect <- R
#>         effect ~ add(effectSd)
#>     })
#> }

readModelDb("PK_2cmt_no_depot") |> addIndirectLin(inhib="in")
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>      lcl      lvc      lvp       lq   propSd     lkin    lkout      lIk 
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
#>         lIk <- c(-Inf, 0.1, 1)
#>         label("linear inhibition constant (Ik)")
#>         effectSd <- c(0, 0.1)
#>         label("additive error for effect")
#>     })
#>     model({
#>         kin <- exp(lkin)
#>         kout <- exp(lkout)
#>         Ik <- exp(lIk)
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
#>         d/dt(R) <- kin * (1 - Ik * Cc) - kout * R
#>         effect <- R
#>         effect ~ add(effectSd)
#>     })
#> }

readModelDb("PK_2cmt_no_depot") |> addIndirectLin(inhib="out")
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>      lcl      lvc      lvp       lq   propSd     lkin    lkout      lIk 
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
#>         lIk <- c(-Inf, 0.1, 1)
#>         label("linear inhibition constant (Ik)")
#>         effectSd <- c(0, 0.1)
#>         label("additive error for effect")
#>     })
#>     model({
#>         kin <- exp(lkin)
#>         kout <- exp(lkout)
#>         Ik <- exp(lIk)
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
#>         d/dt(R) <- kin - kout * R * (1 - Ik * Cc)
#>         effect <- R
#>         effect ~ add(effectSd)
#>     })
#> }
```
