# Convert linear effect to Emax-Hill effect

Convert linear effect to Emax-Hill effect

## Usage

``` r
convertEmaxHill(
  ui,
  emax = "Emax",
  ec50 = "EC50",
  g = "g",
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

- g:

  hill coefficient

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

## Author

Matthew L. Fidler

## Examples

``` r
readModelDb("PK_2cmt_no_depot") |>
  addIndirectLin(stim="in") |>
  convertEmaxHill()
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>       lcl       lvc       lvp        lq    propSd      lkin     lkout  effectSd 
#>  1.000000  3.000000  5.000000  0.100000  0.500000  0.100000  0.100000  0.100000 
#>     lEmax     lEC50       lgg 
#>  0.100000  0.100000 -2.302585 
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
#>         lgg <- -2.30258509299405
#>         label("logit-constrained Hill coefficient g")
#>     })
#>     model({
#>         Emax <- exp(lEmax)
#>         EC50 <- exp(lEC50)
#>         g <- expit(lgg, 0.1, 10)
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
#>         d/dt(R) <- kin * (1 + Emax * Cc^g/(Cc^g + EC50^g)) - 
#>             kout * R
#>         effect <- R
#>         effect ~ add(effectSd)
#>     })
#> }

# can also specify as emax=1

readModelDb("PK_2cmt_no_depot") |>
  addIndirectLin(stim="in") |>
  convertEmaxHill(emax=1)
#>  
#>  
#>  ── rxode2-based free-form 3-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>       lcl       lvc       lvp        lq    propSd      lkin     lkout  effectSd 
#>  1.000000  3.000000  5.000000  0.100000  0.500000  0.100000  0.100000  0.100000 
#>     lEC50       lgg 
#>  0.100000 -2.302585 
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
#>         lgg <- -2.30258509299405
#>         label("logit-constrained Hill coefficient g")
#>     })
#>     model({
#>         EC50 <- exp(lEC50)
#>         g <- expit(lgg, 0.1, 10)
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
#>         d/dt(R) <- kin * (1 + Cc^g/(Cc^g + EC50^g)) - kout * 
#>             R
#>         effect <- R
#>         effect ~ add(effectSd)
#>     })
#> }
```
