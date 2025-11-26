# Convert a kin/kout indirect response model to R0 and kout

This replaces the kin/kout parameterization to the R0 and kout
parameterization

## Usage

``` r
convertKinR0(ui, kin = "kin", kout = "kout", R = "R", R0 = "R0")
```

## Arguments

- ui:

  a rxode2 user function

- kin:

  the kin variable (by default is "kin")

- kout:

  the kout variable (by default is "kout")

- R:

  the compartment variable (by default is "R")

- R0:

  the R0 variable (by default is "R0")

## Value

a model where the estimated kin is changed to the estimated R0

## See also

Other Indirect response:
[`addIndirect()`](https://nlmixr2.github.io/nlmixr2lib/reference/addIndirect.md),
[`addIndirectLin()`](https://nlmixr2.github.io/nlmixr2lib/reference/addIndirectLin.md)

## Author

Matthew L. Fidler

## Examples

``` r
addIndirect(stim="in") |> convertKinR0()
#>  
#>  
#>  ── rxode2-based free-form 1-cmt ODE model ────────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lkout effectSd    lEmax    lEC50      uR0 
#>      0.1      0.1      0.1      0.1      0.1 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name
#> 1                  1                R
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     ini({
#>         lkout <- 0.1
#>         label("first order rate response loss (kout)")
#>         effectSd <- c(0, 0.1)
#>         label("additive error for effect")
#>         lEmax <- 0.1
#>         label("Maximum effect (Emax)")
#>         lEC50 <- 0.1
#>         label("Concentration of 50% Emax (Emax)")
#>         uR0 <- 0.1
#>         label("untransformed baseline (R0)")
#>     })
#>     model({
#>         R0 <- uR0
#>         Emax <- exp(lEmax)
#>         EC50 <- exp(lEC50)
#>         kout <- exp(lkout)
#>         R(0) <- R0
#>         d/dt(R) <- kout * R0 * (1 + Emax * Cc/(Cc + EC50)) - 
#>             kout * R
#>         effect <- R
#>         effect ~ add(effectSd)
#>     })
#> }
```
