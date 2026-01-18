# Get the model from the model library

This function gets a model from the available model library

## Usage

``` r
modellib(name = NULL, eta = NULL, reserr = NULL)
```

## Arguments

- name:

  character with the name of the model to load (if `NULL`, lists all
  available base models)

- eta:

  vector with the parameters to add random effects (sometimes referred
  to as inter-individual variability, IIV) on

- reserr:

  The type or types of residual error (currently `"addSd"`, `"propSd"`,
  and `"lnormSd"` are accepted)

## Value

The function returns a function the model code (or `NULL` if the
`model = NULL`)

## Details

This is a very first draft just to look at the proof of concept

## Examples

``` r
modellib(name = "PK_1cmt")
#> function() {
#>   description <- "One compartment PK model with linear clearance"
#>   ini({
#>     lka <- 0.45 ; label("Absorption rate (Ka)")
#>     lcl <- 1 ; label("Clearance (CL)")
#>     lvc  <- 3.45 ; label("Central volume of distribution (V)")
#>     propSd <- 0.5 ; label("Proportional residual error (fraction)")
#>   })
#>   model({
#>     ka <- exp(lka)
#>     cl <- exp(lcl)
#>     vc  <- exp(lvc)
#> 
#>     Cc <- linCmt()
#>     Cc ~ prop(propSd)
#>   })
#> }
#> <environment: 0x55b9fb19d4c0>
modellib(name = "PK_1cmt", eta = c("ka", "vc"), reserr = "addSd")
#>  
#>  
#>  
#>  
#> → Adding eta to lka instead of ka due to mu-referencing
#>  
#>  
#> → Adding eta to lvc instead of vc due to mu-referencing
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> ℹ promote `etaKa` to between subject variability with initial estimate 0.1
#> ℹ change initial estimate of `etaKa` to `0.1`
#> ℹ promote `etaVc` to between subject variability with initial estimate 0.1
#> ℹ change initial estimate of `etaVc` to `0.1`
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> ! remove population parameter `propSd`
#> ℹ add residual parameter `CcAddSd` and set estimate to 1
#> ℹ change initial estimate of `CcAddSd` to `1`
#>  ── rxode2-based solved PK 1-compartment model ────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>     lka     lcl     lvc CcAddSd 
#>    0.45    1.00    3.45    1.00 
#> 
#> Omega ($omega): 
#>       etaKa etaVc
#> etaKa   0.1   0.0
#> etaVc   0.0   0.1
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name Rate   Off Internal #
#> 1                  1            depot TRUE FALSE          1
#> 2                  2          central TRUE FALSE          2
#>  ── μ-referencing ($muRefTable): ──  
#>   theta   eta level
#> 1   lka etaKa    id
#> 2   lvc etaVc    id
#> 
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     description <- "One compartment PK model with linear clearance"
#>     ini({
#>         lka <- 0.45
#>         label("Absorption rate (Ka)")
#>         lcl <- 1
#>         label("Clearance (CL)")
#>         lvc <- 3.45
#>         label("Central volume of distribution (V)")
#>         CcAddSd <- c(0, 1)
#>         etaKa ~ 0.1
#>         etaVc ~ 0.1
#>     })
#>     model({
#>         ka <- exp(lka + etaKa)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc + etaVc)
#>         Cc <- linCmt()
#>         Cc ~ add(CcAddSd)
#>     })
#> }
modellib(name = "PK_1cmt", reserr = "addSd")
#>  
#>  
#>  
#>  
#> ! remove population parameter `propSd`
#> ℹ add residual parameter `CcAddSd` and set estimate to 1
#> ℹ change initial estimate of `CcAddSd` to `1`
#>  ── rxode2-based solved PK 1-compartment model ────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>     lka     lcl     lvc CcAddSd 
#>    0.45    1.00    3.45    1.00 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name Rate   Off Internal #
#> 1                  1            depot TRUE FALSE          1
#> 2                  2          central TRUE FALSE          2
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     description <- "One compartment PK model with linear clearance"
#>     ini({
#>         lka <- 0.45
#>         label("Absorption rate (Ka)")
#>         lcl <- 1
#>         label("Clearance (CL)")
#>         lvc <- 3.45
#>         label("Central volume of distribution (V)")
#>         CcAddSd <- c(0, 1)
#>     })
#>     model({
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         Cc <- linCmt()
#>         Cc ~ add(CcAddSd)
#>     })
#> }
```
