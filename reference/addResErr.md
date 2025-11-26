# Add residual error to a model

Add residual error to a model

## Usage

``` r
addResErr(ui, reserr, endpoint)
```

## Arguments

- ui:

  The model as a function

- reserr:

  The type or types of residual error (currently `"addSd"`, `"propSd"`,
  and `"lnormSd"` are accepted)

- endpoint:

  the endpoint to apply the error; will default to the first error in
  the model

## Value

The model with residual error modified

## Details

For `reserr`, the parameter will be named with the dependent variable
from the model as a prefix. For example, if the dependent variable in
the model is `Cc`, the parameter name for `propSd` will become
`CcpropSd`.

## Examples

``` r
library(rxode2)
readModelDb("PK_1cmt") |> addResErr("addSd")
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
readModelDb("PK_1cmt") |> addResErr("lnormSd")
#>  
#>  
#>  
#>  
#> ! remove population parameter `propSd`
#> ℹ add residual parameter `CcLnormSd` and set estimate to 1
#> ℹ change initial estimate of `CcLnormSd` to `0.5`
#>  ── rxode2-based solved PK 1-compartment model ────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>       lka       lcl       lvc CcLnormSd 
#>      0.45      1.00      3.45      0.50 
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
#>         CcLnormSd <- c(0, 0.5)
#>     })
#>     model({
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         Cc <- linCmt()
#>         Cc ~ lnorm(CcLnormSd)
#>     })
#> }
readModelDb("PK_1cmt") |> addResErr(c("addSd", "propSd"))
#>  
#>  
#>  
#>  
#> ! remove population parameter `propSd`
#> ℹ add residual parameter `CcAddSd` and set estimate to 1
#> ℹ add residual parameter `CcPropSd` and set estimate to 1
#> ℹ change initial estimate of `CcAddSd` to `1`
#> ℹ change initial estimate of `CcPropSd` to `0.5`
#>  ── rxode2-based solved PK 1-compartment model ────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>      lka      lcl      lvc  CcAddSd CcPropSd 
#>     0.45     1.00     3.45     1.00     0.50 
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
#>         CcPropSd <- c(0, 0.5)
#>     })
#>     model({
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         Cc <- linCmt()
#>         Cc ~ add(CcAddSd) + prop(CcPropSd)
#>     })
#> }
```
