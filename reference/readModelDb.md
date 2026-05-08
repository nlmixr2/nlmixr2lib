# Read a model from the nlmixr2 model database

Read a model from the nlmixr2 model database

## Usage

``` r
readModelDb(name)
```

## Arguments

- name:

  The name of the model (must be one of `modeldb$name`)

## Value

The model as a function

## Examples

``` r
readModelDb("PK_1cmt")
#> function() {
#>   description <- "One compartment PK model with linear clearance"
#>   reference <- "nlmixr2lib template"
#>   units <- list(time = "time_unit", dosing = "dose_unit", concentration = "conc_unit/vol_unit")
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
#> <environment: 0x55931bb85c90>
```
