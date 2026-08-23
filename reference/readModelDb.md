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
#>   # Issue #482: what each ODE state holds, in what amount units, in what
#>   # biological matrix. analyte/specimen proposed by a local model from the
#>   # model description; units derived from the units block. verified = FALSE
#>   # means NOT checked against the source paper.
#>   compartmentData <- list(
#>     depot   = list(analyte = "drug", units = NA_character_, specimen = "administration site", verified = FALSE),
#>     central = list(analyte = "drug", units = NA_character_, specimen = "plasma", verified = FALSE)
#>   )
#> 
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
#> <environment: 0x55be221d49a0>
```
