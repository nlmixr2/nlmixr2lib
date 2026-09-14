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
#> <environment: 0x558a8a64e030>
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
#> Error in rbind(deparse.level, ...): numbers of columns of arguments do not match
modellib(name = "PK_1cmt", reserr = "addSd")
#>  
#>  
#>  
#>  
#> ! remove population parameter `propSd`
#> ℹ add residual parameter `CcAddSd` and set estimate to 1
#> Error in rbind(deparse.level, ...): numbers of columns of arguments do not match
```
