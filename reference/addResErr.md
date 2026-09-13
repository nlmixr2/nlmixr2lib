# Add residual error to a model

Add residual error to a model

## Usage

``` r
addResErr(ui, reserr, endpoint, model)
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

- model:

  Deprecated alias for `ui`. Supplying `model` instead of `ui` still
  works but emits a deprecation warning.

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
#> Error in rbind(deparse.level, ...): numbers of columns of arguments do not match
readModelDb("PK_1cmt") |> addResErr("lnormSd")
#>  
#>  
#>  
#>  
#> ! remove population parameter `propSd`
#> ℹ add residual parameter `CcLnormSd` and set estimate to 1
#> Error in rbind(deparse.level, ...): numbers of columns of arguments do not match
readModelDb("PK_1cmt") |> addResErr(c("addSd", "propSd"))
#>  
#>  
#>  
#>  
#> ! remove population parameter `propSd`
#> ℹ add residual parameter `CcAddSd` and set estimate to 1
#> Error in rbind(deparse.level, ...): numbers of columns of arguments do not match
```
