<!-- badges: start -->
[![R-CMD-check](https://github.com/nlmixr2/nlmixr2lib/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nlmixr2/nlmixr2lib/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# nlmixr2lib

This is a very early draft to see if a model library can be set up for `nlmixr2` using a separate R package.
The package currently contains one function `modellib` that loads the base models using the 'modellib.rds', both available in the 'inst' directory.
The 'modellib.rds' file is a saved version of a dataframe including all meta data for the available models.
When running the function without arguments, it will list all available models. When a modelname is given, the function will return a vector with the model code.
The function can be easily adapted to return other types as well, e.g.:

-  to console: `cat(modr,sep="\n)`
-  to file: `writeLines(modr,paste0(model,".r"))`
-  directly evaluated: `eval(parse(text=modr))`
-  to Rstudio's current script: `rstudioapi::insertText(paste(modr,"\n"))`

Within the function, the parameters for IIV or the type of residual error can be provided. In this case the model is extended using the available pipe syntax.
Some examples for usage:

```
modellib(model="PK_1cmt")
modellib(model="PK_1cmt",iiv = c("ka","v"),reserr = "add")
modellib(model="PK_1cmt",reserr = "add")
```

# Assumptions and remarks

Please take into account the following when evaluating the package:

- This package currently only provides a proof of concept
- The dataframe including all meta data contains the bare minimum,
  other variables to include should be considered (e.g. modeltype=
  PK,PD,PKPD,etc.)
  
- Added DV variable in rds dataset, this is currently used for
  residual error (coded as cp~prop(prop.err) or
  linCmt()~prop(prop.err))

- Piping is now used to add IIV and residual error. For IIV the
  parameters are used/checked from the rds dataset. For residual
  error, now a default of proportional error is always assumed and the
  user can choose from pre-defined options: add, prop, add+prop

# Possible extensions

The `modellib` function is set-up in way that it can be easily
extended and used in other applications.  A possible extension could
be implementation in a shiny app. An app can be created to easily add
new models to the model library database (curated?), and directly make
these models available for other users.  I believe there can be added
value in having a base model library that can be easily extended by
the community this way.
