# Model library for nlmixr2

This is a data frame of the available models in nlmixr2lib, it is
generated with the package. Custom modeldb may be used.

## Usage

``` r
modeldb
```

## Format

A data frame with 57 rows and 9 columns

- name:

  Model name that can be used to extract the model from the model
  library

- description:

  Model description in free from text; in model itself

- parameters:

  A comma separated string listing either the parameter in the model
  defined by population/individual effects or a population effect
  parameter

- DV:

  The definition of the dependent variable(s)

- linCmt:

  Logical flag indicating if solved models are used (TRUE) or not
  (FALSE)

- algebraic:

  Logical flag indicating if the model is purely algebraic: TRUE no
  linCmt() and no ODEs; FALSE otherwise

- dosing:

  A comma separated string of identified dosing compartments

- depends:

  A comma separated string of objects the model depends on

- filename:

  Filename of the model. By default these are installed in the model
  library and read on demand
