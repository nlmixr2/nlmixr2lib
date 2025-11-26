# Search within a model to replace part of the model

Search within a model to replace part of the model

## Usage

``` r
searchReplace(object, find, replace)

searchReplaceHelper(object, find, replace)
```

## Arguments

- object:

  function specifying the nlmixr2 model

- find, replace:

  Character scalars of parts of the model to replace

## Value

`object` with `find` replaced with `replace`

## Functions

- `searchReplaceHelper()`: A helper function for searchReplace (not
  intended for users to use directly)
