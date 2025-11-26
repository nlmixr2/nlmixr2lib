# Change the default combine type for the package

Change the default combine type for the package

## Usage

``` r
setCombineType(combineType = c("snake", "camel", "dot", "blank"))
```

## Arguments

- combineType:

  this is the default combine type: - `"default"`: default combine -
  `"snake"`: snake_case combine - `"camel"`: camelCase combine -
  `"dot"`: dot combine (i.e. "a.b") - `"blank"`: no separator (i.e.
  "ab")

## Author

Matthew L. Fidler

## Examples

``` r
# Change to the popular snake_case
setCombineType("snake")
defaultCombine("a", "b")
#> [1] "a_b"

# Change back to nlmixr2/rxode2 default camelCase

setCombineType("camel")
defaultCombine("a", "b")
#> [1] "aB"

# This is used to change the naming convention for parameters
# produced by this package
```
