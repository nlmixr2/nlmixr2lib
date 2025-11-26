# Combine two strings using a naming convention

Combine two in a manner similar to \`paste()\` strings using the default
combine type

## Usage

``` r
combinePaste2(
  a,
  b,
  combineType = c("default", "snake", "camel", "dot", "blank")
)
```

## Arguments

- a:

  first string to combine

- b:

  second string to combine

- combineType:

  is the type of combination; can be:

  \- `"default"`: default combine (set with \`defaultCombine()\`)

  \- `"camel"`: camelCase combine

  \- `"snake"`: snake_case combine

  \- `"dot"`: dot combine (i.e. "a.b")

  \- `"blank"`: no separator (i.e. "ab")

## Value

Combined strings separated with \`defaultCombine()\`

## Author

Matthew L. Fidler

## Examples

``` r
combinePaste2("f", "depot")
#> [1] "fDepot"

combinePaste2("f", "depot", "snake")
#> [1] "f_depot"

combinePaste2("f", "depot", "dot")
#> [1] "f.depot"

combinePaste2("f", "depot", "blank")
#> [1] "fdepot"
```
