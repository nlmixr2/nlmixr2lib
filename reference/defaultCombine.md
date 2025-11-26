# Default combine strings

Default combine strings

## Usage

``` r
defaultCombine(...)

snakeCombine(...)

camelCombine(...)

dotCombine(...)

blankCombine(...)
```

## Arguments

- ...:

  uses default to combine strings

## Value

combined strings

## Functions

- `snakeCombine()`: use snake_case to combine 2 strings

- `camelCombine()`: use camelCase to combine strings

- `dotCombine()`: use the default method for combining two strings

- `blankCombine()`: combine using a blank separator

## Author

Matthew L. Fidler

## Examples

``` r
# default combine

defaultCombine("f", "depot")
#> [1] "fDepot"

defaultCombine(list(c("a", "funny", "c")))
#> [1] "aFunnyC"

defaultCombine(c("a", "funny", "c"))
#> [1] "aFunnyC"

# snake combine

snakeCombine("f", "depot")
#> [1] "f_depot"

snakeCombine(list(c("a", "funny", "c")))
#> [1] "a_funny_c"

snakeCombine(c("a", "funny", "c"))
#> [1] "a_funny_c"

# dot combine

dotCombine("f", "depot")
#> [1] "f.depot"

dotCombine(list(c("a", "funny", "c")))
#> [1] "a.funny.c"

dotCombine(c("a", "funny", "c"))
#> [1] "a.funny.c"

# blank combine

blankCombine("f", "depot")
#> [1] "fdepot"

blankCombine(list(c("a", "funny", "c")))
#> [1] "afunnyc"

blankCombine(c("a", "funny", "c"))
#> [1] "afunnyc"
```
