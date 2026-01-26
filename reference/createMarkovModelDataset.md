# Create a Markov model dataset for use in estimating Markov models

Create a Markov model dataset for use in estimating Markov models

## Usage

``` r
createMarkovModelDataset(x, ...)

# Default S3 method
createMarkovModelDataset(
  x,
  colCur,
  colPrev = x,
  prefixPrev = "prev",
  prefixCur = "cur",
  ...
)

# S3 method for class 'data.frame'
createMarkovModelDataset(x, colPrev, colCur, ...)
```

## Arguments

- x:

  The object to create a Markov dataset within

- ...:

  Passed to methods

- colPrev, colCur:

  Column names in the dataset for the previous and current states (may
  be any R object that can be coerced to a character string, ordered
  objects are usually preferred)

- prefixPrev, prefixCur:

  Column name prefixes for the previous and current states

## Value

The data.frame modified by adding one-hot columns for the previous and
current states

The data.frame modified by adding one-hot columns for the previous and
current states

## Methods (by class)

- `createMarkovModelDataset(default)`: Create a Markov dataset from a
  pair of vectors

- `createMarkovModelDataset(data.frame)`: Create a Markov dataset from a
  data.frame

## See also

Other Markov models:
[`createMarkovModel()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovModel.md),
[`createMarkovModelFromSingleState()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovModelFromSingleState.md),
[`createMarkovTransitionMatrix()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovTransitionMatrix.md),
[`simMarkov()`](https://nlmixr2.github.io/nlmixr2lib/reference/simMarkov.md),
[`simMarkovId()`](https://nlmixr2.github.io/nlmixr2lib/reference/simMarkovId.md)

## Examples

``` r
createMarkovModelDataset(c(1, 2, 1), colCur = c(1, 2, 3))
#>   prevX1 curX1 prevX2 curX2 prevX3 curX3
#> 1   TRUE  TRUE  FALSE FALSE  FALSE FALSE
#> 2  FALSE FALSE   TRUE  TRUE  FALSE FALSE
#> 3   TRUE FALSE  FALSE FALSE  FALSE  TRUE
```
