# Simulate a Markov model broken down by ID

Simulate a Markov model broken down by ID

## Usage

``` r
simMarkovId(data, initialState, prCols)
```

## Arguments

- data:

  A data.frame for the individual with columns for each of the state
  transition probabilities

- initialState:

  The initial Markov state

- prCols:

  A named list of probability columns. List names are the previous
  state, and list values are a character vector of probability columns.

## Value

A data.frame with two columns named "prev" and "cur"

## See also

Other Markov models:
[`createMarkovModel()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovModel.md),
[`createMarkovModelDataset()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovModelDataset.md),
[`createMarkovModelFromSingleState()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovModelFromSingleState.md),
[`createMarkovTransitionMatrix()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovTransitionMatrix.md),
[`simMarkov()`](https://nlmixr2.github.io/nlmixr2lib/reference/simMarkov.md)
