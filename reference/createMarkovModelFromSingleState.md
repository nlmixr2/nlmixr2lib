# Create the parts of a Markov model for transitioning from a single state

Create the parts of a Markov model for transitioning from a single state

## Usage

``` r
createMarkovModelFromSingleState(transitionRow, stateNames)
```

## Arguments

- transitionRow:

  a single-element named list with a named vector of all transitions

- stateNames:

  a named vector of state names where the name is the name for use in
  the model parameters

## Value

A list with two elements, "ini" and "model", where each element is a
character vector of lines of code for the model

## See also

Other Markov models:
[`createMarkovModel()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovModel.md),
[`createMarkovModelDataset()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovModelDataset.md),
[`createMarkovTransitionMatrix()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovTransitionMatrix.md),
[`simMarkov()`](https://nlmixr2.github.io/nlmixr2lib/reference/simMarkov.md),
[`simMarkovId()`](https://nlmixr2.github.io/nlmixr2lib/reference/simMarkovId.md)
