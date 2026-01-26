# Create a Markov transition matrix with probabilities of transitioning between every state

Create a Markov transition matrix with probabilities of transitioning
between every state

## Usage

``` r
createMarkovTransitionMatrix(
  colPrev,
  colCur,
  estimateZeroTransitions = FALSE,
  estimateZeroTransitionsInitial = FALSE,
  ...
)
```

## Arguments

- colPrev, colCur:

  Column names in the dataset for the previous and current states (may
  be any R object that can be coerced to a character string, ordered
  objects are usually preferred)

- estimateZeroTransitions:

  Should transitions that have zero occurrences be estimated? This is
  done by setting the state to have a single transition.

- estimateZeroTransitionsInitial:

  Should transitions that are only initial states be estimated (ignored
  if \`estimateZeroTransitions = FALSE\`)

- ...:

  Ignored

## Value

A square matrix with row and column names for each state where rows are
the prior state and columns are the current state.

## See also

Other Markov models:
[`createMarkovModel()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovModel.md),
[`createMarkovModelDataset()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovModelDataset.md),
[`createMarkovModelFromSingleState()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovModelFromSingleState.md),
[`simMarkov()`](https://nlmixr2.github.io/nlmixr2lib/reference/simMarkov.md),
[`simMarkovId()`](https://nlmixr2.github.io/nlmixr2lib/reference/simMarkovId.md)
