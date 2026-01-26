# Create a Markov model for \`nlmixr2\`

Create a Markov model for \`nlmixr2\`

## Usage

``` r
createMarkovModel(..., ignoreProbLt = 0, ignoreProbGt = 1, transitions = NULL)
```

## Arguments

- ...:

  Passed to \`createMarkovTransitionMatrix()\`

- ignoreProbLt:

  Do not estimate probabilities less than the given probability. Set to
  0 to estimate all probabilities with any transition chance, use
  \`estimateZeroTransitions\` to give an initial estimate when there is
  zero probability of transition into the state.

- ignoreProbGt:

  Do not estimate probabilities on a row where any probability is
  greater than the given probability (treating it as a collector row).
  Set to 1 to estimate all probabilities with any transition chance away
  from the state, use \`estimateZeroTransitions\` to give an initial
  estimate when there is zero probability of transition from the state.

- transitions:

  Use this manually-created transition matrix rather than an
  automatically-creating one

## Value

A template \`nlmixr2\` model as a character string with \`ini()\` and
\`model()\` blocks for the Markov model

## See also

Other Markov models:
[`createMarkovModelDataset()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovModelDataset.md),
[`createMarkovModelFromSingleState()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovModelFromSingleState.md),
[`createMarkovTransitionMatrix()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovTransitionMatrix.md),
[`simMarkov()`](https://nlmixr2.github.io/nlmixr2lib/reference/simMarkov.md),
[`simMarkovId()`](https://nlmixr2.github.io/nlmixr2lib/reference/simMarkovId.md)
