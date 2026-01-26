# Convert nlmixr2 simulations to Markov states

Simulations are performed per sim.id (if applicable) and per subject id
within the simulation. The \`sim.id\` column must be named exactly that,
and the \`id\` column is detected case-insensitvely.

## Usage

``` r
simMarkov(
  ui,
  initialState,
  states,
  colPrev = "previous",
  colCur = "current",
  ...
)
```

## Arguments

- ui:

  a model fit or simulation object to simulate from or a data.frame
  containing the transition probability columns

- initialState:

  The initial Markov state

- states:

  All states in the model as a character vector

- colPrev, colCur:

  Column names in the dataset for the previous and current states (may
  be any R object that can be coerced to a character string, ordered
  objects are usually preferred)

- ...:

  Passed to methods

## Value

A data.frame with the original data along with the new state columns
(\`colPrev\` and \`colCur\`) added

## Details

See the Markov vignette
([`vignette("markov", package = "nlmixr2lib")`](https://nlmixr2.github.io/nlmixr2lib/articles/markov.md))
for an example.

## See also

Other Markov models:
[`createMarkovModel()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovModel.md),
[`createMarkovModelDataset()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovModelDataset.md),
[`createMarkovModelFromSingleState()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovModelFromSingleState.md),
[`createMarkovTransitionMatrix()`](https://nlmixr2.github.io/nlmixr2lib/reference/createMarkovTransitionMatrix.md),
[`simMarkovId()`](https://nlmixr2.github.io/nlmixr2lib/reference/simMarkovId.md)
