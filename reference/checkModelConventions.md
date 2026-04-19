# Check a model against nlmixr2lib conventions

Parses a model and reports deviations from the nlmixr2lib conventions
documented in \`vignettes/create-model-library.Rmd\` and the
\`extract-literature-model\` skill references (especially
\`naming-conventions.md\` and \`inst/references/covariate-columns.md\`).
The checker inspects: file-level metadata (description, reference,
units, covariateData); parameter names (log-prefix PK params,
\`eta\`-prefix IIV, \`propSd\`/\`addSd\` residual error); parameter
labels; covariates (canonical register, units, declared aliases);
compartment names; the observation variable (\`Cc\`); and a syntactic
dosing-vs-concentration unit cross-check.

## Usage

``` r
checkModelConventions(model, verbose = TRUE)
```

## Arguments

- model:

  Model to check. Accepted forms: - character scalar: resolved via
  \[readModelDb()\]. - function: evaluated via
  \`nlmixr2est::nlmixr()\`. - \`rxUi\` object: used directly. - missing
  / \`NULL\`: iterate over every model in \`modeldb\`.

- verbose:

  Logical; if \`TRUE\` (default) print the per-category report to the
  console via \`cli\`.

## Value

Invisibly, a \`data.frame\` with one row per issue and columns
\`model\`, \`category\`, \`severity\`, \`name\`, \`message\`,
\`suggestion\`.

## Details

When any issue of severity \`"error"\` or \`"warning"\` is found, a
\`warning()\` is emitted whose message directs the user to call this
function for the full report.

## Examples

``` r
if (FALSE) { # \dontrun{
checkModelConventions("PK_1cmt_des")
checkModelConventions()            # check every model in modeldb
} # }
```
