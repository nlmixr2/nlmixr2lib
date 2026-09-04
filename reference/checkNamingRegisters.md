# Check the naming registers for structural defects

Validates \`inst/references/\*.md\`, the registers that define canonical
covariate, parameter and compartment names, rather than the models that
use them. Complements \[checkModelConventions()\], which checks the
reverse direction.

## Usage

``` r
checkNamingRegisters(root = NULL)
```

## Arguments

- root:

  Package root to check. \`NULL\` (default) uses the installed package.

## Value

A data.frame with one row per issue and columns \`register\`, \`check\`,
\`name\`, \`line\`, \`detail\`. Zero rows means clean.

## Details

Checks performed:

- duplicate canonical entries (globally for covariates; per section for
  compartments, where a token may legitimately be both a bare
  compartment and a suffix);

- citations of model files that do not exist on disk;

- canonicals that no model uses;

- \`\[\[cross-references\]\]\` with no matching entry;

- entries with no \`Example models:\` line, outside an explicit
  exemption list of universal parameters;

- entries with no \`Type:\` line, or whose \`Type:\` is outside the
  known vocabulary. A missing \`Type:\` is the quietest defect of the
  set: the entry drops out of the canonical name list
  \[checkModelConventions()\] builds, so the model-facing check stops
  recognising the name while every other register check stays green.
