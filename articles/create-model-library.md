# Creating a model library

Model libraries are useful to have consistent high-quality basic models
that can be used as a model itself or as a building block for other
models.

## Model library conventions within nlmixr2lib

Compartment and parameter names should be all lower case when on their
own and should use snakeCase when combined in some way.

Compartment and parameter names are selected to align with those used by
`rxode2::linCmt()` which are described in the vignette:
vignette(“rxode2-model-types”, package = “rxode2”).

### Compartment naming

Compartment naming follows compartment names with the `linCmt()` with
augmentation for other compartments:

- `depot`: The extravascular dosing compartment (for example, the gut
  for oral dosing or subcutaneous space for subcutaneous dosing)
- `central`: The intravascular compartment used for intravenous dosing
  or for typical pharmacokinetic (PK) model sampling of the drug
- `peripheral1`, `peripheral2`: The first and second peripheral
  compartments for 2- and 3-compartment PK models
- `effect`: The compartment for effect compartment models
- Therapeutic-area-specific models should use consistent compartment and
  parameter naming. When adding a new therapeutic area model to the
  library, please discuss naming first in a new GitHub issue.

### Estimated parameter naming

To enable more consistent cross-model compatibility, the following
conventions should be used unless there is a strong reason for an
exception:

- Pharmacokinetic concentrations in the central compartment should be
  named `Cc`. `Cc` should be used even when using a `linCmt()` model (in
  which case `Cc <- linCmt()` should be used and the residual error
  should be applied to the `Cc` parameter).
- Therapeutic-area-specific models should use consistent compartment and
  parameter naming. When adding a new therapeutic area model to the
  library, please discuss naming first in a new GitHub issue.

### Parameter naming

PK models should use the following parameter naming conventions:

- `ka`: absorption rate
- `cl`: clearance
- `q`: intercompartmental clearance (`central` to and from `peripheral1`
  compartments)
- `q2`: second intercompartmental clearance (`central` to and from
  `peripheral2` compartments)
- `vc`: central volume of distribution
- `vp`, `vp2`: first and second peripheral compartment volumes

When micro-constants are used, they should use the following naming
conventions:

- `kel` elimination rate (`cl/vc`)
- `k12`, `k21`, `k13`, `k31`: intercompartmental transit rates (`q/vc`,
  `q/vp`, `q2/vc`, and `q2/vp2`, respectively)

### Parameter transforms

Parameters are often estimated on a transformed scale. For instance, a
natural logarithm transform is often used for parameters that must be
positive, and a logit transform is often used when a parameter must
remain within a specific range.

Transformed parameters should be prefixed with an indicator of the
transformation. Preferred transformation prefixes are:

- `l` (lower case L): natural log transform
- `logit`: logit transform
- `probit`: probit transform
- Any other transform should similarly include the full transform as the
  prefix

Generally, for any transform other than natural logarithm, include the
full name as a prefix. For example, natural logarithm-transformed `ka`
would be `lka` and logit-transformed `emax` would be `logitemax`.

### Random effects

Random effects are estimates as part of a distribution varying by some
grouping factor. The grouping factor is often a subject in a clinical
trial. (For NONMEM users, random effects are often referred to as
inter-individual variability.)

Random effect parameters should prefix the (transformed) parameter name
with `eta`. For example, a random effect on log-transformed clearance
would be named `etalcl`.

### Drug effects

Different drug effects may be investigated during model building. And,
multiple drug effect styles (linear, E_(max), threshold, etc.) may be
investigated by the user.

To enable simpler changes to drug effects and to minimize the chance of
parameter name collisions when combining models, the following rules are
strongly recommended:

- Drug effects should be calculated on a model by themselves to enable
  changing the type of drug effect (e.g. linear to E_(max)).
- The parameter name for the drug effect should be called `drugEffect`
  followed by the name of the part of the model that is most closely
  associated with the drug effect. For example, in the Simeoni 2004
  model, the drug effect is called `drugEffectCyclingCells`.

## Model files

Files in a model library should have the following characteristics:

- The first line inside the function should have a description
  assignment. That is
  `description <- "This is the description of the model"` right inside
  the `function()` before the `ini({})` block.

- If the model has a literature reference associated with it, then the
  second line of inside the function should have the reference, for
  example,
  `reference <- "Richard Hooijmaijers, Matthew Fidler, William S. Denney (2022). nlmixr2lib: A Model Library for 'nlmixr2'. https://nlmixr2.github.io/nlmixr2lib/"`

- If it would be helpful to give the user some information about the
  model on load, it can be added as meta-data as a `"message"` attribute
  to the model. Note that in that case, you must give the function name
  as the last line of the model to ensure that it is the returned value
  from evaluation of the file. (See `oncology_xenograft_simeoni_2004.R`
  for an example of adding a message.)

- If the model is to be combined with other models and it expects
  certain objects will be defined a depends value should be specified.
  For example if there is a tumor growth model that is driven by the
  drug concentration in the central compartment, then the following
  could be used: `depends <"Cc"`

- It can also be helpful to specify the compartments where dosing is
  expected. This can be done in the following manner:
  `dosing <- c("central", "depot")`

- Units used in the model can be specified using a list
  `units <- list(dosing= "mg/kg", time="hr")` To add more fields to this
  list please discuss first in a GitHub issue.

- The remainder of the file should be an nlmixr2 model in a function
  with a typical
  [`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html) and
  [`model()`](https://nlmixr2.github.io/rxode2/reference/model.html)
  block.

- The name of the file should match the name of the model within the
  file.

If a function to modify, self-start, or otherwise help the user would
make sense, add it as a new file in the `R/` directory with the file
name and function name `updateModelName()` using the word update
followed by the model name in camelCase
(e.g. `updateOncologyXenograftSimeoni2004`). If such a function is
added, please add it in the `messages` described above, as well. Update
functions must be able to take in a function, an rxUi object, or an
nlmixr2fitCore object and should usually return an rxUi object.

For examples, see the package installation directory.
