# compartmentData backfill TODO (issue #482, PR fix/issues-474-479 partial scope)

## Status

Issue #482 asked for a machine-readable specimen type per compartment, so that
"restrict to blood, serum or plasma" is a filter rather than a judgement call made
by reading state names. The operator extended the requirement: each entry also
carries the **analyte** and the **state units**, because otherwise that information
is not recoverable for a cross-model analysis either.

**Done in `fix/issues-474-479`:**

- `conventions$specimenVocabulary` (23 terms) and `conventions$compartmentDataFields`
  in `R/conventions.R`.
- `.checkCompartmentData()` in `checkModelConventions()`, which `buildModelDb()`
  calls, so every build validates: one entry per `d/dt()` state, no entries for
  states that do not exist, all four fields present, and `specimen` drawn from the
  vocabulary.
- The schema is in the extraction skill's `model-file-template.md` and
  `verification-checklist.md`, so **new models are required to carry it**.
- Worked examples on the three models issue #482 called out as hardest:
  `PerezRuixo_2025_posdinemab` (serum / CSF / brain ISF, 8 states),
  `HuttonSmith_2018_ranibizumab` (retina / vitreous / aqueous humour, 12 states),
  `Le_2015_lampalizumab_cyno` (vitreous + serum + tissue, 6 states).

**Not yet done:** the backfill. 1,403 models have ODE states and 6,174 `d/dt()`
declarations. The check is registered at **warning** severity for a missing block
so the database still builds; a *malformed* block is already an error.

**When the backfill reaches zero, change the missing-block branch of
`.checkCompartmentData()` from `"warning"` to `"error"`.**

## Why this is not being auto-generated

`specimen` could largely be derived from the state name, and `units` from the
model's `units` list. `analyte` cannot: it is the one field that genuinely
requires the source paper for anything beyond a single-drug model. Two models in
issue #479's own list looked mechanically identical and turned out to differ once
the paper was read, and the same trap applies here -- a confidently wrong
`analyte` is worse than an absent one, because a cross-model analysis would
silently pool the wrong species.

That is what the `verified` field is for. `verified = FALSE` means the entry was
derived from the state name and still needs checking; `verified = TRUE` means
analyte and specimen were confirmed against the paper. **Do not set `TRUE` to
silence the check.**

## Suggested order of work

Sorted by payoff, because the state distribution is very skewed -- 4 names cover
44% of all declarations and 1,408 names appear exactly once.

1. **761 models whose states are all canonical** (`central`, `peripheral1..3`,
   `depot`, `transit*`, `urine`, organ names). `specimen` and `units` are
   mechanical here; only `analyte` needs a look, and for a single-drug PK model it
   is the drug named in the filename. Highest models-per-hour.
2. **202 models that already declare `paper_specific_compartments`.** That list is
   exactly the set whose states are *known* to be outside the canonical register,
   so it is a ready-made work queue for the cases that need the paper.
3. **The remaining ~440 models** with non-canonical states that have not declared
   them.

## Vocabulary notes

Two entries are deliberately not biological matrices:

- `"administration site"` -- depot and transit states, holding drug that has not
  reached a matrix yet.
- `"not applicable"` -- latent, PD and bookkeeping states. A scan found these in
  **642 of the 1,403** models with ODEs (`effect` 89, `precursor1-4` 108 combined,
  `cumhaz` 25, `total_target` 25, `circ` 20, `complex` 16, bacterial
  subpopulations, tumour-size states). Forcing them into a specimen would put
  false data in the field that the whole issue exists to make trustworthy.

If a matrix is genuinely missing from the vocabulary, add it to
`conventions$specimenVocabulary` in the same change that first uses it, rather
than approximating with a near neighbour.

## Related

- Issue #482 and PR `fix/issues-474-479`.
- `inst/references/compartment-names.md` -- the 430-entry canonical compartment
  register. Its `Role:` field already describes the matrix in prose for the
  canonical states, which is a good starting point for step 1 above.
- `inst/references/fixed-provenance-followup.md` -- same shape of document for
  issue #479.
