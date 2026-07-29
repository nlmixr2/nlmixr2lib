# Replicate the author's modeling structure

The default policy for `extract-literature-model` is to build the model **exactly as the authors built it**. Almost always allow the original authors' added complexity; the library exists to faithfully reproduce the literature, and collapsing or splitting structure the authors chose is a loss of fidelity.

This file documents the policy with worked examples. The decision logic itself lives in `SKILL.md` Phase 1 step 10 and the "Standing default" block in `pre-flight-checklist.md`.

## Decision logic

| Author's structure | Extraction shape | Vignette |
|---|---|---|
| N independent models on different cohorts / endpoints / species | N `.R` files | 1 vignette per paper |
| One joint coupled model (parent + metabolite simultaneously, N species sharing parameters, single multi-output) | 1 `.R` file | 1 vignette |
| Base + final in a model-development paper | 1 `.R` file (final only) | 1 vignette |
| Sensitivity analyses / robustness checks the authors did not report as final | exclude | n/a |

The number of `.R` files mirrors the author's modeling design. The vignette is always **one per paper**, regardless of N, because it walks the paper's narrative as a unit.

## When to sidecar

Only two narrow conditions:

1. **Infeasible** — replicating the structure requires NONMEM (or other-platform) features that rxode2 / nlmixr2 cannot express (e.g. a specific mixture-model construct, a censoring scheme without a clean nlmixr2 analogue). Sidecar with the specific feature + a proposed simplification.
2. **Genuinely ambiguous from the paper text** — the text leaves it unclear whether a second parameterisation is intended as a final result or as a robustness check / obsolete iteration. Sidecar with the textual evidence for each interpretation.

Do NOT sidecar to ask "the paper has multiple models, which do you want?" That question is answered by the policy above. The recommended answer is almost always "extract all of them as the authors built them, in N files, with one vignette tying them together."

## Worked example: N-file extraction

**de Vries Schultink 2018** (J Pharmacokinet Pharmacodyn) reported two independent PD models on different cohorts:

- Anthracycline → troponin T K-PD model (one cohort, one drug, one biomarker).
- Trastuzumab → LVEF effect-compartment turnover model (a different cohort, a different drug, a different biomarker).

The authors did NOT fit a joint cardiotoxicity model — these were two parallel non-hierarchical fits sharing only the cardiotoxicity theme of the paper.

Per this policy, the extraction is:

```
inst/modeldb/specificDrugs/de_vries_schultink_2018_anthracycline.R
inst/modeldb/specificDrugs/de_vries_schultink_2018_trastuzumab.R
vignettes/articles/de_vries_schultink_2018_anthracycline_trastuzumab_cardiotoxicity.Rmd
```

Both `.R` files set `vignette <- "de_vries_schultink_2018_anthracycline_trastuzumab_cardiotoxicity"` so they point to the same rendered article. The vignette's narrative walks the paper section by section, calling each `modellib()` at the point in the narrative where that model is introduced.

Wrong shape for this paper:

- One file collapsing both models into a single `de_vries_schultink_2018_cardiotoxicity.R` — discards the author's deliberate independent-fit structure.
- Two vignettes (`de_vries_schultink_2018_anthracycline.Rmd` + `de_vries_schultink_2018_trastuzumab.Rmd`) — fragments the reviewer's view of the paper as a unit and duplicates the population / methods narrative.

## Worked example: one-file extraction

**Yoshida 2021 ipatasertib** (joint parent + metabolite model). The authors fit a single coupled model in which the parent compartment feeds the metabolite compartment through metabolic clearance, and the parent + metabolite parameters were estimated simultaneously in one NONMEM run with shared between-subject covariance.

Per this policy, the extraction is:

```
inst/modeldb/specificDrugs/Yoshida_2021_ipatasertib.R
vignettes/articles/Yoshida_2021_ipatasertib.Rmd
```

One file, one vignette. The single `.R` file contains both `central` and `metabolite` compartments and their shared `ini()`. Splitting parent and metabolite into two files would lose the coupling and the shared random-effects structure that defines the model.

Wrong shape for this paper:

- Two files `Yoshida_2021_ipatasertib_parent.R` + `Yoshida_2021_ipatasertib_metabolite.R` — destroys the joint fit; the metabolite parameters are only identifiable in the context of the parent in this paper.

## Filename rules for N-file extractions

When the policy yields N files for one paper, disambiguate the filename suffix by the dimension the authors split along:

| Split dimension | Filename pattern |
|---|---|
| Per drug (independent models on different drugs) | `Author_Year_drug1.R`, `Author_Year_drug2.R` |
| Per endpoint (independent models on different biomarkers) | `Author_Year_<drug>_<endpoint>.R` |
| Per species (preclinical + human-scaled projection) | `Author_Year_<drug>_rat.R`, `Author_Year_<drug>_human.R` |
| Per subpopulation (independent fits on different patient groups) | `Author_Year_<drug>_<population>.R` |

All N files in the set point to the same vignette via `vignette <- "Author_Year_<paper-topic>"`. The vignette's `\VignetteIndexEntry{...}` and YAML `title:` name the paper's overall subject (e.g. `Anthracycline + trastuzumab cardiotoxicity (de Vries Schultink 2018)`), not any single one of the contributed drugs.

## Why this matters

- The paper's choice to separate or couple is itself information about the modeling problem (identifiability, scientific framing, computational tractability). Erasing it erases the modeling-decision context.
- A user who looks up `modellib("Author_Year_drug")` expects the model that paper reported, not a unification or fragmentation by the librarian.
- Splits are easy to combine downstream; combining splits loses are harder to undo. Defaulting to the author's structure is the conservative choice.

## Cross-references

- `SKILL.md` Phase 1 step 10 — the policy in the workflow.
- `pre-flight-checklist.md` "Standing default" block — sidecar triggers.
- `vignette-template.md` — header conventions including the multi-drug human-form title.
- Memory: `build-models-as-authors-built.md` in the operator's memory store.
