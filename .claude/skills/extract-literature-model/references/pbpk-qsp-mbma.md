# PBPK / QSP / MBMA discipline

Load this reference when the source paper describes a whole-body PBPK, semi-mechanistic PBPK, quantitative systems pharmacology (QSP), or model-based meta-analysis (MBMA) model. These mechanism papers are first-class members of nlmixr2lib when the source provides enough information to reproduce the model end-to-end, but the discipline for sourcing parameters and recording provenance is **stricter** than for popPK extractions.

Standard popPK conventions still apply. This file adds the extra rules.

## `description` prefix and filename suffix

Prepend the modeling framework to `description` so users filtering by mechanism class see it in `modellib()` listings. Examples:

- `"PBPK (whole-body, SimCYP V12.1). APAP disposition in pregnancy ..."`
- `"QSP. 5-lipoxygenase pathway with 33 ODEs, 64 rate laws, 113 parameters ..."`
- `"MBMA. Trastuzumab + T-DXd from 103 published trials; between-study variability only ..."`

When a popPK extraction already exists for the same drug, suffix the filename so the mechanism class is visible at file-tree scan and the new file does not collide with its sibling: `Author_Year_<drug>_pbpk.R`, `_qsp.R`, or `_mbma.R`.

## Parameter sourcing — NEVER substitute from another paper unless that paper is on disk

PBPK and QSP models often have 50–200+ parameters, many of which live in supplements rather than the main text. If a parameter is not reported in any on-disk source (paper text, on-disk supplement, on-disk upstream citation), DO NOT fill in the gap from training-data knowledge, class-typical values, default SimCYP / GastroPlus / Open Systems Pharmacology library entries, or "what a similar PBPK paper used for the same drug class." Instead, fire the Phase 4 missing-parameter sidecar and let the operator decide between (A) author correspondence, (B) operator drops the supplement, or (C) skip the paper.

This rule is **stricter** for PBPK / QSP / MBMA than for popPK because the literature is full of mechanism papers that cite each other circularly — substituting from "a representative PBPK paper for the same drug" introduces silent provenance ambiguity (which paper's parameters? which body-composition assumptions? which permeability model? which between-study covariance prior?) that a downstream user cannot audit. The acceptable substitution path is exactly one: the current paper explicitly states "we used the <Author Year> body-composition / metabolism / permeability model" AND that paper is on disk in the source directory; the inherited parameters then get an inline `# inherited from <Author Year> Table N` comment pointing to the upstream paper.

Class-typical training-data defaults are NOT an acceptable substitute. "kdeg for trastuzumab is typically 0.5/day in the literature" is not a citation; it's a guess. Sidecar instead.

## Reproducibility check before drafting

Walk the paper for these provenance markers:

- Are all structural parameters (Vmax, Km, kcat, organ volumes, blood flows, permeabilities, partition coefficients, between-study covariance, etc.) tabulated explicitly with units?
- Are the ODEs written out, or only described in prose ("we used the standard SimCYP whole-body model")?
- Are the dosing-event handling and observation-mapping rules unambiguous?
- For MBMA: are the per-study weights / variance estimates explicit, or only summary parameter point estimates?

If any of these are insufficient, sidecar-ask the operator before drafting. Do not "fill in the gaps" from PBPK / QSP / MBMA class knowledge.

## PBPK / QSP platforms vs nlmixr2

SimCYP, GastroPlus, PK-Sim/MoBi, OSP, and similar commercial platforms come with built-in whole-body models that users parameterise. nlmixr2 does not have those platforms — when extracting a SimCYP-published model, you must encode the platform's whole-body ODEs explicitly in `model()`. If the paper only states "we used SimCYP's standard whole-body model" without writing out the ODEs and parameter values, the model is not reproducible from on-disk sources and the appropriate action is sidecar-ask or skip — not fill in from training-data knowledge of how SimCYP's standard model is typically parameterised.

## MBMA `between-study variability` ≠ popPK BSV

MBMA papers typically report between-study variance (a study-level random effect on summary metrics) rather than between-subject variance (individual-level η on θ). When extracting an MBMA model:

- Encode between-study variance as a study-level eta (e.g., `eta_study_lcl ~ <var>`) clearly labelled to distinguish from the popPK pattern; do NOT silently relabel it as `etalcl` (between-subject).
- The model is suitable for simulating study-level summary outcomes (per-arm mean response) but NOT individual-level concentrations.
- Document the simulation scope prominently in `description` and in the vignette's Assumptions and deviations section.
