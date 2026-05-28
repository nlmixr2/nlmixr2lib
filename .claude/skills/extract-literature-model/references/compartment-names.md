# Compartment naming for nlmixr2lib models

Authoritative source: `vignettes/create-model-library.Rmd` and `R/conventions.R`. When in doubt, prefer this file; if this file conflicts with `create-model-library.Rmd`, raise the conflict with the user.

**Stop-and-ask gate (Phase 1 pre-flight + Phase 3 drafting):** If the model you are extracting needs a compartment role that is NOT in this document, file a sidecar BEFORE writing the model file. Propose the canonical role-based name, its source-paper local name(s), and any cross-precedent in existing registered models. Wait for operator approval. Never introduce numbered placeholders like `cmt1` / `compartment_3` / `c1` silently — those are red flags that the role hasn't been canonicalised. Trivial casing differences (the paper's `Central` → canonical `central`) translate silently and do NOT need a sidecar.

## Compartments

Lower case. Snake case only when combining concepts.

- `depot` — extravascular dosing compartment (oral, SC, IM).
- `central` — IV / central sampling compartment.
- `peripheral1`, `peripheral2` — peripheral compartments for 2- and 3-compartment models.
- `effect` — effect compartment for PK/PD models.
- `transit1`, `transit2`, … — transit-compartment absorption chains.
- `precursor1`, `precursor2`, … — precursor / maturation chains (e.g., platelet maturation in Petrov 2024 romiplostim).
- `lat0`, `lat1`, `lat2`, … — latency chains (e.g., CSF/BBB delay in LeTilly 2021 trastuzumab).
- `target` — free (unbound) target species in explicit-binding TMDD models.
- `complex` — drug–target complex in explicit-binding TMDD models (Mager & Jusko 2001).
- `total_target` — total (free + bound) target in QSS / MM TMDD approximations where the bound species is not carried as a separate state (Gibiansky et al. 2008).
- `liver`, `renal_cortex`, `csf`, `isf`, `cumhaz` — semi-physiologic / hazard compartments registered in `R/conventions.R`.
- Therapeutic-area or mechanism-specific compartments: open a GitHub issue before adding new names.

The `target` / `complex` / `total_target` names follow the convention proposed by @iamstein in review of PR #60 and are standard in the TMDD literature.

### Numbered compartments and chain prefixes

Generic numbered compartments (`central1`, `central2`, `cmt1`, `comp_a`) are
**forbidden** — use the canonical `central` / `peripheral1` / `peripheral2`
names. The only numeric suffixes the conventions allow are:

1. The blessed chain prefixes: `transit<n>`, `effect<n>`, `precursor<n>`, `lat<n>`, `depot<n>` (for parallel-absorption / IOV models with multiple absorption routes — `depot` + `depot2` is the common idiom).
2. DAR-numbered ADC isoforms: `dar<n>_<base>` (e.g., `dar0_central`, `dar3_peripheral1`). Each `dar<n>` represents a biologically distinct drug-antibody-ratio isoform of the parent ADC, not arbitrary numbering.
3. Metabolite-suffixed compartments: `<canonical>_<metab>` (see § Parent drug + metabolite below).
4. Target species in physiologic compartments: `target_<location>` / `complex_<location>` where `<location>` is `csf` or `isf` (the `targetLocationRegex` in `R/conventions.R`).

### Parent drug + metabolite compartments (ADCs and similar)

For models that track a parent drug AND one or more metabolites
(antibody-drug conjugates with payload, mAbs that release a small
molecule, etc.), the **parent uses the canonical names** and the
**metabolite carries a lowercase suffix** drawn from the registered
metabolite list:

```
# Parent ADC
central
peripheral1
peripheral2

# MMAE metabolite (suffix at end)
central_mmae
peripheral1_mmae
```

Registered metabolite suffixes are defined in `R/conventions.R::registeredMetabolites` (`mmae`, `dxd`, `sn38`, `dm4`, `medm4`, `mcmmaf`, `complex`, `ige`, `tab`, `nab`, `dar0..dar8`, and many more). Add a new metabolite by editing `R/conventions.R::registeredMetabolites`; do not introduce ad hoc prefixes (`adc_central`, `mmae_central`).

### Do not declare compartment order

Do **not** add explicit `cmt(depot)`, `cmt(central)`, … declarations at the top
of `model()`. rxode2 / nlmixr2 infer compartment order from the order the
`d/dt(...)` lines appear, and user data can set the compartment via the `cmt`
column on dose rows. Explicit `cmt()` declarations are almost never needed; they
clutter the model and can mask ordering bugs. Only add them when the source
paper *requires* a specific numeric order that isn't what the ODE declarations
would imply (very rare) and flag the reason in a comment.

## Observation variable

- `Cc` — concentration in the central compartment. Use even with `linCmt()`: `Cc <- linCmt()`.
- Multi-output PK metabolite outputs use `Cc_<metab>` (e.g., `Cc_mmae`, `Cc_dxd`, `Cc_sn38`, `Cc_tab`, `Cc_dar0`, `Cc_dm4`). Do **not** use the deprecated `C<metab>` form (`Cmmae`, `Cdxd`).
- Non-PK multi-output variables (e.g., `tumorSize`, `freeIgE`, `totalIgE`, `Cbrain_cerebellum`, `Ccsf`) are paper-named and exempt from the `Cc_<metab>` rule.
- Apply residual error to each output by name. Multi-output residual error uses parameter-name-then-output-suffix: `propSd_<output>` / `addSd_<output>` / `expSd_<output>` (parent `Cc` uses the bare suffix-free form `propSd` / `addSd` / `expSd`).
