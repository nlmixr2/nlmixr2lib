# fixed() provenance follow-up TODO (issue #479, PR fix/issues-474-479 partial scope)

## Status

Issue #479 established that `iniDf$fix` is the only machine-readable signal that a
parameter was not estimated from the source study's data, and that several models
say "fixed" / "assumed" / "taken from the literature" in the label while leaving the
parameter estimable. Any consumer that trusts `fix` therefore reads an assumption as
evidence.

PR `fix/issues-474-479` fixed every parameter #479 enumerated by hand, and added
`.checkFixedLabelAgreement()` to `checkModelConventions()` so the class cannot
silently grow. Because `buildModelDb()` calls `checkModelConventions()`, the check
runs on every database build.

**Not yet done (follow-up work):** the 3 parameters in 2 models listed below. They
were found by the new check, not by hand, and were left for a separate change so the
PR stayed reviewable. The check is registered at **warning** severity precisely
because of this backlog -- erroring would fail the build for pre-existing models
rather than for the new mistake the rule exists to catch.

**When this list reaches zero, promote the check to `"error"`** in
`.checkFixedLabelAgreement()` (`R/checkModelConventions.R`) and update the comment
above the severity argument.

## How the check decides

A label claims *this* parameter was not estimated only when the claim opens its own
clause. The first draft of the rule matched the bare word "fixed" anywhere in the
label and flagged 50 parameters across 40 models -- roughly two thirds of them false
positives, because model labels routinely describe a **different** parameter as fixed:

| Label | Why it is not a claim about this parameter |
|---|---|
| `Table IV omega_Ka = 1.41 (Ka pop is fixed but IIV is estimated)` | the *typical value* was fixed; the variance was estimated |
| `IC50 for inhibition of precursor synthesis (Imax fixed to 1)` | `Imax` was fixed, not `lec50` |
| `Apparent DM4 clearance CL_DM4 (L/day; V_DM4 fixed to 1 L)` | `V_DM4` was fixed -- #479 itself lists this one under "not bugs" |
| `Scaling factor K relating eta_V/F to eta_CL/F (correlation fixed to 1)` | the *correlation* was fixed |
| `CMS nonrenal clearance CL_NRCMS (L/h); CL_RCMS structurally fixed at 0` | `CL_RCMS` was fixed |
| `Absorption lag time for Group 1 (h; Group 2 lag is fixed at 0)` | Group 2's lag was fixed |
| `Apparent formation clearance CL_met/F (L/h); fixed effect, no IIV` | "fixed effect" is the mixed-models term; Lehr 2010 reports RSE 6.9% |

Three narrowing rules bring the count to the 15 below, and each is locked in by a
test in `tests/testthat/test-checkModelConventions.R`:

1. **Variance terms are exempt.** A label on an `eta*` parameter almost always
   discusses the fixed status of the corresponding typical value.
2. **An explicit "was estimated" disclaimer wins** over any claim in the same label.
3. **The claim must open its clause** (after `^`, `;`, `(`, `)`, `,` or a dash), and
   `fixed effect(s)` is excluded.

Anything that survives all three is a genuine disagreement between the label and
`fix`. **If a new false-positive class appears, tighten the pattern and add a case to
that test rather than adding an exemption at the call site.**

## Remaining work (3 parameters, 2 models)

The 12 parameters whose provenance could be established were cleared in
PR `fix/issues-474-479`. `Archary_2019_abacavir` / `Archary_2019_lamivudine`
("Allometric exponents were fixed to 0.75 for CL/F and 1 for apparent volume of
distribution"), `Aksenov_2018_uricAcid` ("The value of Rmax for febuxostat was fixed
at 1") and `Caldes_2009_ganciclovir` (CL reported as `7.49 * (CLCR/57)`, a linear form
whose exponent is 1 by construction) were confirmed against their source papers;
`Bergmann_2014_tacrolimus`, `Boer-Perez_2026_piperacillin` and `Dunlap_2025_tacrolimus`
were taken on their self-declaring labels plus the exact 0.75 / 1.0 values;
`Cao_2013_mab7E3` carries over the 0.950 already confirmed for its two siblings.

What remains is only the residual-error group, which is a value decision rather than a
flag decision -- see the group heading below.

| Model | Parameter | Value | Label |
|---|---|---|---|
| `Belldina_2003_cysteamine` | `addSd_cystine` | 0.1 | Additive residual error on WBC cystine output (nmol cystine / mg protein; assumed; paper does not specify a residual error structure for cystine) |
| `Belldina_2003_cysteamine` | `propSd` | 0.15 | Proportional residual error on plasma cysteamine Cc (fraction; assumed; paper Methods specify proportional structure but SD value is not reported) |
| `Toffoli_2001_etoposide` | `propSd` | 0.1 | Proportional residual error (fraction) - ASSUMED from assay validation; not reported in Toffoli 2001 |

## Per-group guidance

### Residual-error parameters assumed by the encoder (3 parameters, 2 models)

`Belldina_2003_cysteamine` (`addSd_cystine`, `propSd`), `Toffoli_2001_etoposide`
(`propSd`).

**Do not wrap these without deciding the policy question first.** These values were
not fixed by the source paper -- they were assumed by the encoder because the paper
reports no residual-error value. The standing library policy for that situation is to
encode unreported residual error as `fixed(0)` plus a vignette erratum, *not* to
substitute a class-typical or assay-precision-derived placeholder. `Toffoli_2001_etoposide`
says so in the label itself ("ASSUMED from assay validation"). Wrapping the current
placeholder in `fixed()` would make a value the paper never reported look deliberate.

The right change for this group is therefore a *value* change (to `fixed(0)` with an
erratum), not a flag change -- which is why it is separated here.

This group is also the case #479 raised in its "Suggested fix" section: the library
cannot currently distinguish **fixed by the source paper** from **assumed by the
encoder**, since both must be spelled `fixed()`. A structured provenance field would
let consumers exclude assumed values from summaries without parsing English. That is
a design change beyond this follow-up, but this group is the evidence for it.

## Related

- Issue #479 and PR `fix/issues-474-479`.
- `inst/references/parameter-names.md`, "Retired names" section.
- `Marcantonio_2022_amivantamab` / `Marcantonio_2022_panitumumab` record
  `propSd = 0.01` with "not paper-derived" in the **units** string rather than the
  label, so the check does not see them. Either move the note into the label or
  extend the check to read `units` as well.
