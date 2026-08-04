# compartmentData backfill (issue #482)

## Status

**Backfilled 2026-08-03: 1,463 of 1,518 models, 6,219 compartment entries.**
The remaining models either have no ODE states (nothing to annotate) or are the
two exceptions listed below.

Every backfilled entry carries `verified = FALSE`. **That is not a formality.**
None has been checked against its source paper. `verified = TRUE` appears only
on the three worked examples that were hand-checked
(`PerezRuixo_2025_posdinemab`, `HuttonSmith_2018_ranibizumab`,
`Le_2015_lampalizumab_cyno`).

## How it was produced

| Field | Source | Why |
|---|---|---|
| `units` | the model's own `units` block (amount basis) | The local model got this wrong 20/20 in testing -- it returned the *concentration* string for states that hold *amounts*. It is derivable, so it is never asked for. |
| `specimen` | curated map first, local model only for unknown states | The map is authoritative and suffix-aware, so `dar8_peripheral1` and `total_target_peripheral1` inherit the base compartment's matrix. |
| `analyte` | filename for pure-PK models; local model otherwise | The one field that genuinely needs the paper. |

Two tiers:

- **769 models fully deterministic.** Every state is a drug-disposition state,
  so one analyte is provably correct across all of them. No model call.
- **689 models via the local model**, reconciled against the authoritative
  state list -- its state *names* are advisory only, because it autocorrects
  (`nows` -> `NOWS`, `aroff` -> `arooff`, `nmta_f` -> `nm_ta_f`). 685 resolved;
  the other 4 were written by decoding their own systematic naming schemes.

## What the spot-checks caught

Recorded because each would have written plausible, schema-valid, build-passing
data that was wrong:

1. **`Friberg_2002_paclitaxel`.** `circ` holds circulating neutrophils and
   `precursor1-4` progenitor cells, but `circ` maps to a real matrix (whole
   blood), so under the first rule it inherited "paclitaxel" as its analyte and
   umol as its units. Both wrong, across 62 models. Fixed by restricting the
   deterministic tier to models whose states are *all* drug-disposition states.
2. **`peripheral1`.** The local model answered "tissue"; the curated map, issue
   #482's own example table and the hand-written models all say "plasma". Same
   state, different answer depending on which tier processed it -- valid, and
   useless for filtering. Fixed by making the curated map authoritative.
3. **The audit itself.** A first audit asked the model "do all these states hold
   the same analyte?" and reported 0 flags in 40 models. It then failed its own
   validation: it said `Friberg_2002_paclitaxel` was single-analyte. The clean
   result was worthless. Rebuilt to ask for per-state analytes (which the model
   does well) and compare them in code; the rebuilt audit flags Friberg
   correctly and found 1 real issue in 45 sampled models.

**The lesson worth keeping: on this database, a mechanical rule over names needs
a sampled read before it is trusted, and an auditor needs a known-bad case
before its clean result means anything.**

## Not done

- `Anbari_2023_atezolizumab_cibisatamab_qsp` (90 states) and
  `Ippolito_2024_pacmilimab_qsp` (240 states). More than half their states are
  immunological-synapse and antigen-processing constructs
  (`q_syn_T_C1_PD1_PDLN`, `q_A_e_M1p0`). A bespoke decoder for two models risks
  hundreds of wrong entries; a missing block is a warning by design, which is
  exactly this case.
- **Verification against source papers for all 6,219 entries.** This is the
  real remaining work. The `verified` flag makes it enumerable:
  `sum(!verified)` is the size of the backlog.

**When `verified = FALSE` reaches zero, change the missing-block branch of
`.checkCompartmentData()` from `"warning"` to `"error"`.**
