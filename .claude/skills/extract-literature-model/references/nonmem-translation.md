# NONMEM → nlmixr2 syntax translation

Load this reference only when the source paper supplies a NONMEM control stream (`.mod` / `.ctl`, supplement, regulatory review with `$PK`/`$DES` blocks). For popPK papers that report parameter tables without exposing NONMEM syntax, this file is not needed.

## Parameter symbols

| NONMEM | nlmixr2 | Notes |
|--------|---------|-------|
| `THETA(i)` | named log-transformed parameter (`lka`, `lcl`, `lvc`, …) | Always log-transform positive structural parameters; the bare value (`exp(lcl)`) reappears inside `model()`. Choose the canonical name from `parameter-names.md` § "Structural PK parameters"; never carry over `THETA1`, `THETA2`, … as parameter names. |
| `ETA(i)` | `eta` + transformed parameter (`etalcl`, `etalvc`, …) | One `eta*` per `ETA(i)` slot referenced in `$PK` / `$DES`. Block `$OMEGA` translates to a correlated-IIV form (see below). |
| `EPS(i)` | `propSd` / `addSd` / `expSd` (or `propSd_<output>` / `addSd_<output>` / `expSd_<output>` for multi-output) | Map `EPS(i)` slots to residual-error parameters by their role in `$ERROR` (multiplicative vs additive vs combined vs log-normal). |

## `$OMEGA` blocks

- **Diagonal** (no `BLOCK`): each line is an independent `eta* ~ var` ini line.
  ```
  $OMEGA  0.09        ; IIV CL
          0.16        ; IIV VC
  ```
  →
  ```r
  etalcl ~ 0.09  # NONMEM $OMEGA line 1, IIV CL
  etalvc ~ 0.16  # NONMEM $OMEGA line 2, IIV VC
  ```
- **Block** (`$OMEGA BLOCK(n)`): each block becomes a correlated `eta1 + eta2 + ... ~ c(...)` line. NONMEM stores the lower triangle row-by-row; nlmixr2 expects the same lower-triangle order in `c(...)`.
  ```
  $OMEGA BLOCK(2)
   0.09           ; IIV CL
   0.05  0.16     ; cov(CL,VC), IIV VC
  ```
  →
  ```r
  etalcl + etalvc ~ c(0.09, 0.05, 0.16)  # NONMEM $OMEGA BLOCK(2)
  ```
- `$OMEGA BLOCK(n) FIXED` → wrap the whole `c(...)` in `fixed(...)`: `etalcl + etalvc ~ fixed(c(0.09, 0.05, 0.16))`.
- `$OMEGA BLOCK(n) SAME` (block re-uses the previous block's values) → spell out the values; nlmixr2 has no `SAME` shortcut.

## CV%, variance, log-vs-linear

- `omega²` in `$OMEGA` is variance on the **internal** (log) scale for log-normal parameters.
- Convert paper-reported CV%: `omega² = log(1 + CV²)`. Always carry the conversion in a comment, e.g. `etalcl ~ 0.0905  # log(1 + 0.30^2) per Table 3 CV% = 30%`.
- A `THETA` is itself the back-transformed value (e.g., `THETA(1) = 0.225`); wrap it in `log(...)`: `lcl <- log(0.225)`. **Never** write `lcl <- 0.225` because that interprets 0.225 as `log(CL)` rather than `CL` itself.
- `$THETA  0.225 FIXED` → `lcl <- fixed(log(0.225))` in `ini()`.

## `$SUBROUTINE ADVAN…` mapping

| ADVAN | nlmixr2 path | Notes |
|-------|--------------|-------|
| `ADVAN1` (1-cmt IV) | `linCmt()` with `cl` + `vc` | TRANS1 / TRANS2 unified — both pick the `cl, vc` parameterization in nlmixr2. |
| `ADVAN2` (1-cmt oral, 1st-order absorption) | `linCmt()` with `ka` + `cl` + `vc` | |
| `ADVAN3` (2-cmt IV) | `linCmt()` with `cl, vc, q, vp` | TRANS3 / TRANS4 / TRANS5 / TRANS6 — pick the parameterization the source uses; `linCmt()` accepts micro-constants (`k12`, `k21`) too. |
| `ADVAN4` (2-cmt oral) | `linCmt()` with `ka, cl, vc, q, vp` | |
| `ADVAN11` (3-cmt IV) | `linCmt()` with `cl, vc, q, vp, q2, vp2` | |
| `ADVAN12` (3-cmt oral) | `linCmt()` with `ka, cl, vc, q, vp, q2, vp2` | |
| `ADVAN5` (general linear, user `$MODEL`) | explicit `d/dt(...)` ODEs | Map each `$MODEL` `COMP=` to a named compartment (canonical names where possible). |
| `ADVAN6` / `ADVAN8` / `ADVAN9` / `ADVAN13` (general non-linear ODE solvers) | explicit `d/dt(...)` ODEs in `model()` | Translate each `$DES` line into a `d/dt(<state>)` assignment. Solver tolerance differences are not material at the model-file level. |

## `$PK` and `$DES` blocks

- `$PK` block lines computing individual parameters → reproduced inside `model()` after derived covariate terms but before `d/dt(...)`. Drop NONMEM's verbose `TVCL` / `CL` distinction in favor of the nlmixr2 canonical pattern: compute `cl <- exp(lcl + etalcl) * <covariate effects>`.
- `$DES` block → `d/dt(<state>)` lines in `model()`. The state names come from `$MODEL COMP=` (canonical compartment names) or, for `ADVAN6` without explicit `$MODEL`, from the index implied by `A(1)`, `A(2)`, … (operator must declare them).

## `$ERROR` block patterns

| NONMEM | nlmixr2 |
|--------|---------|
| `Y = F + EPS(1)` | `Cc ~ add(addSd)` |
| `Y = F * (1 + EPS(1))` | `Cc ~ prop(propSd)` |
| `Y = F * (1 + EPS(1)) + EPS(2)` | `Cc ~ prop(propSd) + add(addSd)` |
| `Y = LOG(F) + EPS(1)` (log-transform-both-sides) | `Cc ~ prop(propSd)` (NONMEM "additive on log scale" ≡ proportional in linear space) — **or** `Cc ~ lnorm(expSd)` when the paper reports the residual as the log-normal SD directly. |
| Conditional `IF (CMT.EQ.1) Y = … ELSE Y = …` (multi-output) | per-output residual lines: `Cc_<metab> ~ prop(propSd_<metab>)` etc. |
| `IPRED=F`, `IRES=DV-IPRED`, `IWRES=IRES/W` | not needed — derived in `model()` only when the vignette consumes them; otherwise omit. |

## Scale factors and bioavailability

- `S1 = V/1000` (or similar) is a NONMEM unit-rescaling shorthand: dose in mg, volume in L, concentration in mg/L is the canonical case. Translate by **declaring the units explicitly in `units$concentration`** (e.g., `"mg/L"`) and dropping the `S1`. Add a unit comment if the source's unit choice is non-obvious.
- `F1 = THETA(k)` (bioavailability on compartment 1) → `f(depot) <- exp(lfdepot)` (or `f(<compartment>) <- ...` for a non-depot target). Confirm which compartment NONMEM is targeting by reading `$MODEL`.
- `ALAG1 = THETA(k)` (lag time on compartment 1) → `lag(depot) <- exp(llag)`.
- `D1 = THETA(k)` (zero-order infusion duration on compartment 1) → `dur(depot) <- exp(ldur)`.
- `R1 = THETA(k)` (zero-order infusion rate on compartment 1) → `rate(depot) <- exp(lrate)`.

## Covariate effects

- `TVCL = THETA(1) * (WT/70)**THETA(2)` → power-form effect in `model()`:
  ```r
  cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
  ```
  with `lcl <- log(<THETA1>)` and `e_wt_cl <- <THETA2>` in `ini()`.
- `TVCL = THETA(1) * EXP(THETA(2) * (AGE - 40))` → exponential form: `cl <- exp(lcl + etalcl) * exp(e_age_cl * (AGE - 40))` with `e_age_cl <- <THETA2>`.
- `IF (SEX.EQ.1) TVCL = TVCL * THETA(2)` → multiplicative categorical effect: `cl <- exp(lcl + etalcl) * (1 + e_sex_cl * SEX)` (binary `SEX`) or `cl <- exp(lcl + etalcl) * exp(log(<THETA2>) * SEX)` if the paper uses the multiplicative-from-1 form. Match the paper's convention; mark the reference category in `covariateData[[<COV>]]$reference_category`.
- Hill / sigmoidal covariate forms (e.g., maturation): translate the symbolic form directly; do not "simplify" to a generic shape that loses parameter identifiability.

## NONMEM constants and shortcuts

- `A(i)` — `i`-th compartment amount; replaced by the canonical compartment state in `d/dt(<state>)`. The `A(1) = ...` initialization in `$DES` translates to `<state>(0) <- ...` in `model()`.
- `T` — current time; nlmixr2 exposes `t` (lowercase). Source `IF(T.GE.X)` becomes `if (t >= X)`.
- `MIXEST` (mixture model index) — nlmixr2 mixture support is more limited; sidecar-ask the operator before extracting a mixture model.
- `$MODEL TOL=` — NONMEM solver tolerance, ignored at the nlmixr2 model-file level.
- `$ESTIMATION METHOD=…` — irrelevant for the model file (estimation method is supplied at fit time, not in the model definition).

When in doubt about a NONMEM construct that is not in this table, sidecar-ask the operator with the verbatim source line; do not invent translations.
