# Which covariate interpolation should articles use? (locf vs nocb)

Study date 2026-11-13. Measured on rxode2 5.1.8.9000 (e4cb67f), nlmixr2est 7.2.1 and nonmem2rx 0.1.12.9000, against the articles as of `origin/main` d8a308afa plus the two open branches `fix/lagged-dose-covariate-notes` and `investigate/tsirizani-advan4-ka-switch` (five articles differ from `origin/main`), with nlmixr2lib built from the same tree. This is a study and a proposal; no article or model was changed.

## Question

Many models read a covariate that can change within a subject. Examples are an occasion (`OCC`, between-occasion variability), postnatal age, weight, a per-dose dose level, or a flag such as active haemodialysis. Between two records, the solver needs a value for such a covariate:
- `locf` (rxode2's default): each record's value holds from that record until the next one.
- `nocb`: each record's value holds back to the previous record.

Most library models were fitted in NONMEM. This study asks which rule the articles should use, and what each choice would cost.

## Recommendation

1. **Use `locf` for every article, NONMEM-sourced models included.** Build event tables so that each record carries the value in force from that record onward: a dose carries the occasion it opens, and a visit carries the weight measured at it.
   - This is rxode2's default and the rule nlmixr2 estimates with.
   - It is the convention pattern 14 of the extraction skill already requires for lagged doses (`addlKeepsCov = FALSE`).
   - It changes no current result.
2. **Use `nocb` only for a table built NONMEM's way**, where each record carries the value for the interval that ends at it: a replayed source dataset, or a NONMEM test case. Say so beside the call. Royston 2025 letermovir already does this correctly.
3. **Never solve a table built for one rule with the other.** That mismatch is what this study measures. It moves the median subject's AUC by up to 43% (Jonsson 2011 ethambutol).
4. **State the rule where the result depends on it.** Where the two rules differ by 1% or more in median subject AUC on an article's own table (25 articles today), the article should:
   - pass `covsInterpolation` explicitly;
   - say in its assumptions which rule the table follows.
5. **For covariates that change continuously** (postnatal or postmenstrual age, weight growth, a concentration fed in as a covariate), a step rule is a discretisation. Either:
   - put records densely enough that `locf` and `nocb` agree to the article's tolerance, or
   - use `covsInterpolation = "linear"` where the values are defined at the records, as 10 affected articles already do for some of their calls.

## The two rules and who uses which

- **NONMEM** advances the system from one record to the next using the PK parameters computed at the later record, so a record's covariates apply to the interval that ends at it.
  - `?rxSolve` describes `covsInterpolation = "nocb"` as "the same method that NONMEM uses".
  - `nonmem = TRUE` evaluates `$PK`-like statements "with `time` set to the time of the data record that ends the interval being integrated".
  - No NONMEM installation was available here, so this rests on the nlmixr2 ecosystem's documentation and tools, not on running NONMEM.
- **nonmem2rx**, which imports NONMEM runs and checks them against NONMEM's own predictions:
  - Its `rxSolve()` method defaults to `covsInterpolation = "nocb"` ("using nocb interpolation like NONMEM") and to `nonmem = TRUE`.
  - Its validation against NONMEM's IPRED and PRED solves with `covsInterpolation = "nocb"`.
- **nlmixr2** estimates with `locf`. `foceiControl()` and `saemControl()` both carry `rxControl(covsInterpolation = "locf")` by default, so an nlmixr2 fit already assumes rxode2's default.

The two rules only differ between a record where a covariate changes and the record before it. Two tables describe the same design if each is built for its rule. A NONMEM dataset puts an occasion's value on the dose that opens it, often with a pre-dose sample at the same time, and NONMEM then applies the value back to that sample. A `locf` table puts it on the dose and applies it forward.

So when a NONMEM model is simulated for a new design, what matters is the design, and a `locf` table states the design directly. The one case where `nocb` is the faithful rule is replaying NONMEM's own records.

A NONMEM fit on sparse data may have applied a new occasion's values back across a long record-free stretch, for example to the absorption tail of the previous dose. The fitted parameters then carry that, and no publication reports enough to undo it. On data with a record at or near each change, the two rules agree, as the measurements below show.

## Method

- **Enumeration.** 3,237 articles were scanned in two steps.
  - A static pass over the R chunks kept the 289 that set `covsInterpolation`, or that assign a candidate covariate from an expression involving time, the record, the dose number or the occasion. The candidates are `OCC`, PNA, PMA, PAGE, GA, WT, dose levels, `CONMED_*`, `MULTI_DOSE*`, `CYCLE`, `FED` and similar.
  - Each candidate was then run with every `rxSolve()` call routed through a wrapper (`design/covariate-interpolation-harness.R`). It counts as affected if, in at least one call, a covariate the model reads varies within a subject in the events: 164 articles.
- **Measurement.**
  - The wrapper leaves the article's own solve untouched.
  - For each call that varies a covariate, it takes up to 20 subjects in which the covariate changes. It solves them twice on the article's own solver path and settings, once with `covsInterpolation = "locf"` and once with `"nocb"`.
  - Each subject keeps the etas the article drew, on a `zeroRe()` model, so the two solves differ only in the rule.
  - It records the change in each subject's trapezoidal AUC of the concentration output, and the largest pointwise difference in any output relative to that output's maximum.
  - It also records the time from each covariate change back to the record before it: the stretch over which the rules disagree.
- **What the number means.**
  - Every table was built for `locf`, the default all but one affected article solves with. So the measured change is what the article would lose if its table were solved with the other rule. It shows how much a result depends on the convention, not which rule is "right".
  - The change includes both effects of the rule: on parameters that act between records (`ka`, CL, V), and on which record a lagged dose reads its `f()`, `dur()` and `rate()` from.
- **Exact check.** For Tsirizani 2025 ritonavir, three children were also integrated exactly by matrix exponential under each rule.
  - The ODE path matches the exact `locf` solution to 1e-10.
  - Exact `nocb` differs from exact `locf` by up to 104% at single points and 10.5% in AUC on that table.

## Which articles are affected

164 articles (the full list is `design/covariate-interpolation-articles.csv`), by what varies:

| What varies within a subject | Articles |
|---|---|
| Occasion (`OCC`, or an occasion eta supplied as data) | 83 |
| Physiology or time since an event (PNA, PAGE, WT, HT, AGE, POD, CYCLE, pregnancy, albumin, CrCl, body temperature) | 32 |
| Dose, regimen or exposure input (dose level, cumulative dose, dialysis or wash flags, a perpetrator concentration, exposure metrics) | 49 |

- Source software named in the model files of the 164: NONMEM in 99 (3 also name Monolix and 1 Phoenix); Monolix alone in 5, Phoenix in 5 and nlmixr2 in 1. 31 name none, and 23 articles have no model file of the same name.
- Already explicit about the rule:
  - 30 articles pass it on every compared call: 21 `"locf"`, 8 `"linear"`, and 1 a mix of the two.
  - 2 more pass it on some calls: one `"linear"`, and Royston 2025 letermovir `"nocb"`, on a table built NONMEM's way.
  - The other 132 rely on the default `locf`.
- Not measured:
  - 10 candidate articles had calls the wrapper could not inspect: events it could not identify, or per-study random effects. These are Barnett 2018 coproporphyrin I, Courlet 2023, deWit 2016 everolimus, Gall 2023, Huh 2024, Krishna 2011, Leding 2026, Nguyen 2026 dupilumab, Stachanow 2023 and Zou 2026 pembrolizumab QoL.
  - Neroutsos 2022 busulfan stopped partway: the wrapper's comparison solves reset the random seed, and a later assertion failed. Its one compared call is included.
  - Articles outside the 289 static candidates were not run. An article that varies a covariate within a subject in a way the static pass does not recognise is missed.

## How much the rule moves results

Largest median-subject AUC change over each article's calls:

| Class | < 0.01% | 0.01-0.1% | 0.1-1% | 1-5% | >= 5% | no concentration output |
|---|---|---|---|---|---|---|
| Occasion (83) | 41 | 12 | 12 | 6 | 11 | 1 |
| Physiology / time since event (32) | 4 | 6 | 10 | 6 | 0 | 6 |
| Dose, regimen or exposure input (49) | 17 | 2 | 8 | 1 | 1 | 20 |

The size follows the time from a covariate change back to the previous record. Occasion articles by their largest gap:

| Gap | < 0.1% | 0.1-1% | >= 1% |
|---|---|---|---|
| up to 1 h | 26 | 5 | 0 |
| 1-6 h | 2 | 0 | 1 |
| 6-12 h | 5 | 3 | 2 |
| over 12 h | 20 | 4 | 13 |

So articles that sample near each change are insensitive. The sensitive ones change occasion or dose level across dose-only stretches, typically a run-in with no observations, where the two rules disagree for a whole dosing interval.

Articles where the median subject's AUC moves by 1% or more:

| Article | Class | Covariates that vary | Gap (h) | Median subject AUC | Largest subject AUC | Source |
|---|---|---|---|---|---|---|
| Jonsson_2011_ethambutol | occasion | OCC | 24 | 43.1% | 384.0% | NONMEM |
| Riglet_2020_mycophenolic_acid | occasion | CRCL, OCC | 696 | 12.8% | 68.4% | unstated |
| Royston_2025_letermovir | occasion | OCC | NA | 12.3% | 44.7% | NONMEM+Monolix |
| Francis_2020_lumefantrine | occasion | OCC | 16 | 11.5% | 23.4% | NONMEM |
| Sawe_2025_levofloxacin | occasion | OCC | 24 | 8.2% | 1140.0% | NONMEM |
| Tsirizani_2024_darunavir | occasion | OCC | 24 | 7.2% | 18.9% | NONMEM |
| Tsirizani_2025_ritonavir | occasion | OCC | 12 | 6.9% | 66.4% | NONMEM |
| Kengo_2025_atazanavir_ritonavir_rifampicin | occasion | OCC | 24 | 6.9% | 16.1% | no model file matched |
| Umpierrez_2021_cyclosporine | occasion | OCC | 6 | 6.9% | 16.3% | Monolix |
| Dickinson_2021_dolutegravir | occasion | OCC, PREG | 23.9 | 5.5% | 16.1% | NONMEM |
| Barnett_2018_rosuvastatin | dose, regimen or exposure input | CP_RIF_UM | 0.25 | 5.4% | 8.7% | NONMEM |
| Vet_2016_midazolam | occasion | OCC | 15.5 | 5.3% | 12.1% | NONMEM |
| Upton_2025_cycloserine_clofazimine | occasion | OCC | 24 | 4.5% | 14.6% | no model file matched |
| Miano_2024_tacrolimus | physiology / time since event | POD | 11 | 4.2% | 13.3% | NONMEM |
| Punyawudho_2025_dolutegravir | occasion | OCC | 24 | 4.0% | 41.1% | NONMEM |
| Mockeliunas_2022_linezolid | occasion | OCC | 24 | 3.5% | 8.8% | NONMEM |
| Itohara_2019_tacrolimus | physiology / time since event | POD | 60 | 3.5% | 3.6% | NONMEM |
| Yeung_2025_indomethacin | physiology / time since event | PNA | 4 | 2.7% | 4.0% | unstated |
| Alihodzic_2022_ciprofloxacin | occasion | OCC | 12 | 2.4% | 8.8% | NONMEM |
| Wallender_2021_piperaquine | occasion | OCC, PAGE, WT | 56 | 2.3% | 2.9% | NONMEM |
| Kim_2021_YH12852 | occasion | OCC | 24 | 1.7% | 15.0% | NONMEM |
| Wu_2023_paracetamol | physiology / time since event | PNA | 6 | 1.6% | 1.8% | NONMEM |
| Chen_2020_ciclosporin | physiology / time since event | POD | 12 | 1.1% | 4.0% | unstated |
| Guo_2016_taurocholate_schh | dose, regimen or exposure input | WASH_ACTIVE | 0.25 | 1.1% | 1.1% | Phoenix |
| Verrest_2023_miltefosine_paromomycin | physiology / time since event | DOSE_MF_CUM_MGKG, NEUT | 24 | 1.1% | 1.7% | no model file matched |

"Gap" is the longest time from a covariate change back to the previous record in the measured calls. Royston 2025 has no gap entry because its time column is named differently. Its 12% is the cost of solving its NONMEM-style table with `locf`, the reverse of every other row.

## Cost of each policy

**Recommended (`locf`, with the rule stated where it matters):**
- No result changes. Every affected article except Royston 2025 already solves its step covariates with `locf` (10 use `"linear"` for continuous ones), and Royston's `nocb` is consistent with its table.
- Edits: for the 25 articles in the table above (24 excluding Royston, which is already explicit), add `covsInterpolation = "locf"` to the calls that vary a covariate and one sentence in the assumptions. These are prose and argument changes with identical numbers.
- The extraction skill gains the rule (a short section beside pattern 14), and `vignettes/create-model-library.Rmd` gains a line.
- Optional: a static check, for example in the pairing lint, that flags a time-varying covariate changing across a dose-only stretch without an explicit `covsInterpolation`.

**Alternative (`nocb` for NONMEM-sourced models):**
- Solving today's tables with `nocb` moves the median subject's AUC by 1% or more in 25 articles, and by 5% or more in 12. Keeping their designs would mean rebuilding each table in NONMEM's convention: at least 99 NONMEM-sourced articles reviewed and rebuilt.
- Monolix- and nlmixr2-sourced models would keep `locf`, so the library would carry two conventions, chosen per model.
- The gain is fidelity only to how NONMEM treated the source's own records, and no article simulates those.

## Limitations

- The static pass decides which articles were run, so coverage of covariates it does not recognise is incomplete (see "Not measured").
- Up to 20 subjects per call were compared, so the largest subject change is a floor, not the cohort maximum.
- Comparisons run on the article's own path. Two affected articles use `linCmt()`. That path keeps the previous value of a covariate across a lagged dose's record (rxode2 PR 1625), so for them part of the difference may come from `linCmt()` rather than from the rule.
- A NONMEM fit's parameters reflect how NONMEM applied the source's records, which the publications do not report. No simulation rule can recover that.

## Files

- `design/covariate-interpolation-harness.R`: the wrapper. Run it with `NLMIXR2LIB_DIR` set to the repository and `NLMIXR2LIB_LIB` set to a library holding nlmixr2lib built from it. CI does not run it.
- `design/covariate-interpolation-articles.csv`: one row per affected article, with class, covariates, source software, current `covsInterpolation`, largest gap, median-subject and largest-subject AUC change, largest pointwise change, and whether it uses `linCmt()`.
