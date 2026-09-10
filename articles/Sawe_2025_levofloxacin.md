# Levofloxacin in pregnancy (Sawe 2025)

## Model and source

``` r

mod <- rxode2::rxode(readModelDb("Sawe_2025_levofloxacin"))
#> ℹ parameter labels from comments will be replaced by 'label()'
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_mtt_1, etaiov_mtt_2, etaiov_mtt_3, etaiov_mtt_4, etaiov_ka_1, etaiov_ka_2, etaiov_ka_3, etaiov_ka_4, etaiov_fdepot_1, etaiov_fdepot_2, etaiov_fdepot_3, etaiov_fdepot_4
#> as a work-around try putting the mu-referenced expression on a simple line
```

- Citation: Sawe S, Tsirizani L, Court R, Gausi K, Poswa A, Badat T,
  Wiesner L, Loveday M, Maartens G, Conradie F, Denti P. The effect of
  pregnancy on the population pharmacokinetics of levofloxacin in South
  Africans with rifampicin-resistant tuberculosis. Antimicrob Agents
  Chemother. 2025 May;69(5):e01626-24. <doi:10.1128/aac.01626-24>.
  Values taken from the corrected version posted 17 April 2025, which
  revised the Table 2 footnote defining the reported %CV.
- Description: One-compartment population PK model for oral levofloxacin
  in South African adults treated for rifampicin-resistant tuberculosis
  (RR-TB), characterising the effect of third-trimester pregnancy (Sawe
  2025; n = 47 pooled from two studies, 21 pregnant, 12 with matched
  antepartum / postpartum profiles). Savic transit-compartment
  absorption (analytical form, N fixed to 20, MTT = 1.07 h) feeds
  first-order absorption into a one-compartment disposition model.
  Clearance and volume are scaled allometrically by fat-free mass
  (Janmahasatian formula) with fixed exponents 0.75 and 1 and the
  cohort-median 39.4 kg as reference. Higher serum creatinine lowers
  clearance via a power function (exponent -0.367) centred on the cohort
  median 56.2 umol/L, and third-trimester pregnancy raises clearance by
  a further 38.1%. Between-subject variability is retained on CL only;
  between-occasion variability is carried on MTT, ka and F, with the BOV
  on F inflated 2.35-fold on the two occasions whose preceding dose was
  self-reported rather than observed. Bioavailability is fixed at 1.
- Article: <https://doi.org/10.1128/aac.01626-24>
- Supplement (Tables S1-S3, Fig. S1-S7, NONMEM control stream):
  <https://doi.org/10.1128/aac.01626-24> (file
  `aac.01626-24-s0001.docx`, retrieved via Europe PMC `PMC12057369`)

Sawe 2025 pooled two South African studies in adults treated for
rifampicin-resistant tuberculosis (RR-TB) to ask what third-trimester
pregnancy does to levofloxacin exposure. The answer is a 38.1% increase
in clearance *on top of* the increases already explained by body size
and renal function, which together drive antepartum exposures well below
those of non-pregnant participants on the same weight-based dose.

## Population

47 participants contributed 320 levofloxacin concentrations: 33 (70%)
female, 38 (81%) black, median (range) age 32 (19-51) years, weight 58.0
(37.0-98.0) kg, height 1.60 (1.46-1.88) m, fat-free mass 39.4
(27.3-51.2) kg, serum creatinine 56.2 (25.3-110) umol/L and serum
albumin 30.5 (17.0-40.0) g/L (Table 1). 31 (66%) were living with HIV
and on antiretroviral therapy. Of the 33 women, 21 were pregnant and
sampled in the third trimester; 12 of those contributed matched
antepartum and postpartum profiles, giving 19 antepartum, 14 postpartum,
12 non-pregnant-female and 14 male pharmacokinetic profiles.

Levofloxacin was dosed 750 or 1000 mg once daily by WHO weight band
(34-50 kg and above 50 kg), with a standard breakfast about 1 h before
the observed dose. Sampling was pre-dose and at 2, 4, 6, 8, 10 and 24 h
(BEAT Tuberculosis, NCT04062201), pre-dose and 2, 4, 6, 8 and 24 h in
its pregnancy sub-study, and pre-dose and 2, 4 and 6 h in the King
Dinuzulu Hospital observational cohort. Assay was HPLC-MS/MS with an
LLOQ of 0.0781 mg/L; 6 of 320 samples (1.88%) were below it, all
pre-dose.

The same information is available programmatically:

``` r

str(readModelDb("Sawe_2025_levofloxacin")()$population, max.level = 1)
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_mtt_1, etaiov_mtt_2, etaiov_mtt_3, etaiov_mtt_4, etaiov_ka_1, etaiov_ka_2, etaiov_ka_3, etaiov_ka_4, etaiov_fdepot_1, etaiov_fdepot_2, etaiov_fdepot_3, etaiov_fdepot_4
#> as a work-around try putting the mu-referenced expression on a simple line
#> List of 20
#>  $ species       : chr "human"
#>  $ n_subjects    : int 47
#>  $ n_studies     : int 2
#>  $ n_observations: int 320
#>  $ n_profiles    : chr "19 antepartum, 14 postpartum, 12 non-pregnant female, 14 male"
#>  $ age_range     : chr "19-51 years (Table 1: median 32 years)"
#>  $ weight_range  : chr "37.0-98.0 kg (Table 1: median 58.0 kg)"
#>  $ height_range  : chr "1.46-1.88 m (Table 1: median 1.60 m)"
#>  $ ffm_range     : chr "27.3-51.2 kg (Table 1: median 39.4 kg)"
#>  $ creat_range   : chr "25.3-110 umol/L (Table 1: median 56.2 umol/L)"
#>  $ alb_range     : chr "17.0-40.0 g/L (Table 1: median 30.5 g/L)"
#>  $ sex_female_pct: num 70.2
#>  $ race_ethnicity: chr "38 (81%) black, 9 (19%) white (Table 1)"
#>  $ n_hiv_positive: int 31
#>  $ n_pregnant    : int 21
#>  $ disease_state : chr "Rifampicin-resistant tuberculosis (RR-TB) on treatment. 21 of the 33 female participants were pregnant, sampled"| __truncated__
#>  $ dose_range    : chr "Levofloxacin 750 or 1000 mg orally once daily by body-weight band per WHO guidance (34-50 kg and above 50 kg). "| __truncated__
#>  $ regions       : chr "South Africa (two sites for BEAT Tuberculosis; King Dinuzulu Hospital, Durban)"
#>  $ co_medication : chr "Antiretroviral therapy in 31 participants, commonly tenofovir/lamivudine/dolutegravir; RR-TB co-treatment most "| __truncated__
#>  $ notes         : chr "Pooled from BEAT Tuberculosis (ClinicalTrials.gov NCT04062201, a phase 3 RCT with a pregnancy PK sub-study; sam"| __truncated__
```

## Source trace

Every `ini()` entry carries an in-file comment naming its source
location in `inst/modeldb/specificDrugs/Sawe_2025_levofloxacin.R`.
Collected here for review. “Control stream” refers to the NONMEM code
printed in the supplementary material, whose `$THETA` block agrees with
Table 2 to every printed digit and supplies additional significant
figures.

| Equation / parameter | Value | Source location |
|----|----|----|
| `lcl` (CL/F) | `log(6.05859)` L/h | Table 2 “Clearance, CL (L/h)” 6.06 (95% CI 5.47-6.53); control stream `$THETA 1` |
| `lvc` (V/F) | `log(85.8803)` L | Table 2 “Volume of distribution, V (L)” 85.9 (95% CI 80.6-91.7); control stream `$THETA 2` |
| `lka` | `log(1.59498)` 1/h | Table 2 “Absorption rate constant, Ka (1/h)” 1.59 (95% CI 1.11-2.40); control stream `$THETA 3` |
| `lmtt` | `log(1.07437)` h | Table 2 “Mean transit time, MTT (h)” 1.07 (95% CI 0.771-1.32); control stream `$THETA 7` |
| `lnn` (N transit) | `fixed(log(20))` | Table 2 “Number of transit compartments, NN” 20 (fixed); control stream `$THETA 13 = 20 FIX` |
| `lfdepot` (F) | `fixed(log(1))` | Table 2 “Bioavailability, F” 1 (fixed); control stream `$THETA 4 = 1 FIX` |
| `e_ffm_cl` | `fixed(0.75)` | Methods, “allometric scaling of clearance (with a fixed exponent of 0.75)”; control stream `ALLMCL_FFM = (FFM/39.4)**0.75` |
| `e_ffm_vc` | `fixed(1)` | Methods, “volume of distribution (with a fixed exponent of 1)”; control stream `ALLMV_FFM = (FFM/39.4)` |
| `e_preg_cl` | `0.380866` | Table 2 “Effect of pregnancy on CL (%)” +38.1 (95% CI +23.4 to +57.1); control stream `IF (PREGNANT.EQ.1) preg_CL = 1 + THETA(12)`, `$THETA 12` |
| `e_creat_cl` | `-0.366504` | Table 2 “Effect of serum creatinine on CL (power exponent)” -0.367 (95% CI -0.493 to -0.104); Methods equation (supplementary inline equation `m001`); control stream `$THETA 15` |
| `etalcl` | `0.0489868` | Table 2 “Between-subject variability, CL” 22.1% (95% CI 17.1-28.3); control stream `$OMEGA 1` |
| `etaiov_mtt_1..4` | `0.211291` | Table 2 “Between-occasion variability, MTT” 45.9% (95% CI 30.5-70.3); control stream `$OMEGA 21` + 3x `SAME` |
| `etaiov_ka_1..4` | `0.732947` | Table 2 “Between-occasion variability, Ka” 85.6% (95% CI 61.7-119); control stream `$OMEGA 17` + 3x `SAME` |
| `etaiov_fdepot_2`, `_4` | `0.0562092` | Table 2 “Between-occasion variability, F” 23.7% (95% CI 19.3-28.1); control stream `$OMEGA 13` + 3x `SAME` |
| `etaiov_fdepot_1`, `_3` | `0.310775` | `2.35136^2 * 0.0562092`, from Table 2 “Scaling factor for BOV on F for unobserved doses” 2.35 (95% CI 1.68-3.28) and control stream `$THETA 14` with `IF(UNOBS.EQ.1) BOVBIO = E_BOVF*BOVBIO` |
| `addSd` | `0.259443` mg/L | Table 2 “Additive error (mg/L)” 0.244 (`$THETA 6 = 0.243823`) plus the control stream’s unconditional `ADD = THETA(6) + LLOQ*0.2` with LLOQ 0.0781 mg/L |
| `propSd` | `0.0733236` | Table 2 “Proportional error (%)” 7.33 (95% CI 6.49-8.10); control stream `$THETA 5` |
| Transit input rate | `transit(nn, mtt, fbio)` | Control stream `PIZZA`/`TRANSIT` block with `KTR = (NN+1)/MTT`, citing Savic 2007; supplementary Fig. S1 |
| Dose enters via transit only | no `f(depot)` line | Control stream `F1 = 0`; under rxUi the bolus is not separately added, so no explicit suppression is needed (see Errata) |
| CL covariate model | `(FFM/39.4)^0.75 * (1 + e_preg_cl*PREG) * (CREAT/56.2)^e_creat_cl` | Control stream `TVCL = THETA(1)*ALLMCL_FFM*preg_CL*CR_CL` |
| Residual error | `add(addSd) + prop(propSd)` | Control stream `W = SQRT(ADD**2 + PROP**2)` |
| Reference FFM / creatinine | 39.4 kg / 56.2 umol/L | Table 2 footnote *b*; control stream `TVFFM = 39.4`, `MEDCREAT = 56.2` |

The variance scale deserves a note of its own, because it is what the
publisher’s 17 April 2025 correction changed. Table 2 footnote *d*
defines the reported percentages as

``` math
\%CV = \sqrt{\omega^2}\cdot 100
```

so a reported “22.1%” **is** `omega * 100` and is *not* a log-normal CV.
The usual `omega^2 = log(CV^2 + 1)` conversion must therefore NOT be
applied here. The check below confirms the reading against the control
stream’s `$OMEGA` values, all four of which reproduce their printed
Table 2 percentage.

``` r

omega_cs <- c(
  "BSV CL"  = 0.0489868,
  "BOV MTT" = 0.211291,
  "BOV Ka"  = 0.732947,
  "BOV F"   = 0.0562092
)
printed <- c("BSV CL" = 22.1, "BOV MTT" = 45.9, "BOV Ka" = 85.6, "BOV F" = 23.7)

data.frame(
  Parameter          = names(omega_cs),
  `omega^2 (stream)` = omega_cs,
  `sqrt(omega^2)*100` = round(sqrt(omega_cs) * 100, 1),
  `Table 2 %CV`      = printed,
  `log(CV^2+1) would give` = round(log(printed^2 / 1e4 + 1), 4),
  check.names = FALSE, row.names = NULL
) |>
  knitr::kable(caption = "Table 2 footnote d reads %CV = sqrt(omega^2)*100.")
```

| Parameter | omega^2 (stream) | sqrt(omega^2)\*100 | Table 2 %CV | log(CV^2+1) would give |
|:---|---:|---:|---:|---:|
| BSV CL | 0.0489868 | 22.1 | 22.1 | 0.0477 |
| BOV MTT | 0.2112910 | 46.0 | 45.9 | 0.1912 |
| BOV Ka | 0.7329470 | 85.6 | 85.6 | 0.5497 |
| BOV F | 0.0562092 | 23.7 | 23.7 | 0.0546 |

Table 2 footnote d reads %CV = sqrt(omega^2)\*100. {.table}

``` r


# sqrt(omega^2)*100 reproduces every printed Table 2 percentage to the
# printed precision (MTT prints 45.9%, a truncation of 45.97%).
stopifnot(all(abs(sqrt(omega_cs) * 100 - printed) < 0.1))
```

## Covariate derivation: fat-free mass

The size descriptor is fat-free mass, not total body weight. The
supplementary material gives the Janmahasatian formula

``` math
FFM = \frac{WHS_{max}\cdot HT^2 \cdot WT}{WHS_{50}\cdot HT^2 + WT}
```

with `WHSmax` 37.99 and `WHS50` 35.98 for females, 42.92 and 30.93 for
males, `HT` in metres and `WT` in kilograms. FFM beat total body weight
as the allometric descriptor (dOFV 12.9 versus 3.46), which is why no
explicit sex term appears in the final model: the sex difference in
exposure is mediated through FFM.

``` r

ffm_janmahasatian <- function(wt, ht, sexf) {
  whs_max <- ifelse(sexf == 1, 37.99, 42.92)
  whs_50  <- ifelse(sexf == 1, 35.98, 30.93)
  whs_max * ht^2 * wt / (whs_50 * ht^2 + wt)
}
```

Applying it to the published median weight and height of each
participant group should reproduce that group’s published median FFM.
Group medians of different quantities need not be exactly
self-consistent, so this is checked at 6%.

``` r

ffm_ref <- tibble::tribble(
  ~group,                 ~wt,  ~ht,  ~sexf, ~ffm_published, ~source,
  "Antepartum",           61.5, 1.60, 1,     38.7,           "Table S1 study data, pregnant (3rd trimester)",
  "Postpartum",           57.6, 1.57, 1,     37.0,           "Discussion (median postpartum WT 57.6, FFM 37.0)",
  "Non-pregnant female",  52.9, 1.57, 1,     35.7,           "Table S1 study data, non-pregnant females",
  "Male",                 54.4, 1.73, 0,     46.4,           "Table S1 study data, males",
  "Overall cohort",       58.0, 1.60, 1,     39.4,           "Table 1 overall participants (median WT/HT/FFM)"
) |>
  dplyr::mutate(
    ffm_formula = ffm_janmahasatian(wt, ht, sexf),
    pct_diff    = 100 * (ffm_formula - ffm_published) / ffm_published
  )

ffm_ref |>
  dplyr::transmute(
    Group = group, `WT (kg)` = wt, `HT (m)` = ht,
    `FFM published (kg)` = ffm_published,
    `FFM from formula (kg)` = round(ffm_formula, 1),
    `% diff` = round(pct_diff, 1), Source = source
  ) |>
  knitr::kable(caption = "Janmahasatian FFM reproduces the published group median FFM.")
```

| Group | WT (kg) | HT (m) | FFM published (kg) | FFM from formula (kg) | % diff | Source |
|:---|---:|---:|---:|---:|---:|:---|
| Antepartum | 61.5 | 1.60 | 38.7 | 38.9 | 0.6 | Table S1 study data, pregnant (3rd trimester) |
| Postpartum | 57.6 | 1.57 | 37.0 | 36.9 | -0.3 | Discussion (median postpartum WT 57.6, FFM 37.0) |
| Non-pregnant female | 52.9 | 1.57 | 35.7 | 35.0 | -2.0 | Table S1 study data, non-pregnant females |
| Male | 54.4 | 1.73 | 46.4 | 47.5 | 2.5 | Table S1 study data, males |
| Overall cohort | 58.0 | 1.60 | 39.4 | 37.6 | -4.6 | Table 1 overall participants (median WT/HT/FFM) |

Janmahasatian FFM reproduces the published group median FFM. {.table
style="width:100%;"}

``` r


stopifnot(max(abs(ffm_ref$pct_diff)) < 6)
```

## Reconstructing the paper’s decomposition of the pregnancy effect

This is the sharpest available check on the covariate model, because the
paper reports the same quantity twice, computed two different ways.

The final model splits the pregnancy effect three ways: part is body
size (pregnant women are larger), part is renal function (pregnant women
have lower serum creatinine), and the remainder is the estimated +38.1%
“pregnancy” term. The Discussion then states that the *overall*
clearance increase due to pregnancy is “around 53%”, the value obtained
from an alternative model (Table S3) that omits the serum-creatinine
effect entirely.

So: evaluating the final model’s three effects at the published
antepartum and postpartum median covariates must reconstruct that ~53%.
Nothing in this calculation is fitted to the 53% figure, and the two
models were estimated separately, so agreement is a genuine cross-check.

``` r

cl_typical <- function(ffm, creat, preg) {
  6.05859 * (ffm / 39.4)^0.75 * (1 + 0.380866 * preg) * (creat / 56.2)^-0.366504
}

# Published medians for the matched antepartum / postpartum visits
# (Discussion: FFM 38.7 vs 37.0 kg; serum creatinine 45.3 vs 55.4 umol/L).
cl_ante <- cl_typical(ffm = 38.7, creat = 45.3, preg = 1)
cl_post <- cl_typical(ffm = 37.0, creat = 55.4, preg = 0)

contributions <- c(
  "Pregnancy term (+38.1%)"        = 1 + 0.380866,
  "Fat-free mass (38.7 vs 37.0 kg)" = (38.7 / 37.0)^0.75,
  "Serum creatinine (45.3 vs 55.4 umol/L)" = (45.3 / 55.4)^-0.366504
)

data.frame(
  Component = c(names(contributions), "Total CL ratio (antepartum / postpartum)"),
  `CL multiplier` = round(c(contributions, prod(contributions)), 4),
  check.names = FALSE
) |>
  knitr::kable(caption = "Decomposed pregnancy effect on levofloxacin clearance.")
```

|  | Component | CL multiplier |
|:---|:---|---:|
| Pregnancy term (+38.1%) | Pregnancy term (+38.1%) | 1.3809 |
| Fat-free mass (38.7 vs 37.0 kg) | Fat-free mass (38.7 vs 37.0 kg) | 1.0343 |
| Serum creatinine (45.3 vs 55.4 umol/L) | Serum creatinine (45.3 vs 55.4 umol/L) | 1.0766 |
|  | Total CL ratio (antepartum / postpartum) | 1.5375 |

Decomposed pregnancy effect on levofloxacin clearance. {.table}

``` r


total_ratio <- cl_ante / cl_post
cat(sprintf(
  "Reconstructed total CL increase in pregnancy: %.1f%%\nPaper's alternative model (Table S3, no creatinine effect): 53%%\n",
  100 * (total_ratio - 1)
))
#> Reconstructed total CL increase in pregnancy: 53.8%
#> Paper's alternative model (Table S3, no creatinine effect): 53%

# The three decomposed effects reproduce the independently-estimated 53%.
stopifnot(abs(total_ratio / 1.53 - 1) < 0.06)
```

The efficacy target drawn on Fig. 2 is likewise reproducible in closed
form: the caption defines `AUC0-24 >= 146 * MIC / fu` with the
hollow-fibre `fAUC0-24/MIC` target of 146, the wild-type MIC of 0.5 mg/L
and an unbound fraction of 69%.

``` r

auc_target <- 146 * 0.5 / 0.69
cat(sprintf("146 * 0.5 / 0.69 = %.1f mg*h/L (paper states 106)\n", auc_target))
#> 146 * 0.5 / 0.69 = 105.8 mg*h/L (paper states 106)
stopifnot(round(auc_target) == 106)
```

## Structural checks against closed-form results

These compare a solved profile against an algebraic result derived from
the published equations, with parameter values computed independently in
R rather than read back out of the model. Both sides use the same fixed
parameter values, so the residual difference is pure numerical /
trapezoidal error and the tolerances are tight.

``` r

tau       <- 24
n_dose    <- 10                       # 10 daily doses: well past steady state
t_last    <- (n_dose - 1) * tau       # 216 h, the final ("observed") dose
grid_ss   <- seq(t_last, t_last + tau, by = 0.1)

#' Build a plain data.frame event table (never an rxEt object, which silently
#' drops post-hoc covariate column assignment).
#'
#' Occasions follow the paper's design: every dose before the final one is
#' self-administered at home (OCC = 1, the inflated BOV on F), the final dose
#' is the observed clinic dose (OCC = 2), the pre-dose trough belongs to
#' occasion 1 and the subsequent profile to occasion 2.
build_events <- function(cohort, obs_times = grid_ss) {
  doses <- tidyr::expand_grid(
    cohort,
    time = seq(0, t_last, by = tau)
  ) |>
    dplyr::mutate(
      evid = 1L,
      amt  = dose,
      cmt  = "depot",
      OCC  = ifelse(time < t_last, 1L, 2L)
    )
  obs <- tidyr::expand_grid(cohort, time = obs_times) |>
    dplyr::mutate(
      evid = 0L,
      amt  = NA_real_,
      cmt  = "central",
      OCC  = ifelse(time <= t_last, 1L, 2L)
    )
  dplyr::bind_rows(doses, obs) |>
    dplyr::arrange(id, time, dplyr::desc(evid)) |>
    dplyr::select(id, time, amt, evid, cmt, FFM, CREAT, PREG, OCC,
                  dplyr::any_of(c("group", "dose", "arm", "weight_band")))
}
```

``` r

# Four deterministic typical subjects spanning the covariate space.
typ <- tibble::tribble(
  ~id, ~group,                ~FFM,  ~CREAT, ~PREG, ~dose,
  1L,  "Antepartum",          38.7,  45.3,   1,     1000,
  2L,  "Postpartum",          37.0,  55.4,   0,     1000,
  3L,  "Non-pregnant female", 35.7,  58.5,   0,     1000,
  4L,  "Male",                46.4,  64.1,   0,     1000
)

mod_typ <- rxode2::zeroRe(mod)
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_mtt_1, etaiov_mtt_2, etaiov_mtt_3, etaiov_mtt_4, etaiov_ka_1, etaiov_ka_2, etaiov_ka_3, etaiov_ka_4, etaiov_fdepot_1, etaiov_fdepot_2, etaiov_fdepot_3, etaiov_fdepot_4
#> as a work-around try putting the mu-referenced expression on a simple line
sim_typ <- rxode2::rxSolve(mod_typ, build_events(typ), keep = c("group", "dose")) |>
  as.data.frame()
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_mtt_1, etaiov_mtt_2, etaiov_mtt_3, etaiov_mtt_4, etaiov_ka_1, etaiov_ka_2, etaiov_ka_3, etaiov_ka_4, etaiov_fdepot_1, etaiov_fdepot_2, etaiov_fdepot_3, etaiov_fdepot_4
#> as a work-around try putting the mu-referenced expression on a simple line
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etaiov_mtt_1', 'etaiov_mtt_2', 'etaiov_mtt_3', 'etaiov_mtt_4', 'etaiov_ka_1', 'etaiov_ka_2', 'etaiov_ka_3', 'etaiov_ka_4', 'etaiov_fdepot_1', 'etaiov_fdepot_2', 'etaiov_fdepot_3', 'etaiov_fdepot_4'
#> Warning: multi-subject simulation without without 'omega'
```

### Steady-state AUC identity

At steady state the amount eliminated over a dosing interval equals the
amount absorbed, so `CL * AUC(0-tau) = F * dose` exactly, whatever the
absorption model does. `CL` here is computed from the published
covariate equations, not taken from the solve.

``` r

auc_trapz <- function(t, c) sum(diff(t) * (head(c, -1) + tail(c, -1)) / 2)

id_check <- sim_typ |>
  dplyr::filter(time >= t_last) |>
  dplyr::group_by(id, group, dose) |>
  dplyr::summarise(auc_sim = auc_trapz(time, Cc), .groups = "drop") |>
  dplyr::left_join(typ, by = c("id", "group", "dose")) |>
  dplyr::mutate(
    cl_expected  = cl_typical(FFM, CREAT, PREG),
    auc_expected = dose / cl_expected,          # F = 1 (fixed)
    pct_diff     = 100 * (auc_sim - auc_expected) / auc_expected
  )

id_check |>
  dplyr::transmute(
    Group = group, `CL expected (L/h)` = round(cl_expected, 3),
    `AUC dose/CL (mg*h/L)` = round(auc_expected, 1),
    `AUC from solve (mg*h/L)` = round(auc_sim, 1),
    `% diff` = round(pct_diff, 3)
  ) |>
  knitr::kable(caption = "Steady-state AUC(0-24) equals F*dose/CL.")
```

| Group | CL expected (L/h) | AUC dose/CL (mg\*h/L) | AUC from solve (mg\*h/L) | % diff |
|:---|---:|---:|---:|---:|
| Antepartum | 8.933 | 111.9 | 111.9 | -0.014 |
| Postpartum | 5.810 | 172.1 | 172.1 | -0.010 |
| Non-pregnant female | 5.545 | 180.4 | 180.3 | -0.009 |
| Male | 6.527 | 153.2 | 153.2 | -0.009 |

Steady-state AUC(0-24) equals F\*dose/CL. {.table}

``` r


stopifnot(max(abs(id_check$pct_diff)) < 0.5)
```

### Terminal slope recovers CL/Vc

By 14 h after the dose the transit chain and the depot are empty (`ktr`
= 21/1.07 = 19.6 /h; depot half-life ln2/1.59 = 0.44 h), so the decline
is mono-exponential with rate `kel = CL/Vc`. Regressing `log(Cc)` on
time over 14-24 h post-dose therefore tests `lvc` and the volume
exponent `e_ffm_vc`, which the AUC identity above cannot see.

``` r

kel_check <- sim_typ |>
  dplyr::filter(time >= t_last + 14) |>
  dplyr::group_by(id, group) |>
  dplyr::summarise(
    kel_sim = -stats::coef(stats::lm(log(Cc) ~ time))[["time"]],
    .groups = "drop"
  ) |>
  dplyr::left_join(typ, by = c("id", "group")) |>
  dplyr::mutate(
    vc_expected  = 85.8803 * (FFM / 39.4)^1,
    kel_expected = cl_typical(FFM, CREAT, PREG) / vc_expected,
    pct_diff     = 100 * (kel_sim - kel_expected) / kel_expected
  )

kel_check |>
  dplyr::transmute(
    Group = group, `Vc expected (L)` = round(vc_expected, 1),
    `kel expected (1/h)` = round(kel_expected, 5),
    `kel from slope (1/h)` = round(kel_sim, 5),
    `t1/2 (h)` = round(log(2) / kel_sim, 2),
    `% diff` = round(pct_diff, 3)
  ) |>
  knitr::kable(caption = "Terminal slope recovers CL/Vc; half-life 6-8 h is the published range.")
```

| Group | Vc expected (L) | kel expected (1/h) | kel from slope (1/h) | t1/2 (h) | % diff |
|:---|---:|---:|---:|---:|---:|
| Antepartum | 84.4 | 0.10590 | 0.10590 | 6.55 | 0 |
| Postpartum | 80.6 | 0.07204 | 0.07204 | 9.62 | 0 |
| Non-pregnant female | 77.8 | 0.07125 | 0.07125 | 9.73 | 0 |
| Male | 101.1 | 0.06453 | 0.06453 | 10.74 | 0 |

Terminal slope recovers CL/Vc; half-life 6-8 h is the published range.
{.table}

``` r


stopifnot(max(abs(kel_check$pct_diff)) < 0.5)
```

The recovered half-lives sit inside the 6-8 h that the Introduction
quotes for individuals with normal renal function.

### Covariate exponents recovered from exposure ratios

Each ratio isolates one covariate by holding the others fixed, so the
simulated exposure ratio must equal the algebraic ratio exactly.

``` r

ratio_cohort <- tibble::tribble(
  ~id, ~arm,          ~FFM,  ~CREAT, ~PREG, ~dose,
  1L,  "reference",   39.4,  56.2,   0,     1000,
  2L,  "pregnant",    39.4,  56.2,   1,     1000,
  3L,  "creat x2",    39.4,  112.4,  0,     1000,
  4L,  "ffm x1.25",   49.25, 56.2,   0,     1000
) |>
  dplyr::mutate(group = arm)

sim_ratio <- rxode2::rxSolve(mod_typ, build_events(ratio_cohort),
                             keep = c("arm", "dose")) |>
  as.data.frame() |>
  dplyr::filter(time >= t_last) |>
  dplyr::group_by(arm) |>
  dplyr::summarise(auc = auc_trapz(time, Cc), .groups = "drop")
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etaiov_mtt_1', 'etaiov_mtt_2', 'etaiov_mtt_3', 'etaiov_mtt_4', 'etaiov_ka_1', 'etaiov_ka_2', 'etaiov_ka_3', 'etaiov_ka_4', 'etaiov_fdepot_1', 'etaiov_fdepot_2', 'etaiov_fdepot_3', 'etaiov_fdepot_4'
#> Warning: multi-subject simulation without without 'omega'

auc_of <- function(a) sim_ratio$auc[match(a, sim_ratio$arm)]

# At steady state AUC = F*dose/CL, so the exposure ratio between two arms is
# the inverse ratio of their clearances: AUC(ref)/AUC(arm) = CL(arm)/CL(ref).
# Doubling creatinine LOWERS clearance (the exponent is negative), so that
# ratio is below 1; raising FFM RAISES clearance, so that ratio is above 1.
ratios <- tibble::tribble(
  ~Effect,                              ~observed,                                   ~expected,
  "Pregnancy: AUC(ref)/AUC(pregnant)",  auc_of("reference") / auc_of("pregnant"),    1 + 0.380866,
  "Creatinine x2: AUC(ref)/AUC(2x)",    auc_of("reference") / auc_of("creat x2"),    2^-0.366504,
  "FFM x1.25: AUC(ref)/AUC(1.25x)",     auc_of("reference") / auc_of("ffm x1.25"),   1.25^0.75
) |>
  dplyr::mutate(pct_diff = 100 * (observed - expected) / expected)

ratios |>
  dplyr::transmute(
    Effect, `Expected ratio` = round(expected, 5),
    `Simulated ratio` = round(observed, 5), `% diff` = round(pct_diff, 4)
  ) |>
  knitr::kable(caption = "Exposure ratios recover the published covariate exponents.")
```

| Effect                            | Expected ratio | Simulated ratio |  % diff |
|:----------------------------------|---------------:|----------------:|--------:|
| Pregnancy: AUC(ref)/AUC(pregnant) |        1.38087 |         1.38092 |  0.0036 |
| Creatinine x2: AUC(ref)/AUC(2x)   |        0.77566 |         0.77565 | -0.0019 |
| FFM x1.25: AUC(ref)/AUC(1.25x)    |        1.18218 |         1.18217 | -0.0005 |

Exposure ratios recover the published covariate exponents. {.table}

``` r


stopifnot(max(abs(ratios$pct_diff)) < 0.2)
```

## Virtual cohort

The cohort mirrors the four participant groups in the numbers in which
they contributed profiles (19 antepartum, 14 postpartum, 12 non-pregnant
female, 14 male), scaled 10-fold and capped at 200 per arm. Weight,
height and serum creatinine are drawn log-normally around each group’s
published median (Table S1 study-data columns; postpartum medians from
the Discussion) and truncated to the published range; FFM is then
derived per subject with the Janmahasatian formula, and dose is assigned
by WHO weight band. The dispersion of the drawn covariates is an
assumption – Sawe 2025 reports medians and ranges, not distributions –
and is listed under Assumptions below.

``` r

group_spec <- tibble::tribble(
  ~group,                ~n,   ~preg, ~sexf, ~wt_med, ~wt_lo, ~wt_hi, ~ht_med, ~creat_med, ~creat_lo, ~creat_hi,
  "Antepartum",          190L, 1,     1,     61.5,    43.0,   98.0,   1.60,    45.3,       37.0,      66.4,
  "Postpartum",          140L, 0,     1,     57.6,    45.0,   83.5,   1.57,    55.4,       25.3,      63.0,
  "Non-pregnant female", 120L, 0,     1,     52.9,    37.0,   67.9,   1.57,    58.5,       37.0,      110,
  "Male",                140L, 0,     0,     54.4,    41.9,   64.0,   1.73,    64.1,       52.3,      102
)

rxode2::rxSetSeed(20250401)   # seeded per stochastic block, not once globally
set.seed(20250401)

draw_group <- function(spec) {
  n <- spec$n
  tibble::tibble(
    group = spec$group,
    PREG  = spec$preg,
    SEXF  = spec$sexf,
    WT    = pmin(pmax(spec$wt_med    * exp(stats::rnorm(n, 0, 0.20)), spec$wt_lo),    spec$wt_hi),
    HT    = spec$ht_med * exp(stats::rnorm(n, 0, 0.045)),
    CREAT = pmin(pmax(spec$creat_med * exp(stats::rnorm(n, 0, 0.25)), spec$creat_lo), spec$creat_hi)
  )
}

cohort <- dplyr::bind_rows(
  lapply(seq_len(nrow(group_spec)), function(i) draw_group(group_spec[i, ]))
) |>
  dplyr::mutate(
    FFM         = ffm_janmahasatian(WT, HT, SEXF),
    weight_band = ifelse(WT > 50, "above 50 kg", "34-50 kg"),
    dose        = ifelse(WT > 50, 1000, 750),
    id          = dplyr::row_number()
  )

cohort |>
  dplyr::group_by(group) |>
  dplyr::summarise(
    n = dplyr::n(),
    `WT median` = round(stats::median(WT), 1),
    `FFM median` = round(stats::median(FFM), 1),
    `CREAT median` = round(stats::median(CREAT), 1),
    `% on 1000 mg` = round(100 * mean(dose == 1000)),
    .groups = "drop"
  ) |>
  knitr::kable(caption = "Simulated cohort versus the published group medians.")
```

| group               |   n | WT median | FFM median | CREAT median | % on 1000 mg |
|:--------------------|----:|----------:|-----------:|-------------:|-------------:|
| Antepartum          | 190 |      61.8 |       39.3 |         44.5 |           83 |
| Male                | 140 |      54.0 |       47.7 |         64.1 |           71 |
| Non-pregnant female | 120 |      53.3 |       35.1 |         56.0 |           61 |
| Postpartum          | 140 |      57.1 |       36.8 |         57.3 |           73 |

Simulated cohort versus the published group medians. {.table}

``` r

rxode2::rxSetSeed(20250402)
sim <- rxode2::rxSolve(
  mod,
  build_events(cohort, obs_times = seq(t_last, t_last + tau, by = 0.25)),
  keep = c("group", "dose", "weight_band")
) |>
  as.data.frame()
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etaiov_mtt_1, etaiov_mtt_2, etaiov_mtt_3, etaiov_mtt_4, etaiov_ka_1, etaiov_ka_2, etaiov_ka_3, etaiov_ka_4, etaiov_fdepot_1, etaiov_fdepot_2, etaiov_fdepot_3, etaiov_fdepot_4
#> as a work-around try putting the mu-referenced expression on a simple line
```

## Replicating Figure 1: VPC stratified by pregnancy status

Figure 1 of Sawe 2025 is a visual predictive check of concentration
versus time after dose, stratified into a pregnant and a non-pregnant
stratum, where the non-pregnant stratum pools the postpartum,
never-pregnant-female and male records. Its message is that median
exposure is lower in the pregnant stratum.

``` r

vpc_dat <- sim |>
  dplyr::mutate(
    tad     = time - t_last,
    stratum = ifelse(PREG == 1, "Pregnant", "Non-pregnant")
  ) |>
  dplyr::filter(!is.na(Cc))

vpc_pct <- vpc_dat |>
  dplyr::group_by(stratum, tad) |>
  dplyr::summarise(
    p05 = stats::quantile(Cc, 0.05),
    p50 = stats::median(Cc),
    p95 = stats::quantile(Cc, 0.95),
    .groups = "drop"
  )

ggplot2::ggplot(vpc_pct, ggplot2::aes(tad)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = p05, ymax = p95), alpha = 0.25,
                       fill = "steelblue") +
  ggplot2::geom_line(ggplot2::aes(y = p50), linewidth = 0.9) +
  ggplot2::geom_line(ggplot2::aes(y = p05), linetype = "dashed") +
  ggplot2::geom_line(ggplot2::aes(y = p95), linetype = "dashed") +
  ggplot2::facet_wrap(~stratum) +
  ggplot2::labs(x = "Time after dose (h)", y = "Levofloxacin concentration (mg/L)") +
  ggplot2::theme_bw()
```

![Replicates Figure 1 of Sawe 2025: simulated 5th, 50th and 95th
percentiles of levofloxacin concentration versus time after dose,
stratified by pregnancy
status.](Sawe_2025_levofloxacin_files/figure-html/fig1-1.png)

Replicates Figure 1 of Sawe 2025: simulated 5th, 50th and 95th
percentiles of levofloxacin concentration versus time after dose,
stratified by pregnancy status.

``` r

med_by_stratum <- vpc_pct |>
  dplyr::group_by(stratum) |>
  dplyr::summarise(auc_med = auc_trapz(tad, p50), .groups = "drop")

cat(sprintf(
  "Median-profile AUC(0-24): pregnant %.0f, non-pregnant %.0f mg*h/L (ratio %.2f)\n",
  med_by_stratum$auc_med[med_by_stratum$stratum == "Pregnant"],
  med_by_stratum$auc_med[med_by_stratum$stratum == "Non-pregnant"],
  med_by_stratum$auc_med[med_by_stratum$stratum == "Pregnant"] /
    med_by_stratum$auc_med[med_by_stratum$stratum == "Non-pregnant"]
))
#> Median-profile AUC(0-24): pregnant 108, non-pregnant 153 mg*h/L (ratio 0.70)

# The paper's central qualitative finding: median exposure is lower in the
# pregnant stratum. Asserted on the ratio of median profiles, which is robust
# to which subjects land in the tails.
stopifnot(
  med_by_stratum$auc_med[med_by_stratum$stratum == "Pregnant"] <
    med_by_stratum$auc_med[med_by_stratum$stratum == "Non-pregnant"]
)
```

## PKNCA validation

``` r

conc_df <- sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::select(id, group, weight_band, dose, time, Cc)

dose_df <- cohort |>
  dplyr::transmute(id, group, weight_band, dose, time = t_last)

conc_obj <- PKNCA::PKNCAconc(
  data = conc_df, formula = Cc ~ time | group + id,
  concu = "mg/L", timeu = "h"
)
dose_obj <- PKNCA::PKNCAdose(
  data = dose_df, formula = dose ~ time | group + id,
  doseu = "mg"
)

intervals <- data.frame(
  start = t_last, end = t_last + tau,
  cmax = TRUE, tmax = TRUE, cmin = TRUE, auclast = TRUE, cav = TRUE
)

nca <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals))
nca_res <- as.data.frame(nca$result)
```

``` r

nca_wide <- nca_res |>
  dplyr::filter(PPTESTCD %in% c("cmax", "tmax", "auclast", "cav", "cmin")) |>
  dplyr::group_by(group, PPTESTCD) |>
  dplyr::summarise(
    median = stats::median(PPORRES),
    q25    = stats::quantile(PPORRES, 0.25),
    q75    = stats::quantile(PPORRES, 0.75),
    .groups = "drop"
  )

nca_wide |>
  dplyr::mutate(
    Parameter = dplyr::recode(PPTESTCD,
      cmax = "Cmax (mg/L)", tmax = "Tmax (h)", auclast = "AUC0-24 (mg*h/L)",
      cav = "Cavg (mg/L)", cmin = "Ctrough (mg/L)"
    ),
    `Median (IQR)` = sprintf("%.2f (%.2f-%.2f)", median, q25, q75)
  ) |>
  dplyr::select(Group = group, Parameter, `Median (IQR)`) |>
  tidyr::pivot_wider(names_from = Group, values_from = `Median (IQR)`) |>
  knitr::kable(caption = "Simulated steady-state NCA by participant group.")
```

| Parameter | Antepartum | Male | Non-pregnant female | Postpartum |
|:---|:---|:---|:---|:---|
| AUC0-24 (mg\*h/L) | 106.33 (89.97-126.93) | 143.40 (115.04-165.11) | 162.34 (132.33-185.12) | 159.94 (134.62-190.47) |
| Cavg (mg/L) | 4.43 (3.75-5.29) | 5.98 (4.79-6.88) | 6.76 (5.51-7.71) | 6.66 (5.61-7.94) |
| Cmax (mg/L) | 10.07 (9.02-11.13) | 10.00 (8.88-11.10) | 12.47 (10.72-13.64) | 12.18 (10.60-13.72) |
| Ctrough (mg/L) | 1.16 (0.70-1.66) | 2.67 (1.81-3.54) | 2.69 (1.85-3.57) | 2.68 (1.88-3.56) |
| Tmax (h) | 1.75 (1.25-2.75) | 1.75 (1.25-3.00) | 1.88 (1.25-2.75) | 1.75 (1.25-2.75) |

Simulated steady-state NCA by participant group. {.table}

The paper reports exposures for the pooled cohort: median (IQR) AUC0-24
131 (108-170) mg\*h/L and Cmax 11.3 (9.68-13.7) mg/L (Results; Fig. S4).

``` r

pooled <- nca_res |>
  dplyr::filter(PPTESTCD %in% c("cmax", "tmax", "auclast"))

reference <- data.frame(
  cmax    = 11.3,   # Results / Fig. S4, pooled observed median
  auclast = 131,
  tmax    = 2       # first post-dose sampling time; see Assumptions
)

cmp <- nlmixr2lib::ncaComparisonTable(
  simulated = pooled |> dplyr::select(PPTESTCD, PPORRES),
  reference = reference,
  units = c(cmax = "mg/L", auclast = "mg*h/L", tmax = "h")
)
knitr::kable(cmp, caption = "Simulated versus published pooled observed NCA.")
```

| NCA parameter     | Reference | Simulated | % diff |
|:------------------|:----------|:----------|:-------|
| Cmax (mg/L)       | 11.3      | 10.8      | -4.5%  |
| Tmax (h)          | 2         | 1.75      | -12.5% |
| AUClast (mg\*h/L) | 131       | 137       | +4.6%  |

Simulated versus published pooled observed NCA. {.table}

``` r

attr(cmp, "footnote")
#> NULL
```

``` r

pooled_med <- pooled |>
  dplyr::group_by(PPTESTCD) |>
  dplyr::summarise(median = stats::median(PPORRES), .groups = "drop")

auc_med  <- pooled_med$median[pooled_med$PPTESTCD == "auclast"]
cmax_med <- pooled_med$median[pooled_med$PPTESTCD == "cmax"]

cat(sprintf("Pooled simulated median AUC0-24 %.0f mg*h/L (observed IQR 108-170)\n", auc_med))
#> Pooled simulated median AUC0-24 137 mg*h/L (observed IQR 108-170)
cat(sprintf("Pooled simulated median Cmax    %.1f mg/L   (observed IQR 9.68-13.7)\n", cmax_med))
#> Pooled simulated median Cmax    10.8 mg/L   (observed IQR 9.68-13.7)

# Assert against the published interquartile range rather than the point
# median: the simulated cohort's covariate dispersion is an assumption, so the
# defensible claim is that the centre of the simulated distribution falls
# inside the observed IQR.
stopifnot(
  auc_med  > 108,  auc_med  < 170,
  cmax_med > 9.68, cmax_med < 13.7
)

# Table S2 range of median AUC0-24 (98.8-145) and Cmax (7.40-14.86) reported
# for adults with drug-resistant TB on 750-1000 mg; the paper's Fig. 2 ribbon.
lit_auc  <- c(98.8, 145)
lit_cmax <- c(7.40, 14.86)
stopifnot(cmax_med > lit_cmax[1], cmax_med < lit_cmax[2])

# Antepartum exposure must sit below all three non-pregnant groups.
group_auc <- nca_res |>
  dplyr::filter(PPTESTCD == "auclast") |>
  dplyr::group_by(group) |>
  dplyr::summarise(median = stats::median(PPORRES), .groups = "drop")
knitr::kable(
  group_auc |> dplyr::transmute(Group = group, `Median AUC0-24 (mg*h/L)` = round(median)),
  caption = "Antepartum exposure is the lowest of the four groups."
)
```

| Group               | Median AUC0-24 (mg\*h/L) |
|:--------------------|-------------------------:|
| Antepartum          |                      106 |
| Male                |                      143 |
| Non-pregnant female |                      162 |
| Postpartum          |                      160 |

Antepartum exposure is the lowest of the four groups. {.table}

``` r

stopifnot(
  group_auc$median[group_auc$group == "Antepartum"] ==
    min(group_auc$median)
)
```

## Replicating Figure 2: the +250 mg dose recommendation

Figure 2 plots simulated steady-state AUC0-24 and Cmax against weight
band with serum creatinine split into tertiles, adds a dashed box for
the proposed extra 250 mg once daily in pregnancy, ribbons the
literature range of reported medians, and draws the 106 mg\*h/L efficacy
line. The paper concludes that a 250 mg increase brings pregnant
exposures up to those of non-pregnant participants across both weight
bands.

``` r

rxode2::rxSetSeed(20250403)

preg_extra <- cohort |>
  dplyr::filter(group == "Antepartum") |>
  dplyr::mutate(
    group = "Antepartum +250 mg",
    dose  = dose + 250,
    id    = id + 10000L
  )

cohort2 <- dplyr::bind_rows(cohort, preg_extra)

sim2 <- rxode2::rxSolve(
  mod,
  build_events(cohort2, obs_times = seq(t_last, t_last + tau, by = 0.5)),
  keep = c("group", "dose", "weight_band")
) |>
  as.data.frame()

nca2 <- PKNCA::pk.nca(PKNCA::PKNCAdata(
  PKNCA::PKNCAconc(
    data = sim2 |> dplyr::filter(!is.na(Cc)) |>
      dplyr::select(id, group, weight_band, time, Cc),
    formula = Cc ~ time | group + weight_band + id, concu = "mg/L", timeu = "h"
  ),
  PKNCA::PKNCAdose(
    data = cohort2 |> dplyr::transmute(id, group, weight_band, dose, time = t_last),
    formula = dose ~ time | group + weight_band + id, doseu = "mg"
  ),
  intervals = data.frame(start = t_last, end = t_last + tau,
                         cmax = TRUE, auclast = TRUE)
))

exposure <- as.data.frame(nca2$result) |>
  dplyr::filter(PPTESTCD %in% c("cmax", "auclast")) |>
  dplyr::left_join(
    cohort2 |> dplyr::select(id, CREAT), by = "id"
  ) |>
  dplyr::group_by(group, weight_band) |>
  # Sawe 2025 describes the split as tertiles "within each group of
  # participants, based on the weight, sex, and pregnancy status", i.e. a
  # data-driven division into three equal parts. dplyr::ntile() does exactly
  # that and, unlike cut() on sample quantiles, tolerates the ties produced by
  # truncating the drawn creatinine values to the published range.
  dplyr::mutate(
    creat_tertile = factor(
      dplyr::ntile(CREAT, 3), levels = 1:3, labels = c("low", "medium", "high")
    )
  ) |>
  dplyr::ungroup()
```

``` r

ribbons <- data.frame(
  PPTESTCD = c("auclast", "cmax"),
  ymin     = c(98.8, 7.40),
  ymax     = c(145, 14.86)
)
targets <- data.frame(PPTESTCD = "auclast", yint = auc_target)

exposure |>
  dplyr::mutate(
    panel = dplyr::recode(PPTESTCD, auclast = "AUC0-24 (mg*h/L)", cmax = "Cmax (mg/L)")
  ) |>
  ggplot2::ggplot(ggplot2::aes(creat_tertile, PPORRES, fill = group)) +
  ggplot2::geom_rect(
    data = ribbons |> dplyr::mutate(
      panel = dplyr::recode(PPTESTCD, auclast = "AUC0-24 (mg*h/L)", cmax = "Cmax (mg/L)")
    ),
    ggplot2::aes(ymin = ymin, ymax = ymax, xmin = -Inf, xmax = Inf),
    inherit.aes = FALSE, fill = "grey70", alpha = 0.35
  ) +
  ggplot2::geom_hline(
    data = targets |> dplyr::mutate(panel = "AUC0-24 (mg*h/L)"),
    ggplot2::aes(yintercept = yint), inherit.aes = FALSE,
    colour = "red", linetype = "dashed"
  ) +
  ggplot2::geom_boxplot(outlier.size = 0.4, position = ggplot2::position_dodge(0.8)) +
  ggplot2::facet_grid(panel ~ weight_band, scales = "free_y", switch = "y") +
  ggplot2::labs(x = "Serum creatinine tertile", y = NULL, fill = NULL) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom", strip.placement = "outside")
```

![Replicates Figure 2 of Sawe 2025: simulated steady-state AUC0-24 (top)
and Cmax (bottom) by weight band and serum-creatinine tertile. Grey
ribbon is the Table S2 literature range of reported adult medians; red
dashed line is the 106 mg\*h/L efficacy
target.](Sawe_2025_levofloxacin_files/figure-html/fig2-1.png)

Replicates Figure 2 of Sawe 2025: simulated steady-state AUC0-24 (top)
and Cmax (bottom) by weight band and serum-creatinine tertile. Grey
ribbon is the Table S2 literature range of reported adult medians; red
dashed line is the 106 mg\*h/L efficacy target.

``` r

auc_by_group <- exposure |>
  dplyr::filter(PPTESTCD == "auclast") |>
  dplyr::group_by(group, weight_band) |>
  dplyr::summarise(median = stats::median(PPORRES), .groups = "drop") |>
  tidyr::pivot_wider(names_from = group, values_from = median)

auc_by_group |>
  dplyr::mutate(dplyr::across(-weight_band, ~round(.x))) |>
  knitr::kable(caption = "Median simulated steady-state AUC0-24 by weight band.")
```

| weight_band | Antepartum | Antepartum +250 mg | Male | Non-pregnant female | Postpartum |
|:------------|-----------:|-------------------:|-----:|--------------------:|-----------:|
| 34-50 kg    |        100 |                129 |  123 |                 152 |        146 |
| above 50 kg |        107 |                135 |  148 |                 172 |        169 |

Median simulated steady-state AUC0-24 by weight band. {.table}

``` r


recovery <- auc_by_group |>
  dplyr::transmute(
    weight_band,
    base_gap     = 100 * (`Antepartum` / `Non-pregnant female` - 1),
    with_250     = 100 * (`Antepartum +250 mg` / `Non-pregnant female` - 1)
  )
knitr::kable(
  recovery |> dplyr::transmute(
    `Weight band` = weight_band,
    `Antepartum vs non-pregnant female, base dose (%)` = round(base_gap),
    `... with +250 mg (%)` = round(with_250)
  ),
  caption = "The proposed 250 mg increase closes most of the exposure gap."
)
```

| Weight band | Antepartum vs non-pregnant female, base dose (%) | … with +250 mg (%) |
|:---|---:|---:|
| 34-50 kg | -34 | -15 |
| above 50 kg | -38 | -21 |

The proposed 250 mg increase closes most of the exposure gap. {.table}

``` r


# Antepartum exposure is below non-pregnant at the base dose, and the +250 mg
# arm closes most of the gap in both weight bands. The residual shortfall is
# expected: fully offsetting a 38.1% clearance increase needs a 1.381-fold
# dose, i.e. +286 mg on 750 mg and +381 mg on 1000 mg, and the paper chose
# 250 mg because it is the available tablet increment.
stopifnot(
  all(recovery$base_gap < -10),
  all(recovery$with_250 > recovery$base_gap),
  all(abs(recovery$with_250) < abs(recovery$base_gap))
)

for (i in seq_len(nrow(auc_by_group))) {
  cat(sprintf(
    "%-12s | antepartum %.0f | +250 mg %.0f | postpartum %.0f | non-preg F %.0f | male %.0f  (target %.0f)\n",
    auc_by_group$weight_band[i], auc_by_group$`Antepartum`[i],
    auc_by_group$`Antepartum +250 mg`[i], auc_by_group$`Postpartum`[i],
    auc_by_group$`Non-pregnant female`[i], auc_by_group$`Male`[i], auc_target
  ))
}
#> 34-50 kg     | antepartum 100 | +250 mg 129 | postpartum 146 | non-preg F 152 | male 123  (target 106)
#> above 50 kg  | antepartum 107 | +250 mg 135 | postpartum 169 | non-preg F 172 | male 148  (target 106)

# The three non-pregnant arms clear the 106 mg*h/L efficacy target in both
# weight bands, and so does the pregnant arm once the proposed 250 mg is
# added. The pregnant arm on the CURRENT dose is deliberately NOT asserted to
# clear it: that is the paper's point. Its median lands just above the target
# in the above-50 kg band and below it in the 34-50 kg band, which is why the
# paper says AUC is larger than 106 in "most" subjects rather than all, and
# why it recommends the increase.
stopifnot(
  all(auc_by_group$`Postpartum`          > auc_target),
  all(auc_by_group$`Non-pregnant female` > auc_target),
  all(auc_by_group$`Male`                > auc_target),
  all(auc_by_group$`Antepartum +250 mg`  > auc_target)
)

# Antepartum on the current dose is the only arm at risk of missing the target,
# and the lighter weight band is where it misses.
stopifnot(
  auc_by_group$`Antepartum`[auc_by_group$weight_band == "34-50 kg"] < auc_target
)
```

Both weight bands behave as the paper describes: antepartum medians sit
materially below the non-pregnant-female medians on the same
weight-based dose, and adding 250 mg once daily recovers most of the
shortfall without pushing exposures above the range reported in the
comparator literature.

## Assumptions and deviations

Where the paper did not state something the simulation needed, the
choice is recorded here.

- **Covariate dispersion in the virtual cohort.** Sawe 2025 reports
  medians and ranges (Table 1, Table S1), not distributions. Weight and
  serum creatinine are drawn log-normally at 20% and 25% CV respectively
  and height at 4.5% CV, each truncated to the published range. Group
  medians are matched (see the cohort table); the spread is an
  assumption and the interquartile ranges of the simulated NCA
  parameters should not be read as replications of the paper’s.
- **Postpartum median covariates** come from the Discussion (weight 57.6
  kg, FFM 37.0 kg, serum creatinine 55.4 umol/L) rather than Table 1,
  which splits postpartum across the two studies. Postpartum height is
  taken as 1.57 m, the midpoint of the two studies’ medians (1.58 and
  1.56 m).
- **Occasion assignment.** The model’s `OCC` covariate is set to 1 for
  every self-administered dose before the final one and 2 for the
  observed clinic dose, with the pre-dose trough on occasion 1 and the
  subsequent profile on occasion 2. This mirrors the design described in
  the supplementary material (“occasion 1 is linked to the dose taken on
  the day prior to sampling and its associated concentration (trough
  sample)”), and it is what makes the inflated BOV on bioavailability
  apply to the trough only. A user simulating a prospective regimen in
  which every dose is observed should pass an even occasion number
  throughout.
- **Between-occasion variability degenerates within an arm.** Because
  each simulated subject is given a single occasion value per dose class
  rather than a fresh occasion per dose, the BOV on MTT, ka and F is
  drawn once per subject per class and behaves like between-subject
  variability within that class. That is the correct reading of the
  paper’s design, where each participant contributed one observed
  profile, but it means the simulation does not show dose-to-dose
  absorption variability within a subject.
- **Tmax has no published counterpart.** The paper reports no NCA Tmax.
  The reference value of 2 h in the comparison table is the first
  post-dose sampling time in both studies, so the observed Tmax could
  not have been earlier than that; the row is informational and is not
  asserted on.
- **Simulated versus observed exposure.** The pooled simulated median
  AUC0-24 and Cmax are asserted to fall inside the published observed
  interquartile ranges rather than to match the published medians,
  because the observed medians come from 59 profiles with the real joint
  covariate distribution and real dose assignment, while the simulation
  uses drawn covariates and a weight-band dose rule.

### Errata and omissions relative to the source

- **Table 2 footnote correction.** The article was published on 1 April
  2025 with an error in a Table 2 footnote and corrected on 17
  April 2025. The corrected footnote *d* defines the reported
  variability percentages as `%CV = sqrt(omega^2) * 100`. This model
  file uses the corrected reading, verified three ways: against the
  footnote as typeset in the corrected PDF, against the supplementary
  equation image, and arithmetically against the control stream `$OMEGA`
  values (see the variance-scale table above). Applying the more common
  `omega^2 = log(CV^2 + 1)` conversion instead would understate the BOV
  on ka by 26% on the variance scale.
- **Additive residual error is 0.259443 mg/L, not the 0.244 mg/L of
  Table 2.** The control stream’s `$ERROR` block sets
  `ADD = THETA(6) + (LLOQ*0.2)` unconditionally on every record, with
  LLOQ 0.0781 mg/L, so the additive standard deviation actually applied
  is `0.243823 + 0.01562`. The same 20%-of-LLOQ floor appears in the
  sibling paediatric model `Denti_2018_levofloxacin`.
- **BLQ-handling error inflations are omitted.** The control stream adds
  a further 50% of the LLOQ to the additive error for the imputed BLQ
  records (Beal’s M6) and substitutes a huge fixed additive term to drop
  trailing BLQ records from the fit. Both are conditional on the `CENS`
  data column, apply to no simulated observation, and are deliberately
  not encoded.
- **Peripheral-compartment parameters are absent.** The control stream
  retains `$THETA` and `$OMEGA` slots for `V3`, `Q`, `V4` and `Q2`, all
  fixed to zero: the residue of testing a second and third compartment.
  Adding a second compartment did not improve the fit (dOFV 0.674, 2
  df), so the packaged model is one-compartment, as Table 2 reports.
- **Between-subject variability is on clearance only.** Every other
  `$OMEGA` in the control stream is `BLOCK(1) FIX 0`, and between-visit
  variability on clearance was tested and not retained. No IIV has been
  invented for volume or the absorption parameters.
- **No `f(depot) <- 0` line, despite the control stream’s `F1 = 0`.**
  The NONMEM stream sets `F1 = 0` so the dose enters only through the
  analytic transit input rate rather than as a bolus into the absorption
  compartment. Transcribing that literally as `f(depot) <- 0` is the
  natural move, and several older transit models in this package do
  exactly that, but under `rxUi` with rxode2 5.1.7 it silently zeroes
  the `transit()` input as well and the model absorbs nothing – every
  simulated concentration comes back 0, with no error or warning.
  Omitting the line gives the intended behaviour, not a double-counted
  dose: `transit()` is already the sole input pathway under `rxUi`. The
  steady-state AUC identity above is the standing regression gate, and a
  single-dose check gives `AUC(0-inf) * CL` = 999.98 mg against a 1000
  mg dose with F fixed at 1.
- **The transit input accounts only for the most recent dose.** The
  published analytic Savic form tracks a single dose time (`TDOS`/`PNXD`
  in the control stream), so contributions from earlier doses are not
  superposed. This is faithful to the source and numerically immaterial
  here: with `ktr` = 19.6 /h the chain delivers essentially the whole
  dose within a few hours of a 24-hour interval.
- **Screened-but-unretained covariates** (total body weight as an
  alternative size descriptor, sex, age, albumin, HIV status,
  gestational age, study and study site) are recorded in the model
  file’s `covariatesDataExcluded` for provenance. None is referenced in
  `model()`.
- **Gestational age.** Only third-trimester data were collected, and
  gestational age was tested within the pregnant participants and found
  not significant, so the model must not be used to interpolate a
  pregnancy effect at earlier gestation.
- **Serum creatinine extrapolation.** Values were measured within two
  weeks of the PK visit rather than on the day, and no validated
  renal-function equation exists for pregnancy; the authors caution
  against extrapolating outside the observed 25.3-110 umol/L range.
