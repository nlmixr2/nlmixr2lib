# Pregnancy and postpartum physiology (Berton 2024)

## Model and source

    #> ℹ parameter labels from comments will be replaced by 'label()'

- Citation: Berton M, Stader F, Bettonte S, Battegay M, Marzolini C.
  Ideal Time to Conduct a Pharmacokinetic Investigation After Delivery
  to Fully Capture the Effect of Pregnancy on Drug Exposure. Open Forum
  Infect Dis. 2024;11(10):ofae585. <doi:10.1093/ofid/ofae585>.
  Equations: Supplementary Table 2 (supplementary data file 1).
- Article: <https://doi.org/10.1093/ofid/ofae585>
- PubMed Central:
  <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11495486/>

Berton 2024 asks when, after delivery, a pharmacokinetic study should be
run so that the postpartum reference fully captures the effect of
pregnancy on drug exposure. To answer it, the authors extended the Basel
in-house Matlab whole-body PBPK framework (Stader 2019) with a pregnant
and postpartum virtual population. They pooled 224 articles, normalised
every anatomical, physiological and enzyme-abundance measurement to its
non-pregnant value, and regressed the fold change on gestational age
(Supplementary Table 2). The population was then verified with seven
antiretrovirals and used to show that drug exposure is mostly back to
its pre-pregnancy level 4 to 5 weeks after delivery.

The packaged model is that **system layer**: the virtual pregnant /
postpartum woman. It has no drug, no dose, no compartments and no ODEs.
Every output is an algebraic function of the `rxode2` time variable,
which the model interprets as **gestational age in weeks**. Delivery is
at 40 weeks for every woman, so `time > 40` is postpartum and
`time - 40` is the week after delivery; the equations were fitted and
the population generated over 0 to 66 weeks.

> PBPK system-parameter (physiology) model of healthy pregnancy and the
> postpartum period (Berton 2024 Open Forum Infect Dis), the pregnant /
> postpartum virtual population of the Basel in-house Matlab whole-body
> PBPK framework. This is the SYSTEM layer only: it has no drug, no
> dosing, no compartments and no ODEs, and every output is an algebraic
> function of the rxode2 time variable, which the model interprets as
> gestational age in WEEKS, with delivery at 40 weeks and the postpartum
> period beyond it (valid 0 to 66 weeks). Each woman’s non-pregnant
> anatomy and physiology come from the framework’s age- and
> height-dependent regressions with their population CVs (Supplementary
> Table 2, citing Stader 2019); the paper’s own gestational-age
> fold-change regressions then scale body weight, adipose / brain /
> heart / kidney / uterus / mammary / plasma / red-cell weights, cardiac
> output, regional blood flows, GFR, haematocrit, albumin, alpha-1-acid
> glycoprotein and the hepatic and intestinal abundance of CYP1A2,
> CYP2B6, CYP2C9, CYP2C19, CYP2D6, CYP3A4, CYP3A5, UGT1A1 and UGT1A4.
> Fetal, placental and amniotic-fluid weights and the fetoplacental
> blood flow are added for the pregnancy. The seven antiretroviral drug
> PBPK models the paper verifies against this population (ritonavir,
> efavirenz, rilpivirine, doravirine, dolutegravir, bictegravir,
> raltegravir) are NOT included: their drug parameters are inherited
> from earlier publications and not printed here, and the way the
> fetoplacental unit is wired into the PBPK ODEs is not described.

### Why the drug layer is not included

The seven drug PBPK models the paper verifies against this population
(ritonavir, efavirenz, rilpivirine, doravirine, dolutegravir,
bictegravir and raltegravir) “were developed previously” in the authors’
obese and ageing papers (Berton 2024 Clin Infect Dis; Stader 2021 Br J
Clin Pharmacol) and “were used without changing any drug parameters”.
None of their drug parameters is printed in Berton 2024, and the paper
does not describe how the fetoplacental unit is connected to the PBPK
ODEs (its partition coefficient, whether it is a perfusion-limited
organ, how it is fed from the uterine flow). The drug layer can be added
once those inputs are available. A bictegravir PBPK on the 2017 version
of the same framework is already in the library
(`Stader_2021_bictegravir_pbpk`); its organ-weight and flow regressions
are an earlier version of the non-pregnant equations used here, so the
two models cannot be chained without re-deriving that PBPK on this
population.

## Population

The development and verification data come from 224 published studies of
pregnant and postpartum women (Supplementary Table 1). The authors
generated 1000 virtual women across 0 to 66 weeks of gestational age to
compare against the verification data set (Figure 1; Supplementary
Figures 1 and 2). The verification trials (Supplementary Table 3)
enrolled pregnant women living with HIV aged 18 to 44 years, sampled at
gestational ages of 13 to 40 weeks and 41 to 80 weeks (postpartum).

The non-pregnant baseline of every woman comes from the framework’s
age-, height- and sex-dependent regressions (Supplementary Table 2,
citing the Stader 2019 ageing repository of White adults), each with a
population CV. Age is the only covariate; sex is fixed to female.

## Source trace

Every gestational-age coefficient is a point estimate of a regression
the authors fitted in Matlab 2020a; **no standard error or confidence
interval is reported**, so every `ini()` entry is `fixed()`. Equations
inherited from other sources (non-pregnant regressions from Stader 2019;
uterus, mammary gland, fetus, placenta and amniotic fluid from
Abduljalil 2012; fetoplacental flow in the Gaohua 2012 form) are written
as constants in `model()`, each with its Supplementary Table 2 row cited
in an in-file comment.

| Output (units) | Supplementary Table 2 rows | Form |
|----|----|----|
| `ht` (cm), `wt_np` (kg) | Body height; Body weight, Non pregnant | age (and height) regressions, CV 3.8% / 15.2% |
| `wt` (kg), `bmi`, `bsa` | Body weight FC; BMI; BSA | quadratic FC in pregnancy and postpartum; DuBois BSA |
| `w_adipose`, `w_brain`, `w_heart`, `w_kidney` (kg) | organ weight + organ weight FC | linear FC, own postpartum line, capped at 1 where printed |
| `w_lung`, `w_bone`, `w_gonads`, `w_skin`, `w_thymus`, `w_gut`, `w_spleen`, `w_pancreas`, `w_liver` (kg) | organ weight; FC = 1 | non-pregnant regression only |
| `w_uterus`, `w_mammary`, `w_muscle` (kg) | Uterus / Mammary gland (FC); Total muscle weight | Abduljalil 2012; both inside muscle |
| `w_plasma`, `w_rbc`, `w_blood` (kg) | Plasma / Red blood cells weight (FC); Blood weight | see Assumptions for the blood-weight row |
| `w_fetus`, `w_placenta`, `w_amniotic`, `w_fetoplacental` (kg) | Foetal / Placental / Amniotic fluid weight | Gompertz, cubic, quintic in GA; 0 after delivery |
| `co` (L/h) | Cardiac output; Cardiac output FC | quadratic FC to 50.9 weeks |
| `fq_*` (% of CO), `q_*` (L/h) | regional flows; Blood flows FC; Brain / Kidney / Uterus blood flow FC; Muscle blood flow; Foetal-placental blood flow | see the flow-budget check below |
| `gfr` (mL/min) | GFR; GFR FC | quadratic FC to 67.5 weeks |
| `hct` (fraction) | HCT; HCT FC | quadratic FC to 65 weeks |
| `hsa`, `agp` (g/L) | Albumin (FC); AAG (FC) | quadratic FC; albumin has its own postpartum line |
| `fc_cyp1a2` … `fc_ugt1a4` | CYP / UGT FC rows | quadratic FC, 1 after the printed knot; CYP3A5 = CYP3A4 |

Gestational-age fold-change coefficients as packaged (all fixed; Berton
2024 Supplementary Table 2). {.table}

| Parameter | Value | Fixed | Label |
|:---|---:|:---|:---|
| bw_preg_c2 | 4.85e-05 | yes | Body weight FC, pregnancy, quadratic GA coefficient (1/week^2) |
| bw_preg_c1 | 4.05e-03 | yes | Body weight FC, pregnancy, linear GA coefficient (1/week) |
| bw_pp_c2 | 3.40e-06 | yes | Body weight FC, postpartum, quadratic GA coefficient (1/week^2) |
| bw_pp_c1 | -1.28e-03 | yes | Body weight FC, postpartum, linear GA coefficient (1/week) |
| bw_pp_c0 | 1.11e+00 | yes | Body weight FC, postpartum, intercept (unitless) |
| adipose_preg_c1 | 6.75e-03 | yes | Adipose weight FC, pregnancy, linear GA coefficient (1/week) |
| adipose_pp_c1 | -1.88e-03 | yes | Adipose weight FC, postpartum, linear GA coefficient (1/week) |
| adipose_pp_c0 | 1.24e+00 | yes | Adipose weight FC, postpartum, intercept (unitless) |
| brain_preg_c1 | -1.49e-03 | yes | Brain weight FC, pregnancy, linear GA coefficient (1/week) |
| brain_pp_c1 | 1.22e-03 | yes | Brain weight FC, postpartum (40 to 42 weeks), linear GA coefficient (1/week) |
| brain_pp_c0 | 9.50e-01 | yes | Brain weight FC, postpartum (40 to 42 weeks), intercept (unitless) |
| heart_preg_c1 | 6.74e-03 | yes | Heart weight FC, pregnancy, linear GA coefficient (1/week) |
| heart_pp_c1 | -1.35e-03 | yes | Heart weight FC, postpartum (40 to 86 weeks), linear GA coefficient (1/week) |
| heart_pp_c0 | 1.12e+00 | yes | Heart weight FC, postpartum (40 to 86 weeks), intercept (unitless) |
| kidney_preg_c1 | 7.32e-03 | yes | Kidney weight FC, pregnancy, linear GA coefficient (1/week) |
| kidney_pp_c1 | -2.93e-02 | yes | Kidney weight FC, postpartum (40 to 50 weeks), linear GA coefficient (1/week) |
| kidney_pp_c0 | 2.46e+00 | yes | Kidney weight FC, postpartum (40 to 50 weeks), intercept (unitless) |
| plasma_preg_c1 | 1.43e-02 | yes | Plasma weight FC, pregnancy, linear GA coefficient (1/week) |
| plasma_pp_c1 | -4.16e-03 | yes | Plasma weight FC, postpartum (40 to 57.7 weeks), linear GA coefficient (1/week) |
| plasma_pp_c0 | 1.24e+00 | yes | Plasma weight FC, postpartum (40 to 57.7 weeks), intercept (unitless) |
| rbc_preg_c2 | 2.47e-04 | yes | Red-blood-cell weight FC, pregnancy, quadratic GA coefficient (1/week^2) |
| rbc_preg_c1 | 2.49e-04 | yes | Red-blood-cell weight FC, pregnancy, linear GA coefficient (1/week) |
| rbc_pp_c2 | 2.07e-03 | yes | Red-blood-cell weight FC, postpartum (40 to 46.5 weeks), quadratic GA coefficient (1/week^2) |
| rbc_pp_c1 | -1.95e-01 | yes | Red-blood-cell weight FC, postpartum (40 to 46.5 weeks), linear GA coefficient (1/week) |
| rbc_pp_c0 | 5.64e+00 | yes | Red-blood-cell weight FC, postpartum (40 to 46.5 weeks), intercept (unitless) |
| co_c2 | -4.77e-04 | yes | Cardiac output FC (0 to 50.9 weeks), quadratic GA coefficient (1/week^2) |
| co_c1 | 2.43e-02 | yes | Cardiac output FC (0 to 50.9 weeks), linear GA coefficient (1/week) |
| qfrac_preg_c1 | -9.90e-03 | yes | FC on the percent-of-cardiac-output flows of organs with unchanged absolute flow, pregnancy, linear GA coefficient (1/week) |
| qfrac_pp_c1 | 7.30e-03 | yes | FC on the percent-of-cardiac-output flows, postpartum (40 to 55 weeks), linear GA coefficient (1/week) |
| qfrac_pp_c0 | 5.95e-01 | yes | FC on the percent-of-cardiac-output flows, postpartum (40 to 55 weeks), intercept (unitless) |
| qbrain_c2 | -1.24e-04 | yes | Absolute brain blood flow FC (0 to 82.4 weeks), quadratic GA coefficient (1/week^2) |
| qbrain_c1 | 1.13e-02 | yes | Absolute brain blood flow FC (0 to 82.4 weeks), linear GA coefficient (1/week) |
| qkidney_c2 | -1.03e-03 | yes | Absolute kidney blood flow FC (0 to 40.7 weeks), quadratic GA coefficient (1/week^2) |
| qkidney_c1 | 4.19e-02 | yes | Absolute kidney blood flow FC (0 to 40.7 weeks), linear GA coefficient (1/week) |
| quterus_preg_c2 | 2.62e-03 | yes | Uterus blood flow FC, pregnancy, quadratic GA coefficient (1/week^2) |
| quterus_preg_c1 | 1.83e-01 | yes | Uterus blood flow FC, pregnancy, linear GA coefficient (1/week) |
| quterus_pp_c1 | -6.60e-01 | yes | Uterus blood flow FC, postpartum (40 to 53 weeks), linear GA coefficient (1/week) |
| quterus_pp_c0 | 3.60e+01 | yes | Uterus blood flow FC, postpartum (40 to 53 weeks), intercept (unitless) |
| gfr_c2 | -2.33e-04 | yes | GFR FC (0 to 67.5 weeks), quadratic GA coefficient (1/week^2) |
| gfr_c1 | 1.57e-02 | yes | GFR FC (0 to 67.5 weeks), linear GA coefficient (1/week) |
| hct_c2 | 8.31e-05 | yes | Haematocrit FC (0 to 65 weeks), quadratic GA coefficient (1/week^2) |
| hct_c1 | -5.37e-03 | yes | Haematocrit FC (0 to 65 weeks), linear GA coefficient (1/week) |
| hsa_preg_c2 | 1.75e-04 | yes | Albumin FC, pregnancy, quadratic GA coefficient (1/week^2) |
| hsa_preg_c1 | -1.32e-02 | yes | Albumin FC, pregnancy, linear GA coefficient (1/week) |
| hsa_pp_c1 | 5.00e-02 | yes | Albumin FC, postpartum (40 to 45 weeks), linear GA coefficient (1/week) |
| hsa_pp_c0 | -1.23e+00 | yes | Albumin FC, postpartum (40 to 45 weeks), intercept (unitless) |
| agp_c2 | 2.07e-04 | yes | Alpha-1-acid glycoprotein FC (0 to 63.7 weeks), quadratic GA coefficient (1/week^2) |
| agp_c1 | -1.32e-02 | yes | Alpha-1-acid glycoprotein FC (0 to 63.7 weeks), linear GA coefficient (1/week) |
| cyp1a2_c2 | 9.54e-04 | yes | CYP1A2 abundance FC (0 to 51.8 weeks), quadratic GA coefficient (1/week^2) |
| cyp1a2_c1 | -4.94e-02 | yes | CYP1A2 abundance FC (0 to 51.8 weeks), linear GA coefficient (1/week) |
| cyp2b6_c2 | -7.36e-04 | yes | CYP2B6 abundance FC (0 to 50.8 weeks), quadratic GA coefficient (1/week^2) |
| cyp2b6_c1 | 3.73e-02 | yes | CYP2B6 abundance FC (0 to 50.8 weeks), linear GA coefficient (1/week) |
| cyp2c9_c2 | -5.26e-03 | yes | CYP2C9 abundance FC (0 to 42 weeks), quadratic GA coefficient (1/week^2) |
| cyp2c9_c1 | 2.18e-01 | yes | CYP2C9 abundance FC (0 to 42 weeks), linear GA coefficient (1/week) |
| cyp2c19_c2 | -1.75e-03 | yes | CYP2C19 abundance FC (0 to 43.5 weeks), quadratic GA coefficient (1/week^2) |
| cyp2c19_c1 | 7.61e-02 | yes | CYP2C19 abundance FC (0 to 43.5 weeks), linear GA coefficient (1/week) |
| cyp2d6_c2 | -1.17e-03 | yes | CYP2D6 abundance FC (0 to 51.4 weeks), quadratic GA coefficient (1/week^2) |
| cyp2d6_c1 | 5.99e-02 | yes | CYP2D6 abundance FC (0 to 51.4 weeks), linear GA coefficient (1/week) |
| cyp3a4_preg_c2 | -1.11e-03 | yes | CYP3A4 abundance FC, pregnancy, quadratic GA coefficient (1/week^2) |
| cyp3a4_preg_c1 | 6.23e-02 | yes | CYP3A4 abundance FC, pregnancy, linear GA coefficient (1/week) |
| cyp3a4_pp_c1 | -1.43e-01 | yes | CYP3A4 abundance FC, postpartum (40 to 45 weeks), linear GA coefficient (1/week) |
| cyp3a4_pp_c0 | 7.44e+00 | yes | CYP3A4 abundance FC, postpartum (40 to 45 weeks), intercept (unitless) |
| ugt1a1_c2 | -4.14e-04 | yes | UGT1A1 abundance FC (0 to 46.5 weeks), quadratic GA coefficient (1/week^2) |
| ugt1a1_c1 | 1.93e-02 | yes | UGT1A1 abundance FC (0 to 46.5 weeks), linear GA coefficient (1/week) |
| ugt1a4_c2 | -2.10e-03 | yes | UGT1A4 abundance FC (0 to 52.1 weeks), quadratic GA coefficient (1/week^2) |
| ugt1a4_c1 | 1.09e-01 | yes | UGT1A4 abundance FC (0 to 52.1 weeks), linear GA coefficient (1/week) |

Virtual-population draws: Mean x (1 + eta), eta ~ N(0, (CV/100)^2), CVs
from Supplementary Table 2. {.table}

| Draw          |  Variance | CV (%) |
|:--------------|----------:|-------:|
| etaht         | 0.0014440 |   3.80 |
| etawt         | 0.0231040 |  15.20 |
| etaw_adipose  | 0.0876160 |  29.60 |
| etaw_bone     | 0.0174240 |  13.20 |
| etaw_brain    | 0.0081000 |   9.00 |
| etaw_gonads   | 0.1211040 |  34.80 |
| etaw_heart    | 0.0412090 |  20.30 |
| etaw_kidney   | 0.0538240 |  23.20 |
| etaw_uterus   | 0.0219040 |  14.80 |
| etaw_mammary  | 0.2667720 |  51.65 |
| etaw_muscle   | 0.0139240 |  11.80 |
| etaw_skin     | 0.0068890 |   8.30 |
| etaw_thymus   | 0.2007040 |  44.80 |
| etaw_gut      | 0.0053290 |   7.30 |
| etaw_spleen   | 0.2672890 |  51.70 |
| etaw_pancreas | 0.0772840 |  27.80 |
| etaw_liver    | 0.0561690 |  23.70 |
| etaw_blood    | 0.0108160 |  10.40 |
| etaw_fetus    | 0.0756800 |  27.51 |
| etaw_placenta | 0.0320410 |  17.90 |
| etaw_amniotic | 0.0895206 |  29.92 |
| etaco         | 0.0445210 |  21.10 |
| etagfr        | 0.0216090 |  14.70 |
| etahct        | 0.0207360 |  14.40 |
| etahsa        | 0.0062410 |   7.90 |
| etaagp        | 0.0590490 |  24.30 |

## Simulation

There is nothing to dose: a solve needs only an observation grid and
each woman’s age.

``` r

mod <- readModelDb("Berton_2024_pregnancy_physiology")
mod_typ <- rxode2::zeroRe(mod)
#> ℹ parameter labels from comments will be replaced by 'label()'
#> Warning: No sigma parameters in the model

ga_grid <- seq(0, 66, by = 0.5)

# Typical 30-year-old woman (all etas zero).
d_typ <- data.frame(id = 1L, time = ga_grid, evid = 0L, AGE = 30)
sim_typ <- rxode2::rxSolve(mod_typ, d_typ, returnType = "data.frame")
#> ℹ omega/sigma items treated as zero: 'etaht', 'etawt', 'etaw_adipose', 'etaw_bone', 'etaw_brain', 'etaw_gonads', 'etaw_heart', 'etaw_kidney', 'etaw_uterus', 'etaw_mammary', 'etaw_muscle', 'etaw_skin', 'etaw_thymus', 'etaw_gut', 'etaw_spleen', 'etaw_pancreas', 'etaw_liver', 'etaw_blood', 'etaw_fetus', 'etaw_placenta', 'etaw_amniotic', 'etaco', 'etagfr', 'etahct', 'etahsa', 'etaagp'
stopifnot(nrow(sim_typ) == length(ga_grid))

# Typical value at given gestational ages, matched by value rather than by row
# position.
typical_at <- function(ga) {
  idx <- match(ga, sim_typ$time)
  stopifnot(!anyNA(idx))
  sim_typ[idx, , drop = FALSE]
}

# Virtual cohort of 200 women (the per-arm cap for these vignettes), ages drawn
# uniformly over the 18 to 44 years of the verification trials, each followed
# from conception to 26 weeks after delivery.
rxode2::rxSetSeed(20241018)
set.seed(20241018)
n_women <- 200L
ages <- data.frame(id = seq_len(n_women), AGE = runif(n_women, 18, 44))
d_pop <- tidyr::expand_grid(id = seq_len(n_women), time = seq(0, 66, by = 1)) |>
  dplyr::left_join(ages, by = "id") |>
  dplyr::mutate(evid = 0L)
sim_pop <- rxode2::rxSolve(mod, d_pop, returnType = "data.frame")
#> ℹ parameter labels from comments will be replaced by 'label()'
stopifnot(dplyr::n_distinct(sim_pop$id) == n_women)
```

## Validation 1: changes at the third trimester

The Results report the third-trimester change relative to non-pregnant
women for several organ weights, cardiac output and the plasma proteins.
The typical woman reproduces them.

``` r

base <- typical_at(0)
term <- typical_at(40)
pct <- function(x) round(100 * (term[[x]] / base[[x]] - 1), 1)

tab_t3 <- data.frame(
  Quantity = c("Adipose weight", "Heart weight", "Kidney weight", "Blood weight",
               "Brain weight", "Cardiac output", "Albumin", "Alpha-1-acid glycoprotein"),
  Output = c("w_adipose", "w_heart", "w_kidney", "w_blood", "w_brain", "co", "hsa", "agp"),
  Published = c(25, 30, 30, 47, -7, 20, -25, -20)
)
tab_t3$Model <- vapply(tab_t3$Output, pct, numeric(1))
tab_t3 |>
  dplyr::rename(
    `Published change, % (Results)` = Published,
    `Model at 40 weeks, %` = Model
  ) |>
  knitr::kable(caption = "Typical 30-year-old woman at 40 weeks vs non-pregnant; published values from the Results text.")
```

| Quantity | Output | Published change, % (Results) | Model at 40 weeks, % |
|:---|:---|---:|---:|
| Adipose weight | w_adipose | 25 | 27.0 |
| Heart weight | w_heart | 30 | 27.0 |
| Kidney weight | w_kidney | 30 | 29.3 |
| Blood weight | w_blood | 47 | 51.0 |
| Brain weight | w_brain | -7 | -6.0 |
| Cardiac output | co | 20 | 20.9 |
| Albumin | hsa | -25 | -24.8 |
| Alpha-1-acid glycoprotein | agp | -20 | -19.7 |

Typical 30-year-old woman at 40 weeks vs non-pregnant; published values
from the Results text. {.table}

``` r


# The published values are rounded figures from the text; every modelled change
# must be within 5 percentage points of them and in the same direction.
stopifnot(
  all(sign(tab_t3$Model) == sign(tab_t3$Published)),
  all(abs(tab_t3$Model - tab_t3$Published) <= 5)
)
```

Cardiac output “returned to baseline 10 weeks after delivery” (Results):

``` r

co_fc <- sim_typ$co / base$co
first_baseline <- min(sim_typ$time[sim_typ$time > 40 & abs(co_fc - 1) < 1e-8])
first_baseline - 40
#> [1] 11
stopifnot(abs((first_baseline - 40) - 10) <= 1)
```

## Validation 2: Table 1, enzyme abundance and albumin

Table 1 gives, for each enzyme and for albumin, the change at 37 to 40
weeks and the postpartum week by which it is back to baseline.

``` r

fc_cols <- c(CYP1A2 = "fc_cyp1a2", CYP2B6 = "fc_cyp2b6", CYP2C9 = "fc_cyp2c9",
             CYP2C19 = "fc_cyp2c19", CYP2D6 = "fc_cyp2d6", CYP3A4 = "fc_cyp3a4",
             UGT1A1 = "fc_ugt1a1", UGT1A4 = "fc_ugt1a4")
sim_typ$fc_hsa_out <- sim_typ$hsa / base$hsa
fc_cols <- c(fc_cols, Albumin = "fc_hsa_out")

t1_change <- c(CYP1A2 = -44, CYP2B6 = 30, CYP2C9 = 56, CYP2C19 = 24, CYP2D6 = 44,
               CYP3A4 = 69, UGT1A1 = 23, UGT1A4 = 116, Albumin = -25)
t1_weeks <- c(CYP1A2 = 12, CYP2B6 = 11, CYP2C9 = 2, CYP2C19 = 4, CYP2D6 = 11,
              CYP3A4 = 5, UGT1A1 = 6, UGT1A4 = 12, Albumin = 5)

w3740 <- sim_typ$time >= 37 & sim_typ$time <= 40
tab_t1 <- data.frame(
  Variable = names(fc_cols),
  pub_change = t1_change[names(fc_cols)],
  mod_mean = vapply(fc_cols, function(x) 100 * (mean(sim_typ[[x]][w3740]) - 1), numeric(1)),
  mod_40 = vapply(fc_cols, function(x) 100 * (sim_typ[[x]][sim_typ$time == 40] - 1), numeric(1)),
  pub_weeks = t1_weeks[names(fc_cols)],
  mod_weeks = vapply(fc_cols, function(x) {
    post <- sim_typ$time > 40
    # First postpartum time from which the fold change stays at 1.
    off <- sim_typ$time[post][abs(sim_typ[[x]][post] - 1) > 1e-8]
    (if (length(off)) max(off) + 0.5 else 40.5) - 40
  }, numeric(1)),
  row.names = NULL
)
tab_t1 |>
  dplyr::mutate(dplyr::across(c(mod_mean, mod_40), \(x) round(x, 1))) |>
  dplyr::rename(
    `Table 1 change, %` = pub_change,
    `Model mean 37-40 wk, %` = mod_mean,
    `Model at 40 wk, %` = mod_40,
    `Table 1 back to baseline, pp wk` = pub_weeks,
    `Model back to baseline, pp wk` = mod_weeks
  ) |>
  knitr::kable(caption = "Table 1 of Berton 2024 against the typical woman.")
```

| Variable | Table 1 change, % | Model mean 37-40 wk, % | Model at 40 wk, % | Table 1 back to baseline, pp wk | Model back to baseline, pp wk |
|:---|---:|---:|---:|---:|---:|
| CYP1A2 | -44 | -48.7 | -45.0 | 12 | 12.0 |
| CYP2B6 | 30 | 34.4 | 31.4 | 11 | 11.0 |
| CYP2C9 | 56 | 59.1 | 30.4 | 2 | 2.5 |
| CYP2C19 | 24 | 33.4 | 24.4 | 4 | 4.0 |
| CYP2D6 | 44 | 57.1 | 52.4 | 11 | 11.5 |
| CYP3A4 | 69 | 75.2 | 71.6 | 5 | 5.5 |
| UGT1A1 | 23 | 12.9 | 11.0 | 6 | 7.0 |
| UGT1A4 | 116 | 108.2 | 100.0 | 12 | 12.5 |
| Albumin | -25 | -24.9 | -24.8 | 5 | 5.5 |

Table 1 of Berton 2024 against the typical woman. {.table}

``` r


stopifnot(
  # Direction of every change matches Table 1.
  all(sign(tab_t1$mod_mean) == sign(tab_t1$pub_change)),
  # The return-to-baseline week follows from the printed knots and agrees with
  # Table 1 to within one week for all nine variables.
  all(abs(tab_t1$mod_weeks - tab_t1$pub_weeks) <= 1)
)
```

The return-to-baseline weeks agree for all nine variables. The size of
the change is less consistent, and Table 1 does not use one rule for it:
the albumin, CYP1A2, CYP2B6, CYP2C19 and CYP3A4 entries match the
equations at 40 weeks to within 3 percentage points, while the CYP2C9
entry (+56%) matches the 37 to 40 week mean (+59%) and is twice the
40-week value. The CYP2D6 (+44%) and UGT1A1 (+23%) entries match
neither; UGT1A1’s +23% is the equation’s peak at about 23 weeks, and
UGT1A4’s +116% is its value at 37 weeks. Table 1 is most likely
summarised from the observed clinical data or the simulated population
rather than read off the fitted curves; the packaged model follows the
Supplementary Table 2 equations.

## Validation 3: the regional blood-flow budget

Supplementary Table 2 gives regional flows as percent of cardiac output.
Most organs are assumed to keep their absolute flow, so a single “Blood
flows FC” scales their percentages while the uterus and the
fetoplacental unit take a growing share. If the fetoplacental flow is
read as percent of cardiac output, as its place under the “% of cardiac
output” heading implies, the listed percentages add to a near-constant
total throughout pregnancy and after delivery: the shared FC was
calibrated to keep the budget. Read as L/h (the unit of the Gaohua 2012
equation it is modelled on) the total would fall from 78% to about 68%
at term.

``` r

budget <- sim_typ |>
  dplyr::mutate(
    total = fq_adipose + fq_bone + fq_gonads + fq_heart + fq_skin + fq_thymus +
      fq_liver + fq_kidney + fq_muscle + fq_fetoplacental,
    total_if_lh = total - fq_fetoplacental + 100 * fq_fetoplacental / co
  )
budget |>
  dplyr::filter(time %in% c(0, 10, 20, 30, 39, 41, 50, 60)) |>
  dplyr::select(time, total, total_if_lh, fq_uterus, fq_fetoplacental) |>
  dplyr::mutate(dplyr::across(-time, \(x) round(x, 1))) |>
  dplyr::rename(
    `GA (wk)` = time,
    `Sum of listed flows, % CO` = total,
    `Same, fetoplacental read as L/h` = total_if_lh,
    `Uterus, % CO` = fq_uterus,
    `Fetoplacental, % CO` = fq_fetoplacental
  ) |>
  knitr::kable(caption = "Percent of cardiac output taken by the organs listed in Supplementary Table 2 (typical woman).")
```

| GA (wk) | Sum of listed flows, % CO | Same, fetoplacental read as L/h | Uterus, % CO | Fetoplacental, % CO |
|---:|---:|---:|---:|---:|
| 0 | 78.4 | 78.4 | 0.9 | 0.0 |
| 10 | 78.1 | 76.8 | 2.9 | 1.8 |
| 20 | 79.6 | 74.4 | 5.4 | 6.6 |
| 30 | 80.9 | 71.4 | 8.4 | 12.3 |
| 39 | 79.6 | 67.5 | 11.5 | 15.9 |
| 41 | 79.7 | 79.7 | 8.5 | 0.0 |
| 50 | 77.9 | 77.9 | 2.8 | 0.0 |
| 60 | 78.4 | 78.4 | 0.9 | 0.0 |

Percent of cardiac output taken by the organs listed in Supplementary
Table 2 (typical woman). {.table}

``` r


stopifnot(
  # Percent reading: the budget stays within 3 percentage points of baseline.
  max(abs(budget$total - budget$total[budget$time == 0])) < 3,
  # L/h reading: the budget would lose more than 8 points near term.
  min(budget$total_if_lh) < budget$total[budget$time == 0] - 8
)
```

## Validation 4: the virtual population against Figure 1 and Supplementary Figure 1

``` r

pop_summ <- function(var) {
  sim_pop |>
    dplyr::group_by(time) |>
    dplyr::summarise(
      p05 = quantile(.data[[var]], 0.05),
      p50 = median(.data[[var]]),
      p95 = quantile(.data[[var]], 0.95),
      .groups = "drop"
    )
}
bw <- pop_summ("wt")
alb <- pop_summ("hsa")
co_p <- pop_summ("co")
qk <- pop_summ("q_kidney")

# Centres of the grey simulated clouds, read by eye. Figure 1A: body weight
# about 65 kg before pregnancy and about 75 kg at 36 to 40 weeks; Figure 1B:
# albumin falls from about 45 to about 34 g/L at term and is back
# near 45 g/L from about 45 weeks; Supplementary Figure 1c: cardiac output about
# 330 L/h before pregnancy and about 420 L/h through the second and third
# trimesters; Supplementary Figure 1d: renal blood flow about 60 L/h before and
# about 85 L/h at mid-pregnancy.
checks <- data.frame(
  Check = c("Body weight median, GA 0 (kg)", "Body weight median, GA 38 (kg)",
            "Albumin median, GA 0 (g/L)", "Albumin median, GA 38 (g/L)",
            "Albumin median, GA 50 (g/L)", "Cardiac output median, GA 0 (L/h)",
            "Cardiac output median, GA 30 (L/h)", "Renal blood flow median, GA 0 (L/h)",
            "Renal blood flow median, GA 20 (L/h)"),
  Model = c(bw$p50[bw$time == 0], bw$p50[bw$time == 38],
            alb$p50[alb$time == 0], alb$p50[alb$time == 38], alb$p50[alb$time == 50],
            co_p$p50[co_p$time == 0], co_p$p50[co_p$time == 30],
            qk$p50[qk$time == 0], qk$p50[qk$time == 20]),
  Figure = c(65, 75, 45, 34, 45, 330, 420, 60, 85)
)
checks$Model <- round(checks$Model, 1)
checks |>
  dplyr::rename(`Model (200 women)` = Model, `Read from figure` = Figure) |>
  knitr::kable(caption = "Cohort medians against values read from the published figures.")
```

| Check                                | Model (200 women) | Read from figure |
|:-------------------------------------|------------------:|-----------------:|
| Body weight median, GA 0 (kg)        |              61.6 |               65 |
| Body weight median, GA 38 (kg)       |              75.4 |               75 |
| Albumin median, GA 0 (g/L)           |              45.2 |               45 |
| Albumin median, GA 38 (g/L)          |              34.0 |               34 |
| Albumin median, GA 50 (g/L)          |              45.2 |               45 |
| Cardiac output median, GA 0 (L/h)    |             329.2 |              330 |
| Cardiac output median, GA 30 (L/h)   |             427.8 |              420 |
| Renal blood flow median, GA 0 (L/h)  |              57.8 |               60 |
| Renal blood flow median, GA 20 (L/h) |              82.5 |               85 |

Cohort medians against values read from the published figures. {.table}

``` r


# Medians are robust to which women land in the tails; 15% covers reading the
# centre of a scatter cloud by eye.
stopifnot(all(abs(checks$Model / checks$Figure - 1) < 0.15))
```

## Replicate published figures

``` r

dplyr::bind_rows(
  dplyr::mutate(bw, var = "Body weight (kg)"),
  dplyr::mutate(alb, var = "Albumin (g/L)")
) |>
  ggplot(aes(time, p50)) +
  geom_ribbon(aes(ymin = p05, ymax = p95), fill = "grey80") +
  geom_line() +
  geom_vline(xintercept = 40, linetype = "dashed") +
  facet_wrap(~var, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, 66, by = 6)) +
  labs(x = "Gestational age (weeks)", y = NULL)
```

![Replicates Figure 1A and 1B of Berton 2024: body weight and albumin
against gestational age (200 virtual women; median and 5th to 95th
percentiles).](Berton_2024_pregnancy_physiology_files/figure-html/fig1-1.png)

Replicates Figure 1A and 1B of Berton 2024: body weight and albumin
against gestational age (200 virtual women; median and 5th to 95th
percentiles).

``` r

sim_typ |>
  dplyr::select(time, dplyr::all_of(unname(fc_cols[names(fc_cols) != "Albumin"]))) |>
  tidyr::pivot_longer(-time, names_to = "enzyme", values_to = "fc") |>
  dplyr::mutate(enzyme = toupper(sub("fc_", "", enzyme))) |>
  ggplot(aes(time, fc)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_vline(xintercept = 40, linetype = "dashed") +
  facet_wrap(~enzyme, ncol = 4) +
  labs(x = "Gestational age (weeks)", y = "Fold change vs non-pregnant")
```

![Replicates the fold-change shapes behind Figure 1C-1D and
Supplementary Figure 2 of Berton 2024: liver enzyme abundance relative
to the non-pregnant woman (typical
values).](Berton_2024_pregnancy_physiology_files/figure-html/fig-enzymes-1.png)

Replicates the fold-change shapes behind Figure 1C-1D and Supplementary
Figure 2 of Berton 2024: liver enzyme abundance relative to the
non-pregnant woman (typical values).

``` r

dplyr::bind_rows(
  dplyr::mutate(co_p, var = "Cardiac output (L/h)"),
  dplyr::mutate(qk, var = "Renal blood flow (L/h)"),
  dplyr::mutate(pop_summ("gfr"), var = "GFR (mL/min)")
) |>
  ggplot(aes(time, p50)) +
  geom_ribbon(aes(ymin = p05, ymax = p95), fill = "grey80") +
  geom_line() +
  geom_vline(xintercept = 40, linetype = "dashed") +
  facet_wrap(~var, scales = "free_y", ncol = 3) +
  labs(x = "Gestational age (weeks)", y = NULL)
```

![Replicates Supplementary Figure 1c-1e of Berton 2024: cardiac output,
renal blood flow and GFR (200 virtual women; median and 5th to 95th
percentiles).](Berton_2024_pregnancy_physiology_files/figure-html/fig-supp1-1.png)

Replicates Supplementary Figure 1c-1e of Berton 2024: cardiac output,
renal blood flow and GFR (200 virtual women; median and 5th to 95th
percentiles).

Figures 2 to 4 and Tables 2 and 3 of the paper are drug-exposure results
of the antiretroviral PBPK models, which are not part of this model.

## Assumptions and deviations

- **Blood weight in pregnancy.** The Supplementary Table 2 row prints
  `Plasma weight * Plasma weight FC / (1 - 0.91*HCT) + RBC weight * RBC weight FC`.
  The first term is already a whole-blood weight, so the red cells are
  counted twice and blood weight more than doubles at term (+105% at 38
  weeks). The Results state that blood weight is 47% higher in the third
  trimester. The model uses the sum of the two components, plasma x its
  FC plus red cells x their FC, which gives +48% at 38 weeks. Dividing
  the pregnant plasma weight by `1 - 0.91 * pregnant HCT` instead gives
  practically the same curve.
- **Fetoplacental flow units.** The row gives no unit. It is listed
  under the “Regional blood flows are expressed as % of cardiac output”
  heading and is read as percent of cardiac output; Validation 3 shows
  that this reading keeps the flow budget constant, which the L/h
  reading of the Gaohua 2012 equation it is modelled on does not. The
  coefficients differ from Gaohua 2012 (refitted). The flow is zero
  before 2.23 weeks, where the cubic is negative, and after delivery.
  The row’s printed range is “\> 2.23 and \< 40”; the model closes it at
  40 weeks like every other pregnancy row, so that the flow and the
  fetoplacental weight end together at delivery.
- **Non-pregnant body surface area.** BSA appears both in the
  non-pregnant regressions (heart, muscle, spleen, liver, blood weights,
  cardiac output, brain and kidney flows, GFR) and as a pregnancy-wide
  row. The model evaluates the non-pregnant regressions with the
  non-pregnant BSA, so that gestational weight gain acts only through
  the printed fold changes, and reports `bsa` at the current body weight
  separately.
- **Distribution of the population CVs.** The paper gives a CV per
  quantity but not the distribution. The model uses the
  `normrnd(Mean, (CV/100) * Mean)` draw of the framework’s published
  population generator (Stader 2019), so `Mean x (1 + eta)` with
  `eta ~ N(0, (CV/100)^2)`. For the widest CVs (mammary gland 51.65%,
  spleen 51.7%, thymus 44.8%) about 2% of draws are negative; the
  framework’s handling of such draws in the 2020a generator is not
  published, and no truncation is applied here. The etas are
  independent.
- **Height and weight are generated, not covariates.** As in the paper’s
  virtual population, height and non-pregnant body weight are drawn from
  their age regressions (CV 3.8% and 15.2%). Age is the only covariate;
  sex is fixed to female (`Sex = 1` in the framework convention).
- **Discontinuities at the knots are as printed.** Several fold changes
  jump at delivery (body weight, adipose, heart, plasma and uterus
  weights, uterine flow, the shared blood-flow FC) or at the week where
  they are capped at 1 (red cells at 46.5 weeks, brain flow at 82.4
  weeks). These reflect the printed equations. Two printed curves do not
  return to 1 inside 0 to 66 weeks (postpartum body weight and adipose
  weight), and the brain-flow FC is still 1.23 at 60 weeks; they are
  used as printed.
- **Age range of the virtual cohort.** The paper draws age from a
  uniform distribution without stating the range for the development
  population; this vignette uses 18 to 44 years, the range of the
  verification trials (Supplementary Table 3).
- **Absolute enzyme abundances are not part of the model.** The paper
  reports only fold changes; the non-pregnant abundances in Figure 1C-1D
  come from the framework and are not printed. The fold changes apply
  equally to the liver and to the intestinal enzymes.
- **Organs without a printed row.** Lymph node and remaining-tissue
  weights and flows, lymph flow, and the hepatic-artery / portal split
  are part of the framework but have no Supplementary Table 2 row and
  are not modelled.
- **No PKNCA validation.** There is no drug and no concentration-time
  profile, so NCA does not apply. Validation uses the published
  third-trimester changes, Table 1, the blood-flow budget and the
  published figures instead.
