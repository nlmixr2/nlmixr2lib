# Joint population PK / PANSS PD model for RP5063 (brilaroxazine), a multimodal
# dopamine-serotonin stabilizer, in 175 adults with an acute exacerbation of
# schizophrenia or schizoaffective disorder (Cantillon 2018, Eur J Drug Metab
# Pharmacokinet 43(5):573-585; doi:10.1007/s13318-018-0472-z).

Cantillon_2018_brilaroxazine <- function() {
  description <- paste(
    "Sequential population pharmacokinetic / pharmacodynamic model for RP5063",
    "(brilaroxazine), a multimodal dopamine (D2/3/4) and serotonin",
    "(5-HT1A/2A/2B/2C/7) receptor stabilizer, in 175 adults with an acute",
    "exacerbation of schizophrenia or schizoaffective disorder dosed 15, 30 or",
    "50 mg orally once daily for 28 days in the phase 2 REFRESH trial",
    "(NCT01490086). The PK layer is a one-compartment model with first-order",
    "absorption and an absorption lag time; a two-compartment model was fitted",
    "but abandoned because most parameter correlations exceeded 0.95 and",
    "standard errors could not be computed, and the authors note the initial",
    "distribution phase carries under 10 percent of the total AUC. Body mass",
    "index is the single retained PK covariate, entering the apparent central",
    "volume as a power function centred on the cohort mean of 23.01 kg/m2;",
    "clearance carries no covariate, so BMI does not shift average",
    "steady-state exposure. The PD layer is an Emax model in which total",
    "Positive and Negative Syndrome Scale (PANSS) score is driven by CUMULATIVE",
    "AUC from the first dose rather than by instantaneous concentration - the",
    "cumulative-exposure predictor beat plasma concentration, effect-compartment",
    "concentration and daily average concentration on objective function. A",
    "placebo-effect term was tested and rejected (the objective function",
    "increased), so none is carried here. The single retained PD covariate is a",
    "geographic-site indicator for the Moldova sites, which flips the sign of",
    "Emax from -31.6 to +29.4 PANSS units; the authors attribute this to a",
    "site-level PANSS rating artefact rather than to pharmacology, and report",
    "no site effect on the PK. Random effects on E0 and Emax are ADDITIVE, so",
    "an individual may either improve or worsen with exposure, and the",
    "bookkeeping state auc_central integrates plasma concentration to supply",
    "the exposure driver.",
    sep = " "
  )
  reference <- paste(
    "Cantillon M, Ings R, Prakash A, Bhat L (2018).",
    "A population pharmacokinetic and pharmacodynamic analysis of RP5063",
    "phase 2 study data in patients with schizophrenia or schizoaffective",
    "disorder. European Journal of Drug Metabolism and Pharmacokinetics",
    "43(5):573-585. doi:10.1007/s13318-018-0472-z.",
    sep = " "
  )
  vignette <- "Cantillon_2018_brilaroxazine"

  # Bookkeeping state that integrates Cc so the Emax layer can read cumulative
  # AUC (the paper's Eq. 11 predictor) off the solve. Same idiom as
  # auc_central in Beguin_2024_carboplatin_dog.R and
  # Assmus_2025_benznidazole_qpcr.R; not a biological compartment.
  paper_specific_compartments <- c("auc_central")

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  compartmentData <- list(
    depot = list(analyte = "brilaroxazine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "brilaroxazine", units = "mg", specimen = "plasma", verified = TRUE),
    auc_central = list(analyte = "brilaroxazine", units = "ug*h/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    BMI = list(
      description = "Baseline body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the apparent central volume as a power function centred on",
        "23.01 kg/m2 (Cantillon 2018 Eq. 18, Vj = V (BMIj/23.01)^c1). The",
        "centring constant is the cohort mean, consistent with the paper's",
        "general covariate form Eq. 3 which divides by mean(cov); Table 1",
        "gives per-arm BMI means of 23.2, 22.3 and 23.4 kg/m2, pooling to",
        "roughly 23.0. Time-fixed at baseline. BMI was the ONLY covariate",
        "retained after backward elimination against the Bonferroni-corrected",
        "critical value of 10.86, with an objective-function drop of 26.8",
        "points; it acts on volume only, so it does not change average",
        "steady-state plasma levels."
      ),
      source_name = "BMI"
    ),
    REGION_MOLDOVA = list(
      description = "Study site located in Moldova (1 = Moldova sites, 0 = all other sites)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (study sites in the USA, India, the Philippines or Malaysia)",
      notes = paste(
        "The paper's 'Geographic Area 5' indicator, GEOG5. Multiplies the",
        "typical Emax, Emaxj = Emax (1 + c1 GEOG5) with c1 = -1.93, so the",
        "Moldova sites carry Emax = -31.6 * (1 - 1.93) = +29.4 PANSS units -",
        "i.e. predicted PANSS RISES with cumulative exposure there. Cantillon",
        "2018 Sect. 3.4 and Sect. 4 read this as an artefact of the PANSS",
        "measurements at that single site ('such data would be considered as",
        "an outlier'), not as pharmacology: there was no geographic-site",
        "effect on the pharmacokinetics. Retained because it was the only",
        "covariate surviving backward elimination against the",
        "Bonferroni-corrected critical value of 13.53, dropping the objective",
        "function by 40 points. Set to 0 to obtain the base dose-response",
        "relationship (the solid line of Fig. 6)."
      ),
      source_name = "GEOG 5"
    )
  )

  # Screened in the covariate analysis but NOT retained in either final model.
  # Documented here so the paper's covariate screen is preserved without
  # declaring covariates that model() never references.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Flagged by the generalized-additive-model screen as influencing",
        "Cl/F, and the stepwise addition ranked sex on Vc third, but it",
        "dropped only 3.1 objective-function points against the",
        "Bonferroni-corrected critical value of 10.86 and was eliminated.",
        "No point estimate is reported. Cohort was 80 percent male."
      ),
      source_name = "sex"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      notes = paste(
        "Flagged by the GAM screen as influencing Vc/F; did not survive",
        "stepwise addition/backward elimination. No point estimate reported.",
        "Cohort mean 36 years, inclusion range 18-65 years."
      ),
      source_name = "age"
    ),
    CRCL = list(
      description = "Creatinine clearance by Cockcroft-Gault, used as a surrogate for glomerular filtration rate",
      units = "mL/min",
      type = "continuous",
      notes = paste(
        "Flagged by the GAM screen as influencing Cl/F; did not survive",
        "stepwise addition/backward elimination. No point estimate reported.",
        "RP5063 is eliminated mainly by CYP3A4 (64 percent) and CYP2D6",
        "(17 percent) metabolism, so a renal covariate would not be expected."
      ),
      source_name = "creatinine clearance (Cockcroft Gault)"
    ),
    SMOKER = list(
      description = "Current smoker indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Flagged by the GAM screen as influencing Cl/F; did not survive",
        "stepwise addition/backward elimination. No point estimate reported."
      ),
      source_name = "smoking"
    ),
    CONMED = list(
      description = "Any concomitant medication use indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Carried into the stepwise search alongside the five GAM-screened",
        "covariates; no influence detected. Concomitant drugs were mainly",
        "benzodiazepines (lorazepam 24 patients, zolpidem 22 patients) plus",
        "occasional antihypertensives, antiepileptics and antibiotics",
        "(Supplemental data Appendix D, which was not available when this model was built). Cantillon 2018 Sect. 4",
        "treats the null covariate result as its drug-drug interaction",
        "assessment."
      ),
      source_name = "concomitant drug use"
    ),
    WT = list(
      description = "Body weight at baseline",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Collected as a covariate (Sect. 2.1) but not reported as reaching",
        "the GAM screen; BMI is the size descriptor the final model uses.",
        "No point estimate reported."
      ),
      source_name = "body weight"
    ),
    RACE_ASIAN = list(
      description = "Asian / Indian race indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Race/ethnicity was collected and screened; the GAM screen found no",
        "detectable effect on the empirical Bayes estimates and it was not",
        "retained. Cantillon 2018 Sect. 4 argues the null result is real",
        "rather than a power artefact, citing the phase 1 single-dose study",
        "that found comparable PK between Japanese and Caucasian subjects.",
        "Cohort was 89 percent Asian/Indian. No point estimate reported."
      ),
      source_name = "race/ethnicity"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 175L,
    n_studies = 1L,
    age_range = "18-65 years (inclusion criterion)",
    age_median = "mean 36 years (per-arm means 36, 37, 35)",
    sex_female_pct = 24,
    race_ethnicity = c(Asian = 89, Black = 5, White = 5, Other = 1),
    disease_state = "Acute exacerbation of schizophrenia (96 percent) or schizoaffective disorder (4 percent); mean duration of illness 9 years; mean baseline total PANSS 87.4",
    dose_range = "15, 30 or 50 mg orally once daily for 28 days, dosed after an overnight fast and 1 h before breakfast",
    regions = "USA, India, Philippines, Malaysia, Moldova",
    bmi_range = "per-arm means 23.2 (SD 4.2), 22.3 (SD 4.6) and 23.4 (SD 3.3) kg/m2",
    notes = paste(
      "Baseline characteristics from Cantillon 2018 Table 1 (REFRESH phase 2,",
      "NCT01490086). 234 patients were randomised to RP5063 (15/30/50 mg),",
      "aripiprazole 15 mg or placebo in a 3:3:3:1:2 ratio; only the 175 who",
      "received RP5063 entered the population analysis. Sect. 3.1 states 80",
      "percent male / 20 percent female, but the Table 1 counts give 133 male",
      "and 42 female, i.e. 24 percent female; the Table 1 counts are used",
      "here. Five plasma samples per patient (pre-dose baseline and on days 1,",
      "8, 22 and 28, in four time-blocks) out to 220 h after the last dose,",
      "spanning non-steady-state and steady-state. Total PANSS was collected",
      "pre-dose on day 1 and at least 2 h post-dose on days 4, 8, 15, 22 and 28."
    )
  )

  ini({
    # ---------------- Pharmacokinetics (Cantillon 2018 Table 2) --------------
    # Final one-compartment model with first-order absorption and lag time,
    # after stepwise covariate addition/backward elimination. Cl and Vc are
    # apparent (oral) values, Cl/F and Vc/F; the paper reports no separate F.
    lcl <- log(5.11); label("Apparent oral clearance, Cl/F (L/h)") # Table 2: 5.11 (SE 0.11); bootstrap 5.11 (SE 0.17)
    lvc <- log(328); label("Apparent oral central volume of distribution at BMI 23.01 kg/m2, Vc/F (L)") # Table 2: 328.00 (SE 31.40); bootstrap 329.00 (SE 2.05)
    lka <- log(0.42); label("First-order absorption rate constant, ka (1/h)") # Table 2: 0.42 (SE 0.17); bootstrap 0.45 (SE 0.12)
    ltlag <- log(0.41); label("Absorption lag time, t lag (h)") # Table 2: 0.41 (SE 0.02); bootstrap 0.47 (SE 0.11)

    # Power exponent of BMI on Vc, centred on the cohort mean 23.01 kg/m2.
    e_bmi_vc <- 0.90; label("Power exponent of BMI/23.01 on Vc/F (unitless)") # Table 2 row 'c': 0.90 (SE 0.36); bootstrap 0.84 (SE 0.34); Eq. 18

    # IIV. Table 2 reports a variance-covariance matrix: the rows are labelled
    # Var(.) and Cov(.,.), so the tabulated numbers are VARIANCES on the log
    # scale, not SDs. Corroborated two ways: (a) the covariance row implies a
    # correlation of 0.04/sqrt(0.16*0.28) = 0.19 as variances, against 0.89
    # under an SD reading; (b) the day-averaged empirical-Bayes Cl/F of
    # Sect. 3.3, 5.17 +/- 0.24 L/h over 58 subjects, is a mean +/- standard
    # error, implying an EBE SD of 0.24*sqrt(58) = 1.83 L/h and a CV of 35
    # percent - consistent with the 42 percent CV of Var = 0.16 after
    # shrinkage, but impossible under the 16 percent CV of an SD reading,
    # since EBE spread cannot exceed the true IIV. See vignette Errata for the
    # conflicting base-model percentages quoted in Sect. 3.2.
    etalcl + etalvc ~ c(0.16, 0.04, 0.28) # Table 2, rows 'Var ( g 1 )' = 0.16, 'Cov ( g 1, g 2 )' = 0.04, 'Var ( g 2 )' = 0.28; Eq. 1 makes g1/g2 correlated
    etalka ~ 2.09 # Table 2, row 'Var ( g 3 )' = 2.09 (SE 6.63); Eq. 1 makes g3 independent; CV 266 percent, the paper's poorly-identified parameter

    # No IIV on t lag: Eq. 1 writes t-lagj = tlag with no random effect.

    # Residual error. Sect. 3.2 calls the intra-subject variability a
    # coefficient of variation, so it is proportional; Table 2 reports it as
    # the variance sigma_1^2.
    propSd <- 0.2646; label("Proportional residual error on plasma concentration (fraction)") # Table 2, row 'r 1 2' (sigma_1^2) = 0.07 (SE 0.02); SD = sqrt(0.07) = 0.2646

    # ---------------- Pharmacodynamics (Cantillon 2018 Table 3) -------------
    # Emax model for total PANSS driven by CUMULATIVE AUC (Eq. 12 with the
    # Eq. 11 predictor). E0 and Emax are NOT log-transformed: Eqs. 14-15 give
    # them ADDITIVE random effects, and Emax is negative.
    e0 <- 87.3; label("Baseline total PANSS score at zero cumulative exposure (PANSS units)") # Table 3: 87.3 (SE 0.711); matches the observed mean baseline PANSS of 87.4 (Sect. 3.1)
    emax <- -31.6; label("Maximal change in total PANSS at infinite cumulative exposure (PANSS units)") # Table 3: -31.6 (SE 4.05); negative = improvement
    lauc50 <- log(89.6); label("Cumulative AUC producing half the maximal PANSS change, AUC50 (ug*h/mL)") # Table 3: 89.6 (SE 30.1)

    # Geographic-site effect on Emax, Emaxj = Emax (1 + c1 GEOG5) (Sect. 3.4,
    # Eq. 4 indicator form). -1.93 flips the sign: -31.6 * (1 - 1.93) = +29.4,
    # the positive Emax the paper quotes for that site.
    e_region_moldova_emax <- -1.93; label("Fractional change in Emax at the Moldova study sites (unitless)") # Table 3 row 'c': -1.93 (SE 0.535)

    # PD IIV. Table 3 uses the same Var(.) labelling as Table 2, so these are
    # variances. Independently confirmed on the E0 row: Var = 164 gives an SD
    # of 12.8 PANSS units, against the observed per-arm baseline PANSS SDs of
    # 13.3, 13.4 and 14.9 in Table 1. An SD reading would put 164 PANSS units
    # of spread on a scale bounded at 30-210, which is impossible.
    etae0 ~ 164 # Table 3, row 'Var ( g 1 )' = 164 (SE 632); Eq. 14 E0j = E0 + g1j, additive
    etaemax ~ 464 # Table 3, row 'Var ( g 2 )' = 464 (SE 308); Eq. 15 Emaxj = Emax + g2j, additive; SD 21.5 pu is what lets an individual Emax turn positive
    etalauc50 ~ 0.476 # Table 3, row 'Var ( g 3 )' = 0.476 (SE 0.878); Eq. 16 x50j = x50 exp(g3j), log-normal

    addSd_PANSS <- 6.863; label("Additive residual error on total PANSS (PANSS units)") # Table 3, row 'r 2' (sigma^2) = 47.1 (SE 6.04); SD = sqrt(47.1) = 6.863
  })

  model({
    # ---------------- 1. Individual pharmacokinetic parameters --------------
    cl <- exp(lcl + etalcl)
    # Eq. 18: Vj = V (BMIj / 23.01)^c1. 23.01 kg/m2 is the cohort mean BMI.
    vc <- exp(lvc + etalvc) * (BMI / 23.01)^e_bmi_vc
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)

    kel <- cl / vc

    # ---------------- 2. ODE system -----------------------------------------
    # One compartment, first-order absorption from an oral depot (Sect. 3.2).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and vc in L give central/vc in mg/L = ug/mL, so integrating
    # it accumulates ug*h/mL - the units of the Table 3 AUC50. This is the
    # Eq. 11 predictor, cumulative AUC from the first dose to the current
    # time, NOT a dosing-interval or steady-state AUC.
    d/dt(auc_central) <- central / vc

    # ---------------- 3. Absorption lag -------------------------------------
    alag(depot) <- tlag

    # ---------------- 4. Observations ---------------------------------------
    # Plasma concentration in ng/mL, matching the assay range of 1.00-500
    # ng/mL used for the dependent variable (Sect. 2.1-2.2).
    Cc <- 1000 * central / vc

    # Emax model for total PANSS, Eq. 12: E(x) = E0 + Emax * x / (x + x50),
    # with x the cumulative AUC. Random effects on E0 and Emax are additive
    # (Eqs. 14-15); the one on AUC50 is log-normal (Eq. 16). The site
    # covariate multiplies the TYPICAL Emax, reproducing the paper's quoted
    # average of +29.4 PANSS units at the Moldova sites.
    e0_i <- e0 + etae0
    emax_i <- emax * (1 + e_region_moldova_emax * REGION_MOLDOVA) + etaemax
    auc50_i <- exp(lauc50 + etalauc50)

    PANSS <- e0_i + emax_i * auc_central / (auc_central + auc50_i)

    # No placebo term: Eq. 13's alpha*t was tested and rejected (Sect. 3.4,
    # "Inclusion of a placebo effect did not improve fit"), and no alpha
    # estimate is reported.

    Cc ~ prop(propSd)
    PANSS ~ add(addSd_PANSS)
  })
}
