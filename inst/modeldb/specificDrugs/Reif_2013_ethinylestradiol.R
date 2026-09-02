Reif_2013_ethinylestradiol <- function() {
  description <- paste(
    "Three-compartment population PK model with first-order absorption and",
    "an absorption lag for ethinylestradiol (EE) in young healthy women",
    "taking an EE 20 ug / drospirenone 3 mg combined oral contraceptive in",
    "conventional 24/4-day, fixed extended 120/4-day and flexible extended",
    "24-120/4-day regimens. Apparent oral clearance carries body-weight and",
    "log-age effects; relative bioavailability is 8.15% higher at the Week 27",
    "sampling occasion than at Week 3. Companion model:",
    "modellib('Reif_2013_drospirenone')."
  )
  reference <- paste(
    "Reif S, Snelder N, Blode H. Characterisation of the pharmacokinetics of",
    "ethinylestradiol and drospirenone in extended-cycle regimens: population",
    "pharmacokinetic analysis from a randomised Phase III study.",
    "J Fam Plann Reprod Health Care. 2013;39(2):e1.",
    "doi:10.1136/jfprhc-2012-100397"
  )
  vignette <- "Reif_2013_ethinylestradiol_drospirenone"
  units <- list(
    time = "h",
    # Doses are entered in ng (the 20 ug tablet is 20000 ng) so that
    # central / vc lands in ng/L, which is identical to the pg/mL the paper
    # reports EE serum concentrations and AUC0-24h,ss in.
    dosing = "ng",
    concentration = "pg/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters CL/F as a linear fractional deviation from the cohort median",
        "of 62 kg: CO2 = 1 + CL_BW * (BW - 62). Recorded at both the Week 3",
        "and Week 27 visits; the paper reports only minor within-subject",
        "change between visits. Cohort median 62 kg, 5th-95th percentile",
        "51-79.8 kg (Table 2)."
      ),
      source_name        = "BW"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters CL/F on the natural-log scale as a linear fractional",
        "deviation from the cohort median of 24 years:",
        "CO1 = 1 + CL_AGE * (LOG(AGE) - LOG(24)). Cohort median 24 years,",
        "5th-95th percentile 19-34 years (Table 2)."
      ),
      source_name        = "AGE"
    ),
    OCC = list(
      description        = "Sampling occasion",
      units              = "(count)",
      type               = "categorical",
      reference_category = "1 (Week 3, days 15-21 of the first cycle)",
      notes              = paste(
        "Two occasions: OCC = 1 is the Week 3 visit and OCC = 2 is the",
        "Week 27 visit after about 6 months of treatment. Decomposed inside",
        "model() into the binary indicator oc2, which is the paper's OCA",
        "term (Appendix 1 footnote: 'week3: OCA = 0, other visits: OCA = 1').",
        "It gates a typical-value step in relative bioavailability, not an",
        "inter-occasion random effect; the EE model has no IOV."
      ),
      source_name        = "OCA"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "ethinylestradiol", units = "ng",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "ethinylestradiol", units = "ng",
      specimen = "serum", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "ethinylestradiol", units = "ng",
      specimen = "serum", verified = TRUE
    ),
    peripheral2 = list(
      analyte = "ethinylestradiol", units = "ng",
      specimen = "serum", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1109,
    n_studies      = 1,
    n_observations = 4218,
    age_range      = "19-34 years (5th-95th percentile)",
    age_median     = "24 years",
    weight_range   = "51-79.8 kg (5th-95th percentile)",
    weight_median  = "62 kg",
    bmi_median     = "22 kg/m^2",
    sex_female_pct = 100,
    race_ethnicity = c(Caucasian = 98),
    disease_state  = "healthy young women using combined oral contraception",
    dose_range     = paste(
      "ethinylestradiol 20 ug / drospirenone 3 mg once daily by mouth;",
      "conventional 24/4-day, fixed extended 120/4-day, or flexible extended",
      "24-120/4-day regimens over 1 year"
    ),
    regions        = "multicentre (study 308683 / NCT00266032)",
    notes          = paste(
      "Baseline demographics are Table 2 of Reif 2013 (the ethinylestradiol",
      "PK dataset, n = 1109 of the 1134 subjects who received study drug).",
      "Sparse sampling: two serum samples 45-120 minutes apart on each of two",
      "occasions (Week 3 and Week 27); 99% of samples fell 0-36 h after the",
      "preceding dose. More than 98% of subjects were Caucasian, so ethnic",
      "group was not tested as a covariate. Smoking and concomitant",
      "medication occurred in fewer than 10% of subjects and were likewise",
      "not tested."
    )
  )

  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Pre-selected for the covariate analysis but not retained in the",
        "final model; strongly correlated with body weight, which was",
        "retained instead."
      )
    ),
    ALCOHOL = list(
      description = "Alcohol consumption category",
      units       = "(category)",
      type        = "categorical",
      notes       = paste(
        "Pre-selected for the covariate analysis (never 12.1%, seldom 51.7%,",
        "occasionally 30.2%, regularly 6.0%; Table 2) but not retained in the",
        "final ethinylestradiol model."
      )
    )
  )

  ini({
    # Structural parameters - Appendix 1 (final model for ethinylestradiol).
    # All disposition parameters are apparent (/F) oral values.
    lka <- log(0.295)
    label("First-order absorption rate constant (1/h)")  # Appendix 1: ka = 0.295 1/hour (RSE 6.98%)
    lcl <- log(25.3)
    label("Apparent oral clearance for a typical subject aged 24 years weighing 62 kg (CL/F, L/h)")  # Appendix 1: TVCL/F = 25.3 l/hour (RSE 1.24%)
    lvc <- log(23.9)
    label("Apparent central volume of distribution (V2/F, L)")  # Appendix 1: V2/F = 23.9 l (RSE 13.6%)
    lvp <- log(1330)
    label("Apparent deep peripheral volume of distribution (V3/F, L)")  # Appendix 1: V3/F = 1,330 l (RSE 3.62%)
    lq <- log(52.9)
    label("Apparent inter-compartmental clearance to the deep peripheral compartment (Q3/F, L/h)")  # Appendix 1: Q3/F = 52.9 l/hour (RSE 7.01%)
    lq2 <- log(8.49)
    label("Apparent inter-compartmental clearance to the shallow peripheral compartment (Q4/F, L/h)")  # Appendix 1: Q4/F = 8.49 l/hour (RSE 34.3%)
    ltlag <- log(0.353)
    label("Absorption lag time (h)")  # Appendix 1: ALAG = 0.353 hour (RSE 2.78%)
    # V4/F is not a separate ini() parameter: Appendix 1 footnote reads 'For
    # the model it was assumed that V2F = V4F', i.e. one estimated volume is
    # used for both. It is set to vc inside model().
    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability at the Week 3 occasion (unitless)")  # Appendix 1: F_week3 = 1 (structural anchor, no RSE reported)

    # Occasion effect on relative bioavailability - Appendix 1 footnote:
    # F = F_week3 * (1 + F_week27 * OCA); (week3: OCA = 0, other visits: OCA = 1)
    e_occ2_fdepot <- 0.0815
    label("Fractional change in relative oral bioavailability at Week 27 relative to Week 3 (unitless)")  # Appendix 1: Diff_F_week27 = 8.15 % difference from Week 3 (RSE 11.0%)

    # Covariate effects on CL/F - Appendix 1 footnote:
    # CL = TVCL*EXP(ETA_CL)*CO2*CO1;
    # CO1 = (1 + CL_AGE * (LOG(AGE)-LOG(24))), CO2 = (1 + CL_BW*(BW-62))
    e_wt_cl <- 0.00591
    label("Fractional change in CL/F per kg of body weight away from 62 kg (1/kg)")  # Appendix 1: CL_BW = 0.591 %/kg (RSE 20.1%)
    e_age_cl <- 0.208
    label("Fractional change in CL/F per unit of ln(age in years) away from ln(24) (unitless)")  # Appendix 1: CL_AGE = 20.8 %/LN(year) (RSE 29.1%)

    # IIV - Appendix 1 reports 33.4 CV% on clearance. Appendix 1 footnote
    # defines CV% = sqrt(exp(OMEGA) - 1) * 100, so
    # omega^2 = log(CV^2 + 1) = log(0.334^2 + 1) = 0.105761.
    etalcl ~ 0.105761

    # Residual error - Appendix 1 reports a 24.4 CV% proportional error and
    # marks it with the same footnote, so sigma = sqrt(log(0.244^2 + 1))
    # = 0.240481 rather than 0.244 read off directly.
    propSd <- 0.240481
    label("Proportional residual error (fraction)")  # Appendix 1: proportional error = 24.4 CV% (RSE 1.38%)
  })

  model({
    # 1. Occasion indicator. OCC = 1 is Week 3 (the paper's OCA = 0) and
    #    OCC = 2 is Week 27 (OCA = 1).
    oc2 <- (OCC == 2)

    # 2. Individual parameters. The covariate model is multiplicative in two
    #    linear deviation terms, exactly as printed in the Appendix 1 footnote.
    cl <- exp(lcl + etalcl) *
      (1 + e_wt_cl * (WT - 62)) *
      (1 + e_age_cl * (log(AGE) - log(24)))
    vc <- exp(lvc)
    vp <- exp(lvp)
    # V4/F assumed equal to V2/F (Appendix 1 footnote).
    vp2 <- vc
    q <- exp(lq)
    q2 <- exp(lq2)
    ka <- exp(lka)
    tlag <- exp(ltlag)

    # 3. Micro-constants. peripheral1 is the paper's deep compartment (V3/F,
    #    Q3/F); peripheral2 is the paper's shallow compartment (V4/F, Q4/F).
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 4. ODE system
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # 5. Bioavailability and lag time
    f(depot) <- exp(lfdepot) * (1 + e_occ2_fdepot * oc2)
    alag(depot) <- tlag

    # 6. Observation and error
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
