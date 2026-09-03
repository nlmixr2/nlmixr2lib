Reif_2013_drospirenone <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption and an",
    "absorption lag for drospirenone (DRSP) in young healthy women taking an",
    "ethinylestradiol 20 ug / DRSP 3 mg combined oral contraceptive in",
    "conventional 24/4-day, fixed extended 120/4-day and flexible extended",
    "24-120/4-day regimens. Apparent oral clearance carries a body-weight",
    "effect and is 6.55% lower at the Week 27 sampling occasion than at",
    "Week 3; relative bioavailability carries correlated inter-individual and",
    "inter-occasion variability. Companion model:",
    "modellib('Reif_2013_ethinylestradiol')."
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
    # Doses are entered in ug (the 3 mg tablet is 3000 ug) so that
    # central / vc lands in ug/L, which is identical to the ng/mL the paper
    # reports DRSP serum concentrations and AUC0-24h,ss in.
    dosing = "ug",
    concentration = "ng/mL"
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
        "51-79.8 kg (Table 2, reported for the ethinylestradiol dataset; the",
        "paper states the drospirenone dataset was similar)."
      ),
      source_name        = "BW"
    ),
    OCC = list(
      description        = "Sampling occasion",
      units              = "(count)",
      type               = "categorical",
      reference_category = "1 (Week 3, days 15-21 of the first cycle)",
      notes              = paste(
        "Two occasions: OCC = 1 is the Week 3 visit and OCC = 2 is the",
        "Week 27 visit after about 6 months of treatment. Decomposed inside",
        "model() into the binary indicators oc1 and oc2; oc2 is the paper's",
        "OCA term. It does double duty here: it gates the typical-value step",
        "in CL/F between visits, and it selects which per-occasion",
        "bioavailability IOV eta applies."
      ),
      source_name        = "OCA"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "drospirenone", units = "ug",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "drospirenone", units = "ug",
      specimen = "serum", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "drospirenone", units = "ug",
      specimen = "serum", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1096,
    n_studies      = 1,
    n_observations = 4042,
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
      "The drospirenone PK dataset held 4042 serum concentrations from 1096",
      "subjects (Table 1). Demographics are quoted from Table 2, which",
      "tabulates the ethinylestradiol dataset (n = 1109); the paper states",
      "'similar results were observed in the DRSP PK dataset (data not",
      "shown)'. Sparse sampling: two serum samples 45-120 minutes apart on",
      "each of two occasions (Week 3 and Week 27); 99% of samples fell 0-36 h",
      "after the preceding dose. Records judged not to represent steady-state",
      "concentrations were excluded from the drospirenone dataset."
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Pre-selected for the covariate analysis but not retained in the",
        "final drospirenone model; Table 3 tabulates drospirenone exposure",
        "under the 'Age NA' heading. Age IS retained in the companion",
        "ethinylestradiol model."
      )
    ),
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
        "Orienting covariate analyses indicated alcohol consumption might",
        "affect the central volume V2, but the paper reports the effect on",
        "drospirenone exposure as negligible and it was not retained in the",
        "final model."
      )
    )
  )

  ini({
    # Structural parameters - Appendix 2 (final model for drospirenone).
    # All disposition parameters are apparent (/F) oral values.
    lka <- log(2.18)
    label("First-order absorption rate constant (1/h)")  # Appendix 2: ka = 2.18 1/hour (RSE 8.26%)
    lcl <- log(3.52)
    label("Apparent oral clearance at the Week 3 occasion for a typical subject weighing 62 kg (CL/F, L/h)")  # Appendix 2: TVCL_week3/F = 3.52 l/hour (RSE 0.98%)
    lvc <- log(51.6)
    label("Apparent central volume of distribution (V2/F, L)")  # Appendix 2: V2/F = 51.6 L (RSE 3.64%)
    lvp <- log(204)
    label("Apparent peripheral volume of distribution (V3/F, L)")  # Appendix 2: V3/F = 204 L (RSE 7.99%)
    lq <- log(17.5)
    label("Apparent inter-compartmental clearance (Q3/F, L/h)")  # Appendix 2: Q3/F = 17.5 l/hour (RSE 4.62%)
    ltlag <- log(0.372)
    label("Absorption lag time (h)")  # Appendix 2: ALAG = 0.372 hour (RSE 2.54%)
    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability (unitless)")  # Appendix 2: F = 1 (structural anchor, no RSE reported)

    # Occasion effect on CL/F - Appendix 2 footnote:
    # CL = TVCL_week3 * (1 + TVCL_visit5*OCA) * EXP(ETA(1)) * CO2
    e_occ2_cl <- -0.0655
    label("Fractional change in CL/F at Week 27 relative to Week 3 (unitless)")  # Appendix 2: Diff_CL_week27/F = -6.55 % difference vs Week 3 (RSE 13.0%); 3.52 -> 3.29 L/h

    # Covariate effect on CL/F - Appendix 2 footnote: CO2 = (1 + CL_BW*(BW-62))
    e_wt_cl <- 0.00672
    label("Fractional change in CL/F per kg of body weight away from 62 kg (1/kg)")  # Appendix 2: CL_BW = 0.672 %/kg (RSE 14.5%)

    # IIV on CL/F and F, correlated. Appendix 2 reports 55.4 CV% on CL/F,
    # 47.6 CV% on F and a correlation coefficient of 0.91. The Appendix 2
    # footnote defines CV% = sqrt(exp(OMEGA) - 1) * 100, so
    #   var(etalcl)     = log(0.554^2 + 1) = 0.267670
    #   var(etalfdepot) = log(0.476^2 + 1) = 0.204227
    #   cov             = 0.91 * sqrt(0.267670 * 0.204227) = 0.212764
    etalcl + etalfdepot ~ c(
      0.267670,
      0.212764, 0.204227
    )

    # IOV on F across the two sampling occasions. Appendix 2 reports one
    # shared magnitude of 19.5 CV%, so var = log(0.195^2 + 1) = 0.037320 and
    # the second occasion repeats it (NONMEM $OMEGA BLOCK(1) SAME idiom).
    etaiov_fdepot_1 ~ 0.037320
    etaiov_fdepot_2 ~ fixed(0.037320)

    # Residual error - Appendix 2 reports a 9.61 CV% proportional error and
    # marks it with the same footnote, so sigma = sqrt(log(0.0961^2 + 1))
    # = 0.095879 rather than 0.0961 read off directly. The additive term is
    # reported directly in ug/l (= ng/mL) and carries no such footnote.
    propSd <- 0.095879
    label("Proportional residual error (fraction)")  # Appendix 2: proportional error = 9.61 CV% (RSE 1.88%)
    addSd <- 0.95
    label("Additive residual error (ng/mL)")  # Appendix 2: additive error = 0.95 ug/l (RSE 26.0%)
  })

  model({
    # 1. Occasion indicators. OCC = 1 is Week 3 (the paper's OCA = 0) and
    #    OCC = 2 is Week 27 (OCA = 1).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2

    # 2. Individual parameters, exactly as printed in the Appendix 2 footnote.
    cl <- exp(lcl + etalcl) *
      (1 + e_occ2_cl * oc2) *
      (1 + e_wt_cl * (WT - 62))
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)
    ka <- exp(lka)
    tlag <- exp(ltlag)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central -
      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Bioavailability and lag time
    f(depot) <- exp(lfdepot + etalfdepot + iov_fdepot)
    alag(depot) <- tlag

    # 6. Observation and error
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
