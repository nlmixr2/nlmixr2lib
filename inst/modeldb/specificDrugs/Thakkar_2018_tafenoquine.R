Thakkar_2018_tafenoquine <- function() {
  description <- "Two-compartment first-order-absorption population PK model with an absorption lag time for oral tafenoquine, pooled across five phase 1 to phase 3 studies in healthy volunteers and Plasmodium vivax malaria patients (Thakkar 2018). Allometric body-weight scaling on CL/F, V2/F, Q/F and V3/F; capsule-versus-tablet formulation effects on relative bioavailability and on the absorption rate constant; and health status (healthy volunteer versus patient) on both apparent volumes of distribution. Interindividual variability is carried on CL/F and V2/F as a correlated block, plus Ka, the lag time, and the residual error magnitude."
  reference <- paste(
    "Thakkar N, Green JA, Koh GCKW, Duparc S, Tenero D, Goyal N (2018).",
    "Population pharmacokinetics of tafenoquine, a novel antimalarial.",
    "Antimicrob Agents Chemother 62(11):e00711-18.",
    "doi:10.1128/AAC.00711-18.",
    sep = " "
  )
  vignette <- "Thakkar_2018_tafenoquine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed (study baseline). Enters CL/F, V2/F, Q/F and V3/F as allometric power terms with exponents fixed at 0.75 (clearances) and 1 (volumes) per Thakkar 2018 Methods 'Covariate analysis': 'Fixed exponents of 0.75 and 1 were applied for the clearance and volume parameters, respectively'. The paper does not state the normalising weight; 70 kg is used here as the conventional allometric reference. The analysis-set median weight is 69.3 kg (Table 2), so the choice shifts the typical clearances by only (69.3/70)^0.75 = 0.993 and the typical volumes by 69.3/70 = 0.990. See vignette Assumptions and deviations.",
      source_name = "WT"
    ),
    FORM_CAPSULE = list(
      description = "Capsule formulation indicator (1 = capsule, 0 = tablet)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (tablet). Thakkar 2018 Table 3 footnote a: 'The tablet formulation was considered the reference, i.e., F1tablet = 1.' The tablet arm therefore carries the structural F = 1 anchor (lfdepot fixed) and the capsule arm carries the estimated relative bioavailability 0.863.",
      notes = "Time-fixed per subject; each study used a single formulation (Table 1: DDI, SIL, DETECTIVE part 2 and GATHER used tablets; TQT and DETECTIVE part 1 used capsules). Acts on two parameters: relative bioavailability (0.863) and the absorption rate constant (0.924), both as multiplicative factors raised to the indicator per Thakkar 2018 equation 2. Capsule was retained on Ka despite backward elimination favouring its removal (dOBJFV 20.4) because 'use of the formulation-specific Ka allowed for a better approximation of the formulation-specific Cmax' (Results, 'Covariate analysis').",
      source_name = "formulation status"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer indicator (1 = healthy volunteer, 0 = P. vivax malaria patient)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (P. vivax malaria patient). The typical CL/F of 2.96 L/h and V2/F of 915 L quoted in the Thakkar 2018 abstract are explicitly 'in P. vivax-infected subjects', so the patient cohort is the reference and the healthy-volunteer ratios multiply onto it.",
      notes = "Time-fixed per subject and determined by study (Table 1): the DDI, SIL and TQT studies enrolled healthy volunteers (193 subjects, 28.6%), the DETECTIVE part 1 and part 2 studies enrolled patients (482 subjects, 71.4%). Acts multiplicatively on both apparent volumes per equation 2: V2/F ratio (healthy/patient) = 1.35 and V3/F ratio (healthy/patient) = 0.347 (Table 3). The authors attribute the higher central volume in healthy volunteers to the absence of background chloroquine and to dehydration in acute malaria, and the higher peripheral volume in patients to infection-related vascular leak (Discussion).",
      source_name = "health status"
    )
  )

  # Screened during covariate model building but not retained in the final
  # model (Thakkar 2018 Results, 'Covariate analysis': 'other covariates, such
  # as age and ethnicity, were also evaluated and did not demonstrate any
  # relevant impact on tafenoquine PK'; Discussion: 'The model-building
  # exercise demonstrated a lack of an effect of demographics, such as age,
  # gender, and ethnicity, on tafenoquine PK'). The paper reports no point
  # estimate for any of them, so none can be encoded.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      notes = "Listed among the covariates evaluated (Methods, 'Covariate analysis'); no relevant impact and not retained. Analysis-set median 35.0 years (range 15.0-79.0; Table 2)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Gender was one of the categorical covariates tested under equation 2 (Methods, 'Covariate analysis'); not retained. Analysis set was 171 female (25.3%) and 504 male (74.7%) (Table 2)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Ethnicity was evaluated and showed no relevant impact (Results, 'Covariate analysis'). Recorded here as the representative member of the race indicator family; the analysis set was 14.1% Caucasian, 18.2% African American, 23.7% Asian, 28.6% American Indian/Alaska native, 0.1% Other and 15.3% Multiple (Table 2). No per-race point estimate is reported."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tafenoquine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tafenoquine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tafenoquine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 675L,
    n_studies = 5L,
    age_range = "15.0-79.0 years",
    age_median = "35.0 years",
    weight_range = "37.2-138 kg",
    weight_median = "69.3 kg",
    sex_female_pct = 25.3,
    race_ethnicity = c(
      Caucasian = 14.1,
      `African American` = 18.2,
      Asian = 23.7,
      `American Indian/Alaska native` = 28.6,
      Multiple = 15.3,
      Other = 0.1
    ),
    disease_state = "Pooled healthy volunteers (193 subjects, 28.6%) and patients with acute Plasmodium vivax malaria (482 subjects, 71.4%)",
    dose_range = "50, 100, 300, 600 and 1,200 mg tafenoquine orally. Single dose in every study except the supratherapeutic 1,200 mg thorough-QTc cohort, which received 400 mg once daily for 3 days. Patients received tafenoquine with background chloroquine.",
    regions = "Multiregional phase 1 to phase 3 programme (ClinicalTrials.gov NCT02184637, NCT02751294, NCT01928914, NCT01376167)",
    notes = "Demographics from Thakkar 2018 Table 2. The parameter-estimation data set comprised 5,286 tafenoquine plasma observations from 675 subjects across 5 studies (Table 1): 200951 (DDI, phase 1, tablet, healthy), 201780 (SIL, phase 1, tablet, healthy), TAF114582 (TQT, phase 1, capsule, healthy), TAF112582 DETECTIVE part 1 (phase 2B, capsule, patients) and TAF112582 DETECTIVE part 2 (phase 3, tablet, patients). A sixth study, TAF116564 (GATHER, phase 3, 166 patients, tablet, 1,001 samples), was deliberately held out of estimation and used only for external validation; its demographics are in Table S1 (median weight 64.8 kg, 31.3% female). Plasma tafenoquine was assayed by LC-MS-MS with an LLQ of 2 ng/mL (0.5 ng/mL in the SIL study); fewer than 4% of observations were below the quantification limit and were excluded. NONMEM 7.3.0, FOCE-I, log-transformed concentrations. Formulation split: 297 subjects (44.0%) tablet, 378 (56.0%) capsule."
  )

  ini({
    # ----------------------------------------------------------------------
    # Structural parameters. Thakkar 2018 Table 3, column 'Final run'.
    # These are the typical values in the reference cells: P. vivax patients
    # (DIS_HEALTHY = 0) receiving the tablet formulation (FORM_CAPSULE = 0),
    # at the 70 kg allometric reference weight.
    # ----------------------------------------------------------------------
    lka <- log(0.252); label("Absorption rate constant, tablet (1/h)") # Table 3 'K a (h-1)' = 0.252 (bootstrap median 0.254, 90% CI 0.226-0.296)
    lcl <- log(2.96); label("Apparent oral clearance CL/F (L/h)") # Table 3 'CL/F (liters/h)' = 2.96 (bootstrap median 2.96, 90% CI 2.87-3.05)
    lvc <- log(915); label("Apparent central volume V2/F (L)") # Table 3 'V 2/F (liters)' = 915 (bootstrap median 913, 90% CI 879-956)
    lq <- log(5.09); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 'Q/F (liters/h)' = 5.09 (bootstrap median 5.10, 90% CI 4.76-5.43)
    lvp <- log(664); label("Apparent peripheral volume V3/F (L)") # Table 3 'V 3/F (liters)' = 664 (bootstrap median 665, 90% CI 634-692)
    ltlag <- log(0.908); label("Absorption lag time (h)") # Table 3 'Absorption lag time (h)' = 0.908 (bootstrap median 0.930, 90% CI 0.904-0.950)

    # Tablet is the structural bioavailability anchor, not an estimate.
    lfdepot <- fixed(log(1)); label("Relative bioavailability of the tablet formulation (unitless)") # Table 3 footnote a: 'The tablet formulation was considered the reference, i.e., F 1tablet = 1.'

    # ----------------------------------------------------------------------
    # Allometric body-weight scaling. Thakkar 2018 Methods, 'Covariate
    # analysis': 'Fixed exponents of 0.75 and 1 were applied for the clearance
    # and volume parameters, respectively'. One exponent is shared by both
    # clearances and one by both volumes, hence the two-parameter shared-
    # exponent naming. Normalised to 70 kg (see covariateData$WT$notes).
    # ----------------------------------------------------------------------
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent on WT/70 shared by CL/F and Q/F (unitless)") # Methods, 'Covariate analysis'
    e_wt_vc_vp <- fixed(1); label("Allometric exponent on WT/70 shared by V2/F and V3/F (unitless)") # Methods, 'Covariate analysis'

    # ----------------------------------------------------------------------
    # Categorical covariate effects. Thakkar 2018 equation 2,
    # P_ij = theta_pop,j * theta_cov^cat, so each value below is the
    # multiplicative factor applied when the indicator is 1 and the reference
    # cell (indicator 0) is left unchanged.
    # ----------------------------------------------------------------------
    e_form_capsule_fdepot <- 0.863; label("Relative bioavailability of capsule versus tablet (unitless)") # Table 3 'Relative bioavailability (capsule)' = 0.863 (bootstrap median 0.866, 90% CI 0.833-0.900)
    e_form_capsule_ka <- 0.924; label("Multiplicative effect of the capsule formulation on Ka (unitless)") # Table 3 'Capsule effect on K a' = 0.924 (bootstrap median 0.914, 90% CI 0.805-1.03)
    e_dis_healthy_vc <- 1.35; label("V2/F ratio, healthy volunteers versus patients (unitless)") # Table 3 'V 2/F ratio (healthy volunteers/patients)' = 1.35 (bootstrap median 1.35, 90% CI 1.30-1.41)
    e_dis_healthy_vp <- 0.347; label("V3/F ratio, healthy volunteers versus patients (unitless)") # Table 3 'V 3/F ratio (healthy volunteers/patients)' = 0.347 (bootstrap median 0.340, 90% CI 0.295-0.396)

    # ----------------------------------------------------------------------
    # Interindividual variability. Table 3 footnote a states the IIV rows are
    # 'expressed as the percent coefficient of variation', so each diagonal is
    # converted to the internal log-scale variance via omega^2 = log(CV^2 + 1):
    #   CL/F 32.1% -> log(1 + 0.321^2) = 0.09807
    #   V2/F 34.4% -> log(1 + 0.344^2) = 0.11184
    #   Ka   40.4% -> log(1 + 0.404^2) = 0.15119
    #   ALAG1 44.3% -> log(1 + 0.443^2) = 0.17919
    #
    # The 'IIV CL-V2 block' row (33.3) is the OMEGA BLOCK(2) off-diagonal,
    # reported as a CORRELATION in percent rather than as a %CV. Reading it
    # literally as a %CV like the diagonals is arithmetically impossible:
    # log(1 + 0.333^2) = 0.10516, which divided by sqrt(0.09807 * 0.11184) =
    # 0.10473 gives a correlation of 1.004 -- outside [-1, 1], so the implied
    # matrix is not positive definite. With rho = 0.333 the covariance is
    # 0.333 * 0.10473 = 0.03488, which is what is used here. See vignette
    # Assumptions and deviations.
    #
    # IIV on Q/F and V3/F was tested and rejected (Methods, 'Population PK
    # model development': 'The interindividual variability (IIV) parameter was
    # evaluated for other population parameters (e.g., Q/F, V3/F) without any
    # significant improvement in model fit or a drop in the objective function
    # value ... Thus, it was not included in the model'), so neither carries an
    # eta here.
    # ----------------------------------------------------------------------
    etalcl + etalvc ~ c(
      0.09807,
      0.03488, 0.11184
    ) # Table 3 rows 'IIV CL/F' = 32.1, 'IIV V 2/F' = 34.4, 'IIV CL-V 2 block' = 33.3 (correlation)
    etalka ~ 0.15119 # Table 3 row 'IIV K a' = 40.4
    etaltlag ~ 0.17919 # Table 3 row 'IIV ALAG1' = 44.3; ETA shrinkage 54.3% (Results, 'Final model')

    # ----------------------------------------------------------------------
    # Residual variability. Methods, 'Population PK model development': 'An
    # additive residual error with IIV was used to describe the residual
    # variability. The additive error with the log-transformed data reflected
    # an exponential residual error model.' Additive-on-log-scale is exactly
    # nlmixr2's lnorm() error structure.
    #
    # The magnitude itself carries an eta: Table 3 reports both a 'Random
    # residual variability (% CV)' of 15.0 and an 'IIV error' of 33.0, and
    # Results ('Final model') quotes a separate 'ETA shrinkage for residual
    # variability' of 8.1%. So the per-subject residual SD is
    # expSd * exp(etaexpSd), applied in model().
    #   sigma  15.0% -> sqrt(log(1 + 0.150^2)) = 0.14917
    #   IIV    33.0% -> log(1 + 0.330^2)       = 0.10337
    # ----------------------------------------------------------------------
    expSd <- 0.14917; label("Residual SD on the natural-log concentration scale (log ng/mL)") # Table 3 'Random residual variability (% CV)' = 15.0 (bootstrap median 14.9, 90% CI 14.3-15.8)
    etaexpSd ~ 0.10337 # Table 3 row 'IIV error' = 33.0; ETA shrinkage 8.1% (Results, 'Final model')
  })

  model({
    # ------------------------------------------------------------------
    # Individual parameters. Allometry uses equation-1 power form on
    # WT/70; the categorical terms use equation-2 form theta_cov^cat, so a
    # subject in the reference cell (tablet, patient) multiplies by 1.
    # ------------------------------------------------------------------
    ka <- exp(lka + etalka) * e_form_capsule_ka^FORM_CAPSULE
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp * e_dis_healthy_vc^DIS_HEALTHY
    q <- exp(lq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp) * (WT / 70)^e_wt_vc_vp * e_dis_healthy_vp^DIS_HEALTHY
    tlag <- exp(ltlag + etaltlag)
    fdepot <- exp(lfdepot) * e_form_capsule_fdepot^FORM_CAPSULE

    # Micro-constants for the two-compartment disposition.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment model with first-order absorption and elimination
    # (Results, 'Population PK model development').
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Relative bioavailability and the absorption lag act on the depot.
    f(depot) <- fdepot
    alag(depot) <- tlag

    # Plasma concentration. Dose is in mg and volumes are in L, so
    # central / vc is mg/L; multiply by 1000 to match the ng/mL used
    # throughout the paper (assay LLQ 2 ng/mL, observed range 2.13 to
    # 1,013 ng/mL per Methods, 'Population PK model development').
    Cc <- central / vc * 1000

    # Per-subject residual SD (see the ini() residual block).
    expSdi <- expSd * exp(etaexpSd)
    Cc ~ lnorm(expSdi)
  })
}
