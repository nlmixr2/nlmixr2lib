Srimani_2022_ixazomib_diarrhea <- function() {
  description <- "Discrete-time Markov model of weekly diarrhea grade (0-3) for oral ixazomib plus lenalidomide-dexamethasone (LenDex) in relapsed/refractory multiple myeloma from the phase III TOURMALINE-MM1 trial (Srimani 2022). Ixazomib plasma concentrations come from the three-compartment population PK model of Gupta 2017 (fixed) and are integrated to the weekly AUC. For each current grade, the next-week grade follows a proportional-odds cumulative-logit model with a per-grade random intercept. The weekly ixazomib AUC raises all transitions out of grade 0; a first-week onset term raises the grade-0-to-1 transition; a slowly rising time effect makes recovery from grade 1 to grade 0 less likely; IMiD-naive patients leave grade 3 faster. The 16 transition probabilities are model outputs; the Markov chain is advanced week by week outside the solve, as in the vignette."
  reference <- paste(
    "Srimani JK, Diderichsen PM, Hanley MJ, Venkatakrishnan K, Labotka R, Gupta N.",
    "Population pharmacokinetic/pharmacodynamic joint modeling of ixazomib efficacy",
    "and safety using data from the pivotal phase III TOURMALINE-MM1 study in",
    "multiple myeloma patients. CPT Pharmacometrics Syst Pharmacol.",
    "2022;11(8):1085-1099. doi:10.1002/psp4.12815.",
    "Ixazomib PK layer: Gupta N, Diderichsen PM, Hanley MJ, et al. Clin Pharmacokinet.",
    "2017;56(11):1355-1368. doi:10.1007/s40262-017-0526-4",
    "(also available as modellib('Gupta_2017_ixazomib'))."
  )
  vignette <- "Srimani_2022_ixazomib"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL",
    outcome = "weekly diarrhea grade transition probabilities (fraction)"
  )

  compartmentData <- list(
    depot = list(analyte = "ixazomib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ixazomib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ixazomib", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "ixazomib", units = "mg", specimen = "tissue", verified = TRUE),
    auc_central = list(analyte = "ixazomib", units = "ng*h/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area at baseline.",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only by the Gupta 2017 ixazomib PK layer: power covariate on the second peripheral volume, reference 1.87 m^2, exponent 2.06.",
      source_name = "BSA"
    ),
    PRIOR_IMID = list(
      description = "Prior immunomodulatory-drug (IMiD: thalidomide, lenalidomide, pomalidomide) therapy indicator (1 = exposed, 0 = IMiD-naive).",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (IMiD-exposed); the paper's effect is on the naive group.",
      notes = "Additive shift of -2.87 on the cumulative logits out of grade 3 for IMiD-naive patients (Table 3 row 'B3x (PIMID)'), so IMiD-naive patients leave grade 3 diarrhea faster. The printed transition matrix (Supplementary Equation S3) is the IMiD-exposed one. 55.0% of the safety dataset was IMiD-exposed (Supplementary Table 4).",
      source_name = "PIMID"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 720L,
    n_studies = 1L,
    age_range = "30-91 years (median 66)",
    age_median = "66 years",
    sex_female_pct = 43.3,
    race_ethnicity = c(White = 85.0, Black = 1.8, Asian = 8.9, `Not reported` = 2.6, Other = 1.7),
    disease_state = "Relapsed and/or refractory multiple myeloma after 1-3 prior lines of therapy.",
    dose_range = "Ixazomib 4 mg or matching placebo orally on days 1, 8 and 15 of 28-day cycles, both arms with lenalidomide 25 mg (10 mg for reduced creatinine clearance) on days 1-21 and dexamethasone 40 mg on days 1, 8, 15 and 22.",
    regions = "Global (TOURMALINE-MM1, C16010).",
    biomarkers = "Worst diarrhea grade: 0 in 58.2%, 1 in 22.8%, 2 in 14.6%, 3 in 4.4%; grade 1 at baseline in 1.2% (Supplementary Table 4).",
    notes = "Safety population: all 720 treated patients (361 ixazomib, 359 placebo). Adverse events were recorded from the first dose to 30 days after the last dose. Demographics from Srimani 2022 Supplementary Table 4."
  )

  ini({
    # Ixazomib PK -- Gupta 2017 Table 3, fixed (Srimani 2022 Methods).
    lka <- fixed(log(0.34)); label("Ixazomib first-order absorption rate constant (1/h)") # Gupta 2017 Table 3: Ka = 0.34/h
    lcl <- fixed(log(1.86)); label("Ixazomib clearance (L/h)") # Gupta 2017 Table 3: CL = 1.86 L/h
    lvc <- fixed(log(13.7)); label("Ixazomib central volume (L)") # Gupta 2017 Table 3: V2 = 13.7 L
    lfdepot <- fixed(log(0.58)); label("Ixazomib oral bioavailability (fraction)") # Gupta 2017 Table 3: F = 58%
    lq <- fixed(log(5.18)); label("Ixazomib intercompartmental clearance to peripheral1 (L/h)") # Gupta 2017 Table 3: Q3 = 5.18 L/h
    lvp <- fixed(log(309)); label("Ixazomib first peripheral volume (L)") # Gupta 2017 Table 3: V3 = 309 L
    lq2 <- fixed(log(26.1)); label("Ixazomib intercompartmental clearance to peripheral2 (L/h)") # Gupta 2017 Table 3: Q4 = 26.1 L/h
    lvp2 <- fixed(log(205)); label("Ixazomib second peripheral volume (L)") # Gupta 2017 Table 3: V4 = 205 L
    ltlag <- fixed(log(13 / 60)); label("Ixazomib absorption lag time (h)") # Gupta 2017 Table 3: TLAG = 13 min
    e_bsa_vp2 <- fixed(2.06); label("Power exponent of BSA on the second peripheral volume (unitless)") # Gupta 2017 Table 3: V4[BSA] = 2.06
    etalcl + etalfdepot ~ fixed(c(0.17697, 0.22550, 0.42726)) # Gupta 2017 Table 3: CL 44%CV, F 73%CV, rho 0.82
    etalvp2 ~ fixed(0.48492) # Gupta 2017 Table 3: V4 79%CV

    # Cumulative-logit transition parameters -- Srimani 2022 Table 3.
    # b<i>1 is the logit of moving from grade i to grade >= 1; b<i>2 and b<i>3
    # are log-decrements, so logit(>= 2) = b<i>1 - exp(b<i>2) and
    # logit(>= 3) = logit(>= 2) - exp(b<i>3) (Equation 14). With no drug,
    # time or covariate terms these reproduce Supplementary Equation S3.
    b01 <- -5.28; label("Logit of moving from diarrhea grade 0 to grade >= 1") # Table 3 B01 = -5.28
    b02 <- 0.211; label("Log decrement from the grade >= 1 to the grade >= 2 logit, from grade 0") # Table 3 B02 = 0.211
    b03 <- 0.657; label("Log decrement from the grade >= 2 to the grade >= 3 logit, from grade 0") # Table 3 B03 = 0.657
    b11 <- 0.677; label("Logit of moving from diarrhea grade 1 to grade >= 1") # Table 3 B11 = 0.677
    b12 <- 1.99; label("Log decrement from the grade >= 1 to the grade >= 2 logit, from grade 1") # Table 3 B12 = 1.99
    b13 <- 0.541; label("Log decrement from the grade >= 2 to the grade >= 3 logit, from grade 1") # Table 3 B13 = 0.541
    b21 <- 1.94; label("Logit of moving from diarrhea grade 2 to grade >= 1") # Table 3 B21 = 1.94
    b22 <- -1.92; label("Log decrement from the grade >= 1 to the grade >= 2 logit, from grade 2") # Table 3 B22 = -1.92
    b23 <- 2.21; label("Log decrement from the grade >= 2 to the grade >= 3 logit, from grade 2") # Table 3 B23 = 2.21
    b31 <- 3.63; label("Logit of moving from diarrhea grade 3 to grade >= 1 (IMiD-exposed)") # Table 3 B31 = 3.63
    b32 <- -2.15; label("Log decrement from the grade >= 1 to the grade >= 2 logit, from grade 3") # Table 3 B32 = -2.15
    b33 <- -1.78; label("Log decrement from the grade >= 2 to the grade >= 3 logit, from grade 3") # Table 3 B33 = -1.78

    # Explicit time effects. Table 3 prints the slow
    # time-effect magnitude on the log scale (exp(2.25) = 9.49 log OR) and
    # the rate constant in log 1/h (exp(-10.4) * 168 = 0.0051/week).
    lemax_time_g10 <- 2.25; label("Log maximal slow time effect slowing recovery from grade 1 to grade 0 (log of log-odds units)") # Table 3 P1|0T = 2.25 (9.52 log OR)
    lk_time_g10 <- -10.4; label("Log rate constant of the slow time effect on recovery from grade 1 (log 1/h)") # Table 3 K1|0T = -10.4 (0.00532/week)
    e_week1_g01 <- 0.933; label("First-week shift of the grade 0 to 1 logit (log-odds units)") # Table 3 P0|1I = 0.933 (rapid diarrhea onset in week 1)

    # Ixazomib exposure and covariate effects.
    slp_aucwk_g0 <- 0.000715; label("Linear effect of weekly ixazomib AUC on the logits out of grade 0 (per ng*h/mL)") # Table 3 SLP0 = 0.000715 (0.715 per ug*h/mL)
    e_imidnaive_b3 <- -2.87; label("Shift of the logits out of grade 3 for IMiD-naive patients (log-odds units)") # Table 3 B3x (PIMID) = -2.87

    # Random intercepts on b<i>1, diagonal; Table 3 prints the variance
    # (sqrt(1.86) = 1.364 -> 136%CV).
    etab01 ~ 1.86 # Table 3 eta0 = 1.86 (136%CV)
    etab11 ~ 3.72 # Table 3 eta1 = 3.72 (193%CV)
    etab21 ~ 3.24 # Table 3 eta2 = 3.24 (180%CV)
    etab31 ~ 2.74 # Table 3 eta3 = 2.74 (165%CV)
  })

  model({
    # Ixazomib PK (Gupta 2017)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2 + etalvp2) * (BSA / 1.87)^e_bsa_vp2
    fdepot <- exp(lfdepot + etalfdepot)
    tlag <- exp(ltlag)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (cl / vc) * central -
      (q / vc) * central + (q / vp) * peripheral1 -
      (q2 / vc) * central + (q2 / vp2) * peripheral2
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1
    d/dt(peripheral2) <- (q2 / vc) * central - (q2 / vp2) * peripheral2
    f(depot) <- fdepot
    alag(depot) <- tlag

    # mg / L * 1000 = ng/mL
    Cc <- 1000 * central / vc

    # Weekly AUC (ng*h/mL) over the 168 h ending at the evaluation time, which
    # is the week whose transition the probabilities below describe.
    # delay() forces the non-stiff dense dop853 solver; use atol = rtol = 1e-6
    # in rxSolve() if an occasional subject fails at the default tolerances.
    d/dt(auc_central) <- Cc
    aucwk <- auc_central - delay(auc_central, 168)

    # Explicit time effects. Appendix S1 writes p_i|j for the move from grade
    # i to grade j, so P0|1I ("rapid diarrhea onset in week 1") acts on the
    # grade 0 -> 1 logit and P1|0T on the grade 1 -> 0 recovery: added to the
    # grade >= 1 logit from grade 1, it slows recovery and raises the grade 1
    # prevalence over time. The 1 - exp(-k * t) form is not printed; with
    # this placement it reproduces the paper's 1-year prevalences (see the
    # vignette).
    time_g10 <- exp(lemax_time_g10) * (1 - exp(-exp(lk_time_g10) * t))
    week1_g01 <- e_week1_g01 * (t <= 168)

    # Cumulative logits by current grade (Equation 14)
    lg0 <- b01 + etab01 + slp_aucwk_g0 * aucwk
    lg0_ge1 <- lg0 + week1_g01
    lg0_ge2 <- lg0 - exp(b02)
    lg0_ge3 <- lg0_ge2 - exp(b03)

    lg1 <- b11 + etab11
    lg1_ge2 <- lg1 - exp(b12)
    lg1_ge3 <- lg1_ge2 - exp(b13)

    lg2 <- b21 + etab21
    lg2_ge2 <- lg2 - exp(b22)
    lg2_ge3 <- lg2_ge2 - exp(b23)

    lg3 <- b31 + etab31 + e_imidnaive_b3 * (1 - PRIOR_IMID)
    lg3_ge2 <- lg3 - exp(b32)
    lg3_ge3 <- lg3_ge2 - exp(b33)

    # Transition probabilities p<from><to>; each row sums to 1.
    p00 <- 1 - expit(lg0_ge1)
    p01 <- expit(lg0_ge1) - expit(lg0_ge2)
    p02 <- expit(lg0_ge2) - expit(lg0_ge3)
    p03 <- expit(lg0_ge3)

    p10 <- 1 - expit(lg1 + time_g10)
    p11 <- expit(lg1 + time_g10) - expit(lg1_ge2)
    p12 <- expit(lg1_ge2) - expit(lg1_ge3)
    p13 <- expit(lg1_ge3)

    p20 <- 1 - expit(lg2)
    p21 <- expit(lg2) - expit(lg2_ge2)
    p22 <- expit(lg2_ge2) - expit(lg2_ge3)
    p23 <- expit(lg2_ge3)

    p30 <- 1 - expit(lg3)
    p31 <- expit(lg3) - expit(lg3_ge2)
    p32 <- expit(lg3_ge2) - expit(lg3_ge3)
    p33 <- expit(lg3_ge3)
  })
}
