Srimani_2022_ixazomib_rash <- function() {
  description <- "Discrete-time Markov model of weekly rash grade (0-3) for oral ixazomib plus lenalidomide-dexamethasone (LenDex) in relapsed/refractory multiple myeloma from the phase III TOURMALINE-MM1 trial (Srimani 2022). Ixazomib plasma concentrations come from the three-compartment population PK model of Gupta 2017 (fixed) and are integrated to the weekly AUC. For each current grade, the next-week grade follows a proportional-odds cumulative-logit model with a per-grade random intercept. The weekly ixazomib AUC, Asian race and a transient early term that decays over the first months raise all transitions out of grade 0; a slowly rising time effect makes recovery from grade 1 to grade 0 less likely. Grade 3 carries no random effect. The 16 transition probabilities are model outputs; the Markov chain is advanced week by week outside the solve, as in the vignette."
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
    outcome = "weekly rash grade transition probabilities (fraction)"
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
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = non-Asian).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian).",
      notes = "Additive shift of +1.19 on the cumulative logits out of grade 0 (Supplementary Table 5 row 'B0x (RACE)', Asian vs non-Asian), so Asian patients move to every non-zero rash grade more often (Supplementary Figure S12). 8.9% of the safety dataset was Asian (Supplementary Table 4).",
      source_name = "RACE"
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
    biomarkers = "Worst rash grade: 0 in 70.6%, 1 in 18.2%, 2 in 7.9%, 3 in 3.3%; grade 1 or 2 at baseline in 1.1% (Supplementary Table 4).",
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

    # Cumulative-logit transition parameters -- Srimani 2022 Supplementary
    # Table 5. b<i>1 is the logit of moving from grade i to grade >= 1; b<i>2
    # and b<i>3 are log-decrements, so logit(>= 2) = b<i>1 - exp(b<i>2) and
    # logit(>= 3) = logit(>= 2) - exp(b<i>3) (main-text Equation 14). With no
    # drug, time or covariate terms these reproduce Supplementary Equation 14.
    b01 <- -7.37; label("Logit of moving from rash grade 0 to grade >= 1") # Suppl. Table 5 B01 = -7.37
    b02 <- 0.269; label("Log decrement from the grade >= 1 to the grade >= 2 logit, from grade 0") # Suppl. Table 5 B02 = 0.269
    b03 <- 0.436; label("Log decrement from the grade >= 2 to the grade >= 3 logit, from grade 0") # Suppl. Table 5 B03 = 0.436
    b11 <- 0.402; label("Logit of moving from rash grade 1 to grade >= 1") # Suppl. Table 5 B11 = 0.402
    b12 <- 1.87; label("Log decrement from the grade >= 1 to the grade >= 2 logit, from grade 1") # Suppl. Table 5 B12 = 1.87
    b13 <- -0.527; label("Log decrement from the grade >= 2 to the grade >= 3 logit, from grade 1") # Suppl. Table 5 B13 = -0.527
    b21 <- 1.48; label("Logit of moving from rash grade 2 to grade >= 1") # Suppl. Table 5 B21 = 1.48
    b22 <- -2.52; label("Log decrement from the grade >= 1 to the grade >= 2 logit, from grade 2") # Suppl. Table 5 B22 = -2.52
    b23 <- 1.99; label("Log decrement from the grade >= 2 to the grade >= 3 logit, from grade 2") # Suppl. Table 5 B23 = 1.99
    b31 <- 1.85; label("Logit of moving from rash grade 3 to grade >= 1") # Suppl. Table 5 B31 = 1.85
    b32 <- -2.82; label("Log decrement from the grade >= 1 to the grade >= 2 logit, from grade 3") # Suppl. Table 5 B32 = -2.82
    b33 <- -2.19; label("Log decrement from the grade >= 2 to the grade >= 3 logit, from grade 3") # Suppl. Table 5 B33 = -2.19

    # Explicit time effects. The slow effect magnitude is printed on the log
    # scale (exp(0.921) = 2.51 log OR); rate constants are in log 1/h
    # (exp(-7.83) * 168 = 0.0667/week).
    lemax_time_g10 <- 0.921; label("Log maximal slow time effect slowing recovery from grade 1 to grade 0 (log of log-odds units)") # Suppl. Table 5 P1|0T = 0.921 (2.51 log OR)
    lk_time_g10 <- -7.83; label("Log rate constant of the slow time effect on recovery from grade 1 (log 1/h)") # Suppl. Table 5 K1|0T = -7.83 (0.0666/wk)
    e_transient_g0 <- 3.81; label("Initial transient shift of the logits out of grade 0 (log-odds units)") # Suppl. Table 5 ETIME = 3.81
    lk_transient_g0 <- -7.13; label("Log decay rate constant of the transient grade 0 term (log 1/h)") # Suppl. Table 5 KTIME = -7.13 (0.134/wk)

    # Ixazomib exposure and covariate effects.
    slp_aucwk_g0 <- 0.000929; label("Linear effect of weekly ixazomib AUC on the logits out of grade 0 (per ng*h/mL)") # Suppl. Table 5 SLP0 = 0.000929 (0.929 per ug*h/mL)
    e_asian_b0 <- 1.19; label("Shift of the logits out of grade 0 for Asian patients (log-odds units)") # Suppl. Table 5 B0x (RACE) = 1.19

    # Random intercepts on b01, b11 and b21 (none on b31); Supplementary
    # Table 5 prints the variance (sqrt(1.80) = 1.342 -> 134%cv).
    etab01 ~ 1.80 # Suppl. Table 5 eta0 = 1.80 (134%cv)
    etab11 ~ 0.685 # Suppl. Table 5 eta1 = 0.685 (82.8%cv)
    etab21 ~ 0.523 # Suppl. Table 5 eta2 = 0.523 (72.3%cv)
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

    # Appendix S1 writes p_i|j for the move from grade i to grade j, so the
    # slow effect P1|0T acts on the grade 1 -> 0 recovery: added to the
    # grade >= 1 logit from grade 1, it slows recovery and raises the grade 1
    # prevalence over time. The 1 - exp(-k * t) form is not printed; with
    # this placement it reproduces the supplement's 1-year rash prevalences
    # (see the vignette). The transient term, like the drug and race effects,
    # shifts every logit out of grade 0.
    time_g10 <- exp(lemax_time_g10) * (1 - exp(-exp(lk_time_g10) * t))
    transient_g0 <- e_transient_g0 * exp(-exp(lk_transient_g0) * t)

    # Cumulative logits by current grade (main-text Equation 14)
    lg0 <- b01 + etab01 + slp_aucwk_g0 * aucwk + e_asian_b0 * RACE_ASIAN + transient_g0
    lg0_ge1 <- lg0
    lg0_ge2 <- lg0 - exp(b02)
    lg0_ge3 <- lg0_ge2 - exp(b03)

    lg1 <- b11 + etab11
    lg1_ge2 <- lg1 - exp(b12)
    lg1_ge3 <- lg1_ge2 - exp(b13)

    lg2 <- b21 + etab21
    lg2_ge2 <- lg2 - exp(b22)
    lg2_ge3 <- lg2_ge2 - exp(b23)

    lg3 <- b31
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
