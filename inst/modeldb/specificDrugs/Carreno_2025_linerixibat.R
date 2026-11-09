Carreno_2025_linerixibat <- function() {
  description <- "Kinetic-pharmacodynamic (K-PD) model of serum 7-alpha-hydroxy-4-cholesten-3-one (C4), a marker of bile acid synthesis, after oral linerixibat (an ileal bile acid transporter inhibitor) in healthy volunteers and patients with primary biliary cholangitis, linked to a longitudinal proportional-odds model of the patient-reported GSRS diarrhea-bother score (Question 11, scores 1-7) in patients with primary biliary cholangitis and pruritus (Carreno 2025). Linerixibat is minimally absorbed, so the dose enters a virtual compartment that empties first-order at a fixed rate; the dose rate leaving it stimulates the zero-order C4 synthesis rate through an Emax function. C4 synthesis follows a dual-cosine diurnal rhythm (8-h and 12-h periods) whose first amplitude falls with the individual baseline C4. The cumulative logits of the GSRS score carry an Emax function of the model-predicted C4 concentration and one between-subject random effect shared by all thresholds."
  reference <- paste(
    "Carreno F, Mehta R, Sherman AR, Collins J, Swift B.",
    "Analysis of C4 Concentrations to Predict Impact of Patient-Reported Diarrhea",
    "Associated With the Ileal Bile Acid Transporter Inhibitor Linerixibat.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(3):596-605.",
    "doi:10.1002/psp4.13300"
  )
  vignette <- "Carreno_2025_linerixibat"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL",
    outcome = "GSRS Question 11 diarrhea-bother score (ordinal, 1-7)"
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened by stepwise covariate modeling on the K-PD parameters; not retained (Carreno 2025 Methods and Discussion).",
      source_name = "body weight"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened by stepwise covariate modeling on the K-PD parameters; not retained.",
      source_name = "age"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened by stepwise covariate modeling on the K-PD parameters; not retained.",
      source_name = "body mass index"
    ),
    DIS_PBC = list(
      description = "Participant status, patient with primary biliary cholangitis (1) versus healthy volunteer (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy volunteer)",
      notes = "Screened by stepwise covariate modeling on the K-PD parameters; not retained, so healthy volunteers and patients share one parameter set.",
      source_name = "participant status"
    ),
    RACE_ASIAN = list(
      description = "Race indicator (1 = Asian)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "Race was screened on the K-PD and on the proportional-odds models (38 Japanese versus 103 Caucasian participants, Discussion) and not retained.",
      source_name = "race"
    )
  )

  compartmentData <- list(
    depot_kpd = list(analyte = "linerixibat", units = "mg", specimen = "administration site", verified = TRUE),
    c4 = list(analyte = "7-alpha-hydroxy-4-cholesten-3-one (C4)", units = "ng/mL", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 255L,
    n_studies = 5L,
    age_range = "Geometric mean age 31 to 55 years across the five studies",
    weight_range = "Geometric mean body weight 62.8 to 82.2 kg across the five studies",
    sex_female_pct = 63,
    race_ethnicity = "Mixed; the data included 38 Japanese and 103 Caucasian participants (Discussion).",
    disease_state = "86 healthy volunteers (three phase 1 studies) and 169 patients with primary biliary cholangitis and pruritus (phase 2a NCT01899703 and phase 2b GLIMMER, NCT02966834). The proportional-odds GSRS model used the 147 GLIMMER patients only.",
    dose_range = "Linerixibat 0.1-180 mg once daily or 3-90 mg twice daily (oral solution or tablet), single dose to 12 weeks.",
    regions = "Multinational (phase 1 studies in the US, UK and Japan; GLIMMER multinational).",
    biomarkers = "Baseline serum C4 geometric mean 7.4-14.1 ng/mL across studies (Table 1); C4 assay range 2-200 ng/mL (LC-MS/MS).",
    notes = "3292 C4 concentrations from 255 participants (Table 1). GSRS dataset: 147 GLIMMER patients, 713 observations in the 4-week single-blind placebo run-in and 1640 in the 12-week treatment period."
  )

  ini({
    # K-PD virtual compartment (Supplementary Methods S1; Table 2).
    lkel <- fixed(log(0.110))
    label("First-order elimination rate from the virtual dosing compartment KDE (1/h)") # Table 2: KDE = 0.110 1/h FIX
    lemax <- log(5.55)
    label("Maximal fractional stimulation of C4 synthesis Emax (unitless)") # Table 2: Emax = 5.55 (RSE 10.4%)
    ledk50 <- log(3.03)
    label("Dose rate leaving the virtual compartment giving half-maximal stimulation EDK50 (mg/h)") # Table 2: EDK50 = 3.03 mg/h (RSE 23.5%)

    # C4 indirect response with dual-cosine synthesis rhythm (Table 2).
    lrbase <- log(11.9)
    label("Baseline serum C4 concentration, Kin/Kout (ng/mL)") # Table 2: Baseline = 11.9 ng/mL (RSE 5.58%)
    lkout <- log(0.135)
    label("First-order C4 elimination rate constant Kout (1/h)") # Table 2: Kout = 0.135 1/h (RSE 6.19%)
    kin_circ_famp1 <- 0.268
    label("Fractional amplitude of the 8-h cosine on C4 synthesis at the typical baseline (unitless)") # Table 2: Amplitude1 = 0.268 (RSE 21.3%)
    tacro1 <- 5.77
    label("Acrophase of the 8-h cosine on C4 synthesis (h)") # Table 2: Acrophase1 = 5.77 h (RSE 0.996%)
    e_rbase_kin_circ_famp1 <- -0.0130
    label("Linear change in the 8-h fractional amplitude per ng/mL of individual baseline C4 above the typical baseline (mL/ng)") # Table 2: Effect of baseline C4 on amplitude = -0.0130 (RSE 61.7%)
    kin_circ_famp2 <- 0.0227
    label("Fractional amplitude of the 12-h cosine on C4 synthesis (unitless)") # Table 2: Amplitude2 = 0.0227 (RSE 38.7%)
    tacro2 <- 21.7
    label("Acrophase of the 12-h cosine on C4 synthesis (h)") # Table 2: Acrophase2 = 21.7 h (RSE 3.86%)

    # Proportional-odds model of GSRS Question 11 (Table 2). The paper
    # parameterises P(score >= m): intercept m = 2 is logit P(score >= 2) and
    # the m = 3..7 rows are the positive decrements to each next threshold.
    # Here logit P(score <= m) = -logit P(score >= m + 1), so logite0 is the
    # sign-flipped m = 2 intercept and the decrements become logitstep1-5.
    logite0 <- 4.05
    label("Cumulative logit of P(GSRS score <= 1) without drug, minus the m = 2 intercept (unitless)") # Table 2: Intercept m = 2 = -4.05 (RSE 8%)
    logitstep1 <- 1.26
    label("Increment from the logit of P(score <= 1) to the logit of P(score <= 2) (unitless)") # Table 2: Intercept m = 3 = 1.26 (RSE 3%)
    logitstep2 <- 0.983
    label("Increment from the logit of P(score <= 2) to the logit of P(score <= 3) (unitless)") # Table 2: Intercept m = 4 = 0.983 (RSE 4%)
    logitstep3 <- 1.15
    label("Increment from the logit of P(score <= 3) to the logit of P(score <= 4) (unitless)") # Table 2: Intercept m = 5 = 1.15 (RSE 5%)
    logitstep4 <- 1.11
    label("Increment from the logit of P(score <= 4) to the logit of P(score <= 5) (unitless)") # Table 2: Intercept m = 6 = 1.11 (RSE 6%)
    logitstep5 <- 1.53
    label("Increment from the logit of P(score <= 5) to the logit of P(score <= 6) (unitless)") # Table 2: Intercept m = 7 = 1.53 (RSE 7%)
    lemax_gsrs_diarrhea <- log(11.6)
    label("Maximal C4 effect on the cumulative logits of a higher GSRS diarrhea score Emax (unitless)") # Table 2: PO Emax = 11.6 (RSE 4%)
    lec50_gsrs_diarrhea <- log(42.9)
    label("Serum C4 concentration giving the half-maximal effect on the GSRS logits EC50 (ng/mL)") # Table 2: PO EC50 = 42.9 ng/mL (RSE 10%)

    # IIV. Table 2 footnote: %CV = sqrt(exp(omega^2) - 1) for exponential
    # random effects, so omega^2 = log(1 + CV^2); the logit random effect is
    # additive and the column reports its SD, so omega^2 = SD^2.
    etalrbase ~ log(1 + 0.892^2) # Table 2: Baseline IIV 89.2% CV
    etalemax ~ log(1 + 0.558^2) # Table 2: Emax IIV 55.8% CV
    etaledk50 ~ log(1 + 3.58^2) # Table 2: EDK50 IIV 358% CV
    etalogite0 ~ 5.44^2 # Table 2: intercept m = 2 IIV SD 5.44

    # Residual error on C4 (Table 2), combined additive + proportional.
    addSd_c4 <- fixed(0.000492)
    label("Additive residual SD of serum C4 (ng/mL)") # Table 2: Additive error (SD) = 0.000492 FIX
    propSd_c4 <- 0.448
    label("Proportional residual SD of serum C4 (fraction)") # Table 2: Proportional error = 44.8 CV% (RSE 2.05%)
  })

  model({
    # Individual parameters
    kel <- exp(lkel)
    emax <- exp(lemax + etalemax)
    edk50 <- exp(ledk50 + etaledk50)
    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout)
    kin <- rbase * kout

    # Linerixibat virtual compartment; dr is the dose rate (mg/h) leaving it.
    d/dt(depot_kpd) <- -kel * depot_kpd
    dr <- kel * depot_kpd
    eff <- emax * dr / (edk50 + dr)

    # Dual-cosine diurnal rhythm on C4 synthesis. t is the time axis of the
    # source data (time 0 = first dose on Day 1). The 8-h amplitude changes
    # linearly with the individual baseline, centred on the typical baseline;
    # the paper names this effect but does not print its equation.
    famp1 <- kin_circ_famp1 + e_rbase_kin_circ_famp1 * (rbase - exp(lrbase))
    cyr1 <- famp1 * cos((t - tacro1) * 2 * pi / 8)
    cyr2 <- kin_circ_famp2 * cos((t - tacro2) * 2 * pi / 12)

    c4(0) <- rbase
    d/dt(c4) <- kin * (1 + cyr1 + cyr2) * (1 + eff) - kout * c4

    # Proportional-odds GSRS model driven by the predicted C4 concentration.
    # logit P(score <= m) for m = 1..6; a higher C4 lowers each one.
    eff_gsrs <- exp(lemax_gsrs_diarrhea) * c4 / (exp(lec50_gsrs_diarrhea) + c4)
    lp <- logite0 + etalogite0 - eff_gsrs
    pcum1 <- expit(lp)
    pcum2 <- expit(lp + logitstep1)
    pcum3 <- expit(lp + logitstep1 + logitstep2)
    pcum4 <- expit(lp + logitstep1 + logitstep2 + logitstep3)
    pcum5 <- expit(lp + logitstep1 + logitstep2 + logitstep3 + logitstep4)
    pcum6 <- expit(lp + logitstep1 + logitstep2 + logitstep3 + logitstep4 + logitstep5)
    p1 <- pcum1
    p2 <- pcum2 - pcum1
    p3 <- pcum3 - pcum2
    p4 <- pcum4 - pcum3
    p5 <- pcum5 - pcum4
    p6 <- pcum6 - pcum5
    # Probability of moderate-to-very-severe bother (score >= 4).
    p_ge4 <- 1 - pcum3

    c4 ~ add(addSd_c4) + prop(propSd_c4)
    # rxode2 ordinal endpoint: simulated category k is GSRS score k (1-7).
    gsrs_diarrhea ~ c(p1, p2, p3, p4, p5, p6)
  })
}
