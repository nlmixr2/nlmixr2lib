Hanke_2026_iclepertin <- function() {
  description <- paste(
    "Two-compartment population PK model with lagged first-order absorption for oral",
    "iclepertin (BI 425809, GlyT1 inhibitor) tablets in healthy adults and in patients with",
    "cognitive impairment associated with schizophrenia (CIAS) or Alzheimer's disease.",
    "Body-weight power effects on all disposition parameters; age, race and smoking on",
    "apparent clearance; food and diagnosis on the absorption rate; food and dose level on",
    "relative bioavailability. Supplies the individual steady-state AUC that drives the",
    "companion hemoglobin model Hanke_2026_iclepertin_hemoglobin."
  )
  reference <- paste(
    "Hanke N, Callisto SP, Nagy P, Desch M, Abdelwahab MT.",
    "Pharmacokinetic/pharmacodynamic analysis to characterize the effect of long-term",
    "GlyT1 inhibitor iclepertin exposure on hemoglobin levels.",
    "Front Pharmacol. 2026;17:1770980. doi:10.3389/fphar.2026.1770980.",
    "The popPK parameter estimates are in the article's Supplementary Material, Table S2."
  )
  vignette <- "Hanke_2026_iclepertin"

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effects on CL/F (estimated), V2/F (estimated), Q/F (0.75 fixed) and V3/F",
        "(1.00 fixed) per Table S2. The reference weight is not reported; 70 kg is assumed",
        "(see the vignette's Assumptions and deviations section)."
      ),
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL/F, exponent -1.64 (Table S2 'CL AGE'). Neither the functional form",
        "nor the reference age is reported; a normalized power model (AGE / 40)^theta is",
        "assumed, the form the paper states for its continuous PD covariates, with 40 years as",
        "a round adult reference (see the vignette's Assumptions and deviations section)."
      ),
      source_name = "AGE"
    ),
    RACE_BLACK = list(
      description = "Black race indicator (1 = Black, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "White (all race indicators 0)",
      notes = "Multiplicative factor 1.35 on CL/F (Table S2 'CL BLACK').",
      source_name = "RACE"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "White (all race indicators 0)",
      notes = "Multiplicative factor 1.00 on CL/F (Table S2 'CL ASIAN', 95% CDI 0.910-1.10).",
      source_name = "RACE"
    ),
    RACE_OTHER = list(
      description = "Other race indicator (1 = race other than White, Black or Asian)",
      units = "(binary)",
      type = "binary",
      reference_category = "White (all race indicators 0)",
      notes = "Multiplicative factor 0.963 on CL/F (Table S2 'CL OTHER').",
      source_name = "RACE"
    ),
    SMOKE = list(
      description = "Current smoker indicator (1 = smoker, 0 = non-smoker)",
      units = "(binary)",
      type = "binary",
      reference_category = "non-smoker (SMOKE = 0)",
      notes = "Multiplicative factor 1.02 on CL/F (Table S2 'CL SMOKE').",
      source_name = "SMOKE"
    ),
    FED = list(
      description = "Dose taken with food (1 = fed, 0 = fasted)",
      units = "(binary)",
      type = "binary",
      reference_category = "fasted (FED = 0)",
      notes = "Multiplicative factors 1.68 on KA and 1.15 on relative bioavailability (Table S2 'KA FED', 'F FED').",
      source_name = "FED"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer indicator (1 = healthy volunteer, 0 = patient with CIAS or Alzheimer's disease)",
      units = "(binary)",
      type = "binary",
      reference_category = "patient (DIS_HEALTHY = 0)",
      notes = paste(
        "Table S2 reports 'KA DIAG, Effect of diagnosis on KA' = 2.31 without stating the",
        "coding. It is applied here as a factor on KA for subjects WITH a diagnosis (the Phase",
        "II CIAS and Alzheimer's disease patients), i.e. ka x 2.31^(1 - DIS_HEALTHY); see the",
        "vignette's Assumptions and deviations section. Affects only the absorption rate, so",
        "the steady-state AUC is unaffected."
      ),
      source_name = "DIAG"
    ),
    DOSE = list(
      description = "Administered iclepertin dose level",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Selects the dose-level relative bioavailability: 0.931 for 10 mg and 0.823 for 25 mg",
        "doses (Table S2 'F1-10mg', 'F1-25mg'), 1 for doses below 10 mg. Doses of 50 mg and",
        "above (Phase I only) are assigned the 25 mg value; doses between the reported levels",
        "take the band below them (see the vignette's Assumptions and deviations section)."
      ),
      source_name = "DOSE"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "iclepertin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "iclepertin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "iclepertin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1031,
    n_studies = 9,
    disease_state = paste(
      "Healthy volunteers (seven Phase I studies: single and multiple rising dose, proof of",
      "mechanism, Japanese/Chinese single rising dose, itraconazole DDI, absolute",
      "bioavailability, human ADME) and patients with cognitive impairment associated with",
      "schizophrenia (Phase II, n = 325) or Alzheimer's disease (Phase II, n = 477)"
    ),
    dose_range = "0.5-150 mg single oral doses; 2-75 mg once daily for 12 days to 12 weeks",
    regions = "Not reported (Japanese and Chinese single-rising-dose study included)",
    notes = paste(
      "Supplementary Table S1 lists the nine studies and per-dose-group subject counts",
      "(SRD 65, MRD 75, proof of mechanism 25, Japanese/Chinese SRD 36, itraconazole DDI 16,",
      "absolute BA 6, human ADME 6, Phase II CIAS 325, Phase II AD 477; placebo excluded).",
      "Demographic summaries for the popPK dataset are not reported. Geometric mean",
      "steady-state exposures for the Phase II CIAS population are in Table S3",
      "(10 mg once daily: AUCtau,ss 4020 nmol*h/L, 64.5% CV). NONMEM 7.4.3, FOCE-I."
    )
  )

  ini({
    # Typical values: White, non-smoking, fasted, dose < 10 mg, at the (assumed)
    # reference weight 70 kg and age 40 years. KA is the healthy-volunteer value
    # (the diagnosis factor multiplies it for patients).
    lka <- log(0.383); label("Absorption rate constant KA (1/h)") # Table S2 KA = 0.383 (95% CDI 0.318-0.462)
    lcl <- log(3.98); label("Apparent clearance CL/F (L/h)") # Table S2 CL/F = 3.98 (3.70-4.28)
    lvc <- log(120); label("Apparent central volume V2/F (L)") # Table S2 V2/F = 120 (109-133)
    lq <- log(6.99); label("Apparent intercompartmental clearance Q/F (L/h)") # Table S2 Q/F = 6.99 (5.91-8.14)
    lvp <- log(86.2); label("Apparent peripheral volume V3/F (L)") # Table S2 V3/F = 86.2 (77.5-95.3)
    ltlag <- log(0.233); label("Absorption lag time ALAG (h)") # Table S2 ALAG = 0.233 (0.230-0.235)
    lfdepot <- fixed(log(1)); label("Relative bioavailability for doses below 10 mg, fasted (reference, unitless)") # Table S2 reports F1 only as dose-level factors, so the lowest doses are the reference

    # Body-weight power exponents on (WT / 70)
    e_wt_cl <- 0.164; label("Power exponent of body weight on CL/F (unitless)") # Table S2 CL WT = 0.164 (0.0120-0.321)
    e_wt_vc <- 0.864; label("Power exponent of body weight on V2/F (unitless)") # Table S2 V2 WT = 0.864 (0.673-1.05)
    e_wt_q <- fixed(0.75); label("Power exponent of body weight on Q/F (unitless)") # Table S2 Q WT = 0.750 Fixed
    e_wt_vp <- fixed(1); label("Power exponent of body weight on V3/F (unitless)") # Table S2 V3 WT = 1.00 Fixed

    # Covariate effects on CL/F
    e_age_cl <- -1.64; label("Power exponent of age on CL/F (unitless)") # Table S2 CL AGE = -1.64 (-2.07 to -1.21)
    e_race_black_cl <- 1.35; label("CL/F multiplier for Black race (unitless)") # Table S2 CL BLACK = 1.35 (1.19-1.53)
    e_race_asian_cl <- 1.00; label("CL/F multiplier for Asian race (unitless)") # Table S2 CL ASIAN = 1.00 (0.910-1.10)
    e_race_other_cl <- 0.963; label("CL/F multiplier for other race (unitless)") # Table S2 CL OTHER = 0.963 (0.793-1.17)
    e_smoke_cl <- 1.02; label("CL/F multiplier for current smokers (unitless)") # Table S2 CL SMOKE = 1.02 (0.935-1.12)

    # Covariate effects on absorption and bioavailability
    e_fed_ka <- 1.68; label("KA multiplier when dosed with food (unitless)") # Table S2 KA FED = 1.68 (1.48-1.90)
    e_dis_healthy_ka <- 2.31; label("KA multiplier for patients with a diagnosis vs healthy volunteers (unitless)") # Table S2 KA DIAG = 2.31 (1.69-3.29)
    e_fed_fdepot <- 1.15; label("Relative bioavailability multiplier when dosed with food (unitless)") # Table S2 F FED = 1.15 (1.10-1.20)
    e_dose_10mg_fdepot <- 0.931; label("Relative bioavailability of 10 mg doses (unitless)") # Table S2 F1-10mg = 0.931 (0.872-0.998)
    e_dose_25mg_fdepot <- 0.823; label("Relative bioavailability of 25 mg doses (unitless)") # Table S2 F1-25mg = 0.823 (0.775-0.873)

    # IIV: Table S2 variances and covariances (correlations in parentheses there)
    etalcl + etalvc + etalka ~ c(
      0.242,
      -0.0371, 0.118,
      -0.171, 0.227, 1.03
    ) # Table S2 omega2 CL/F 0.242 (52.4%CV), COV CL-V2 -0.0371, omega2 V2/F 0.118 (35.4%CV), COV CL-KA -0.171, COV V2-KA 0.227, omega2 KA 1.03 (134%CV)

    # Residual error: sigma2 PROP 0.0919 (30.3%CV) -> SD = sqrt(0.0919)
    propSd <- 0.30315; label("Proportional residual error (fraction)") # Table S2 sigma2 PROP = 0.0919, sqrt = 0.30315 (30.3%CV)
  })

  model({
    # Covariate multipliers. The reference weight (70 kg) and age (40 years)
    # are not reported in the source and are assumptions.
    raceCl <- e_race_black_cl^RACE_BLACK * e_race_asian_cl^RACE_ASIAN *
      e_race_other_cl^RACE_OTHER

    # Dose-level relative bioavailability: < 10 mg reference, 10 mg, >= 25 mg.
    fDose <- 1
    if (DOSE >= 10) fDose <- e_dose_10mg_fdepot
    if (DOSE >= 25) fDose <- e_dose_25mg_fdepot

    ka <- exp(lka + etalka) * e_fed_ka^FED * e_dis_healthy_ka^(1 - DIS_HEALTHY)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (AGE / 40)^e_age_cl * raceCl *
      e_smoke_cl^SMOKE
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_q
    vp <- exp(lvp) * (WT / 70)^e_wt_vp
    tlag <- exp(ltlag)
    fdepot <- exp(lfdepot) * fDose * e_fed_fdepot^FED

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    # Dose in mg and volumes in L give mg/L; x 1000 converts to ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
