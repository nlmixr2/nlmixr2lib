Williams_2026_cabozantinib <- function() {
  description <- "Two-compartment population PK model for oral cabozantinib tablets in healthy volunteers and adults with advanced solid tumors, including castration-resistant prostate cancer and patients on cabozantinib plus atezolizumab in COSMIC-021 (Williams 2026, n=2757 across 8 studies). Dual-phase absorption: a fraction F1 of the dose is released zero-order (duration D1, lag ALAG1) into a gut depot and absorbed first-order (Ka) into central, while the remaining 1 - F1 enters central directly by a zero-order input (duration D2) after a long lag (ALAG2 about 22 h). Linear elimination from central. Body weight (power, reference 78 kg) acts on CL/F and on V2/F and V3/F with a shared exponent; female sex lowers CL/F by about 20%. IIV on Ka, D1 and F1 applies to healthy volunteers only, and the combined proportional plus additive residual error has separate magnitudes for healthy volunteers and patients."
  reference <- paste(
    "Williams JH, Li H, Othman AA.",
    "Population Pharmacokinetics of Cabozantinib: Updated Analyses Including",
    "Data in Combination With Atezolizumab in Patients With Castration-Resistant",
    "Prostate Cancer and Other Tumor Types.",
    "Clin Transl Sci. 2026;19(10):e70738.",
    "doi:10.1111/cts.70738.",
    sep = " "
  )
  vignette <- "Williams_2026_cabozantinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Each oral administration is two dose records with the same amt and
  # rate = -2: one into depot (fraction F1, zero-order release over D1 after
  # ALAG1) and one into central (fraction 1 - F1, zero-order input over D2
  # after ALAG2). f(), alag() and dur() split and time them, mirroring the
  # source control stream's ADVAN4 F1/F2, ALAG1/ALAG2 and D1/D2.
  dosing <- c("depot", "central")

  # The canonical residual-error matcher recognises only the bare names
  # propSd / addSd / expSd, so the four cohort-stratified residual SDs are
  # declared here (same pattern as Ozdin_2025_dexamethasone.R).
  paper_specific_residual_sds <- c(
    "propSdHealthy",
    "addSdHealthy",
    "propSdPatient",
    "addSdPatient"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot = list(analyte = "cabozantinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cabozantinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cabozantinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline. Power effect (WT / 78)^0.0885 on CL/F and (WT / 78)^1.04 on both V2/F and V3/F (one shared exponent). Reference 78 kg from the Supplementary Material control stream (BWT = 78; also the missing-weight imputation) and Figure 4 (typical patient 78 kg).",
      source_name = "WEIGHT"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Time-fixed. Linear fractional effect CL/F x (1 + (-0.201) x SEXF). Source control stream: IF(SEXN.EQ.1) CLSEXN = (1 + THETA(17)); the Results state female patients had about 20% lower CL/F, so SEXN = 1 is female and maps directly to SEXF.",
      source_name = "SEXN"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer indicator, 1 = healthy volunteer, 0 = patient with cancer",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (patient with cancer)",
      notes = "Time-fixed. Source control stream: HV = 1 when POP = 0. Switches on the IIV of Ka, D1 and F1 (ETA x HV) and selects the healthy-volunteer residual error magnitudes. Only Study XL184-020 (63 healthy volunteers) has DIS_HEALTHY = 1.",
      source_name = "POP (HV = 1 when POP = 0)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      notes = "Tested on CL/F, V2/F and V3/F in the stepwise covariate search; not retained (Results, Discussion)."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested on CL/F; not retained."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested on CL/F; not retained."
    ),
    TBILI = list(
      description = "Baseline total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Tested on CL/F; not retained."
    ),
    CREAT = list(
      description = "Baseline serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Tested on CL/F; not retained."
    ),
    CONMED_ATEZOLIZUMAB = list(
      description = "Combination treatment with atezolizumab 1200 mg IV q3w, 1 = yes",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL/F; no significant effect (Abstract, Results). Race, tumor type (CRPC, RCC, HCC, DTC) and NCI-ODWG hepatic function category were also tested on CL/F (race also on volumes) and not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2757,
    n_studies = 8,
    n_observations = 7513,
    age_range = "19-95 years (median 66)",
    weight_range = "35.0-191 kg (median 77.7)",
    sex_female_pct = 21.4,
    race_ethnicity = c(
      White = 74.9,
      Black = 2.9,
      Asian = 8.9,
      `American Indian/Native American` = 0.3,
      Other = 2.8,
      Missing = 10.3
    ),
    disease_state = "Healthy volunteers (n=63, Study XL184-020) and adults with advanced solid tumors: CRPC 31.8%, RCC 26.4%, HCC 17.5%, NSCLC 5.9%, DTC 5.3%, urothelial carcinoma 4.1%, and colorectal, esophageal, gastric, head and neck, ovarian and triple-negative breast cancers (about 1% each, COSMIC-021 only).",
    dose_range = "20, 40 or 60 mg tablet once daily (single 20-60 mg doses in healthy volunteers); 27.9% received atezolizumab 1200 mg IV q3w (COSMIC-021).",
    hepatic_function = "NCI-ODWG normal 63.2%, mild 19.8%, moderate 0.6%, severe 1 subject, missing 16.3%.",
    studies = "XL184-020, COSMIC-021 (XL184-021), XL184-306, XL184-307 (COMET-1), XL184-308 (METEOR), XL184-309 (CELESTIAL), XL184-311 (COSMIC-311), CA2099ER (CheckMate 9ER).",
    notes = "Baseline demographics from Williams 2026 Table 1 (overall column); per-study demographics in Supplementary Table S3."
  )

  ini({
    # Disposition (Williams 2026 Table 2); typical values for a 78 kg male
    lcl <- log(2.06); label("Apparent clearance CL/F (L/h)") # Table 2 CL/F = 2.06 L/h (RSE 1.05%)
    lvc <- log(103); label("Apparent central volume V2/F (L)") # Table 2 V2/F = 103 L (RSE 3.48%)
    lq <- log(19.3); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2 Q/F = 19.3 L/h (RSE 4.88%)
    lvp <- log(192); label("Apparent peripheral volume V3/F (L)") # Table 2 V3/F = 192 L (RSE 1.43%)

    # Absorption, held at the Study XL184-020 healthy-volunteer estimates
    # (Table 2 footnote: absorption parameters fixed to the XL184-020-alone fit)
    lka <- fixed(log(0.826)); label("First-order absorption rate constant Ka from depot (1/h)") # Table 2 Ka = 0.826 1/h (fix)
    ltlag <- fixed(log(0.359)); label("Lag time ALAG1 of the depot (first-order) pathway (h)") # Table 2 ALAG1 = 0.359 h (fix)
    ld1 <- fixed(log(0.385)); label("Duration D1 of the zero-order release into depot (h)") # Table 2 D1 = 0.385 h (fix)
    ltlag2 <- fixed(log(22.2)); label("Lag time ALAG2 of the direct zero-order input into central (h)") # Table 2 ALAG2 = 22.2 h (fix)
    ld2 <- fixed(log(5.94)); label("Duration D2 of the direct zero-order input into central (h)") # Table 2 D2 = 5.94 h (fix)
    logitffo <- fixed(logit(0.708)); label("Logit of F1, the fraction of the dose on the depot pathway (logit)") # Table 2 F1 = 0.708 (fix); control stream LGF = LOG(THETA(10)/(1-THETA(10)))

    # Covariate effects (Table 2)
    e_wt_cl <- 0.0885; label("Power exponent of WT/78 on CL/F (unitless)") # Table 2 'Weight on CL/F' = 0.0885 (RSE 49.1%)
    e_wt_vc_vp <- 1.04; label("Power exponent of WT/78 on V2/F and V3/F (unitless)") # Table 2 'Weight on V2/F, weight on V3/F' = 1.04 (RSE 4.69%)
    e_sexf_cl <- -0.201; label("Fractional change in CL/F for females (unitless)") # Table 2 'Sex on CL/F' = -0.201 (RSE 9.22%)

    # IIV (Table 2 omega^2). Ka, D1 and F1 etas act only in healthy volunteers
    # (control stream ETA(4), ETA(5), ETA(6) multiplied by HV).
    etalcl + etalvc ~ c(
      0.169,
      0.130, 0.576
    ) # Table 2 omega2 CL/F = 0.169, omega2 V2/F x CL/F = 0.130, omega2 V2/F = 0.576
    etalka ~ 0.373 # Table 2 omega2 Ka = 0.373
    etald1 ~ 1.59 # Table 2 omega2 D1 = 1.59
    etalogitffo ~ 0.424 # Table 2 omega2 F1 = 0.424 (logit scale, control stream LGF)

    # Residual error (Table 2): combined proportional + additive, by cohort
    propSdHealthy <- 0.217; label("Proportional residual SD, healthy volunteers (fraction)") # Table 2 sigma proportional HV = 21.7%
    addSdHealthy <- 0.351; label("Additive residual SD, healthy volunteers (ng/mL)") # Table 2 sigma additive HV = 0.351 ng/mL
    propSdPatient <- 0.302; label("Proportional residual SD, patients (fraction)") # Table 2 sigma proportional patient = 30.2%
    addSdPatient <- 53.7; label("Additive residual SD, patients (ng/mL)") # Table 2 sigma additive patient = 53.7 ng/mL
  })

  model({
    # Covariate effects (Supplementary Material control stream $PK)
    cl <- exp(lcl + etalcl) * (1 + e_sexf_cl * SEXF) * (WT / 78)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 78)^e_wt_vc_vp
    q <- exp(lq)
    vp <- exp(lvp) * (WT / 78)^e_wt_vc_vp

    # Absorption: IIV switched on for healthy volunteers only
    ka <- exp(lka + etalka * DIS_HEALTHY)
    d1 <- exp(ld1 + etald1 * DIS_HEALTHY)
    tlag <- exp(ltlag)
    tlag2 <- exp(ltlag2)
    d2 <- exp(ld2)
    ffo <- expit(logitffo + etalogitffo * DIS_HEALTHY)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- ffo
    alag(depot) <- tlag
    dur(depot) <- d1
    f(central) <- 1 - ffo
    alag(central) <- tlag2
    dur(central) <- d2

    # mg / L x 1000 = ng/mL
    Cc <- central / vc * 1000

    propSd <- propSdHealthy * DIS_HEALTHY + propSdPatient * (1 - DIS_HEALTHY)
    addSd <- addSdHealthy * DIS_HEALTHY + addSdPatient * (1 - DIS_HEALTHY)
    Cc ~ add(addSd) + prop(propSd)
  })
}
