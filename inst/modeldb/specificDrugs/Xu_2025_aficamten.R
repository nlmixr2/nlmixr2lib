Xu_2025_aficamten <- function() {
  description <- "Two-compartment population PK model for aficamten with first-order absorption and a formulation-dependent lag time in healthy participants and participants with obstructive hypertrophic cardiomyopathy (Xu 2025)"
  reference <- "Xu D, Li H, Heitner SB, Jacoby DL, Kupfer S, German P, Lutz J. Toward a Quantitative Understanding of Aficamten Clinical Pharmacology: Population Pharmacokinetic Modeling. CPT Pharmacometrics Syst Pharmacol. 2025;14(12):1982-1992. doi:10.1002/psp4.70099"
  vignette <- "Xu_2025_aficamten"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Confirmed against the NONMEM control stream in
  # Supporting Information Data S1 (ADVAN4 TRANS4: depot + central +
  # peripheral, S2 = V2/1000 so that a dose in mg gives ng/mL).
  compartmentData <- list(
    depot       = list(analyte = "aficamten", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "aficamten", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "aficamten", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (baseline value). Power model normalised to a standard 80 kg, applied with a shared exponent across CL/F and Q/F and a second shared exponent across Vc/F and Vp/F (Xu 2025 Figure S3).",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Source NONMEM column SEX is coded 1 = male, 2 = female; SEXF = 1 for SEX == 2. Proportional effect on CL/F and Vp/F relative to the male reference (Xu 2025 Table 1).",
      source_name        = "SEX"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant indicator (vs participant with obstructive hypertrophic cardiomyopathy)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (participant with oHCM)",
      notes              = "Source NONMEM column PTYPE; PTYPE == 1 identifies healthy participants. Proportional effect on CL/F and Vp/F. The oHCM cohort is the reference, so the typical values in ini() describe a male participant with oHCM weighing 80 kg (Xu 2025 Section 3.1.3).",
      source_name        = "PTYPE"
    ),
    FORM_CAPSULE = list(
      description        = "Capsule formulation indicator (vs tablet)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tablet: phase 2, phase 3, or commercial tablet)",
      notes              = "Source NONMEM column FORM (<= 1 capsule; 2 phase 2 tablet; 3 phase 3 tablet; 4 commercial tablet). Selects between two separately estimated absorption lag times; the three tablet levels share a single lag time. Formulation affects only Tlag in the final model -- no formulation effect on relative bioavailability or Ka was retained (Xu 2025 Section 3.1.3).",
      source_name        = "FORM"
    ),
    FED = list(
      description        = "Fed-state-at-dosing indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (fed: standard meal, high-fat meal, or unspecified 'with or without food')",
      notes              = "Source NONMEM column FOOD (0 fasted; 1 high-fat meal; 2 standard meal); FED = 0 only when FOOD == 0. The fasted effect on Ka is applied as (1 - FED), so the reference arm is the standard meal. Phase 2 / phase 3 oHCM records, which the source records as 'with or without food' and which set none of the NONMEM FASTED / HIGHFAT / STDMEAL indicators, therefore carry FED = 1 to place them on the reference arm (Xu 2025 Table S5 and Data S1 control stream).",
      source_name        = "FOOD"
    ),
    FED_HIGHFAT = list(
      description        = "High-fat-meal-at-dosing indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted, standard meal, or unspecified 'with or without food')",
      notes              = "Source NONMEM column FOOD == 1. Carries all three high-fat-meal effects (relative bioavailability, Tlag, Ka). High-fat-meal records also carry FED = 1. Evaluated only in the phase 1 food-effect study, in which a single 20 mg dose was given fasted or fed (Xu 2025 Table S1).",
      source_name        = "FOOD"
    )
  )

  # Covariates the source screened but did not retain. Each of these carries an
  # explicitly zero-fixed coefficient in the Data S1 control stream, so the
  # final model is numerically identical to one that omits them; they are
  # recorded here to preserve the provenance of the paper's covariate screen.
  # Disease-related covariates (LVEF, NYHA class, NT-proBNP) and concomitant
  # medications (beta-blockers, calcium channel blockers, disopyramide,
  # amiodarone) were also screened and dropped (Xu 2025 Section 3.1.2); they
  # have no coefficient at all in the control stream and no canonical column in
  # inst/references/covariate-columns.md, so they are described in
  # population$notes rather than given entries here.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units       = "years",
      type        = "continuous",
      notes       = "Significant on CL/F and Vp/F in univariate screening (p < 0.01) but removed in backward elimination (p > 0.001). Data S1 control stream retains the power terms (AGE/60)^THETA(11) on CL/F and (AGE/60)^THETA(23) on Vp/F with both thetas '0 FIX' (Xu 2025 Section 3.1.2)."
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units       = "g/dL (source convention; canonical column units are g/L)",
      type        = "continuous",
      notes       = "Significant on CL/F in univariate screening but removed in backward elimination. Data S1 control stream retains (ALB/4.5)^THETA(12) with THETA(12) '0 FIX'; the 4.5 normalisation constant confirms the source reported albumin in g/dL (Xu 2025 Section 3.1.2)."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Significant on CL/F in univariate screening but removed in backward elimination. Data S1 control stream retains (AST/20)^THETA(13) with THETA(13) '0 FIX' (Xu 2025 Section 3.1.2)."
    ),
    CRCL = list(
      description = "Baseline creatinine clearance",
      units       = "mL/min (source convention; canonical column units are mL/min/1.73 m^2)",
      type        = "continuous",
      notes       = "Significant on CL/F in univariate screening but removed in backward elimination; mild and moderate renal impairment had no significant effect in the final model. Data S1 control stream retains (CRCL/100)^THETA(14) with THETA(14) '0 FIX' (Xu 2025 Sections 3.1.2 and 3.1.3)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 447,
    n_studies      = 9,
    age_range      = "18-83 years",
    age_median     = "43 years",
    weight_range   = "47.0-156 kg",
    weight_median  = "77.2 kg",
    sex_female_pct = 39.4,
    race_ethnicity = c(
      White = 71.1, Black = 11.2, Asian = 13.4,
      AmericanIndianAlaskaNative = 0.9, NativeHawaiianPacificIslander = 0.2,
      Multiple = 2.5, Other = 0.7
    ),
    disease_state  = "264 healthy participants (7 phase 1 studies) and 183 participants with symptomatic obstructive hypertrophic cardiomyopathy (phase 2 REDWOOD-HCM, NCT04219826; phase 3 SEQUOIA-HCM, NCT05186818)",
    dose_range     = "single oral doses of 1-75 mg; multiple once-daily oral doses of 5-30 mg",
    renal_function = "Normal 83.7%, mild impairment 13.4%, moderate impairment 2.9% (by CrCL; no severe impairment)",
    n_observations = 9963,
    notes          = "Baseline demographics per Xu 2025 Table S5; study designs and sampling schedules per Table S1; sample accounting per Table S3. Samples below the 1 ng/mL LLOQ (n = 373) were excluded from the analysis. Participants with hepatic impairment and with non-obstructive HCM were excluded from the population PK dataset. Covariates screened but not retained in the final model (Xu 2025 Section 3.1.2): age, serum albumin, AST, and creatinine clearance on CL/F, age on Vp/F, and race on CL/F and Vc/F -- all significant in univariate screening (p < 0.01) but eliminated at the stricter backward-elimination threshold (p > 0.001), and all retained in the Data S1 control stream with zero-fixed coefficients (see covariatesDataExcluded; race enters as a 3-indicator contrast against a reference level via THETA(17)-THETA(19) on CL/F and THETA(20)-THETA(22) on Vc/F, and the source does not state which numeric RACE level is the reference). Baseline disease-related covariates (left ventricular ejection fraction, NYHA class, N-terminal prohormone of brain natriuretic peptide) and concomitant medications (beta-blockers, calcium channel blockers, disopyramide, amiodarone) were screened after the demographic covariates were incorporated and showed no significant impact on clearance (p >= 0.01); these carry no coefficient in the control stream."
  )

  ini({
    # Structural parameters. Typical values describe a MALE participant with
    # oHCM weighing 80 kg, dosed as a tablet without a high-fat meal
    # (Xu 2025 Section 3.1.3 and Figure S3).
    lcl <- log(2.62); label("Apparent clearance (CL/F, L/h)")                            # Xu 2025 Table 1
    lvc <- log(18.1); label("Apparent central volume of distribution (Vc/F, L)")         # Xu 2025 Table 1
    lq  <- log(57.6); label("Apparent intercompartmental clearance (Q/F, L/h)")          # Xu 2025 Table 1
    lvp <- log(295);  label("Apparent peripheral volume of distribution (Vp/F, L)")      # Xu 2025 Table 1

    # Ka: Table 1 reports 0.337 1/h; Figure S3 prints 0.373 1/h for the same
    # quantity. Table 1 is used -- see the vignette Errata section for the full
    # reconciliation. Table 1 lists the unit as "L/h", which is a typographical
    # error; Ka is a first-order rate constant (1/h) per the Data S1 control
    # stream comment ";5 KA 1/hr".
    lka <- log(0.337); label("First-order absorption rate constant (Ka, 1/h)")           # Xu 2025 Table 1

    # Absorption lag time: separately estimated for the capsule (phase 1) and
    # the pooled tablet (phase 2 / phase 3 / commercial) formulations.
    ltlag_tab <- log(0.229); label("Absorption lag time, tablet formulations (Tlag, h)") # Xu 2025 Table 1
    ltlag_cap <- log(0.248); label("Absorption lag time, capsule formulation (Tlag, h)") # Xu 2025 Table 1

    # Relative bioavailability is anchored at 1 for the reference arm
    # (Data S1 control stream: TVF1 = 1 * (1 + FASTED*THETA(29)) * ...,
    # with the fasted and standard-meal coefficients fixed to 0).
    lfdepot <- fixed(log(1)); label("Relative bioavailability of the reference arm (F1, unitless)") # Xu 2025 Data S1 control stream

    # Body-weight power exponents, normalised to a standard 80 kg. One exponent
    # is shared by CL/F and Q/F, a second by Vc/F and Vp/F.
    e_wt_cl_q  <- 0.586; label("Body-weight power exponent on CL/F and Q/F (unitless)")   # Xu 2025 Table 1
    e_wt_vc_vp <- 0.882; label("Body-weight power exponent on Vc/F and Vp/F (unitless)")  # Xu 2025 Table 1

    # Categorical covariate effects, all proportional: P = theta * (1 + e * Q)
    # (Xu 2025 Figure S1).
    e_dis_healthy_cl <-  0.356; label("Effect of healthy status (vs oHCM) on CL/F (fraction)")  # Xu 2025 Table 1
    e_sexf_cl        <- -0.128; label("Effect of female sex (vs male) on CL/F (fraction)")      # Xu 2025 Table 1
    e_dis_healthy_vp <-  0.266; label("Effect of healthy status (vs oHCM) on Vp/F (fraction)")  # Xu 2025 Table 1
    e_sexf_vp        <- -0.193; label("Effect of female sex (vs male) on Vp/F (fraction)")      # Xu 2025 Table 1

    # Food effects. e_fed_ka is the effect of the FASTED state and is applied to
    # (1 - FED); the fed / unspecified arm is the reference.
    e_fed_ka             <-  0.170;  label("Effect of fasted dosing (FED = 0, vs fed reference) on Ka (fraction)")  # Xu 2025 Table 1 "Fasted on Ka"
    e_fed_highfat_ka     <- -0.365;  label("Effect of a high-fat meal on Ka (fraction)")                            # Xu 2025 Table 1
    e_fed_highfat_tlag   <-  0.115;  label("Effect of a high-fat meal on Tlag (fraction)")                          # Xu 2025 Table 1
    e_fed_highfat_fdepot <-  0.0687; label("Effect of a high-fat meal on relative bioavailability F1 (fraction)")   # Xu 2025 Table 1

    # Inter-individual variability, reported as omega^2 (variance) on the log
    # scale; the source's %CV is sqrt(exp(omega^2) - 1) * 100 (Xu 2025 Table 1
    # note). IIV on Q/F was fixed to 0 in the source control stream
    # ($OMEGA 0 FIX ; BSVQ/F) and is therefore not carried here.
    etalcl ~ 0.0793 # 28.7% CV  # Xu 2025 Table 1
    etalvc ~ 2.35   # 308% CV   # Xu 2025 Table 1
    etalvp ~ 0.041  # 20.5% CV  # Xu 2025 Table 1
    etalka ~ 0.329  # 62.4% CV  # Xu 2025 Table 1

    # Residual error: additive on the natural-log-transformed concentration
    # (NONMEM Y = LOG(F) + EPS(1)), i.e. log-normal on the linear scale.
    # expSd = sqrt(0.0414) = 0.2035, which the source reports as 20.3% CV.
    expSd <- 0.2035; label("Log-normal residual error SD (log ng/mL)")  # Xu 2025 Table 1
  })
  model({
    # 1. Derived covariate multipliers (proportional form, Xu 2025 Figure S1/S3)
    cl_cov <- (1 + e_dis_healthy_cl * DIS_HEALTHY) * (1 + e_sexf_cl * SEXF)
    vp_cov <- (1 + e_dis_healthy_vp * DIS_HEALTHY) * (1 + e_sexf_vp * SEXF)
    ka_cov <- (1 + e_fed_ka * (1 - FED)) * (1 + e_fed_highfat_ka * FED_HIGHFAT)

    # 2. Individual parameters (body weight normalised to a standard 80 kg)
    cl <- exp(lcl + etalcl) * (WT / 80)^e_wt_cl_q  * cl_cov
    vc <- exp(lvc + etalvc) * (WT / 80)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 80)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 80)^e_wt_vc_vp * vp_cov
    ka <- exp(lka + etalka) * ka_cov

    # Formulation selects the lag time; a high-fat meal lengthens it
    tlag <- (exp(ltlag_tab) * (1 - FORM_CAPSULE) + exp(ltlag_cap) * FORM_CAPSULE) *
      (1 + e_fed_highfat_tlag * FED_HIGHFAT)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. ODE system
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Bioavailability and lag time
    f(depot)    <- exp(lfdepot) * (1 + e_fed_highfat_fdepot * FED_HIGHFAT)
    alag(depot) <- tlag

    # 6. Observation. Dose in mg and vc in L give mg/L = ug/mL; the factor of
    # 1000 converts to ng/mL, matching the source scaling S2 = V2/1000.
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
