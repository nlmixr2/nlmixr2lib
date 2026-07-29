Hong_2006_amphotericinB_liposomal <- function() {
  description <- "Two-compartment population PK model for liposomal amphotericin B (AmBisome) in 39 pediatric oncology patients receiving 1-h IV infusions (Hong 2006). Clearance and central volume scale exponentially with body weight centered at the cohort-median 21 kg; the paper additionally reports substantial between-occasion variability on CL and V1 that is encoded here as IIV on fixed-at-1 multiplicative anchors (Bellanti 2015 IOV-as-IIV pattern)."
  reference <- "Hong Y, Shaw PJ, Nath CE, Yadav SP, Stephen KR, Earl JW, McLachlan AJ. Population pharmacokinetics of liposomal amphotericin B in pediatric patients with malignant diseases. Antimicrob Agents Chemother. 2006 Mar;50(3):935-942. doi:10.1128/AAC.50.3.935-942.2006"
  vignette <- "Hong_2006_amphotericinB_liposomal"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying across dosing courses. Reference value 21 kg (cohort median, Hong 2006 Table 1). Covariate enters as an exponential effect on CL and V1: parameter = typical * exp(coefficient * (WT - 21)). The cohort included patients between 6.1 and 84.1 kg.",
      source_name        = "WT"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 39L,
    n_studies        = 1L,
    age_range        = "0.17-17 years",
    age_median       = "6.5 years (mean 7.1, SD 5.1)",
    weight_range     = "6.1-84.1 kg",
    weight_median    = "21.1 kg (mean 28.8, SD 19.8)",
    sex_female_pct   = 33,
    race_ethnicity   = "Not reported (single-centre paediatric oncology cohort in Sydney, Australia)",
    disease_state    = "Paediatric oncology patients receiving liposomal amphotericin B (AmBisome) for antifungal therapy or prophylaxis; all neutropenic and febrile. Diagnoses: 15 acute lymphoblastic leukaemia, 4 acute myeloid leukaemia, 20 other; 16 of 39 had received a bone marrow transplant (12 allogeneic, 4 autologous).",
    dose_range       = "0.8-5.9 mg/kg/day liposomal amphotericin B administered as a 1-h IV infusion through a central venous catheter, once daily. Typical prophylaxis 1-3 mg/kg/day; escalated to 5-6 mg/kg/day on evidence of invasive fungal infection.",
    regions          = "Australia (single centre, the Children's Hospital at Westmead, Sydney, NSW)",
    n_concentrations = 637L,
    n_courses        = 48L,
    notes            = "Baseline characteristics from Hong 2006 Table 1 (mean (SD), median, range): age 7.1 yr (5.1), 6.5, 0.17-17; weight 28.8 kg (19.8), 21.1, 6.1-84.1; height 119.2 cm (33.7), 118.5, 61.5-190; BSA 0.96 m^2 (0.46), 0.82, 0.32-2.11; creatinine 74 umol/L (108), 47, 19-659; GGT 105 U/L (211), 39, 11-1398. 637 plasma concentration observations from 48 dosing courses across 39 patients (mean 13.3 samples per course, range 3-30). NONMEM 5.1.1 first-order conditional estimation; bootstrap n=1000 in Wings for NONMEM. Renal-function markers were not screened as covariates because renal elimination is a minor pathway for L-AmB; sex did not significantly affect PK by Mann-Whitney U; weight and height were highly correlated (r2 = 0.93) and weight was selected as the size descriptor over height; age was tested on V1 but did not improve the fit beyond weight. Hong 2006 reports between-occasion variability (BOV, 46% on CL and 56% on V1) that exceeds the within-subject IIV on the same parameters. nlmixr2lib has no canonical pattern for occasion-keyed IOV without an OCC column, so this BOV is encoded here as additional log-normal IIV on fixed-at-1 multipliers (lcl_bov_mult, lvc_bov_mult) following the Bellanti 2015 deferoxamine precedent; see vignette Assumptions and deviations for the loss of within-subject occasion drift."
  )

  ini({
    # Structural parameters (Hong 2006 Table 3, final/full model). Reference subject is 21 kg.
    lcl <- log(0.44); label("Clearance at WT=21 kg (L/h)")          # Hong 2006 Table 3 row CL (RSE 27%)
    lvc <- log(3.12); label("Central volume V1 at WT=21 kg (L)")    # Hong 2006 Table 3 row V1 (RSE 40%)
    lq  <- log(0.73); label("Intercompartmental clearance Q (L/h)") # Hong 2006 Table 3 row Q (RSE 18%)
    lvp <- log(18.0); label("Peripheral volume V2 (L)")             # Hong 2006 Table 3 row V2 (RSE 40%)

    # Covariate effects: CL and V1 follow exp(coefficient * (WT - 21)) per the full
    # model expressions in Hong 2006 (text following Table 3 and the "TABLE 2 model 9"
    # entry: CL = theta_1 * EXP(theta_5 * (WT - 21)) and V1 = theta_2 * EXP(theta_6 *
    # (WT - 21))). The coefficient carries units of 1/kg.
    e_wt_cl <- 0.0152; label("Exponential coefficient on (WT-21) for CL (1/kg)")  # Hong 2006 Table 3 row theta_5 (RSE 12%)
    e_wt_vc <- 0.0241; label("Exponential coefficient on (WT-21) for V1 (1/kg)")  # Hong 2006 Table 3 row theta_6 (RSE 29%)

    # Inter-individual variability. The "%" values in Hong 2006 Table 3 are the
    # log-scale SD (omega) expressed as a percent (so omega^2 ~ ("%"/100)^2);
    # the variance scale matches the omega^2 values in Hong 2006 Table 4.
    etalcl ~ 0.0108  # Hong 2006 Table 3 CL IIV 10%; matches Table 4 omega^2 CL = 0.0108
    # IIV on V1 failed to converge in the original data set (Hong 2006 Table 3 leaves the
    # IIV-on-V1 column blank; the Discussion notes that NONMEM struggled to estimate it
    # because BOV greatly exceeded IIV). The bootstrap mean from 1000 successful runs
    # (Table 4 omega^2 V1 = 0.0450) is used here as the best available point estimate.
    etalvc ~ 0.0450  # Hong 2006 Table 4 bootstrap omega^2 V1 = 0.0450
    etalq  ~ 0.593   # Hong 2006 Table 3 Q IIV 77%; matches Table 4 omega^2 Q = 0.593
    etalvp ~ 0.546   # Hong 2006 Table 3 V2 IIV 74%; matches Table 4 omega^2 V2 = 0.546

    # Between-occasion variability (BOV) encoded as IIV on fixed-at-1 multipliers
    # (Bellanti 2015 deferoxamine pattern). Hong 2006 estimated BOV on CL and V1 via
    # NONMEM's BLOCK(2) SAME (single magnitude shared across all occasions). The
    # multiplicative anchors lcl_bov_mult / lvc_bov_mult have typical value
    # exp(log(1)) = 1 and carry the paper's BOV variance as ordinary IIV; this
    # preserves the population-level CL / V1 spread but loses the within-subject
    # between-occasion drift the paper estimates. See vignette Assumptions and
    # deviations.
    lcl_bov_mult <- fixed(log(1))
    label("Multiplicative CL factor anchor (log scale, fixed; carries IIV originally reported as BOV)")
    etalcl_bov_mult ~ 0.209  # Hong 2006 Table 3 CL BOV 46%; matches Table 4 pi^2 CL = 0.209

    lvc_bov_mult <- fixed(log(1))
    label("Multiplicative V1 factor anchor (log scale, fixed; carries IIV originally reported as BOV)")
    etalvc_bov_mult ~ 0.319  # Hong 2006 Table 3 V1 BOV 56%; matches Table 4 pi^2 V1 = 0.319

    # Combined additive + proportional residual error (Hong 2006 Table 3).
    propSd <- 0.27;  label("Proportional residual error (fraction)")  # Hong 2006 Table 3 proportional 27% (RSE 7%); matches Table 4 sigma^2 prop = 0.0709
    addSd  <- 0.034; label("Additive residual error (mg/L)")           # Hong 2006 Table 3 additive 0.034 mg/L (RSE 24%); matches sqrt of Table 4 sigma^2 add = 0.0011
  })
  model({
    # Multiplicative BOV-anchor factors (typical value 1; carry the paper's BOV
    # magnitude as IIV). Total log-scale random effect on CL is etalcl + etalcl_bov_mult
    # with variances 0.0108 + 0.209 = 0.220 (CL) and 0.045 + 0.319 = 0.364 (V1).
    cl_bov_mult <- exp(lcl_bov_mult + etalcl_bov_mult)
    vc_bov_mult <- exp(lvc_bov_mult + etalvc_bov_mult)

    # Individual PK parameters. CL and V1 scale exponentially with (WT - 21) per the
    # Hong 2006 full model. Q and V2 carry no covariates (the paper found no
    # physiological / clinical rationale to screen them).
    cl <- exp(lcl + e_wt_cl * (WT - 21) + etalcl) * cl_bov_mult
    vc <- exp(lvc + e_wt_vc * (WT - 21) + etalvc) * vc_bov_mult
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV with zero-order input. The 1-h infusion is supplied via the
    # dataset RATE / DUR fields; the model treats the dose as entering the central
    # compartment directly (no depot).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> central/vc has units mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
