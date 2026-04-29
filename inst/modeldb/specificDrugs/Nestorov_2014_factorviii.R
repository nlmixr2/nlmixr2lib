Nestorov_2014_factorviii <- function() {
  description <- "Two-compartment population PK model for recombinant factor VIII Fc fusion protein (rFVIIIFc, efmoroctocog alfa) in previously treated patients with severe hemophilia A (Nestorov 2014; final covariate model with VWF on CL and WT and HCT on V1)"
  reference <- "Nestorov I, Neelakantan S, Ludden TM, Li S, Jiang H, Rogge M. Population pharmacokinetics of recombinant factor VIII Fc fusion protein. Clin Pharmacol Drug Dev. 2015;4(3):163-174. doi:10.1002/cpdd.167"
  vignette <- "Nestorov_2014_factorviii"
  units <- list(time = "h", dosing = "IU", concentration = "IU/dL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power scaling on V1 with reference weight 73 kg (Nestorov 2014 Eq. 2). Treated as baseline body weight in this model.",
      source_name        = "WT"
    ),
    VWF = list(
      description        = "Plasma von Willebrand factor (VWF) antigen concentration; FVIII-protective carrier protein.",
      units              = "IU/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL with reference 118 IU/dL (Nestorov 2014 Eq. 2). Negative exponent: higher VWF protects FVIII activity from clearance. Time-varying within an individual; the model treats VWF as the covariate value at the observation time; for simulation a typical baseline value (e.g., 100-118 IU/dL) is normally used because the within-subject VWF time course is not characterized in this paper.",
      source_name        = "VWF"
    ),
    HCT = list(
      description        = "Hematocrit",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on V1 with reference 45 % (Nestorov 2014 Eq. 2). Negative exponent: a higher hematocrit reduces plasma fraction of total volume, which decreases the plasma-restricted V1 of FVIII activity.",
      source_name        = "HCT"
    )
  )

  population <- list(
    n_subjects     = 180L,
    n_studies      = 2L,
    age_range      = "phase 3 enrolled subjects >=12 years; phase 1/2a enrolled adults",
    age_median     = "not reported in source; mixed adolescent / adult cohort",
    weight_range   = "not reported as a range; reference weight 73 kg used in V1 allometry (Nestorov 2014 Eq. 2)",
    weight_median  = "not reported; allometric reference 73 kg",
    sex_female_pct = 0,
    race_ethnicity = "not reported in source",
    disease_state  = "Severe hemophilia A (FVIII activity <1 IU/dL); previously treated patients",
    dose_range     = "Single intravenous injection: 25, 50, or 65 IU/kg in phase 1/2a (Cohorts A and B) and phase 3 (Arms 1, 2, 3); Arm 1 prophylaxis individualized over 25-65 IU/kg every 3-5 days; Arm 2 fixed 65 IU/kg weekly; Arm 3 episodic 10-50 IU/kg",
    regions        = "phase 3 was multinational, multicenter; specific regions not enumerated in source",
    co_medication  = "rFVIII (Advate) administered as comparator in phase 1/2a and phase 3 Arm 1 sequential PK subgroup; analysed in a separate population PK model not implemented in this file",
    notes          = "n = 16 in phase 1/2a (Cohort A 25 IU/kg, Cohort B 65 IU/kg) and n = 164 in phase 3 (Arms 1-3). Hemophilia A is X-linked recessive and the registered phase 1/2a and phase 3 studies enrolled male patients; sex_female_pct = 0 reflects this. Baseline demographics summarised from the Methods of Nestorov 2014; the publication does not present a tabulated summary of baseline demographics for the population PK analysis cohort."
  )

  ini({
    # Structural parameters at the paper's reference covariate values (WT = 73 kg,
    # VWF = 118 IU/dL, HCT = 45 %). Volumes in dL and clearances in dL/h match the
    # paper's units (Table 1, "Final Model rFVIIIFc, Without BLQ values" column).
    lcl <- log(1.65);   label("Typical clearance TVCL at reference VWF (dL/h)")                              # Nestorov 2014 Table 1: CL = 1.65 dL/h (final model, without BLQ)
    lvc <- log(37.5);   label("Typical central volume TVV1 at reference WT and HCT (dL)")                    # Nestorov 2014 Table 1: V1 = 37.5 dL
    lq  <- log(0.0746); label("Intercompartmental clearance Q (dL/h)")                                       # Nestorov 2014 Table 1: Q = 0.0746 dL/h
    lvp <- log(6.92);   label("Peripheral volume V2 (dL)")                                                   # Nestorov 2014 Table 1: V2 = 6.92 dL

    # Covariate effects (power form; Nestorov 2014 Eq. 2)
    e_vwf_cl <- -0.343; label("Power exponent of VWF on CL (unitless)")  # Nestorov 2014 Table 1: Exponent on VWF = -0.343
    e_wt_vc  <-  0.382; label("Allometric exponent of WT on V1 (unitless)")  # Nestorov 2014 Table 1: Allometric exponent on V1 = 0.382
    e_hct_vc <- -0.419; label("Power exponent of HCT on V1 (unitless)")  # Nestorov 2014 Table 1: Exponent on HCT = -0.419

    # Inter-individual variability: BLOCK(2) on CL and V1 (Nestorov 2014 Table 1).
    # Variances: omega^2 = (CV%/100)^2 reproduces the paper's reported NONMEM
    # %CV (24.3% for CL, 13.4% for V1) under the convention CV% = sqrt(omega^2)*100.
    # Covariance v12 = 0.0179 is the paper's reported population mean; corresponds
    # to a correlation of 0.55, consistent with the Table 1 value of 0.548.
    # Inter-occasion variability (also reported in Table 1, 20.6% on CL and 12.0%
    # on V1 with correlation 0.639) is not implemented; see vignette deviations.
    etalcl + etalvc ~ c(0.059049,
                        0.0179,   0.017956)  # Nestorov 2014 Table 1: IIV CL 24.3% (omega^2 = 0.0590), v12 = 0.0179, IIV V1 13.4% (omega^2 = 0.0180)

    # Residual error (final rFVIIIFc model). The paper estimates a common
    # proportional error and study-specific additive errors; the phase 3 additive
    # value (0.208 IU/dL) is used here as the representative residual because the
    # phase 3 study contributed the majority (164 / 180) of subjects. The phase
    # 1/2a additive value of 0.421 IU/dL is documented in the vignette deviations.
    propSd <- 0.136; label("Proportional residual error (fraction)")    # Nestorov 2014 Table 1: proportional error = 13.6%
    addSd  <- 0.208; label("Additive residual error (IU/dL)")           # Nestorov 2014 Table 1: phase 3 additive error = 0.208 IU/dL
  })
  model({
    # Individual PK parameters with covariate effects (Nestorov 2014 Eq. 2).
    cl <- exp(lcl + etalcl) * (VWF / 118)^e_vwf_cl
    vc <- exp(lvc + etalvc) * (WT  /  73)^e_wt_vc * (HCT / 45)^e_hct_vc
    q  <- exp(lq)
    vp <- exp(lvp)

    # Micro-constants for explicit two-compartment ODEs.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODEs: rFVIIIFc is administered as an intravenous bolus / short infusion;
    # doses enter the central compartment directly. Activity (FVIII:C) circulates
    # in plasma so V1 approximates plasma volume.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Observation: FVIII activity in IU/dL (1 IU/dL = 1% of pooled normal plasma).
    # Dose in IU and V1 in dL give central / vc in IU/dL directly.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
