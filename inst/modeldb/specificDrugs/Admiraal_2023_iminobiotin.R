Admiraal_2023_iminobiotin <- function() {
  description <- "Two-compartment IV population PK model for 2-iminobiotin (2-IB, a selective neuronal nitric oxide synthase inhibitor) in adults after out-of-hospital cardiac arrest, with a power-model eGFR-on-clearance covariate effect."
  reference <- "Admiraal MM, Velseboer DC, Tjabbes H, Vis P, Peeters-Scholte C, Horn J. Neuroprotection after cardiac arrest with 2-iminobiotin: a single center phase IIa study on safety, tolerability, and pharmacokinetics. Front Neurol. 2023;14:1136046. doi:10.3389/fneur.2023.1136046"
  vignette <- "Admiraal_2023_iminobiotin"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Estimated glomerular filtration rate on hospital admission, by the Modification of Diet in Renal Disease (MDRD) equation: 186 * (Creat/88.4)^-1.154 * Age^-0.203 * (0.742 if female) * (1.210 if Black).",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-model effect on clearance with reference 90 mL/min/1.73 m^2. The reference value is NOT reported in Admiraal 2023 or its supplement; 90 mL/min/1.73 m^2 was set per the nlmixr2lib registry's most common adult-MDRD/CKD-EPI convention (Bajaj 2017, Li 2019), with operator approval (task 171 sidecar 001, response value 'A'). At eGFR = 90 the model returns the published typical CL = 12.10 L/h; behaviour at any non-reference eGFR depends on this REF choice (raised to the 1.03 power exponent). Source column 'eGFR' (MDRD-based) maps to the canonical CRCL.",
      source_name        = "eGFR"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 21L,
    n_studies        = 1L,
    age_range        = "median 60.5 yr (cohort A IQR 59.5-67); 64 yr (cohort B IQR 63-66); 74 yr (cohort C IQR 69-76)",
    age_median       = "65 years (pooled, approximated from cohort medians)",
    weight_range     = "median 87.5 kg (cohort A IQR 82-112); 76.5 kg (cohort B IQR 67.5-96.5); 83 kg (cohort C IQR 77-84)",
    weight_median    = "83 kg (pooled, approximated from cohort medians)",
    sex_female_pct   = 24,
    race_ethnicity   = "Not reported in detail",
    disease_state    = "Adult survivors of out-of-hospital cardiac arrest (OHCA) admitted to the intensive care unit after return of spontaneous circulation; treated with targeted temperature management (36 C for 24 h) and sedated with propofol +/- remifentanil during study drug administration.",
    dose_range       = "Three dose-escalation cohorts: cohort A 0.055 mg/kg/dose (body-weight dosing); cohorts B and C eGFR-based bins targeting AUC0-24h 2,100-3,300 ng*h/mL (cohort B) and 7,200-8,400 ng*h/mL (cohort C). All cohorts received six 15-minute IV infusions, every 4 h, over 24 h.",
    regions          = "Single centre, Amsterdam UMC, Netherlands",
    renal_function   = "MDRD eGFR on hospital admission spanned the dosing-table range of 0-220 mL/min/1.73 m^2 (Supplementary Table S1). Population median eGFR not reported in the paper.",
    notes            = "Demographics from Admiraal 2023 Table 1. N = 21 (cohort A 8, cohort B 8, cohort C 5). One patient (patient 12) was screen-failed post hoc and excluded from PK analysis; the modelled N is therefore 21, not the protocol's planned N = 24. Median time from OHCA to first dose was 5.3 h (IQR 4.8-5.6). NCT02836340 / EUDRACT 2015-003902-17."
  )

  ini({
    # Structural parameters - Admiraal 2023 Supplementary Table S3 (final
    # population PK model fit to pooled cohorts A + B + C, n = 21).
    lcl       <- log(12.10);  label("Typical clearance at reference eGFR (L/h)")    # Supplementary Table S3: Clearance = 12.10 L/h (SE 1.31)
    lvc       <- log(10.1);   label("Central volume of distribution (L)")            # Supplementary Table S3: Vcentral = 10.1 L (SE 1.17)
    lq        <- log(14.6);   label("Intercompartmental clearance (L/h)")            # Supplementary Table S3: Q = 14.6 L/h (SE 1.15)
    lvp       <- log(10.8);   label("Peripheral volume of distribution (L)")          # Supplementary Table S3: Vperipheral = 10.8 L (SE 0.72)

    # Covariate effect: power exponent of eGFR on CL. Reference eGFR is not
    # reported in the paper; set to 90 mL/min/1.73 m^2 in model() per the
    # registry's adult-MDRD/CKD-EPI convention (Bajaj 2017, Li 2019). See
    # covariateData[[CRCL]]$notes and the vignette Assumptions section.
    e_crcl_cl <- 1.03;        label("Power exponent of CRCL (MDRD eGFR) on CL (unitless)")  # Supplementary Table S3: eGFR0 on Clearance = 1.03 (SE 0.43)

    # Inter-individual variability (variances on the log scale, NONMEM $OMEGA).
    # Admiraal 2023 Supplementary Table S3 reports omega^2 = 0.225 on CL and
    # 0.27 on Vcentral; no off-diagonal covariance reported.
    etalcl ~ 0.225                                                                    # Supplementary Table S3: IIV (CL) variance = 0.225 (SE 0.09)
    etalvc ~ 0.27                                                                     # Supplementary Table S3: IIV Vcentral variance = 0.27 (SE 0.10)

    # Residual error. Admiraal 2023 Supplementary Table S3 reports a single
    # residual variance term of 0.06 (SE 0.02), interpreted as the NONMEM
    # $SIGMA variance of the proportional (log-scale additive) error. The
    # nlmixr2 proportional SD is therefore sqrt(0.06) ~ 0.245 (~24.5% CV).
    propSd <- sqrt(0.06); label("Proportional residual error (fraction)")             # Supplementary Table S3: Residual Error variance = 0.06 (SE 0.02)
  })

  model({
    # Individual PK parameters with eGFR-on-CL power-model covariate effect
    # (reference eGFR = 90 mL/min/1.73 m^2; see covariateData[[CRCL]]$notes).
    cl <- exp(lcl + etalcl) * (CRCL / 90)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV model (dosing is into central; no absorption phase)
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # central is in mg, vc in L, so central/vc is in mg/L = 1000 * ng/mL.
    # Multiply by 1000 to express Cc in ng/mL (matches the paper's units).
    Cc <- (central / vc) * 1000

    Cc ~ prop(propSd)
  })
}
