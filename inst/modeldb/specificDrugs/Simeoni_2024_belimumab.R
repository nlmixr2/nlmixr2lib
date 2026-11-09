Simeoni_2024_belimumab <- function() {
  description <- "Linear two-compartment IV population PK model for belimumab in adults with active lupus nephritis (BLISS-LN), with baseline fat-free mass allometry on all disposition parameters and time-varying proteinuria (urine protein/creatinine ratio) and serum albumin on clearance (Simeoni 2024); used by Dimelow 2024 to extrapolate exposure to children aged 5-17 years with lupus nephritis"
  reference <- paste(
    "Simeoni M, Yang S, Tompson DJ, Dimelow R. Longitudinal modeling of efficacy",
    "response in patients with lupus nephritis receiving belimumab.",
    "J Pharmacokinet Pharmacodyn. 2024;51(3):289-301. doi:10.1007/s10928-024-09907-w",
    "(PK model parameters: Online Resource 7).",
    "Application: Dimelow R, Liefaard L, Green Y, Tomlinson R. Extrapolation of the",
    "Efficacy and Pharmacokinetics of Belimumab to Support its Use in Children with",
    "Lupus Nephritis. Clin Pharmacokinet. 2024;63(9):1313-1326.",
    "doi:10.1007/s40262-024-01422-y",
    sep = " "
  )
  vignette <- "Dimelow_2024_belimumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "belimumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "belimumab", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    FFM = list(
      description = "Baseline fat-free mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed at baseline (FFMBL in Simeoni 2024 Online Resource 7). Power effects (FFM / 38.4)^0.624 on CL and Q and (FFM / 38.4)^0.723 on V1 and V2; 38.4 kg is the BLISS-LN median (Simeoni 2024 Table 1). Derived from body weight, BMI and sex with the Janmahasatian et al. (Clin Pharmacokinet 2005;44:1051-1065) formula (Dimelow 2024 Methods, ref. 25).",
      source_name = "FFMBL"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "TIME-VARYING. Power effect (ALB / 42)^(theta1 / (1 + theta2 * UPCR)) on CL, with an exponent that shrinks as proteinuria rises (Simeoni 2024 Online Resource 7). BLISS-LN baseline median 31 g/L (Simeoni 2024 Table 1).",
      source_name = "ALB"
    ),
    UPCR = list(
      description = "Urine protein-to-creatinine ratio (proteinuria)",
      units = "g/g",
      type = "continuous",
      reference_category = NULL,
      notes = "TIME-VARYING, measured every 4 weeks in BLISS-LN. Linear effect (1 + 0.0663 * UPCR) on CL and modifier of the albumin exponent (Simeoni 2024 Online Resource 7, where it is called PROT, 'proteinuria level (g/g)'). Dimelow 2024 ESM Section S5 labels PROT 'g/day'; the primary source's g/g (the uPCR that BLISS-LN measured) is used. BLISS-LN baseline median 2.61 g/g in the belimumab arm (Simeoni 2024 Table 1).",
      source_name = "PROT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 224,
    n_studies = 1,
    age_range = "18-63 years",
    age_median = "31 years",
    weight_range = "36.9-136.9 kg",
    weight_median = "59.0 kg",
    ffm_median = "38.1 kg (range 26.3-77.4)",
    sex_female_pct = 88,
    race_ethnicity = "Asian 51%, White 33%, Black African ancestry 13%, American Indian 2%, multiple 1%",
    disease_state = "Active, biopsy-confirmed lupus nephritis (adults), on standard therapy (mycophenolate mofetil or cyclophosphamide induction + glucocorticoids)",
    dose_range = "Belimumab 10 mg/kg IV on days 0, 14, 28 and every 28 days thereafter, for 104 weeks",
    regions = "Multinational",
    renal_function = "Baseline eGFR median 99 (range 14-208) mL/min/1.73 m^2; uPCR median 2.61 (0.18-16.57) g/g; serum albumin median 31 (16-42) g/L",
    studies = "BEL114054 (BLISS-LN; NCT01639339), belimumab arm",
    notes = "Belimumab-arm demographics from Simeoni 2024 Table 1. PK samples on days 0, 3, 14, 28, 56, 168, 171, 364 and 728 (Simeoni 2024 Methods). Dimelow 2024 applies the model unchanged to children aged 5-17 years, sampling body size from NHANES (90% female, weight and BMI scaled to 86% and 94%) and proteinuria / albumin time courses from adult BLISS-LN patients."
  )

  ini({
    # Structural parameters -- Simeoni 2024 Online Resource 7, published in
    # mL/day and mL and divided by 1000 here to give L/day and L.
    lcl <- log(0.175); label("Clearance at FFM 38.4 kg, albumin 42 g/L and zero proteinuria (L/day)") # Online Resource 7: CL = 175 mL/day (RSE 2.7%)
    lvc <- log(2.728); label("Central volume of distribution at FFM 38.4 kg (L)") # Online Resource 7: V1 = 2728 mL (RSE 1.8%)
    lq <- log(0.487); label("Intercompartmental clearance at FFM 38.4 kg (L/day)") # Online Resource 7: Q = 487 mL/day (RSE 6.3%)
    lvp <- log(1.992); label("Peripheral volume of distribution at FFM 38.4 kg (L)") # Online Resource 7: V2 = 1992 mL (RSE 4.99%)

    # Covariate effects -- Simeoni 2024 Online Resource 7
    e_ffm_cl <- 0.624; label("Power exponent of baseline FFM on CL and Q (unitless)") # Online Resource 7: x (FFMBL/38.4)^theta on CL = 0.624; Q row 'See FFMBL on CL'
    e_ffm_vc <- 0.723; label("Power exponent of baseline FFM on V1 and V2 (unitless)") # Online Resource 7: x (FFMBL/38.4)^theta on V1 = 0.723; V2 row 'See FFMBL on V1'
    e_alb_cl <- -1.74; label("Power exponent of serum albumin on CL at zero proteinuria (unitless)") # Online Resource 7: theta1 = -1.74 in (ALB/42)^(theta1/(1+theta2*PROT))
    e_upcr_alb_cl <- 0.0832; label("Attenuation of the albumin exponent on CL per unit uPCR (g/g)^-1") # Online Resource 7: theta2 = 0.0832 in (ALB/42)^(theta1/(1+theta2*PROT))
    e_upcr_cl <- 0.0663; label("Linear increase in CL per unit uPCR (g/g)^-1") # Online Resource 7: x (1 + theta*PROT) = 0.0663

    # IIV -- log-scale variances, Simeoni 2024 Online Resource 7 (no
    # covariances reported; none on Q)
    etalcl ~ 0.0593 # Online Resource 7: w2 CL = 0.0593 (CV 24.7%)
    etalvc ~ 0.0322 # Online Resource 7: w2 V1 = 0.0322 (CV 18.1%)
    etalvp ~ 0.133 # Online Resource 7: w2 V2 = 0.133 (CV 37.7%)
    etapropSd ~ 0.346 # Online Resource 7: w2 PROP = 0.346 (CV 64.3%), between-patient variance on the proportional residual error

    # Residual error -- Simeoni 2024 Online Resource 7
    propSd <- 0.241; label("Proportional residual error before the per-patient modifier (fraction)") # Online Resource 7: PROP = 0.241 (RSE 4.5%)
    addSd <- fixed(0.1); label("Additive residual error (mg/L)") # Online Resource 7: ADD = 0.1 ug/mL (Fixed)
  })
  model({
    # Covariate effects (Simeoni 2024 Online Resource 7; also printed in
    # Dimelow 2024 ESM Section S5)
    ffm_cl <- (FFM / 38.4)^e_ffm_cl
    ffm_vc <- (FFM / 38.4)^e_ffm_vc
    alb_cl <- (ALB / 42)^(e_alb_cl / (1 + e_upcr_alb_cl * UPCR))
    upcr_cl <- 1 + e_upcr_cl * UPCR

    cl <- exp(lcl + etalcl) * ffm_cl * alb_cl * upcr_cl
    vc <- exp(lvc + etalvc) * ffm_vc
    q <- exp(lq) * ffm_cl
    vp <- exp(lvp + etalvp) * ffm_vc

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Linear two-compartment IV model
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Serum belimumab concentration: mg / L = ug/mL
    Cc <- central / vc

    # Combined proportional + additive error; the proportional SD carries
    # between-patient variability (w2 PROP in Online Resource 7)
    propSd_i <- propSd * exp(etapropSd)
    Cc ~ add(addSd) + prop(propSd_i) + combined2()
  })
}
