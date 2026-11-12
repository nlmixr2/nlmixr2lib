Zhou_2026_gabapentin <- function() {
  description <- "Joint population PK/PD model of oral immediate-release gabapentin and Visual Analog Scale (VAS, 0-10) pain score in adults with chronic neuropathic pain (Zhou 2026). Two-compartment PK with first-order absorption, an absorption lag time and linear elimination; apparent clearance scales with CKD-EPI 2021 eGFR as a power function centred on 84.8532 mL/min/1.73 m^2. Pain score follows a full inhibitory Imax model (Imax = 1) driven by an effect compartment with separate first-order input (k1e, from plasma) and output (ke0) rate constants, dCe/dt = k1e * Cc - ke0 * Ce; carriers of the SLC22A2 (OCT2) c.808G>T variant (heterozygous GT) have an 11.5-fold lower k1e. Baseline pain is logit-normal on (0, 10.1). Concentrations in mg/L (= ug/mL) for doses in mg."
  reference <- "Zhou L, Yamamoto PA, Walker M, Costa ACC, Lauretti GR, Dach F, Schmidt S, de Moraes NV. Gabapentin CNS exposure and analgesic response are modulated by OCT2 genotype in patients with chronic neuropathic pain. Front Pharmacol. 2026;17:1760901. doi:10.3389/fphar.2026.1760901"
  vignette <- "Zhou_2026_gabapentin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description = "Estimated glomerular filtration rate by the 2021 CKD-EPI creatinine equation, BSA-normalized",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL/F: (CRCL / 84.8532)^1.34. The reference 84.8532 mL/min/1.73 m^2 is the observation-weighted mean eGFR (Zhou 2026 Methods quotes 84.85; the Simulx code in Supplementary Data 2 uses 84.8532). Missing eGFR (12% of subjects) was imputed with the median during estimation. Cohort mean (SD) 90.29 (22.79) mL/min/1.73 m^2 (Supplementary Table S1).",
      source_name = "EGFR_CKD"
    ),
    SNP_SLC22A2_808GT = list(
      description = "SLC22A2 (OCT2) c.808G>T variant carrier indicator, 1 = GT, 0 = GG",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (homozygous wild-type 808GG)",
      notes = "Time-fixed germline genotype. Enters the effect-compartment input rate constant multiplicatively as exp(-2.44 * SNP_SLC22A2_808GT) (Zhou 2026 Table 1 'OCT2 on ke1'; Methods: GG coded 0, GT coded 1). The cohort contained 79 GG and 15 GT subjects and no 808TT homozygotes (Supplementary Table S1), so the effect is estimated for heterozygotes only; applying it to TT homozygotes is an extrapolation.",
      source_name = "OCT2"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened (Zhou 2026 Methods, Population pharmacokinetic analysis); not retained.",
      source_name = "SEX"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained.",
      source_name = "AGE"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained.",
      source_name = "WT"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained.",
      source_name = "BMI"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained (eGFR was retained instead).",
      source_name = "SCR"
    ),
    SNP_SLC22A4_1507CT = list(
      description = "SLC22A4 (OCTN1) c.1507C>T genotype",
      units = "(categorical)",
      type = "categorical",
      reference_category = "0 (CC)",
      notes = "Screened on the PK and PD parameters, coded 0 = CC, 1 = CT, 2 = TT (Zhou 2026 Methods); not retained. Cohort 42 CC, 39 CT, 13 TT (Supplementary Table S1).",
      source_name = "OCTN1"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "gabapentin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "gabapentin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "gabapentin", units = "mg", specimen = "plasma", verified = TRUE),
    effect = list(
      analyte = "gabapentin, effect site (scaled by k1e / ke0 relative to plasma)",
      units = "mg/L",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 94L,
    n_studies = 2L,
    n_observations = "568 plasma gabapentin concentrations (9.15% below the 0.2 ug/mL LLOQ, handled as left-censored) and 392 VAS pain scores",
    age_mean = "53 years (SD 11)",
    weight_mean = "81.97 kg (SD 17.34)",
    bmi_mean = "29.9 kg/m^2 (SD 6.15)",
    sex_female_pct = 56.4,
    race_ethnicity = "Brazilian (Latin American) cohort; race not tabulated",
    disease_state = "Chronic neuropathic pain (diabetic neuropathy or compressive etiologies such as cervical or lumbar disc herniation; VAS above 4 in Study 1)",
    renal_function = "Mean eGFR (CKD-EPI 2021) 90.29 (SD 22.79) mL/min/1.73 m^2: 49 normal or high, 22 mild, 11 mild-to-severe, 1 severe, 11 missing",
    genotype = "SLC22A2 c.808G>T: 79 GG, 15 GT, 0 TT; SLC22A4 c.1507C>T: 42 CC, 39 CT, 13 TT",
    dose_range = "Study 1: single oral 300 mg IR capsule (n = 29, rich sampling); Study 2: multiple oral IR doses of 600-3600 mg/day (n = 65, sparse sampling)",
    regions = "Brazil (Clinical Hospital of the University of Sao Paulo, Ribeirao Preto)",
    notes = "Demographics from Zhou 2026 Methods 'Study population' and Supplementary Table S1. Joint PK/PD estimation in Monolix 2024R1 (SAEM)."
  )

  ini({
    # PK (Zhou 2026 Table 1, final joint PK/PD estimates)
    ltlag <- log(0.34); label("Absorption lag time, tlag (h)") # Table 1 'tlag (h)' = 0.34 (RSE 8.03%)
    lka <- log(0.14); label("First-order absorption rate constant, ka (1/h)") # Table 1 'ka (h-1)' = 0.14 (RSE 7.04%)
    lcl <- log(10.16); label("Apparent clearance at eGFR 84.8532 mL/min/1.73 m^2, CL/F (L/h)") # Table 1 'CL/F (L/h)' = 10.16 (RSE 12.15%)
    e_crcl_cl <- 1.34; label("Power exponent of eGFR / 84.8532 on CL/F (unitless)") # Table 1 'logt (eGFR/84.85) on CL/F' = 1.34 (RSE 24.38%)
    lvc <- log(18.16); label("Apparent central volume, V1/F (L)") # Table 1 'V1/F (L)' = 18.16 (RSE 15.13%)
    lq <- log(6.58); label("Apparent intercompartmental clearance, Q/F (L/h)") # Table 1 'Q (L/h)' = 6.58 (RSE 18.68%)
    lvp <- log(357.67); label("Apparent peripheral volume, V2/F (L)") # Table 1 'V2/F (L)' = 357.67 (RSE 31.27%)

    # PD: effect compartment and Imax model of the VAS pain score
    lk1e <- log(0.53); label("Effect-compartment input rate constant from plasma, ke1, in 808GG (1/h)") # Table 1 'ke1 (h-1)' = 0.53 (RSE 55.76%)
    e_snp_oct2_k1e <- -2.44; label("Log-scale effect of SLC22A2 c.808GT on ke1 (unitless)") # Table 1 'OCT2 on ke1' = -2.44 (RSE 41.81%); exp(-2.44) = 0.087
    lke0 <- fixed(log(1)); label("Effect-compartment output rate constant, ke2 (1/h)") # Table 1 'ke2 (h-1)' = 1 (fixed)
    logite0 <- qlogis(7.31 / 10.1); label("Logit of the baseline VAS pain score E0 on (0, 10.1) (unitless)") # Table 1 'E0' = 7.31 (RSE 5.86%); Supplementary Data 2 E0 logitNormal max = 10.1
    lic50 <- fixed(log(0.26311)); label("Effect-site concentration at half-maximal pain inhibition, IC50 (mg/L)") # Table 1 'IC50 (ng/mL)' = 263.11 (fixed) = 0.26311 mg/L
    limax <- fixed(log(1)); label("Maximum fractional inhibition of the pain score, Imax (unitless)") # Table 1 'Imax' = 1 (fixed)

    # IIV: Table 1 reports omega as the SD of the random effects; variances are SD^2
    etaltlag ~ 0.2704 # Table 1 'omega Tlag' = 0.52 (SD); 0.52^2
    etalka ~ 0.0225 # Table 1 'omega ka' = 0.15 (SD); 0.15^2
    etalcl ~ 0.2116 # Table 1 'omega CL' = 0.46 (SD); 0.46^2
    etalvc ~ 0.2916 # Table 1 'omega V1' = 0.54 (SD); 0.54^2
    etalq ~ 0.4096 # Table 1 'omega Q' = 0.64 (SD); 0.64^2
    etalk1e ~ 4.4521 # Table 1 'omega ke1' = 2.11 (SD); 2.11^2
    etalogite0 ~ 0.9604 # Table 1 'omega E0' = 0.98 (SD, logit scale); 0.98^2

    # Residual error. PK: Monolix logNormal distribution with constant error,
    # i.e. additive on log(Cc). PD: Monolix logitNormal(min -0.2, max 10.2) with
    # combined1 error, i.e. logit(y) = logit(E) + (a + b * logit(E)) * eps.
    expSd <- 0.19; label("Additive residual SD of log plasma concentration (unitless)") # Table 1 PK 'Additive (a)' = 0.19 (RSE 4.44%)
    addSd_score <- 0.52; label("Additive residual SD of the pain score on the logit scale (unitless)") # Table 1 PD 'Additive (a)' = 0.52 (RSE 5.81%)
    propSd_score <- 0.30; label("Proportional residual SD of the pain score on the logit scale (unitless)") # Table 1 PD 'Proportional (b)' = 0.30 (RSE 9.63%)
  })

  model({
    tlag <- exp(ltlag + etaltlag)
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (CRCL / 84.8532)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp)

    k1e <- exp(lk1e + etalk1e + e_snp_oct2_k1e * SNP_SLC22A2_808GT)
    ke0 <- exp(lke0)
    e0 <- expit(logite0 + etalogite0, 0, 10.1)
    ic50 <- exp(lic50)
    imax <- exp(limax)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag

    Cc <- central / vc

    # Supplementary Data 1: ddt_Ce = ke1 * Cc - ke2 * Ce
    d/dt(effect) <- k1e * Cc - ke0 * effect

    # Supplementary Data 1: E = E0 * (1 - Ce / (Ce + IC50)); Imax = 1
    score <- e0 * (1 - imax * effect / (effect + ic50))

    Cc ~ lnorm(expSd)
    score ~ logitNorm(addSd_score, -0.2, 10.2) + propT(propSd_score) + combined1()
  })
}
