Aoyama_2012_sepantronium <- function() {
  description <- "One-compartment IV population PK model for sepantronium bromide (YM155), a small-molecule survivin suppressant administered as a 7-day continuous IV infusion every 21 days, with power-form covariate effects of creatinine clearance and alanine aminotransferase and proportional cancer-type effects (hormone-refractory prostate cancer and melanoma vs non-small cell lung cancer) on clearance, in adults with NSCLC, HRPC, or unresectable stage III/IV melanoma (Aoyama 2012)"
  reference <- "Aoyama Y, Kaibara A, Takada A, Nishimura T, Katashima M, Sawamoto T. Population pharmacokinetic modeling of sepantronium bromide (YM155), a small molecule survivin suppressant, in patients with non-small cell lung cancer, hormone refractory prostate cancer, or unresectable stage III or IV melanoma. Invest New Drugs. 2013;31(2):443-451. doi:10.1007/s10637-012-9867-x"
  vignette <- "Aoyama_2012_sepantronium"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLCR. Computed by the Cockcroft-Gault equation (paper Eq. 5) in raw mL/min (NOT BSA-normalized to mL/min/1.73 m^2). Stored under the canonical CRCL column per inst/references/covariate-columns.md, following the raw-Cockcroft-Gault pattern of Delattre 2010 amikacin. Reference value 79 mL/min (population median, Aoyama 2012 Table 2). Effect form: (CRCL / 79)^0.425 multiplicative on CL.",
      source_name        = "CLCR"
    ),
    ALT = list(
      description        = "Alanine aminotransferase activity at baseline",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference value 20 U/L (population median, Aoyama 2012 Table 2). Effect form: (ALT / 20)^0.124 multiplicative on CL. The paper notes the exponent is close to 0 so the practical influence is small; the covariate was nonetheless retained in the final model per the backward-elimination likelihood-ratio test.",
      source_name        = "ALT"
    ),
    TUMTP_HRPC = list(
      description        = "Hormone-refractory prostate cancer tumor-type indicator (historical term; modern equivalent castration-resistant prostate cancer / CRPC)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-HRPC; NSCLC is the model's implicit reference when paired with TUMTP_MEL = 0)",
      notes              = "Proportional change on CL: ratio 0.955 (i.e. -4.5% vs NSCLC) per Aoyama 2012 Table 3. The paper uses the power form THETA_HRPC^I_HRPC; for a binary indicator this equals the proportional form (1 + (THETA_HRPC - 1) * I_HRPC) with coefficient -0.045. Cancer type was tested as a three-level categorical (NSCLC reference, HRPC, MM = melanoma); the implicit NSCLC reference holds when both TUMTP_HRPC = 0 and TUMTP_MEL = 0. Decompose a source TUMTP column into TUMTP_HRPC = as.integer(TUMTP == 'HRPC').",
      source_name        = "TUMTP (HRPC level)"
    ),
    TUMTP_MEL = list(
      description        = "Melanoma tumor-type indicator (unresectable stage III/IV cutaneous melanoma)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-melanoma; NSCLC is the model's implicit reference when paired with TUMTP_HRPC = 0)",
      notes              = "Proportional change on CL: ratio 1.24 (i.e. +24% vs NSCLC) per Aoyama 2012 Table 3. The paper's 'MM' label denotes malignant melanoma (unresectable stage III/IV cutaneous melanoma) and is NOT multiple myeloma; do not confuse with the canonical MM register entry. The paper uses the power form THETA_MM^I_MM; for a binary indicator this equals the proportional form (1 + (THETA_MM - 1) * I_MM) with coefficient +0.24. Decompose a source TUMTP column into TUMTP_MEL = as.integer(TUMTP == 'MM').",
      source_name        = "TUMTP (MM = melanoma level)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 96L,
    n_studies       = 3L,
    age_range       = "29-90 years",
    age_median      = "64 years",
    weight_range    = "50-114 kg",
    weight_median   = "81 kg",
    sex_female_pct  = 17.7,
    race_ethnicity  = c(
      `African American` = 4.2,
      `Caucasian (unspecified Hispanic status)` = 45.8,
      `Caucasian Hispanic or Latino` = 8.3,
      `Caucasian non-Hispanic or Latino` = 41.7
    ),
    disease_state   = "Advanced solid tumors: non-small cell lung cancer (NSCLC) 34.4%, hormone-refractory prostate cancer (HRPC; modern term castration-resistant prostate cancer / CRPC) 35.4%, unresectable stage III or IV melanoma (MM = melanoma) 30.2%",
    dose_range      = "4.8 mg/m^2/day continuous IV infusion (CIVI) over 7 days (168 h) every 21 days; first cycle's dose computed using actual body surface area; concentrations and doses expressed as the cationic moiety of sepantronium bromide",
    regions         = "NSCLC cohort enrolled in Europe (LUCY trial); HRPC and melanoma cohorts enrolled in North America (PACY and MACY trials)",
    ecog_distribution = "ECOG 0 (asymptomatic) 39.6%, ECOG 1 (symptomatic) 55.2%, ECOG 2 (ambulatory <50%) 5.2%; melanoma trial accepted only ECOG 0-1, NSCLC and HRPC accepted ECOG 0-2",
    renal_function  = "Cockcroft-Gault creatinine clearance 31-180 mL/min (median 79); raw mL/min, NOT BSA-normalized",
    hepatic_function = "Baseline ALT 6-185 U/L (median 20); AST 12-92 U/L (median 25); no patients with severe hepatic impairment enrolled per Discussion",
    notes           = "Baseline demographics per Aoyama 2012 Table 2 (N = 96; 578 plasma sepantronium concentrations across cycles 1-6). Three open-label multicenter Phase 2 trials pooled: LUCY (NSCLC, 33 patients, Europe), PACY (HRPC, 34 patients, North America), MACY (melanoma, 29 patients, North America). All HRPC patients were male; all NSCLC patients were Caucasian. Median values of alpha-1-acid glycoprotein and AST differed by cancer type (AAG 22 umol/L in melanoma vs 35 umol/L in NSCLC/HRPC; AST 34 U/L in HRPC vs 21-24 U/L in NSCLC/melanoma); most other demographics did not differ significantly across cancer types (P > 0.05 for age, height, serum creatinine, and CLCR). Bioanalysis: LC-MS/MS at PPD Central Laboratory, LLOQ 0.05 ng/mL. PK samples drawn during the 7-day CIVI of cycle 1 plus 0.5-4 h and 6-24 h after stop of CIVI, and during the CIVI of subsequent cycles."
  )

  ini({
    # Structural parameters at the typical patient: CRCL 79 mL/min (median),
    # ALT 20 U/L (median), NSCLC cancer type (TUMTP_HRPC = 0 and TUMTP_MEL = 0).
    # Values from Aoyama 2012 Table 3 final-model 'Estimate' column.
    lcl <- log(42.1); label("Clearance (L/h)")                             # Aoyama 2012 Table 3 final model: CL = 42.1 L/h
    lvc <- log(319);  label("Central volume of distribution (L)")          # Aoyama 2012 Table 3 final model: V = 319 L

    # Continuous covariate effects on CL - power-form on (cov / ref).
    e_crcl_cl <- 0.425; label("Power exponent of CRCL on CL (unitless)")   # Aoyama 2012 Table 3 final model: CLCR (power) = 0.425
    e_alt_cl  <- 0.124; label("Power exponent of ALT on CL (unitless)")    # Aoyama 2012 Table 3 final model: ALT (power) = 0.124

    # Categorical covariate effects on CL - proportional-change form
    # (1 + e * I). The paper uses the power form THETA^I; for a binary
    # indicator this is mathematically equivalent to a proportional change
    # with coefficient (THETA - 1). NSCLC is the implicit reference when
    # both TUMTP_HRPC and TUMTP_MEL are 0.
    e_hrpc_cl <- -0.045; label("Proportional change of HRPC vs NSCLC on CL (unitless)") # Aoyama 2012 Table 3 final model: HRPC ratio = 0.955; proportional = (0.955 - 1) = -0.045
    e_mel_cl  <-  0.24;  label("Proportional change of melanoma vs NSCLC on CL (unitless)") # Aoyama 2012 Table 3 final model: MM (melanoma) ratio = 1.24; proportional = (1.24 - 1) = +0.24

    # Inter-individual variability on CL only. The paper reports the variance
    # of an exponential IIV (CL_i = TVCL * exp(eta_CL)) as omega^2 = 0.0385,
    # corresponding to 19.6% CV by the small-omega approximation sqrt(0.0385).
    # No IIV was kept on V because the available terminal-phase data (only
    # 0.5-24 h post-CIVI sampling) made V and CL inter-individual variability
    # difficult to separate (paper Discussion).
    etalcl ~ 0.0385                                                         # Aoyama 2012 Table 3 final model: omega^2_CL = 0.0385 (CV approx 19.6%)

    # Residual error (proportional). Aoyama 2012 fit two proportional residual
    # errors: a typical-patient error (sigma^2 = 0.0934; SD = sqrt(0.0934) =
    # 0.3056, i.e. 30.6% CV) and a much larger outlier-patient error
    # (sigma^2 = 31.6) applied to 11 of 96 patients (~11.5%) who had at least
    # one plasma concentration exceeding the post-hoc IQR-based threshold
    # Q3 + 3*IQR = 23.13 ng/mL during the CIVI. The outlier classification is
    # data-driven (a post-hoc identification of patients with anomalously high
    # observed concentrations) and is not a population-level predictor that a
    # simulation of new patients could mimic; only the typical-patient error
    # is encoded here. See the vignette 'Assumptions and deviations' section
    # for the discussion of the second residual error.
    propSd <- sqrt(0.0934); label("Proportional residual error (fraction)") # Aoyama 2012 Table 3 final model: sigma^2 w/o outliers = 0.0934; SD = sqrt(0.0934) = 0.3056 (30.6% CV)
  })

  model({
    # Individual PK parameters at the typical patient (NSCLC, CRCL 79 mL/min,
    # ALT 20 U/L). Continuous covariates apply as power factors on (cov / ref);
    # categorical covariates apply as proportional changes (1 + e * I).
    # NSCLC is the implicit reference category for cancer type.
    cl <- exp(lcl + etalcl) *
      (CRCL / 79)^e_crcl_cl *
      (ALT  / 20)^e_alt_cl *
      (1 + e_hrpc_cl * TUMTP_HRPC) *
      (1 + e_mel_cl  * TUMTP_MEL)

    vc <- exp(lvc)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> central / vc has units mg/L. Convert mg/L to
    # ng/mL by multiplying by 1000 (1 mg/L = 1000 ng/mL).
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
