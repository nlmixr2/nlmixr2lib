Kloprogge_2020_hrze_bacterialload <- function() {
  description <- "Sputum bacillary-load PKPD model for Malawian adults with drug-sensitive pulmonary tuberculosis on standard isoniazid-rifampicin-pyrazinamide-ethambutol therapy (Kloprogge 2020). A single ODE gives first-order decline of sputum CFU from an estimated baseline, with a kill rate LAM that falls by a fraction BETA with half-time T1/2 (biphasic decline). LAM increases with steady-state isoniazid AUC0-24 and baseline bilirubin and is lower with alcohol consumption; T1/2 lengthens with steady-state rifampicin AUC0-24. Exposures enter as per-subject covariates (umol*h/L), not as dynamic concentrations. Full 4 x 4 IIV block; additive residual error on log10 CFU/mL."
  reference <- paste(
    "Kloprogge F, Mwandumba HC, Banda G, Kamdolozi M, Shani D, Corbett EL,",
    "Kontogianni N, Ward S, Khoo SH, Davies GR, Sloan DJ. (2020).",
    "Longitudinal pharmacokinetic-pharmacodynamic biomarkers correlate with",
    "treatment outcome in drug-sensitive pulmonary tuberculosis: a population",
    "pharmacokinetic-pharmacodynamic analysis.",
    "Open Forum Infect Dis 7(7):ofaa218. doi:10.1093/ofid/ofaa218.",
    sep = " "
  )
  vignette <- "Kloprogge_2020_tuberculosis"
  units <- list(time = "h", dosing = "mg", concentration = "CFU/mL (cfu state); log10 CFU/mL (log_cfu output)")

  covariateData <- list(
    AUC_INH = list(
      description = "Steady-state isoniazid AUC0-24 (individual empirical Bayes estimate from the isoniazid popPK model).",
      units = "umol*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Exponential effect on LAM, exp(theta * (AUC_INH - median)) (Supplementary Materials section 4). Molar units per the Figure 4 caption ('54.9-515 hr x umol/L' for isoniazid), consistent with the umol/L PK data (Figure 1 y-axis). The median of the 102-patient PKPD cohort is not printed; the model centres on the 154-patient PK-cohort median AUC0-24 of 18.83 mg*h/L (Results 'Pharmacokinetics') = 18.83 / 137.14 g/mol x 1000 = 137.3 umol*h/L. Compute from Kloprogge_2020_isoniazid as dose / CL at steady state, converted with MW 137.14 g/mol.",
      source_name = "AUCinh"
    ),
    AUC_RIF = list(
      description = "Steady-state rifampicin AUC0-24 (individual empirical Bayes estimate from the rifampicin popPK model).",
      units = "umol*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Exponential effect on T1/2, exp(theta * (AUC_RIF - median)) (Supplementary Materials section 4). Molar units per the Figure 4 caption ('19.9-145 hr x umol/L' for rifampicin). Centred on the PK-cohort median AUC0-24 of 29.10 mg*h/L (Results 'Pharmacokinetics') = 29.10 / 822.94 g/mol x 1000 = 35.36 umol*h/L. Compute from Kloprogge_2020_rifampicin as dose / CL at steady state, converted with MW 822.94 g/mol.",
      source_name = "AUCrif"
    ),
    TBILI = list(
      description = "Baseline total serum bilirubin.",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline value, time-fixed per subject. Exponential effect on LAM, exp(theta * (TBILI - median)), centred on the PKPD-cohort median 7.5 umol/L (Table 1 'Pharmacodynamic Data' column, range 1-32 umol/L).",
      source_name = "BLbilirubin"
    ),
    ALCOHOL_USE = list(
      description = "Current alcohol consumption (any beer or spirits) at study recruitment, 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no alcohol consumption)",
      notes = "Self-reported current practice at recruitment (Methods 'Study Population and Study Design': 'Alcohol consumption (any beer or spirits) and smoking were reported as binary covariates based on practice at the time of study recruitment'); 32% of the PKPD cohort (Table 1). Proportional effect on LAM, (1 + theta * ALCOHOL_USE) (Supplementary Materials section 4 'Categorical variables were evaluated using a proportional ((1 + theta) equation').",
      source_name = "alcohol"
    )
  )

  covariatesDataExcluded <- list(
    AUC_INH_MIC = list(
      description = "Isoniazid AUC0-24 / MIC.",
      units = "h",
      type = "continuous",
      notes = "Tested in place of AUC_INH on LAM in the 50 patients with MIC data (Supplementary S3 Table); not part of the final full-cohort model and its coefficient is not reported."
    ),
    AUC_RIF_MIC = list(
      description = "Rifampicin AUC0-24 / MIC.",
      units = "h",
      type = "continuous",
      notes = "Tested in place of AUC_RIF on T1/2 in the MIC subset; the model did not converge (Supplementary S3 Table)."
    )
  )

  compartmentData <- list(
    cfu = list(analyte = "Mycobacterium tuberculosis", units = "CFU/mL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 102L,
    n_studies = 1L,
    age_range = "17-60 years",
    age_median = "30 years",
    weight_range = "35-74 kg",
    weight_median = "53 kg",
    sex_female_pct = 26,
    disease_state = "Smear-positive, drug-sensitive pulmonary tuberculosis on standard first-line therapy (RZHE intensive phase); 58% HIV co-infected; 32% reported alcohol consumption.",
    dose_range = "Daily fixed-dose-combination RZHE tablets (rifampicin 150 mg, isoniazid 75 mg, pyrazinamide 400 mg, ethambutol 275 mg) by weight band (2-5 tablets).",
    regions = "Malawi (Queen Elizabeth Central Hospital, Blantyre)",
    notes = "Serial sputum colony counting (SSCC, log10 CFU/mL on solid culture) on days 0, 2, 4, 7, 14, 28, 49 and 56; patients with at least 2 bacterial-load measurements were included. Counts below the limit of quantification were handled with the M3 method. Demographics from Table 1 'Pharmacodynamic Data' column."
  )

  ini({
    # Supplementary S2 Table 'Summary of population pharmacokinetic-
    # pharmacodynamic parameter estimates' (Estimate column).
    e0_cfu <- 6.3867; label("Baseline sputum bacillary load (log10 CFU/mL)") # S2 Table 'BL (log10 CFU/ml)' = 6.3867; normally distributed (Supplementary section 4)
    llam <- log(0.0389); label("Initial bacterial kill rate LAM (1/h)")      # S2 Table 'LAM (CFU/ml per hour)' = 0.0389
    lt12 <- log(149.58); label("Half-time of the fall in kill rate T1/2 (h)") # S2 Table 'T1/2 (hours)' = 149.5800
    logitbeta <- logit(0.6445); label("Maximum fractional fall in kill rate BETA (fraction, logit scale)") # S2 Table 'BETA [0,1]' = 0.6445

    # Covariate effects (Supplementary section 4 forms: exponential for
    # continuous covariates centred on the median, (1 + theta) for binary).
    e_auc_inh_lam <- 0.0005; label("Exponential effect of isoniazid AUC0-24 on LAM (per umol*h/L)")    # S2 Table 'LAM_AUCinh' = 0.0005
    e_tbili_lam <- 0.0029; label("Exponential effect of baseline bilirubin on LAM (per umol/L)")         # S2 Table 'LAM_BLbilirubin' = 0.0029
    e_alcohol_lam <- -0.0377; label("Proportional effect of alcohol consumption on LAM (fraction)")      # S2 Table 'LAM_alcohol' = -0.0377
    e_auc_rif_t12 <- 0.0335; label("Exponential effect of rifampicin AUC0-24 on T1/2 (per umol*h/L)")  # S2 Table 'T12_AUCrif' = 0.0335

    # IIV: full block in S2 Table order BL, LAM, T1/2, BETA ('Between patient
    # variability was estimated in a block with OMEGA's being the
    # off-diagnals and parameter names with _var being the diagnoal
    # estimates'). BL IIV is additive on log10 CFU/mL; LAM and T1/2 are
    # log-normal; the BETA IIV sits on the logit scale (see model()).
    etae0_cfu + etallam + etalt12 + etalogitbeta ~ c(
      0.1135,
      0.0028, 0.0383,
      0, -0.0037, 0.0526,
      0.3846, 0.0102, -0.0085, 1.3149
    ) # S2 Table BL_var 0.1135; OMEGA.2.1. 0.0028; LAM_var 0.0383; OMEGA.3.1. 0.0000; OMEGA.3.2. -0.0037; T1/2_var 0.0526; OMEGA.4.1. 0.3846; OMEGA.4.2. 0.0102; OMEGA.4.3. -0.0085; BETA_var 1.3149

    # S2 Table 'RUV: additive residual variability on log 10 transformed
    # data' = 2.5832, a NONMEM variance: SD = sqrt(2.5832) = 1.607. The
    # Figure 3 VPC baseline 2.5th-97.5th percentile spread (about 3.2 to 9.6
    # log10 CFU/mL) reproduces with the variance and not with SD = 2.58.
    addSd <- 1.607; label("Additive residual error (log10 CFU/mL)") # S2 Table 'RUV' = 2.5832 (variance); sqrt(2.5832) = 1.607
  })
  model({
    # Individual parameters. Covariate centring values: isoniazid and
    # rifampicin PK-cohort median AUC0-24 (18.83 and 29.10 mg*h/L, Results
    # 'Pharmacokinetics') in umol*h/L (MW 137.14 and 822.94 g/mol), and the
    # PKPD-cohort median bilirubin 7.5 umol/L (Table 1).
    e0 <- e0_cfu + etae0_cfu
    lam <- exp(llam + etallam) *
      exp(e_auc_inh_lam * (AUC_INH - 137.3)) *
      exp(e_tbili_lam * (TBILI - 7.5)) *
      (1 + e_alcohol_lam * ALCOHOL_USE)
    t12 <- exp(lt12 + etalt12) * exp(e_auc_rif_t12 * (AUC_RIF - 35.36))
    # Supplementary section 4 states that all PKPD parameters other than
    # the baseline were log-normal, but S2 Table bounds BETA to [0, 1] and a
    # log-normal BETA with variance 1.31 exceeds 1 in 35% of patients,
    # giving regrowth that the Figure 3 VPC (95% of samples below the
    # limit of quantification by 1100 h) rules out. The logit-normal form
    # reproduces Figure 3.
    beta <- expit(logitbeta + etalogitbeta)

    # Supplementary section 4: dB/dt = -knet * B with
    # knet = LAM * (1 - BETA * (1 - exp(-time * ln(2) / T1/2))). The
    # printed knet carries a leading minus sign that, combined with the
    # minus in dB/dt, would make the load grow; the decline described
    # throughout the paper is used. `time` is time since the start of
    # tuberculosis treatment (h).
    knet <- lam * (1 - beta * (1 - exp(-time * log(2) / t12)))
    cfu(0) <- 10^e0
    d/dt(cfu) <- -knet * cfu

    log_cfu <- log10(cfu)
    log_cfu ~ add(addSd)
  })
}
