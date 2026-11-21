Jin_2024_parkinson_progression <- function() {
  description <- paste0(
    "Long-term (about 40-year) disease-progression model of Parkinson's ",
    "disease fit with SReFT (Statistical Restoration of Fragmented Time ",
    "course) to four biomarkers of de novo and LRRK2-genetic-cohort PD ",
    "patients in the Parkinson's Progression Markers Initiative (PPMI): ",
    "MDS-UPDRS Part I-III total score, SCOPA-AUT score, DaTSCAN mean-striatum ",
    "specific binding ratio and CSF alpha-synuclein. No drug input. Each ",
    "biomarker follows a log-Gompertz curve on its log-standardized scale, ",
    "Y = alpha + beta / gamma * (exp(gamma * dT * t) - 1), in disease time t ",
    "(years since onset, defined as MDS-UPDRS total = 10). A single per-subject ",
    "progression factor dT = 1 + effects of female sex, LRRK2 G2019S ",
    "(rs34637584) carriage and LRRK2 rs76904798 carriage scales disease time ",
    "in every biomarker (dT < 1 is slower progression). Additive IIV on ",
    "alpha and gamma (SCOPA-AUT, DaTSCAN, alpha-synuclein) and on beta ",
    "(MDS-UPDRS); additive residual error on each log-standardized scale. ",
    "Model time is disease time, not study time: SReFT solves each subject's ",
    "onset-to-entry offset by an inverse-Gompertz average rather than as a ",
    "random effect, so a user with study-time data adds the subject's offset ",
    "to time. Raw-scale outputs (points, SBR ratio, pg/mL) are algebraic ",
    "back-transforms of the fitted log-standardized predictions."
  )
  reference <- paste(
    "Jin R, Yoshioka H, Sato H, Hisaka A.",
    "Data-driven disease progression model of Parkinson's disease and effect",
    "of sex and genetic variants.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(4):649.",
    "doi:10.1002/psp4.13112.",
    sep = " "
  )
  vignette <- "Jin_2024_parkinson_progression"

  # SReFT places the random effects additively on the Gompertz shape
  # parameters (eta ~ N(0, omega^2) added to alpha / beta / gamma; Methods
  # Equation 1 and the authors' pred_sreft_sim_final.f90, a = theta + eta).
  # The shape parameters are signed and not log-transformed, so the
  # eta + l<param> pairing does not apply.
  paper_specific_etas <- c(
    "etagomp_beta_updrs",
    "etagomp_alpha_scopaaut",
    "etagomp_alpha_datscan",
    "etagomp_alpha_asyn",
    "etagomp_gamma_scopaaut",
    "etagomp_gamma_datscan",
    "etagomp_gamma_asyn"
  )

  units <- list(
    time = "year (disease time since onset; onset = MDS-UPDRS Part I-III total of 10)",
    dosing = "(none; disease-progression model with no drug input)",
    concentration = paste(
      "updrs_z, scopaaut_z, datscan_z, asyn_z: log-standardized biomarker",
      "(unitless z-score of the log-transformed observation);",
      "updrs, scopaaut: score points; datscan: specific binding ratio",
      "(unitless); asyn: pg/mL"
    )
  )

  covariateData <- list(
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male).",
      notes = "Enters the disease-time progression factor additively: ftdis = 1 + e_sexf_ftdis * SEXF + ...; women progress 21.6% slower (Jin 2024 Results; Table 2 column 'Sex: Women').",
      source_name = "Sex"
    ),
    SNP_LRRK2_RS34637584 = list(
      description = "LRRK2 p.G2019S (dbSNP rs34637584) carrier indicator: 1 = at least one mutant allele, 0 = none.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no mutant allele).",
      notes = "The paper dichotomized the genotype because few participants were homozygous (Methods 'SNP selection'; Table S1: 62 heterozygous, 1 homozygous). Carriers progress 25.4% slower in disease time (Results 'dT,G2019S = -0.254').",
      source_name = "G2019S"
    ),
    SNP_LRRK2_RS76904798 = list(
      description = "LRRK2-region rs76904798 carrier indicator: 1 = at least one mutant allele, 0 = none.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no mutant allele).",
      notes = "Dichotomized as for G2019S (Table S1: 124 heterozygous, 5 homozygous). Carriers progress 10.4% slower; the bootstrap 95% CI crosses zero (Results).",
      source_name = "rs76904798"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 481L,
    n_studies = 1L,
    age_range = "mean (SD) 61.52 (9.67) years at study inclusion (Table 1); PPMI required age >= 30 years at diagnosis",
    age_median = NA_character_,
    weight_range = NA_character_,
    weight_median = NA_character_,
    sex_female_pct = 38.0,
    race_ethnicity = c(White = 94.6, Black = 1.46, Asian = 1.66, Other = 2.29),
    disease_state = "Parkinson's disease: PPMI de novo PD cohort plus the genetic cohort and genetic registry participants with PD (GenPD, RegPD). Hoehn and Yahr stage I in 40.1% and II in 57.8% at baseline; mean disease duration since diagnosis 0.90 years.",
    dose_range = "n/a (no drug input; PD medication could start at any time after enrollment and was not modelled)",
    regions = "PPMI international multicentre study (sites not itemized in the paper)",
    genotype = "LRRK2 G2019S carriers 63 (13.1%), rs76904798 carriers 129 (26.8%), both 11 (2.3%) (Results; Table 1)",
    n_observations = c(updrs = 5916L, scopaaut = 3521L, datscan = 1529L, asyn = 1685L),
    notes = "Observation periods up to 8 (MDS-UPDRS), 9 (SCOPA-AUT), 4 (DaTSCAN) and 3 (alpha-synuclein) years; SReFT aligns them on a common disease-time axis spanning about 30 years. MDS-UPDRS Part III scores in the ON-medication state and alpha-synuclein > 8000 pg/mL were excluded (Methods 'Biomarker selection')."
  )

  ini({
    # Typical values: authors' final-model simulation program
    # (Data S1 Simultation_Finalmodel.R, $THETA). They round to Jin 2024
    # Table 2 except gomp_beta_asyn (Table 2 -7.72e-5) and the third digits
    # of gomp_gamma_scopaaut / gomp_gamma_asyn and e_snp_lrrk2_rs34637584_ftdis;
    # see the vignette Assumptions and deviations section.
    gomp_alpha_updrs <- fixed(-2.4517); label("Gompertz alpha, MDS-UPDRS total: log-standardized level at disease time 0 (the value that puts the score at 10 at onset)") # Data S1 THETA(1) -2.4517 FIXED; Table 2 alpha MDS-UPDRS -2.45 (fixed)
    gomp_alpha_scopaaut <- -1.6802; label("Gompertz alpha, SCOPA-AUT: log-standardized level at disease time 0") # Data S1 THETA(2); Table 2 -1.68
    gomp_alpha_datscan <- 1.5773; label("Gompertz alpha, DaTSCAN SBR: log-standardized level at disease time 0") # Data S1 THETA(3); Table 2 1.58
    gomp_alpha_asyn <- 0.015194; label("Gompertz alpha, CSF alpha-synuclein: log-standardized level at disease time 0") # Data S1 THETA(4); Table 2 0.0153

    gomp_beta_updrs <- 0.21741; label("Gompertz beta, MDS-UPDRS total: initial slope (1/year)") # Data S1 THETA(5); Table 2 0.217
    gomp_beta_scopaaut <- 0.16796; label("Gompertz beta, SCOPA-AUT: initial slope (1/year)") # Data S1 THETA(6); Table 2 0.168
    gomp_beta_datscan <- -0.091738; label("Gompertz beta, DaTSCAN SBR: initial slope (1/year)") # Data S1 THETA(7); Table 2 -0.0917
    gomp_beta_asyn <- -4.0855e-05; label("Gompertz beta, CSF alpha-synuclein: initial slope (1/year)") # Data S1 THETA(8); Table 2 prints -7.72e-5

    gomp_gamma_updrs <- -0.0010217; label("Gompertz gamma, MDS-UPDRS total: curvature rate (1/year)") # Data S1 THETA(9); Table 2 -0.00102
    gomp_gamma_scopaaut <- -0.026625; label("Gompertz gamma, SCOPA-AUT: curvature rate (1/year)") # Data S1 THETA(10); Table 2 -0.0265
    gomp_gamma_datscan <- 0.069864; label("Gompertz gamma, DaTSCAN SBR: curvature rate (1/year)") # Data S1 THETA(11); Table 2 0.0698
    gomp_gamma_asyn <- -0.079492; label("Gompertz gamma, CSF alpha-synuclein: curvature rate (1/year)") # Data S1 THETA(12); Table 2 -0.0797

    e_sexf_ftdis <- -0.216; label("Additive effect of female sex on the disease-time progression factor (unitless)") # Data S1 THETA(13) -2.16E-01; Results dT,Sex = -0.216
    e_snp_lrrk2_rs34637584_ftdis <- -0.254; label("Additive effect of LRRK2 G2019S carriage on the disease-time progression factor (unitless)") # Data S1 THETA(14) -2.54E-01; Results dT,G2019S = -0.254 (Table 2 -0.255)
    e_snp_lrrk2_rs76904798_ftdis <- -0.104; label("Additive effect of LRRK2 rs76904798 carriage on the disease-time progression factor (unitless)") # Data S1 THETA(15) -1.04E-01; Results dT,rs76904798 = -0.104

    # Residual error: $SIGMA holds variances on the log-standardized scale
    # (Data S1 0.29943 0.2489 0.13736 0.21598; Table 2 'Residuals' 0.299,
    # 0.249, 0.137, 0.216). Converted to SD.
    addSd_updrs_z <- sqrt(0.29943); label("Additive residual SD, log-standardized MDS-UPDRS total") # Data S1 SIGMA(1) variance 0.29943
    addSd_scopaaut_z <- sqrt(0.2489); label("Additive residual SD, log-standardized SCOPA-AUT") # Data S1 SIGMA(2) variance 0.2489
    addSd_datscan_z <- sqrt(0.13736); label("Additive residual SD, log-standardized DaTSCAN SBR") # Data S1 SIGMA(3) variance 0.13736
    addSd_asyn_z <- sqrt(0.21598); label("Additive residual SD, log-standardized CSF alpha-synuclein") # Data S1 SIGMA(4) variance 0.21598

    # IIV (additive, variances). alpha_updrs, gamma_updrs and the beta of the
    # other three biomarkers carry no eta (0 FIXED in Data S1 $OMEGA).
    etagomp_beta_updrs ~ 1.87e-09 # Data S1 OMEGA(5) variance; Table 2 prints 2.15e-9
    etagomp_alpha_scopaaut ~ 0.49882 # Data S1 OMEGA(2); Table 2 0.498
    etagomp_alpha_datscan ~ 0.63096 # Data S1 OMEGA(3); Table 2 0.631
    etagomp_alpha_asyn ~ 0.76947 # Data S1 OMEGA(4); Table 2 0.769
    etagomp_gamma_scopaaut ~ 0.00021902 # Data S1 OMEGA(10); Table 2 2.19e-4
    etagomp_gamma_datscan ~ 0.0012104 # Data S1 OMEGA(11); Table 2 0.00121
    etagomp_gamma_asyn ~ 0.00056341 # Data S1 OMEGA(12); Table 2 5.63e-4
  })
  model({
    # Disease-time progression factor dT (Methods, Equation 1; Data S1
    # covt = 1 + THETA(13)*Sex + THETA(14)*G2019S + THETA(15)*rs76904798).
    ftdis <- 1 + e_sexf_ftdis * SEXF +
      e_snp_lrrk2_rs34637584_ftdis * SNP_LRRK2_RS34637584 +
      e_snp_lrrk2_rs76904798_ftdis * SNP_LRRK2_RS76904798
    # time is disease time (years since onset); dT rescales it.
    tdis_scaled <- ftdis * time

    # Individual Gompertz shape parameters (additive etas)
    alpha_updrs <- gomp_alpha_updrs
    beta_updrs <- gomp_beta_updrs + etagomp_beta_updrs
    gamma_updrs <- gomp_gamma_updrs

    alpha_scopaaut <- gomp_alpha_scopaaut + etagomp_alpha_scopaaut
    beta_scopaaut <- gomp_beta_scopaaut
    gamma_scopaaut <- gomp_gamma_scopaaut + etagomp_gamma_scopaaut

    alpha_datscan <- gomp_alpha_datscan + etagomp_alpha_datscan
    beta_datscan <- gomp_beta_datscan
    gamma_datscan <- gomp_gamma_datscan + etagomp_gamma_datscan

    alpha_asyn <- gomp_alpha_asyn + etagomp_alpha_asyn
    beta_asyn <- gomp_beta_asyn
    gamma_asyn <- gomp_gamma_asyn + etagomp_gamma_asyn

    # Log-Gompertz: Y = alpha + beta / gamma * (exp(gamma * dT * t) - 1),
    # with the linear limit alpha + beta * dT * t at gamma = 0 (the
    # c(bm) == 0 branch of the authors' pred_sreft_sim_final.f90).
    gfun_updrs <- tdis_scaled
    if (gamma_updrs != 0) gfun_updrs <- (exp(gamma_updrs * tdis_scaled) - 1) / gamma_updrs
    gfun_scopaaut <- tdis_scaled
    if (gamma_scopaaut != 0) gfun_scopaaut <- (exp(gamma_scopaaut * tdis_scaled) - 1) / gamma_scopaaut
    gfun_datscan <- tdis_scaled
    if (gamma_datscan != 0) gfun_datscan <- (exp(gamma_datscan * tdis_scaled) - 1) / gamma_datscan
    gfun_asyn <- tdis_scaled
    if (gamma_asyn != 0) gfun_asyn <- (exp(gamma_asyn * tdis_scaled) - 1) / gamma_asyn

    updrs_z <- alpha_updrs + beta_updrs * gfun_updrs
    scopaaut_z <- alpha_scopaaut + beta_scopaaut * gfun_scopaaut
    datscan_z <- alpha_datscan + beta_datscan * gfun_datscan
    asyn_z <- alpha_asyn + beta_asyn * gfun_asyn

    # Back-transform to the raw scale (Data S1 DVToDVraw()). The
    # standardization constants (mean and SD of log(y) in the analysis
    # data) are AVEBM = -1.6702, -1.6641, 0.13179, 7.1978 and
    # SDBM = 0.57152, 0.83675, 0.36341, 0.41332. The two questionnaire
    # scores were first mapped onto (0, 1) by Equations 3-4,
    # X = (x + 0.5) / (Xmax + 1), y = X / (1 - X), with Xmax = 236
    # (MDS-UPDRS I-III total) and 69 (SCOPA-AUT); DaTSCAN and
    # alpha-synuclein were log-transformed directly.
    y_updrs <- exp(updrs_z * 0.57152 - 1.6702)
    updrs <- y_updrs / (1 + y_updrs) * (236 + 1) - 0.5
    y_scopaaut <- exp(scopaaut_z * 0.83675 - 1.6641)
    scopaaut <- y_scopaaut / (1 + y_scopaaut) * (69 + 1) - 0.5
    datscan <- exp(datscan_z * 0.36341 + 0.13179)
    asyn <- exp(asyn_z * 0.41332 + 7.1978)

    updrs_z ~ add(addSd_updrs_z)
    scopaaut_z ~ add(addSd_scopaaut_z)
    datscan_z ~ add(addSd_datscan_z)
    asyn_z ~ add(addSd_asyn_z)
  })
}
