Aoyama_2017_meloxicam <- function() {
  description <- "Two-compartment population PK model for oral meloxicam in healthy Japanese, Chinese, Korean and white adult men, with parallel zero-order (into central) and lagged first-order absorption, CYP2C9 *2 / *3 allele-count effects on CL/F and a lean-body-mass power effect on Vc/F, coupled to the indirect-response serum thromboxane B2 (TXB2) model with Hill-type inhibition of TXB2 production used for the paper's steady-state covariate simulations (Aoyama 2017)"
  reference <- paste(
    "Aoyama T, Hasunuma T, Yamaguchi T, Matsumoto Y. Pharmacokinetics and",
    "Pharmacodynamics of Meloxicam in East Asian Populations: The Role of",
    "Ethnicity on Drug Response. CPT Pharmacometrics Syst Pharmacol.",
    "2017;6(12):823-832. doi:10.1002/psp4.12259.",
    sep = " "
  )
  vignette <- "Aoyama_2017_meloxicam"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")
  # Unit note: doses are in mg and volumes in L, so central / vc is in mg/L
  # (= ug/mL); the factor 1000 converts to ng/mL, the unit of every reported
  # concentration and of the IC50 (Table 3). TXB2 is a percentage of the
  # subject's own basal serum TXB2 generation (Eq. 7), not a concentration.

  covariateData <- list(
    CYP2C9_S2_COUNT = list(
      description = "Number of CYP2C9*2 alleles (0, 1 or 2)",
      units = "(count, 0/1/2 alleles per subject)",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL/F additively per allele, CL/F = 0.391 * (1 + n*2 * (-0.147) + n*3 * (-0.400)) (Eq. 5; Table 2 footnote). The control stream (Supplementary Files) codes the two alleles as CYPL in {1, 2} and CYPR in {1, 2, 3}, so the allele counts are the number of 2's (for this column) and 3's (for CYP2C9_S3_COUNT) across the two. Cohort: 5 *1/*2 and 1 *2/*2 (white) plus 1 *1/*2 (Chinese) of 119 subjects.",
      source_name = "CYPL / CYPR (n*2)"
    ),
    CYP2C9_S3_COUNT = list(
      description = "Number of CYP2C9*3 alleles (0, 1 or 2)",
      units = "(count, 0/1/2 alleles per subject)",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL/F additively per allele (Eq. 5). Cohort: 15 *1/*3 subjects (4 Japanese, 3 Chinese, 4 Korean, 4 white); no *2/*3 or *3/*3 subjects were observed, so the *2/*3 and *3/*3 CL/F values (0.45 and 0.20 x the *1/*1 value) are extrapolations of the per-allele assumption, as the paper notes in its Discussion.",
      source_name = "CYPL / CYPR (n*3)"
    ),
    LBM = list(
      description = "Lean body mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "James formula (Methods, ref. 24): for men LBM = 1.1 * WT - 128 * (WT / height)^2 with WT in kg and height in cm; all subjects were men. Power effect on Vc/F normalised to the cohort median 55.0 kg (Table 2 footnote: Vcpop = 7.79 * (LBM / 55.0)^1.05; control stream COV1 = LBM/55). Cohort range 43.6-75.2 kg (Table 1, Methods 'Simulation').",
      source_name = "LBM"
    )
  )

  covariatesDataExcluded <- list(
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Entered the full model on Vc/F in the forward step (Results; Supplementary Table S2) but was removed in backward elimination; not in the final model."
    ),
    RACE_CHINESE = list(
      description = "Chinese ethnicity indicator (1 = Chinese, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Japanese, the reference of Eq. 6)",
      notes = "Ethnic differences (n_ch, n_ko, n_ca vs Japanese; Eq. 6) were tested last in the forward step on CL/F, Vc/F, F and Ka (Supplementary Table S1); none was statistically significant (Results). Not retained."
    ),
    RACE_KOREAN = list(
      description = "Korean ethnicity indicator (1 = Korean, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Japanese, the reference of Eq. 6)",
      notes = "Tested as n_ko in Eq. 6; not statistically significant, not retained."
    ),
    RACE_WHITE = list(
      description = "White ethnicity indicator (1 = white, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Japanese, the reference of Eq. 6)",
      notes = "Tested as n_ca in Eq. 6; not statistically significant, not retained."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "meloxicam", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "meloxicam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "meloxicam", units = "mg", specimen = "plasma", verified = TRUE),
    txb2 = list(analyte = "thromboxane B2", units = "% of basal value", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 119,
    n_studies = 1,
    age_range = "21-35 years (group medians 24-31)",
    weight_range = "51.0-100 kg (group medians 63.5-74.5)",
    sex_female_pct = 0,
    race_ethnicity = c(Japanese = 30, Chinese = 30, Korean = 29, White = 30),
    disease_state = "Healthy volunteers",
    dose_range = "Single oral 7.5 mg meloxicam (PK model); PD layer from single oral 15 mg in 21 healthy Korean men (Bae 2007)",
    regions = "Japan (single-site bioanalysis of a Japanese / Chinese / Korean / white trial, UMIN000004173)",
    cyp2c9_genotype = c(`*1/*1` = 97, `*1/*2` = 6, `*2/*2` = 1, `*1/*3` = 15),
    lbm_range = "43.6-75.2 kg (median 55.0)",
    notes = "PK: Hasunuma et al. trial with a single drug lot, matched meals and drinking water; plasma sampled pre-dose and 1, 2, 3, 4, 5, 6, 8, 12, 24, 36, 48, 60 and 72 h after a single 7.5 mg oral dose (Methods; Table 1). PD (TXB2) parameters are from a separate sequential PK/PD fit to digitised mean data from Bae et al. 2007 (12 CYP2C9*1/*1 and 9 *1/*13 Korean men, 15 mg single dose; Table 3), which the paper combines with this population PK model for the steady-state simulations of Figure 4 and Results 'Simulation'."
  )

  ini({
    # Structural PK (Table 2 'Original dataset' estimates; NONMEM 7.2 FOCE-I).
    lcl <- log(0.391); label("Apparent clearance CL/F for CYP2C9*1/*1 (L/h)") # Table 2 'CL (L/h)' 0.391 (0.375 to 0.407)
    lvc <- log(7.79); label("Apparent central volume Vc/F at LBM 55.0 kg (L)") # Table 2 'Vc (L)' 7.79 (7.24 to 8.34)
    lq <- log(1.24); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2 'Q (L/h)' 1.24 (0.948 to 1.53)
    lvp <- log(2.73); label("Apparent peripheral volume Vp/F (L)") # Table 2 'Vp (L)' 2.73 (2.20 to 3.26)
    lka <- log(2.00); label("First-order absorption rate constant Ka (1/h)") # Table 2 'Ka (/h)' 2.00 (1.38 to 2.62)
    ld1 <- log(1.91); label("Duration of zero-order input into central DT, also the first-order lag time (h)") # Table 2 'DT (h)' 1.91 (1.86 to 1.96); footnote 'tlag = DT'
    # Table 2 reports F, the fraction of the dose absorbed by the ZERO-order
    # process. The canonical logitffo is the complementary first-order
    # fraction, so logitffo = logit(1 - F) = -logit(F).
    logitffo <- qlogis(1 - 0.425); label("Logit of the fraction of the dose absorbed by the lagged first-order route, 1 - F (unitless)") # Table 2 'F' 0.425 (0.367 to 0.483)

    # Covariate effects (Eq. 5; Table 2 footnote equations).
    e_cyp2c9s2_cl <- -0.147; label("Fractional change in CL/F per CYP2C9*2 allele (unitless)") # Table 2 'CYP2C9 *2 on CL' -0.147 (-0.234 to -0.0604)
    e_cyp2c9s3_cl <- -0.400; label("Fractional change in CL/F per CYP2C9*3 allele (unitless)") # Table 2 'CYP2C9 *3 on CL' -0.400 (-0.488 to -0.312)
    e_lbm_vc <- 1.05; label("Power exponent of (LBM/55.0) on Vc/F (unitless)") # Table 2 'LBM on Vd' 1.05 (0.695 to 1.40)

    # PD: indirect response with inhibition of TXB2 production (Eq. 7),
    # from the sequential PK/PD fit (Table 3, shared by both genotype groups).
    lkout <- log(0.912); label("First-order loss rate constant of percent serum TXB2 generation kout (1/h)") # Table 3 'kout (/h)' 0.912 (RSE 17.2%)
    lic50 <- log(1390); label("Plasma meloxicam concentration halving TXB2 production IC50 (ng/mL)") # Table 3 'IC50 (ng/mL)' 1390 (RSE 10.2%)
    lhill <- log(0.739); label("Sigmoidicity of the inhibition gamma (unitless)") # Table 3 'gamma' 0.739 (RSE 9.62%)

    # Between-subject variability. Table 2 prints omegaCL, omegaVd and omegaKa
    # as CV% = 100 * sqrt(exp(omega^2) - 1): the printed +/-1.96 SE interval
    # of omegaKa (68.1 to 201%) is symmetric about the estimate only on that
    # scale (variance interval 0.383 to 1.621, midpoint 1.002 vs 0.999),
    # whereas the 100 * sqrt(omega^2) reading gives an asymmetric one
    # (midpoint 2.25 vs 1.72). Variances are therefore log(1 + CV^2).
    etalcl ~ 0.044371 # Table 2 'omegaCL (CV%)' 21.3 -> log(1 + 0.213^2)
    etalvc ~ 0.029153 # Table 2 'omegaVd (CV%)' 17.2 -> log(1 + 0.172^2)
    # Ka and logit(F) etas are correlated (Table 2 'omegaKa-omegaF' correlation
    # 0.243). omega^2 Ka = log(1 + 1.31^2) = 0.999203; omega^2 logit(F) = 2.03
    # as printed. Because logitffo = -logit(F), the eta on logitffo is the
    # negated eta on logit(F): same variance, covariance sign flipped:
    # -0.243 * sqrt(0.999203 * 2.03) = -0.346081.
    etalka + etalogitffo ~ c(0.999203, -0.346081, 2.03) # Table 2 'omegaKa (CV%)' 131, 'omega2 F' 2.03, 'omegaKa-omegaF' 0.243

    propSd <- 0.123; label("Proportional residual error (fraction)") # Table 2 'sigma (CV%)' 12.3 (11.3 to 13.2); $ERROR Y = IPRED*(1+EPS(1))
  })

  model({
    # Individual PK parameters (Eq. 5 and Table 2 footnote; control stream $PK).
    cl <- exp(lcl + etalcl) * (1 + e_cyp2c9s2_cl * CYP2C9_S2_COUNT + e_cyp2c9s3_cl * CYP2C9_S3_COUNT)
    vc <- exp(lvc + etalvc) * (LBM / 55.0)^e_lbm_vc
    q <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka + etalka)
    d1 <- exp(ld1)
    ffo <- expit(logitffo + etalogitffo)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # PD parameters (Eq. 7): basal percent TXB2 generation is 100%, so
    # kin = 100 * kout (Methods; Table 3 footnote).
    kout <- exp(lkout)
    kin <- 100 * kout
    ic50 <- exp(lic50)
    hill <- exp(lhill)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Parallel absorption (Figure 2; control stream ADVAN6 with F1 = 1 - F2,
    # ALAG1 = D2): each oral dose is TWO records with the same amt -- one on
    # cmt = "central" with rate = -2 (the zero-order fraction F = 1 - ffo,
    # infused over d1) and one on cmt = "depot" (the first-order fraction ffo,
    # lagged by tlag = DT so it starts as the zero-order input ends).
    f(central) <- 1 - ffo
    dur(central) <- d1
    f(depot) <- ffo
    alag(depot) <- d1

    Cc <- 1000 * central / vc

    inh <- Cc^hill / (ic50^hill + Cc^hill)
    d/dt(txb2) <- kin * (1 - inh) - kout * txb2
    txb2(0) <- 100
    TXB2 <- txb2

    Cc ~ prop(propSd)
  })
}
