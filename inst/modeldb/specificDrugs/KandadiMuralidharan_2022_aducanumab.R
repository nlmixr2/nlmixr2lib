KandadiMuralidharan_2022_aducanumab <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order elimination for",
    "intravenous aducanumab in patients with early Alzheimer's disease,",
    "sequentially linked to an indirect-response model in which serum",
    "aducanumab stimulates (Emax) the first-order loss of the composite",
    "florbetapir amyloid PET standard uptake value ratio (SUVR)"
  )
  reference <- paste(
    "Kandadi Muralidharan K, Tong X, Kowalski KG, Rajagovindan R, Lin L,",
    "Budd Haberlain S, Nestorov I. Population pharmacokinetics and standard",
    "uptake value ratio of aducanumab, an amyloid plaque-removing agent, in",
    "patients with Alzheimer's disease. CPT Pharmacometrics Syst Pharmacol.",
    "2022;11(1):7-19. doi:10.1002/psp4.12728"
  )
  vignette <- "KandadiMuralidharan_2022_aducanumab"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Amyloid PET SUVR turnover state of the PK-PD layer; a paper-specific PD
  # readout with no registered compartment canonical.
  paper_specific_compartments <- c("suvr")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline (time-fixed). Power effects on CL, V1 and V2 normalised to",
        "70.9 kg (Model S1 control stream, CLWTBL / V1WTBL / V2WTBL) and on",
        "Kout normalised to 72 kg (Model S2 control stream, KOUTWTBL). The",
        "Results text's 'weighing 71.9 kg' reference patient is a typo for the",
        "70.9 kg used in the control stream and in the Figure 2 / Figure 3",
        "captions."
      ),
      source_name = "WTBL"
    ),
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline (time-fixed). Power effects normalised to 71 years on V2",
        "(Model S1, V2AGE), on baseline SUVR (Model S2, BSLAGE) and on Emax",
        "(Model S2, EMAXAGE)."
      ),
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 = female (the paper's reference patient is female)",
      notes = paste(
        "The source column SEXN is 1 = male, 0 = female; SEXF = 1 - SEXN.",
        "The paper's male effects enter as (1 + theta) multipliers when male,",
        "encoded here as (1 + e_male_<param> * (1 - SEXF)); coefficient signs",
        "are unchanged."
      ),
      source_name = "SEXN"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = White or Other)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = White + Other (pooled because of the small 'other race' category)",
      notes = paste(
        "Model S1 recodes RACEN = 1 if RACE == 2 (Asian), else 0; the paper",
        "tested race on PK as White + Other versus Asian (Results, PopPK",
        "model)."
      ),
      source_name = "RACEN"
    ),
    SCORE_MMSE = list(
      description = "Baseline Mini-Mental State Examination score",
      units = "(SCORE_MMSE units, 0-30 score)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline (time-fixed). Power effects normalised to 26 on V2",
        "(Model S1, V2MMSEBL) and on baseline SUVR (Model S2, BSLMMSEBL)."
      ),
      source_name = "MMSEBL"
    ),
    APOE4_HET = list(
      description = "APOE-epsilon4 heterozygote indicator (1 = one epsilon4 allele)",
      units = "(binary)",
      type = "binary",
      reference_category = paste(
        "Model reference genotype is the HETEROZYGOTE (APOENCAT = 1, 'Most",
        "common' in Model S2): APOE4_HET = 1, APOE4_HOM = 0 gives a baseline",
        "SUVR multiplier of 1."
      ),
      notes = paste(
        "Derived from the three-level APOENCAT column (0 = non-carrier,",
        "1 = one copy, 2 = two copies): APOE4_HET = as.integer(APOENCAT == 1).",
        "A non-carrier has APOE4_HET = APOE4_HOM = 0 and receives the",
        "e_apoe4non_rbase multiplier; the indicators are mutually exclusive."
      ),
      source_name = "APOENCAT"
    ),
    APOE4_HOM = list(
      description = "APOE-epsilon4 homozygote indicator (1 = two epsilon4 alleles)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = not homozygous (see APOE4_HET for the heterozygote model reference)",
      notes = "Derived from APOENCAT: APOE4_HOM = as.integer(APOENCAT == 2).",
      source_name = "APOENCAT"
    )
  )

  compartmentData <- list(
    central = list(analyte = "aducanumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "aducanumab", units = "mg", specimen = "tissue", verified = TRUE),
    suvr = list(
      analyte = "amyloid-beta plaque (composite cortical florbetapir PET SUVR, whole-cerebellum reference)",
      units = "(ratio)",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2961,
    n_studies = 5,
    age_range = "50-91 years",
    age_mean = "70.3 years (PopPK data set); median 71.3 years",
    weight_range = "35.6-162 kg",
    weight_mean = "71.6 kg (PopPK data set)",
    sex_female_pct = 52,
    race_ethnicity = c(White = 79, NonWhite = 21),
    disease_state = paste(
      "Early Alzheimer's disease (80% mild cognitive impairment due to AD) in",
      "the phase III ENGAGE (221AD301) and EMERGE (221AD302) trials; prodromal",
      "or mild AD in PRIME (221AD103); mild to moderate AD in 221AD101 and in",
      "the Japanese PROPEL (221AD104) study"
    ),
    dose_range = paste(
      "IV infusion: single doses 0.3-60 mg/kg; 1, 3, 6 or 10 mg/kg Q4W; and",
      "titration regimens to maintenance 3, 6 or 10 mg/kg Q4W"
    ),
    regions = "Multinational (phase III); Japan (PROPEL)",
    apoe4_carrier_pct = 69,
    baseline_mmse = "mean 26.0 (range 14-30)",
    n_observations = "50,306 serum concentrations (PK); 3655 SUVR measurements",
    pd_population = paste(
      "PopPK-PD data set: 1125 patients from PRIME, ENGAGE and EMERGE; mean",
      "age 70.4 (50-91) years, mean weight 73.2 (40.8-143) kg, 52% female,",
      "88% Caucasian, 70% ApoE epsilon4 carriers (52% one copy) (Table S2)"
    ),
    notes = paste(
      "Baseline demographics from Supplementary Table S2 and the Results",
      "'Analysis population' section. Anti-drug antibody incidence was 1% and",
      "was not tested as a covariate."
    )
  )

  ini({
    # --- PopPK model (Table 1; functional forms from Model S1) ---
    # Reference patient: female, White (non-Asian), 70.9 kg; V2 additionally
    # referenced to age 71 years and baseline MMSE 26.
    lcl <- log(0.0159); label("Clearance CL for the reference patient (L/h)") # Table 1 CL = 0.0159 L/h
    lvc <- log(3.59); label("Central volume of distribution V1 for the reference patient (L)") # Table 1 V1 = 3.59 L
    lvp <- log(6.04); label("Peripheral volume of distribution V2 for the reference patient (L)") # Table 1 V2 = 6.04 L
    lq <- log(0.0194); label("Intercompartmental clearance Q (L/h)") # Table 1 Q = 0.0194 L/h

    e_wt_cl <- 0.561; label("Power exponent of (WT/70.9) on CL (unitless)") # Table 1 WT on CL = 0.561
    e_wt_vc <- 0.506; label("Power exponent of (WT/70.9) on V1 (unitless)") # Table 1 WT on V1 = 0.506
    e_wt_vp <- 0.320; label("Power exponent of (WT/70.9) on V2 (unitless)") # Table 1 WT on V2 = 0.320
    e_age_vp <- 0.207; label("Power exponent of (AGE/71) on V2 (unitless)") # Table 1 Age on V2 = 0.207
    e_mmse_vp <- 0.182; label("Power exponent of (SCORE_MMSE/26) on V2 (unitless)") # Table 1 MMSEBL on V2 = 0.182
    e_race_asian_cl <- 0.125; label("Fractional change in CL for Asian vs White + Other (unitless)") # Table 1 race Asian CL = 0.125; form (1 + theta), Model S1 CLRACEN
    e_race_asian_vc <- -0.044; label("Fractional change in V1 for Asian vs White + Other (unitless)") # Table 1 race Asian V1 = -0.044
    e_race_asian_vp <- -0.148; label("Fractional change in V2 for Asian vs White + Other (unitless)") # Table 1 race Asian V2 = -0.148
    e_male_cl <- 0.134; label("Fractional change in CL for male vs female (unitless)") # Table 1 sex male CL = 0.134; form (1 + theta), Model S1 CLSEXN
    e_male_vc <- 0.146; label("Fractional change in V1 for male vs female (unitless)") # Table 1 sex male V1 = 0.146
    e_male_vp <- 0.129; label("Fractional change in V2 for male vs female (unitless)") # Table 1 sex male V2 = 0.129

    # Table 1 %CV = 100*sqrt(exp(omega^2) - 1) (footnote b), so
    # omega^2 = log(1 + CV^2): CL 0.216 -> 0.045600, V1 0.148 -> 0.021668,
    # V2 0.170 -> 0.028490. Covariances = rho * sqrt(omega_i^2 * omega_j^2):
    # rho(CL,V1) 0.378 -> 0.011882, rho(CL,V2) -0.407 -> -0.014670,
    # rho(V1,V2) 0.392 -> 0.009740. (Model S1 $OMEGA BLOCK(3) carries
    # 0.0457 / 0.0120 / 0.0217 / -0.0147 / 0.00973 / 0.0285, which reproduce
    # the Table 1 %CV and correlations to the printed precision.)
    etalcl + etalvc + etalvp ~ c(
      0.045600,
      0.011882, 0.021668,
      -0.014670, 0.009740, 0.028490
    )
    # Table 1 'Weighting on residual error' 34.9% CV -> log(1 + 0.349^2) = 0.114935;
    # Model S1 W2 = TVW2 * EXP(ETA(4)), $OMEGA 0.115
    etaruv ~ 0.114935

    propSd <- 0.148; label("Proportional residual error on serum aducanumab (fraction)") # Table 1 proportional error = 14.8%
    addSd <- 0.202; label("Additive residual error on serum aducanumab (mg/L)") # Table 1 additive error = 0.202 mg/L

    # --- PopPK-PD (SUVR) model (Table 2; final estimates and functional forms
    # from Model S2, whose $THETA / $OMEGA carry the final estimates at full
    # precision). Reference patient: ApoE4 heterozygous, age 71 years,
    # baseline MMSE 26, weight 72 kg (Kout).
    lrbase <- log(1.40005); label("Baseline composite SUVR for the reference patient (ratio)") # Table 2 BL = 1.40; Model S2 THETA(1) = 1.40005
    lkout <- log(8.51197e-05); label("First-order SUVR elimination rate constant Kout (1/h)") # Table 2 Kout = 8.51E-05 1/h; Model S2 THETA(2)
    lemax <- log(0.702378); label("Maximum fold change in SUVR elimination Emax (unitless)") # Table 2 Emax = 0.702; Model S2 THETA(3)
    lec50 <- log(46.4175); label("Serum aducanumab concentration giving half-maximal stimulation EC50 (mg/L)") # Table 2 EC50 = 46.4 mg/L; Model S2 THETA(4)

    e_wt_kout <- 0.413748; label("Power exponent of (WT/72) on Kout (unitless)") # Table 2 WT on Kout = 0.414; Model S2 THETA(12)
    e_age_rbase <- 0.100713; label("Power exponent of (AGE/71) on baseline SUVR (unitless)") # Table 2 age on BL = 0.101; Model S2 THETA(7)
    e_age_emax <- 1.94016; label("Power exponent of (AGE/71) on Emax (unitless)") # Model S2 THETA(11) EMAXAGE = 1.94016 (row absent from Table 2)
    e_mmse_rbase <- -0.186335; label("Power exponent of (SCORE_MMSE/26) on baseline SUVR (unitless)") # Table 2 MMSEBL on BL = -0.186; Model S2 THETA(10)
    e_apoe4non_rbase <- -0.0403713; label("Fractional change in baseline SUVR for ApoE4 non-carriers vs heterozygotes (unitless)") # Table 2 noncarrier on BL = -0.0404; Model S2 THETA(8)
    e_apoe4hom_rbase <- 0.00895844; label("Fractional change in baseline SUVR for ApoE4 homozygotes vs heterozygotes (unitless)") # Table 2 carrier (2 copies) on BL = 0.00896; Model S2 THETA(9)

    # Model S2 $OMEGA BLOCK(4), eta order Kout, Emax, BSL, EC50, banded with
    # the Kout-BSL, Kout-EC50 and Emax-EC50 covariances fixed to 0 (Results:
    # 'banded covariance structure'). Table 2 prints %CV = 100*sqrt(omega^2)
    # (44.2 / 25.3 / 12.8 / 84.1) and rho = -0.872 / 0.439 / -0.406, which
    # these values reproduce; its footnote b formula does not.
    etalkout + etalemax + etalrbase + etalec50 ~ c(
      0.195612,
      -0.0976648, 0.0641612,
      fixed(0), 0.0141831, 0.0163087,
      fixed(0), fixed(0), -0.04365, 0.707621
    )

    propSd_suvr <- 0.0403608; label("Proportional residual error on composite SUVR (fraction)") # Table 2 proportional error = 4.04%; Model S2 THETA(5); additive THETA(6) 0 FIX (dropped, Table S5 run 75)
  })

  model({
    # Covariate multipliers (Model S1)
    male <- 1 - SEXF
    cl_cov <- (WT / 70.9)^e_wt_cl *
      (1 + e_male_cl * male) *
      (1 + e_race_asian_cl * RACE_ASIAN)
    vc_cov <- (WT / 70.9)^e_wt_vc *
      (1 + e_male_vc * male) *
      (1 + e_race_asian_vc * RACE_ASIAN)
    vp_cov <- (WT / 70.9)^e_wt_vp *
      (AGE / 71)^e_age_vp *
      (SCORE_MMSE / 26)^e_mmse_vp *
      (1 + e_male_vp * male) *
      (1 + e_race_asian_vp * RACE_ASIAN)

    cl <- exp(lcl + etalcl) * cl_cov
    vc <- exp(lvc + etalvc) * vc_cov
    vp <- exp(lvp + etalvp) * vp_cov
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Covariate multipliers (Model S2). The ApoE4 heterozygote is the
    # reference genotype.
    apoe4_non <- 1 - APOE4_HET - APOE4_HOM
    rbase_cov <- (AGE / 71)^e_age_rbase *
      (SCORE_MMSE / 26)^e_mmse_rbase *
      (1 + e_apoe4non_rbase * apoe4_non + e_apoe4hom_rbase * APOE4_HOM)
    rbase <- exp(lrbase + etalrbase) * rbase_cov
    kout <- exp(lkout + etalkout) * (WT / 72)^e_wt_kout
    emax <- exp(lemax + etalemax) * (AGE / 71)^e_age_emax
    ec50 <- exp(lec50 + etalec50)
    kin <- rbase * kout

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc

    # Equation 1: dSUVR/dt = Kin - Kout * (1 + Emax * C / (EC50 + C)) * SUVR
    d/dt(suvr) <- kin - kout * (1 + emax * Cc / (ec50 + Cc)) * suvr
    suvr(0) <- rbase

    # Model S1 $ERROR: W2 = sqrt((THETA(5)*IPRED)^2 + THETA(6)^2) * exp(ETA(4)),
    # i.e. combined2() error with the whole residual SD scaled per subject.
    propSdCc <- propSd * exp(etaruv)
    addSdCc <- addSd * exp(etaruv)
    Cc ~ add(addSdCc) + prop(propSdCc) + combined2()
    suvr ~ prop(propSd_suvr)
  })
}
