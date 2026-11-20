DeCarlo_2026_carboplatin <- function() {
  description <- paste(
    "Population PK-PD model of carboplatin-induced myelosuppression in adult",
    "cancer patients (De Carlo 2026): a fixed two-compartment carboplatin PK",
    "model taken from Ekhart 2006 (clearance from age, weight, sex and serum",
    "creatinine) driving two parallel Friberg myelosuppression chains for",
    "platelets and neutrophils (proliferating pool, three transit compartments",
    "and a circulating pool each, linear drug effect on proliferation), with",
    "serum albumin, CKD-EPI eGFR, paclitaxel and pemetrexed co-medication,",
    "oesophageal/gastric cancer and treatment cycle as covariates."
  )
  reference <- paste(
    "De Carlo A, Tosca EM, Crul M, Schutte T, van Zuijlen L, Bahce I, Said MM,",
    "Buter J, Huls H, Magni P, Bartelink I. Model-informed precision dosing of",
    "carboplatin in cancer patients by leveraging myelosuppression data from",
    "electronic health records. Br J Clin Pharmacol. 2026;92(5):1457.",
    "doi:10.1002/bcp.70413 (PMID 41431252; PMCID PMC13122292).",
    "PK layer: Ekhart C, de Jonge ME, Huitema ADR, Schellens JHM, Rodenhuis S,",
    "Beijnen JH. Flat dosing of carboplatin is justified in adult patients",
    "with normal renal function. Clin Cancer Res. 2006;12(21):6502-6508.",
    "doi:10.1158/1078-0432.CCR-05-1076."
  )
  vignette <- "DeCarlo_2026_carboplatin"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "mg/L",
    PLT = "10^9/L",
    ANC = "10^9/L"
  )

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters the fixed Ekhart 2006 carboplatin clearance equation only. Training-set median 68 years (95% range 39-84; De Carlo 2026 Table 1).",
      source_name = "Age"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying in the source. Enters the fixed Ekhart 2006 carboplatin clearance equation only. Training-set median 73 kg (95% range 50-116; De Carlo 2026 Table 1).",
      source_name = "Weight"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Ekhart 2006 / De Carlo 2026 Eq. S1.1 'Gender = 1 if female'; the female factor 0.736 multiplies the creatinine-dependent clearance arm. 49% female in the training set (Table 1).",
      source_name = "Gender"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying in the source. Enters the fixed Ekhart 2006 carboplatin clearance equation in umol/L. Training-set median 72 umol/L (95% range 45-139; Table 1).",
      source_name = "Creatinine"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying in the source. Power effect (ALB / 44)^beta on the",
        "platelet and neutrophil baselines and on the platelet drug-effect",
        "slope; 44 g/L is the reference chosen by the authors in the middle",
        "of the normal range (Table 2 footnote). Training-set median 37 g/L",
        "(95% range 25-47; Table 1)."
      ),
      source_name = "ALBUMIN"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate by the CKD-EPI creatinine equation",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying in the source. Power effect (CRCL / 100)^beta on the",
        "neutrophil drug-effect slope; 100 is the reference chosen by the",
        "authors within the normal range (Table 2 footnote). The paper labels",
        "the unit 'mL/min' but the CKD-EPI equation returns mL/min/1.73 m^2;",
        "supply the CKD-EPI value as the equation computes it. The main text",
        "cites the 2009 equation, the Table S1.1 footnote the 2021 equation.",
        "Training-set median 85.89 (95% range 37.80-115.91; Table 1). Not",
        "used for the carboplatin PK, which uses CREAT directly."
      ),
      source_name = "GFR CKDEPI"
    ),
    CONMED_PACLITAXEL = list(
      description = "Concomitant paclitaxel indicator (1 = receiving paclitaxel with carboplatin, 0 = not)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant paclitaxel)",
      notes = "Log-linear effects on the platelet mean maturation time and the platelet drug-effect slope. 51% of the training set (Table 1).",
      source_name = "PTX"
    ),
    CONMED_PEMETREXED = list(
      description = "Concomitant pemetrexed indicator (1 = receiving pemetrexed with carboplatin, 0 = not)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant pemetrexed)",
      notes = "Log-linear effect on the neutrophil baseline. 70% of the training set (Table 1).",
      source_name = "Pemetrexed"
    ),
    TUMTP_ESOPHAGOGASTRIC = list(
      description = "Oesophageal or gastric cancer indicator (1 = oesophageal/gastric neoplasm, 0 = lung or gynaecological cancer)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (lung or gynaecological cancer)",
      notes = "Log-linear effect on the platelet baseline. The source pools oesophageal and gastric neoplasms of any histology into one diagnosis group (32% of the training set, Table 1); Eq. S2.1 names the indicator 'Stomach'.",
      source_name = "Stomach"
    ),
    CYCLE = list(
      description = "Carboplatin treatment-cycle number (1 = first cycle)",
      units = "(count)",
      type = "count",
      reference_category = NULL,
      notes = paste(
        "Time-varying. Enters the platelet drug-effect slope as the uncentred",
        "linear factor (1 + 0.756 * CYCLE) of Eq. S2.1, so the tabulated",
        "SlopePLT is the value extrapolated to CYCLE = 0. The paper does not",
        "state the origin of its 'Treat.Cycle' counter; the register's",
        "1-based counter is used (see the vignette's assumptions section)."
      ),
      source_name = "Treat.Cycle"
    )
  )

  covariatesDataExcluded <- list(
    LDH = list(
      description = "Serum lactate dehydrogenase",
      units = "U/L",
      type = "continuous",
      notes = "Screened (Methods 2.3) but not retained in the final model."
    ),
    CONMED_CAPECITABINE = list(
      description = "Concomitant capecitabine indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Table 1) but not retained; 2% of the training set."
    ),
    CONMED_ETOPOSIDE = list(
      description = "Concomitant etoposide indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Table 1) but not retained; 11% of the training set."
    ),
    RADIOTHERAPY = list(
      description = "Concomitant radiotherapy indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Table 1) but not retained; 31% of the training set."
    )
  )

  compartmentData <- list(
    central = list(analyte = "carboplatin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "carboplatin", units = "mg", specimen = "plasma", verified = TRUE),
    prol_plt = list(analyte = "platelets", units = "10^9/L", specimen = "not applicable", verified = TRUE),
    transit1_plt = list(analyte = "platelets", units = "10^9/L", specimen = "not applicable", verified = TRUE),
    transit2_plt = list(analyte = "platelets", units = "10^9/L", specimen = "not applicable", verified = TRUE),
    transit3_plt = list(analyte = "platelets", units = "10^9/L", specimen = "not applicable", verified = TRUE),
    circ_plt = list(analyte = "platelets", units = "10^9/L", specimen = "whole blood", verified = TRUE),
    prol_anc = list(analyte = "neutrophils", units = "10^9/L", specimen = "not applicable", verified = TRUE),
    transit1_anc = list(analyte = "neutrophils", units = "10^9/L", specimen = "not applicable", verified = TRUE),
    transit2_anc = list(analyte = "neutrophils", units = "10^9/L", specimen = "not applicable", verified = TRUE),
    transit3_anc = list(analyte = "neutrophils", units = "10^9/L", specimen = "not applicable", verified = TRUE),
    circ_anc = list(analyte = "neutrophils", units = "10^9/L", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 580L,
    n_studies = 1L,
    age_range = "95% range 39-84 years",
    age_median = "68 years",
    weight_range = "95% range 50-116 kg",
    weight_median = "73 kg",
    sex_female_pct = 49,
    disease_state = paste(
      "Adult cancer patients treated with carboplatin-based chemotherapy for",
      "lung (52%), gynaecological (16%) or oesophageal/gastric (32%) cancer;",
      "61% curative intent, 31% concomitant radiotherapy."
    ),
    dose_range = "Carboplatin IV per the Calvert formula (routine clinical dosing, e.g. AUC 5 every 3 weeks or AUC 2 weekly with radiotherapy).",
    regions = "Netherlands (Amsterdam UMC)",
    co_medication = "Paclitaxel 51%, pemetrexed 70%, etoposide 11%, capecitabine 2%; immunotherapy in 29% (Table 1).",
    notes = paste(
      "Retrospective electronic-health-record data, January 2019 - June 2022",
      "(training set, 580 patients, 3635 platelet and 2726 neutrophil",
      "observations; Table 1). A separate test set of 210 patients treated",
      "June 2022 - January 2024 was used for Bayesian forecasting evaluation",
      "(Table S1.1). No carboplatin concentrations were measured; the PK is",
      "fixed to Ekhart 2006 (178 Dutch patients with creatinine clearance",
      "> 50 mL/min). Monolix 2023R1."
    )
  )

  ini({
    # ---- Carboplatin PK: fixed to Ekhart 2006 (De Carlo 2026 Methods 2.2:
    # 'individual PK parameters were fixed to the population values found in
    # the original work and adjusted by patient covariates'). No PK IIV. ----
    lcl_nonren <- fixed(log(115.5)); label("Creatinine-independent carboplatin clearance intercept (mL/min)") # De Carlo 2026 Eq. S1.1 '115.5'; Ekhart 2006 Results equation 'Cl_carbo (mL/min) = 115.5 + ...' (Table 3 model 4)
    cl_renal_fraction <- fixed(0.212); label("Coefficient of (140 - AGE) * WT / CREAT in carboplatin clearance (mL/min per year.kg/(umol/L))") # De Carlo 2026 Eq. S1.1 '0.212'; Ekhart 2006 Results equation
    e_sexf_cl_renal <- fixed(0.736); label("Multiplicative factor on the creatinine-dependent clearance arm for women (unitless)") # De Carlo 2026 Eq. S1.1 '0.736^Gender'; Ekhart 2006 Results equation
    lvc <- fixed(log(16.3)); label("Central volume of distribution V (L)") # Ekhart 2006 Table 2 'Volume of distribution (L)' = 16.3
    lk12 <- fixed(log(0.104)); label("Distribution rate constant k12, central to peripheral (1/h)") # Ekhart 2006 Table 2 'Distribution microconstant k12' = 0.104 1/h
    lk21 <- fixed(log(0.171)); label("Distribution rate constant k21, peripheral to central (1/h)") # Ekhart 2006 Table 2 'Distribution microconstant k21' = 0.171 1/h

    # ---- Platelet Friberg chain (De Carlo 2026 Table 2) ----
    lcirc0_plt <- log(293); label("Baseline platelet count Circ0PLT (10^9/L)") # Table 2 'Circ0PLT' = 293
    e_tumtp_esophagogastric_circ0_plt <- -0.16; label("Log-linear effect of oesophageal/gastric cancer on Circ0PLT (unitless)") # Table 2 'beta gastric/esophageal on Circ0PLT' = -0.16
    e_alb_circ0_plt <- -0.503; label("Power exponent of (ALB / 44 g/L) on Circ0PLT (unitless)") # Table 2 'beta Albumin on Circ0PLT' = -0.503
    lmtt_plt <- log(6.27); label("Mean maturation time of platelets MTTPLT (day)") # Table 2 'MTTPLT' = 6.27 days
    e_conmed_paclitaxel_mtt_plt <- 0.306; label("Log-linear effect of paclitaxel co-medication on MTTPLT (unitless)") # Table 2 'beta PTX on MTTPLT' = 0.306
    lgamma_plt <- log(0.188); label("Platelet feedback exponent gammaPLT (unitless)") # Table 2 'gammaPLT' = 0.188
    lslope_plt <- log(0.282); label("Linear carboplatin effect on platelet proliferation SlopePLT at CYCLE = 0 (L/mg)") # Table 2 'SlopePLT' = 0.282 L/mg
    e_cycle_slope_plt <- 0.756; label("Linear effect of treatment cycle on SlopePLT, factor (1 + beta * CYCLE) (unitless)") # Table 2 'beta Treatment on SlopePLT' = 0.756; Eq. S2.1
    e_alb_slope_plt <- -0.26; label("Power exponent of (ALB / 44 g/L) on SlopePLT (unitless)") # Table 2 'beta Albumin on SlopePLT' = -0.26
    e_conmed_paclitaxel_slope_plt <- -0.562; label("Log-linear effect of paclitaxel co-medication on SlopePLT (unitless)") # Table 2 'beta PTX on SlopePLT' = -0.562

    # ---- Neutrophil Friberg chain (De Carlo 2026 Table 2) ----
    lcirc0_anc <- log(5.16); label("Baseline neutrophil count Circ0NT (10^9/L)") # Table 2 'Circ0NT' = 5.16
    e_alb_circ0_anc <- -0.471; label("Power exponent of (ALB / 44 g/L) on Circ0NT (unitless)") # Table 2 'beta Albumin on Circ0NT' = -0.471
    e_conmed_pemetrexed_circ0_anc <- 0.355; label("Log-linear effect of pemetrexed co-medication on Circ0NT (unitless)") # Table 2 'beta Pemetrexed on Circ0NT' = 0.355
    lmtt_anc <- log(6.24); label("Mean maturation time of neutrophils MTTNT (day)") # Table 2 'MTTNT' = 6.24 days
    lgamma_anc <- log(0.074); label("Neutrophil feedback exponent gammaNT (unitless)") # Table 2 'gammaNT' = 0.074
    lslope_anc <- log(0.509); label("Linear carboplatin effect on neutrophil proliferation SlopeNT (L/mg)") # Table 2 'SlopeNT' = 0.509 L/mg
    e_crcl_slope_anc <- -0.602; label("Power exponent of (CRCL / 100) on SlopeNT (unitless)") # Table 2 'beta GFR CKDEPI on SlopeNT' = -0.602

    # ---- IIV: Table 2 reports omega as SDs; variance = omega^2. The second
    # block of four omega rows is mislabelled 'PLT' in the published table;
    # it is the NT block in the same order as the PLT block (Circ0, MTT,
    # gamma, Slope), matching the eight random effects of Eq. S2.1. ----
    etalcirc0_plt + etalcirc0_anc ~ c(0.071824, 0.0384419, 0.106276) # Table 2 omega Circ0PLT = 0.268, omega Circ0NT = 0.326 (2nd block row 1), rho = 0.440; cov = 0.440 * 0.268 * 0.326
    etalmtt_plt ~ 0.0361 # Table 2 omega MTTPLT = 0.190
    etalgamma_plt ~ 0.128164 # Table 2 omega gammaPLT = 0.358
    etalslope_plt ~ 0.299209 # Table 2 omega SlopePLT = 0.547
    etalmtt_anc ~ 0.103684 # Table 2 omega MTTNT = 0.322 (2nd block row 2)
    etalgamma_anc + etalslope_anc ~ c(0.391876, 0.283202, 0.389376) # Table 2 omega gammaNT = 0.626 (2nd block row 3), omega SlopeNT = 0.624 (2nd block row 4), rho(SlopeNT, gammaNT) = 0.725; cov = 0.725 * 0.626 * 0.624

    # ---- Residual error (Table 2) ----
    addSd_PLT <- 5.92; label("Additive residual error on platelets (10^9/L)") # Table 2 'sigma PLT, Additive' = 5.92
    propSd_PLT <- 0.218; label("Proportional residual error on platelets (fraction)") # Table 2 'sigma PLT, Proportional' = 0.218
    propSd_ANC <- 0.321; label("Proportional residual error on neutrophils (fraction)") # Table 2 'sigma NT, Proportional' = 0.321
  })

  model({
    # Carboplatin clearance, Ekhart 2006 model 4 as reproduced in De Carlo
    # 2026 Eq. S1.1, in mL/min with CREAT in umol/L; converted to L/day
    # (x 1.44). The rate constants are reported per hour (x 24 to 1/day).
    cl_mlmin <- exp(lcl_nonren) + cl_renal_fraction * (140 - AGE) * WT * e_sexf_cl_renal^SEXF / CREAT
    cl <- cl_mlmin * 1.44
    vc <- exp(lvc)
    k12 <- exp(lk12) * 24
    k21 <- exp(lk21) * 24
    kel <- cl / vc

    # Platelet chain parameters (Eq. S2.1)
    circ0_plt <- exp(lcirc0_plt + e_tumtp_esophagogastric_circ0_plt * TUMTP_ESOPHAGOGASTRIC + e_alb_circ0_plt * log(ALB / 44) + etalcirc0_plt)
    mtt_plt <- exp(lmtt_plt + e_conmed_paclitaxel_mtt_plt * CONMED_PACLITAXEL + etalmtt_plt)
    gamma_plt <- exp(lgamma_plt + etalgamma_plt)
    slope_plt <- exp(lslope_plt + e_alb_slope_plt * log(ALB / 44) + e_conmed_paclitaxel_slope_plt * CONMED_PACLITAXEL + etalslope_plt) * (1 + e_cycle_slope_plt * CYCLE)
    ktr_plt <- 4 / mtt_plt # (N + 1) / MTT with N = 3 transit compartments

    # Neutrophil chain parameters (Eq. S2.1)
    circ0_anc <- exp(lcirc0_anc + e_conmed_pemetrexed_circ0_anc * CONMED_PEMETREXED + e_alb_circ0_anc * log(ALB / 44) + etalcirc0_anc)
    mtt_anc <- exp(lmtt_anc + etalmtt_anc)
    gamma_anc <- exp(lgamma_anc + etalgamma_anc)
    slope_anc <- exp(lslope_anc + e_crcl_slope_anc * log(CRCL / 100) + etalslope_anc)
    ktr_anc <- 4 / mtt_anc

    # Carboplatin PK (Eq. S1.1)
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    Cc <- central / vc

    # Linear drug effect on proliferation (Eq. S1.2), kprol = ktr = kcirc
    edrug_plt <- slope_plt * Cc
    edrug_anc <- slope_anc * Cc

    d/dt(prol_plt) <- ktr_plt * prol_plt * (1 - edrug_plt) * (circ0_plt / circ_plt)^gamma_plt - ktr_plt * prol_plt
    d/dt(transit1_plt) <- ktr_plt * prol_plt - ktr_plt * transit1_plt
    d/dt(transit2_plt) <- ktr_plt * transit1_plt - ktr_plt * transit2_plt
    d/dt(transit3_plt) <- ktr_plt * transit2_plt - ktr_plt * transit3_plt
    d/dt(circ_plt) <- ktr_plt * transit3_plt - ktr_plt * circ_plt

    d/dt(prol_anc) <- ktr_anc * prol_anc * (1 - edrug_anc) * (circ0_anc / circ_anc)^gamma_anc - ktr_anc * prol_anc
    d/dt(transit1_anc) <- ktr_anc * prol_anc - ktr_anc * transit1_anc
    d/dt(transit2_anc) <- ktr_anc * transit1_anc - ktr_anc * transit2_anc
    d/dt(transit3_anc) <- ktr_anc * transit2_anc - ktr_anc * transit3_anc
    d/dt(circ_anc) <- ktr_anc * transit3_anc - ktr_anc * circ_anc

    prol_plt(0) <- circ0_plt
    transit1_plt(0) <- circ0_plt
    transit2_plt(0) <- circ0_plt
    transit3_plt(0) <- circ0_plt
    circ_plt(0) <- circ0_plt

    prol_anc(0) <- circ0_anc
    transit1_anc(0) <- circ0_anc
    transit2_anc(0) <- circ0_anc
    transit3_anc(0) <- circ0_anc
    circ_anc(0) <- circ0_anc

    PLT <- circ_plt
    ANC <- circ_anc

    PLT ~ add(addSd_PLT) + prop(propSd_PLT) + combined1()
    ANC ~ prop(propSd_ANC)
  })
}
