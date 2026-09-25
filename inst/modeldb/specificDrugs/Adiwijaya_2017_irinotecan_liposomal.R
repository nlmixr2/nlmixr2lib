Adiwijaya_2017_irinotecan_liposomal <- function() {
  description <- paste(
    "Coupled two-analyte population PK model for nanoliposomal irinotecan",
    "(nal-IRI, MM-398/PEP02; output Cc = total irinotecan) and total SN-38",
    "(output Cc_sn38) in adults with metastatic pancreatic, gastric/GEJ,",
    "colorectal and other solid tumours (Adiwijaya 2017, pooled from six",
    "phase I-III studies including NAPOLI-1). Total irinotecan is a",
    "two-compartment model with first-order elimination. Total SN-38 is the",
    "sum of encapsulated SN-38 (a time-invariant mass fraction of total",
    "irinotecan, i.e. a co-encapsulated manufacturing contaminant) and",
    "unencapsulated SN-38 (a one-compartment model formed from total",
    "irinotecan by the first-order rate constant kmet and sharing the",
    "irinotecan central volume). The two sub-models are sequential: the",
    "irinotecan model is independent, and the SN-38 model depends on it,",
    "including through the individual irinotecan CL and V1 acting as",
    "covariates on the formation rate.",
    sep = " "
  )
  reference <- paste(
    "Adiwijaya BS, Kim J, Lang I, Csoszi T, Cubillo A, Chen JS, Wong M,",
    "Park JO, Kim JS, Rau KM, Melichar B, Gallego JB, Fitzgerald J,",
    "Belanger B, Molnar I, Ma WW. Population Pharmacokinetics of Liposomal",
    "Irinotecan in Patients With Cancer. Clin Pharmacol Ther.",
    "2017;102(6):997-1005. doi:10.1002/cpt.720. PMCID: PMC5697569.",
    sep = " "
  )
  vignette <- "Adiwijaya_2017_irinotecan_liposomal"
  units <- list(time = "week", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    BSA = list(
      description = "Body surface area; nal-IRI is dosed per BSA",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters V1 and the SN-38 formation rate kmet as a median-centred",
        "linear term inside exp() per supplement Equation S1, i.e.",
        "exp(theta * (BSA - 1.70)). Centring constant 1.70 m^2 is the cohort",
        "median from Table 1. The BSA formula used is not stated in the source."
      ),
      source_name = "BSA"
    ),
    RACE_ASIAN = list(
      description = "East Asian race indicator (1 = East Asian, 0 = Caucasian or Other)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian: Caucasian 52% or Other 6% of the cohort)",
      notes = paste(
        "The source's race categories are Caucasian / East Asian / Others",
        "(Table 1) and the model theta is labelled 'race == Asian' with the",
        "reference 'race == non-Asian' (Tables S6, S7), so the two non-Asian",
        "categories are pooled into the reference. Race was the single most",
        "significant baseline factor for both analytes: East Asian patients",
        "had higher irinotecan CL and lower SN-38 CL."
      ),
      source_name = "race"
    ),
    LMET = list(
      description = "Baseline presence of liver metastases",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no liver metastases)",
      notes = paste(
        "Recorded only in NAPOLI-1 (66% positive, Table 1); the source imputed",
        "'No' for the other five studies and reports that a sensitivity",
        "analysis evaluated that imputation (Methods, PK data)."
      ),
      source_name = "Liver metastasis"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Log10-transformed and median-centred: exp(theta * (log10(ALT) -",
        "log10(25))), median 25 U/L from Table 1. The source log-transformed",
        "all laboratory covariates because log-normal distributions were",
        "observed (Methods, PK data). Because the term is centred, the",
        "coefficient is invariant to the concentration unit."
      ),
      source_name = "ALT"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Log10-transformed and median-centred: exp(theta * (log10(ALB) -",
        "log10(40))), median 40 g/L from Table 1 (already SI). Because the",
        "term is centred, the coefficient is invariant to the unit, so no",
        "g/L -> g/dL conversion is required. The source notes the albumin -",
        "total-irinotecan association runs opposite to the direction expected",
        "from hepatic impairment and is unlikely to be clinically relevant."
      ),
      source_name = "albumin"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Log10-transformed and median-centred: exp(theta * (log10(TBILI) -",
        "log10(7))), median 7 umol/L from Table 1 (already SI). Because the",
        "term is centred, the coefficient is invariant to the unit, so no",
        "umol/L -> mg/dL conversion is required. Bilirubin is the strongest",
        "laboratory predictor of SN-38: patients with bilirubin >= 1 mg/dL",
        "(17.1 umol/L) had 43% higher uSN38 Cavg and 35% higher uSN38 Cmax.",
        "Only 20 of 353 patients had bilirubin >= 1 mg/dL."
      ),
      source_name = "bilirubin"
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters as a median-centred LINEAR (not log) term inside exp():",
        "exp(theta * (CRCL - 81.6)); the Tables S6/S7 theta unit 'min/mL'",
        "fixes the covariate scale as mL/min. Table 1 reports the median as",
        "1.36e-3 L/s = 81.6 mL/min. The source does not state whether its",
        "creatinine clearance was BSA-normalised; because the cohort median",
        "BSA is 1.70 m^2, the absolute and 1.73 m^2-normalised values coincide",
        "to within 2%, so the distinction does not materially affect the",
        "covariate term. Creatinine clearance was not a significant predictor",
        "of SN-38 after adjusting for BSA."
      ),
      source_name = "CrCl"
    ),
    STUDY_NAPOLI1 = list(
      description = "NAPOLI-1 study / NAPOLI-1 drug-product manufacturing site indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (PharmaEngine (PEI) phase I-II studies / PEI manufacturing site)",
      notes = paste(
        "Table S3 names this covariate 'manufacturing site' with the two",
        "levels 'NAPOLI' and 'PEI', tested because nal-IRI used in the phase",
        "III study was manufactured differently from the material used in the",
        "phase I-II studies. The manufacturing site is perfectly collinear",
        "with the study (NAPOLI-1 = 258 of 353 patients, Table 1), so the",
        "canonical study indicator carries it. It acts on V1, irinotecan CL,",
        "the SN-38 formation rate kmet, and the encapsulated-SN-38 fraction."
      ),
      source_name = "mfg"
    ),
    CONMED_FLUOROURACIL = list(
      description = "Concomitant fluorouracil (5-FU, given with leucovorin) coadministration indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (nal-IRI monotherapy)",
      notes = paste(
        "Tested for a drug-drug interaction (Table S3). In NAPOLI-1, 116 of",
        "258 patients (45%) received nal-IRI + 5-FU/LV (Table 1). The",
        "estimated effect on both clearances is small."
      ),
      source_name = "Treatment contains 5FU"
    ),
    UGT1A1_STAR28_HOM = list(
      description = "UGT1A1*28 7/7 homozygous-variant genotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (UGT1A1*28 non-7/7: 6/6 or 6/7)",
      notes = paste(
        "Genotyped in NAPOLI-1 only; 14 of 258 patients (5%) were 7/7",
        "homozygous (Table 1). The estimated effect on SN-38 clearance is",
        "-1.46e-05 on the log scale, i.e. numerically indistinguishable from",
        "zero: the source states the 7/7 clearance is 1.000-times the non-7/7",
        "clearance (0.0% difference) and concludes UGT1A1*28 is not a",
        "significant covariate with nal-IRI, in contrast to nonliposomal",
        "irinotecan. The coefficient is retained here because the source's",
        "full-covariate approach retained it in the final model."
      ),
      source_name = "UGT1A1*28"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Pre-specified as a covariate on both clearances in Table S3, but no",
        "coefficient is reported for age in either final-model parameter",
        "table (Table S6, Table S7), so it cannot be encoded. Cohort median",
        "63 years (5th-95th percentile 39.8-79.2), Table 1."
      )
    ),
    SEXF = list(
      description = "Sex",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Pre-specified as a covariate on both clearances in Table S3, but no",
        "coefficient is reported in Table S6 or Table S7. The Results state",
        "that sex was not significantly associated with SN-38 after adjusting",
        "for BSA. Cohort 44% female (Table 1)."
      )
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = paste(
        "Pre-specified as a hepatic-function covariate on both clearances in",
        "Table S3, but no coefficient is reported in Table S6 or Table S7",
        "(ALT, albumin and bilirubin are the hepatic markers that are).",
        "Cohort median 29 U/L (5th-95th percentile 14.7-81.9), Table 1."
      )
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "irinotecan",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "irinotecan",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_sn38 = list(
      analyte = "SN-38",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 353,
    n_studies = 6,
    age_median = "63 years",
    age_range = "39.8-79.2 years (5th-95th percentile)",
    sex_female_pct = 44,
    race_ethnicity = c(Caucasian = 52, `East Asian` = 42, Other = 6),
    disease_state = paste(
      "Advanced solid tumours: metastatic pancreatic cancer 73%,",
      "gastric and gastroesophageal-junction cancer 10%, solid tumour 11%,",
      "colorectal cancer 5%"
    ),
    dose_range = paste(
      "50-150 mg/m^2 irinotecan free base (60-180 mg/m^2 irinotecan",
      "hydrochloride trihydrate salt) intravenously; most patients received",
      "100 mg/m^2 Q3W (53%) or 70 mg/m^2 Q2W (40%)"
    ),
    bsa_median = "1.7 m^2 (5th-95th percentile 1.3-2.2)",
    hepatic_function = paste(
      "Albumin 40 g/L (29-47), ALT 25 U/L (8.9-96.3), AST 29 U/L (14.7-81.9),",
      "total bilirubin 7 umol/L (3-19), all median (5th-95th percentile).",
      "Hepatic impairment was an exclusion criterion, but 20 patients were",
      "enrolled with bilirubin >= 1 mg/dL."
    ),
    renal_function = paste(
      "Creatinine clearance 1.36e-3 L/s = 81.6 mL/min",
      "(5th-95th percentile 0.66-2.53e-3 L/s = 39.6-151.8 mL/min).",
      "Renal impairment was an exclusion criterion."
    ),
    co_medication = "45% of NAPOLI-1 patients received nal-IRI with 5-fluorouracil and leucovorin",
    regions = "Europe, East Asia, Australia, North America",
    notes = paste(
      "Baseline characteristics from Table 1; the six contributing studies",
      "are listed in Table S1. PK samples were collected during cycle 1 only",
      "(intensive sampling in the five phase I-II studies, sparse sampling in",
      "phase III NAPOLI-1): 1,792 total-irinotecan and 1,765 total-SN-38",
      "concentrations. Doses throughout this model file are expressed as",
      "irinotecan free base, matching the source's convention."
    )
  )

  ini({
    # --- Total irinotecan (tIRI): structural, Table S 6 'Estimated values (final model)'
    lvc <- log(4.60)
    label("Central volume of total irinotecan V1 (L)") # Adiwijaya 2017 Table S 6, Volume (V1) = 4.60 L
    lcl <- log(13.6)
    label("Clearance of total irinotecan CL (L/week)") # Adiwijaya 2017 Table S 6, Clearance (CL) = 13.6 L/week
    lq <- log(0.471)
    label("Intercompartmental clearance of total irinotecan Q (L/week)") # Adiwijaya 2017 Table S 6, Q = 0.471 L/week
    lvp <- log(48.7)
    label("Peripheral volume of total irinotecan V2 (L)") # Adiwijaya 2017 Table S 6, V2 = 48.7 L

    # --- Total irinotecan: covariate effects (supplement Equation S1; log-scale,
    #     continuous covariates median-centred, categorical as 0/1 indicators)
    e_bsa_vc <- 0.416
    label("Effect of BSA on V1 (per m^2, median-centred at 1.70 m^2)") # Adiwijaya 2017 Table S 6, theta{V1,BSA} = 0.416 1/m2
    e_study_napoli1_vc <- -0.172
    label("Effect of NAPOLI-1 study / manufacturing site on V1") # Adiwijaya 2017 Table S 6, theta{V1,mfg==NAPOLI} = -0.172 (relative to mfg==PEI)
    e_race_asian_cl <- 0.647
    label("Effect of East Asian race on irinotecan CL") # Adiwijaya 2017 Table S 6, theta{CL,race==Asian} = 0.647 (relative to race==non-Asian)
    e_conmed_fluorouracil_cl <- 0.075
    label("Effect of 5-FU coadministration on irinotecan CL") # Adiwijaya 2017 Table S 6, theta{CL,treatment contain 5FU} = 0.075
    e_study_napoli1_cl <- -0.189
    label("Effect of NAPOLI-1 study / manufacturing site on irinotecan CL") # Adiwijaya 2017 Table S 6, theta{CL,mfg==NAPOLI} = -0.189
    e_lmet_cl <- -0.075
    label("Effect of baseline liver metastases on irinotecan CL") # Adiwijaya 2017 Table S 6, theta{CL,liver metastasis} = -0.075
    e_alt_cl <- 0.016
    label("Effect of log10 ALT on irinotecan CL (median-centred at 25 U/L)") # Adiwijaya 2017 Table S 6, theta{CL,ALT} = 0.016 per unit change of log10 ALT
    e_alb_cl <- -1.79
    label("Effect of log10 albumin on irinotecan CL (median-centred at 40 g/L)") # Adiwijaya 2017 Table S 6, theta{CL,albumin} = -1.79 per unit change of log10 albumin
    e_tbili_cl <- -0.0670
    label("Effect of log10 total bilirubin on irinotecan CL (median-centred at 7 umol/L)") # Adiwijaya 2017 Table S 6, theta{CL,bilirubin} = -0.0670 per unit change of log10 bilirubin
    e_crcl_cl <- 0.003
    label("Effect of creatinine clearance on irinotecan CL (per mL/min, median-centred at 81.6 mL/min)") # Adiwijaya 2017 Table S 6, theta{CL,creatinine clearance} = 0.003 min/mL

    # --- SN-38: structural, Table S 7 'Estimated Values (Final Model)'
    lcl_sn38 <- log(14.2)
    label("Clearance of unencapsulated SN-38 CL_SN38 (L/week)") # Adiwijaya 2017 Table S 7, Clearance (CL_SN38) = 14.2 L/week
    lkmet <- log(0.00072)
    label("First-order formation rate constant of unencapsulated SN-38 from total irinotecan Kcov (1/week)") # Adiwijaya 2017 Table S 7, Conversion flux rate from irinotecan (Kcov) = 0.00072 1/week
    lfsn38 <- log(9.0e-5)
    label("Encapsulated-SN-38 to total-irinotecan mass ratio fSN38 (unitless, w/w)") # Adiwijaya 2017 Table S 7, fSN38 = 0.090 ng (SN-38)/ug (irinotecan) = 9.0e-5 w/w; consistent with the Methods statement that the eSN38 fraction of tIRI was estimated to be 0.01% (in vitro 0.015%)

    # --- SN-38: covariate effects (supplement Equation S1, as above)
    e_race_asian_cl_sn38 <- -0.161
    label("Effect of East Asian race on SN-38 CL") # Adiwijaya 2017 Table S 7, theta{CL_SN38,race==Asian} = -0.161 (relative to race == non Asian)
    e_ugt1a1_star28_hom_cl_sn38 <- -1.46e-05
    label("Effect of UGT1A1*28 7/7 homozygosity on SN-38 CL") # Adiwijaya 2017 Table S 7, theta{CL_SN38,UGT1A1*28==homozygous} = -1.46E-05 (relative to non-homozygous)
    e_conmed_fluorouracil_cl_sn38 <- -1.53e-04
    label("Effect of 5-FU coadministration on SN-38 CL") # Adiwijaya 2017 Table S 7, theta{CL_SN38,treatment contains 5FU} = -1.53E-04
    e_lmet_cl_sn38 <- -0.002
    label("Effect of baseline liver metastases on SN-38 CL") # Adiwijaya 2017 Table S 7, theta{CL_SN38,liver metastasis = YES} = -0.002
    e_alt_cl_sn38 <- -1.09e-05
    label("Effect of log10 ALT on SN-38 CL (median-centred at 25 U/L)") # Adiwijaya 2017 Table S 7, theta{CL_SN38,ALT} = -1.09E-05 per unit change of log10 ALT
    e_alb_cl_sn38 <- -0.207
    label("Effect of log10 albumin on SN-38 CL (median-centred at 40 g/L)") # Adiwijaya 2017 Table S 7, theta{CL_SN38,albumin} = -0.207 per unit change of log10 albumin
    e_tbili_cl_sn38 <- -0.852
    label("Effect of log10 total bilirubin on SN-38 CL (median-centred at 7 umol/L)") # Adiwijaya 2017 Table S 7, theta{CL_SN38,bilirubin} = -0.852 per unit change of log10 bilirubin
    e_crcl_cl_sn38 <- -5.64e-05
    label("Effect of creatinine clearance on SN-38 CL (per mL/min, median-centred at 81.6 mL/min)") # Adiwijaya 2017 Table S 7, theta{CL_SN38,CRCL} = -5.64E-05 min/mL
    e_study_napoli1_kmet <- 3.79e-05
    label("Effect of NAPOLI-1 study / manufacturing site on the SN-38 formation rate") # Adiwijaya 2017 Table S 7, theta{Kcov,mfg==NAPOLI} = 3.79E-05
    e_cl_kmet <- 2.095
    label("Effect of log10 individual irinotecan CL on the SN-38 formation rate (centred at the typical 13.6 L/week)") # Adiwijaya 2017 Table S 7, theta{Kcov,tIRI_logCL} = 2.095 per unit change of log10 tIRI CL
    e_vc_kmet <- -0.867
    label("Effect of log10 individual irinotecan V1 on the SN-38 formation rate (centred at the typical 4.60 L)") # Adiwijaya 2017 Table S 7, theta{Kcov,tIRI_logV1} = -0.867 per unit change of log10 tIRI V1
    e_bsa_kmet <- -1.121
    label("Effect of BSA on the SN-38 formation rate (per m^2, median-centred at 1.70 m^2)") # Adiwijaya 2017 Table S 7, theta{Kcov,BSA} = -1.121 1/m2
    e_study_napoli1_fsn38 <- -0.615
    label("Effect of NAPOLI-1 study / manufacturing site on the encapsulated-SN-38 fraction") # Adiwijaya 2017 Table S 7, theta{fSN38,mfg==NAPOLI} = -0.615

    # --- Inter-individual variability (exponential, per Equation S1)
    # Table S 6 reports the V1-CL pair as a 2x2 block; the implied correlation is
    # 0.184 / sqrt(0.068 * 0.843) = 0.77, which is admissible, confirming the
    # 'unitless (variance)' column label is to be read as variances.
    etalvc + etalcl ~ c(
      0.068,
      0.184, 0.843
    ) # Adiwijaya 2017 Table S 6, Random effects: sigma2(V1) = 0.068, sigma2(V1-CL)(off-diagonal) = 0.184, sigma2(CL) = 0.843
    etalcl_sn38 ~ 0.155 # Adiwijaya 2017 Table S 7, Random effects: sigma2(CL_SN38) = 0.155
    etalkmet ~ 0.184 # Adiwijaya 2017 Table S 7, Random effects: sigma2(Kcov) = 0.184
    etalfsn38 ~ 0.500 # Adiwijaya 2017 Table S 7, Random effects: sigma2(fSN38) = 0.500

    # --- Residual error
    # The source modelled residual variability as additive on the log10 scale
    # (Laplacian estimation with M3 handling of values below the limit of
    # quantification). Converting the reported log10-scale variance to a
    # natural-log SD for rxode2's lnorm() error model:
    # sqrt(0.038) * log(10) = 0.4488 and sqrt(0.021) * log(10) = 0.3337.
    expSd <- 0.4488
    label("Log-normal residual SD for total irinotecan (natural-log scale)") # Adiwijaya 2017 Table S 6, Residuals: sigma2 (in log10 concentration) = 0.038; sqrt(0.038) * log(10) = 0.4488
    expSd_sn38 <- 0.3337
    label("Log-normal residual SD for total SN-38 (natural-log scale)") # Adiwijaya 2017 Table S 7, Residuals: sigma2 (in log10 concentration) = 0.021; sqrt(0.021) * log(10) = 0.3337
  })

  model({
    # 1. Derived covariate terms. Supplement Equation S1 centres every
    #    continuous covariate on its population median and enters categorical
    #    covariates as 0/1 indicators, all inside a single exp().
    #    Centring constants are the Table 1 medians. The three log10 terms are
    #    differences of logarithms, so they are invariant to the unit in which
    #    the laboratory value is supplied.
    bsa_c <- BSA - 1.70
    alt_c <- log10(ALT) - log10(25)
    alb_c <- log10(ALB) - log10(40)
    tbili_c <- log10(TBILI) - log10(7)
    crcl_c <- CRCL - 81.6

    # 2. Individual parameters, total irinotecan
    vc <- exp(
      lvc + etalvc +
        e_bsa_vc * bsa_c +
        e_study_napoli1_vc * STUDY_NAPOLI1
    )
    cl <- exp(
      lcl + etalcl +
        e_race_asian_cl * RACE_ASIAN +
        e_conmed_fluorouracil_cl * CONMED_FLUOROURACIL +
        e_study_napoli1_cl * STUDY_NAPOLI1 +
        e_lmet_cl * LMET +
        e_alt_cl * alt_c +
        e_alb_cl * alb_c +
        e_tbili_cl * tbili_c +
        e_crcl_cl * crcl_c
    )
    vp <- exp(lvp)
    q <- exp(lq)

    # 3. Individual parameters, SN-38. The source used the individual
    #    irinotecan CL and V1 as covariates on the formation rate, because
    #    irinotecan clearance is a surrogate for mononuclear-phagocyte-system
    #    activity, which is hypothesised also to drive release of irinotecan
    #    from the liposome and its conversion to SN-38 (Table S3 rationale,
    #    supplement 'Final PK model of total SN-38'). They are centred here on
    #    the typical values so that the printed Kcov is the typical-subject
    #    value, which is the reading that reproduces the published uSN38 Cmax.
    cl_sn38 <- exp(
      lcl_sn38 + etalcl_sn38 +
        e_race_asian_cl_sn38 * RACE_ASIAN +
        e_ugt1a1_star28_hom_cl_sn38 * UGT1A1_STAR28_HOM +
        e_conmed_fluorouracil_cl_sn38 * CONMED_FLUOROURACIL +
        e_lmet_cl_sn38 * LMET +
        e_alt_cl_sn38 * alt_c +
        e_alb_cl_sn38 * alb_c +
        e_tbili_cl_sn38 * tbili_c +
        e_crcl_cl_sn38 * crcl_c
    )
    kmet <- exp(
      lkmet + etalkmet +
        e_study_napoli1_kmet * STUDY_NAPOLI1 +
        e_cl_kmet * (log10(cl) - log10(13.6)) +
        e_vc_kmet * (log10(vc) - log10(4.60)) +
        e_bsa_kmet * bsa_c
    )
    fsn38 <- exp(
      lfsn38 + etalfsn38 +
        e_study_napoli1_fsn38 * STUDY_NAPOLI1
    )

    # 4. Micro-constants. The unencapsulated-SN-38 compartment shares the
    #    irinotecan central volume: 'Parameter V is obtained from the V1
    #    estimates of the tIRI model' (supplement, final SN-38 model).
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_sn38 <- cl_sn38 / vc

    # 5. ODE system (Figure S1). Distribution of unencapsulated irinotecan to
    #    peripheral tissues is assumed negligible, and release from the
    #    liposome plus conversion to SN-38 are lumped into the single
    #    first-order step kmet.
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_sn38) <- kmet * central - kel_sn38 * central_sn38

    # 6. Observations. All concentrations are in ug/mL (= mg/L), the unit in
    #    which Table 2 reports total irinotecan. The source reports both SN-38
    #    quantities in ng/mL, so multiply Cc_sn38, Cu_sn38 and Ce_sn38 by 1000
    #    when comparing against Table 2 or Figure 1.
    #    Cu_sn38 is the unencapsulated SN-38 that drives the paper's
    #    exposure-response analyses (Cmax for neutropenia, time above a
    #    0.03 ng/mL threshold for survival); Ce_sn38 is the encapsulated
    #    contaminant, which is measured but not eliminated or glucuronidated.
    Cc <- central / vc
    Cu_sn38 <- central_sn38 / vc
    Ce_sn38 <- fsn38 * Cc
    Cc_sn38 <- Cu_sn38 + Ce_sn38

    Cc ~ lnorm(expSd)
    Cc_sn38 ~ lnorm(expSd_sn38)
  })
}
