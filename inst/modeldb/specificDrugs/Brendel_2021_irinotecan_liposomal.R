Brendel_2021_irinotecan_liposomal <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for liposomal irinotecan",
    "(nal-IRI; output Cc = total, i.e. encapsulated plus unencapsulated,",
    "irinotecan) and SN-38 (output Cc_sn38) in adults with solid tumours,",
    "mainly metastatic pancreatic ductal adenocarcinoma (Brendel 2021,",
    "N = 440 pooled from seven phase I-III studies including NAPOLI-1).",
    "Total irinotecan is a two-compartment model with first-order",
    "elimination. SN-38 is a one-compartment model sharing the irinotecan",
    "central volume and is formed from the irinotecan central compartment by",
    "two first-order pathways: a direct one (about 9% of irinotecan",
    "elimination) and a delayed one through a single transit compartment",
    "(about 35%). Covariates: Asian race, sex, oxaliplatin",
    "co-administration and drug-product manufacturing site on irinotecan",
    "clearance; body surface area, sex and manufacturing site on the",
    "irinotecan central volume; manufacturing site on the delayed-pathway",
    "fraction; bilirubin, creatinine clearance, sex and oxaliplatin on SN-38",
    "clearance. Time is in weeks.",
    sep = " "
  )
  reference <- paste(
    "Brendel K, Bekaii-Saab T, Boland PM, Dayyani F, Dean A, Macarulla T,",
    "Maxwell F, Mody K, Pedret-Dunn A, Wainberg ZA, Zhang B. Population",
    "pharmacokinetics of liposomal irinotecan in patients with cancer and",
    "exposure-safety analyses in patients with metastatic pancreatic cancer.",
    "CPT Pharmacometrics Syst Pharmacol. 2021;10:1550-1563.",
    "doi:10.1002/psp4.12725. PMCID: PMC8674005.",
    "Parameter values from Table 2 (identical to Supplementary Table S3);",
    "model structure from the final NONMEM control stream in the",
    "Supplementary Material.",
    sep = " "
  )
  vignette <- "Brendel_2021_irinotecan_liposomal"
  units <- list(
    time = "week",
    dosing = "mg (irinotecan free base)",
    concentration = "ug/mL (Cc, total irinotecan); ng/mL (Cc_sn38, SN-38)"
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on the irinotecan central volume, (BSA/1.71)^0.573;",
        "1.71 m^2 is the cohort median (Supplementary Table S1: median 1.71,",
        "range 1.29-2.48). The BSA formula is not stated in the source."
      ),
      source_name = "BSA"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male; the 'most common' category in the control stream)",
      notes = paste(
        "Control stream codes SEX = 1 as the effect category, and the Results",
        "state that clearances of total irinotecan and SN-38 were 20% lower in",
        "women than in men, so SEX = 1 is female and maps to SEXF with no",
        "transformation. 48.6% of the 440 patients were female",
        "(Supplementary Table S1; 226/440 = 51% men in the Results)."
      ),
      source_name = "SEX"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = any other race)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "35.2% of the cohort was Asian (Supplementary Table S1).",
      source_name = "ASIAN"
    ),
    CONMED_OXALIPLATIN = list(
      description = "Oxaliplatin co-administration indicator (1 = liposomal irinotecan given with oxaliplatin, 0 = without)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no oxaliplatin; liposomal irinotecan alone or with 5-FU/LV)",
      notes = paste(
        "Only the phase I/II first-line mPDAC study (NCT02551991) gave",
        "oxaliplatin, always together with 5-FU/LV (the NALIRIFOX regimen);",
        "12.7% of the cohort (Supplementary Table S1). Raises irinotecan",
        "clearance by 33.9% and lowers SN-38 clearance by 34.4%. The",
        "control-stream column is TRTOXA."
      ),
      source_name = "TRTOXA"
    ),
    FORM_NALIRI_PREVSITE = list(
      description = "Liposomal irinotecan drug-product manufacturing site (1 = previous site, 0 = current commercial site)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (current site; used in all pivotal studies including NAPOLI-1 and NCT02551991)",
      notes = paste(
        "Control stream column MFG: MFG = 1 is the current site (the 'most",
        "common' category, 81.4% of patients) and MFG = 0 the previous site",
        "(18.6%), so FORM_NALIRI_PREVSITE = 1 - MFG. Material from the",
        "previous site was used in phase I studies conducted before 2012",
        "(Results). The previous-site product has a 51.5% higher irinotecan",
        "clearance, a 12.8% lower irinotecan central volume and a 37.6%",
        "higher delayed-pathway fraction ratio. Set to 0 to simulate the",
        "commercial product."
      ),
      source_name = "MFG"
    ),
    TBILI = list(
      description = "Total serum bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The source used mg/dL with reference 0.41 mg/dL, the cohort median",
        "(Supplementary Table S1: median 0.41, range 0.12-2.11 mg/dL). The",
        "model converts the canonical umol/L column back to mg/dL",
        "(TBILI / 17.1) before applying (BIL/0.41)^-0.266 to SN-38 clearance.",
        "The control stream sets the effect to 1 when bilirubin is missing",
        "(BIL = -99); impute the median to reproduce that."
      ),
      source_name = "BIL"
    ),
    CRCL = list(
      description = "Creatinine clearance, absolute (not normalised to body surface area)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect (CrCL/85.04)^0.25 on SN-38 clearance; 85.04 mL/min is",
        "the cohort median (Supplementary Table S1: median 85, range",
        "27-177 mL/min). The source reports mL/min without a /1.73 m^2",
        "normalisation and does not name the estimating equation. The",
        "control stream sets the effect to 1 when CRCL is missing (-99);",
        "impute 85.04 mL/min to reproduce that."
      ),
      source_name = "CRCL"
    )
  )

  covariatesDataExcluded <- list(
    UGT1A1_STAR28_HOM = list(
      description = "UGT1A1*28 homozygous 7/7 genotype indicator (1 = 7/7, 0 = other or unknown)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Tested on SN-38 clearance in the stepwise covariate search and not",
        "significant by the likelihood ratio test (Results; Figure 2). 6.1% of",
        "the cohort was 7/7 homozygous."
      )
    ),
    CONMED_FLUOROURACIL = list(
      description = "5-FU/leucovorin co-administration indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Concomitant therapy (5-FU/LV and oxaliplatin) was screened (Methods,",
        "Covariate model building); only oxaliplatin was retained. 42.7% of",
        "patients received 5-FU/LV (Supplementary Table S1)."
      )
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (Methods, Covariate model building) and not retained. Median 62 years (range 28-87)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened as a body-size covariate alongside height, BMI and BSA; only BSA was retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = "Screened (liver function tests) and not retained. Median 24 IU/L (range 4-202)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = "Screened (liver function tests) and not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Screened (liver function tests) and not retained. Median 4 g/dL (range 2.1-5.1)."
    )
  )

  compartmentData <- list(
    central = list(
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
    ),
    peripheral1 = list(
      analyte = "irinotecan",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "irinotecan",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 440L,
    n_studies = 7L,
    age_median = "62 years",
    age_range = "28-87 years",
    sex_female_pct = 48.6,
    race_ethnicity = c(Asian = 35.2, `non-Asian` = 64.8),
    disease_state = paste(
      "Advanced solid tumours: metastatic pancreatic ductal adenocarcinoma",
      "(NAPOLI-1, N = 260; NCT02551991, N = 56), metastatic colorectal cancer",
      "(PIST-CRC-01, N = 18) and various tumour types (four phase I/II",
      "studies, N = 106)"
    ),
    dose_range = paste(
      "35-156 mg/m^2 irinotecan free base (40-180 mg/m^2 irinotecan",
      "hydrochloride trihydrate) as a 90-minute intravenous infusion every 2",
      "or 3 weeks, alone, with 5-FU/LV, or with 5-FU/LV plus oxaliplatin"
    ),
    bsa_median = "1.71 m^2 (range 1.29-2.48)",
    hepatic_function = "Bilirubin median 0.41 mg/dL (range 0.12-2.11); ALT median 24 IU/L (4-202); albumin median 4 g/dL (2.1-5.1)",
    renal_function = "Creatinine clearance median 85 mL/min (range 27-177)",
    co_medication = "5-FU/LV 42.7%; oxaliplatin (with 5-FU/LV) 12.7%",
    notes = paste(
      "Baseline covariates from Supplementary Table S1; study designs from",
      "Table 1. 1887 total-irinotecan and 1827 SN-38 concentrations above the",
      "lower limit of quantification were modelled (M1 method: 23% and 25% of",
      "samples were below it and excluded). UGT1A1*28 7/7 homozygous: 6.1%.",
      "Liposomal irinotecan manufactured at the previous site: 18.6%.",
      "Doses in this model are irinotecan free base, the source's convention."
    )
  )

  ini({
    # --- Total irinotecan, structural (Table 2 / Supplementary Table S3)
    lcl <- log(17.9)
    label("Irinotecan total clearance (L/week)") # Table 2, 'Irinotecan total clearance, L/week' = 17.9 (THETA(1))
    lvc <- log(4.09)
    label("Irinotecan central volume (L)") # Table 2, 'Irinotecan central volume, L' = 4.09 (THETA(2))
    lq <- log(1.35)
    label("Irinotecan intercompartmental clearance (L/week)") # Table 2, 'Irinotecan inter-compartmental clearance, L/week' = 1.35 (THETA(3))
    lvp <- log(0.421)
    label("Irinotecan peripheral volume (L)") # Table 2, 'Irinotecan peripheral volume, L' = 0.421 (THETA(4))

    # --- SN-38 formation: FR1 and FR2 are the direct- and delayed-pathway
    #     ratios; the control stream normalises them as FM = FR / (1 + FR1 + FR2)
    lfr1 <- log(0.152)
    label("Direct SN-38 formation pathway ratio FR1 (unitless)") # Table 2, 'Fraction of direct irinotecan total rate of elimination' = 0.152 (THETA(5), TVFR1)
    lfr2 <- log(0.629)
    label("Delayed (transit) SN-38 formation pathway ratio FR2 (unitless)") # Table 2, 'Fraction of delayed irinotecan total rate of elimination' = 0.629 (THETA(6), TVFR2)
    lktr <- log(2)
    label("Rate constant of SN-38 formation out of the transit compartment (1/week)") # Table 2, 'Rate of transformation after delay, 1/week' = 2 (THETA(7), KFM)
    lcl_sn38 <- log(19800)
    label("SN-38 total clearance (L/week)") # Table 2, 'SN-38 total clearance, L/week' = 19,800 (THETA(8))

    # --- Covariate effects. The control stream codes every categorical effect
    #     as a factor (1 + THETA); Table 2 prints that factor, so each value
    #     below is the printed factor minus 1.
    e_race_asian_cl <- 0.204
    label("Fractional change in irinotecan clearance for Asian patients (unitless)") # Table 2, CL 'Asian race' = 1.204 (THETA(15))
    e_form_naliri_prevsite_cl <- 0.515
    label("Fractional change in irinotecan clearance for previous-site drug product (unitless)") # Table 2, CL 'Manufacturing site' = 1.515 (THETA(16))
    e_sexf_cl <- -0.201
    label("Fractional change in irinotecan clearance for women (unitless)") # Table 2, CL 'Gender' = 0.799 (THETA(17))
    e_conmed_oxaliplatin_cl <- 0.339
    label("Fractional change in irinotecan clearance with oxaliplatin (unitless)") # Table 2, CL 'Oxaliplatin administration' = 1.339 (THETA(18))
    e_bsa_vc <- 0.573
    label("Power exponent of BSA/1.71 on irinotecan central volume (unitless)") # Table 2, V 'Body surface area' = (BSA/1.71)^0.573 (THETA(20))
    e_form_naliri_prevsite_vc <- -0.128
    label("Fractional change in irinotecan central volume for previous-site drug product (unitless)") # Table 2, V 'Manufacturing site' = 0.872 (THETA(21))
    e_sexf_vc <- -0.114
    label("Fractional change in irinotecan central volume for women (unitless)") # Table 2, V 'Gender' = 0.886 (THETA(22))
    e_form_naliri_prevsite_fr2 <- 0.376
    label("Fractional change in the delayed-pathway ratio FR2 for previous-site drug product (unitless)") # Table 2, delayed fraction 'Manufacturing site' = 1.376 (THETA(19))
    e_tbili_cl_sn38 <- -0.266
    label("Power exponent of bilirubin/0.41 mg/dL on SN-38 clearance (unitless)") # Table 2, SN-38 CL 'Bilirubin' = (BIL/0.41)^-0.266 (THETA(11))
    e_crcl_cl_sn38 <- 0.25
    label("Power exponent of CrCL/85.04 mL/min on SN-38 clearance (unitless)") # Table 2, SN-38 CL 'Creatinine clearance' = (CrCL/85.04)^0.25 (THETA(12))
    e_sexf_cl_sn38 <- -0.198
    label("Fractional change in SN-38 clearance for women (unitless)") # Table 2, SN-38 CL 'Gender' = 0.802 (THETA(13))
    e_conmed_oxaliplatin_cl_sn38 <- -0.344
    label("Fractional change in SN-38 clearance with oxaliplatin (unitless)") # Table 2, SN-38 CL 'Oxaliplatin administration' = 0.656 (THETA(14))

    # --- Inter-individual variability (variances; Table 2 'IIV (% CV)' column,
    #     CV = sqrt(exp(omega^2) - 1) reproduces every printed %CV). Block on
    #     CL, V and FR1 per the three printed covariances.
    etalcl + etalvc + etalfr1 ~ c(
      0.545,
      0.117, 0.066,
      -0.558, -0.103, 0.928
    ) # Table 2: omega^2 CL 0.545, V 0.066, FR1 0.928; cov CL-V 0.117, CL-FR1 -0.558, V-FR1 -0.103
    etalfr2 ~ 0.188 # Table 2, delayed-fraction IIV 0.188 (45.4% CV)
    etalktr ~ 0.135 # Table 2, transformation-rate IIV 0.135 (38% CV)
    etalcl_sn38 ~ 0.126 # Table 2, SN-38 clearance IIV 0.126 (36.6% CV)

    # --- Residual error (control stream $ERROR: Y = IPRED + THETA*IPRED*EPS)
    propSd <- 0.243
    label("Proportional residual error, total irinotecan (fraction)") # Table 2, 'Proportional error on irinotecan' = 0.243 (THETA(9))
    propSd_sn38 <- 0.291
    label("Proportional residual error, SN-38 (fraction)") # Table 2, 'Proportional error on SN-38' = 0.291 (THETA(10))
  })

  model({
    # 1. Derived covariate terms. Bilirubin is supplied in umol/L and converted
    #    to the source's mg/dL before centring on the 0.41 mg/dL median.
    tbili_mgdl <- TBILI / 17.1

    # 2. Individual parameters (control stream $PK)
    cl <- exp(lcl + etalcl) *
      (1 + e_race_asian_cl * RACE_ASIAN) *
      (1 + e_form_naliri_prevsite_cl * FORM_NALIRI_PREVSITE) *
      (1 + e_sexf_cl * SEXF) *
      (1 + e_conmed_oxaliplatin_cl * CONMED_OXALIPLATIN)
    vc <- exp(lvc + etalvc) *
      (BSA / 1.71)^e_bsa_vc *
      (1 + e_form_naliri_prevsite_vc * FORM_NALIRI_PREVSITE) *
      (1 + e_sexf_vc * SEXF)
    q <- exp(lq)
    vp <- exp(lvp)
    fr1 <- exp(lfr1 + etalfr1)
    fr2 <- exp(lfr2 + etalfr2) *
      (1 + e_form_naliri_prevsite_fr2 * FORM_NALIRI_PREVSITE)
    ktr <- exp(lktr + etalktr)
    cl_sn38 <- exp(lcl_sn38 + etalcl_sn38) *
      (tbili_mgdl / 0.41)^e_tbili_cl_sn38 *
      (CRCL / 85.04)^e_crcl_cl_sn38 *
      (1 + e_sexf_cl_sn38 * SEXF) *
      (1 + e_conmed_oxaliplatin_cl_sn38 * CONMED_OXALIPLATIN)

    # SN-38 shares the irinotecan central volume ('VCM = VCP'); the Methods
    # state this was needed to prevent identifiability issues.
    vc_sn38 <- vc

    # 3. Micro-constants. fm1 / fm2 are the fractions of irinotecan total
    #    elimination forming SN-38 directly and via the transit compartment
    #    (Figure 1: FM1 about 9%, FM2 about 35% at the typical values).
    kel <- cl / vc
    fm1 <- fr1 / (1 + fr1 + fr2)
    fm2 <- fr2 / (1 + fr1 + fr2)
    k12 <- q / vc
    k21 <- q / vp
    kel_sn38 <- cl_sn38 / vc_sn38

    # Irinotecan is converted to SN-38 one mole for one mole. The source
    # modelled both analytes in molar units, so the formation flux (mg of
    # irinotecan) is multiplied by the molecular-weight ratio to give mg of
    # SN-38. Irinotecan free base 586.678 g/mol is from the Methods; SN-38
    # (C22H20N2O5) 392.41 g/mol is not printed in the source and was computed
    # from standard atomic weights.
    mw_ratio_sn38 <- 392.41 / 586.678

    # 4. ODE system (control stream compartments PARENT, MET1, PERIP, TRANSIT1)
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(central_sn38) <- mw_ratio_sn38 * (fm1 * kel * central + ktr * transit1) -
      kel_sn38 * central_sn38
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(transit1) <- fm2 * kel * central - ktr * transit1

    # 6. Observations. Irinotecan in ug/mL (= mg/L); SN-38 in ng/mL, the units
    #    of the source's figures (control-stream scale S2 = VCP / 1000).
    Cc <- central / vc
    Cc_sn38 <- 1000 * central_sn38 / vc_sn38

    Cc ~ prop(propSd)
    Cc_sn38 ~ prop(propSd_sn38)
  })
}
