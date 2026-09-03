Cao_2026_lecanemab <- function() {
  description <- paste(
    "QSP. Neuro-Dynamic quantitative systems pharmacology model of",
    "Alzheimer's disease pathophysiology and anti-amyloid treatment",
    "effects, with lecanemab as the reference drug. Eleven ODEs across",
    "three interlinked modules. Amyloid module: a four-stage reversible",
    "A-beta aggregation cascade (monomer A -> oligomer O -> protofibril F",
    "-> plaque P) in which F and P cooperatively accelerate every forward",
    "aggregation step. Tau module: five neuron states (N healthy, I0",
    "minimal misfolded tau, I1 misfolded tau oligomers, I2 tau tangles,",
    "DN irreversible neurodegeneration) driven by spontaneous misfolding,",
    "tau-seed (TS) infection of healthy neurons, and amyloid promotion,",
    "plus a phosphorylated-tau pool (TP); transitions are reversible",
    "except into DN. Cognitive module: CDR-SB and ADAS-Cog are power",
    "functions of a weighted composite of I1, I2 and DN neuron fractions",
    "plus a direct amyloid-plaque effect. Six observables are fitted",
    "(amyloid PET Centiloid, medial-temporal tau PET SUVR, plasma",
    "A-beta 42/40 ratio, plasma p-tau181, CDR-SB, ADAS-Cog) and two are",
    "predicted but not fitted (CSF A-beta 42, CSF p-tau181). Lecanemab",
    "acts as additive first-order clearance on protofibrils and plaques,",
    "driven by the steady-state average serum concentration covariate",
    "CSS_LEC; the estimated protofibril effect is about 2.3-fold the",
    "plaque effect at the cohort median age. Model time is YEARS SINCE AD",
    "PATHOLOGICAL ONSET (not study time): the analysis dataset was aligned",
    "to a subject-specific onset using the GRACE method. Fitted by SAEM in",
    "Monolix 2021R1 to 4056 subjects from lecanemab Study 201, Study 301",
    "(Clarity AD) Core and OLE, and the ADNI natural-history cohort.",
    "Of 74 parameters, 35 are fixed, 32 are estimated with IIV and 7 with",
    "fixed effects only.",
    sep = " "
  )
  reference <- paste(
    "Cao Y, Willis BA, Horie K, Wildsmith KR, Koyama A, Sachdev P,",
    "Penner N, Charil A, Irizarry M, Reyderman L (2026).",
    "Neuro-Dynamic Quantitative Systems Pharmacology (QSP) model",
    "describing Alzheimer's disease pathophysiology and treatment",
    "effects.",
    "npj Systems Biology and Applications.",
    "doi:10.1038/s41540-026-00677-4.",
    sep = " "
  )
  vignette <- "Cao_2026_lecanemab"

  # Every state in this model is a paper-mechanistic species or neuron
  # population; none maps onto a canonical PK compartment name. Names follow
  # the paper's own symbols (Eq 1-2 and the "Full ODE system" section of
  # Methods) so the source trace is one-to-one.
  paper_specific_compartments <- c(
    "A", "O", "F", "P", "N", "I0", "I1", "I2", "DN", "TS", "TP"
  )

  units <- list(
    time          = "year",
    dosing        = paste(
      "(this model has no dosing events; anti-amyloid drug exposure enters",
      "through the time-varying covariate CSS_LEC, the steady-state average",
      "serum concentration in ug/mL. See covariateData.)",
      sep = " "
    ),
    concentration = paste(
      "A, O, F, P, TS and TP are concentrations in pg/mL; N, I0, I1, I2 and",
      "DN are neuron densities in neurons/mL. Observables carry their own",
      "scales: CENTILOID (Centiloid), TAUSUVR and ABSUVR (SUVR), PAB",
      "(dimensionless ratio), PPTAU and CSFPTAU (pg/mL), CSFAB (pg/mL),",
      "CDRSB (0-18 points), ADASCOG (0-90 points).",
      sep = " "
    )
  )

  covariateData <- list(
    AGE = list(
      description        = paste(
        "Baseline age in years. Enters three log-scale multiplicative",
        "effects, all UNCENTRED (the raw age in years multiplies the",
        "coefficient): on the plasma A-beta 42/40 ratio floor PABFLOOR",
        "(Table S14, beta_PABFLOOR_AGE = -0.0019), on the ADAS-Cog shape",
        "parameter (beta_shape_adas_AGE = -0.021), and on the lecanemab",
        "plaque drug effect DELECAP (beta_DELECAP_AGE = 0.041). The paper",
        "does not print the covariate equation, but the uncentred form is",
        "confirmed by three independent back-calculations against the",
        "individual-parameter medians in Supplementary Table S15 at the",
        "cohort mean age 71.3 years: PABFLOOR 0.091 * exp(-0.0019 * 71.3)",
        "= 0.0795 vs S15 median 0.079; shape_adas 6.23 * exp(-0.021 *",
        "71.3) = 1.40 vs S15 median 1.3; DELECAP 0.0015 * exp(0.041 *",
        "71.3) = 0.0277 vs S15 median 0.03 (and the main-text statement",
        "that the median plaque effect is 0.03). A centred form reproduces",
        "none of the three.",
        sep = " "
      ),
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "UNCENTRED - use the raw age, not age minus a reference. Cohort",
        "means (Supplementary Table A1): ADNI 73.4 (SD 7.7), Study 201",
        "71.3 (8.2), Study 301 71.3 (7.8); overall range 50-91 years.",
        sep = " "
      ),
      source_name        = "AGE"
    ),
    APOE4_HET = list(
      description        = paste(
        "1 = subject carries exactly one APOE-epsilon4 allele",
        "(heterozygous); 0 = otherwise (non-carrier or homozygote).",
        "Paired with APOE4_HOM to encode the three-level APOE4 genotype",
        "with non-carrier as the reference. Multiplicative log-scale",
        "effect on the A-beta monomer aggregation rate kAggAO",
        "(Supplementary Table S14, beta_kAggAO_APOEGEN_Hetero = 0.15).",
        sep = " "
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = APOE-epsilon4 non-carrier (the reference genotype in the paper's covariate model)",
      notes              = paste(
        "APOE4_HET and APOE4_HOM are mutually exclusive; a non-carrier has",
        "both set to 0. Cohort frequencies (Table A1): heterozygous 41.2%",
        "(ADNI), 54.9% (Study 201), 53.3% (Study 301). The paper reports",
        "distinct, non-additive heterozygote and homozygote coefficients",
        "(0.15 vs 0.18), so the additive allele-count encoding",
        "APOE4_COUNT cannot reproduce this model.",
        sep = " "
      ),
      source_name        = "APOEGEN (Monolix categorical covariate, level 'Hetero')"
    ),
    APOE4_HOM = list(
      description        = paste(
        "1 = subject carries two APOE-epsilon4 alleles (homozygous);",
        "0 = otherwise (non-carrier or heterozygote). Multiplicative",
        "log-scale effect on kAggAO (Supplementary Table S14,",
        "beta_kAggAO_APOEGEN_Homo = 0.18).",
        sep = " "
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = APOE-epsilon4 non-carrier",
      notes              = paste(
        "Cohort frequencies (Table A1): homozygous 13.4% (ADNI), 16.5%",
        "(Study 201), 15.3% (Study 301). See APOE4_HET.",
        sep = " "
      ),
      source_name        = "APOEGEN (Monolix categorical covariate, level 'Homo')"
    ),
    APOE4_CARRIER = list(
      description        = paste(
        "1 = subject carries at least one APOE-epsilon4 allele",
        "(heterozygous or homozygous); 0 = non-carrier. The paper uses",
        "this binary collapse - separately from the three-level genotype",
        "used on kAggAO - for the effects on the CDR-SB and ADAS-Cog",
        "scaling factors (Supplementary Table S14,",
        "beta_alpha_tAPOECarrier = -0.21 on alphaCdr and",
        "beta_alpha_adas_tAPOECarrier = -0.097 on alphaAdas).",
        sep = " "
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = APOE-epsilon4 non-carrier",
      notes              = paste(
        "Derivable as APOE4_CARRIER = max(APOE4_HET, APOE4_HOM); carried",
        "as its own column because the paper fits it as a distinct",
        "covariate from the three-level genotype. Both codings are",
        "required to reproduce the published model.",
        sep = " "
      ),
      source_name        = "tAPOECarrier (Monolix transformed binary carrier covariate)"
    ),
    RACE_ASIAN = list(
      description        = paste(
        "1 = Asian race (including the Japanese, Korean, Chinese and",
        "'Other Asian' categories reported separately in Supplementary",
        "Table A1); 0 = otherwise. Multiplicative log-scale effects on",
        "the plasma A-beta 42/40 ratio floor PABFLOOR",
        "(beta_PABFLOOR_RACE_Asian = 0.038) and on the ADAS-Cog scaling",
        "factor alphaAdas (beta_alpha_adas_RACE_Asian = 0.14);",
        "Supplementary Table S14.",
        sep = " "
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = non-Asian (predominantly White: 92.9% ADNI, 90.6% Study 201, 76.9% Study 301)",
      notes              = paste(
        "Study 301 enrolled a substantial East Asian cohort (Japanese",
        "8.5%, Korean 7.2%, Chinese 0.7%, Other Asian 0.3%, Asian 0.2% =",
        "about 17% in total); Study 201 about 6.5%; ADNI 1.8%.",
        sep = " "
      ),
      source_name        = "RACE (Monolix categorical covariate, level 'Asian')"
    ),
    DIS_AD_LMCI = list(
      description        = paste(
        "1 = baseline diagnosis of late mild cognitive impairment (LMCI,",
        "an ADNI diagnostic category); 0 = otherwise. One of four",
        "indicators encoding the pooled baseline-diagnosis covariate with",
        "EMCI (early mild cognitive impairment) as the reference level.",
        "Multiplicative log-scale effect on alphaAdas",
        "(beta_alpha_adas_bDIAG_LMCI = 0.13; Supplementary Table S14).",
        sep = " "
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0; the reference level of the baseline-diagnosis covariate is EMCI (all four DIS_AD_* indicators set to 0)",
      notes              = paste(
        "The pooled dataset carries study-specific diagnostic vocabularies",
        "(Supplementary Table A1): ADNI recorded EMCI / LMCI / Mild AD",
        "dementia, while Studies 201 and 301 recorded MCI / Mild Dementia",
        "due to AD. The five levels appearing in Supplementary Table S14",
        "are therefore EMCI (reference), LMCI, MCI, Mild_AD and AD. See",
        "the vignette Errata for the residual ambiguity about which",
        "levels map to which source category.",
        sep = " "
      ),
      source_name        = "bDIAG (Monolix categorical covariate, level 'LMCI')"
    ),
    DIS_AD_MCI = list(
      description        = paste(
        "1 = baseline diagnosis of mild cognitive impairment (MCI, the",
        "Study 201 / Study 301 category); 0 = otherwise. Multiplicative",
        "log-scale effects on the A-beta monomer aggregation rate kAggAO",
        "(beta_kAggAO_DIAG_MCI = 0.12) and on alphaAdas",
        "(beta_alpha_adas_bDIAG_MCI = 0.14); Supplementary Table S14.",
        sep = " "
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0; reference level EMCI",
      notes              = paste(
        "Cohort frequencies (Table A1): MCI 64.1% (Study 201), 61.9%",
        "(Study 301). See DIS_AD_LMCI.",
        sep = " "
      ),
      source_name        = "DIAG / bDIAG (Monolix categorical covariate, level 'MCI')"
    ),
    DIS_AD_MILD = list(
      description        = paste(
        "1 = baseline diagnosis of mild dementia due to Alzheimer's",
        "disease (the Study 201 / Study 301 category 'Mild Dementia due",
        "to AD'); 0 = otherwise. Multiplicative log-scale effects on",
        "kAggAO (beta_kAggAO_DIAG_Mild_AD = 0.032), on alphaCdr",
        "(beta_alpha_bDIAG_Mild_AD = 0.18) and on alphaAdas",
        "(beta_alpha_adas_bDIAG_Mild_AD = 0.16); Supplementary Table S14.",
        sep = " "
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0; reference level EMCI",
      notes              = paste(
        "Cohort frequencies (Table A1): mild dementia due to AD 35.9%",
        "(Study 201), 38.1% (Study 301). See DIS_AD_LMCI.",
        sep = " "
      ),
      source_name        = "DIAG / bDIAG (Monolix categorical covariate, level 'Mild_AD')"
    ),
    DIS_AD_DEMENTIA = list(
      description        = paste(
        "1 = baseline diagnosis of Alzheimer's dementia recorded under the",
        "ADNI vocabulary ('Mild AD dementia', the level printed as 'AD' in",
        "Supplementary Table S14); 0 = otherwise. Multiplicative log-scale",
        "effects on alphaCdr (beta_alpha_bDIAG_AD = 0.39) and on alphaAdas",
        "(beta_alpha_adas_bDIAG_AD = 0.39).",
        sep = " "
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0; reference level EMCI",
      notes              = paste(
        "Cohort frequency (Table A1): Mild AD dementia 27.6% of the ADNI",
        "cohort. Kept distinct from DIS_AD_MILD because Supplementary",
        "Table S14 estimates separate coefficients for the 'AD' and",
        "'Mild_AD' levels (0.39 vs 0.18 on alphaCdr). See DIS_AD_LMCI.",
        sep = " "
      ),
      source_name        = "bDIAG (Monolix categorical covariate, level 'AD')"
    ),
    CSS_LEC = list(
      description        = paste(
        "Steady-state average serum concentration of lecanemab, generated",
        "externally by the lecanemab population PK model (Hayato et al.,",
        "references 22 and 24 of the paper) and fed into the QSP model as",
        "the drug-exposure driver. It enters two additive clearance-",
        "enhancement terms on the amyloid states: Eff_LECAF = 1 + DELECAF",
        "* CSS_LEC on the protofibril clearance cF, and Eff_LECAP = 1 +",
        "DELECAP * CSS_LEC on the plaque clearance cP (Methods, 'Full ODE",
        "system').",
        sep = " "
      ),
      units              = "ug/mL (micrograms per millilitre; DELECAF and DELECAP carry the reciprocal units mL/ug per Supplementary Table S13)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TIME-VARYING and equal to 0 whenever no anti-amyloid antibody is",
        "on board. The paper writes the drug effect as an if / elseif /",
        "else block over two treatment windows (Core and OLE) that sets",
        "Eff_LECAF = Eff_LECAP = 1 outside treatment; setting CSS_LEC = 0",
        "off treatment is algebraically identical and covers any number of",
        "treatment and washout windows. Reference value for the approved",
        "10 mg/kg IV Q2W regimen: the Vpop301 Css,av distribution",
        "(Supplementary Figure E21) is log-normal with a mode near",
        "115-130 ug/mL and a median near 135 ug/mL, spanning about 30-370",
        "ug/mL; the figure is the only place the value is reported, so",
        "135 ug/mL is a digitised median. Member of the CSS_<drug>",
        "canonical family (compare CSS_RBV, CSS_ALBIFN, CSS_DFO). The",
        "same column drives donanemab, aducanumab and gantenerumab in the",
        "paper's cross-antibody simulations via drug-specific exposure",
        "multipliers, which the paper does not tabulate.",
        sep = " "
      ),
      source_name        = "CSSAV01 / CSSAV02 (per-window steady-state average concentration in the paper's Methods pseudocode)"
    )
  )

  compartmentData <- list(
    A  = list(analyte = "Amyloid-beta 42 monomer",       units = "pg/mL", specimen = "brain ISF", verified = TRUE),
    O  = list(analyte = "Amyloid-beta oligomer",         units = "pg/mL", specimen = "brain ISF", verified = TRUE),
    F  = list(analyte = "Amyloid-beta protofibril",      units = "pg/mL", specimen = "brain ISF", verified = TRUE),
    P  = list(analyte = "Amyloid-beta plaque",           units = "pg/mL", specimen = "brain ISF", verified = TRUE),
    N  = list(analyte = "Healthy neurons (no tau pathology)",              units = "neurons/mL", specimen = "tissue", verified = TRUE),
    I0 = list(analyte = "Neurons with minimal misfolded tau",              units = "neurons/mL", specimen = "tissue", verified = TRUE),
    I1 = list(analyte = "Neurons with misfolded tau oligomers",            units = "neurons/mL", specimen = "tissue", verified = TRUE),
    I2 = list(analyte = "Neurons with tau neurofibrillary tangles",        units = "neurons/mL", specimen = "tissue", verified = TRUE),
    DN = list(analyte = "Degenerated (irreversibly lost) neurons",         units = "neurons/mL", specimen = "tissue", verified = TRUE),
    TS = list(analyte = "Tau seed",                     units = "pg/mL", specimen = "brain ISF", verified = TRUE),
    TP = list(analyte = "Phosphorylated tau 181",       units = "pg/mL", specimen = "CSF",       verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 4056L,
    n_studies      = 3L,
    age_range      = "50-91 years (ADNI 54.4-91.4, Study 201 50-90, Study 301 50-90; Supplementary Table A1)",
    age_median     = "72 years (Studies 201 and 301); 73.9 years (ADNI)",
    weight_range   = "(body weight was not a covariate in this model and is not tabulated in the source)",
    sex_female_pct = 48.0,
    race_ethnicity = c(White = 83.7, Asian = 9.5, Black = 3.0, Other = 3.8),
    disease_state  = paste(
      "Amyloid-positive early symptomatic Alzheimer's disease: mild",
      "cognitive impairment or mild dementia due to AD. Baseline CDR-SB",
      "mean (SD) 2.3 (1.8) in ADNI, 2.9 (1.4) in Study 201, 3.2 (1.3) in",
      "Study 301. Baseline diagnosis: ADNI EMCI 27.0% / LMCI 45.3% / mild",
      "AD dementia 27.6%; Study 201 MCI 64.1% / mild dementia due to AD",
      "35.9%; Study 301 MCI 61.9% / mild dementia due to AD 38.1%. ADNI",
      "subjects with a baseline diagnosis of cognitively normal or",
      "significant memory concern were excluded. APOE-epsilon4",
      "non-carrier / heterozygous / homozygous: 45.5 / 41.2 / 13.4%",
      "(ADNI), 28.6 / 54.9 / 16.5% (Study 201), 31.4 / 53.3 / 15.3%",
      "(Study 301).",
      sep = " "
    ),
    dose_range     = paste(
      "Lecanemab 2.5, 5 and 10 mg/kg IV every two weeks or monthly",
      "(Study 201 Core), 10 mg/kg IV every two weeks (Study 301 Clarity",
      "AD Core and both open-label extensions). The model itself takes",
      "exposure as the CSS_LEC covariate rather than dose amounts.",
      sep = " "
    ),
    regions        = "North America, Europe, Japan, China, South Korea and other Asia-Pacific sites (Studies 201 and 301); North America (ADNI)",
    notes          = paste(
      "Model time is YEARS SINCE ESTIMATED AD PATHOLOGICAL ONSET, not",
      "study time. Individual onset times were estimated with the GRACE",
      "method (Donohue et al.) from ADNI plus lecanemab data, then shifted",
      "by +20 years so all times are positive, with a further +1 year",
      "adjustment applied to the lecanemab-treated group to remove a",
      "systematic offset caused by using baseline-only data for treated",
      "subjects (Methods, 'GRACE'). At t = 0 the model represents a state",
      "in which A-beta monomer production has begun but no oligomer,",
      "protofibril or plaque has accumulated, no tau pathology has",
      "initiated, amyloid and tau PET are negative and cognition is",
      "normal. Follow-up: about 4 years (Study 301), 10 years (Study 201",
      "including the Gap period) and 15 years (ADNI). Endpoint counts",
      "include 7061 amyloid PET records. Tau PET (MK-6240, medial",
      "temporal region) came only from a Study 301 substudy (N = 355);",
      "ADNI tau PET was Braak-staged and could not be converted, and ADNI",
      "plasma biomarkers were excluded for assay incomparability. CSF",
      "A-beta 42 and CSF p-tau181 were NOT used to train the model and",
      "are pure predictions. Model qualification used goodness-of-fit and",
      "VPCs; all estimated parameters had RSE <= 35% and shrinkage <=",
      "5.34%, with a final Fisher-information condition number of 3492.1.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # AMYLOID PATHWAY MODULE
    # Source: Supplementary Table S3 ("Final QSP Model Amyloid Pathway
    # Module Parameters") unless noted. Fixed values are cross-checked
    # against Supplementary Table S2 ("Fixed Parameters in the Amyloid
    # Pathway Module") and Supplementary Note 2. All rates are per year.
    # Monolix estimated these on the log scale with log-normal IIV, so the
    # log transform below reproduces the published distribution exactly.
    # =====================================================================
    lA0       <- fixed(log(1700))   ; label("Initial CSF A-beta 42 monomer concentration A0 (log-pg/mL)")                  # Table S3 A0 = 1700 (Fixed); Note 2: matches CSF Ab(1-42) 1678 +/- 436 pg/mL in healthy controls (Andreasen 1999)
    lpNA      <- log(0.054)         ; label("A-beta monomer production rate per neuron pNA (log-pg/neuron/year)")          # Table S3 pNA = 0.054 (RSE 0.482%)
    lkAggAO   <- log(20.11)         ; label("A-beta monomer -> oligomer aggregation rate kAggAO (log-1/year)")             # Table S3 kAggAO = 20.11 (RSE 1.00%)
    lkDeAggOA <- log(693.62)        ; label("A-beta oligomer -> monomer deaggregation rate kDeAggOA (log-1/year)")         # Table S3 kDeAggOA = 693.62 (RSE 0.124%)
    lkAggOF   <- fixed(log(283.98)) ; label("A-beta oligomer -> protofibril aggregation rate kAggOF (log-1/year)")         # Table S3 kAggOF = 283.98 (Fixed); Note 2 fixes it at 284/year from Michaels 2020 (9E-6 /s)
    lkDeAggFO <- fixed(log(10))     ; label("A-beta protofibril -> oligomer deaggregation rate kDeAggFO (log-1/year)")     # Table S3 kDeAggFO = 10 (Fixed); Table S2 "Assumed"
    lkAggFP   <- log(0.34)          ; label("A-beta protofibril -> plaque aggregation rate kAggFP (log-1/year)")           # Table S3 kAggFP = 0.34 (RSE 4.50%)
    lkDeAggPF <- fixed(log(0.025))  ; label("A-beta plaque -> protofibril deaggregation rate kDeAggPF (log-1/year)")       # Table S3 kDeAggPF = 0.025 (Fixed fixed-effect, but IIV estimated); Table S2 "Assumed"
    lsA       <- log(1.4e-5)        ; label("Cooperative effect of F and P on monomer aggregation sA (log-mL/pg)")         # Table S3 sA = 1.4E-05 (RSE 0.0501%)
    lsO       <- log(3.4e-4)        ; label("Cooperative effect of F and P on oligomer aggregation sO (log-mL/pg)")        # Table S3 sO = 3.4E-04 (RSE 0.0262%)
    lsF       <- fixed(log(6.9e-8)) ; label("Cooperative effect of F and P on protofibril aggregation sF (log-mL/pg)")     # Table S3 sF = 6.9E-08 (Fixed fixed-effect, but IIV estimated); Table S2 "Calibrated"
    lsP       <- log(5.6e-8)        ; label("Scaling factor of plaque to amyloid PET SUVR sP (log-SUVR*mL/pg)")            # Table S3 sP = 5.6E-08 (RSE 0.0636%)
    lcA       <- fixed(log(1500))   ; label("A-beta monomer clearance rate cA (log-1/year)")                               # Table S3 cA = 1500 (Fixed); Note 2: 3-h monomer half-life implies cA + kAggAO ~ 2024/year (May 2015, Ovod 2017)
    lcO       <- fixed(log(0.5))    ; label("A-beta oligomer clearance rate cO (log-1/year)")                              # Table S3 cO = 0.5 (Fixed); Note 2 derives a near-zero value from Michaels 2020 and fixes it at 0.5
    lcF       <- log(0.12)          ; label("A-beta protofibril clearance rate cF (log-1/year)")                           # Table S3 cF = 0.12 (RSE 0.0941%)
    lcP       <- log(0.26)          ; label("A-beta plaque clearance rate cP (log-1/year)")                                # Table S3 cP = 0.26 (RSE 0.0316%)
    lsC2PAb   <- log(9.9e-6)        ; label("Scaling from CSF A-beta 42 to plasma A-beta 42/40 ratio sC2PAb (log-mL/pg)")  # Table S3 sC2PAb = 9.9E-06 (RSE 2.10%)
    lPABFLOOR <- log(0.091)         ; label("Plasma A-beta 42/40 ratio floor PABFLOOR (log-ratio)")                        # Table S3 PABFLOOR = 0.091 (RSE 1.74%)
    lABSUVR0  <- log(1.06)          ; label("Amyloid PET SUVR floor ABSUVR0 (log-SUVR)")                                   # Table S3 ABSUVR0 = 1.06 (RSE 0.260%); Note 2 quotes 1.05, corresponding to 1.5 Centiloid

    # =====================================================================
    # TAU PATHOLOGY PATHWAY MODULE
    # Source: Supplementary Table S6 ("Final QSP Model Tau Pathway Module
    # Parameters"), cross-checked against Supplementary Table S5 (fixed
    # parameters) and Supplementary Note 3. The paper's Greek symbols are
    # transliterated: eta01/eta12/eta2d -> eFT01/eFT12/eFT2d (renamed to
    # avoid colliding with the "eta" prefix reserved for IIV terms), and
    # beta -> betaTS (renamed to avoid colliding with the "beta_" prefix
    # Monolix uses for covariate coefficients).
    # =====================================================================
    ltN         <- fixed(log(1e-5))    ; label("Spontaneous transition rate N -> I0, tN (log-1/year)")                          # Table S6 tN = 1.0E-05 (Fixed); Table S5 "Calibrated"
    lbetaTS     <- fixed(log(10))      ; label("Tau-seed transmission rate beta (log-mL/pg/year)")                              # Table S6 beta = 10 (Fixed); Table S5 "Calibrated"
    lsN1        <- fixed(log(5e-5))    ; label("Acceleration by F and P of spontaneous tau generation sN1 (log-mL/pg)")         # Table S6 sN1 = 5.0E-05 (Fixed); Table S5 "Calibrated"
    lsN2        <- fixed(log(2e-7))    ; label("Acceleration by F and P of tau-seed transmission sN2 (log-mL/pg)")              # Table S6 sN2 = 2.0E-07 (Fixed); Table S5 "Calibrated"
    lrI0        <- fixed(log(82.95))   ; label("Recovery rate I0 -> N, rI0 (log-1/year)")                                       # Table S6 rI0 = 82.95 (Calibrated)
    lrI1        <- fixed(log(0.23))    ; label("Recovery rate I1 -> I0, rI1 (log-1/year)")                                      # Table S6 rI1 = 0.23 (Calibrated)
    lrI2        <- fixed(log(0.00037)) ; label("Recovery rate I2 -> I1, rI2 (log-1/year)")                                      # Table S6 rI2 = 0.00037 (Calibrated); Note 3 quotes rI2 = 0.00037/year
    ltI01       <- log(0.0051)         ; label("Transition rate I0 -> I1, tI01 (log-1/year)")                                   # Table S6 tI01 = 0.0051 (RSE 0.531%)
    ltI12       <- log(0.01)           ; label("Transition rate I1 -> I2, tI12 (log-1/year)")                                   # Table S6 tI12 = 0.01 (RSE 0.0614%); Note 3 quotes tI12 = 0.01/year
    leFT01      <- fixed(log(1.0e-7))  ; label("Effect of F and P on the I0 -> I1 transition, eta01 (log-mL/pg)")               # Table S6 eta01 = 1.00E-07 (Calibrated)
    leFT12      <- fixed(log(1.2e-8))  ; label("Effect of F and P on the I1 -> I2 transition, eta12 (log-mL/pg)")               # Table S6 eta12 = 1.2E-08 (Calibrated)
    leFT2d      <- fixed(log(8.8e-10)) ; label("Effect of F and P on the I2 -> DN transition, eta2d (log-mL/pg)")               # Table S6 eta2d = 8.8E-10 (Calibrated)
    lrho01      <- fixed(log(7.3e-6))  ; label("Effect of tau seeds on the I0 -> I1 transition, rho01 (log-mL/pg)")             # Table S6 rho01 = 7.3E-06 (Calibrated); Note 3 gives 7.26E-06 and rho01*TS = 58% at a baseline TS of 8E+04 pg/mL
    lrho12      <- fixed(log(1e-14))   ; label("Effect of tau seeds on the I1 -> I2 transition, rho12 (log-mL/pg)")             # Note 3: "the estimated tau seeding effect on transition rates from I1->I2 and I2->neurodegeneration stage is negligible; therefore, they were fixed to rho12 = rho2d = 1E-14"
    lrho2d      <- fixed(log(1e-14))   ; label("Effect of tau seeds on the I2 -> DN transition, rho2d (log-mL/pg)")             # Note 3, same sentence as rho12
    ldI2        <- fixed(log(0.012))   ; label("Neurodegeneration rate I2 -> DN, dI2 (log-1/year)")                             # Table S6 dI2 = 0.012 (Calibrated); Note 3: corresponds to a 5.8-year I2 half-life
    lpSI0       <- log(2.0e-6)         ; label("Tau-seed production rate from I0 neurons pSI0 (log-pg/neuron/year)")            # Table S6 pSI0 = 2.0E-06 (RSE 0.279%); Note 3 quotes 2.04E-06
    lpSI1       <- log(8.2e-5)         ; label("Tau-seed production rate from I1 neurons pSI1 (log-pg/neuron/year)")            # Table S6 pSI1 = 8.2E-05 (RSE 0.111%)
    lpSI2       <- fixed(log(1.65))    ; label("Tau-seed production rate from I2 neurons pSI2 (log-pg/neuron/year)")            # Table S6 pSI2 = 1.65 (Calibrated); Note 3 quotes pSI2 = 1.65 pg/neuron/year
    lpPI0       <- log(1.8e-6)         ; label("p-tau181 production rate from I0 neurons pPI0 (log-pg/neuron/year)")            # Table S6 pPI0 = 1.8E-06 (RSE 4.87%)
    lpPI1       <- log(4.9e-5)         ; label("p-tau181 production rate from I1 neurons pPI1 (log-pg/neuron/year)")            # Table S6 pPI1 = 4.9E-05 (RSE 4.11%)
    lpPI2       <- fixed(log(4.0e-5))  ; label("p-tau181 production rate from I2 neurons pPI2 (log-pg/neuron/year)")            # Table S6 pPI2 = 4.0E-05 (Calibrated)
    lpPN        <- fixed(log(1.0e-6))  ; label("p-tau181 production rate from healthy neurons pPN (log-pg/neuron/year)")        # Table S6 pPN = 1.0E-06 (Calibrated)
    lcTS        <- fixed(log(8.6))     ; label("Tau-seed clearance rate cTS (log-1/year)")                                      # Table S6 cTS = 8.6 (Fixed); Table S5 "Calculated from Sato et al. (2018)"
    lcTP        <- fixed(log(8.6))     ; label("CSF p-tau181 clearance rate cTP (log-1/year)")                                  # Table S6 cTP = 8.6 (Fixed); Table S5 "Calculated from Sato et al. (2018)"; Note 3 quotes a 30.7-day tau half-life
    lsTauS      <- log(8.8e-8)         ; label("Acceleration by F and P of tau-seed production sTauS (log-mL/pg)")              # Table S6 sTauS = 8.8E-08 (RSE 0.0821%)
    lsTauP      <- log(8.2e-8)         ; label("Acceleration by F and P of p-tau181 production sTauP (log-mL/pg)")              # Table S6 sTauP = 8.2E-08 (RSE 0.410%)
    lsC2PPtau   <- log(0.078)          ; label("Scaling from CSF to plasma p-tau181 sC2PPtau (log-dimensionless)")              # Table S6 sC2PPtau = 0.078 (RSE 3.88%)
    lPPTAUFLOOR <- log(0.32)           ; label("Plasma p-tau181 floor PPTAUFLOOR (log-pg/mL)")                                  # Table S6 PPTAUFLOOR = 0.32 (RSE 3.25%)
    lTAUSUVR0   <- log(0.69)           ; label("Medial-temporal tau PET SUVR floor TAUSUVR0 (log-SUVR)")                        # Table S6 TAUSUVR0 = 0.69 (RSE 0.232%)
    lsI2TAUPET  <- log(6.0e-7)         ; label("Scaling from I1 and I2 neurons to tau PET SUVR sI2TAUPET (log-SUVR*mL/neuron)") # Table S6 sI2TAUPET = 6.0E-07 (RSE 1.96%)
    lfI1T       <- fixed(log(0.47))    ; label("Relative contribution of I1 neurons to tau PET SUVR fI1T (log-dimensionless)")  # Table S6 fI1T = 0.47 (Calibrated)
    lShapeTau   <- fixed(log(1.5))     ; label("Shape (power) parameter for tau PET SUVR, ShapeTau (log-dimensionless)")        # Table S6 ShapeTau = 1.5 (Calibrated)

    # =====================================================================
    # COGNITIVE FUNCTION MODULE
    # Source: Supplementary Table S9 ("Final QSP Model Cognitive Function
    # Module Parameters"), cross-checked against Supplementary Table S8
    # (fixed parameters) and Supplementary Note 4. N0 comes from
    # Supplementary Note 2.
    # =====================================================================
    lN0        <- fixed(log(6.14e7)) ; label("Initial healthy-neuron density N0 (log-neurons/mL)")                        # Note 2: 86 billion neurons / 1400 mL brain volume (Herculano-Houzel 2009) = 6.14E+07 neurons/mL
    lCDR0      <- log(0.69)          ; label("CDR-SB floor effect CDR0 (log-CDR-SB points)")                              # Table S9 CDR0 = 0.69 (RSE 2.12%)
    lalphaCdr  <- log(3.58)          ; label("Scaling from pathological neurons to CDR-SB, alpha_cdr (log-points*mL/neuron)")   # Table S9 alpha_cdr = 3.58 (RSE 1.97%)
    lehcCdr    <- log(2.4e-8)        ; label("Direct effect of A-beta plaque on CDR-SB, ehc_cdr (log-mL/pg)")             # Table S9 ehc_cdr = 2.4E-08 (RSE 1.47%)
    lshapeCdr  <- log(2.92)          ; label("CDR-SB shape (power) parameter shape_cdr (log-dimensionless)")              # Table S9 shape_cdr = 2.92 (RSE 3.77%)
    ldcdnCdr   <- log(1075.4)        ; label("Contribution of degenerated neurons to CDR-SB, dcdn_cdr (log-mL/neuron)")   # Table S9 dcdn_cdr = 1075.4 (RSE 0.295%)
    ldc2Cdr    <- fixed(log(20))     ; label("Contribution of I2 neurons to CDR-SB, dc2_cdr (log-mL/neuron)")             # Table S9 dc2_cdr = 20 (Fixed); Table S8 "Calibrated and fixed"
    ldc1Cdr    <- fixed(log(2.5))    ; label("Contribution of I1 neurons to CDR-SB, dc1_cdr (log-mL/neuron)")             # Table S9 dc1_cdr = 2.5 (Fixed); Table S8 "Calibrated and fixed"
    lADAS0     <- log(7.54)          ; label("ADAS-Cog floor effect ADAS0 (log-ADAS-Cog points)")                         # Table S9 ADAS0 = 7.54 (RSE 1.30%)
    lalphaAdas <- log(16.9)          ; label("Scaling from pathological neurons to ADAS-Cog, alpha_adas (log-points*mL/neuron)") # Table S9 alpha_adas = 16.9 (RSE 2.78%)
    lehcAdas   <- log(2.0e-8)        ; label("Direct effect of A-beta plaque on ADAS-Cog, ehc_adas (log-mL/pg)")          # Table S9 ehc_adas = 2.0E-08 (RSE 0.594%)
    lshapeAdas <- log(6.23)          ; label("ADAS-Cog shape (power) parameter shape_adas (log-dimensionless)")           # Table S9 shape_adas = 6.23 (RSE 13.3%)
    ldcdnAdas  <- log(1323.24)       ; label("Contribution of degenerated neurons to ADAS-Cog, dcdn_adas (log-mL/neuron)")# Table S9 dcdn_adas = 1323.24 (RSE 0.325%)
    ldc2Adas   <- fixed(log(2.5))    ; label("Contribution of I2 neurons to ADAS-Cog, dc2_adas (log-mL/neuron)")          # Table S9 dc2_adas = 2.5 (Fixed); Table S8 "Calibrated and fixed"
    ldc1Adas   <- fixed(log(1.6))    ; label("Contribution of I1 neurons to ADAS-Cog, dc1_adas (log-mL/neuron)")          # Table S9 dc1_adas = 1.6 (Fixed); Table S8 "Calibrated and fixed"

    # =====================================================================
    # INTER-MODULE CONTRIBUTION FACTORS (protofibril vs plaque potency)
    # Source: Supplementary Table S12 ("Final QSP Model Protofibril/Plaque
    # Impact Parameters"), cross-checked against Supplementary Table S11.
    # fPPT = 0.39 is the paper's headline finding: for the same mass,
    # plaque exerts only about 39% of the protofibril effect on tau.
    # =====================================================================
    lfPFA <- fixed(log(1))      ; label("Relative contribution of protofibril to amyloid aggregation fPFA (log-dimensionless)") # Table S12 fPFA = 1 (Fixed)
    lfPFT <- fixed(log(1))      ; label("Relative contribution of protofibril to the tau pathway fPFT (log-dimensionless)")     # Table S12 fPFT = 1 (Fixed)
    lfPPT <- log(0.39)          ; label("Relative contribution of plaque to the tau pathway fPPT (log-dimensionless)")          # Table S12 fPPT = 0.39 (RSE 0.193%)
    lfPFP <- fixed(log(0.0001)) ; label("Relative contribution of protofibril to amyloid PET SUVR fPFP (log-dimensionless)")    # Table S12 fPFP = 0.0001 (Fixed)
    lfPFC <- fixed(log(0.001))  ; label("Relative contribution of protofibril to the direct cognitive effect fPFC (log-dimensionless)") # Table S12 fPFC = 0.001 (Calibrated)

    # =====================================================================
    # LECANEMAB DRUG EFFECT
    # Source: Supplementary Table S13 ("Parameters of Drug Effect Model").
    # Both act as additive clearance enhancement driven by CSS_LEC.
    # =====================================================================
    lDELECAF <- log(0.068)  ; label("Lecanemab drug effect on protofibril clearance DELECAF (log-mL/ug)") # Table S13 DELECAF = 0.068 (RSE 10.1%)
    lDELECAP <- log(0.0015) ; label("Lecanemab drug effect on plaque clearance DELECAP (log-mL/ug)")      # Table S13 DELECAP = 0.0015 (RSE 24.7%); this is the value at AGE = 0 -- see e_AGE_DELECAP

    # =====================================================================
    # COVARIATE EFFECTS -- Supplementary Table S14 ("Final QSP Model
    # Covariates and Parameter Values"). All are Monolix beta coefficients
    # on the log scale of a log-normally distributed parameter, i.e.
    # param_i = param_pop * exp(sum(beta_k * cov_k)) * exp(eta_i).
    # Continuous AGE enters UNCENTRED (see covariateData$AGE for the
    # three back-calculations against Supplementary Table S15 that
    # establish this).
    # =====================================================================
    e_APOE4HET_kAggAO   <- 0.15    ; label("Effect of APOE-e4 heterozygosity on kAggAO (log-scale)")            # Table S14 beta_kAggAO_APOEGEN_Hetero = 0.15 (RSE 5.22%)
    e_APOE4HOM_kAggAO   <- 0.18    ; label("Effect of APOE-e4 homozygosity on kAggAO (log-scale)")              # Table S14 beta_kAggAO_APOEGEN_Homo = 0.18 (RSE 6.93%)
    e_MCI_kAggAO        <- 0.12    ; label("Effect of baseline MCI diagnosis on kAggAO (log-scale)")            # Table S14 beta_kAggAO_DIAG_MCI = 0.12 (RSE 9.14%)
    e_MILD_kAggAO       <- 0.032   ; label("Effect of baseline mild-AD-dementia diagnosis on kAggAO (log-scale)") # Table S14 beta_kAggAO_DIAG_Mild_AD = 0.032 (RSE 34.4%)
    e_AGE_PABFLOOR      <- -0.0019 ; label("Effect of age on PABFLOOR (log-scale per year, uncentred)")         # Table S14 beta_PABFLOOR_AGE = -0.0019 (RSE 13.0%)
    e_ASIAN_PABFLOOR    <- 0.038   ; label("Effect of Asian race on PABFLOOR (log-scale)")                      # Table S14 beta_PABFLOOR_RACE_Asian = 0.038 (RSE 13.4%)
    e_APOE4CAR_alphaCdr <- -0.21   ; label("Effect of APOE-e4 carriage on alpha_cdr (log-scale)")               # Table S14 beta_alpha_tAPOECarrier = -0.21 (RSE 9.00%)
    e_DEM_alphaCdr      <- 0.39    ; label("Effect of the ADNI AD-dementia diagnosis on alpha_cdr (log-scale)") # Table S14 beta_alpha_bDIAG_AD = 0.39 (RSE 7.37%)
    e_MILD_alphaCdr     <- 0.18    ; label("Effect of mild dementia due to AD on alpha_cdr (log-scale)")        # Table S14 beta_alpha_bDIAG_Mild_AD = 0.18 (RSE 10.8%)
    e_APOE4CAR_alphaAdas<- -0.097  ; label("Effect of APOE-e4 carriage on alpha_adas (log-scale)")              # Table S14 beta_alpha_adas_tAPOECarrier = -0.097 (RSE 12.1%)
    e_ASIAN_alphaAdas   <- 0.14    ; label("Effect of Asian race on alpha_adas (log-scale)")                    # Table S14 beta_alpha_adas_RACE_Asian = 0.14 (RSE 11.2%)
    e_DEM_alphaAdas     <- 0.39    ; label("Effect of the ADNI AD-dementia diagnosis on alpha_adas (log-scale)")# Table S14 beta_alpha_adas_bDIAG_AD = 0.39 (RSE 8.11%)
    e_LMCI_alphaAdas    <- 0.13    ; label("Effect of baseline LMCI diagnosis on alpha_adas (log-scale)")       # Table S14 beta_alpha_adas_bDIAG_LMCI = 0.13 (RSE 22.4%)
    e_MCI_alphaAdas     <- 0.14    ; label("Effect of baseline MCI diagnosis on alpha_adas (log-scale)")        # Table S14 beta_alpha_adas_bDIAG_MCI = 0.14 (RSE 20.0%)
    e_MILD_alphaAdas    <- 0.16    ; label("Effect of mild dementia due to AD on alpha_adas (log-scale)")       # Table S14 beta_alpha_adas_bDIAG_Mild_AD = 0.16 (RSE 17.8%)
    e_AGE_shapeAdas     <- -0.021  ; label("Effect of age on shape_adas (log-scale per year, uncentred)")       # Table S14 beta_shape_adas_AGE = -0.021 (RSE 8.66%)
    e_AGE_DELECAP       <- 0.041   ; label("Effect of age on DELECAP (log-scale per year, uncentred)")          # Table S14 beta_DELECAP_AGE = 0.041 (RSE 8.14%)

    # =====================================================================
    # INTER-INDIVIDUAL VARIABILITY -- 32 parameters, matching the paper's
    # statement that "32 were estimated with inter-subject variability".
    # Values are the "Std dev of the random effects" columns of
    # Supplementary Tables S3, S6, S9 and S13, squared to give variances
    # on the log scale (Monolix log-normal parameter distributions).
    #
    # The correlated block reproduces the two off-diagonal correlations in
    # Supplementary Table S14: corr(kAggAO, alpha_cdr) = -0.24 and
    # corr(pNA, alpha_cdr) = -0.32. No correlation is reported between
    # pNA and kAggAO, so that covariance is fixed to zero. The block
    # covariances are corr * sd_i * sd_j.
    # =====================================================================
    etalpNA + etalkAggAO + etalalphaCdr ~
      c(0.13^2,
        fixed(0),                 0.093^2,
        -0.32 * 0.13 * 0.3, -0.24 * 0.093 * 0.3, 0.3^2)   # Table S3 pNA SD 0.13, kAggAO SD 0.093; Table S9 alpha_cdr SD 0.3; Table S14 corr_kAggAO_alpha_cdr = -0.24 (RSE 25.5%), corr_pNA_alpha_cdr = -0.32 (RSE 12.0%)

    etalkDeAggOA   ~ 0.01^2    # Table S3 kDeAggOA SD 0.01 (RSE 8.76%)
    etalkAggFP     ~ 0.78^2    # Table S3 kAggFP SD 0.78 (RSE 4.06%)
    etalkDeAggPF   ~ 0.26^2    # Table S3 kDeAggPF SD 0.26 (RSE 7.92%) -- IIV estimated even though the fixed effect was held at 0.025
    etalsA         ~ 0.0031^2  # Table S3 sA SD 0.0031 (RSE 19.3%)
    etalsO         ~ 0.0025^2  # Table S3 sO SD 0.0025 (RSE 7.92%)
    etalsF         ~ 0.15^2    # Table S3 sF SD 0.15 (RSE 10.7%) -- IIV estimated even though the fixed effect was held at 6.9E-08
    etalsP         ~ 0.0047^2  # Table S3 sP SD 0.0047 (RSE 9.45%)
    etalcF         ~ 0.0057^2  # Table S3 cF SD 0.0057 (RSE 10.2%)
    etalcP         ~ 0.0022^2  # Table S3 cP SD 0.0022 (RSE 8.18%)
    etalABSUVR0    ~ 0.07^2    # Table S3 ABSUVR0 SD 0.07 (RSE 2.80%)
    etalPABFLOOR   ~ 0.072^2   # Table S3 PABFLOOR SD 0.072 (RSE 2.18%)
    etalsTauS      ~ 0.0047^2  # Table S6 sTauS SD 0.0047 (RSE 10.5%)
    etalsTauP      ~ 0.021^2   # Table S6 sTauP SD 0.021 (RSE 21.1%)
    etaltI01       ~ 0.15^2    # Table S6 tI01 SD 0.15 (RSE 2.07%)
    etaltI12       ~ 0.0048^2  # Table S6 tI12 SD 0.0048 (RSE 7.21%)
    etalpSI0       ~ 0.019^2   # Table S6 pSI0 SD 0.019 (RSE 10.9%)
    etalpSI1       ~ 0.0077^2  # Table S6 pSI1 SD 0.0077 (RSE 16.0%)
    etalpPI0       ~ 0.53^2    # Table S6 pPI0 SD 0.53 (RSE 3.79%)
    etalpPI1       ~ 0.56^2    # Table S6 pPI1 SD 0.56 (RSE 3.44%)
    etalPPTAUFLOOR ~ 0.28^2    # Table S6 PPTAUFLOOR SD 0.28 (RSE 10.3%)
    etalTAUSUVR0   ~ 0.022^2   # Table S6 TAUSUVR0 SD 0.022 (RSE 8.73%)
    etalsI2TAUPET  ~ 0.37^2    # Table S6 sI2TAUPET SD 0.37 (RSE 4.24%)
    etalshapeCdr   ~ 0.85^2    # Table S9 shape_cdr SD 0.85 (RSE 1.77%)
    etalehcCdr     ~ 0.12^2    # Table S9 ehc_cdr SD 0.12 (RSE 7.49%)
    etalalphaAdas  ~ 0.2^2     # Table S9 alpha_adas SD 0.2 (RSE 2.39%)
    etalshapeAdas  ~ 0.73^2    # Table S9 shape_adas SD 0.73 (RSE 1.78%)
    etalehcAdas    ~ 0.053^2   # Table S9 ehc_adas SD 0.053 (RSE 11.3%)
    etalDELECAF    ~ 1.89^2    # Table S13 DELECAF SD 1.89 (RSE 5.22%)
    etalDELECAP    ~ 0.53^2    # Table S13 DELECAP SD 0.53 (RSE 4.73%)

    # =====================================================================
    # RESIDUAL ERROR -- Supplementary Table S14, "Final QSP Error Model
    # Parameters". Monolix reports one a (additive) and/or b
    # (proportional) term per endpoint. The paper does not state whether
    # the two-term endpoints used Monolix's combined1 (sd = a + b*f) or
    # combined2 (sd = sqrt(a^2 + (b*f)^2)) form; nlmixr2's add() + prop()
    # is combined2. See the vignette Errata.
    # =====================================================================
    propSd_PAB       <- 0.081 ; label("Proportional residual SD, plasma A-beta 42/40 ratio")  # Table S14 "Plasma Ab42/40 Ratio Proportional (bpabr)" = 0.081 (RSE 0.716%)
    addSd_CENTILOID  <- 4.47  ; label("Additive residual SD, amyloid PET (Centiloid)")        # Table S14 "Amyloid Centiloid Additive (acent)" = 4.47 (RSE 2.69%)
    propSd_CENTILOID <- 0.091 ; label("Proportional residual SD, amyloid PET (Centiloid)")    # Table S14 "Amyloid Centiloid Proportional (bcent)" = 0.091 (RSE 3.03%)
    propSd_PPTAU     <- 0.19  ; label("Proportional residual SD, plasma p-tau181")           # Table S14 "Plasma pTau181 Proportional (bpptau)" = 0.19 (RSE 0.920%)
    propSd_TAUSUVR   <- 0.064 ; label("Proportional residual SD, medial-temporal tau PET SUVR") # Table S14 "Tau PET Proportional (b10)" = 0.064 (RSE 3.11%)
    addSd_CDRSB      <- 0.52  ; label("Additive residual SD, CDR-SB (points)")                # Table S14 "CDR-SB Additive (acdr)" = 0.52 (RSE 1.08%)
    propSd_CDRSB     <- 0.16  ; label("Proportional residual SD, CDR-SB")                     # Table S14 "CDR-SB Proportional (bcdr)" = 0.16 (RSE 1.09%)
    addSd_ADASCOG    <- 2.98  ; label("Additive residual SD, ADAS-Cog (points)")              # Table S14 "ADAS-Cog Additive (aadas)" = 2.98 (RSE 0.923%)
    propSd_ADASCOG   <- 0.07  ; label("Proportional residual SD, ADAS-Cog")                   # Table S14 "ADAS-Cog Proportional (badas)" = 0.07 (RSE 2.46%)
  })

  model({
    # ---------------------------------------------------------------------
    # 1. Individual parameters. Monolix log-normal distributions, so every
    #    parameter is exp(log-theta + covariate terms + eta). Covariate
    #    coefficients follow Supplementary Table S14; AGE is uncentred.
    # ---------------------------------------------------------------------
    A0       <- exp(lA0)
    N0       <- exp(lN0)

    pNA      <- exp(lpNA + etalpNA)
    kAggAO   <- exp(lkAggAO +
                      e_APOE4HET_kAggAO * APOE4_HET +
                      e_APOE4HOM_kAggAO * APOE4_HOM +
                      e_MCI_kAggAO      * DIS_AD_MCI +
                      e_MILD_kAggAO     * DIS_AD_MILD +
                      etalkAggAO)
    kDeAggOA <- exp(lkDeAggOA + etalkDeAggOA)
    kAggOF   <- exp(lkAggOF)
    kDeAggFO <- exp(lkDeAggFO)
    kAggFP   <- exp(lkAggFP + etalkAggFP)
    kDeAggPF <- exp(lkDeAggPF + etalkDeAggPF)
    sA       <- exp(lsA + etalsA)
    sO       <- exp(lsO + etalsO)
    sF       <- exp(lsF + etalsF)
    sP       <- exp(lsP + etalsP)
    cA       <- exp(lcA)
    cO       <- exp(lcO)
    cF       <- exp(lcF + etalcF)
    cP       <- exp(lcP + etalcP)
    sC2PAb   <- exp(lsC2PAb)
    PABFLOOR <- exp(lPABFLOOR +
                      e_AGE_PABFLOOR   * AGE +
                      e_ASIAN_PABFLOOR * RACE_ASIAN +
                      etalPABFLOOR)
    ABSUVR0  <- exp(lABSUVR0 + etalABSUVR0)

    tN         <- exp(ltN)
    betaTS     <- exp(lbetaTS)
    sN1        <- exp(lsN1)
    sN2        <- exp(lsN2)
    rI0        <- exp(lrI0)
    rI1        <- exp(lrI1)
    rI2        <- exp(lrI2)
    tI01       <- exp(ltI01 + etaltI01)
    tI12       <- exp(ltI12 + etaltI12)
    eFT01      <- exp(leFT01)
    eFT12      <- exp(leFT12)
    eFT2d      <- exp(leFT2d)
    rho01      <- exp(lrho01)
    rho12      <- exp(lrho12)
    rho2d      <- exp(lrho2d)
    dI2        <- exp(ldI2)
    pSI0       <- exp(lpSI0 + etalpSI0)
    pSI1       <- exp(lpSI1 + etalpSI1)
    pSI2       <- exp(lpSI2)
    pPI0       <- exp(lpPI0 + etalpPI0)
    pPI1       <- exp(lpPI1 + etalpPI1)
    pPI2       <- exp(lpPI2)
    pPN        <- exp(lpPN)
    cTS        <- exp(lcTS)
    cTP        <- exp(lcTP)
    sTauS      <- exp(lsTauS + etalsTauS)
    sTauP      <- exp(lsTauP + etalsTauP)
    sC2PPtau   <- exp(lsC2PPtau)
    PPTAUFLOOR <- exp(lPPTAUFLOOR + etalPPTAUFLOOR)
    TAUSUVR0   <- exp(lTAUSUVR0 + etalTAUSUVR0)
    sI2TAUPET  <- exp(lsI2TAUPET + etalsI2TAUPET)
    fI1T       <- exp(lfI1T)
    ShapeTau   <- exp(lShapeTau)

    CDR0      <- exp(lCDR0)
    alphaCdr  <- exp(lalphaCdr +
                       e_APOE4CAR_alphaCdr * APOE4_CARRIER +
                       e_DEM_alphaCdr      * DIS_AD_DEMENTIA +
                       e_MILD_alphaCdr     * DIS_AD_MILD +
                       etalalphaCdr)
    ehcCdr    <- exp(lehcCdr + etalehcCdr)
    shapeCdr  <- exp(lshapeCdr + etalshapeCdr)
    dcdnCdr   <- exp(ldcdnCdr)
    dc2Cdr    <- exp(ldc2Cdr)
    dc1Cdr    <- exp(ldc1Cdr)
    ADAS0     <- exp(lADAS0)
    alphaAdas <- exp(lalphaAdas +
                       e_APOE4CAR_alphaAdas * APOE4_CARRIER +
                       e_ASIAN_alphaAdas    * RACE_ASIAN +
                       e_DEM_alphaAdas      * DIS_AD_DEMENTIA +
                       e_LMCI_alphaAdas     * DIS_AD_LMCI +
                       e_MCI_alphaAdas      * DIS_AD_MCI +
                       e_MILD_alphaAdas     * DIS_AD_MILD +
                       etalalphaAdas)
    ehcAdas   <- exp(lehcAdas + etalehcAdas)
    shapeAdas <- exp(lshapeAdas + e_AGE_shapeAdas * AGE + etalshapeAdas)
    dcdnAdas  <- exp(ldcdnAdas)
    dc2Adas   <- exp(ldc2Adas)
    dc1Adas   <- exp(ldc1Adas)

    fPFA <- exp(lfPFA)
    fPFT <- exp(lfPFT)
    fPPT <- exp(lfPPT)
    fPFP <- exp(lfPFP)
    fPFC <- exp(lfPFC)

    DELECAF <- exp(lDELECAF + etalDELECAF)
    DELECAP <- exp(lDELECAP + e_AGE_DELECAP * AGE + etalDELECAP)

    # ---------------------------------------------------------------------
    # 2. Anti-amyloid drug effect (Methods, "Full ODE system"). The paper
    #    writes an if / elseif / else block over two treatment windows that
    #    sets Eff_LECAF = 1 + DELECAF * CSSAV and Eff_LECAP = 1 + DELECAP *
    #    CSSAV while treatment is on board and 1 otherwise. Driving both
    #    from a single time-varying CSS_LEC covariate that is 0 off
    #    treatment is algebraically identical and generalises to any
    #    number of treatment / washout windows.
    # ---------------------------------------------------------------------
    effLecaF <- 1 + DELECAF * CSS_LEC
    effLecaP <- 1 + DELECAP * CSS_LEC

    # ---------------------------------------------------------------------
    # 3. Cooperativity terms. ampA is the protofibril-plus-plaque driver of
    #    amyloid aggregation (Eq 1, "kAggAO * A * (1 + sA * (fPFA * F +
    #    P))"); ampT is the corresponding driver of tau pathology (Eq 2,
    #    "(1 + eta01 * (fPFT * F + fPPT * P) + rho01 * TS)"). fPFA and fPFT
    #    are both fixed to 1 and fPPT = 0.39 -- the paper's finding that
    #    plaque is about 39% as neurotoxic as protofibril per unit mass.
    # ---------------------------------------------------------------------
    ampA <- fPFA * F + P
    ampT <- fPFT * F + fPPT * P

    # ---------------------------------------------------------------------
    # 4. Amyloid pathway ODEs -- Eq 1 / Methods "Full ODE system".
    # ---------------------------------------------------------------------
    d/dt(A) <- pNA * (N + I0 + I1 + I2) -
      kAggAO * A * (1 + sA * ampA) + kDeAggOA * O - cA * A
    d/dt(O) <- kAggAO * A * (1 + sA * ampA) -
      kAggOF * O * (1 + sO * ampA) + kDeAggFO * F - kDeAggOA * O - cO * O
    d/dt(F) <- kAggOF * O * (1 + sO * ampA) -
      kAggFP * F * (1 + sF * ampA) + kDeAggPF * P - kDeAggFO * F -
      effLecaF * cF * F
    d/dt(P) <- kAggFP * F * (1 + sF * ampA) - kDeAggPF * P -
      effLecaP * cP * P

    # ---------------------------------------------------------------------
    # 5. Tau pathway ODEs -- Eq 2 / Methods "Full ODE system". Named flux
    #    intermediates keep each transition written exactly once, so the
    #    "out of state i" and "into state i+1" terms cannot drift apart.
    # ---------------------------------------------------------------------
    fluxN01 <- tN * N * (1 + sN1 * ampT) +
      betaTS * N * TS * (1 + sN2 * ampT)
    flux01  <- tI01 * I0 * (1 + eFT01 * ampT + rho01 * TS)
    flux12  <- tI12 * I1 * (1 + eFT12 * ampT + rho12 * TS)
    flux2d  <- dI2  * I2 * (1 + eFT2d * ampT + rho2d * TS)

    d/dt(N)  <- -fluxN01 + rI0 * I0
    d/dt(I0) <- fluxN01 - flux01 - rI0 * I0 + rI1 * I1
    d/dt(I1) <- flux01 - flux12 + rI2 * I2 - rI1 * I1
    d/dt(I2) <- flux12 - flux2d - rI2 * I2
    d/dt(DN) <- flux2d
    d/dt(TS) <- (pSI0 * I0 + pSI1 * I1 + pSI2 * I2) * (1 + sTauS * ampT) -
      cTS * TS
    d/dt(TP) <- (pPI0 * I0 + pPI1 * I1 + pPI2 * I2) * (1 + sTauP * ampT) +
      pPN * N - cTP * TP

    # ---------------------------------------------------------------------
    # 6. Initial conditions -- Methods, "Base model development stages":
    #    at disease onset A-beta monomer production has started but no
    #    aggregate has formed, and all neurons are healthy, so only A and N
    #    are non-zero.
    # ---------------------------------------------------------------------
    A(0) <- A0
    N(0) <- N0

    # ---------------------------------------------------------------------
    # 7. Observation models -- Eq 3 and Methods "Algebraic equation to
    #    derive clinical outcomes".
    # ---------------------------------------------------------------------
    # Cognitive composites and scores. The paper caps CDR-SB at 18 points
    # and ADAS-Cog at 90 points, the maxima of the two scales.
    compEffCdr  <- (dcdnCdr  * DN + dc2Cdr  * I2 + dc1Cdr  * I1) / N0
    compEffAdas <- (dcdnAdas * DN + dc2Adas * I2 + dc1Adas * I1) / N0
    CDRSB   <- min(CDR0  + alphaCdr  * (1 + ehcCdr  * (fPFC * F + P)) *
                     compEffCdr^shapeCdr, 18)
    ADASCOG <- min(ADAS0 + alphaAdas * (1 + ehcAdas * (fPFC * F + P)) *
                     compEffAdas^shapeAdas, 90)

    # Amyloid PET. The SUVR -> Centiloid conversion constants 236.6 and
    # 246.9 are the paper's own (Methods); the separate ADNI conversion
    # Centiloid = 196.9 * ABSUVR - 196.03 was used only to harmonise the
    # ADNI florbetapir data into the modelling dataset, not in the model.
    ABSUVR    <- sP * (fPFP * F + P) + ABSUVR0
    CENTILOID <- ABSUVR * 236.6 - 246.9

    # Medial-temporal tau PET SUVR: I1 and I2 neurons carry the tracer
    # signal, I0 and DN contribute negligibly.
    TAUSUVR <- (sI2TAUPET * (fI1T * I1 + I2))^ShapeTau + TAUSUVR0

    # Fluid biomarkers. CSFAB and CSFPTAU are predicted, NOT fitted -- the
    # paper held CSF A-beta 42 and CSF p-tau181 out of the training set
    # and overlaid the predictions post hoc (Supplementary Figures S3 and
    # S4, which use the further display scalings CSFAB42 = 0.5 * A + 120
    # and CSFPTAU181 = 2.1 * TP). Because they carry no residual-error
    # model they are plain observables here.
    CSFAB   <- A
    PAB     <- sC2PAb * CSFAB + PABFLOOR
    CSFPTAU <- TP
    PPTAU   <- sC2PPtau * CSFPTAU + PPTAUFLOOR

    # ---------------------------------------------------------------------
    # 8. Residual error, one endpoint per fitted observable
    #    (Supplementary Table S14, "Final QSP Error Model Parameters").
    # ---------------------------------------------------------------------
    PAB       ~ prop(propSd_PAB)
    CENTILOID ~ add(addSd_CENTILOID) + prop(propSd_CENTILOID)
    TAUSUVR   ~ prop(propSd_TAUSUVR)
    PPTAU     ~ prop(propSd_PPTAU)
    CDRSB     ~ add(addSd_CDRSB) + prop(propSd_CDRSB)
    ADASCOG   ~ add(addSd_ADASCOG) + prop(propSd_ADASCOG)
  })
}
