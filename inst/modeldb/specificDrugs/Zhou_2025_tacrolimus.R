Zhou_2025_tacrolimus <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and",
    "elimination for oral tacrolimus in Chinese adult lung transplant",
    "recipients during the first three months post-transplantation",
    "(Zhou 2025 Eur J Clin Pharmacol, Phoenix NLME 8.2). Apparent oral",
    "clearance CL/F carries four covariate effects: median-normalised power",
    "effects of haematocrit and postoperative day, an exponential",
    "three-level CYP3A5*3 (rs776746) genotype effect with the *3/*3 poor",
    "metaboliser as reference, and an exponential azole-antifungal effect",
    "with a separate coefficient for each of voriconazole, posaconazole and",
    "itraconazole. Absorption rate constant fixed at 4.48 1/h from the",
    "literature because only trough concentrations were available; no",
    "covariate was retained on Vd/F. Exponential IIV on both CL/F and",
    "Vd/F, combined proportional-plus-additive residual error."
  )
  reference <- paste(
    "Zhou S, Lian Q, Luo H, Xie H, Guan Y, He J, Wei L, Ju C (2025).",
    "Population pharmacokinetic characteristics of tacrolimus in Chinese",
    "lung transplant recipients and optimisation of dosing regimen during",
    "the early post-transplantation phase.",
    "Eur J Clin Pharmacol 81(12):1841-1852.",
    "doi:10.1007/s00228-025-03920-9."
  )
  vignette <- "Zhou_2025_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Zhou 2025 Methods "Blood concentration detection":
  # "the whole blood concentration of tacrolimus was detected by
  # chemiluminescence microparticle immunoassay (CMIA)".
  compartmentData <- list(
    depot   = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    HCT = list(
      description        = "Haematocrit",
      units              = "% (volume fraction times 100)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying laboratory value collected throughout the treatment",
        "period (Methods 'General clinical data': 'laboratory indexes",
        "during treatment: ... haematocrit (HCT)'). Enters CL/F as the",
        "median-normalised power form (HCT / 26.5)^-0.71 per Zhou 2025",
        "Eq. 4.",
        "",
        "SCALE WARNING. Zhou 2025 reports haematocrit as a volume FRACTION",
        "(Table 1 cohort mean +/- SD 0.28 +/- 0.05; Eq. 4 normalising",
        "median 0.265), whereas the canonical HCT register column is on the",
        "PERCENT scale. The paper's median 0.265 is therefore recorded here",
        "as the reference 26.5 %, and a dataset supplied to this model must",
        "put HCT on the percent scale (e.g. 28, not 0.28). Because the",
        "effect is a ratio, mixing the two scales between the column and",
        "the reference would misscale CL/F by 100^-0.71 = 1/26.3, i.e. a",
        "26-fold error - the two must always be checked together.",
        "",
        "Mechanism (Discussion): tacrolimus is extensively bound to red",
        "blood cells, so 'Recipients with low levels of HCT have more free",
        "tacrolimus in the plasma, resulting in a corresponding increase in",
        "clearance', giving the negative exponent on a whole-blood-assay",
        "clearance. Reported sensitivity: 'When HCT decreased from 0.30 to",
        "0.20, tacrolimus clearance elevated by 33.36%' - reproduced",
        "exactly by (0.20 / 0.30)^-0.71 = 1.3336, which is the check that",
        "confirms the power (rather than exponential) reading of Eq. 4 for",
        "this covariate. The cohort is markedly anaemic relative to a",
        "healthy population (Hb 90.73 +/- 17.86 g/L, RBC 2.96 +/- 0.63e9),",
        "which is why the normalising median (26.5 %) sits far below the",
        "45 % used by Nestorov 2014."
      ),
      source_name        = "HCT"
    ),
    POD = list(
      description        = "Postoperative day (days elapsed since lung transplantation)",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying within subject. Enters CL/F as the median-normalised",
        "power form (POD / 49)^0.14 per Zhou 2025 Eq. 4. Table 1 cohort",
        "mean +/- SD 52.76 +/- 17.13 days; the Eq. 4 normalising median is",
        "49 days.",
        "",
        "The exponent is POSITIVE, i.e. apparent oral clearance RISES with",
        "time post-transplant - the opposite direction to the kidney",
        "transplant model Bergmann_2014_tacrolimus.R. Discussion: 'a",
        "significant elevation of 16.63% in tacrolimus clearance between",
        "the 30-day and 90-day periods after transplantation. It may be",
        "attributed to the gradual restoration of gastrointestinal motility",
        "and CYP3A enzyme activity during the early postoperative phase.'",
        "That figure is reproduced exactly by (90 / 30)^0.14 = 1.1663,",
        "which is the check that confirms the power reading of Eq. 4 for",
        "this covariate.",
        "",
        "SIMULATION FLOOR. The power form is undefined at POD = 0:",
        "(0 / 49)^0.14 = 0, which collapses CL/F to zero and removes all",
        "elimination. Zhou 2025 sampled only after tacrolimus reached",
        "steady state ('approximately 48 h after first tacrolimus dose')",
        "and within three months of surgery, so keep POD >= 1 in any",
        "simulation and treat the model as calibrated over roughly",
        "POD 2-90 days."
      ),
      source_name        = "POD"
    ),
    CYP3A5_STAR1_HET = list(
      description        = "CYP3A5*1/*3 heterozygote indicator (one functional CYP3A5*1 allele)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5*3/*3 poor metaboliser, when paired with CYP3A5_STAR1_HOM = 0)",
      notes              = paste(
        "Time-fixed per subject (germline genotype). 1 = subject is",
        "CYP3A5*1/*3 heterozygote (genotype GA at rs776746, one functional",
        "*1 allele); 0 = otherwise (the union of *3/*3 poor metabolisers",
        "and *1/*1 rapid metabolisers; the paired indicator",
        "CYP3A5_STAR1_HOM flags the homozygous group). Zhou 2025 genotyped",
        "the CYP3A5 6986 A>G variant by PCR-RFLP (Methods 'Genotype",
        "detection'); Table 1 counts AA (*1/*1) 13 (9.15%), GA (*1/*3) 56",
        "(39.44%), GG (*3/*3) 73 (51.41%). The paired HET / HOM",
        "decomposition applies because Table 2 reports a distinct",
        "coefficient for each of the three strata (0 / 0.63 / 1.00), which",
        "the Discussion highlights as the paper's novel contribution:",
        "'this study ... further expands the study of the fast-metabolic",
        "population'. Note the allele orientation: at rs776746 the G",
        "allele is CYP3A5*3 (the loss-of-function variant), so GG is the",
        "poor-metaboliser reference and A-allele carriage confers",
        "function - the *1-carrier-equals-1 orientation of these",
        "canonicals."
      ),
      source_name        = "CYP3A5 *1/*3"
    ),
    CYP3A5_STAR1_HOM = list(
      description        = "CYP3A5*1/*1 homozygote indicator (two functional CYP3A5*1 alleles)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5*3/*3 poor metaboliser, when paired with CYP3A5_STAR1_HET = 0)",
      notes              = paste(
        "Time-fixed per subject (germline genotype). 1 = subject is",
        "CYP3A5*1/*1 rapid metaboliser (genotype AA at rs776746); 0 =",
        "otherwise. Only 13 of 142 recipients (9.15%) - see the",
        "CYP3A5_STAR1_HET notes for the rs776746 allele orientation and the",
        "full genotype counts. Zhou 2025 Table 2 reports the coefficient as",
        "1.00 with RSE 9.23% (bootstrap median 1.00, 95% CI 0.80-1.19),",
        "back-transformed to the 2.72-fold clearance ratio quoted in the",
        "Abstract and Eq. 4. The genotype distribution satisfied",
        "Hardy-Weinberg equilibrium (Methods 'Data statistics')."
      ),
      source_name        = "CYP3A5 *1/*1"
    ),
    CONMED_VORICONAZOLE = list(
      description        = "Concomitant voriconazole coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant azole antifungal)",
      notes              = paste(
        "Time-varying. 1 = voriconazole coadministered at the time of the",
        "trough measurement; 0 = no azole antifungal. Table 1 counts 741 of",
        "988 trough records on voriconazole, making it by far the most",
        "common azole in the cohort. Reduces tacrolimus CL/F by 38.21%",
        "(Discussion), which closely reproduces the 36.2% reduction",
        "reported by Cai et al. in an independent lung-transplant cohort.",
        "The three azole indicators are mutually exclusive in the Zhou 2025",
        "dataset (no recipient received two mould-active azoles",
        "concurrently), and the reference category for all three is the",
        "no-AFD state. Methods 'Administration regimens' records that no",
        "enrolled patient received caspofungin, Wuzhi capsule, macrolide",
        "antibiotics or calcium channel blockers, so those competing",
        "CYP3A / P-gp perpetrators do not confound the azole coefficients."
      ),
      source_name        = "Voriconazole"
    ),
    CONMED_POSACONAZOLE = list(
      description        = "Concomitant posaconazole coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant azole antifungal)",
      notes              = paste(
        "Time-varying. 1 = posaconazole coadministered at the time of the",
        "trough measurement; 0 = no azole antifungal. Table 1 counts 157 of",
        "988 trough records. Reduces tacrolimus CL/F by 26.30%",
        "(Discussion) - the smallest of the three azole effects. The",
        "Discussion attributes the small magnitude partly to formulation:",
        "'the patients enrolled in our study mainly used posaconazole as",
        "suspension, which shows smaller impact than delayed-release",
        "tablets', so this coefficient should not be transferred to a",
        "delayed-release-tablet cohort without adjustment. See the",
        "CONMED_VORICONAZOLE notes for the mutual-exclusivity and",
        "competing-perpetrator context."
      ),
      source_name        = "Posaconazole"
    ),
    CONMED_ITRACONAZOLE = list(
      description        = "Concomitant itraconazole coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant azole antifungal)",
      notes              = paste(
        "Time-varying. 1 = itraconazole coadministered at the time of the",
        "trough measurement; 0 = no azole antifungal. Table 1 counts only",
        "38 of 988 trough records, the smallest azole subgroup, yet the",
        "coefficient is the most precisely estimated of the three (RSE",
        "-6.41%). Reduces tacrolimus CL/F by 57.98% (Discussion) - the",
        "largest of the three effects, consistent with the potency ranking",
        "the Discussion cites ('At equimolar concentrations, ketoconazole",
        "is the most potent CYP3A4 inhibitor, followed by itraconazole,",
        "voriconazole, and fluconazole') and with itraconazole being 'the",
        "most lipophilic azole which is prone to enzymatic and",
        "transporter-mediated interactions'. See the CONMED_VORICONAZOLE",
        "notes for the mutual-exclusivity and competing-perpetrator",
        "context."
      ),
      source_name        = "Itraconazole"
    )
  )

  # Covariates that Zhou 2025 collected and screened but did NOT retain in
  # the final model (Methods 'General clinical data' lists the full screening
  # set; Table 3 shows the stepwise forward-inclusion / backward-elimination
  # path, which retained only CYP3A5, AFD, HCT and POD). Documented here so
  # the paper's covariate screen is preserved without carrying
  # declared-but-unreferenced convention warnings.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Collected (Methods 'General clinical data') but not retained.",
        "Table 1 cohort mean +/- SD 50.26 +/- 11.13 kg. The Discussion",
        "invokes the low body weight to explain why this cohort's CL/F",
        "(7.58 L/h) is lower than in published cohorts ('our study",
        "subjects had a comparatively lower body weight, which means that",
        "tacrolimus was distributed in a smaller body surface area,",
        "resulting in a lower clearance'), and closes by naming weight as",
        "unfinished business: 'In addition, weight, age, and renal",
        "function indexes might affect the drug concentration of",
        "tacrolimus, which needs to be further discussed.' No allometric",
        "scaling is present in the final model."
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Collected but not retained. Table 1 cohort mean +/- SD 57.10 +/-",
        "11.46 years. Like weight, age is used narratively to explain the",
        "low cohort clearance ('the subjects in our study were generally",
        "older patients, with a mean age of 57 years ... there is a prone",
        "for tacrolimus to be easily accumulated in the elderly",
        "population') and flagged for further study, but it carries no",
        "coefficient in Table 2."
      )
    ),
    SEXF = list(
      description = "Sex, female indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Collected but not retained. Table 1 reports 116 male / 26 female",
        "(18.3% female)."
      )
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Collected but not retained. Table 1 cohort mean +/- SD 57.05 +/-",
        "20.85 mL/min (serum creatinine 99.29 +/- 46.45 mmol/L as printed).",
        "Recipients with creatinine clearance below 15 mL/min were excluded",
        "at enrolment (Methods exclusion criterion 2). Named in the",
        "Discussion's further-work list alongside weight and age."
      )
    ),
    HGB = list(
      description = "Haemoglobin",
      units       = "g/L",
      type        = "continuous",
      notes       = paste(
        "Collected but not retained. Table 1 cohort mean +/- SD 90.73 +/-",
        "17.86 g/L. Strongly collinear with the retained HCT covariate;",
        "the Discussion frames the haematological effect entirely through",
        "HCT and red-cell binding."
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Collected but not retained. Table 1 cohort mean +/- SD 37.33 +/- 4.08 g/L."
    ),
    TBILI = list(
      description = "Total serum bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = paste(
        "Collected but not retained. Table 1 cohort mean +/- SD 13.35 +/-",
        "23.51 umol/L. Recipients with Child-Pugh C liver dysfunction were",
        "excluded at enrolment (Methods exclusion criterion 2). Contrast",
        "Pei_2023_tacrolimus.R, where total bilirubin IS retained on CL/F",
        "in a heart transplant cohort."
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Collected but not retained. Table 1 cohort mean +/- SD 27.38 +/- 33.32 IU/L."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Collected but not retained. Table 1 cohort mean +/- SD 27.93 +/- 21.52 IU/L."
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor use",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Collected but not retained. Table 1 counts 446 of 988 trough",
        "records on a proton pump inhibitor."
      )
    ),
    CONMED_STEROID = list(
      description = "Systemic corticosteroid administration indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Collected but not retained. Table 1 counts 952 of 988 trough",
        "records on a glucocorticoid (prednisone or methylprednisolone) -",
        "i.e. 96% of records, so the covariate is nearly constant across",
        "the dataset and cannot be identified."
      )
    ),
    CONMED_MPA = list(
      description = "Concomitant mycophenolic acid coadministration indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Collected but not retained. Table 1 counts 578 trough records on",
        "mycophenolate mofetil and 67 on enteric-coated mycophenolate",
        "sodium, out of 988. Part of the triple-immunotherapy backbone",
        "(Methods 'Administration regimens')."
      )
    ),
    TX_LUNG_BILATERAL = list(
      description = "Bilateral (double) lung transplantation indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Collected but not retained. Table 1 counts 74 single-lung and 68",
        "double-lung transplantations."
      )
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 142L,
    n_studies        = 1L,
    n_concentrations = 988L,
    age_range        = "adults; Table 1 mean +/- SD 57.10 +/- 11.46 years",
    weight_range     = "Table 1 mean +/- SD 50.26 +/- 11.13 kg",
    height_range     = "Table 1 mean +/- SD 165.07 +/- 6.76 cm",
    sex_female_pct   = 18.3,
    race_ethnicity   = c(Asian = 100),
    disease_state    = paste(
      "Chinese adult first-time lung transplant recipients at the First",
      "Affiliated Hospital of Guangzhou Medical University, January 2019 to",
      "December 2021, on tacrolimus-based triple immunosuppression",
      "(tacrolimus + mycophenolate + glucocorticoid). Primary indication:",
      "interstitial lung disease 58, chronic obstructive pulmonary disease",
      "44, bronchiectasis 15, other 25. Procedure: single lung 74,",
      "double lung 68. Excluded were combined-organ transplant recipients,",
      "those with severe systemic complications, haemodynamic dysfunction,",
      "Child-Pugh C liver dysfunction or creatinine clearance below",
      "15 mL/min, and those with poor compliance or incomplete records.",
      "The cohort is notably anaemic (haematocrit 0.28 +/- 0.05,",
      "haemoglobin 90.73 +/- 17.86 g/L) and of low body weight, both of",
      "which the Discussion invokes to explain the low apparent clearance",
      "relative to published cohorts."
    ),
    dose_range       = paste(
      "Oral tacrolimus capsules (Prograf, Astellas). Initial dose",
      "50-150 ug/kg/day administered every 12 h, then titrated on trough",
      "concentration to a target of 10-15 ng/mL. Table 1 mean +/- SD daily",
      "dose 2.51 +/- 1.37 mg/day; observed trough 12.82 +/- 4.92 ng/mL."
    ),
    regions          = "China (Guangzhou)",
    sampling_design  = paste(
      "988 steady-state whole-blood trough concentrations from 142",
      "recipients, collected by routine therapeutic drug monitoring within",
      "three months of transplantation. Sampling began after tacrolimus",
      "reached steady state (approximately 48 h after the first dose);",
      "each sample was drawn 30 min before the next dose from EDTA",
      "anticoagulated whole blood and assayed by chemiluminescence",
      "microparticle immunoassay (CMIA) on an Abbott ARCHITECT i1000SR",
      "with the Abbott 1L77 tacrolimus kit. Trough-only sampling is why",
      "Ka was fixed rather than estimated and why no absorption-phase or",
      "peripheral-compartment parameters are identifiable."
    ),
    cyp3a5_genotype  = c(`*3/*3 (GG)` = 51.41, `*1/*3 (GA)` = 39.44, `*1/*1 (AA)` = 9.15),
    notes            = paste(
      "Software: Phoenix NLME 8.2 (Certara); R 4.0.4 for diagnostics.",
      "Structural model selected on OFV (-2LL) and goodness-of-fit plots;",
      "one- and two-compartment models with several residual structures",
      "were compared and the one-compartment model with first-order",
      "absorption and elimination was retained. Covariates screened by",
      "univariate analysis then stepwise regression with forward inclusion",
      "(-2LL drop > 3.84) and backward elimination (-2LL rise < 6.64);",
      "Table 3 records the path (base OFV 5829.11 -> final 5647.30, with",
      "CYP3A5 contributing the largest single drop, -84.92). Evaluated by",
      "goodness-of-fit plots (Figure 1), a prediction- and",
      "variability-corrected VPC (Figure 2), and a 1000-replicate",
      "bootstrap (Table 2) in which every final estimate fell inside the",
      "bootstrap 95% CI with bias within +/- 20%. Ethics: First Affiliated",
      "Hospital of Guangzhou Medical University ES-2023-K010."
    )
  )

  ini({
    # ---- Structural parameters (Zhou 2025 Table 2, final model) ----------
    # Ka was FIXED at a literature value, not estimated: Table 2 prints
    # 4.48 with "-" in both RSE columns, and Methods "Establishment of
    # basic model" states "According to a comprehensive literature report
    # [13], the absorption rate constant (Ka) is fixed at 4.48 h-1". The
    # Discussion gives the reason: "the model had only trough concentration
    # data and no sampling in the absorption stage. Therefore, based on the
    # literature report, the absorption rate constant (Ka) was fixed at
    # 4.48 h-1, and 74% of relevant studies selected fixed parameters to
    # build the model".
    lka <- fixed(log(4.48)); label("First-order absorption rate constant (1/h), literature value")  # Zhou 2025 Table 2, theta_ka = 4.48 (no RSE reported); Methods "Establishment of basic model"

    # Apparent oral clearance and volume for the REFERENCE subject:
    # CYP3A5*3/*3 poor metaboliser, no azole antifungal, haematocrit at the
    # cohort median (0.265 = 26.5%) and postoperative day at the cohort
    # median (49 days). Both include the oral bioavailability F, which is
    # not separately identifiable from trough-only data (Methods: "the
    # clearance (CL/F) and apparent distribution volume (Vd/F) include the
    # bioavailability of oral preparations").
    lcl <- log(7.58);   label("Apparent oral clearance CL/F for the reference subject (L/h)")  # Zhou 2025 Table 2, theta_CL/F = 7.58 (RSE 10.60%, bootstrap median 7.55, 95% CI 6.14-8.83)
    lvc <- log(701.39); label("Apparent volume of distribution Vd/F (L)")                      # Zhou 2025 Table 2, theta_Vd/F = 701.39 (RSE 15.43%, bootstrap median 684.78, 95% CI 531.16-952.14)

    # ---- Continuous covariate effects on CL/F (Zhou 2025 Table 2) --------
    # Median-normalised power form (Eq. 2: theta_i = theta * (COV_i /
    # COV_median)^theta_COV), with the normalising medians given in Eq. 4.
    #
    # Haematocrit, reference 26.5% (the paper's 0.265 volume fraction
    # rescaled to the canonical percent column - see covariateData). The
    # power reading is confirmed by the paper's own sensitivity statement:
    # (0.20/0.30)^-0.71 = 1.3336, matching "When HCT decreased from 0.30 to
    # 0.20, tacrolimus clearance elevated by 33.36%" to four digits.
    e_hct_cl <- -0.71; label("Power exponent of (HCT / 26.5%) on CL/F (unitless)")  # Zhou 2025 Table 2, theta_HCT = -0.71 (RSE -19.03%, bootstrap median -0.70, 95% CI -0.88 to -0.35)

    # Postoperative day, reference 49 days. The power reading is confirmed
    # by (90/30)^0.14 = 1.1663, matching the Discussion's "significant
    # elevation of 16.63% in tacrolimus clearance between the 30-day and
    # 90-day periods after transplantation" to four digits.
    e_pod_cl <-  0.14; label("Power exponent of (POD / 49 days) on CL/F (unitless)")  # Zhou 2025 Table 2, theta_POD = 0.14 (RSE 34.90%, bootstrap median 0.13, 95% CI 0.02-0.22)

    # ---- Categorical covariate effects on CL/F (Zhou 2025 Table 2) -------
    # Exponential form (Eq. 3: theta_i = theta * exp(theta_cov)) on each 0/1
    # indicator. Table 2 reports the coefficients on the log scale; Eq. 4
    # prints their back-transformed multipliers, which is the cross-check
    # that fixes the exponential reading:
    #   exp(0.63) = 1.878 vs Eq. 4 "1.87" for CYP3A5 *1/*3
    #   exp(1.00) = 2.718 vs Eq. 4 "2.72" for CYP3A5 *1/*1
    #   exp(-0.48) = 0.619 vs Eq. 4 "0.62" for voriconazole
    #   exp(-0.31) = 0.733 vs Eq. 4 "0.74" for posaconazole
    #   exp(-0.87) = 0.419 vs Eq. 4 "0.42" for itraconazole
    # The Table 2 log-scale values are used here because they are the
    # estimated quantities carrying RSEs and bootstrap CIs; see the
    # vignette's Assumptions and deviations section for the one place the
    # rounding differs materially (posaconazole, 0.733 vs 0.74).

    # CYP3A5*3 (rs776746) three-level genotype; *3/*3 poor metaboliser is
    # the reference with both indicators = 0 and Table 2 coefficient 0.
    e_cyp3a5_het_cl <- 0.63; label("Exponential coefficient of CYP3A5*1/*3 on CL/F (unitless)")  # Zhou 2025 Table 2, theta_CYP3A5 *1/*3 = 0.63 (RSE 11.21%, bootstrap median 0.62, 95% CI 0.51-0.73)
    e_cyp3a5_hom_cl <- 1.00; label("Exponential coefficient of CYP3A5*1/*1 on CL/F (unitless)")  # Zhou 2025 Table 2, theta_CYP3A5 *1/*1 = 1.00 (RSE 9.23%, bootstrap median 1.00, 95% CI 0.80-1.19)

    # Azole antifungal drugs, one coefficient per agent; no-AFD is the
    # reference with all three indicators = 0.
    e_vori_cl <- -0.48; label("Exponential coefficient of concomitant voriconazole on CL/F (unitless)")   # Zhou 2025 Table 2, theta_AFD Voriconazole = -0.48 (RSE -22.31%, bootstrap median -0.46, 95% CI -0.57 to -0.25)
    e_posa_cl <- -0.31; label("Exponential coefficient of concomitant posaconazole on CL/F (unitless)")   # Zhou 2025 Table 2, theta_AFD Posaconazole = -0.31 (RSE -31.83%, bootstrap median -0.32, 95% CI -0.55 to -0.09)
    e_itra_cl <- -0.87; label("Exponential coefficient of concomitant itraconazole on CL/F (unitless)")   # Zhou 2025 Table 2, theta_AFD Itraconazole = -0.87 (RSE -6.41%, bootstrap median -0.86, 95% CI -1.05 to -0.54)

    # ---- Inter-individual variability (Zhou 2025 Table 2) ----------------
    # Eq. 1: P_i = theta_i * exp(eta_i), eta ~ N(0, omega^2). Table 2's
    # rows are labelled "omega^2 Vd" and "omega^2 CL", so the reported
    # numbers are variances already and no CV-to-variance conversion is
    # needed. Note the unusually large variance on Vd/F: omega^2 = 0.80
    # corresponds to sqrt(exp(0.80) - 1) = 111% CV, which is expected for a
    # volume identified from trough-only data.
    etalcl ~ 0.13   # Zhou 2025 Table 2, omega^2 CL = 0.13 (RSE 18.49%, bootstrap median 0.13)
    etalvc ~ 0.80   # Zhou 2025 Table 2, omega^2 Vd = 0.80 (RSE 21.22%, bootstrap median 0.68)

    # ---- Residual error (Zhou 2025 Table 2) ------------------------------
    # Mixed proportional-plus-additive ("the inter-individual and
    # intra-individual variation adopted exponential and a mixed
    # proportional-additive model, respectively"). Table 2 reports
    # "Proportional residual error (%) -21.58" and "Additive residual error
    # (ng/ml) 2.11". The proportional term is read as a magnitude of
    # 21.58% = 0.2158: the printed minus sign is a Phoenix NLME sign
    # convention artefact on a standard deviation (the RSE column for the
    # same row is also negative, -10.57%, and an RSE cannot be negative
    # either), and a residual SD is defined only up to sign.
    propSd <- 0.2158; label("Proportional residual error (fraction)")  # Zhou 2025 Table 2, epsilon_1 = -21.58% (RSE -10.57%, bootstrap median -22.54, 95% CI -26.91 to -18.37)
    addSd  <- 2.11;   label("Additive residual error (ng/mL)")         # Zhou 2025 Table 2, epsilon_2 = 2.11 ng/mL (RSE 17.34%, bootstrap median 1.78, 95% CI 0.05-2.53)
  })

  model({
    # ---- 1. Covariate factors on CL/F ------------------------------------
    # Zhou 2025 Eq. 4:
    #   CL/F (L/h) = 7.58 * (HCT/0.265)^-0.71 * (POD/49)^0.14
    #                     * CYP3A5 * AFDs
    # where CYP3A5 is 2.72 (*1/*1), 1.87 (*1/*3), 1 (*3/*3) and AFDs is
    # 0.62 (voriconazole), 0.74 (posaconazole), 0.42 (itraconazole), 1 (no
    # AFD). The HCT reference is written as 26.5 here because the canonical
    # HCT column is on the percent scale while the paper prints volume
    # fractions; the ratio is unchanged.
    hct_factor <- (HCT / 26.5)^e_hct_cl
    pod_factor <- (POD / 49)^e_pod_cl

    # Three-level CYP3A5*3 genotype, exponential form on the paired
    # mutually exclusive indicators. The *3/*3 reference has both
    # indicators = 0, so the factor collapses to 1.
    cyp3a5_factor <- exp(e_cyp3a5_het_cl * CYP3A5_STAR1_HET +
                         e_cyp3a5_hom_cl * CYP3A5_STAR1_HOM)

    # Azole antifungals, exponential form with one coefficient per agent.
    # The indicators are mutually exclusive in the source dataset, so at
    # most one term is nonzero; the no-AFD reference collapses the factor
    # to 1.
    afd_factor <- exp(e_vori_cl * CONMED_VORICONAZOLE +
                      e_posa_cl * CONMED_POSACONAZOLE +
                      e_itra_cl * CONMED_ITRACONAZOLE)

    # ---- 2. Individual parameters ----------------------------------------
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * hct_factor * pod_factor *
          cyp3a5_factor * afd_factor
    vc <- exp(lvc + etalvc)

    # ---- 3. ODE system (one compartment, first-order absorption) ---------
    kel <- cl / vc
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ---- 4. Observation and error ----------------------------------------
    # central is in mg and vc in L, so central / vc is mg/L = ug/mL; the
    # factor 1000 converts to the ng/mL reported by Zhou 2025.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
