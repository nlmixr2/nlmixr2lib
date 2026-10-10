Cao_2022_pemetrexed <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order elimination for",
    "intravenous pemetrexed (15-minute infusion, 500 mg/m^2) in Chinese adults",
    "with primary advanced non-small cell lung carcinoma (Cao 2022; 116",
    "patients, 192 plasma concentrations). Clearance (8.29 L/h at the cohort",
    "mean) scales as a power function of Cockcroft-Gault creatinine clearance",
    "normalised to 93.6 mL/min (exponent 0.58); intercompartmental clearance",
    "is raised in ERCC1 rs3212986 C/C homozygotes and in CYP3A5 rs776746 T/C",
    "(*1/*3) heterozygotes (exponential effects). Exponential interindividual",
    "variability on clearance and intercompartmental clearance (none on the",
    "volumes, dropped for >99% shrinkage) and a proportional residual error.",
    sep = " "
  )
  reference <- paste(
    "Cao P, Guo W, Wang J, Wu S, Huang Y, Wang Y, Liu Y, Zhang Y (2022).",
    "Population pharmacokinetic study of pemetrexed in chinese primary",
    "advanced non-small cell lung carcinoma patients. Frontiers in",
    "Pharmacology 13:954242. doi:10.3389/fphar.2022.954242.",
    "Final parameter estimates from Table 2 and the final-model equation in",
    "the Results ('Final PPK model validation').",
    sep = " "
  )
  vignette <- "Cao_2022_pemetrexed"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "pemetrexed", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "pemetrexed", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation (raw, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column CrCl. Calculated by the Cockcroft-Gault equation",
        "(Methods 'Dosage regimen and pharmacokinetic sampling'), raw mL/min,",
        "NOT BSA-normalized -- the raw Cockcroft-Gault variant of the CRCL",
        "canonical (Delattre_2010_amikacin.R / Gijsen_2022_meropenem.R",
        "precedents). Enters clearance as the power function (CRCL / 93.6)^0.58",
        "per the final-model equation; CrClMean = 93.6 mL/min is the cohort",
        "mean (Table 1 mean 93.6 (SD 26.5), median 89.5, range 47.5-179.7).",
        "Strongly correlated with CL (Figure 2A Spearman r = 0.6187,",
        "p < 0.0001)."
      ),
      source_name = "CrCl"
    ),
    SNP_ERCC1_RS3212986_CC = list(
      description = "ERCC1 rs3212986 (C8092A) homozygous wild-type C/C genotype indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (A/C or A/A at rs3212986)",
      notes = paste(
        "1 = C/C homozygote at ERCC1 rs3212986; 0 = A/C or A/A. Per the",
        "final-model text: 'If ERCC1 phenotype (rs3212986) = C/C, theta_ERCC1",
        "= theta_2 (0.83); If ERCC1 phenotype (rs3212986) = A/C or A/A,",
        "theta_ERCC1 = 0.' C/C homozygotes exhibited higher Q than the two",
        "variant genotypes (Figure 2B). Genotype counts A/A:A/C:C/C =",
        "15:56:45 (Supplementary Table S2). A plain mutant-allele-presence",
        "indicator would group A/C with A/A and flag the reference stratum,",
        "so the wild-type-homozygote-specific C/C indicator is used."
      ),
      source_name = "ERCC1 (rs3212986) C/C"
    ),
    CYP3A5_STAR1_HET = list(
      description = "CYP3A5 rs776746 (*3) heterozygote (*1/*3) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (C/C or T/T at rs776746, i.e. *3/*3 nonexpresser or *1/*1 homozygous expresser)",
      notes = paste(
        "1 = T/C heterozygote at CYP3A5 rs776746; 0 = C/C or T/T. Per the",
        "final-model text: 'If CYP3A5 phenotype (rs776746) = T/C,",
        "theta_CYP3A5 = theta_3 (0.62); If CYP3A5 phenotype (rs776746) = C/C",
        "or T/T, theta_CYP3A5 = 0.' rs776746 is CYP3A5*3; on the strand the",
        "paper reports, T = A = *1 (functional) and C = G = *3 (nonfunctional),",
        "so T/C = *1/*3 heterozygote, which is exactly the CYP3A5_STAR1_HET",
        "canonical (1 = *1/*3; 0 = union of *3/*3 and *1/*1). The paper's",
        "effect is heterozygote-specific (T/C showed higher Q than BOTH",
        "homozygotes, Figure 2C), matching the het-vs-both-homozygotes",
        "reference of this canonical. Genotype counts C/C:T/C:T/T = 60:49:7",
        "(Supplementary Table S2)."
      ),
      source_name = "CYP3A5 (rs776746) T/C"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Carried into the stepwise covariate search after collinearity screening and not retained. Cohort median 57.0 years (range 27.0-73.0) (Table 1)."
    ),
    BSA = list(
      description = "Body surface area (Mosteller equation)",
      units = "m^2",
      type = "continuous",
      notes = "Carried into the stepwise covariate search and not retained. Cohort median 1.7 m^2 (range 1.4-2.0) (Table 1). BSA drives the 500 mg/m^2 dose but was not a covariate on any PK parameter."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Carried into the stepwise covariate search and not retained. Cohort median 37.7 g/L (range 27.4-50.9) (Table 1)."
    ),
    WBC = list(
      description = "White blood cell count",
      units = "10^9/L",
      type = "continuous",
      notes = "Carried into the stepwise covariate search and not retained. Cohort median 7.0 (range 3.8-39.4) (Table 1)."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      notes = "Carried into the stepwise covariate search and not retained. Cohort median 91.0 U/L (range 31.0-730.0) (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 116L,
    n_studies = 1L,
    n_observations = 192L,
    age_median = "57.0 years (range 27.0-73.0)",
    sex_female_pct = 100 * 47 / 116,
    disease_state = paste(
      "Chinese adults with histologically diagnosed primary advanced",
      "non-small cell lung carcinoma (NSCLC), receiving at least two cycles",
      "of pemetrexed plus platinum (cisplatin/carboplatin/nedaplatin",
      "43:37:31) as primary treatment."
    ),
    race_ethnicity = "Chinese (Han, Wuhan).",
    renal_function = paste(
      "Cockcroft-Gault creatinine clearance median 89.5 mL/min (mean 93.6,",
      "SD 26.5, range 47.5-179.7) (Table 1). Serum creatinine median 66.0",
      "umol/L (range 36.8-123.5)."
    ),
    weight_range = "Body surface area median 1.7 m^2 (range 1.4-2.0); body weight not tabulated.",
    dose_range = paste(
      "Pemetrexed 500 mg/m^2 as a 15-minute intravenous infusion in 100 mL",
      "saline at the first chemotherapy cycle (about 850 mg at the median 1.7",
      "m^2 BSA)."
    ),
    regions = "China (Wuhan Union Hospital, single centre).",
    notes = paste(
      "Prospective single-centre cohort, February 2018 to December 2019.",
      "One or two plasma samples per patient (1-3 total) drawn at 0.5, 1, 3,",
      "5, 7, 24, 48 or 72 h after infusion (192 samples). UPLC-ESI-MS/MS",
      "assay, lower limit of quantification 2.50 ng/mL. All patients",
      "genotyped on the CBT-PMRA array for 17 continuous and 45 categorical",
      "covariates. Model built in Phoenix NLME 8.2 (FOCE-ELS). Random effects",
      "on V1 and V2 were dropped for high shrinkage (99.51% and 99.97%).",
      "Demographics from Table 1; parameters from Table 2."
    )
  )

  ini({
    # Structural parameters -- Table 2 final-model estimates. The final-model
    # equation parameterises in V1, V2, CL and Q with first-order elimination.
    lcl <- log(8.29); label("Clearance at CRCL = 93.6 mL/min (L/h)")  # Table 2 CL = 8.29 L/h (RSE 4.43%)
    lvc <- log(18.94); label("Central volume of distribution V1 (L)")  # Table 2 V1 = 18.94 L (RSE 5.94%)
    lq <- log(0.10); label("Intercompartmental clearance Q at reference genotype (L/h)")  # Table 2 Q = 0.10 L/h (RSE 18.03%)
    lvp <- log(5.12); label("Peripheral volume of distribution V2 (L)")  # Table 2 V2 = 5.12 L (RSE 18.85%)

    # Covariate effects (Table 2 theta_1/theta_2/theta_3; final-model equation).
    # CL = theta_CL * (CrCl/CrClMean)^theta_1 * exp(etaCL), CrClMean = 93.6.
    e_crcl_cl <- 0.58; label("Power exponent of CRCL on clearance (unitless)")  # Table 2 theta_1 = 0.58 (RSE 26.90%)
    # Q = theta_Q * exp(theta_2 * ERCC1_CC) * exp(theta_3 * CYP3A5_het) * exp(etaQ).
    e_ercc1_q <- 0.83; label("ERCC1 rs3212986 C/C exponential effect on Q")  # Table 2 theta_2 = 0.83 (RSE 29.66%)
    e_cyp3a5_q <- 0.62; label("CYP3A5 rs776746 *1/*3 exponential effect on Q")  # Table 2 theta_3 = 0.62 (RSE 36.45%)

    # Interindividual variability. Exponential model Pi = theta * exp(eta)
    # (Methods). Table 2 'IIV(omega%)' 5.61 (CL) and 24.17 (Q) are the Phoenix
    # Omega diagonal (variances) x 100, not SDs or CVs: a 5.6% CV on CL cannot
    # give the reported 13.89% eta-shrinkage from 1-3 samples with a 24.7%
    # residual, nor the ~+/-30% EBE scatter of CL about the CrCl trend in
    # Figure 2A, nor the Table 3 AUC 95% CI (e.g. 92.08-248.72 around 165.41,
    # a log SD of ~0.25 = sqrt(0.0561)), nor the 10th-90th percentile bands
    # of Figure 6. See the vignette variance and Figure 6 checks.
    # No IIV on V1/V2 (dropped for 99.51%/99.97% shrinkage).
    etalcl ~ 0.0561  # Table 2 IIV CL omega% = 5.61 -> variance 0.0561
    etalq ~ 0.2417   # Table 2 IIV Q omega% = 24.17 -> variance 0.2417

    # Residual error. A proportional model Y = IPRED * (1 + eps) was selected
    # (smallest OFV). Table 2 residual sigma(%) = 24.72 is read as the SD of
    # eps: Phoenix estimates the proportional error as a standard deviation
    # (CEps / stdev0), and the row carries an RSE and bootstrap CI like the
    # structural thetas.
    propSd <- 0.2472; label("Proportional residual error (fraction)")  # Table 2 sigma% = 24.72 (RSE 14.11%)
  })

  model({
    # 1. Individual PK parameters. Clearance scales with Cockcroft-Gault CrCl
    # by a power function normalised to the cohort mean 93.6 mL/min.
    cl <- exp(lcl + etalcl) * (CRCL / 93.6)^e_crcl_cl
    vc <- exp(lvc)
    # Intercompartmental clearance is raised in ERCC1 rs3212986 C/C homozygotes
    # and in CYP3A5 rs776746 *1/*3 heterozygotes (exponential effects).
    q <- exp(lq + etalq) *
      exp(e_ercc1_q * SNP_ERCC1_RS3212986_CC) *
      exp(e_cyp3a5_q * CYP3A5_STAR1_HET)
    vp <- exp(lvp)

    # 2. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. Two-compartment disposition with linear elimination. Pemetrexed was
    # given as a 15-minute intravenous infusion into the central compartment;
    # there is no depot.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Observation. Plasma pemetrexed in mg/L with proportional residual
    # error.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
