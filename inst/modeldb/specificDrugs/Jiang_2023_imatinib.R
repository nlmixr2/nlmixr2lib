Jiang_2023_imatinib <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order absorption and ",
    "first-order elimination for oral adjuvant imatinib in postoperative ",
    "Chinese adults with gastrointestinal stromal tumors (Jiang 2023). ",
    "The absorption rate constant is fixed at 1.22 1/h. Apparent oral ",
    "clearance CL/F carries two covariate effects: a power-form red blood ",
    "cell count effect ((RBC/3.7)^0.49) and a three-level ABCG2 rs2231142 ",
    "genotype effect encoded with paired binary indicators (GG wild-type ",
    "reference = 1, GT heterozygote = 0.879, TT homozygous variant = ",
    "0.976). Inter-individual variability is estimated on CL/F only; the ",
    "apparent volume of distribution V/F is a typical value with no IIV ",
    "because only trough concentrations were collected. Residual error is ",
    "proportional."
  )
  reference <- paste0(
    "Jiang X, Fu Q, Jing Y, Kong Y, Liu H, Peng H, Rexiti K, Wei X. ",
    "Personalized dose of adjuvant imatinib in patients with ",
    "gastrointestinal stromal tumors: results from a population ",
    "pharmacokinetic analysis. Drug Des Devel Ther. 2023;17:809-820. ",
    "doi:10.2147/DDDT.S400986."
  )
  vignette <- "Jiang_2023_imatinib"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    RBC = list(
      description        = "Red blood cell (erythrocyte) count",
      units              = "10^12 cells/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Baseline (per-subject) laboratory value from the electronic ",
        "medical record. Jiang 2023 Table 1 modeling dataset: mean +/- SD ",
        "3.7 +/- 0.6 x 10^12/L (range 2.2-5.6); validation dataset 3.8 +/- ",
        "0.5 (range 3.1-5.0). Enters CL/F as the power function ",
        "(RBC / 3.7)^0.49 (Jiang 2023 Equation 7 applied to the final-model ",
        "equation on page 814), where the reference 3.7 x 10^12/L is the ",
        "cohort central value defined by Jiang 2023 Covariate Model as ",
        "'Ref represents the median of the covariate values'. RBC was the ",
        "strongest single covariate in univariate forward inclusion ",
        "(dOFV = -12.50, p < 0.001) and was retained in backward ",
        "elimination (dOFV = +7.27, p < 0.01; Table 4). Jiang 2023 ",
        "Discussion: 'a reduction in RBC by 50% was connected with a 29% ",
        "drop in imatinib CL/F' -- 0.5^0.49 = 0.712, which confirms the ",
        "power form and the exponent sign. Hemoglobin (HGB) was screened ",
        "but is strongly correlated with RBC and was dropped before the ",
        "stepwise process (Jiang 2023 Covariate Model: 'Of the two ",
        "correlated covariates, only one was included to avoid ",
        "collinearity')."
      ),
      source_name        = "RBC"
    ),
    SNP_ABCG2_RS2231142_HET = list(
      description        = "ABCG2 (BCRP) rs2231142 heterozygote indicator (GT genotype)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ABCG2 rs2231142 GG wild-type homozygote, when paired with SNP_ABCG2_RS2231142_HOM = 0)",
      notes              = paste0(
        "Time-fixed per subject (germline genotype). 1 = subject carries ",
        "exactly one rs2231142 variant allele (genotype GT on the reported ",
        "strand, equivalently 421C/A); 0 = otherwise (the union of GG ",
        "wild-type homozygotes and TT homozygous-variant carriers; the ",
        "paired indicator SNP_ABCG2_RS2231142_HOM flags the TT group). ",
        "Jiang 2023 reports rs2231142 alleles on the G>T strand; dbSNP's ",
        "reference orientation for the same Q141K coding variant is C>A, ",
        "so GG = 421C/C, GT = 421C/A, TT = 421A/A. Jiang 2023 Table 2 ",
        "modeling cohort (n = 85): GG 34 (40%), GT 45 (53%), TT 6 (7%); ",
        "allele frequency G 0.66 / T 0.34; Hardy-Weinberg chi-square 2.99, ",
        "p = 0.22; not in linkage disequilibrium with the other five ",
        "genotyped SNPs (Supplementary Figure 1). The GG wild-type ",
        "reference (both paired indicators = 0) has a fixed CL/F factor of ",
        "1 (Jiang 2023 Table 3, 'rs2231142 GG: 1 (Fixed)')."
      ),
      source_name        = "rs2231142 GT"
    ),
    SNP_ABCG2_RS2231142_HOM = list(
      description        = "ABCG2 (BCRP) rs2231142 homozygous-variant indicator (TT genotype)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ABCG2 rs2231142 GG wild-type homozygote, when paired with SNP_ABCG2_RS2231142_HET = 0)",
      notes              = paste0(
        "Time-fixed per subject (germline genotype). 1 = subject carries ",
        "two rs2231142 variant alleles (genotype TT on the reported ",
        "strand, equivalently 421A/A); 0 = otherwise (the union of GG ",
        "wild-type homozygotes and GT heterozygotes; the paired indicator ",
        "SNP_ABCG2_RS2231142_HET flags the GT group). Jiang 2023 Table 2 ",
        "modeling cohort: TT 6/85 (7%). The column definition is identical ",
        "to the recessive-model use in Ueshima_2018_apixaban.R (1 = ",
        "421A/A); Jiang 2023 differs only in pairing it with a separate ",
        "heterozygote indicator so all three genotype strata carry ",
        "distinct typical-value factors, following the ",
        "CYP3A5_STAR1_HET / CYP3A5_STAR1_HOM and SLCO1B1_HAP15_HET / ",
        "SLCO1B1_HAP15_HOM precedents in ",
        "inst/references/covariate-columns.md. The TT factor (0.976, RSE ",
        "9%) is poorly identified despite the small reported RSE: only 6 ",
        "of 85 subjects were TT, and the bootstrap 95% CI (0.801-1.582) ",
        "is wide and spans 1, so the TT stratum is not distinguishable ",
        "from the GG reference."
      ),
      source_name        = "rs2231142 TT"
    )
  )

  # Covariates Jiang 2023 screened but did NOT retain in the final model.
  # Documentation only: `checkModelConventions()` does not require these to
  # appear in `model()`. Entries whose name is not a registered canonical in
  # inst/references/covariate-columns.md keep the paper's own abbreviation and
  # say so in `notes`; no register entry is created for a covariate that no
  # model uses.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste0(
        "Jiang 2023 Table 1: 58.6 +/- 9.6 kg (range 40.0-82.5). Carried ",
        "into the stepwise covariate search on CL/F and V/F but not ",
        "significant at the forward-inclusion threshold (dOFV < 3.84); ",
        "absent from Table 4 and from the final-model equation."
      )
    ),
    WBC = list(
      description = "White blood cell count",
      units       = "10^9 cells/L",
      type        = "continuous",
      notes       = paste0(
        "Jiang 2023 Table 1: median 4.0 (range 1.4-20.9) x 10^9/L. ",
        "Screened and not retained. Jiang 2023 Discussion explicitly ",
        "revisits this: a prior CML study found a WBC effect on imatinib ",
        "clearance, but 'no such correlation was found in the population ",
        "described in this study, which is consistent with the results ",
        "obtained by Delbaldo et al and Shriyan et al'."
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Jiang 2023 Table 1: median 44.0 g/L (range 30.6-54.8). Screened and not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Jiang 2023 Table 1: median 16.9 U/L (range 2.0-645.0). Screened and not retained."
    ),
    CRCL = list(
      description = "Creatinine clearance rate",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste0(
        "Jiang 2023 Table 1: median 69.7 mL/min (range 28.5-118.8). ",
        "Significant in UNIVARIATE forward inclusion on CL/F (dOFV = ",
        "-6.56, p < 0.05; Table 4) but dropped in the multifactorial step ",
        "(dOFV < 3.84) and therefore never reached backward elimination. ",
        "Jiang 2023 Discussion: imatinib clearance was not correlated with ",
        "creatinine clearance, 'as imatinib and its metabolites are barely ",
        "excreted through the kidneys'."
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste0(
        "Jiang 2023 Table 1: median 57 years (range 27-79). Excluded ",
        "BEFORE the stepwise process by the collinearity screen: the ",
        "correlation coefficient between age and CRCL was -0.67 ",
        "(Supplementary Figure 2), exceeding the |r| >= 0.5 threshold, so ",
        "'age appeared only in the initially considered covariates and was ",
        "excluded in the subsequent analysis' (Jiang 2023 Discussion)."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Jiang 2023 modeling cohort: 46 male / 39 female (45.9% female). ",
        "Listed under the investigated categorical covariates (Jiang 2023 ",
        "Covariate Model) but not among the covariates 'introduced into ",
        "the modeling process' in Results/Covariate Modeling, and absent ",
        "from Table 4 and the final model."
      )
    ),
    HGB = list(
      description = "Hemoglobin",
      units       = "g/L",
      type        = "continuous",
      notes       = paste0(
        "Jiang 2023 Table 1: 115.9 +/- 15.8 g/L (range 71.0-156.0). ",
        "Investigated, but not among the covariates carried into the ",
        "stepwise process -- HGB is collinear with the retained RBC ",
        "covariate and was removed by the |r| >= 0.5 correlation screen ",
        "(Supplementary Figure 2)."
      )
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste0(
        "Jiang 2023 Table 1: median 23.9 U/L (range 11.0-724.8). ",
        "Investigated, but not among the covariates carried into the ",
        "stepwise process (collinear with ALT under the |r| >= 0.5 ",
        "correlation screen, Supplementary Figure 2)."
      )
    ),
    PLT = list(
      description = "Platelet count (paper abbreviation; no canonical register entry -- not used by any model)",
      units       = "10^9 cells/L",
      type        = "continuous",
      notes       = paste0(
        "Jiang 2023 Table 1: median 178.0 (range 54.0-512.0) x 10^9/L. ",
        "Carried into the stepwise covariate search and not retained. No ",
        "canonical name is registered in ",
        "inst/references/covariate-columns.md because no nlmixr2lib model ",
        "uses a platelet-count covariate column; the paper's abbreviation ",
        "is kept here for provenance only."
      )
    ),
    GLB = list(
      description = "Serum globulin (paper abbreviation; no canonical register entry -- not used by any model)",
      units       = "g/L",
      type        = "continuous",
      notes       = paste0(
        "Jiang 2023 Table 1: 23.9 +/- 4.2 g/L (range 12.6-37.8). Carried ",
        "into the stepwise covariate search and not retained. Distinct ",
        "from the registered TPRO (total protein) canonical -- globulin is ",
        "total protein minus albumin, not total protein -- so no canonical ",
        "name is asserted here."
      )
    ),
    SNP_CYP1A2_RS11636419 = list(
      description = "CYP1A2 rs11636419 genotype (paper abbreviation; no canonical register entry -- not used by any model)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste0(
        "Jiang 2023 Table 2: GG 6 (0.07), GA 24 (0.28), AA 55 (0.65); ",
        "allele frequency G 0.21 / A 0.79; HWE chi-square 2.02, p = 0.36. ",
        "Screened on CL/F and V/F and not retained (absent from Table 4)."
      )
    ),
    SNP_SLC22A1_RS1867351 = list(
      description = "SLC22A1 (OCT1) rs1867351 genotype (paper abbreviation; no canonical register entry -- not used by any model)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste0(
        "Jiang 2023 Table 2: CC 13 (0.15), CT 40 (0.47), TT 32 (0.38); ",
        "allele frequency C 0.39 / T 0.61; HWE chi-square 0.01, p = 1.00. ",
        "Screened and not retained."
      )
    ),
    SNP_SLC22A1_RS683369 = list(
      description = "SLC22A1 (OCT1) rs683369 genotype (paper abbreviation; no canonical register entry -- not used by any model)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste0(
        "Jiang 2023 Table 2: CC 65 (0.76), CG 19 (0.22), GG 1 (0.01); ",
        "allele frequency C 0.88 / G 0.12; HWE chi-square 0.09, p = 0.96. ",
        "Screened and not retained."
      )
    ),
    SNP_SLC22A1_RS2282143 = list(
      description = "SLC22A1 (OCT1) rs2282143 genotype (paper abbreviation; no canonical register entry -- not used by any model)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste0(
        "Jiang 2023 Table 2: CC 62 (0.73), CT 22 (0.26), TT 1 (0.01); ",
        "allele frequency C 0.86 / T 0.14; HWE chi-square 0.39, p = 0.82. ",
        "Screened and not retained. Note that Jiang 2023 Table 2 typesets ",
        "rs2282143 under the SLCO1A2 gene row, while the Genotyping ",
        "Analysis methods text assigns rs2282143 to SLC22A1 and rs10841803 ",
        "to SLCO1A2; the methods-text assignment is used here."
      )
    ),
    SNP_SLCO1A2_RS10841803 = list(
      description = "SLCO1A2 (OATP1A2) rs10841803 genotype (paper abbreviation; no canonical register entry -- not used by any model)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = paste0(
        "Jiang 2023 Table 2: GG 29 (0.34), GA 43 (0.51), AA 13 (0.15); ",
        "allele frequency G 0.59 / A 0.41; HWE chi-square 0.20, p = 0.90. ",
        "Screened and not retained."
      )
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 85L,
    n_studies        = 1L,
    n_subjects_external_validation = 25L,
    n_subjects_enrolled = 110L,
    age_range        = "27-79 years (Table 1 median 57)",
    age_median       = "57 years",
    weight_range     = "40.0-82.5 kg",
    weight_mean      = "58.6 +/- 9.6 kg (mean +/- SD, Table 1)",
    sex_female_pct   = 45.9,
    race_ethnicity   = "Chinese (single-centre cohort, most patients recruited from Jiangxi Province)",
    disease_state    = paste0(
      "Adults (age >= 18 years at surgery) with histologically confirmed ",
      "gastrointestinal stromal tumor (GIST) who had undergone surgical ",
      "resection and were at intermediate or high risk of recurrence by ",
      "the modified National Institutes of Health criteria, receiving ",
      "adjuvant imatinib. Approximately half the cohort had comorbidities ",
      "such as hypertension and hyperglycemia (Jiang 2023 Discussion). ",
      "GIST KIT / PDGFRA mutation status was NOT determined (Jiang 2023 ",
      "Limitations)."
    ),
    dose_range       = paste0(
      "Oral imatinib 200-600 mg once daily. Initial-dose distribution in ",
      "the modeling dataset (200/300/400/500/600 mg): 1/2/79/2/1 patients ",
      "-- 92.9% received the standard 400 mg daily (Table 1)."
    ),
    regions          = "China (The First Affiliated Hospital of Nanchang University, Nanchang, Jiangxi Province)",
    enrollment_period = "March 2021 to June 2022",
    sampling_design  = paste0(
      "Sparse therapeutic drug monitoring during routine care. Most blood ",
      "samples were drawn at steady state (>= 30 days of continuous ",
      "imatinib). Only trough concentrations were collected, which is why ",
      "inter-individual variability on V/F was not estimable (Jiang 2023 ",
      "Results/PPK Analysis: 'the estimate of inter-individual variation ",
      "in the distribution volume was <1%'). Assay: two-dimensional liquid ",
      "chromatography; calibration range 13.44-1400.00 ng/mL; LLOQ 6.72 ",
      "ng/mL; inter-/intra-day RSD 0.3-4.0% and 1.1-3.2%."
    ),
    genotype_rs2231142 = c(GG = 40.0, GT = 53.0, TT = 7.0),
    lab_baseline     = paste0(
      "Table 1 modeling dataset: RBC 3.7 +/- 0.6 x 10^12/L (2.2-5.6); WBC ",
      "median 4.0 x 10^9/L (1.4-20.9); platelets median 178.0 x 10^9/L ",
      "(54.0-512.0); hemoglobin 115.9 +/- 15.8 g/L (71.0-156.0); albumin ",
      "median 44.0 g/L (30.6-54.8); globulin 23.9 +/- 4.2 g/L (12.6-37.8); ",
      "ALT median 16.9 U/L (2.0-645.0); AST median 23.9 U/L (11.0-724.8); ",
      "creatinine clearance median 69.7 mL/min (28.5-118.8)."
    ),
    notes            = paste0(
      "Software: NONMEM 7.3 (FOCE with interaction); R 4.0.5 for output ",
      "and diagnostics. Structural model: NONMEM ADVAN2 TRANS2 ",
      "(one-compartment, first-order absorption, CL/V parameterisation). ",
      "Model evaluation: 1000-sample bootstrap (990/1000 successful, ",
      "99.0%); prediction- and variability-corrected VPC (n = 1000); ",
      "external validation on a separate 25-patient dataset giving MPE% = ",
      "-11.9% and MAPE% = 17.9%. Monte Carlo dose simulations (n = 1000) ",
      "targeted a steady-state trough threshold of 1100 ng/mL. This is the ",
      "first imatinib popPK model in postoperative Chinese GIST adults to ",
      "incorporate genetic polymorphisms (Jiang 2023 Discussion)."
    )
  )

  ini({
    # ----- Structural parameters (Jiang 2023 Table 3, final model) -----
    # ka was NOT estimated: Jiang 2023 Base Model states 'the absorption
    # rate constant (ka) was fixed to 1.22 h-1' (citing their reference 26),
    # and Table 3 reports 'Ka (1/h): 1.22 (Fixed)' with no RSE and no
    # bootstrap entry.
    lka <- fixed(log(1.22)); label("First-order absorption rate constant ka (1/h)")  # Jiang 2023 Table 3, Ka = 1.22 (Fixed)

    # Typical apparent oral clearance for the REFERENCE subject: ABCG2
    # rs2231142 GG wild-type homozygote with RBC = 3.7 x 10^12/L. Both
    # covariate factors equal 1 at the reference.
    lcl <- log(9.72); label("Apparent oral clearance CL/F at the reference subject (L/h)")  # Jiang 2023 Table 3, CL/F = 9.72 L/h (RSE 8%; bootstrap median 9.72, 95% CI 8.22-11.84)

    # Apparent volume of distribution. No IIV is estimated on V/F (see the
    # IIV block below).
    lvc <- log(229); label("Apparent volume of distribution V/F (L)")  # Jiang 2023 Table 3, V/F = 229 L (RSE 21%; bootstrap median 226.3, 95% CI 156.4-421.6)

    # ----- Red blood cell count effect on CL/F (Jiang 2023 Table 3) -----
    # Power function (Jiang 2023 Equation 7: P = TV(P) * (Cov/Ref)^theta)
    # with Ref = 3.7 x 10^12/L, per the final-model equation printed on
    # page 814: CL/F = 9.72 * (RBC/3.7)^0.49 * ...
    # Cross-check on the exponent: Jiang 2023 Discussion reports that a 50%
    # reduction in RBC gives a 29% drop in CL/F; 0.5^0.49 = 0.712, i.e. a
    # 28.8% drop. Consistent.
    e_rbc_cl <- 0.49; label("Power exponent of (RBC / 3.7) on CL/F (unitless)")  # Jiang 2023 Table 3, RBC = 0.49 (RSE 38%; bootstrap median 0.48, 95% CI 0.10-0.89)

    # ----- ABCG2 rs2231142 genotype effects on CL/F (Jiang 2023 Table 3) --
    # Multiplicative factors relative to the GG wild-type reference, whose
    # factor is fixed at 1 (Table 3: 'rs2231142 GG: 1 (Fixed)'). Applied in
    # the power-of-binary-indicator form of Jiang 2023 Equation 8:
    #   P = TV(P) * theta1^heterozygous * theta2^homozygous
    # with the paired indicators both zero for GG.
    e_abcg2_het_cl <- 0.879; label("CL/F multiplicative factor for ABCG2 rs2231142 GT heterozygote (unitless)")       # Jiang 2023 Table 3, rs2231142 GT = 0.879 (RSE 5%; bootstrap median 0.877, 95% CI 0.785-0.968)
    e_abcg2_hom_cl <- 0.976; label("CL/F multiplicative factor for ABCG2 rs2231142 TT homozygous variant (unitless)")  # Jiang 2023 Table 3, rs2231142 TT = 0.976 (RSE 9%; bootstrap median 0.976, 95% CI 0.801-1.582)

    # ----- Inter-individual variability (Jiang 2023 Table 3) -----
    # Exponential IIV on CL/F only (Jiang 2023 Equation 1:
    # P_i = TV(P) * exp(eta_i), eta ~ N(0, omega^2)).
    # Table 3 reports the IIV for CL/F as 0.139, which the Results text
    # renders as a percentage: 'The inter-individual variation on CL/F and
    # the proportional residual error was found to be 13.9% and 29.6%'.
    # The tabulated 0.139 is therefore omega (the SD on the log scale,
    # approximately the CV), not omega^2. The final-model equation on page
    # 814 prints the variance directly as the exp() term
    # '... x e^0.0192' (= 0.139^2 = 0.0193 to rounding), and nlmixr2 takes
    # the variance here, so 0.0192 is used.
    # NO IIV on V/F: Jiang 2023 Results/PPK Analysis states 'Inter-individual
    # variation on V/F was not estimated in this analysis because only
    # trough concentrations were collected and the estimate of
    # inter-individual variation in the distribution volume was <1%.'
    etalcl ~ 0.0192  # Jiang 2023 final-model equation p. 814 (exp term e^0.0192); Table 3 IIV CL/F omega = 0.139 (RSE 39%; bootstrap median 0.133, 95% CI 0.061-0.197)

    # ----- Residual unexplained variability (Jiang 2023 Table 3) -----
    # Proportional error model, Jiang 2023 Equation 3: Y = F * (1 + eps1).
    # Additive (Eq. 2) and combined (Eq. 4) models were also evaluated; the
    # proportional model was selected (Results/PPK Analysis).
    propSd <- 0.296; label("Proportional residual error (fraction)")  # Jiang 2023 Table 3, Residual Proportional = 0.296 (RSE 27%; bootstrap median 0.288, 95% CI 0.219-0.371)
  })

  model({
    # ----- 1. Covariate factors on CL/F -----
    # Power-form RBC effect, referenced to the cohort central value
    # 3.7 x 10^12/L (Jiang 2023 final-model equation, p. 814).
    rbc_factor <- (RBC / 3.7)^e_rbc_cl

    # Three-level ABCG2 rs2231142 genotype factor in power-of-binary form.
    # The paired indicators are mutually exclusive, so at most one exponent
    # is 1; for the GG wild-type reference both are 0 and the factor
    # collapses to exactly 1 (matching Table 3's fixed GG value of 1).
    abcg2_factor <-
      e_abcg2_het_cl^SNP_ABCG2_RS2231142_HET *
      e_abcg2_hom_cl^SNP_ABCG2_RS2231142_HOM

    # ----- 2. Individual parameters -----
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * rbc_factor * abcg2_factor
    vc <- exp(lvc)

    # ----- 3. Micro-constants -----
    kel <- cl / vc

    # ----- 4. ODE system (NONMEM ADVAN2 TRANS2 equivalent) -----
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ----- 5. Observation and error -----
    # `central` is an amount in mg and `vc` is in L, so central/vc is mg/L;
    # the factor 1000 converts to the ng/mL units the paper reports
    # concentrations and the 1100 ng/mL target trough in.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
