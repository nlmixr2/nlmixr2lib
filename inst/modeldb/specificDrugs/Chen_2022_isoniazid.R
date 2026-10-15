Chen_2022_isoniazid <- function() {
  description <- paste(
    "Integrated parent-metabolite population PK model for oral isoniazid",
    "(INH) and its metabolite acetylisoniazid (AcINH) in 45 healthy Chinese",
    "adults and 157 Chinese adults with tuberculosis (Chen 2022, Final",
    "Model (2)). One-compartment INH with first-order absorption; a",
    "fraction FM of INH clearance forms AcINH, which has its own",
    "one-compartment disposition with first-order elimination rate K30.",
    "The numbers of NAT2 *5 (341T>C), *6 (590G>A) and *7 (857G>A) alleles",
    "each enter exponentially on INH CL/F and on FM. Concentrations are",
    "in umol/L (doses in mg are converted with the INH molecular weight).",
    "Known deviation: the printed AcINH parameters (K30 x V3/F = 3.7 L/h)",
    "give a typical AcINH exposure about twice the paper's own observed",
    "NCA and VPC in NAT2 *4/*4 subjects; INH reproduces. See the vignette."
  )
  reference <- paste(
    "Chen B, Shi H-Q, Feng MR, Wang X-H, Cao X-M, Cai W-M (2022).",
    "Population Pharmacokinetics and Pharmacodynamics of Isoniazid and its",
    "Metabolite Acetylisoniazid in Chinese Population.",
    "Front Pharmacol 13:932686. doi:10.3389/fphar.2022.932686."
  )
  vignette <- "Chen_2022_isoniazid"
  units <- list(time = "h", dosing = "mg", concentration = "umol/L")

  covariateData <- list(
    SNP_NAT2_RS1801280_C_COUNT = list(
      description = "Number of NAT2 *5 alleles (341T>C, rs1801280 variant C allele): 0, 1 or 2.",
      units = "(count, 0/1/2 alleles per subject)",
      type = "continuous",
      reference_category = "0 (no *5 allele)",
      notes = paste(
        "Chen 2022 Methods 'Covariates' method 3 scores each of the *5, *6",
        "and *7 alleles as 0 (w/w), 1 (m/w) or 2 (m/m); Table 3 row 'M341'.",
        "The paper identifies *5 by the 341 SNP alone (allele-specific PCR),",
        "so the *5 allele count equals the rs1801280 variant-allele count.",
        "The paper writes the SNP as 'C341 -> T'; the NAT2*5 defining change",
        "is 341T>C (variant C). Enters as exp(-0.77 * count) on CL/F and",
        "exp(-0.72 * count) on FM."
      ),
      source_name = "M341"
    ),
    SNP_NAT2_RS1799930_A_COUNT = list(
      description = "Number of NAT2 *6 alleles (590G>A, rs1799930 variant A allele): 0, 1 or 2.",
      units = "(count, 0/1/2 alleles per subject)",
      type = "continuous",
      reference_category = "0 (no *6 allele)",
      notes = paste(
        "Chen 2022 Methods 'Covariates' method 3; Table 3 row 'M590'.",
        "Enters as exp(-0.60 * count) on CL/F and exp(-0.45 * count) on FM."
      ),
      source_name = "M590"
    ),
    SNP_NAT2_RS1799931_A_COUNT = list(
      description = "Number of NAT2 *7 alleles (857G>A, rs1799931 variant A allele): 0, 1 or 2.",
      units = "(count, 0/1/2 alleles per subject)",
      type = "continuous",
      reference_category = "0 (no *7 allele)",
      notes = paste(
        "Chen 2022 Methods 'Covariates' method 3; Table 3 row 'M870' and",
        "the Table 3 footnote 'M803' are both misprints of the 857 position",
        "that defines NAT2*7 (Methods 'Genotyping': 'G857 -> A').",
        "Enters as exp(-0.29 * count) on CL/F and exp(-0.14 * count) on FM."
      ),
      source_name = "M870"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "isoniazid", units = "mg", specimen = "plasma", verified = TRUE),
    central_acinh = list(analyte = "acetylisoniazid", units = "umol", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 202L,
    n_studies = 3L,
    age_range = "19-64 years (healthy 21-29 years; patients 19-64 years)",
    weight_range = "39-78 kg",
    sex_female_pct = 33.7,
    race_ethnicity = "Chinese (healthy subjects all Han)",
    disease_state = paste(
      "45 healthy adult male volunteers (two single-dose studies) and 157",
      "adults with pulmonary tuberculosis on 7-14 days of combination",
      "therapy (isoniazid with rifampicin or rifapentine, pyrazinamide and",
      "ethambutol)."
    ),
    dose_range = paste(
      "Healthy: single oral 300 mg (study 1, n = 24) or 320 mg (study 2,",
      "bioequivalence, n = 21). Patients: daily oral isoniazid, sampled",
      "2 and/or 6 h post-dose."
    ),
    regions = "China",
    nat2_genotype_counts = c(
      "*4/*4" = 91L,
      "*4/*5" = 7L,
      "*4/*6" = 43L,
      "*4/*7" = 29L,
      "*5/*7" = 3L,
      "*6/*6" = 20L,
      "*6/*7" = 5L,
      "*7/*7" = 4L
    ),
    notes = paste(
      "Table 1 demographics; NAT2 genotype counts summed from Results",
      "'INH and AcINH in Relation to NAT2 Genotypes' over the three cohorts",
      "(study 1: 8/0/6/2/1/7/0/0; study 2: 11/1/6/2/0/1/0/0; patients:",
      "72/6/31/25/2/12/5/4 in the order listed). 122 subjects formed the",
      "index group and 80 the validation group; the final model was fit",
      "to both. Female percentage is 68/202 (all healthy volunteers male)."
    )
  )

  ini({
    # Chen 2022 Table 3, 'Final Model (2)' column (NONMEM 6, FOCE).
    # The model was fit to concentrations in umol/L (Figures 2-4 axes).
    lka <- log(3.99); label("INH first-order absorption rate constant Ka (1/h)") # Table 3 Final Model (2) theta1 Ka = 3.99 (SE 0.42)
    lcl <- log(30.2); label("INH apparent clearance CL/F for NAT2 *4/*4 (L/h)") # Table 3 Final Model (2) theta2 CL/F = 30.2 (SE 3.04)
    lkel_acinh <- log(0.30); label("AcINH elimination rate constant K30 (1/h)") # Table 3 Final Model (2) theta3 K30 = 0.30 (SE 0.024)
    lvc <- log(51.4); label("INH apparent volume of distribution V2/F (L)") # Table 3 Final Model (2) theta4 V2/F = 51.4 (SE 3.28)
    lvc_acinh <- log(12.4); label("AcINH apparent volume of distribution V3/F (L)") # Table 3 Final Model (2) theta5 V3/F = 12.4 (SE 0.41)
    lfm <- log(0.86); label("Fraction of INH clearance forming AcINH, FM, for NAT2 *4/*4 (fraction)") # Table 3 Final Model (2) theta6 FM = 0.86 (SE 0.34)

    # Covariate exponents: P = theta * exp(sum(theta_k * allele count)),
    # Table 3 footnote 'Final Model (2)'. Back-transforms match Results:
    # one copy of *5/*6/*7 gives CL/F 46.3/54.9/74.8 and FM 48.7/63.8/86.9
    # of the *4/*4 value.
    e_snp_nat2_rs1801280_cl <- -0.77; label("Exponent of NAT2 *5 allele count on CL/F (per allele)") # Table 3 theta7 M341 (CL) = -0.77 (SE 0.119)
    e_snp_nat2_rs1799930_cl <- -0.60; label("Exponent of NAT2 *6 allele count on CL/F (per allele)") # Table 3 theta8 M590 (CL) = -0.60 (SE 0.066)
    e_snp_nat2_rs1799931_cl <- -0.29; label("Exponent of NAT2 *7 allele count on CL/F (per allele)") # Table 3 theta9 M870 (CL) = -0.29 (SE 0.095)
    e_snp_nat2_rs1801280_fm <- -0.72; label("Exponent of NAT2 *5 allele count on FM (per allele)") # Table 3 theta10 M341 (FM) = -0.72 (SE 0.040)
    e_snp_nat2_rs1799930_fm <- -0.45; label("Exponent of NAT2 *6 allele count on FM (per allele)") # Table 3 theta11 M590 (FM) = -0.45 (SE 0.046)
    e_snp_nat2_rs1799931_fm <- -0.14; label("Exponent of NAT2 *7 allele count on FM (per allele)") # Table 3 theta12 M870 (FM) = -0.14 (SE 0.069)

    # Molecular weight used to express INH in umol; the paper's own
    # mg/L <-> umol/L pairs imply 137.1 (Methods: 15.89 mg/L = 115.9 umol/L;
    # Results: 19.7 ug*h/mL = 143.4 umol*h/L).
    mw_inh <- fixed(137.14); label("Isoniazid molecular weight (g/mol)") # C6H7N3O; consistent with the paper's unit-conversion pairs

    # IIV: exponential model (Eq. 1); Table 3 reports omega as a percentage,
    # i.e. omega x 100 with omega the SD on the log scale, so the variance
    # is (P/100)^2.
    etalka ~ 0.286225 # Table 3 Final Model (2) omega Ka = 53.5 percent; 0.535^2
    etalcl ~ 0.078961 # Table 3 Final Model (2) omega CL/F = 28.1 percent; 0.281^2
    etalkel_acinh ~ 0.048841 # Table 3 Final Model (2) omega K30 = 22.1 percent; 0.221^2
    etalvc ~ 0.039204 # Table 3 Final Model (2) omega V2/F = 19.8 percent; 0.198^2
    etalfm ~ 0.0096432 # Table 3 Final Model (2) omega FM = 9.82 percent; 0.0982^2

    # Proportional residual error (Eq. 2) for each analyte.
    propSd <- 0.339; label("INH proportional residual error (fraction)") # Table 3 Final Model (2) sigma INH = 33.9 percent
    propSd_acinh <- 0.295; label("AcINH proportional residual error (fraction)") # Table 3 Final Model (2) sigma AcINH = 29.5 percent
  })

  model({
    nat2_cl <- exp(
      e_snp_nat2_rs1801280_cl * SNP_NAT2_RS1801280_C_COUNT +
        e_snp_nat2_rs1799930_cl * SNP_NAT2_RS1799930_A_COUNT +
        e_snp_nat2_rs1799931_cl * SNP_NAT2_RS1799931_A_COUNT
    )
    nat2_fm <- exp(
      e_snp_nat2_rs1801280_fm * SNP_NAT2_RS1801280_C_COUNT +
        e_snp_nat2_rs1799930_fm * SNP_NAT2_RS1799930_A_COUNT +
        e_snp_nat2_rs1799931_fm * SNP_NAT2_RS1799931_A_COUNT
    )

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * nat2_cl
    vc <- exp(lvc + etalvc)
    kel_acinh <- exp(lkel_acinh + etalkel_acinh)
    vc_acinh <- exp(lvc_acinh)
    # Exponential IIV as printed (Eq. 1); FM is not bounded above by 1.
    fm <- exp(lfm + etalfm) * nat2_fm

    kel <- cl / vc

    # INH in mg; AcINH formed at FM * CL/F * C(INH) in umol (Figure 1).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(central_acinh) <- fm * kel * central * 1000 / mw_inh - kel_acinh * central_acinh

    Cc <- central / vc * 1000 / mw_inh
    Cc_acinh <- central_acinh / vc_acinh

    Cc ~ prop(propSd)
    Cc_acinh ~ prop(propSd_acinh)
  })
}
