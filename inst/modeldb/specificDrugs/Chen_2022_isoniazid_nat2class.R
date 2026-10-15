Chen_2022_isoniazid_nat2class <- function() {
  description <- paste(
    "Integrated parent-metabolite population PK model for oral isoniazid",
    "(INH) and its metabolite acetylisoniazid (AcINH) in 45 healthy Chinese",
    "adults and 157 Chinese adults with tuberculosis (Chen 2022, Final",
    "Model (1)). One-compartment INH with first-order absorption; a",
    "fraction FM of INH clearance forms AcINH, which has its own",
    "one-compartment disposition with first-order elimination rate K30.",
    "A three-level NAT2 genotype class score (0 = wt/wt, 1 = m/wt,",
    "2 = m/m; derived from NAT2_RAPID and NAT2_SLOW) enters exponentially",
    "on INH CL/F and on FM. Concentrations are in umol/L (doses in mg are",
    "converted with the INH molecular weight). Alternative final model to",
    "Chen_2022_isoniazid (per-allele NAT2 covariates)."
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
    NAT2_RAPID = list(
      description = paste(
        "NAT2 wild-type (*4/*4) indicator: 1 = no *5, *6 or *7 allele",
        "(the paper's wt/wt group), 0 = otherwise."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0; the paper's reference (score 0) is NAT2_RAPID = 1.",
      notes = paste(
        "Chen 2022 Methods 'Covariates' method 2: wt/wt, m/wt and m/m",
        "scored 0, 1 and 2 (Eq. 7). Only the *5, *6 and *7 SNPs were typed,",
        "so wt/wt = *4/*4. Score = (1 - NAT2_RAPID) * (1 + NAT2_SLOW)."
      ),
      source_name = "Genotype (wt/wt)"
    ),
    NAT2_SLOW = list(
      description = paste(
        "NAT2 homozygous-mutant indicator: 1 = two of the *5, *6 or *7",
        "alleles (the paper's m/m group), 0 = otherwise."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wt/wt or m/wt)",
      notes = paste(
        "Chen 2022 Methods 'Covariates' method 2 m/m group (score 2):",
        "*5/*7, *6/*6, *6/*7 and *7/*7 in this cohort. The joint state",
        "NAT2_RAPID = 0, NAT2_SLOW = 0 is the m/wt group (score 1)."
      ),
      source_name = "Genotype (m/m)"
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
    nat2_class_counts = c("wt/wt" = 91L, "m/wt" = 79L, "m/m" = 32L),
    notes = paste(
      "Table 1 demographics; NAT2 class counts summed from Results",
      "'INH and AcINH in Relation to NAT2 Genotypes' (study 1: 8/8/8;",
      "study 2: 11/9/1; patients: 72/62/23). 122 subjects formed the index",
      "group and 80 the validation group; the final model was fit to both."
    )
  )

  ini({
    # Chen 2022 Table 3, 'Final Model (1)' column (NONMEM 6, FOCE).
    lka <- log(3.91); label("INH first-order absorption rate constant Ka (1/h)") # Table 3 Final Model (1) theta1 Ka = 3.91 (SE 0.44)
    lcl <- log(28.7); label("INH apparent clearance CL/F for NAT2 wt/wt (L/h)") # Table 3 Final Model (1) theta2 CL/F = 28.7 (SE 3.22)
    lkel_acinh <- log(0.41); label("AcINH elimination rate constant K30 (1/h)") # Table 3 Final Model (1) theta3 K30 = 0.41 (SE 0.051)
    lvc <- log(54.1); label("INH apparent volume of distribution V2/F (L)") # Table 3 Final Model (1) theta4 V2/F = 54.1 (SE 12.5)
    lvc_acinh <- log(17.2); label("AcINH apparent volume of distribution V3/F (L)") # Table 3 Final Model (1) theta5 V3/F = 17.2 (SE 3.21)
    lfm <- log(0.88); label("Fraction of INH clearance forming AcINH, FM, for NAT2 wt/wt (fraction)") # Table 3 Final Model (1) theta6 FM = 0.88 (SE 0.21)

    # P = theta * exp(theta_g * score), score 0/1/2 (Eq. 7).
    e_nat2_cl <- -0.55; label("Exponent of the NAT2 class score on CL/F (per score unit)") # Table 3 Final Model (1) theta7 = -0.55 (SE 0.11); footnote CL/F = 28.7*exp(-0.55*Genotype)
    # Table 3 lists theta10 = -0.47 (SE 0.17); its footnote prints
    # FM = 0.88*exp(-0.55*Genotype), repeating the CL exponent. The table
    # value is used: it is a separately estimated theta with its own SE, and
    # it reproduces the observed AcINH AUC of the m/wt and m/m groups
    # (Table 2) more closely than -0.55 does (see the vignette).
    e_nat2_fm <- -0.47; label("Exponent of the NAT2 class score on FM (per score unit)") # Table 3 Final Model (1) theta10 = -0.47 (SE 0.17)

    mw_inh <- fixed(137.14); label("Isoniazid molecular weight (g/mol)") # C6H7N3O; consistent with the paper's unit-conversion pairs (15.89 mg/L = 115.9 umol/L)

    # IIV: exponential model (Eq. 1); omega reported as omega x 100.
    etalka ~ 0.304704 # Table 3 Final Model (1) omega Ka = 55.2 percent; 0.552^2
    etalcl ~ 0.094249 # Table 3 Final Model (1) omega CL/F = 30.7 percent; 0.307^2
    etalkel_acinh ~ 0.045796 # Table 3 Final Model (1) omega K30 = 21.4 percent; 0.214^2
    etalvc ~ 0.037636 # Table 3 Final Model (1) omega V2/F = 19.4 percent; 0.194^2
    etalfm ~ 0.011449 # Table 3 Final Model (1) omega FM = 10.7 percent; 0.107^2

    propSd <- 0.333; label("INH proportional residual error (fraction)") # Table 3 Final Model (1) sigma INH = 33.3 percent
    propSd_acinh <- 0.302; label("AcINH proportional residual error (fraction)") # Table 3 Final Model (1) sigma AcINH = 30.2 percent
  })

  model({
    # wt/wt = 0, m/wt = 1, m/m = 2 (Methods 'Covariates' method 2).
    nat2_score <- (1 - NAT2_RAPID) * (1 + NAT2_SLOW)

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * exp(e_nat2_cl * nat2_score)
    vc <- exp(lvc + etalvc)
    kel_acinh <- exp(lkel_acinh + etalkel_acinh)
    vc_acinh <- exp(lvc_acinh)
    # Exponential IIV as printed (Eq. 1); FM is not bounded above by 1.
    fm <- exp(lfm + etalfm) * exp(e_nat2_fm * nat2_score)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(central_acinh) <- fm * kel * central * 1000 / mw_inh - kel_acinh * central_acinh

    Cc <- central / vc * 1000 / mw_inh
    Cc_acinh <- central_acinh / vc_acinh

    Cc ~ prop(propSd)
    Cc_acinh ~ prop(propSd_acinh)
  })
}
