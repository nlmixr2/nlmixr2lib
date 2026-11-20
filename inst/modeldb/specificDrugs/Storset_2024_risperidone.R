Storset_2024_risperidone <- function() {
  description <- "Joint parent-metabolite population PK model for oral risperidone and 9-hydroxyrisperidone in adult therapeutic-drug-monitoring patients genotyped for CYP2D6 (Storset 2024). One-compartment disposition for each analyte with first-order absorption (ka fixed at 2.01 1/h) and linear elimination; parameters are apparent (CL/F, V/F for the parent; CLmet/(F x fmet), Vmet/(F x fmet) for the metabolite) and all parent clearance is routed to the metabolite. Risperidone CL/F is the SUM of a CYP2D6-independent base clearance and one partial clearance per expressed CYP2D6 allele, each a fraction of the CYP2D6*1 partial clearance (*2-rs5758550G 0.66, *35 0.57, *2-rs5758550A 0.30, *9 0.39, *10 0.32, *41 0.15, deficient *3-*6 0). An NFIB rs28379954 C-allele carrier has 1.41-fold higher CYP2D6-mediated clearance. Both clearances decline linearly with age above an estimated breakpoint (risperidone -0.9 percent per year above 34 years; 9-hydroxyrisperidone -1.3 percent per year above 39 years). Correlated log-normal IIV on both clearances; combined proportional + additive residual error for risperidone and proportional for 9-hydroxyrisperidone."
  reference <- paste(
    "Storset E, Braten LS, Ingelman-Sundberg M, Johansson I, Molden E,",
    "Kringen MK (2024). Impact of CYP2D6*2, CYP2D6*35, rs5758550, and",
    "related haplotypes on risperidone clearance in vivo. European Journal",
    "of Clinical Pharmacology 80(10):1531-1541.",
    "doi:10.1007/s00228-024-03721-6.",
    sep = " "
  )
  vignette <- "Storset_2024_risperidone"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "risperidone", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "risperidone", units = "mg", specimen = "serum", verified = TRUE),
    central_9oh = list(analyte = "9-hydroxyrisperidone", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Piecewise-linear (hinge) effect on both clearances: factor 1 + e_age_cl * max(AGE - age_hinge, 0) with age_hinge = 34 years for risperidone and 39 years for 9-hydroxyrisperidone (Storset 2024 Table 2 footnote). Factor is 1 at or below the breakpoint. Cohort median 37 years (range 18-89; Table 1).",
      source_name = "Age"
    ),
    SNP_NFIB_RS28379954 = list(
      description = "NFIB rs28379954 (T>C) variant-carrier indicator: 1 = CC or CT, 0 = TT",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (NFIB-T homozygote)",
      notes = "Multiplies the CYP2D6-mediated partial clearance (not the base clearance) by 1.41 in carriers (Storset 2024 Table 2 'NFIB-C'; heterozygous and homozygous carriers merged, Methods 'Covariate model'). Six of 1024 alleles had undetermined NFIB genotype (Table 1); the final model has no separate unknown category, so undetermined genotype is coded 0.",
      source_name = "NFIB-C"
    ),
    CYP2D6_STAR1_HET = list(
      description = "CYP2D6*1 heterozygote indicator (exactly one *1 allele)",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = paste(
        "Copies of CYP2D6*1 = HET + 2 * HOM; each copy contributes 1.0 x the CYP2D6*1 partial clearance to risperidone CL/F (Storset 2024 Table 2 'CYP2D6 allele (fraction of CYP2D6*1)').",
        "Includes *1-rs5758550A and the three *1-rs5758550G alleles the paper merged with it (Results 'Subjects and data')."
      ),
      source_name = "CYP2D6*1"
    ),
    CYP2D6_STAR1_HOM = list(
      description = "CYP2D6*1 homozygote indicator (*1/*1)",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Copies of CYP2D6*1 = HET + 2 * HOM; each copy contributes 1.0 x the CYP2D6*1 partial clearance to risperidone CL/F (Storset 2024 Table 2 'CYP2D6 allele (fraction of CYP2D6*1)').",
      source_name = "CYP2D6*1"
    ),
    CYP2D6_STAR2_RS5758550G_HET = list(
      description = "CYP2D6*2-rs5758550G haplotype heterozygote indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Copies of CYP2D6*2-rs5758550G = HET + 2 * HOM; each copy contributes 0.66 x the CYP2D6*1 partial clearance to risperidone CL/F (Storset 2024 Table 2 'CYP2D6 allele (fraction of CYP2D6*1)').",
      source_name = "*2-rs5758550G (*2-G)"
    ),
    CYP2D6_STAR2_RS5758550G_HOM = list(
      description = "CYP2D6*2-rs5758550G haplotype homozygote indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Copies of CYP2D6*2-rs5758550G = HET + 2 * HOM; each copy contributes 0.66 x the CYP2D6*1 partial clearance to risperidone CL/F (Storset 2024 Table 2 'CYP2D6 allele (fraction of CYP2D6*1)').",
      source_name = "*2-rs5758550G (*2-G)"
    ),
    CYP2D6_STAR2_RS5758550A_HET = list(
      description = "CYP2D6*2-rs5758550A haplotype heterozygote indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Copies of CYP2D6*2-rs5758550A = HET + 2 * HOM; each copy contributes 0.30 x the CYP2D6*1 partial clearance to risperidone CL/F (Storset 2024 Table 2 'CYP2D6 allele (fraction of CYP2D6*1)').",
      source_name = "*2-rs5758550A (*2-A)"
    ),
    CYP2D6_STAR2_RS5758550A_HOM = list(
      description = "CYP2D6*2-rs5758550A haplotype homozygote indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Copies of CYP2D6*2-rs5758550A = HET + 2 * HOM; each copy contributes 0.30 x the CYP2D6*1 partial clearance to risperidone CL/F (Storset 2024 Table 2 'CYP2D6 allele (fraction of CYP2D6*1)').",
      source_name = "*2-rs5758550A (*2-A)"
    ),
    CYP2D6_STAR35_HET = list(
      description = "CYP2D6*35 heterozygote indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Copies of CYP2D6*35 = HET + 2 * HOM; each copy contributes 0.57 x the CYP2D6*1 partial clearance to risperidone CL/F (Storset 2024 Table 2 'CYP2D6 allele (fraction of CYP2D6*1)').",
      source_name = "*35"
    ),
    CYP2D6_STAR35_HOM = list(
      description = "CYP2D6*35 homozygote indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Copies of CYP2D6*35 = HET + 2 * HOM; each copy contributes 0.57 x the CYP2D6*1 partial clearance to risperidone CL/F (Storset 2024 Table 2 'CYP2D6 allele (fraction of CYP2D6*1)').",
      source_name = "*35"
    ),
    CYP2D6_STAR9_HET = list(
      description = "CYP2D6*9 heterozygote indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Copies of CYP2D6*9 = HET + 2 * HOM; each copy contributes 0.39 x the CYP2D6*1 partial clearance to risperidone CL/F (Storset 2024 Table 2 'CYP2D6 allele (fraction of CYP2D6*1)').",
      source_name = "*9"
    ),
    CYP2D6_STAR9_HOM = list(
      description = "CYP2D6*9 homozygote indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Copies of CYP2D6*9 = HET + 2 * HOM; each copy contributes 0.39 x the CYP2D6*1 partial clearance to risperidone CL/F (Storset 2024 Table 2 'CYP2D6 allele (fraction of CYP2D6*1)').",
      source_name = "*9"
    ),
    CYP2D6_STAR10_HET = list(
      description = "CYP2D6*10 heterozygote indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = paste(
        "Copies of CYP2D6*10 = HET + 2 * HOM; each copy contributes 0.32 x the CYP2D6*1 partial clearance to risperidone CL/F (Storset 2024 Table 2 'CYP2D6 allele (fraction of CYP2D6*1)').",
        "Here *10 is the star-allele call from the routine panel, in which *4 (which also carries rs1065852 T) is called separately as a deficient allele."
      ),
      source_name = "*10"
    ),
    CYP2D6_STAR10_HOM = list(
      description = "CYP2D6*10 homozygote indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Copies of CYP2D6*10 = HET + 2 * HOM; each copy contributes 0.32 x the CYP2D6*1 partial clearance to risperidone CL/F (Storset 2024 Table 2 'CYP2D6 allele (fraction of CYP2D6*1)').",
      source_name = "*10"
    ),
    CYP2D6_STAR41_HET = list(
      description = "CYP2D6*41 heterozygote indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Copies of CYP2D6*41 = HET + 2 * HOM; each copy contributes 0.15 x the CYP2D6*1 partial clearance to risperidone CL/F (Storset 2024 Table 2 'CYP2D6 allele (fraction of CYP2D6*1)').",
      source_name = "*41"
    ),
    CYP2D6_STAR41_HOM = list(
      description = "CYP2D6*41 homozygote indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Copies of CYP2D6*41 = HET + 2 * HOM; each copy contributes 0.15 x the CYP2D6*1 partial clearance to risperidone CL/F (Storset 2024 Table 2 'CYP2D6 allele (fraction of CYP2D6*1)').",
      source_name = "*41"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Tested on risperidone CL, 9-hydroxyrisperidone CL and both volumes; not retained (p >= 0.23; Storset 2024 Results 'Covariate model').",
      source_name = "Sex"
    ),
    SNP_CYP3A4_RS35599367 = list(
      description = "CYP3A4*22 (rs35599367) carrier indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-carrier)",
      notes = "Tested on the CYP2D6-independent base clearance; not retained (p = 0.31; Storset 2024 Results 'Covariate model').",
      source_name = "CYP3A4*22"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 512L,
    n_observations = "1565 paired risperidone + 9-hydroxyrisperidone serum concentrations (median 2 per subject, range 1-41)",
    n_studies = 1L,
    age_range = "18-89 years",
    age_median = "37 years",
    sex_female_pct = 46.3,
    disease_state = "Adult patients on oral risperidone with routine therapeutic drug monitoring (indications not recorded); CYP2D6-inhibiting comedication (paroxetine, fluoxetine, bupropion), CYP3A inducers and recent long-acting risperidone/paliperidone injections excluded.",
    dose_range = "Oral risperidone 0.5-12 mg/day (median 3 mg/day); once (789 samples), twice (758) or three times (18) daily",
    regions = "Norway (Center for Psychopharmacology, Diakonhjemmet Hospital, Oslo; 2005-2022)",
    sampling = "Trough-type TDM samples 9-30 h after the last dose (median 13 h), assumed at steady state",
    cyp2d6_alleles = "Of 1026 alleles: *1 401, *2-rs5758550G 166, *2-rs5758550A 27, *35 54, *9 20, *10 23, *41 66, deficient (*3, *4, *5, *6) 269 (Table 1); diplotype counts in Supplementary Table S1",
    notes = "No body size or organ function data were available (Discussion)."
  )

  ini({
    lka <- fixed(log(2.01)); label("Absorption rate constant ka (1/h)") # Table 2: ka = 2.01 h-1 (fix), from literature [31]
    lvc <- log(333); label("Apparent volume of distribution of risperidone V/F (L)") # Table 2: V/F = 333 L (bootstrap 95% CI 177, 586)

    lcl_base <- log(4.2); label("CYP2D6-independent base clearance of risperidone CLbase/F (L/h)") # Table 2: Clearance Base = 4.2 L/h (95% CI 3.6, 4.9)
    lcl_cyp2d6_star1 <- log(23.2); label("Partial risperidone clearance per CYP2D6*1 allele (L/h)") # Table 2: Clearance CYP2D6*1 = 23.2 L/h (95% CI 15.4, 31.6)

    f_cyp2d6_star2g <- 0.66; label("Partial clearance of CYP2D6*2-rs5758550G relative to *1 (fraction)") # Table 2: *2-rs5758550G = 0.66 (95% CI 0.50, 0.85)
    f_cyp2d6_star35 <- 0.57; label("Partial clearance of CYP2D6*35 relative to *1 (fraction)") # Table 2: *35 = 0.57 (95% CI 0.28, 0.90)
    f_cyp2d6_star2a <- 0.30; label("Partial clearance of CYP2D6*2-rs5758550A relative to *1 (fraction)") # Table 2: *2-rs5758550A = 0.30 (95% CI 0.17, 0.49)
    f_cyp2d6_star9 <- 0.39; label("Partial clearance of CYP2D6*9 relative to *1 (fraction)") # Table 2: *9 = 0.39 (95% CI 0.20, 0.68)
    f_cyp2d6_star10 <- 0.32; label("Partial clearance of CYP2D6*10 relative to *1 (fraction)") # Table 2: *10 = 0.32 (95% CI 0.15, 0.56)
    f_cyp2d6_star41 <- 0.15; label("Partial clearance of CYP2D6*41 relative to *1 (fraction)") # Table 2: *41 = 0.15 (95% CI 0.08, 0.26)

    e_snp_nfib_rs28379954_cl <- 1.41; label("Multiplier on CYP2D6-mediated clearance for NFIB-C carriers (ratio)") # Table 2: NFIB-C (fraction of NFIB-T) = 1.41 (95% CI 1.15, 1.73)

    lage_hinge <- log(34); label("Age breakpoint for risperidone clearance decline (years)") # Table 2: Age breaking point = 34 years (95% CI 28, 38)
    e_age_cl <- -0.009; label("Fractional change in risperidone CL/F per year above the age breakpoint (1/year)") # Table 2: Reduction per year = -0.9% (95% CI -1.2%, -0.5%)

    lvc_9oh <- log(96); label("Apparent volume of distribution of 9-hydroxyrisperidone Vmet/(F x fmet) (L)") # Table 2: Vmet/F x fmet = 96 (95% CI 77, 127)
    lcl_9oh <- log(8.0); label("Apparent clearance of 9-hydroxyrisperidone CLmet/(F x fmet) (L/h)") # Table 2: CLmet/F x fmet = 8.0 (95% CI 7.5, 8.6)
    lage_hinge_9oh <- log(39); label("Age breakpoint for 9-hydroxyrisperidone clearance decline (years)") # Table 2: Age breaking point = 39 years (95% CI 30, 48)
    e_age_cl_9oh <- -0.013; label("Fractional change in 9-hydroxyrisperidone CL per year above the age breakpoint (1/year)") # Table 2: Reduction per year = -1.3% (95% CI -1.6%, -1.0%)

    # omega^2 = log(1 + CV^2): 86% CV -> 0.5536; 49% CV -> 0.2152;
    # covariance = 0.43 * sqrt(0.5536 * 0.2152) = 0.1484
    etalcl + etalcl_9oh ~ c(
      0.5536,
      0.1484, 0.2152
    ) # Table 2: BSV CL/F 86% CV, BSV CLmet/fmet 49% CV, correlation 0.43

    propSd <- 0.58; label("Proportional residual error, risperidone (fraction)") # Table 2: Residual error, proportional = 58%
    addSd <- 0.06; label("Additive residual error, risperidone (ng/mL)") # Table 2: Residual error, additive = 0.06 ng/mL
    propSd_9oh <- 0.38; label("Proportional residual error, 9-hydroxyrisperidone (fraction)") # Table 2: 9-hydroxyrisperidone residual error, proportional = 38%
  })
  model({
    # Copies of each CYP2D6 allele; the deficient alleles *3-*6 contribute
    # zero clearance (fixed in the paper) and so need no column.
    n_star1 <- CYP2D6_STAR1_HET + 2 * CYP2D6_STAR1_HOM
    n_star2g <- CYP2D6_STAR2_RS5758550G_HET + 2 * CYP2D6_STAR2_RS5758550G_HOM
    n_star2a <- CYP2D6_STAR2_RS5758550A_HET + 2 * CYP2D6_STAR2_RS5758550A_HOM
    n_star35 <- CYP2D6_STAR35_HET + 2 * CYP2D6_STAR35_HOM
    n_star9 <- CYP2D6_STAR9_HET + 2 * CYP2D6_STAR9_HOM
    n_star10 <- CYP2D6_STAR10_HET + 2 * CYP2D6_STAR10_HOM
    n_star41 <- CYP2D6_STAR41_HET + 2 * CYP2D6_STAR41_HOM

    # Diplotype activity in units of one CYP2D6*1 allele
    cyp2d6_act <- n_star1 + f_cyp2d6_star2g * n_star2g +
      f_cyp2d6_star2a * n_star2a + f_cyp2d6_star35 * n_star35 +
      f_cyp2d6_star9 * n_star9 + f_cyp2d6_star10 * n_star10 +
      f_cyp2d6_star41 * n_star41

    # Table 2 footnote: CL = (CLbase + [CL_allele1 + CL_allele2] x CL_NFIB)
    # x (1 - CLage x [Age - AgeBP]) for Age >= AgeBP
    cl_cyp2d6 <- exp(lcl_cyp2d6_star1) * cyp2d6_act *
      e_snp_nfib_rs28379954_cl^SNP_NFIB_RS28379954
    age_hinge <- exp(lage_hinge)
    age_hinge_9oh <- exp(lage_hinge_9oh)
    age_cl <- 1 + e_age_cl * max(AGE - age_hinge, 0)
    age_cl_9oh <- 1 + e_age_cl_9oh * max(AGE - age_hinge_9oh, 0)

    ka <- exp(lka)
    cl <- (exp(lcl_base) + cl_cyp2d6) * age_cl * exp(etalcl)
    vc <- exp(lvc)
    cl_9oh <- exp(lcl_9oh) * age_cl_9oh * exp(etalcl_9oh)
    vc_9oh <- exp(lvc_9oh)

    kel <- cl / vc
    kel_9oh <- cl_9oh / vc_9oh

    # All parent clearance forms metabolite; the unidentifiable fmet (and
    # the 410.5 -> 426.5 g/mol mass change) is absorbed into the apparent
    # metabolite volume and clearance.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(central_9oh) <- kel * central - kel_9oh * central_9oh

    # mg/L -> ng/mL
    Cc <- 1000 * central / vc
    Cc_9oh <- 1000 * central_9oh / vc_9oh

    Cc ~ prop(propSd) + add(addSd)
    Cc_9oh ~ prop(propSd_9oh)
  })
}
