Xia_2025_sertraline <- function() {
  description <- paste(
    "One-compartment first-order absorption population PK model for sertraline in Chinese",
    "adolescent (13-17 years) inpatients with depressive disorders (Xia 2025), built from",
    "221 routine therapeutic-drug-monitoring trough concentrations in 103 patients.",
    "The model carries NO covariates: height, body weight, age, gender and concomitant",
    "quetiapine / olanzapine / alprazolam were all screened by forward inclusion and none",
    "reduced the objective function value by the required 6.63, which the authors attribute",
    "to the narrow 13-17 year age band and the small weight spread of the cohort.",
    "The absorption rate constant is held at ka = 0.5 1/h taken from Poweleit 2023 because",
    "every sample was an elimination-phase trough and the absorption phase was not",
    "identifiable. Typical CL/F = 65.8 L/h and V/F = 1570 L give kel = 0.0419 1/h and a",
    "16.5 h half-life. The paper's purpose is dosing-remediation simulation: it derives",
    "recommended remedial doses for one, two and three consecutive missed 50 / 100 / 200 mg",
    "QD doses as a function of how late the dose is taken.",
    sep = " "
  )
  reference <- paste(
    "Xia H, Deng G, Liang F, Zhang Z, Huang W, Guo Z, Song Q, Wen Y, Shang D, Tan Y.",
    "Investigating Remedial Strategies for Missed or Delayed Dose of Sertraline in Chinese",
    "Adolescent Patients with Depressive Disorders via Population Pharmacokinetics Modeling",
    "and Simulation Approaches. Drug Des Devel Ther. 2025;19:3001-3016.",
    "doi:10.2147/DDDT.S504521. PMID 40260198. PMCID PMC12011033.",
    "The fixed absorption rate constant ka = 0.5 1/h is taken from",
    "Poweleit EA, Taylor ZL, Mizuno T, et al. Escitalopram and sertraline population",
    "pharmacokinetic analysis in pediatric patients. Clin Pharmacokinet. 2023;62(11):1621-1637.",
    "doi:10.1007/s40262-023-01294-8 (Xia 2025 reference 25).",
    "Xia 2025 reference 34 is the same group's earlier adult-and-adolescent sertraline model,",
    "packaged here as modellib('Zhang_2024_sertraline').",
    sep = " "
  )
  vignette <- "Xia_2025_sertraline"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # No covariate was retained in the final model, so covariateData is empty and
  # every screened covariate is documented in covariatesDataExcluded below.
  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age at the therapeutic-drug-monitoring sample",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened on the median-centred linear form of Xia 2025 Eq. (3) and not retained",
        "(Results, Model Establishment: 'no significant covariates were screened').",
        "Cohort median (range) 15 (13-17) years. The Discussion is explicit that this is a",
        "power limitation rather than evidence of no age effect: 'we did not detect a",
        "significant effect of age on sertraline clearance as the age of adolescent patients",
        "was only distributed within a small range of 13-17'. The same group's earlier",
        "adult-plus-adolescent cohort DID retain age on CL/F - see",
        "modellib('Zhang_2024_sertraline'), which encodes",
        "CL/F = 76.1 * [1 - 0.0068 * (AGE - 22)] over an 11-79 year range.",
        sep = " "
      ),
      source_name        = "Age (Xia 2025 Table 1)"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened via the median-centred linear form of Eq. (3) and not retained.",
        "Cohort median (range) 57.5 (39.9-89) kg in males and 50 (34-91) kg in females",
        "(Table 1); 5.80% of weights were missing and were replaced by the within-sex median.",
        "The Discussion contrasts this with Poweleit 2023 (reference 25), where sertraline",
        "CL and V rose by 4.9% and 16.1% per 10 kg, and attributes the null result to the",
        "relatively small weight difference between the adolescents studied here.",
        "No allometric scaling was applied in the final model.",
        sep = " "
      ),
      source_name        = "Weight (Xia 2025 Table 1)"
    ),
    HT = list(
      description        = "Body height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened via the median-centred linear form of Eq. (3) and not retained.",
        "Cohort median (range) 169 (150-176) cm in males and 159.5 (150-172) cm in females",
        "(Table 1). 37.90% of heights were missing and were replaced by the within-sex median,",
        "which further limits the power of this screen.",
        sep = " "
      ),
      source_name        = "Height (Xia 2025 Table 1)"
    ),
    SEXF = list(
      description        = "Sex, 1 = female",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Screened via the categorical-covariate form of Eq. (4) and not retained.",
        "Xia 2025 states directly under Eq. (4) that 'For gender covariates, a COV of 1",
        "represents females and 0 represents males', which matches the SEXF orientation with",
        "no value transformation. Cohort was 30 male (29.00%) / 73 female (71.00%) (Table 1).",
        sep = " "
      ),
      source_name        = "Gender (Xia 2025 Table 1)"
    ),
    CONMED_QUETIAPINE = list(
      description        = "Concomitant quetiapine coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant quetiapine)",
      notes              = paste(
        "Screened via Eq. (4) (COV of 1 represents co-administration, 0 represents no",
        "co-administration) and not retained. The most common comedication in the cohort:",
        "43 patients (41.75%) per Table 1.",
        sep = " "
      ),
      source_name        = "Quetiapine (Xia 2025 Table 1)"
    ),
    CONMED_OLANZAPINE = list(
      description        = "Concomitant olanzapine coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant olanzapine)",
      notes              = paste(
        "Screened via Eq. (4) and not retained. 32 patients (31.07%) per Table 1.",
        sep = " "
      ),
      source_name        = "Olanzapine (Xia 2025 Table 1)"
    ),
    CONMED_ALPRAZOLAM = list(
      description        = "Concomitant alprazolam coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant alprazolam)",
      notes              = paste(
        "Screened via Eq. (4) and not retained. 34 patients (33.00%) per Table 1.",
        "The Discussion notes that no comedication reached significance because of the low",
        "frequency of concomitant use, and that the same three psychiatric drug classes were",
        "likewise non-significant in the group's earlier cohort (reference 34).",
        sep = " "
      ),
      source_name        = "Alprazolam (Xia 2025 Table 1)"
    ),
    CYP2C19_PHENO = list(
      description        = "CYP2C19 metaboliser phenotype group",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = "Normal / extensive metaboliser",
      notes              = paste(
        "NOT screened - no genotype data existed. Xia 2025 Discussion: 'the CYP2C19 genotype",
        "was not tested in patients in this study, and therefore the effect of CYP2C19",
        "genotyping on our PPK model could not be assessed.' Recorded here because the paper",
        "flags it as the single most important unmeasured source of the residual between-subject",
        "variability: it quotes the CPIC guideline that 13% of East Asian individuals are poor",
        "metabolisers, 46% intermediate and 38% normal. N-desmethyl sertraline concentrations",
        "were likewise unavailable and were not screened. No point estimate is published,",
        "so nothing is encoded.",
        sep = " "
      ),
      source_name        = "CYP2C19 (Xia 2025 Discussion)"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "sertraline", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sertraline", units = "mg", specimen = "serum",               verified = TRUE)
  )

  population <- list(
    species          = "human",
    n_subjects       = 103L,
    n_studies        = 1L,
    n_concentrations = 221L,
    age_range        = "13-17 years",
    age_median       = "15 years",
    weight_range     = "39.9-89 kg (male); 34-91 kg (female)",
    weight_median    = "57.5 kg (male); 50 kg (female)",
    height_range     = "150-176 cm (male); 150-172 cm (female)",
    height_median    = "169 cm (male); 159.5 cm (female)",
    sex_female_pct   = 71,
    race_ethnicity   = "Chinese (single-centre Guangzhou cohort; non-Chinese patients were an explicit exclusion criterion, and the cohort is not further stratified by the source)",
    disease_state    = "Hospitalised adolescents with depressive disorders receiving oral sertraline",
    dose_range       = "Daily dose median 150 mg, range 50-300 mg; the simulations use 50, 100 and 200 mg QD",
    regions          = "China (The Affiliated Brain Hospital of Guangzhou Medical University, Guangzhou, Guangdong)",
    co_medication    = "Quetiapine 43 patients (41.75%); alprazolam 34 (33.00%); olanzapine 32 (31.07%)",
    sampling         = "Retrospective therapeutic drug monitoring; all 221 samples are elimination-phase troughs drawn at 6-7 a.m. before the next scheduled dose. Median (range) observed concentration 63.54 (5.07-299.87) ug/L.",
    notes            = paste(
      "Baseline demographics from Xia 2025 Table 1. Retrospective TDM data collected",
      "1 January 2019 to 31 December 2023; IRB approval 2021027. Serum sertraline quantified",
      "by HPLC-MS/MS over a 5-500 ng/mL calibrated range with intra- and inter-day precision",
      "below 15% RSE; the AGNP therapeutic reference range used is 10-150 ng/mL and the",
      "laboratory alert level is 300 ng/mL. Inclusion required at least two sertraline",
      "concentrations per patient within one hospitalisation. Missing heights (37.90%) and",
      "weights (5.80%) were replaced by within-sex medians. Model fitted in NONMEM 7.3.0",
      "with FOCE-I; evaluated by goodness-of-fit plots, NPDE (variance 1.14, mean 0.152) and",
      "a 1000-run bootstrap that converged in 928 runs (92.8%).",
      sep = " "
    )
  )

  ini({
    # ----------------------------------------------------------------------
    # Structural parameters: Xia 2025 Table 2, 'Final Model - Estimate' column.
    # Cross-checked against the paper's own derived quantities: the Discussion
    # quotes CL/V = 0.042 1/h for this model, and 65.8 / 1570 = 0.04191 1/h
    # reproduces it. See the vignette for the two independent structural gates
    # (steady-state tmax and steady-state trough) that these values pass.
    # ----------------------------------------------------------------------
    lka <- fixed(log(0.5))
    label("First-order absorption rate constant ka (1/h); literature value taken from Poweleit 2023")  # Xia 2025 Table 2 'Ka (h-1) = 0.5 fixed'; Methods, Model Development fixes ka to the value of reference 25 because no sample was drawn in the absorption phase
    lcl <- log(65.8)
    label("Apparent oral clearance CL/F (L/h)")  # Xia 2025 Table 2: CL/F = 65.8, RSE 6%; bootstrap median 65.96 (95% CI 59.11-72.82)
    lvc <- log(1570)
    label("Apparent volume of distribution V/F (L)")  # Xia 2025 Table 2: V/F = 1570, RSE 19%; bootstrap median 1652.18 (95% CI 1129.10-2401.14)

    # ----------------------------------------------------------------------
    # Inter-individual variability, exponential on CL/F and V/F per Eq. (1)
    # P_ij = P_tv,j * EXP(eta_i), with eta ~ N(0, omega^2).
    #
    # SCALE OF THE TABLE 2 'IIV (%)' COLUMN. The column is headed "(%)" but
    # holds 0.134 and 0.326, which are not percentage-magnitude numbers. They
    # are read here as the raw NONMEM $OMEGA variances (omega^2), giving
    # CV = sqrt(exp(omega^2) - 1) of 37.9% on CL/F and 62.1% on V/F. Three
    # independent checks support the variance reading over reading them as
    # standard deviations (which would give 13.4% and 32.6%):
    #
    #  1. The residual rows of the SAME table are provably raw variances, and a
    #     single paste of a NONMEM run would not take the square root of $OMEGA
    #     while leaving $SIGMA untransformed. The bootstrap 95% CI on the
    #     proportional term runs down to 0.01; as a variance that is a 10% CV,
    #     but as a standard deviation it would be a 1% CV - below the assay
    #     precision the paper itself reports (intra- and inter-day RSE < 15%),
    #     so the SD reading is impossible. Symmetrically, the additive term's
    #     bootstrap upper bound of 344.63 read as an SD would exceed the largest
    #     concentration in the whole dataset (299.87 ug/L).
    #  2. RSE consistency. The standard error of a typical value is on the order
    #     of omega/sqrt(N). With N = 103 and the variance reading, omega on CL/F
    #     is sqrt(0.134) = 0.366 and omega/sqrt(N) = 3.6%, the same order as the
    #     6% RSE Table 2 reports. Under the SD reading omega would be 0.134 and
    #     the bound would be 1.3%, i.e. the estimator would have to be 4.6 times
    #     more precise than it actually was, for the parameter that trough data
    #     determine most directly.
    #  3. House style. The same group's earlier paper (reference 34, packaged as
    #     modellib('Zhang_2024_sertraline')) headed its column "IIV (CV%)" and
    #     printed percentage-magnitude integers there (10 and 57), while printing
    #     its proportional residual as a raw variance (0.129) in the estimate
    #     column. Xia 2025 prints fraction-magnitude numbers in both places,
    #     i.e. the untransformed NONMEM output throughout.
    #
    # This is an interpretation, not a printed value; see the vignette's
    # Assumptions and deviations section.
    # ----------------------------------------------------------------------
    etalcl ~ 0.134  # Xia 2025 Table 2 'IIV (%)' on CL/F, read as the NONMEM $OMEGA variance (CV 37.9%); see the note above
    etalvc ~ 0.326  # Xia 2025 Table 2 'IIV (%)' on V/F, read as the NONMEM $OMEGA variance (CV 62.1%); see the note above

    # ----------------------------------------------------------------------
    # Residual error. Eq. (2) is the combined form Y = F * (1 + eps1) + eps2,
    # with eps1 ~ N(0, sigma1^2) proportional and eps2 ~ N(0, sigma2^2)
    # additive. Both rows of Table 2 sit in the 'Estimate' column and are read
    # as the raw NONMEM $SIGMA variances, so the nlmixr2 standard deviations are
    # their square roots: sqrt(0.0959) = 0.3097 (31.0% proportional) and
    # sqrt(24.7) = 4.970 ug/L additive. The "(%)" unit printed on the ADD row is
    # a template artefact - an additive residual on a concentration cannot be a
    # percentage - which is itself a sign that the row labels were not curated.
    # See the argument in the IIV note above for why both are variances.
    # ----------------------------------------------------------------------
    propSd <- sqrt(0.0959)
    label("Proportional residual error (fraction)")  # Xia 2025 Table 2 'PRO (%)' = 0.0959 read as the NONMEM $SIGMA variance; bootstrap median 0.09 (95% CI 0.01-0.12)
    addSd <- sqrt(24.7)
    label("Additive residual error (ug/L)")  # Xia 2025 Table 2 'ADD (%)' = 24.7 read as the NONMEM $SIGMA variance; bootstrap median 26.28 (95% CI 2.68-344.63)
  })

  model({
    # Individual parameters. No covariate entered the final model, so these are
    # the typical values perturbed only by the exponential IIV of Eq. (1).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # One-compartment first-order absorption with first-order elimination.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose is in mg and vc in L, so central/vc is mg/L; the assay, the 10-150
    # ng/mL therapeutic reference range and the 300 ng/mL alert level are all in
    # ng/mL (== ug/L), so scale by 1000.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
