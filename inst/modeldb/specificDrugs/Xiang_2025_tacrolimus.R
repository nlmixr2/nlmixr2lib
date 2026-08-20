Xiang_2025_tacrolimus <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral tacrolimus in adult renal transplant recipients, built from twice-daily trough (C0) therapeutic drug monitoring data. Because no absorption-phase samples were collected the absorption rate constant is fixed at a literature value and the residual error is additive. Apparent clearance rises as a power function of post-operative day and is reduced by concomitant Wuzhi capsule and by the CYP3A5*3/*3 non-expresser genotype; no covariate was retained on apparent volume of distribution. This is the population PK step of a three-model paper; the two sequential exposure-response models built on it are Xiang_2025_tacrolimus_fpg (fasting plasma glucose) and Xiang_2025_tacrolimus_egfr (estimated glomerular filtration rate)."
  reference <- paste(
    "Xiang Q, Yang Y, Li G, Chen S, Yang Y, Liu L, Yu X.",
    "Population Pharmacokinetic/Pharmacodynamic Modeling of Tacrolimus in Renal",
    "Transplant Recipients: Impact of CYP3A5 Genotype and Wuzhi Capsule",
    "Co-Medication. Drug Des Devel Ther. 2025;19:8375-8389.",
    "doi:10.2147/DDDT.S542786.",
    "Parameter values are the final-model estimates in Table 2 ('PPK model'",
    "block); the structural equations are Eq. 1, 5 and 6.",
    sep = " "
  )
  vignette <- "Xiang_2025_tacrolimus"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Tacrolimus is assayed in WHOLE BLOOD, not plasma
  # (Xiang 2025 Methods, "Biological Analysis and Data Collection": whole-blood
  # trough concentrations measured by enzyme-amplified immunoassay), so V/F and
  # CL/F are whole-blood apparent parameters.
  compartmentData <- list(
    depot   = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    POD = list(
      description        = "Post-operative day: days elapsed since renal transplantation",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = "TIME-VARYING within subject. Enters as the power scaling (POD / 34)^0.109 on CL/F; 34 days is the median post-operative day of the index cohort (Xiang 2025 Table 1). Apparent clearance therefore RISES with time since transplantation, which is why simulated C0 falls over the first post-operative month at a constant dose (Xiang 2025 Results, 'Simulations'). Unlike the linear centred form used by some transplant popPK models, the power form is undefined at POD = 0: CL/F collapses to zero on the day of surgery. Supply POD >= 1 (tacrolimus was started within 24 h of transplantation, so POD = 1 is the first dosing day). The index cohort spans POD 0-182 days and no upper cap was applied.",
      source_name        = "POD"
    ),
    CONMED_WUZHI = list(
      description        = "Concomitant Wuzhi capsule (Schisandra sphenanthera extract) indicator (1 = receiving Wuzhi capsule, 0 = not)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant Wuzhi capsule)",
      notes              = "Source column WZ, defined in Xiang 2025 immediately below Eq. 6 as 'WZ = 1 if WZ is present; otherwise = 0'. Enters as exp(-0.211 * CONMED_WUZHI) on CL/F, a 19% reduction in apparent clearance. 33 of 100 index-cohort subjects (33%) took Wuzhi capsules. This co-medication is prescribed deliberately to raise tacrolimus exposure and is therefore confounded with genotype in the source cohort: Xiang 2025 Discussion notes it is 'often prescribed for patients with CYP3A5*1/*1 but rarely for patients with CYP3A5*3/*3', with only 11 subjects in the CYP3A5*3/*3 + Wuzhi cell (Supplementary Table 6).",
      source_name        = "WZ"
    ),
    CYP3A5_EXPR = list(
      description        = "CYP3A5 expresser status (1 = carries at least one CYP3A5*1 allele, genotype *1/*1 or *1/*3; 0 = CYP3A5*3/*3 non-expresser)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5*3/*3 non-expresser)",
      notes              = "VALUE INVERSION relative to the source. Xiang 2025 codes the opposite orientation immediately below Eq. 6: 'Genotype = 1 if the genotype is CYP3A5*3/*3; otherwise, it = 0'. The canonical register mandates the expresser-equals-1 orientation, so CYP3A5_EXPR = 1 - Genotype and the model applies the published coefficient to (1 - CYP3A5_EXPR). The direction of the effect is unchanged from the paper: non-expressers have exp(-0.381) = 0.68 times the apparent clearance of expressers, i.e. 32% lower CL/F, matching Xiang 2025 Results ('CYP3A5*3/*3 genotype reduced the tacrolimus CL/F'). Xiang 2025 pooled *1/*1 (8%) and *1/*3 (31%) into a single CYP3A5*1 expresser group because the three-level genotype effect could not be estimated precisely (RSE > 30%; Discussion); 61% of the index cohort are *3/*3 non-expressers.",
      source_name        = "Genotype"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 100L,
    n_studies      = 1L,
    age_range      = "19-65 years (median 38)",
    age_median     = "38 years",
    weight_range   = "31-99 kg (median 60)",
    weight_median  = "60 kg",
    sex_female_pct = 36.0,
    race_ethnicity = c(Han = 92.0, Other = 8.0),
    disease_state  = "Adult renal transplant recipients on a tacrolimus-based immunosuppressive regimen, from the day of transplantation through the first post-operative months (post-operative day 0-182, median 34)",
    dose_range     = "Oral tacrolimus started within 24 h of transplantation at 0.05-0.15 mg/kg/day divided q12h, then dose-adjusted to the trough TDM target (8-12 ng/mL in month 1, 6-10 ng/mL months 1-3, 4-10 ng/mL months 3-12); observed doses 0.5-5.5 mg per administration (median 2.5)",
    regions        = "Single centre, Chongqing, People's Republic of China",
    notes          = "Retrospective cohort of 126 renal transplant recipients (September 2021 to March 2024) randomly split 4:1 into an index group (n = 100) used to build the model and an external validation group (n = 26); the parameter values here are the index-group fit. 2279 tacrolimus trough concentrations. Baseline demographics from Xiang 2025 Table 1 (index dataset column). Genotype distribution CYP3A5*1/*1 8%, *1/*3 31%, *3/*3 61%; 33% took Wuzhi capsules, 73% amlodipine, 65% metoprolol, 99% omeprazole. Estimated in Phoenix NLME 8.3.5 by FOCE. Registered as NCT05872815."
  )

  ini({
    # -------------------------------------------------------------------------
    # Structural parameters. Xiang 2025 Table 2, 'PPK model' block, Final Model
    # column; the equations are Eq. 5 (V/F) and Eq. 6 (CL/F).
    #
    # Only trough concentrations were collected, so no absorption-phase
    # information exists in the data and Ka was fixed at a literature value.
    # Methods states that 3.09 /h and 4.5 /h were both trialled (from the
    # authors' references 16 and 24) and the Discussion states "we therefore
    # fixed Ka at 3.09 h-1 as the final parameter". Table 2 lists "3.09 FIX" in
    # both the base and final model columns. The Results text's "the Ka was
    # fixed to 3.9 h-1" is a typographical error: 3.9 is neither of the two
    # candidate values and contradicts both Table 2 and the Discussion.
    # -------------------------------------------------------------------------
    lka <- fixed(log(3.09)); label("Absorption rate constant Ka (1/h)")                 # Xiang 2025 Table 2 PPK model: Ka = 3.09 h^-1 FIX; literature value, no absorption-phase data
    lvc <- log(2248);        label("Apparent volume of distribution V/F (L)")           # Xiang 2025 Eq. 5 and Table 2: V/F = 2248 L (RSE 10.4%; bootstrap median 2217, 95% CI 1660-2789)
    lcl <- log(33.7);        label("Apparent clearance CL/F at POD 34 d, CYP3A5 expresser, no Wuzhi capsule (L/h)")  # Xiang 2025 Eq. 6 and Table 2: CL/F = 33.7 L/h (RSE 6.9%; bootstrap median 33.4, 95% CI 27.3-39.4)

    # -------------------------------------------------------------------------
    # Covariate effects on CL/F. Xiang 2025 Eq. 6:
    #   CL/F = 33.7 * (POD/34)^0.109 * e^(-0.211*WZ) * e^(-0.381*Genotype) * e^(etaCL)
    # with WZ = 1 when the Wuzhi capsule is present and Genotype = 1 for
    # CYP3A5*3/*3. The genotype coefficient is applied here to
    # (1 - CYP3A5_EXPR) because the canonical covariate carries the opposite
    # (expresser-equals-1) orientation; see covariateData$CYP3A5_EXPR.
    # -------------------------------------------------------------------------
    e_pod_cl          <- 0.109;  label("Power exponent of (POD / 34 days) on CL/F (unitless)")                                   # Xiang 2025 Table 2: POD on CL/F = 0.109 (RSE 4.7%; bootstrap 95% CI 0.052-0.164)
    e_conmed_wuzhi_cl <- -0.211; label("Exponential effect of concomitant Wuzhi capsule on CL/F (unitless)")                     # Xiang 2025 Table 2: WZ on CL/F = -0.211 (RSE 8.5%; bootstrap 95% CI -0.433 to -0.047)
    e_cyp3a5_expr_cl  <- -0.381; label("Exponential effect of CYP3A5 non-expresser status on CL/F, applied to (1 - CYP3A5_EXPR) (unitless)")  # Xiang 2025 Table 2: CYP3A5*3/*3 on CL/F = -0.381 (RSE 21.1%; bootstrap 95% CI -0.554 to -0.188)

    # -------------------------------------------------------------------------
    # Inter-individual variability. Xiang 2025 Eq. 1 declares an exponential
    # random-effect model, P_i = tvP * exp(eta_i), and Table 2 reports the IIV
    # of each parameter as a percentage. The percentages are read here as
    # coefficients of variation of the log-normal parameter distribution and
    # converted to the internal variance scale with omega^2 = log(CV^2 + 1);
    # see the vignette's Assumptions and deviations section.
    # -------------------------------------------------------------------------
    etalcl ~ 0.1009994  # Xiang 2025 Table 2: IIV CL/F = 32.6% (RSE 17.1%; bootstrap 95% CI 27.0-37.3) -> log(1 + 0.326^2)
    etalvc ~ 0.5794162  # Xiang 2025 Table 2: IIV V/F  = 88.6% (RSE 19.3%; bootstrap 95% CI 65.5-122.5) -> log(1 + 0.886^2)

    # -------------------------------------------------------------------------
    # Residual error. Methods tested additive, proportional and mixed models;
    # Results reports "The additive error model best described the residuals of
    # the PK model."
    # -------------------------------------------------------------------------
    addSd <- 2.20; label("Additive residual error for tacrolimus whole-blood concentration (ng/mL)")  # Xiang 2025 Table 2: Residual error Additive = 2.20 ng/mL (RSE 1.0%; bootstrap 95% CI 2.02-2.33)
  })

  model({
    # -----------------------------------------------------------------------
    # Individual PK parameters. Written in the multiplicative form Xiang 2025
    # Eq. 6 prints, rather than folding the covariate terms inside the
    # exponential, so each factor maps one-to-one onto a published term.
    # -----------------------------------------------------------------------
    ka <- exp(lka)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl) *
      (POD / 34)^e_pod_cl *
      exp(e_conmed_wuzhi_cl * CONMED_WUZHI) *
      exp(e_cyp3a5_expr_cl * (1 - CYP3A5_EXPR))

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # -----------------------------------------------------------------------
    # Observation. The dose is in mg and V/F is in L, so central / vc is in
    # mg/L; tacrolimus troughs are reported in ng/mL and 1 mg/L = 1000 ng/mL.
    # -----------------------------------------------------------------------
    Cc <- 1000 * central / vc

    Cc ~ add(addSd)
  })
}
