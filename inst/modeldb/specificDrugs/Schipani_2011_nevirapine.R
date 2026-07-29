Schipani_2011_nevirapine <- function() {
  description <- "One-compartment population PK model for oral nevirapine in HIV-infected adults (Schipani 2011), with CYP2B6 516G>T (rs3745274) and 983T>C (rs28399499) genotype and body-weight covariate effects on CL/F. Covariate effects are additive on linear-scale CL/F per the published equation."
  reference <- paste(
    "Schipani A, Wyen C, Mahungu T, Hendra H, Egan D, Siccardi M, Davies G, Khoo S,",
    "Fatkenheuer G, Rockstroh J, Brockmeyer NH, Johnson MA, Owen A, Back DJ.",
    "Integration of population pharmacokinetics and pharmacogenetics: an aid to",
    "optimal nevirapine dose selection in HIV-infected individuals.",
    "J Antimicrob Chemother. 2011;66(6):1332-1339. doi:10.1093/jac/dkr087."
  )
  vignette <- "Schipani_2011_nevirapine"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (baseline weight). Linear additive shift on CL/F with reference body weight 72.5 kg (Schipani 2011 Table 1 cohort median). Schipani 2011 final-model equation: CL/F increases by 0.018 L/h per 1 kg above 72.5 kg.",
      source_name        = "BW"
    ),
    SNP_CYP2B6_RS3745274_T_COUNT = list(
      description        = "Count of CYP2B6 c.516G>T (rs3745274, p.Q172H) T-alleles per subject (0/1/2). 0 = GG homozygous wild-type, 1 = GT heterozygous, 2 = TT homozygous variant.",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (germline genotype). The source paper encoded two mutually-exclusive binary indicators X_516GT and X_516TT in NONMEM with independently estimated additive shifts on CL/F; the canonical count column reconstructs them as (count == 1) and (count == 2). Schipani 2011 Table 1 distribution: GG 47%, GT 46%, TT 7%.",
      source_name        = "X_516GT and X_516TT (paired indicators)"
    ),
    SNP_CYP2B6_RS28399499_C_COUNT = list(
      description        = "Count of CYP2B6 c.983T>C (rs28399499, p.I328T) C-alleles per subject (0/1/2). 0 = TT homozygous wild-type, 1 = TC heterozygous, 2 = CC homozygous variant.",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed (germline genotype). The source paper estimated a single additive shift for X_983TC heterozygotes; the canonical count column reconstructs the indicator as (count == 1). No 983CC homozygotes were observed in the Schipani 2011 cohort or described elsewhere in the published literature as of 2011, so the homozygous effect is not estimated. Schipani 2011 Table 1 distribution: TT 97%, TC 3%, CC 0%.",
      source_name        = "X_983TC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 272,
    n_studies      = 2,
    age_range      = "18-82 years",
    age_median     = "42 years",
    weight_range   = "47-132 kg",
    weight_median  = "72.5 kg",
    sex_female_pct = 41.5,
    race_ethnicity = c(Caucasian = 66.5, Black = 33.5),
    disease_state  = "HIV-positive adults on nevirapine-based antiretroviral therapy with virological suppression (HIV viral load < 50 copies/mL at the time of sampling), in combination with 2 NRTIs or 1 NRTI + 1 NtRTI; patients on ritonavir-boosted PIs or other major interacting drugs were excluded.",
    dose_range     = "200 mg orally twice daily (n = 237) or 400 mg orally once daily (n = 38), at steady state.",
    regions        = "UK (Royal Free NHS Trust, London) + Germany (KompNet Cohort)",
    cyp2b6_freq    = "516G>T: GG 47%, GT 46% (n = 126), TT 7% (n = 19). 983T>C: TT 97%, TC 3% (n = 9), CC 0% (Schipani 2011 Table 1).",
    notes          = "Pooled across two cohorts (n = 275 enrolled, 272 in the final model after 3 non-adherent subjects with peak concentrations < 1 mg/L were excluded). One cohort had intensive PK sampling (11 patients with 6 samples per dosing interval, 9 patients with two occasions of one random sample); the remaining patients contributed one random TDM sample at steady state. Total of 403 nevirapine concentrations were analyzed."
  )

  ini({
    # ---- Structural PK parameters (Schipani 2011 Table 2 final model) ----
    lka  <- log(1.20)  ; label("First-order absorption rate constant ka (1/h)")             # Schipani 2011 Table 2 final-model ka = 1.20 h^-1 (RSE 22%)
    lcl  <- log(3.51)  ; label("Apparent oral clearance CL/F at reference (L/h)")           # Schipani 2011 Table 2 final-model CL/F = 3.51 L/h (RSE 3.0%); reference = BW 72.5 kg, CYP2B6 wild-type (516GG / 983TT)
    lvc  <- log(150)   ; label("Apparent volume of distribution V/F (L)")                   # Schipani 2011 Table 2 final-model V/F = 150 L (RSE 8.7%)

    # ---- Covariate effects on CL/F (linear additive per Schipani 2011 equation) ----
    # Final-model equation: TVCL = theta0 + theta_BW * (BW - 72.5)
    #                            + theta_516GT * X_516GT + theta_516TT * X_516TT
    #                            + theta_983TC * X_983TC
    e_wt_cl      <-  0.018 ; label("Linear-additive shift on CL/F per 1 kg above 72.5 kg (L/h/kg)")          # Schipani 2011 Table 2: theta_BW = 0.018 (RSE 32.8%); +0.18 L/h per +10 kg, +5%/10 kg
    e_516gt_cl   <- -0.5   ; label("Linear-additive shift on CL/F for CYP2B6 516GT heterozygotes (L/h)")     # Schipani 2011 Table 2: theta_516GT = -0.5 (RSE 27.3%); 14% lower CL/F vs 516GG reference
    e_516tt_cl   <- -1.3   ; label("Linear-additive shift on CL/F for CYP2B6 516TT homozygotes (L/h)")       # Schipani 2011 Table 2: theta_516TT = -1.3 (RSE 16.7%); 37% lower CL/F vs 516GG reference
    e_983tc_cl   <- -1.4   ; label("Linear-additive shift on CL/F for CYP2B6 983TC heterozygotes (L/h)")     # Schipani 2011 Table 2: theta_983TC = -1.4 (RSE 13.2%); 40% lower CL/F vs 983TT reference

    # ---- IIV (only CL/F has IIV in the final model; absorption and V/F do not) ----
    # Schipani 2011 Table 2: IIV CL/F = 31% CV in the final model (RSE 10.8%)
    # omega^2 = log(CV^2 + 1) = log(0.31^2 + 1) = 0.09171
    etalcl ~ 0.09171

    # ---- Residual error (proportional only; additive did not improve the fit) ----
    propSd <- 0.092 ; label("Proportional residual error (fraction)")                       # Schipani 2011 Table 2: proportional residual error = 9.2% CV (RSE 10.8%); sqrt(0.0085) reported in Results
  })

  model({
    # 1. Derive non-additive genotype indicators from the per-allele count columns.
    #    The Schipani 2011 parameterization treats heterozygous and homozygous
    #    variant carriers as separate categorical groups (not a linear per-allele
    #    effect), so we decompose the count into mutually-exclusive indicators.
    het_516 <- (SNP_CYP2B6_RS3745274_T_COUNT == 1)
    hom_516 <- (SNP_CYP2B6_RS3745274_T_COUNT == 2)
    het_983 <- (SNP_CYP2B6_RS28399499_C_COUNT == 1)

    # 2. Typical CL/F: linear-additive covariate model on linear-scale CL,
    #    matching the Schipani 2011 final equation. exp(lcl) recovers the
    #    reference-cohort CL/F (3.51 L/h) at BW = 72.5 kg, 516GG, 983TT.
    tvcl <- exp(lcl) +
            e_wt_cl    * (WT - 72.5) +
            e_516gt_cl * het_516 +
            e_516tt_cl * hom_516 +
            e_983tc_cl * het_983

    # 3. Individual PK parameters
    cl <- tvcl * exp(etalcl)
    vc <- exp(lvc)
    ka <- exp(lka)

    # 4. Micro-constants
    kel <- cl / vc

    # 5. ODE system: one-compartment model with first-order absorption
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 6. Observation and error
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
