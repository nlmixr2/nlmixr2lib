Li_2017_pomalidomide_renalcat <- function() {
  description <- "One-compartment population PK model with first-order absorption and an absorption lag time for orally administered pomalidomide in patients with relapsed or refractory multiple myeloma (rrMM) and various degrees of renal impairment (Li 2017, Table 3). Renal function enters as a four-level CATEGORICAL covariate on apparent clearance, encoded here as three mutually exclusive binary indicators against a normal-renal-function reference: moderate impairment (RENALIMP_MOD, CL/F ratio 1.02), severe impairment not requiring hemodialysis (RENALIMP_SEV, ratio 1.01), and severe impairment requiring hemodialysis (RRT_HEMODIAL_STATUS, ratio 0.724). Interindividual variability on Ka, and a correlated 2x2 block on V/F and CL/F; residual error is additive on the log scale (log-normal). This is one of three renal-covariate parameterisations the paper reports; see Li_2017_pomalidomide_crcl and Li_2017_pomalidomide_egfr for the continuous reverse-hockey-stick alternatives, and Li_2017_pomalidomide_hemodialysis for this model augmented with the observed dialyzer clearance."
  reference <- "Li Y, Wang X, O'Mara E, Dimopoulos MA, Sonneveld P, Weisel KC, Matous J, Siegel DS, Shah JJ, Kueenburg E, Sternas L, Cavanaugh C, Zaki M, Palmisano M, Zhou S. Population pharmacokinetics of pomalidomide in patients with relapsed or refractory multiple myeloma with various degrees of impaired renal function. Clin Pharmacol Adv Appl. 2017;9:133-145. doi:10.2147/CPAA.S144606"
  vignette <- "Li_2017_pomalidomide"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    RENALIMP_MOD = list(
      description = "Moderately impaired renal function (Li 2017 group 2)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal renal function, Li 2017 group 1: CrCl >= 60 mL/min)",
      notes = "Paper-specific classification scheme, NOT the FDA/EMA 30-59 mL/min/1.73 m^2 band: Li 2017 Table 1 defines the moderate group as 30 < eGFR < 45 mL/min/1.73 m^2 (n = 15, all from study CC-4047-MM-013). Mutually exclusive with RENALIMP_SEV and RRT_HEMODIAL_STATUS; all three are 0 for the normal-renal-function reference group.",
      source_name = "renal function group 2"
    ),
    RENALIMP_SEV = list(
      description = "Severely impaired renal function NOT requiring hemodialysis (Li 2017 group 3)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal renal function, Li 2017 group 1: CrCl >= 60 mL/min)",
      notes = "Li 2017 Table 1: CrCl < 30 mL/min or eGFR < 30 mL/min/1.73 m^2, not requiring dialysis (n = 30). Groups 3 and 4 share this renal classification and are distinguished ONLY by whether the patient requires hemodialysis, which is why the group-4 indicator is the dialysis treatment-status flag RRT_HEMODIAL_STATUS rather than a further RENALIMP_* severity band. Mutually exclusive with RENALIMP_MOD and RRT_HEMODIAL_STATUS.",
      source_name = "renal function group 3"
    ),
    RRT_HEMODIAL_STATUS = list(
      description = "Severely impaired renal function REQUIRING hemodialysis (Li 2017 group 4)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not receiving intermittent hemodialysis)",
      notes = "Li 2017 Table 1: CrCl < 30 mL/min or eGFR < 30 mL/min/1.73 m^2 requiring dialysis (n = 10). Subject-level, time-fixed treatment-status flag; it does NOT switch on and off per dialysis session. The estimated CL/F ratio of 0.724 is the clearance on NON-dialysis days (Li 2017 abstract: 'pomalidomide exposure increased approximately 35% in patients with severe renal impairment requiring dialysis on nondialysis days'). For the within-session dialyzer clearance use the per-time-point RRT_HEMODIAL_ACTIVE gate in Li_2017_pomalidomide_hemodialysis. Mutually exclusive with RENALIMP_MOD and RENALIMP_SEV.",
      source_name = "renal function group 4"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "pomalidomide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pomalidomide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 63,
    n_studies = 2,
    age_range = "46-86 years",
    age_median = "69 years",
    weight_range = "39.9-116.6 kg",
    weight_median = "78.2 kg",
    disease_state = "relapsed or refractory multiple myeloma (rrMM) with normal, moderately impaired, or severely impaired renal function (with and without hemodialysis)",
    dose_range = "2-4 mg once daily, oral",
    renal_function = "CrCl median 28.3 mL/min (range 8.7-115.4); eGFR median 27 mL/min/1.73 m^2 (range 5-84). Group sizes (Table 1): normal 8, moderate 15, severe without hemodialysis 30, severe requiring hemodialysis 10.",
    co_medication = "low-dose dexamethasone",
    regions = "North America (USA, Canada) and Europe (France, Germany, Greece, Italy, Netherlands, Spain, UK, Austria)",
    notes = "Pooled intensive and sparse plasma sampling from studies CC-4047-MM-008 (Phase 1, n = 20) and CC-4047-MM-013 (Phase 2, n = 43). Baseline demographics are Li 2017 Table 2; renal-function group sizes are Table 1. Sex distribution and race/ethnicity are not reported in the source."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters. Li 2017 Table 3 ('final model with renal
    # function treated as a categorical variable'), 'Estimate' column.
    # ------------------------------------------------------------------
    lka <- log(0.724); label("Absorption rate constant (1/h)") # Table 3: Ka = 0.724 1/h (90% bootstrap CI 0.533-1.23)
    lvc <- log(58.3); label("Apparent central volume of distribution V/F (L)") # Table 3: V/F = 58.300 L (90% bootstrap CI 51.597-66.055)
    lcl <- log(5.13); label("Apparent clearance CL/F in patients with normal renal function (L/h)") # Table 3: CL/F normal = 5.130 L/h (90% bootstrap CI 4.283-6.069)
    ltlag <- log(0.154); label("Absorption lag time (h)") # Table 3: Tlag = 0.154 h (90% bootstrap CI 0.002-0.43)

    # ------------------------------------------------------------------
    # Renal-function categorical effect on CL/F. Li 2017 Equation 5 writes
    # the categorical covariate model as P = theta * (1 + theta_cov * Z),
    # so theta_cov = ratio - 1; Table 3 tabulates the RATIOS directly, and
    # they are transcribed verbatim here (the model() block subtracts 1).
    # Group 1 (normal renal function) is the reference: all three
    # indicators are 0 and CL/F reduces to exp(lcl).
    # ------------------------------------------------------------------
    e_renalimp_mod_cl <- 1.02; label("CL/F ratio, moderate renal impairment vs normal (unitless)") # Table 3: CL/F ratio (group 2 vs group 1) = 1.020 (90% bootstrap CI 0.848-1.195)
    e_renalimp_sev_cl <- 1.01; label("CL/F ratio, severe renal impairment without hemodialysis vs normal (unitless)") # Table 3: CL/F ratio (group 3 vs group 1) = 1.010 (90% bootstrap CI 0.871-1.169)
    e_rrt_hemodial_status_cl <- 0.724; label("CL/F ratio, severe renal impairment requiring hemodialysis vs normal (unitless)") # Table 3: CL/F ratio (group 4 vs group 1) = 0.724 (90% bootstrap CI 0.606-0.844)

    # ------------------------------------------------------------------
    # Interindividual variability. Li 2017 Equation 1 (P_i = P * exp(eta),
    # eta ~ N(0, omega^2)), so every omega below is a LOG-scale variance
    # and is transcribed from Table 3 without back-transformation.
    #
    # Table 3 reports a 2x2 OMEGA BLOCK on V/F and CL/F. The row labelled
    # 'omega(V/F):omega(CL/F)' = 0.163 is the COVARIANCE, not a
    # correlation; three independent checks agree:
    #   (1) the Results prose says 'The estimated IIV and associated
    #       covariance were reasonably precise';
    #   (2) the table's Notes define only omega^2 as a variance, leaving
    #       the single non-omega^2 row as the block off-diagonal, which is
    #       the NONMEM $OMEGA BLOCK convention;
    #   (3) the bootstrap 90% CI for this row, 0.088-0.236, brackets the
    #       geometric mean of the two variance CIs,
    #       sqrt(0.069*0.121) = 0.091 and sqrt(0.236*0.273) = 0.254 --
    #       i.e. the off-diagonal tracks sqrt(var_V * var_CL) across the
    #       whole bootstrap distribution, which only holds for a
    #       covariance with a near-unity correlation.
    # Implied correlation 0.163 / sqrt(0.139 * 0.198) = 0.982; the block is
    # positive definite (determinant 0.139*0.198 - 0.163^2 = 9.5e-4 > 0,
    # eigenvalues 0.334 and 0.0029). A near-degenerate V/F-CL/F block is
    # the expected signature of an ORAL model in which both parameters
    # carry the same unestimated 1/F, so the shared bioavailability
    # variability dominates both.
    # ------------------------------------------------------------------
    etalka ~ 0.782 # Table 3: omega^2 (Ka) = 0.782 (90% bootstrap CI 0.441-1.387)
    etalvc + etalcl ~ c(
      0.139,
      0.163, 0.198
    ) # Table 3: omega^2 (V/F) = 0.139; omega(V/F):omega(CL/F) = 0.163; omega^2 (CL/F) = 0.198

    # ------------------------------------------------------------------
    # Residual error. Li 2017 Equation 2: ln(C_obs) = ln(C_pred) + eps,
    # i.e. ADDITIVE on the natural-log scale, which is exactly nlmixr2's
    # lnorm() error model (not prop()). Table 3 reports the variance
    # delta^2 = 0.341, so the log-scale SD is sqrt(0.341) = 0.5840.
    # ------------------------------------------------------------------
    expSd <- 0.584; label("Log-scale residual standard deviation (log ng/mL)") # Table 3: delta^2 = 0.341 (90% bootstrap CI 0.194-0.459); sqrt(0.341) = 0.5840
  })

  model({
    # 1. Individual parameters.
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)

    # Renal-function categorical multiplier on CL/F, in the Li 2017
    # Equation 5 form P = theta * (1 + theta_cov * Z) with
    # theta_cov = (tabulated ratio) - 1. The three indicators are mutually
    # exclusive, so at most one term is non-zero and the multiplier equals
    # that group's tabulated CL/F ratio; with all three at 0 the multiplier
    # is 1 and CL/F is the normal-renal-function value.
    cl <- exp(lcl + etalcl) * (1 +
      (e_renalimp_mod_cl - 1) * RENALIMP_MOD +
      (e_renalimp_sev_cl - 1) * RENALIMP_SEV +
      (e_rrt_hemodial_status_cl - 1) * RRT_HEMODIAL_STATUS)

    tlag <- exp(ltlag)

    # 2. Micro-constants.
    kel <- cl / vc

    # 3. One-compartment model with first-order oral absorption.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 4. Absorption lag time on the depot.
    alag(depot) <- tlag

    # 5. Observation. Dose is in mg and vc in L, so central / vc is mg/L;
    #    multiply by 1000 to give the ng/mL in which Li 2017 reports every
    #    concentration (LLOQ 0.25 ng/mL; observed range 0.255-180.569
    #    ng/mL) and every AUC (Table 4, ng/mL*h).
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
