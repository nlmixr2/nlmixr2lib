Li_2017_pomalidomide_hemodialysis <- function() {
  description <- "One-compartment population PK model with first-order absorption and an absorption lag time for orally administered pomalidomide in patients with relapsed or refractory multiple myeloma (rrMM), augmented with an extracorporeal dialyzer clearance arm so that hemodialysis-day profiles can be simulated (Li 2017, Table 3 plus the observed CL_D). The disposition parameters and the four-level categorical renal-function effect on CL/F are those of Li_2017_pomalidomide_renalcat; on top of them an ADDITIVE dialyzer clearance of 12 L/h is gated on by the time-varying RRT_HEMODIAL_ACTIVE covariate for the duration of each session. Reproduces the paper's Figure 5 and Figure 6 scenarios, which show that starting hemodialysis after a dose cuts the dosing-interval exposure to roughly half, whereas completing hemodialysis before the dose preserves it -- the basis for the recommendation to give pomalidomide post-dialysis."
  reference <- "Li Y, Wang X, O'Mara E, Dimopoulos MA, Sonneveld P, Weisel KC, Matous J, Siegel DS, Shah JJ, Kueenburg E, Sternas L, Cavanaugh C, Zaki M, Palmisano M, Zhou S. Population pharmacokinetics of pomalidomide in patients with relapsed or refractory multiple myeloma with various degrees of impaired renal function. Clin Pharmacol Adv Appl. 2017;9:133-145. doi:10.2147/CPAA.S144606"
  vignette <- "Li_2017_pomalidomide"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    RENALIMP_MOD = list(
      description = "Moderately impaired renal function (Li 2017 group 2)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal renal function, Li 2017 group 1: CrCl >= 60 mL/min)",
      notes = "Paper-specific classification scheme: Li 2017 Table 1 defines the moderate group as 30 < eGFR < 45 mL/min/1.73 m^2 (n = 15). Mutually exclusive with RENALIMP_SEV and RRT_HEMODIAL_STATUS.",
      source_name = "renal function group 2"
    ),
    RENALIMP_SEV = list(
      description = "Severely impaired renal function NOT requiring hemodialysis (Li 2017 group 3)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal renal function, Li 2017 group 1: CrCl >= 60 mL/min)",
      notes = "Li 2017 Table 1: CrCl < 30 mL/min or eGFR < 30 mL/min/1.73 m^2, not requiring dialysis (n = 30). Mutually exclusive with RENALIMP_MOD and RRT_HEMODIAL_STATUS.",
      source_name = "renal function group 3"
    ),
    RRT_HEMODIAL_STATUS = list(
      description = "Severely impaired renal function REQUIRING hemodialysis (Li 2017 group 4)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (not receiving intermittent hemodialysis)",
      notes = "Subject-level, time-fixed flag (Li 2017 Table 1, n = 10). Carries the CL/F ratio of 0.724 that applies on NON-dialysis days. It is deliberately distinct from RRT_HEMODIAL_ACTIVE: this one says the patient is a dialysis patient, that one says a session is running right now. A group-4 patient therefore has RRT_HEMODIAL_STATUS = 1 at all times and RRT_HEMODIAL_ACTIVE = 1 only inside each session.",
      source_name = "renal function group 4"
    ),
    RRT_HEMODIAL_ACTIVE = list(
      description = "Hemodialysis session currently running",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no dialysis session running)",
      notes = "Time-varying WITHIN subject: 1 for the duration of each hemodialysis session, 0 in the interdialytic interval and at all times in non-dialysed patients. Gates the additive cl_hemodialysis arm. Li 2017 simulated 4 h sessions (Figures 5 and 6) and reports that a 6 h session gave similar results ('data not shown').",
      source_name = "hemodialysis procedure duration (Figures 5, 6)"
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
    disease_state = "relapsed or refractory multiple myeloma (rrMM); the dialyzer-clearance arm was derived from the 10 patients with severe renal impairment requiring hemodialysis",
    dose_range = "2-4 mg once daily, oral",
    renal_function = "CrCl median 28.3 mL/min (range 8.7-115.4); eGFR median 27 mL/min/1.73 m^2 (range 5-84). CL_D was computed in 7 of the 10 hemodialysis patients.",
    co_medication = "low-dose dexamethasone",
    regions = "North America (USA, Canada) and Europe (France, Germany, Greece, Italy, Netherlands, Spain, UK, Austria)",
    notes = "Disposition parameters from studies CC-4047-MM-008 and CC-4047-MM-013 pooled (Li 2017 Table 2, Table 3). The dialyzer clearance was measured separately from paired arterial-side and venous-side plasma samples drawn across the dialyzer on hemodialysis days."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters and renal-function categorical effect: exactly
    # Li 2017 Table 3, identical to Li_2017_pomalidomide_renalcat. See
    # that file for the full source-trace commentary on the OMEGA block.
    # ------------------------------------------------------------------
    lka <- log(0.724); label("Absorption rate constant (1/h)") # Table 3: Ka = 0.724 1/h (90% bootstrap CI 0.533-1.23)
    lvc <- log(58.3); label("Apparent central volume of distribution V/F (L)") # Table 3: V/F = 58.300 L (90% bootstrap CI 51.597-66.055)
    lcl <- log(5.13); label("Apparent clearance CL/F in patients with normal renal function (L/h)") # Table 3: CL/F normal = 5.130 L/h (90% bootstrap CI 4.283-6.069)
    ltlag <- log(0.154); label("Absorption lag time (h)") # Table 3: Tlag = 0.154 h (90% bootstrap CI 0.002-0.43)

    e_renalimp_mod_cl <- 1.02; label("CL/F ratio, moderate renal impairment vs normal (unitless)") # Table 3: CL/F ratio (group 2 vs group 1) = 1.020 (90% bootstrap CI 0.848-1.195)
    e_renalimp_sev_cl <- 1.01; label("CL/F ratio, severe renal impairment without hemodialysis vs normal (unitless)") # Table 3: CL/F ratio (group 3 vs group 1) = 1.010 (90% bootstrap CI 0.871-1.169)
    e_rrt_hemodial_status_cl <- 0.724; label("CL/F ratio, severe renal impairment requiring hemodialysis vs normal (unitless)") # Table 3: CL/F ratio (group 4 vs group 1) = 0.724 (90% bootstrap CI 0.606-0.844)

    # ------------------------------------------------------------------
    # Dialyzer clearance. Li 2017 Equation 6 defines
    #   CL_D = Q * (C_a - C_v) / C_a
    # from the blood flow rate Q through the dialyzer and the paired
    # arterial-side (entering) and venous-side (exiting) plasma
    # concentrations. Results: 'CL_D was calculated and the median value
    # from seven patients was approximately 12 L/h, which is higher than
    # pomalidomide total body clearance (5 L/h)'.
    #
    # NOT a NONMEM estimate -- it is an observed median from 7 patients
    # that the paper substituted into its simulations ('Simulations based
    # on the final PPK model and observed CL_D'), so it is encoded as
    # fixed(). No uncertainty or IIV is reported for it.
    #
    # The arm is ADDITIVE to the intrinsic body clearance, not a
    # replacement for it. Equation 6 defines an extraction clearance
    # across an extracorporeal circuit operating in parallel with the
    # body's own elimination, and the Results and Discussion both keep
    # CL_D and total body CL as two separate quantities ('the median
    # values of the CL_D were 2-fold higher than pomalidomide total body
    # clearance (12 and 5 L/h, respectively)'). The additive reading is
    # also the one that reproduces the paper's own simulation results:
    # for a typical group-4 patient it gives a scenario-1 dosing-interval
    # exposure of 58-74% of a non-dialysis day against the paper's
    # reported 50-70%, and a scenario-2 exposure of 87% against the
    # paper's 83-91%; a replacement reading (total CL = 12 L/h while
    # dialysing) gives 68-80% and 90%, missing most of the scenario-1
    # range. The abstract's looser phrasing, 'dialysis increased total
    # body pomalidomide clearance from 5 L/h to 12 L/h', should be read
    # against Equation 6 and not as a replacement rule.
    # ------------------------------------------------------------------
    lcl_hemodialysis <- fixed(log(12)); label("Dialyzer clearance CL_D while a hemodialysis session is running (L/h)") # Results, 'Assessment of the extent to which hemodialysis contributes...': median CL_D from 7 patients ~12 L/h, via Equation 6

    # ------------------------------------------------------------------
    # Interindividual variability and residual error: Li 2017 Table 3.
    # The off-diagonal 0.163 is the V/F-CL/F COVARIANCE (implied
    # correlation 0.982); see Li_2017_pomalidomide_renalcat for the three
    # checks that settle the covariance-vs-correlation reading.
    # ------------------------------------------------------------------
    etalka ~ 0.782 # Table 3: omega^2 (Ka) = 0.782 (90% bootstrap CI 0.441-1.387)
    etalvc + etalcl ~ c(
      0.139,
      0.163, 0.198
    ) # Table 3: omega^2 (V/F) = 0.139; omega(V/F):omega(CL/F) = 0.163; omega^2 (CL/F) = 0.198

    expSd <- 0.584; label("Log-scale residual standard deviation (log ng/mL)") # Table 3: delta^2 = 0.341 (90% bootstrap CI 0.194-0.459); sqrt(0.341) = 0.5840
  })

  model({
    # 1. Individual parameters.
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)

    cl_hemodialysis <- exp(lcl_hemodialysis) # L/h

    # 2. TOTAL clearance: the intrinsic body CL/F carrying the renal
    #    categorical effect, PLUS the dialyzer arm while a session runs.
    #
    #    IMPORTANT -- the gated sum MUST be assigned to `cl` itself, never
    #    to a separate `cl_total` variable. For some model shapes rxode2
    #    solves the system with its analytic linear-compartment kernel
    #    driven by variables named `cl` and `vc` and discards the explicit
    #    d/dt() right-hand side; routing the gated total through another
    #    name then leaves the dialysis arm silently INERT, with `cl_total`
    #    switching correctly while the simulated amounts decay at the
    #    interdialytic rate in both states. Same hazard and same fix as
    #    Dohmann_2025_piperacillin.R and Lee_2024_gentamicin_teigen.R;
    #    regression-tested in tests/testthat/test-modeldb-active-gate.R.
    cl <- exp(lcl + etalcl) * (1 +
      (e_renalimp_mod_cl - 1) * RENALIMP_MOD +
      (e_renalimp_sev_cl - 1) * RENALIMP_SEV +
      (e_rrt_hemodial_status_cl - 1) * RRT_HEMODIAL_STATUS) +
      RRT_HEMODIAL_ACTIVE * cl_hemodialysis

    tlag <- exp(ltlag)

    # 3. Micro-constants.
    kel <- cl / vc

    # 4. One-compartment model with first-order oral absorption.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 5. Absorption lag time on the depot.
    alag(depot) <- tlag

    # 6. Observation, in the ng/mL Li 2017 reports (mg/L * 1000).
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
