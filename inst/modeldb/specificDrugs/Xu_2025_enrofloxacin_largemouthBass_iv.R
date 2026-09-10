Xu_2025_enrofloxacin_largemouthBass_iv <- function() {
  description <- paste(
    "Veterinary (largemouth bass, Micropterus salmoides).",
    "Two-compartment population PK model for enrofloxacin after a single",
    "10 mg/kg body-weight intravenous (caudal vein) dose in largemouth bass",
    "held at 25.0 +/- 0.5 degC. Xu 2025 fitted 24 sparsely sampled fish (four",
    "plasma samples each, drawn from a 16-point 0.083-120 h grid) in Phoenix",
    "NLME 8.0 by first-order conditional estimation, extended least squares.",
    "Dose, clearances and volumes are all body-weight normalised in the source",
    "(dose in mg/kg, tvV1 and tvV2 in L/kg, tvCL and tvCL2 in L/h/kg), so the",
    "model is driven with amt in mg/kg and returns Cc directly in ug/mL.",
    "Final model, Table 6 'Without importing covariate of wt' block:",
    "tvV1 0.57 L/kg, tvV2 1.28 L/kg, tvCL 0.012 L/h/kg, tvCL2 1.00 L/h/kg,",
    "with a Phoenix multiplicative residual error, stdev0 = 0.14, encoded as",
    "prop(propSd). These are absolute (not apparent) values; paired with the",
    "companion oral model Xu_2025_enrofloxacin_largemouthBass_oral they imply",
    "an oral bioavailability of 12.24%. Body weight was screened on V1, V2,",
    "CL and CL2 and rejected (Table 5), so no covariate enters the model.",
    "Exponential between-fish variability was fitted on all four structural",
    "parameters but its magnitude is never published, so the four etas are",
    "carried at fixed(0); see the vignette.",
    sep = " "
  )
  reference <- paste(
    "Xu N, Zhou S, Dong J, Li J, Ding Y, Ai X. (2025). Population",
    "Pharmacokinetics of Enrofloxacin in Micropterus salmoides Based on a",
    "Nonlinear Mixed Effect Model After Intravenous and Oral Administration.",
    "Animals 15(10):1362. doi:10.3390/ani15101362.",
    sep = " "
  )
  vignette <- "Xu_2025_enrofloxacin_largemouthBass"
  units <- list(
    time = "h",
    dosing = "mg/kg",
    concentration = "ug/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Both states hold enrofloxacin; the amounts are
  # body-weight normalised (mg/kg) because Xu 2025 reports the dose in mg/kg
  # BW and V1, V2, CL, CL2 per kg BW. verified = TRUE: the analyte is
  # enrofloxacin throughout (Section 2.3.2, caudal-vein injection of EF
  # solution; the metabolite ciprofloxacin was not measured in this study).
  # The sampled matrix is plasma (Section 2.3.2: caudal-vessel blood
  # centrifuged at 1500 x g). `peripheral1` is a mathematical distribution
  # compartment with no identified tissue, so its specimen is recorded as
  # plasma-equivalent drug in the same analyte; no tissue was assayed.
  compartmentData <- list(
    central = list(analyte = "enrofloxacin", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "enrofloxacin", units = "mg/kg", specimen = "not applicable", verified = TRUE)
  )

  # Xu 2025 screened body weight on every structural parameter and rejected
  # it. Section 2.6 lists body weight, gender, blood indicators and urine
  # indicators as candidate covariates but states that only body weight is
  # readily obtainable in fish, so body weight is the only covariate tested.
  # No covariate is referenced in model(), so the screen is recorded here
  # rather than in covariateData.
  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight of the individual fish. Xu 2025 Table 5 screened body weight against V1, V2, CL and CL2 in all sixteen on/off combinations for the intravenous dataset. Section 3.3 selected V2-wt as the lowest -2LL with the fewest parameters (-2LL 270.34, AIC 290.34, 10 parameters, versus 272.70 / 290.70 / 9 for no covariate) and fitted it (Table 6, 'Importing covariate of wt': dV2dwt = -0.21, SE 0.23, RSE -113.84%, 95% CI -0.67 to 0.26). Section 3.3 then rejected it on two grounds: the standard errors and RSEs of the four fixed effects barely moved relative to the covariate-free fit, and dV2dwt was estimated as a negative value.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened, fitted, and then rejected; not retained in the final model. Note that dV2dwt's own 95% CI (-0.67 to 0.26) spans zero and its RSE exceeds 100%, so the estimate is not distinguishable from no effect -- a stronger argument for dropping it than the sign argument the paper makes. As in the companion oral model, weight scaling is nevertheless built into the units rather than fitted: the dose is mg/kg BW and all volumes and clearances are per kg BW, i.e. exponent-1 proportionality is assumed a priori. The rejected covariate test was therefore for a DEPARTURE from that proportionality. The fish used were 243.11 +/- 54.86 g (Section 2.2); no per-fish weights are published, so the dV2dwt centring/reference weight cannot be recovered and the covariate form is not reproducible even if a user wanted it."
    )
  )

  population <- list(
    species            = "largemouth bass (Micropterus salmoides)",
    n_subjects         = 24L,
    n_studies          = 1L,
    weight_mean        = "243.11 +/- 54.86 g (the 60-fish purchase lot from which both route cohorts were drawn; Section 2.2)",
    dose_range         = "single 10 mg/kg body weight intravenous injection into the caudal vein (20 mg/mL enrofloxacin in pure water, delivered from a 1 mL microinjector)",
    disease_state      = "healthy",
    regions            = "China (Yangtze River Fisheries Research Institute, Wuhan; fish sourced from Huazhong Agricultural University)",
    water_temperature  = "25.0 +/- 0.5 degC (air-conditioned)",
    water_quality      = "480 L tanks at 26 L/min; dissolved oxygen, total ammonia nitrogen, nitrite nitrogen and pH checked daily and held in the ranges of Xu 2023 (Animals 13:1749)",
    design             = "24 fish randomised to four groups of six. Fish that bled heavily after needle withdrawal, or in which the needle translocated during injection, were removed and replaced. At least 14 days of acclimation on antibiotic-free feed before dosing (Sections 2.2 and 2.3.2).",
    sampling           = "Sparse: each fish sampled 4 times from the caudal vessels, on the same schedule as the oral cohort (Section 2.3.2 'The planning of sampling time points was consistent with oral administration'). Four repeating schedules of four times each: group 1 (5 min, 1, 8, 48 h), group 2 (10 min, 2, 12, 72 h), group 3 (15 min, 4, 16, 96 h), group 4 (0.5, 6, 24, 120 h), so all 16 nominal times are covered by six fish each.",
    bioanalysis        = "HPLC with fluorescence detection (excitation 280 nm, emission 450 nm), Poroshell 120 EC-C18; LOD 0.003 ug/mL, LOQ 0.01 ug/mL; recovery 83.29-103.12%, intra-day RSD 2.07-3.21%, inter-day RSD 3.01-6.72% (Sections 2.5 and 3.1, Table 1).",
    notes              = "The same 24-fish design was run at 20 mg/kg orally and fitted separately; see Xu_2025_enrofloxacin_largemouthBass_oral. Only one water temperature was studied. Section 4 compares these estimates with conventional PK studies in other species (allogynogenetic silver crucian carp, snakehead fish, brown trout) that all report substantially shorter terminal half-lives and higher clearances. No individual concentration or weight data are published (Data Availability Statement: available on request)."
  )

  ini({
    # =================================================================
    # Structural parameters -- Xu 2025 Table 6, FINAL model
    # ("Without importing covariate of wt" block), NOT the "Importing
    # covariate of wt" block and NOT the bootstrap block. Section 3.3
    # and Section 5 name these four values as the final result.
    #
    # Table 6's caption says "after oral administration"; this is a
    # copy-paste error from Table 4. Its Note, its parameter names
    # (tvV1/tvV2/tvCL2), the Section 3.3 text that introduces it, and
    # the Table 5 scenario list that precedes it all identify it as the
    # INTRAVENOUS model. See the vignette's Errata.
    #
    # All four are absolute (not apparent) values and are body-weight
    # normalised: volumes in L/kg and clearances in L/h/kg, matched to a
    # dose expressed in mg/kg BW.
    #
    # Note that the "Coefficient of Variation (%)" column of Table 6 is
    # the RELATIVE STANDARD ERROR of the estimate, not an inter-animal
    # CV: 0.08/1.28 = 6.3% reproduces the printed 6.55 for tvV2 and
    # 0.12/1.00 = 12.0% the printed 12.47 for tvCL2. It is therefore NOT
    # the CV(%) = 100 x sqrt(exp(omega^2) - 1) quantity defined in
    # Section 2.6, and must not be read as an IIV.
    # =================================================================
    lvc <- log(0.57)
    label("Central volume of distribution V1 (L/kg)")  # Xu 2025 Table 6, "Without importing covariate of wt": tvV1 = 0.57 L/kg (SE 0.05, RSE 9.63%, 95% CI 0.46-0.67). Also stated in Sections 3.3 and 5. Bootstrap (n = 1000) gave 0.54.
    lvp <- log(1.28)
    label("Peripheral volume of distribution V2 (L/kg)")  # Xu 2025 Table 6, "Without importing covariate of wt": tvV2 = 1.28 L/kg (SE 0.08, RSE 6.55%, 95% CI 1.11-1.44). Also stated in Sections 3.3 and 5. Bootstrap (n = 1000) gave 1.27. Section 4 quotes V1 + V2 = 1.85 L/kg as "the value of V".
    lcl <- log(0.012)
    label("Clearance CL (L/h/kg)")  # Xu 2025 Table 6, "Without importing covariate of wt": tvCL = 0.012 L/h/kg (SE reported as 0.00 at 2 dp, RSE 6.12%, 95% CI reported as 0.01-0.01 at 2 dp). Also stated in Sections 3.3, 4 and 5. Bootstrap (n = 1000) gave 0.013.
    lq <- log(1.00)
    label("Intercompartmental clearance CL2 (L/h/kg)")  # Xu 2025 Table 6, "Without importing covariate of wt": tvCL2 = 1.00 L/h/kg (SE 0.12, RSE 12.47%, 95% CI 0.75-1.25). Also stated in Sections 3.3 and 5. Bootstrap (n = 1000) gave 1.06.

    # =================================================================
    # Between-fish (inter-individual) variability
    #
    # Xu 2025 Section 2.6 states that an exponential model,
    # P_i = tvP x exp(eta_Pi) with eta ~ N(0, omega^2), was used for
    # inter-individual and inter-occasion variability, and defines
    # CV(%) = 100 x sqrt(exp(omega^2) - 1). Section 2.6 adds that the
    # intravenous model's "other parameterizations were consistent with
    # the above description", i.e. the same exponential IIV. The omega^2
    # values are NEVER PUBLISHED -- neither Table 6 nor any other table
    # or figure reports them, and the only variability column in Table 6
    # is the RSE of the fixed effect (see the note above).
    #
    # That four etas WERE estimated is nevertheless recoverable from the
    # parameter counts in Table 5: the covariate-free IV model has 9
    # parameters, of which tvV1, tvV2, tvCL, tvCL2 and stdev0 account for
    # 5, leaving exactly 4 -- one omega per structural parameter. Each
    # "-wt" scenario adds exactly one parameter (the four-covariate
    # scenario has 13), confirming the count.
    #
    # Per the standing policy for unreported IIV with structural values
    # present, the four etas are carried at fixed(0) rather than
    # invented. Simulation from this model is therefore deterministic in
    # the structural parameters; a user who wants a stochastic cohort
    # must supply their own omega. See the vignette's "Assumptions and
    # deviations" section.
    # =================================================================
    etalvc ~ fixed(0)  # Xu 2025 Section 2.6 (exponential IIV declared) + Table 5 parameter count = 9; magnitude never reported
    etalvp ~ fixed(0)  # Xu 2025 Section 2.6 (exponential IIV declared) + Table 5 parameter count = 9; magnitude never reported
    etalcl ~ fixed(0)  # Xu 2025 Section 2.6 (exponential IIV declared) + Table 5 parameter count = 9; magnitude never reported
    etalq ~ fixed(0)  # Xu 2025 Section 2.6 (exponential IIV declared) + Table 5 parameter count = 9; magnitude never reported

    # =================================================================
    # Residual unexplained variability
    #
    # Xu 2025 Section 2.6: "The multiplicative model was employed to
    # calculate the residual unexplained variability", with
    #     OBSV = IPCN x (1 + eps),  eps ~ N(0, stdev0^2)
    # i.e. a proportional residual on the linear concentration scale,
    # which is exactly nlmixr2's `~ prop(propSd)`. stdev0 is a fraction
    # and is dimensionless; Table 6's "-" entry in its Units column is
    # consistent with that.
    # =================================================================
    propSd <- 0.14
    label("Proportional residual error (fraction)")  # Xu 2025 Table 6, "Without importing covariate of wt": stdev0 = 0.14 (SE 0.02, RSE 11.98%, 95% CI 0.11-0.17). Bootstrap (n = 1000) gave 0.12.
  })

  model({
    # 1. Individual parameters. Exponential (log-normal) between-fish
    #    variability, matching Phoenix NLME's multiplicative eta
    #    structure of Section 2.6. All four etas are fixed at zero
    #    because their variances are unpublished (see ini()).
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    # 2. Micro-constants. Section 2.6's secondary-parameter equations
    #    K12 = tvCL2/V1 = 1.00/0.57 = 1.754 1/h and
    #    K21 = tvCL2/V2 = 1.00/1.28 = 0.781 1/h are exactly these two.
    #    kel = tvCL/tvV1 = 0.012/0.57 = 0.0211 1/h is NOT a published
    #    quantity: Xu 2025 instead reports a terminal rate constant
    #    Kbeta = tvCL/(tvV1 + tvV2) = 0.012/1.85 = 0.0065 1/h, which is
    #    the Vss-based approximation to the true beta eigenvalue of this
    #    system. The exact eigenvalue of this ODE system is 0.006449 1/h,
    #    0.57% below the Vss form -- immaterial next to the 6.1% RSE on tvCL,
    #    and the resulting terminal half-life (107.5 h) is within 1% of the
    #    published 106.62 h. The approximation is benign only because
    #    distribution is fast here (k21/beta ~ 120); the vignette gates that
    #    condition explicitly rather than assuming it.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system. Xu 2025 Section 4: the intravenous concentration-time
    #    data were best fitted by a two-compartment model without
    #    absorption. The dose is an IV bolus into `central` (Section
    #    2.3.2: injection into the caudal vein from a 1 mL microinjector;
    #    no infusion duration is reported).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Observation. Doses enter as mg/kg body weight and vc is L/kg, so
    #    central/vc is mg/L == ug/mL, the paper's concentration unit.
    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
