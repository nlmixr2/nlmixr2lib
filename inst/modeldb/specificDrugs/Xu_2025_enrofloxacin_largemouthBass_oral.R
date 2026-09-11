Xu_2025_enrofloxacin_largemouthBass_oral <- function() {
  description <- paste(
    "Veterinary (largemouth bass, Micropterus salmoides).",
    "One-compartment population PK model with first-order absorption and no",
    "absorption lag for enrofloxacin after a single 20 mg/kg body-weight oral",
    "(gavage) dose in largemouth bass held at 25.0 +/- 0.5 degC. Xu 2025",
    "fitted 24 sparsely sampled fish (four plasma samples each, drawn from a",
    "16-point 0.083-120 h grid) in Phoenix NLME 8.0 by first-order conditional",
    "estimation, extended least squares. Dose, clearance and volume are all",
    "body-weight normalised in the source (dose in mg/kg, tvV in L/kg, tvCL in",
    "L/h/kg), so the model is driven with amt in mg/kg and returns Cc directly",
    "in ug/mL. Final model, Table 4 'Without importing covariate of wt' block:",
    "tvKa 0.98 1/h, tvV 6.82 L/kg, tvCL 0.098 L/h/kg, with a Phoenix",
    "multiplicative residual error, stdev0 = 0.30, encoded as prop(propSd).",
    "V and CL are APPARENT (per bioavailability); the companion intravenous",
    "model Xu_2025_enrofloxacin_largemouthBass_iv gives the absolute values,",
    "and the two together imply an oral bioavailability of 12.24%. Body",
    "weight was screened on Ka, V and CL and rejected (Table 3), so no",
    "covariate enters the model. Exponential between-fish variability was",
    "fitted on all three structural parameters but its magnitude is never",
    "published, so the three etas are carried at fixed(0); see the vignette.",
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
  # BW and V/F, CL/F per kg BW. verified = TRUE: the analyte is enrofloxacin
  # throughout (Section 2.3.1, oral gavage of EF solution; the metabolite
  # ciprofloxacin was not measured in this study) and the sampled matrix is
  # plasma (Section 2.3.1: caudal-vessel blood centrifuged at 1500 x g).
  compartmentData <- list(
    depot = list(analyte = "enrofloxacin", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "enrofloxacin", units = "mg/kg", specimen = "plasma", verified = TRUE)
  )

  # Xu 2025 screened body weight on every structural parameter and rejected
  # it. Section 2.6 lists body weight, gender, blood indicators and urine
  # indicators as candidate covariates but states that only body weight is
  # readily obtainable in fish, so body weight is the only covariate tested.
  # No covariate is referenced in model(), so the screen is recorded here
  # rather than in covariateData.
  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight of the individual fish. Xu 2025 Table 3 screened body weight against Ka, V and CL in all eight on/off combinations for the oral dataset. The lowest -2LL with the fewest parameters was the V-wt scenario (-2LL 167.01, AIC 183.01, 8 parameters, versus 169.67 / 183.67 / 7 for no covariate), so V-wt was carried forward and fitted (Table 4, 'Importing covariate of wt': dVdwt = -0.39, RSE -43.90%, 95% CI -0.73 to -0.049). Section 3.3 then rejected it because the standard errors and RSEs of tvKa, tvV and tvCL were essentially unchanged relative to the covariate-free fit, so the final reported model carries no covariate.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened, fitted, and then rejected; not retained in the final model. As in the companion IV model, weight scaling is nevertheless built into the units rather than fitted: the dose is mg/kg BW and V/F and CL/F are per kg BW, i.e. exponent-1 proportionality is assumed a priori. The rejected covariate test was therefore for a DEPARTURE from that proportionality, not for the presence of any weight effect at all. The fish used were 243.11 +/- 54.86 g (Section 2.2); no per-fish weights are published, so the dVdwt centring/reference weight cannot be recovered and the covariate form is not reproducible even if a user wanted it."
    )
  )

  population <- list(
    species            = "largemouth bass (Micropterus salmoides)",
    n_subjects         = 24L,
    n_studies          = 1L,
    weight_mean        = "243.11 +/- 54.86 g (the 60-fish purchase lot from which both route cohorts were drawn; Section 2.2)",
    dose_range         = "single 20 mg/kg body weight oral gavage (20 mg/mL enrofloxacin in pure water, delivered by plastic tube into the stomach from a 2.5 mL microinjector)",
    disease_state      = "healthy",
    regions            = "China (Yangtze River Fisheries Research Institute, Wuhan; fish sourced from Huazhong Agricultural University)",
    water_temperature  = "25.0 +/- 0.5 degC (air-conditioned)",
    water_quality      = "480 L tanks at 26 L/min; dissolved oxygen, total ammonia nitrogen, nitrite nitrogen and pH checked daily and held in the ranges of Xu 2023 (Animals 13:1749)",
    design             = "24 fish randomised to four groups of six. Fish that regurgitated the gavage solution were removed and replaced. At least 14 days of acclimation on antibiotic-free feed before dosing (Sections 2.2 and 2.3.1).",
    sampling           = "Sparse: each fish sampled 4 times from the caudal vessels. Four repeating schedules of four times each (Section 2.3.1): group 1 (5 min, 1, 8, 48 h), group 2 (10 min, 2, 12, 72 h), group 3 (15 min, 4, 16, 96 h), group 4 (0.5, 6, 24, 120 h), so all 16 nominal times are covered by six fish each.",
    bioanalysis        = "HPLC with fluorescence detection (excitation 280 nm, emission 450 nm), Poroshell 120 EC-C18; LOD 0.003 ug/mL, LOQ 0.01 ug/mL; recovery 83.29-103.12%, intra-day RSD 2.07-3.21%, inter-day RSD 3.01-6.72% (Sections 2.5 and 3.1, Table 1).",
    notes              = "The same 24-fish design was run at 10 mg/kg intravenously and fitted separately; see Xu_2025_enrofloxacin_largemouthBass_iv. Only one water temperature was studied. Xu 2025 Section 4 contrasts these estimates with a conventional (non-population) PK study of the same species at 28 degC (Shan 2019, J Vet Pharmacol Ther 43:147-152) that reported Ka 10.20 1/h, V 2.21 L/kg, CL 0.017 L/h/kg and AUC 1185.73 ug.h/mL, i.e. very different values; the parameters here should not be extrapolated to other rearing temperatures. No individual concentration or weight data are published (Data Availability Statement: available on request)."
  )

  ini({
    # =================================================================
    # Structural parameters -- Xu 2025 Table 4, FINAL model
    # ("Without importing covariate of wt" block), NOT the "Importing
    # covariate of wt" block and NOT the bootstrap block. Section 3.3
    # names these three values as the final result.
    #
    # All three are apparent (per bioavailability) values and are body-
    # weight normalised: V/F is L/kg and CL/F is L/h/kg, matched to a
    # dose expressed in mg/kg BW. No absorption lag was fitted.
    #
    # Note that the "Coefficient of Variation (%)" column of Table 4 is
    # the RELATIVE STANDARD ERROR of the estimate, not an inter-animal
    # CV: 0.50/6.82 = 7.3% reproduces the printed 7.34 for tvV, and
    # 0.046/0.30 = 15.3% reproduces the printed 15.10 for stdev0. It is
    # therefore NOT the CV(%) = 100 x sqrt(exp(omega^2) - 1) quantity
    # defined in Section 2.6, and must not be read as an IIV.
    # =================================================================
    lka <- log(0.98)
    label("Absorption rate constant Ka (1/h)")  # Xu 2025 Table 4, "Without importing covariate of wt": tvKa = 0.98 1/h (SE 0.14, RSE 14.00%, 95% CI 0.71-1.26). Also stated in Sections 3.3, 4 and 5. Bootstrap (n = 1000) gave 1.01.
    lvc <- log(6.82)
    label("Apparent volume of distribution V/F (L/kg)")  # Xu 2025 Table 4, "Without importing covariate of wt": tvV = 6.82 L/kg (SE 0.50, RSE 7.34%, 95% CI 5.82-7.81). Also stated in Sections 3.3, 4 and 5. Bootstrap (n = 1000) gave 6.90.
    lcl <- log(0.098)
    label("Apparent clearance CL/F (L/h/kg)")  # Xu 2025 Table 4, "Without importing covariate of wt": tvCL = 0.098 L/h/kg (SE 0.0099, RSE 10.16%, 95% CI 0.078-0.12). Also stated in Sections 3.3, 4 and 5. Bootstrap (n = 1000) gave 0.095.

    # =================================================================
    # Between-fish (inter-individual) variability
    #
    # Xu 2025 Section 2.6 states that an exponential model,
    # P_i = tvP x exp(eta_Pi) with eta ~ N(0, omega^2), was used for
    # inter-individual and inter-occasion variability, and defines
    # CV(%) = 100 x sqrt(exp(omega^2) - 1). The omega^2 values are
    # NEVER PUBLISHED -- neither Table 4 nor any other table or figure
    # reports them, and the only variability column in Table 4 is the
    # RSE of the fixed effect (see the note above).
    #
    # That three etas WERE estimated is nevertheless recoverable from
    # the parameter counts in Table 3: the covariate-free oral model has
    # 7 parameters, of which tvKa, tvV, tvCL and stdev0 account for 4,
    # leaving exactly 3 -- one omega per structural parameter. Each
    # "-wt" scenario adds exactly one parameter, confirming the count.
    #
    # Per the standing policy for unreported IIV with structural values
    # present, the three etas are carried at fixed(0) rather than
    # invented. Simulation from this model is therefore deterministic in
    # the structural parameters; a user who wants a stochastic cohort
    # must supply their own omega. See the vignette's "Assumptions and
    # deviations" section.
    # =================================================================
    etalka ~ fixed(0)  # Xu 2025 Section 2.6 (exponential IIV declared) + Table 3 parameter count = 7; magnitude never reported
    etalvc ~ fixed(0)  # Xu 2025 Section 2.6 (exponential IIV declared) + Table 3 parameter count = 7; magnitude never reported
    etalcl ~ fixed(0)  # Xu 2025 Section 2.6 (exponential IIV declared) + Table 3 parameter count = 7; magnitude never reported

    # =================================================================
    # Residual unexplained variability
    #
    # Xu 2025 Section 2.6: "The multiplicative model was employed to
    # calculate the residual unexplained variability", with
    #     OBSV = IPCN x (1 + eps),  eps ~ N(0, stdev0^2)
    # i.e. a proportional residual on the linear concentration scale,
    # which is exactly nlmixr2's `~ prop(propSd)`. stdev0 is a fraction
    # and is dimensionless; Table 4's "-" entry in its Units column is
    # consistent with that.
    # =================================================================
    propSd <- 0.30
    label("Proportional residual error (fraction)")  # Xu 2025 Table 4, "Without importing covariate of wt": stdev0 = 0.30 (SE 0.046, RSE 15.10%, 95% CI 0.21-0.39). Bootstrap (n = 1000) also gave 0.30.
  })

  model({
    # 1. Individual parameters. Exponential (log-normal) between-fish
    #    variability, matching Phoenix NLME's multiplicative eta
    #    structure of Section 2.6. All three etas are fixed at zero
    #    because their variances are unpublished (see ini()).
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)

    # 2. Micro-constant. Section 2.6's secondary-parameter equation
    #    Ke = tvCL/tvV gives 0.098/6.82 = 0.0144 1/h, matching the
    #    published Ke of 0.014 1/h and T1/2Ke = ln2/Ke = 49.50 h.
    kel <- cl / vc

    # 3. ODE system. Xu 2025 Section 4: the oral concentration-time data
    #    were best fitted by a one-compartment model with first-order
    #    absorption; the two- and three-compartment models could not be
    #    fitted at all (their parameter CVs were not computable). No lag
    #    time and no estimated bioavailability -- V and CL are apparent,
    #    i.e. already divided by F.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 4. Observation. Doses enter as mg/kg body weight and vc is L/kg, so
    #    central/vc is mg/L == ug/mL, the paper's concentration unit.
    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
