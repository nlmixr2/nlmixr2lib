Hong_2016_pregabalin <- function() {
  description <- paste(
    "Two-compartment population pharmacokinetic model for oral pregabalin in",
    "healthy Korean male volunteers, with Savic-parameterised transit-compartment",
    "absorption (continuous, non-integer chain length NN = 3.61) feeding a depot",
    "that empties at first-order ka. Mean transit time and ka are estimated",
    "separately for the overnight-fasted and the fed state; apparent clearance",
    "scales with creatinine clearance as a power function.",
    sep = " "
  )
  reference <- paste(
    "Hong T, Han S, Lee J, Jeon S, Yim DS.",
    "Comparison of oral absorption models for pregabalin: usefulness of transit",
    "compartment model. Drug Des Devel Ther. 2016;10:3995-4003.",
    "doi:10.2147/DDDT.S123318.",
    sep = " "
  )
  vignette <- "Hong_2016_pregabalin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance estimated with the Cockcroft-Gault equation, raw",
        "mL/min and NOT BSA-normalized (Hong 2016 'Covariate selection').",
        "Enters CL/F as the mean-normalized power term (CRCL / 120)^0.511;",
        "the reference 120 mL/min is the pooled five-study mean of Table 1",
        "(120.0 +/- 21.6 mL/min), which the Table 5 parameter description",
        "names explicitly as 'exponent of mean-normalized CLCR'.",
        sep = " "
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Observed range in the source cohort is narrow and supranormal: the",
        "per-study means span 110.8-132.8 mL/min in 88 healthy young Korean men",
        "screened to be free of organ-system disease, so the model carries no",
        "information about renal impairment. The paper makes this limitation",
        "the subject of its own Discussion: Bockbrader 2011 and Shoji 2011 both",
        "found pregabalin CL/F to rise proportionally with CLCR only up to a",
        "breakpoint at 107 mL/min and to plateau above it, whereas this cohort",
        "shows no breakpoint (Figure 4). The authors attribute the difference to",
        "the small, homogeneous healthy-volunteer sample rather than to a real",
        "absence of the plateau, so extrapolating this power term to renally",
        "impaired subjects is not supported.",
        sep = " "
      ),
      source_name = "CLCR"
    ),
    FED = list(
      description = paste(
        "Fed-versus-fasted indicator carried on the DOSE record: 1 = dose taken",
        "after a meal, 0 = dose taken after an overnight fast. Selects between",
        "the two separately-estimated absorption parameter sets (MTT and ka).",
        sep = " "
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (overnight fasted)",
      notes = paste(
        "Dose-record level and therefore time-varying within a subject: studies",
        "2 and 3 dosed the same subject fasted in the morning and fed in the",
        "evening (Table 2), so a single individual contributes both levels.",
        "The source pools every fed record into one level. The authors tested",
        "finer prandial distinctions and could not resolve them: 'details such",
        "as dosing time differences after meals (30 min after versus 4 h after",
        "a meal) were not successfully modeled', and the Discussion adds that",
        "'differences in food types (high fat or regular) or time after meal",
        "(0.5 h versus 4 h) were not significantly different or discernible'.",
        "FED_HIGHFAT is therefore deliberately NOT used even though study 4",
        "gave a high-fat breakfast. Relative bioavailability was not allowed to",
        "change with food because the paper's linear-trapezoidal AUC comparison",
        "found no significant meal effect on exposure, so FED alters only the",
        "RATE of absorption, never the extent.",
        sep = " "
      ),
      source_name = "fasted/fed status"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "pregabalin",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "pregabalin",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "pregabalin",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 88,
    n_studies = 5,
    n_observations = 1615,
    age_range = "20-45 years (protocol eligibility criterion)",
    age_median = "27.2 years (pooled mean, SD 5.0)",
    weight_median = "68.3 kg (pooled mean, SD 7.8)",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state = "healthy volunteers",
    dose_range = "150 mg orally, single dose or every 12 h for up to 3 days",
    regions = "Republic of Korea",
    renal_function = "creatinine clearance 120.0 +/- 21.6 mL/min (Cockcroft-Gault); per-study means 110.8-132.8 mL/min",
    notes = paste(
      "Healthy Korean male volunteers enrolled at the clinical trial center of",
      "Seoul St Mary's Hospital across five studies, all within 20% of ideal",
      "body weight (Table 1). Demographics did not differ significantly between",
      "studies (Kruskal-Wallis, all p > 0.05). Prandial conditions differed by",
      "study (Table 2): overnight fasting, 30 min after a regular or high-fat",
      "meal, or 4 h after a regular meal. Plasma pregabalin was assayed by",
      "LC/MS/MS at four different contract research organizations with LLOQs of",
      "30-100 ng/mL (Table 3). The cohort is all-male, so sex was not testable",
      "as a covariate.",
      sep = " "
    )
  )

  ini({
    # --- Structural disposition: two-compartment, first-order elimination ---
    # All disposition parameters are APPARENT (divided by the unknown oral
    # bioavailability F); the paper reports CL/F, V2/F, Q/F and V3/F throughout.
    lcl <- log(6.25)
    label("Typical apparent clearance CLt/F at CRCL = 120 mL/min (L/h)")
    # Hong 2016 Table 5 row 'CLt/F' = 6.25 (RSE 0.80%; bootstrap median 6.25, 95% CI 6.08-6.42)
    lvc <- log(18.0)
    label("Apparent volume of the central compartment V2/F (L)")
    # Hong 2016 Table 5 row 'V2/F' = 18.0 (RSE 1.78%; bootstrap median 17.7, 95% CI 4.65-22.7)
    lq <- log(26.5)
    label("Apparent intercompartmental clearance Q/F (L/h)")
    # Hong 2016 Table 5 row 'Q/F' = 26.5 (RSE 3.63%; bootstrap median 26.5, 95% CI 19.9-43.5)
    lvp <- log(27.0)
    label("Apparent volume of the peripheral compartment V3/F (L)")
    # Hong 2016 Table 5 row 'V3/F' = 27.0 (RSE 1.93%; bootstrap median 27.2, 95% CI 23.3-37.7)

    # --- Covariate effect on apparent clearance ---
    # Table 5 prints the covariate model as its parameter description:
    #   CL/F = CLt/F x (CLCR/120)^theta_CLCR
    # CLCR on CL was the only covariate retained in the final model.
    e_crcl_cl <- 0.511
    label("Exponent of mean-normalized creatinine clearance on CL/F (unitless)")
    # Hong 2016 Table 5 row 'theta CLCR' = 0.511 (RSE 8.96%; bootstrap median 0.509, 95% CI 0.405-0.601)

    # --- Savic transit-compartment absorption ---
    # Reference 12 of the paper is Savic RM et al., J Pharmacokinet Pharmacodyn
    # 2007;34(5):711-726, whose standard parameterisation places (NN + 1)
    # first-order transfers at the shared rate ktr = (NN + 1) / MTT and lets the
    # resulting gamma-density input drive a depot that empties at rate ka.
    # NN was ESTIMATED on a continuous scale (RSE 5.32%), so it is legitimately
    # non-integer and the analytical Erlang/gamma input form is required rather
    # than an integer chain of ODE states.
    lnn <- log(3.61)
    label("Number of Savic transit compartments NN (unitless, non-integer)")
    # Hong 2016 Table 5 row 'nn' = 3.61 (RSE 5.32%; bootstrap median 3.63, 95% CI 3.19-4.26); also named in the Table 4 final-model row 'Transit compartment + two compartments (nn = 3.61)'

    # MTT and ka are estimated SEPARATELY for the two prandial states rather
    # than as a reference value plus a fed effect; the model reproduces that
    # structure literally (see model() below).
    lmtt_fasted <- log(0.494)
    label("Mean transit time MTT in the overnight-fasted state (h, FED = 0)")
    # Hong 2016 Table 5 row 'MTT fast' = 0.494 (RSE 9.94%; bootstrap median 0.499, 95% CI 0.411-0.600)
    lmtt_fed <- log(0.879)
    label("Mean transit time MTT in the fed state (h, FED = 1)")
    # Hong 2016 Table 5 row 'MTT fed' = 0.879 (RSE 6.88%; bootstrap median 0.884, 95% CI 0.751-1.101). Discussion cross-check: food 'prolonged mean transit time by ~0.39 h'; 0.879 - 0.494 = 0.385 h.
    lka_fasted <- log(5.69)
    label("First-order absorption rate constant ka in the overnight-fasted state (1/h, FED = 0)")
    # Hong 2016 Table 5 row 'Ka fast' = 5.69 (RSE 33.57%; bootstrap median 5.73, 95% CI 2.94-9.11)
    lka_fed <- log(0.713)
    label("First-order absorption rate constant ka in the fed state (1/h, FED = 1)")
    # Hong 2016 Table 5 row 'Ka fed' = 0.713 (RSE 1.47%; bootstrap median 0.709, 95% CI 0.617-0.814). Discussion cross-check: food decreased the rate constant 'by 87.5% compared with overnight fasting'; 0.713 / 5.69 = 0.1253, i.e. -87.5%.

    # --- Bioavailability ---
    # F is not identifiable from oral-only data and was never estimated; the
    # paper folds it into the apparent parameters above. Fixed to 1 so that the
    # transit() input delivers the whole nominal dose and CL/F x AUC = Dose.
    lfdepot <- fixed(log(1))
    label("Oral bioavailability F (unitless; absorbed into the apparent CL/F and V/F)")
    # Hong 2016 Table 5 reports only apparent (/F) disposition parameters; F itself is not estimated and relative bioavailability was deliberately not allowed to vary with food (Results: 'The relative bioavailability changes caused by meals were not considered in the model')

    # --- Between-subject variability ---
    # Table 5 heads the random-effect block 'Random effect (CV%)' and the
    # Methods state that BSV was applied exponentially, P_ij = theta_j x
    # exp(eta_ij). Percent CV of a log-normal is therefore back-transformed to
    # the variance as omega^2 = log(1 + CV^2).
    # BSV on Q/F, V3/F and NN is reported as 'not estimated' and is encoded as
    # absent rather than invented.
    etalcl + etalvc ~ c(
      0.01246598,
      0.01873486, 0.07444306
    )
    # Hong 2016 Table 5 OMEGA BLOCK: 'omega CLt/F' = 11.2% CV -> log(1 + 0.112^2) = 0.01246598; 'omega V2/F' = 27.8% CV -> log(1 + 0.278^2) = 0.07444306; 'rho CLt/F - V2/F' = 0.615 -> covariance 0.615 * sqrt(0.01246598 * 0.07444306) = 0.01873486
    etalmtt_fasted ~ 0.17476726
    # Hong 2016 Table 5 'omega MTTfast' = 43.7% CV -> log(1 + 0.437^2) = 0.17476726
    etalmtt_fed ~ 0.51327921
    # Hong 2016 Table 5 'omega MTTfed' = 81.9% CV -> log(1 + 0.819^2) = 0.51327921
    etalka ~ 0.21836093
    # Hong 2016 Table 5 'omega Ka' = 49.4% CV, described as 'BSV of Ka (fasting and fed)', i.e. ONE eta shared by both prandial states -> log(1 + 0.494^2) = 0.21836093

    # --- Residual unexplained variability ---
    propSd <- 0.19
    label("Proportional residual error (fraction)")
    # Hong 2016 Table 5 row 'sigma prop (%)' = 19.0 (RSE 4.76%; bootstrap median 17.7, 95% CI 15.6-19.9). A proportional error model was selected over additive and combined alternatives (Methods, 'Basic PK model').
  })

  model({
    # --- 1. Individual disposition parameters ------------------------------
    # CL/F = CLt/F x (CLCR/120)^theta_CLCR, exactly as printed in the Table 5
    # parameter-name column. Non-circular check on this line: the Discussion
    # states that 'in a typical individual with the same CLCR [107 mL/min], the
    # CL/F was predicted to be 5.89 L/h by our transit compartment model', and
    # 6.25 * (107/120)^0.511 = 5.894 L/h, which confirms the intercept, the
    # exponent and the 120 mL/min centring value simultaneously.
    cl <- exp(lcl + etalcl) * (CRCL / 120)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    # --- 2. Prandial-state-specific absorption -----------------------------
    # The source estimated MTT and ka independently in the fasted and fed
    # states (four separate THETAs) rather than as a reference value carrying a
    # multiplicative food effect, so the two states are selected literally by
    # the dose-record indicator FED. MTT additionally carries a SEPARATE eta
    # per state (43.7% vs 81.9% CV) while ka carries a SINGLE eta shared across
    # both states, which is why the ka expression has one eta outside the
    # switch and MTT is built from two separate per-state quantities.
    # Each MTT eta is kept on its own simple exp(theta + eta) line so that both
    # stay mu-referenced for nlmixr2 estimation; because FED is strictly 0 or 1
    # the linear blend below selects one branch exactly, with no interpolation.
    mtt_fasted <- exp(lmtt_fasted + etalmtt_fasted)
    mtt_fed <- exp(lmtt_fed + etalmtt_fed)
    mtt <- (1 - FED) * mtt_fasted + FED * mtt_fed
    ka <- exp((1 - FED) * lka_fasted + FED * lka_fed + etalka)
    nn <- exp(lnn)
    fdepot <- exp(lfdepot)

    # --- 3. Micro-constants ------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # --- 4. ODE system -----------------------------------------------------
    # The dose is administered to `depot` in the event table. rxode2's builtin
    # transit(nn, mtt, fdepot) returns the Savic 2007 analytical gamma-density
    # input rate, ktr * (ktr * t)^nn * exp(-ktr * t) / gamma(nn + 1) scaled by
    # fdepot * dose, with ktr = (nn + 1) / mtt computed internally; using the
    # closed form rather than an explicit chain is what permits the estimated
    # non-integer nn = 3.61. f(depot) <- 0 below suppresses the ordinary dose
    # bolus so that transit() is the sole input pathway. Disposition is the
    # standard two-compartment system in micro-constant form.
    d/dt(depot) <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # --- 5. Bioavailability ------------------------------------------------
    # Suppress the bolus; the entire dose enters through transit() instead.
    # This idiom is silently broken for some models at rxode2 5.1.7-5.1.8 (it
    # can zero the whole d/dt(depot) right-hand side and return an all-zero
    # solve with no error), so the companion vignette gates it with the
    # steady-state mass balance CL/F x AUCtau = Dose x F, which is the one
    # check that failure mode cannot pass. Verified here: 1.0000000.
    f(depot) <- 0

    # --- 6. Observation and residual error ---------------------------------
    # Dose in mg and volumes in L give Cc in mg/L, numerically identical to
    # ug/mL. The paper tabulates assay LLOQs in ng/mL (Table 3), which are
    # 1000x the model's ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
