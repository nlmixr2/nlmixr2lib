Zurawska_2026_piperacillin <- function() {
  description <- paste(
    "One-compartment population PK model for piperacillin in critically ill",
    "adults with hospital-acquired pneumonia receiving continuous-infusion",
    "piperacillin-tazobactam (Zurawska 2026). Total clearance is a piecewise",
    "sum of four arms: a creatinine-clearance-scaled renal arm and a non-renal",
    "arm, both switched off whenever renal replacement is running, plus an",
    "effluent-flow-scaled continuous-renal-replacement arm and a",
    "literature-fixed intermittent-hemodialysis arm."
  )
  reference <- paste(
    "Zurawska M, Valadez A, Harlan E, Williamson R, Scheetz MH, Neely MN,",
    "Yarnold PR, Kang M, Donnelly HK, Martinez F, Korth E, Medernach RL,",
    "Nozick SH, Hauser AR, Ozer EA, Diaz E, Misharin AV, Wunderink RG,",
    "Rhodes NJ. Pharmacokinetic-pharmacodynamic target attainment with",
    "continuous infusion piperacillin in patients admitted to the ICU with",
    "hospital-acquired pneumonia. Antimicrob Agents Chemother.",
    "2026;70(2):e01760-25. doi:10.1128/aac.01760-25."
  )
  vignette <- "Zurawska_2026_piperacillin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault estimated creatinine clearance",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "RAW Cockcroft-Gault creatinine clearance in mL/min -- NOT",
        "BSA-normalized and NOT weight-standardized. Zurawska 2026 Simulations",
        "cites Cockcroft-Gault (reference 21) and the covariate enters Eq. 5 as",
        "the through-origin ratio (CRCL / 120 mL/min), so the renal clearance",
        "estimate is the value at a reference CrCL of 120 mL/min. Among the 19",
        "patients not requiring renal replacement the mean +/- SD was 78 +/- 68",
        "mL/min (range 9-229; Results, first paragraph); Table 1 reports the",
        "initial value as 77.7 +/- 68.4 mL/min. Cockcroft-Gault CrCL is",
        "conventionally not defined for renal-replacement-dependent subjects,",
        "and Eq. 5 gates the whole renal-plus-non-renal bracket off whenever",
        "either RRT arm is running, so the value supplied for those records",
        "does not affect the prediction. The paper's simulations exercised the",
        "range 25-150 mL/min."
      ),
      source_name        = "CRCL"
    ),
    RRT_CRRT_ACTIVE = list(
      description        = "Continuous renal replacement therapy currently running",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CRRT running)",
      notes              = paste(
        "Zurawska 2026 Eqs. 5-6 indicator CRRT; 'HD and CRRT are indicator",
        "variables (either on or off)' (Covariate model). Modeled as a",
        "time-varying regressor in Monolix ('we chose to use the regressor",
        "approach available in Monolix to model time-varying covariates'),",
        "hence the ACTIVE rather than the STATUS canonical. 15 of 35 patients",
        "(43%) required CRRT (Table 1)."
      ),
      source_name        = "CRRT"
    ),
    RRT_HEMODIAL_ACTIVE = list(
      description        = "Intermittent hemodialysis session currently running",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no hemodialysis session running)",
      notes              = paste(
        "Zurawska 2026 Eqs. 5 and 7 indicator HD. Time-varying regressor, so",
        "the ACTIVE rather than the STATUS canonical. Only 1 of 35 patients",
        "(3%) required intermittent HD (Table 1), which is why the intra-HD",
        "clearance could not be estimated and was fixed from the literature:",
        "'Given the opportunistic nature of blood sampling (typically collected",
        "in the morning) relative to hemodialysis (HD) times (typically",
        "conducted in the late afternoon), HD clearance was fixed to literature",
        "values' (Population PK modeling)."
      ),
      source_name        = "HD"
    ),
    RRT_CRRT_EFFLUENT_FLOW = list(
      description        = "Total effluent flow rate of the CRRT circuit",
      units              = "mL/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Zurawska 2026 Eq. 6 covariate FLOW, entering as the through-origin",
        "ratio (FLOW / 2 L/h), so the CRRT clearance estimate is the value at a",
        "2 L/h effluent flow -- 'Our population mean CRRT clearance estimate",
        "was 3.4 L/h standardized to a 2 L/hr effluent flow rate' (Discussion).",
        "The prose after Eq. 6 states FLOW is 'in mL/h'; this file keeps the",
        "canonical mL/h units and divides by 2000 inside model() so both the",
        "printed prose and the printed reference constant are honoured.",
        "Observed cohort effluent among CRRT patients: highest mean +/- SD",
        "2.75 +/- 1.1 L/h, i.e. 32 +/- 7.8 mL/kg/h (Results, first paragraph).",
        "Weight-normalized prescriptions must be multiplied by body weight",
        "before being supplied here: the paper's simulations used 25 and 35",
        "mL/kg/h at an assumed 70 kg, i.e. 1750 and 2450 mL/h, with a",
        "sensitivity analysis at 91 and 126 kg. Meaningful only while",
        "RRT_CRRT_ACTIVE = 1."
      ),
      source_name        = "FLOW"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened as a planned regressor on Vd and NOT retained. Table S1 Run 8",
        "(V scaled to WT/70) raised BICc from 2495.46 to 2500.73 and Run 9",
        "(allometric WT/70) raised it to 2499.69 with 'high RSE on all",
        "parameters'; the Discussion states 'we did not observe a significant",
        "effect of body weight on Vd'. Cohort 79.6 +/- 24.1 kg (Table 1). Body",
        "weight still matters operationally, because CRRT effluent is",
        "prescribed in mL/kg/h -- but it enters through",
        "RRT_CRRT_EFFLUENT_FLOW, not as a covariate in its own right."
      )
    ),
    AGE = list(
      description = "Age",
      units       = "year",
      type        = "continuous",
      notes       = paste(
        "Planned covariate, screened and not retained (Covariate, regressor,",
        "and error models). Cohort 62 +/- 16 years (Table 1). No coefficient is",
        "reported anywhere in the paper."
      )
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = paste(
        "Planned covariate, screened and not retained. Cohort 1.9 +/- 0.3 m^2",
        "(Table 1). No coefficient is reported anywhere in the paper."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Planned categorical covariate, screened and not retained; Eq. 4 shows",
        "the fractional-change form that was evaluated. Cohort 18 female (51%)",
        "/ 17 male (49%) (Table 1). No coefficient is reported anywhere in the",
        "paper."
      )
    )
  )

  compartmentData <- list(
    central = list(
      analyte  = "piperacillin",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 35,
    n_studies      = 2,
    n_observations = 162,
    age_mean_sd    = "62 +/- 16 years",
    weight_mean_sd = "79.6 +/- 24.1 kg",
    bsa_mean_sd    = "1.9 +/- 0.3 m^2",
    sex_female_pct = 51,
    disease_state  = "hospital-acquired pneumonia requiring intensive care",
    renal_function = paste(
      "Wide range. Among the 19 patients not requiring renal replacement,",
      "Cockcroft-Gault CrCL was 78 +/- 68 mL/min (range 9-229). 16 of 35 (46%)",
      "required renal replacement: 15 (43%) CRRT and 1 (3%) intermittent",
      "hemodialysis. Highest mean CRRT effluent flow 2.75 +/- 1.1 L/h",
      "(32 +/- 7.8 mL/kg/h). Initial serum creatinine 1.7 +/- 1.1 mg/dL in",
      "non-RRT patients."
    ),
    dose_range     = paste(
      "Institutional-protocol piperacillin-tazobactam dosing at the discretion",
      "of the treating team (asp.nm.org), with indication-based selection and",
      "renal dose adjustment. The model was then applied to simulated",
      "continuous-infusion regimens of 3-12 g/day of piperacillin after a 4 g",
      "loading dose."
    ),
    regions        = "United States (single center, Chicago, Illinois)",
    notes          = paste(
      "Baseline demographics per Zurawska 2026 Table 1. Critically ill adults",
      "with hospital-acquired pneumonia admitted to the medical ICU at",
      "Northwestern Memorial Hospital and treated with piperacillin-tazobactam,",
      "enrolled June 2018 - May 2023 in the SCRIPT (U19AI135964) and RX-SCRIPT",
      "(R21AI174159, R01AI158530) studies. Residual clinical blood samples were",
      "captured opportunistically, so sampling is sparse and irregular relative",
      "to dosing and to hemodialysis sessions. Patients requiring concurrent",
      "ECMO were excluded. 162 plasma samples, of which 7 were below the 0.5",
      "mg/L LLOQ and were handled by interval censoring over [0, LLOQ] rather",
      "than being discarded. Total (not free) piperacillin was quantified by a",
      "validated LC-MS assay linear from 0.5 to 80 mg/L. The companion",
      "tazobactam population PK model is reported separately by Williamson et",
      "al. (doi:10.1128/aac.01766-25) and is not part of this model."
    )
  )

  ini({
    # Structural parameters -- Zurawska 2026 Table 2 (Fixed effects), with the
    # renal arm taken from the Results narrative; see the note on lcl_renal.
    lvc <- log(28.1)
    label("Central volume of distribution (L)")  # Table 2, Vd = 28.1 L

    # Zurawska 2026 reports TWO different values for the renal clearance arm.
    # Table 2 lists "CL_CrCL, L/h (reference, 120 mL/min) = 5.1", while the
    # Final population PK model parameters paragraph states "The population
    # typical estimates for PIP Vd, CL_R, and CL_NR were 28.1, 5.77, and 0.65
    # L/h, respectively" -- Vd and CL_NR agree with Table 2 exactly, only the
    # renal arm disagrees. The paper's own downstream simulation settles it in
    # favour of 5.77: replicating Fig. 3 and Table S2 (12 published regimens,
    # each scored on P(C48 < 16), P(C48 > 96) and P(C48 > 160)) reproduces the
    # published percentages with 5.77 and badly overshoots the excessive-
    # exposure tail with 5.1. The sum-of-squared-error against the published
    # percentages is ~5x smaller with 5.77, and the ranking is unchanged
    # whether or not parameter uncertainty, residual error, or the alternative
    # CV-to-omega convention is included. The vignette carries the full
    # replication. 5.77 is used here; 5.1 remains inside the paper's own
    # bootstrap 95% CI of 3.0-8.0, as does 5.77.
    lcl_renal <- log(5.77)
    label("Renal clearance at a Cockcroft-Gault CrCL of 120 mL/min (L/h)")  # Results, Final population PK model parameters

    lcl_nonren <- log(0.65)
    label("Non-renal clearance (L/h)")  # Table 2, Non-renal CL = 0.65 L/h

    lcl_crrt <- log(3.4)
    label("Continuous-renal-replacement clearance at a 2 L/h effluent flow (L/h)")  # Table 2, CL_CRRT = 3.4 L/h

    # Fixed, not estimated: "Intra-HD clearance for PIP was taken from the
    # literature and fixed at 4.16 L/h (19)" (Covariate model). Table 2 shows
    # this row with "Fixed" in every uncertainty column and rounds the estimate
    # to 4.2; the unrounded 4.16 from the prose is used.
    lcl_hemodialysis <- fixed(log(4.16))
    label("Intra-hemodialysis clearance (L/h)")  # Covariate model, fixed from reference 19

    # IIV -- Table 2 reports interindividual variability as CV%, so the
    # variances below are omega^2 = log(1 + CV^2).
    # Vd:       CV 32.9% -> log(1 + 0.329^2) = 0.1028
    # CL_CrCL:  CV 69.6% -> log(1 + 0.696^2) = 0.3950
    # CL_CRRT:  CV 24.5% -> log(1 + 0.245^2) = 0.0583
    etalvc ~ 0.1028        # Table 2, omega_Vd = 32.9% CV
    etalcl_renal ~ 0.3950  # Table 2, omega_CLCrCL = 69.6% CV
    etalcl_crrt ~ 0.0583   # Table 2, omega_CLCRRT = 24.5% CV
    # No IIV on the non-renal arm: "high relative standard errors and CV% for
    # inter-individual variability (IIV) were observed with CL_NR. Removal of
    # IIV on CL_NR (Run 5) retained OFV improvements (-25.1) yielding an
    # acceptable condition number. Run 5 was retained as the final model."
    # No IIV on the hemodialysis arm, which is fixed entirely.

    # Residual error -- Table S1 selects "y1 prop" for piperacillin in every
    # retained run; Run 6 (combined) and Run 7 both raised BICc.
    propSd <- 0.25
    label("Proportional residual error (fraction)")  # Table 2, proportional error 25%
  })

  model({
    # 1. Individual clearance arms, Zurawska 2026 Eqs. 5-7.
    #
    #    Eq 5: CL_nonCRRT,ind = (CL_CRCL,pop * (CRCL / 120 mL/min) + CL_NR,pop)
    #                           * (1 - CRRT) * (1 - HD)
    #    Eq 6: CL_CRRT,ind    = CL_CRRT,pop * (FLOW / 2 L/h) * CRRT
    #    Eq 7: CL_HD,ind      = CL_HD,pop * HD
    #
    #    Eq. 3 supplies the exp(eta) placement; the random effect sits on the
    #    population parameter of each arm, not on the summed clearance, which
    #    is why Table 2 names the random effects omega_CLCrCL and omega_CLCRRT
    #    (per-arm) and why IIV could be removed from CL_NR alone in Run 5.
    #    The CRCL ratio in Eq. 5 carries no exponent -- it is a through-origin
    #    linear scaling, not the power model of Eq. 3.
    cl_renal <- exp(lcl_renal + etalcl_renal) * (CRCL / 120)
    cl_nonren <- exp(lcl_nonren)
    # The effluent-flow ratio is (FLOW / 2 L/h) with FLOW in mL/h, hence 2000.
    cl_crrt <- exp(lcl_crrt + etalcl_crrt) *
      (RRT_CRRT_EFFLUENT_FLOW / 2000) * RRT_CRRT_ACTIVE
    cl_hemodialysis <- exp(lcl_hemodialysis) * RRT_HEMODIAL_ACTIVE

    # 2. Total clearance, Eq. 1: CL_total = CL_nonCRRT + CL_CRRT + CL_HD.
    #    Note that Eq. 5 gates the ENTIRE renal-plus-non-renal bracket off
    #    while either renal-replacement modality is running -- the non-renal
    #    arm is switched off too, not just the renal one. This is the printed
    #    equation and it is also what the paper's simulations behave like:
    #    retaining the non-renal arm during CRRT drops the simulated
    #    probability of exceeding 96 mg/L on the high-dose CRRT regimens from
    #    ~30% to ~6-17%, against published values of roughly 24% and 32%.
    cl <- (cl_renal + cl_nonren) *
      (1 - RRT_CRRT_ACTIVE) * (1 - RRT_HEMODIAL_ACTIVE) +
      cl_crrt + cl_hemodialysis

    vc <- exp(lvc + etalvc)

    # 3. One-compartment disposition with first-order elimination from the
    #    central compartment. A two-compartment base model (Table S1 Run 3)
    #    improved the objective function but was rejected: "the sparse nature
    #    of opportunistic plasma sampling ... led to unstable estimates of
    #    intercompartmental clearance and peripheral volume."
    kel <- cl / vc

    d/dt(central) <- -kel * central

    # 4. Observation. Total (not free) piperacillin in plasma.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
