Williamson_2026_tazobactam <- function() {
  description <- paste(
    "One-compartment population PK model for tazobactam (given with",
    "piperacillin) in critically ill adults with hospital-acquired pneumonia",
    "(Williamson 2026). Total clearance is a piecewise sum of four arms: a",
    "creatinine-clearance-scaled renal arm and a non-renal arm, both switched",
    "off whenever renal replacement is running, plus an effluent-flow-scaled",
    "continuous-renal-replacement arm and a literature-fixed intermittent-",
    "hemodialysis arm. Unbound concentration uses a fixed 70% unbound fraction."
  )
  reference <- paste(
    "Williamson R, Zurawska M, Valadez A, Harlan E, Scheetz MH, Neely MN,",
    "Yarnold PR, Kang M, Donnelly HK, Martinez F, Korth E, Medernach RL,",
    "Nozick SH, Hauser AR, Ozer EA, Diaz E, Misharin AV, Wunderink RG,",
    "Rhodes NJ. Critical threshold target attainment rates for tazobactam",
    "combined with piperacillin among patients admitted to the ICU with",
    "hospital-acquired pneumonia. Antimicrob Agents Chemother.",
    "2026;70(2):e01766-25. doi:10.1128/aac.01766-25."
  )
  vignette <- "Williamson_2026_tazobactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description = "Cockcroft-Gault estimated creatinine clearance",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "RAW Cockcroft-Gault creatinine clearance in mL/min -- NOT",
        "BSA-normalized and NOT weight-standardized (Table 1 footnote:",
        "'CrCL, Cockcroft Gault estimated creatinine clearance'). Enters Eq. 5",
        "as the through-origin ratio (CrCL / 120 mL/min), so the renal",
        "clearance estimate is the value at a CrCL of 120 mL/min. Time-varying",
        "regressor in Monolix. Among patients not requiring renal replacement",
        "the mean +/- SD was 78 +/- 68 mL/min (range 9-229; Results, first",
        "paragraph). Eq. 5 gates the whole renal-plus-non-renal bracket off",
        "whenever either renal-replacement arm is running, so the value supplied",
        "for those records does not affect the prediction. The paper's",
        "simulations used 25, 50, 75 and 150 mL/min."
      ),
      source_name = "CrCL"
    ),
    RRT_CRRT_ACTIVE = list(
      description = "Continuous renal replacement therapy currently running",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no CRRT running)",
      notes = paste(
        "Williamson 2026 Eqs. 5-6 indicator CRRT; 'HD and CRRT are indicator",
        "variables (either on or off)' (Covariate model). Modeled as a",
        "time-varying regressor in Monolix, hence the ACTIVE rather than the",
        "STATUS canonical. 15 of 35 patients required CRRT (Results, first",
        "paragraph)."
      ),
      source_name = "CRRT"
    ),
    RRT_HEMODIAL_ACTIVE = list(
      description = "Intermittent hemodialysis session currently running",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no hemodialysis session running)",
      notes = paste(
        "Williamson 2026 Eqs. 5 and 7 indicator HD. Time-varying regressor, so",
        "the ACTIVE rather than the STATUS canonical. Only 1 of 35 patients",
        "required intermittent HD (Results, first paragraph); 'Given the",
        "opportunistic nature of blood sampling (typically collected in the",
        "morning) relative to hemodialysis (HD) times (typically conducted in",
        "the late afternoon), TAZ HD clearance was fixed to literature values'",
        "(Population PK modeling)."
      ),
      source_name = "HD"
    ),
    RRT_CRRT_EFFLUENT_FLOW = list(
      description = "Total effluent flow rate of the CRRT circuit",
      units = "mL/h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Williamson 2026 Eq. 6 covariate FLOW, entering as the through-origin",
        "ratio (FLOW / 2 L/hr). The prose after Eq. 6 states FLOW is 'in mL/h';",
        "this file keeps the canonical mL/h units and divides by 2000 inside",
        "model(). Observed cohort effluent among CRRT patients: highest mean",
        "+/- SD 2.75 +/- 1.1 L/h, i.e. 32 +/- 7.8 mL/kg/h (Results, first",
        "paragraph). Weight-normalized prescriptions must be multiplied by body",
        "weight before being supplied here: the paper's simulations used 25 and",
        "35 mL/kg/h at an assumed 70 kg (1750 and 2450 mL/h), with a sensitivity",
        "analysis at 91 and 126 kg. Meaningful only while RRT_CRRT_ACTIVE = 1."
      ),
      source_name = "FLOW"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Planned regressor on Vd and NOT retained. Table S1 Run 8 (V scaled to",
        "WT/70) raised BICc from 2495.46 to 2500.73 and Run 9 (allometric",
        "WT/70) raised it to 2499.69 with 'high RSE on all parameters'. Cohort",
        "mean 79.6 kg (Results, first paragraph). Body weight enters only",
        "through the weight-normalized CRRT effluent prescription, i.e. via",
        "RRT_CRRT_EFFLUENT_FLOW."
      )
    ),
    AGE = list(
      description = "Age",
      units = "year",
      type = "continuous",
      notes = paste(
        "Planned covariate, screened and not retained (Covariate, regressor,",
        "and error models). Cohort mean 62 years. No coefficient reported."
      )
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Planned covariate, screened and not retained. No coefficient reported."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Planned categorical covariate, screened and not retained; Eq. 4 shows",
        "the fractional-change form evaluated. Cohort 51% female. No",
        "coefficient reported."
      )
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "tazobactam",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 35,
    n_studies = 2,
    n_observations = 162,
    age_mean = "62 years",
    weight_mean = "79.6 kg",
    sex_female_pct = 51,
    disease_state = "hospital-acquired pneumonia requiring intensive care",
    renal_function = paste(
      "Wide range. Among patients not requiring renal replacement,",
      "Cockcroft-Gault CrCL was 78 +/- 68 mL/min (range 9-229). 16 of 35",
      "required renal replacement: 15 CRRT and 1 intermittent hemodialysis.",
      "Highest mean CRRT effluent flow 2.75 +/- 1.1 L/h (32 +/- 7.8 mL/kg/h)."
    ),
    dose_range = paste(
      "Institutional-protocol piperacillin-tazobactam dosing (asp.nm.org) with",
      "indication-based selection and renal dose adjustment. The model was",
      "then applied to simulated tazobactam regimens after a 500 mg loading",
      "dose: continuous infusion 375-1500 mg/day, extended infusion 375-500 mg",
      "every 8-12 h over 4 h, and intermittent infusion 375-500 mg every 6-8 h",
      "over 30 min (Table S2)."
    ),
    regions = "United States (single center, Chicago, Illinois)",
    notes = paste(
      "Critically ill adults with hospital-acquired pneumonia admitted to the",
      "medical ICU at Northwestern Memorial Hospital and treated with",
      "piperacillin-tazobactam, enrolled June 2018 - May 2023 in the SCRIPT",
      "and RX-SCRIPT studies; full demographics are in the companion",
      "piperacillin paper (Zurawska 2026, doi:10.1128/aac.01760-25), which",
      "used the same cohort and the same joint Monolix model-building sequence",
      "(Table S1). Residual clinical blood samples were captured",
      "opportunistically. Patients requiring concurrent ECMO were excluded.",
      "162 plasma samples, 4 below the 0.5 mg/L LLOQ (interval-censored over",
      "[0, LLOQ] in Monolix). Total tazobactam was quantified by a validated",
      "LC-MS/MS assay linear from 0.5 to 80 mg/L."
    )
  )

  ini({
    # Structural parameters -- Williamson 2026 Table 1 (Fixed effects); the
    # Results paragraph 'Final population PK model parameters' repeats the same
    # values (32.8 L, 4.9 L/h, 0.91 L/h, and 4 L/h on CRRT).
    lvc <- log(32.8)
    label("Central volume of distribution (L)") # Table 1, Vd = 32.8 L

    lcl_renal <- log(4.9)
    label("Renal clearance at a Cockcroft-Gault CrCL of 120 mL/min (L/h)") # Table 1, CL_CrCL = 4.9 L/h/120 mL/min

    lcl_nonren <- log(0.91)
    label("Non-renal clearance (L/h)") # Table 1, Non-renal CL = 0.91 L/h

    lcl_crrt <- log(4.0)
    label("Continuous-renal-replacement clearance at a 2 L/h effluent flow (L/h)") # Table 1, CL_CRRT = 4.0 L/h

    # Not estimated: 'Intra-HD clearance for TAZ was taken from the literature
    # and fixed at 5.68 L/h (22)' (Covariate model); Table 1 marks the row
    # 'Fixed' in every uncertainty column.
    lcl_hemodialysis <- fixed(log(5.68))
    label("Intra-hemodialysis clearance (L/h)") # Table 1, CL_HD = 5.68 L/h, fixed

    # IIV -- Table 1 reports interindividual variability as 'C.V. %' (footnote:
    # 'reported as CV % based on the standard deviation of the random
    # effects'), converted here with omega^2 = log(1 + CV^2), the convention
    # of the companion piperacillin model (Zurawska_2026_piperacillin). Reading
    # the percentages instead as omega SDs (omega^2 = CV^2) widens only the
    # renal arm materially (1.6129 vs 0.9605) and moves the simulated target
    # attainment of the paper's continuous-infusion arms further from Fig. 3;
    # the vignette shows both.
    etalvc ~ log(1 + 0.190^2) # Table 1, omega_Vd = 19.0% CV
    etalcl_renal ~ log(1 + 1.27^2) # Table 1, omega_CLCrCL = 127% CV
    etalcl_crrt ~ log(1 + 0.261^2) # Table 1, omega_CLCRRT = 26.1% CV
    # No IIV on the non-renal arm (Table S1 Run 5, 'no IIV on CLNR', is the
    # final model) and none on the fixed hemodialysis arm.

    # Residual error -- Table S1 selects 'y2 comb1' for tazobactam in the final
    # run: Monolix combined1, SD = a + b * f (linear sum of the two SDs).
    addSd <- 1.6
    label("Additive residual error (mg/L)") # Table 1, Additive = 1.6 mg/L
    propSd <- 0.064
    label("Proportional residual error (fraction)") # Table 1, Proportional error = 6.4%

    # Unbound fraction. Not estimated: the literature value the authors applied
    # to convert total to unbound concentrations for the target-attainment
    # simulations -- 'Simulations evaluated free (i.e., unbound) TAZ
    # concentrations using published protein binding correction of 70% (12)',
    # with the Discussion clarifying 'assumed that 70% of TAZ remains unbound'.
    fu <- fixed(0.7)
    label("Fraction of tazobactam unbound in plasma (unitless)") # Simulations, 70% unbound (reference 12)
  })

  model({
    # 1. Individual clearance arms, Williamson 2026 Eqs. 5-7:
    #
    #    Eq 5: CL_nonCRRT,ind = (CL_CrCL,pop * (CrCL / 120 mL/min) + CL_NR,pop)
    #                           * (1 - CRRT) * (1 - HD)
    #    Eq 6: CL_CRRT,ind    = CL_CRRT,pop * (FLOW / 2 L/hr) * CRRT
    #    Eq 7: CL_HD,ind      = CL_HD,pop * HD
    #
    #    Eq. 3 places exp(eta) on each population parameter; Table 1 names the
    #    random effects per arm (omega_CLCrCL, omega_CLCRRT).
    cl_renal <- exp(lcl_renal + etalcl_renal) * (CRCL / 120)
    cl_nonren <- exp(lcl_nonren)
    # FLOW is in mL/h and the reference is 2 L/h, hence 2000.
    cl_crrt <- exp(lcl_crrt + etalcl_crrt) *
      (RRT_CRRT_EFFLUENT_FLOW / 2000) * RRT_CRRT_ACTIVE
    cl_hemodialysis <- exp(lcl_hemodialysis) * RRT_HEMODIAL_ACTIVE

    # 2. Total clearance, Eq. 1: CL_total = CL_nonCRRT + CL_CRRT + CL_HD. The
    #    Eq. 5 gate switches off the whole renal-plus-non-renal bracket while
    #    either renal-replacement modality runs ('equal to CRRT or HD CL when on
    #    renal replacement', Base model development).
    cl <- (cl_renal + cl_nonren) *
      (1 - RRT_CRRT_ACTIVE) * (1 - RRT_HEMODIAL_ACTIVE) +
      cl_crrt + cl_hemodialysis

    vc <- exp(lvc + etalvc)

    # 3. One-compartment disposition, first-order elimination from central.
    #    The two-compartment base model (Table S1 Run 3) was rejected: 'Unable
    #    to estimate peripheral V'.
    kel <- cl / vc

    d/dt(central) <- -kel * central

    # 4. Observation. The assay measured TOTAL plasma tazobactam, so Cc is the
    #    total concentration; Cu is the unbound concentration on which the
    #    paper's 100% fT > threshold targets are evaluated.
    Cc <- central / vc
    Cu <- fu * Cc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
