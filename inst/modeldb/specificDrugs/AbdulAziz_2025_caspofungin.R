AbdulAziz_2025_caspofungin <- function() {
  description <- "One-compartment population PK model with first-order elimination for intravenous caspofungin in critically ill adults receiving extracorporeal membrane oxygenation (ECMO), from the ASAP ECMO study program. Central volume scales with Janmahasatian fat-free mass (reported by the source as lean body weight) through an estimated power exponent of 1.24 standardised to the cohort median of 58.9 kg; clearance carries no covariate. Every patient in the cohort was cannulated onto ECMO, so ECMO support is a property of the population rather than a covariate -- ECMO duration, mode and flow rate were all screened and none was retained, which is the paper's central negative finding. Estimated by SAEM in Monolix 2021R2 with an additive residual error. Abdul-Aziz 2025, n = 8 patients, 64 plasma samples over a single dosing interval."
  reference <- "Abdul-Aziz MH, Diehl A, Liu X, Cheng V, Corley A, Gilder E, Levkovich B, McGuinness S, Ordonez J, Parke R, Pellegrino V, Wallis SC, Fraser JF, Shekar K, Roberts JA. Population pharmacokinetics of caspofungin in critically ill patients receiving extracorporeal membrane oxygenation-an ASAP ECMO study. Antimicrob Agents Chemother. 2025;69(2):e01435-24. doi:10.1128/aac.01435-24. PMC11823646. Published online 18 December 2024."
  vignette <- "AbdulAziz_2025_caspofungin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Caspofungin was given as an intermittent intravenous
  # infusion and TOTAL plasma caspofungin was assayed by UHPLC-MS/MS with
  # dicloxacillin as the internal standard (Abdul-Aziz 2025, "Caspofungin
  # assay"); the model therefore predicts total, not unbound, concentration.
  compartmentData <- list(
    central = list(analyte = "caspofungin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass computed with the Janmahasatian et al. equation from body weight, height and sex; reported by the source under the label 'lean body weight' (LBW)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The ONLY covariate retained in the final model, and only on the",
        "central volume. Enters as the power term (FFM / 58.9)^1.24; the",
        "58.9 kg standardisation is the cohort median lean body weight",
        "(Table 1: 58.9 kg, range 37.3-81), consistent with the Methods",
        "statement that 'continuous covariates were centered on their median",
        "values'. Registered as FFM rather than LBM because the register",
        "discriminates the two by the estimating FORMULA, not by the paper's",
        "label, and Abdul-Aziz 2025 cites Janmahasatian et al. (Clin",
        "Pharmacokinet 2005;44:1051-1065) -- the fat-free-mass equation. The",
        "source itself makes the identification explicit in the Discussion:",
        "'Although the Janmahasatian et al. equation was developed to",
        "estimate fat-free mass, the estimate is considered to be a",
        "representation of lean body weight and is commonly used",
        "interchangeably.' Same quantity in kg with no value transformation,",
        "so source_name is LBW. This follows the identical LBW -> FFM",
        "precedent set by Rolsma_2025_cefepime.R. Time-fixed: the",
        "anthropometric data were collected once, on the day of",
        "pharmacokinetic sampling (Methods, 'Study procedures'). Must be",
        "strictly positive -- the covariate enters a power term. Note the",
        "exponent 1.24 is ABOVE the theory-based allometric value of 1 for a",
        "volume, and its 95% bootstrap CI is very wide (0.08-1.91), so the",
        "scaling is poorly identified in a cohort of only 8 patients; see the",
        "vignette Assumptions and deviations.",
        "Sex and serum albumin also survived initial screening (sex and LBW",
        "were tested on volume, albumin on clearance) but only LBW on volume",
        "improved the fit (Results, 'Population pharmacokinetic model",
        "building')."
      ),
      source_name        = "LBW"
    )
  )

  # Covariates the source SCREENED but did not retain in the final model.
  # Documentation only -- none of these is referenced in model(). Recording
  # them preserves the provenance of the paper's covariate screen, which
  # matters unusually much here because the absence of any ECMO effect is the
  # paper's headline conclusion. The screened set is listed in Methods,
  # "Covariate screening and model development".
  #
  # Four screened covariates are deliberately NOT given entries below because
  # no canonical column exists for them and inventing one for a covariate that
  # no model actually uses would be an unwarranted register addition:
  # adjusted body weight, the SOFA score, and the ECMO treatment variables
  # (duration, mode, flow rate). They are described in the vignette narrative
  # and in population$notes instead. Note also that ECMO_STATUS would be
  # meaningless here: all 8 patients were on ECMO, so the indicator is
  # constant at 1 and cannot be estimated as a covariate.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened, not retained. Cohort median 36.5 years (range 20.0-62.0; Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened and one of the three covariates that survived initial",
        "screening -- sex was tested on the volume of distribution alongside",
        "lean body weight -- but it did not improve the model fit and was not",
        "retained (Results, 'Population pharmacokinetic model building').",
        "Cohort 5 of 8 male (62.5%), i.e. 37.5% female (Table 1). Note the",
        "source reports male count, so the canonical SEXF requires the value",
        "transformation SEXF = 1 - SEXM; immaterial here because the term is",
        "absent from the final model."
      )
    ),
    IBW = list(
      description = "Ideal body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened, not retained. One of three body-weight descriptors tested against each other; fat-free mass won. Cohort median 68.0 kg (range 55.0-80.0; Table 1)."
    ),
    WT = list(
      description = "Total (actual) body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened, not retained -- the source's 'body weight measures' screen preferred fat-free mass over total body weight. Cohort median 84.0 kg (range 55.0-130.0; Table 1). The Discussion notes that earlier caspofungin studies found total body weight and fat-free mass predictive, and frames the fat-free-mass result as this paper's novel contribution."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened, not retained. Cohort median 25.8 kg/m^2 (range 19.5-38; Table 1)."
    ),
    CRCL = list(
      description = "Creatinine clearance estimated with the Cockcroft-Gault equation",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Screened, not retained. Cohort median 57.0 mL/min (range 31.0-260.0;",
        "Table 1) -- a span from renal impairment to markedly augmented renal",
        "clearance. Reported as RAW Cockcroft-Gault mL/min, NOT BSA-normalized",
        "to the canonical CRCL units of mL/min/1.73 m^2; a user who retained",
        "this covariate would have to normalize first. Recorded with the raw",
        "units here because the entry is documentation only and never reaches",
        "model()."
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = paste(
        "Screened and one of the three covariates that survived initial",
        "screening -- albumin was tested on clearance -- but it did not",
        "improve the fit and was not retained (Results, 'Population",
        "pharmacokinetic model building'). Cohort median 24.0 g/L (range",
        "17.0-33.0; Table 1). The authors flag this as a discrepancy from",
        "prior caspofungin work and explain it in the limitations: 'This could",
        "be due to small variations in albumin concentrations as all patients",
        "demonstrated some degree of hypoalbuminemia (range 17-33 g/L).'",
        "Already in SI units (g/L), so no g/dL conversion is needed."
      )
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened, not retained. Cohort median 9.8 mmol/L (range 3.4-24.3; Table 1)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened, not retained. Cohort median 12.0 umol/L (range 10.0-27.0; Table 1). Patients with bilirubin > 150 umol/L were excluded at enrolment, so the tested range is narrow by design."
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation II score at ICU admission",
      units       = "(score)",
      type        = "continuous",
      notes       = "Screened, not retained. Cohort median 16.5 on admission (range 10-32; Table 1)."
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous / extended renal-replacement-therapy treatment-status indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened as 'receipt of RRT', not retained. 5 of 8 patients (62.5%)",
        "received concomitant renal replacement therapy during",
        "pharmacokinetic sampling -- 4 on CVVHD and 1 on SLED (Table 1), both",
        "continuous / extended modalities, hence RRT_CRRT_STATUS rather than",
        "RRT_HEMODIAL_STATUS. Unlike Valadez_2025_cefepime.R, which excluded",
        "renal-replacement patients specifically to isolate an ECMO effect,",
        "this cohort is majority-RRT, so any ECMO effect here is confounded",
        "with RRT; with n = 8 the two cannot be separated."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 8L,
    n_studies      = 1L,
    n_samples      = 64L,
    age_median     = "36.5 years (range 20.0-62.0; Table 1). Eligibility 18-90 years.",
    weight_median  = "84.0 kg actual body weight (range 55.0-130.0; Table 1)",
    ffm_median     = "58.9 kg lean body weight / Janmahasatian fat-free mass (range 37.3-81; Table 1). This is the 58.9 kg standardisation in the volume covariate model.",
    bmi_median     = "25.8 kg/m^2 (range 19.5-38; Table 1)",
    sex_female_pct = 37.5,
    race_ethnicity = "Not reported. Table 1 tabulates age, weight descriptors, BMI, sex, severity scores, laboratory values, ECMO indication / mode / flow and renal replacement therapy only.",
    disease_state  = paste(
      "Critically ill adults in the intensive care unit receiving caspofungin",
      "while undergoing ECMO therapy for cardiac and/or respiratory failure.",
      "ECMO indication: acute respiratory distress syndrome 3 (37.5%), lung",
      "transplant 3 (37.5%), cardiac arrest 1 (12.5%), cardiogenic shock 1",
      "(12.5%). Severity APACHE II 16.5 on admission (range 10-32) and SOFA 8.0",
      "on the sampling day (range 5.0-15.0). All patients were",
      "hypoalbuminaemic (albumin 17-33 g/L)."
    ),
    ecmo_support   = paste(
      "ALL 8 patients were cannulated onto ECMO -- veno-venous in 6 (75.0%)",
      "and veno-arterial in 2 (25.0%), median flow rate 3.5 (range 2.9-5.7).",
      "Median time to pharmacokinetic sampling after ECMO initiation was 3.5",
      "days (range 1-30). ECMO support is therefore a property of the",
      "population, not a covariate: the indicator is constant at 1 and cannot",
      "be estimated. The ECMO treatment variables that DO vary (duration, mode,",
      "flow rate) were all screened and none was retained, which the authors",
      "read as evidence that 'the use of ECMO had negligible impact on the",
      "pharmacokinetics of caspofungin in critically ill patients'. Because 5",
      "of 8 patients were also on renal replacement therapy, and because there",
      "is no non-ECMO control arm (a limitation the authors state), this is a",
      "null finding within an all-ECMO cohort rather than a quantified",
      "ECMO-vs-no-ECMO contrast."
    ),
    renal_function = "Cockcroft-Gault creatinine clearance median 57.0 mL/min (range 31.0-260.0). 5 of 8 (62.5%) on concomitant renal replacement therapy during sampling (4 CVVHD, 1 SLED). Neither creatinine clearance nor receipt of RRT was retained as a covariate.",
    dose_range     = "All patients received a 70 mg loading dose on day 1 as a 1 h intravenous infusion, followed by a daily maintenance dose of 50 mg (6 patients) or 70 mg (2 patients). Monte Carlo dosing simulations additionally explored loading doses of 50-200 mg with maintenance doses of 50, 70 or 100 mg daily.",
    regions        = "Australia, New Zealand, South Korea and Switzerland (six ICUs; ASAP ECMO study program, November 2012 - November 2019)",
    sampling       = "Serial arterial-line sampling over a SINGLE dosing interval, once the patient had been stabilised on ECMO: 0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 6, 8, 12 and 24 h after commencement of the infusion. 64 total concentration-time points from 8 patients. Assay limit 0.1 mg/L, linear 0.1-20 mg/L, precision and accuracy within 12%.",
    protein_binding = "Not fitted. Caspofungin is highly protein-bound, which is the mechanism by which the authors hypothesised ECMO circuit sequestration; TOTAL plasma concentrations were assayed and the model predicts total drug. The AUC/MIC targets the paper simulates against are likewise defined on total drug ('the target ratio of total drug area under the concentration-time curve').",
    notes          = paste(
      "Structural model selection: one- and two-compartment models with",
      "first-order elimination were both fitted; the two-compartment model was",
      "rejected because it did not improve the fit (no decrease in BICc).",
      "Additive, proportional and combined residual-error models were tested",
      "and the additive form was selected. Estimated by SAEM in Monolix",
      "2021R2. Final model evaluated by goodness-of-fit plots (Fig. 1), a",
      "500-patient visual predictive check (Fig. 2), and a 1,000-run bootstrap",
      "in Rsmlx 4.0.2 / R 4.1.3 whose medians closely matched the point",
      "estimates (Table 2). Observed cohort exposure: median AUC(0-24) 115",
      "mg.h/L (range 81-151) per the Discussion. The dosing-simulation layer",
      "(probability and fractional target attainment against Candida albicans,",
      "C. glabrata and C. parapsilosis, Figs 3-5 and Table 3) is an",
      "application of this model rather than a separate model, and is",
      "reproduced in the validation vignette. The authors' central dosing",
      "conclusion is that the licensed regimen is likely inadequate and that a",
      "100 mg loading dose followed by 100 mg daily is preferable for patients",
      "with a lean body weight of 40-60 kg."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # Structural parameters: the "Estimate (%RSE)" column of Abdul-Aziz 2025
    # Table 2, cross-checked against the Abstract ("The typical volume of
    # distribution and clearance of caspofungin in this cohort were 8.13 L
    # and 0.55 L/h, respectively"). The paired "Bootstrap median (95% CI)"
    # column is reported for confidence in the point estimates and is quoted
    # in the comments below, but is NOT carried into any omega -- it is
    # parameter precision, not between-subject spread.
    # -----------------------------------------------------------------------

    lcl <- log(0.55)
    label("Clearance (L/h)")
    # Table 2: CL = 0.55 L/h (%RSE 11.9); bootstrap median 0.55
    # (95% CI 0.42-0.68). No covariate was retained on clearance -- serum
    # albumin was tested and rejected (Results, 'Population pharmacokinetic
    # model building').

    lvc <- log(8.13)
    label("Volume of distribution at FFM = 58.9 kg (L)")
    # Table 2: V = 8.13 L (%RSE 7.15); bootstrap median 8.27
    # (95% CI 7.07-10.8). This is the typical value at the reference
    # fat-free mass of 58.9 kg, i.e. Vpop in the covariate equation below.

    e_ffm_vc <- 1.24
    label("Power exponent on (FFM / 58.9) for the volume of distribution (unitless)")
    # Table 2 row 'LBW effect on V' = 1.24 (%RSE 23.3); bootstrap median 1.24
    # (95% CI 0.08-1.91). Estimated, not fixed: it carries an RSE and a
    # bootstrap interval, and it is above the theory-based allometric value
    # of 1 for a volume. The functional form is the paper's own displayed
    # equation, typeset in the Results between the covariate paragraph and
    # the Fig. 1 caption:
    #
    #     V = Vpop * (LBW / 58.9)^1.24
    #
    # This is a POWER model on median-centred fat-free mass, not a linear or
    # multiplicative one; the equation is unambiguous in the typeset source.
    # The 58.9 kg divisor is the cohort median lean body weight (Table 1),
    # matching the Methods statement that continuous covariates were centred
    # on their medians.

    # -----------------------------------------------------------------------
    # Between-subject variability. Table 2 reports the BSV rows as
    # percentages -- "CL (%) 25.5" and "V (%) 7.41" -- i.e. as coefficients
    # of variation, which is Monolix's standard summary for a parameter
    # declared log-normal. The Methods fix the distribution: BSV "was
    # described using an exponential model" theta_j = theta_p * exp(eta_j)
    # with eta_j ~ N(0, omega^2), and "individual estimates for
    # pharmacokinetic parameters were assumed to follow a log-normal
    # distribution". nlmixr2 wants the VARIANCE of eta, so each CV% is
    # back-transformed exactly as omega^2 = log(1 + CV^2).
    #
    # The alternative reading -- that the printed percentages are raw omega
    # SDs scaled by 100 rather than CVs -- gives omega_CL = 0.255 vs 0.251
    # and omega_V = 0.0741 vs 0.0740, a difference under 2% relative that no
    # validation gate can resolve. See the vignette Assumptions and
    # deviations. Both etas are poorly estimated in a cohort of 8 patients
    # (%RSE 34.3 on CL and 101 on V), so treat the stochastic layer as
    # indicative.
    # -----------------------------------------------------------------------

    etalcl ~ 0.062996
    # Table 2, BSV CL = 25.5% CV (%RSE 34.3); bootstrap median 21.3%
    # (95% CI 2.86-34.9). log(1 + 0.255^2) = 0.062996, i.e. omega = 0.251.

    etalvc ~ 0.0054758
    # Table 2, BSV V = 7.41% CV (%RSE 101 -- essentially unidentified);
    # bootstrap median 4.94% (95% CI 1.03-22.6).
    # log(1 + 0.0741^2) = 0.0054758, i.e. omega = 0.0740.

    # -----------------------------------------------------------------------
    # Residual error: ADDITIVE only. "Residual unexplained variability was
    # described by an additive error model" (Results, 'Population
    # pharmacokinetic model building'); additive, proportional and combined
    # forms were all tested (Methods, 'Structural model development'). There
    # is no proportional component to carry.
    # -----------------------------------------------------------------------

    addSd <- 1.20
    label("Additive residual error on Cc (mg/L)")
    # Table 2, Residual error, Additive = 1.20 mg/L (%RSE 10.0); bootstrap
    # median 1.17 (95% CI 0.85-1.47). Note this is large relative to the
    # 0.1 mg/L assay limit and to late-interval troughs, so simulated
    # observations near the trough can go negative; the vignette works with
    # the model's Cc prediction (IPRED) for its exposure checks rather than
    # with residual-perturbed observations.
  })

  model({
    # Reference fat-free mass: the cohort median lean body weight, exactly as
    # printed in the source's displayed covariate equation and in Table 1.
    ffm_ref <- 58.9 # kg

    # No covariate on clearance.
    cl <- exp(lcl + etalcl)

    # V = Vpop * (LBW / 58.9)^1.24, the paper's displayed equation. The
    # fat-free-mass term multiplies the typical value, and the exponential
    # between-subject term multiplies it in turn, per the Methods' separate
    # statements of the covariate model and the BSV model.
    vc <- exp(lvc + etalvc) * (FFM / ffm_ref)^e_ffm_vc

    kel <- cl / vc

    # One compartment, first-order elimination, dosed by intravenous
    # infusion directly into the central compartment (no absorption and no
    # bioavailability term: "Caspofungin was reconstituted and administered
    # as an intermittent intravenous infusion").
    d/dt(central) <- -kel * central

    # TOTAL plasma caspofungin, in mg/L.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
