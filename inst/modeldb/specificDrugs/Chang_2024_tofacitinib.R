Chang_2024_tofacitinib <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and an",
    "absorption lag time for oral tofacitinib (tablet or 1 mg/mL oral",
    "solution, weight-based doses twice daily) in 246 children and",
    "adolescents aged 2-17 years with juvenile idiopathic arthritis, pooled",
    "from a phase 1 (NCT01513902), a phase 3 (NCT02592434) and a long-term",
    "extension (NCT01500551) study (Chang 2024). The model is parameterized in",
    "apparent oral clearance (CL/F) and apparent volume of distribution (V/F),",
    "both allometric in time-varying body weight with estimated exponents",
    "(0.31 and 0.54) about the 46.3 kg median. The oral solution increases",
    "the absorption rate constant by a fractional 1.64 over the tablet",
    "(ka x 2.64). V/F has no eta of its own: its deviation is 0.27 times the",
    "CL/F eta. Log-transformed observations with an additive log-scale",
    "residual whose SD switches at 1.08 h after dose and carries its own",
    "inter-individual variability, correlated with the CL/F and ka etas."
  )
  reference <- paste(
    "Chang C, Vong C, Wang X, Hazra A, Diehl A, Nicholas T, Mukherjee A.",
    "Tofacitinib pharmacokinetics in children and",
    "adolescents with juvenile idiopathic arthritis. CPT Pharmacometrics Syst",
    "Pharmacol. 2024;13(4):599-611. doi:10.1002/psp4.13104. PMID: 38298058.",
    "PMCID: PMC11015083."
  )
  vignette <- "Chang_2024_tofacitinib"
  # etaruv scales the residual SD (control stream W = PROP*EXP(ETA(3))); it
  # has no structural typical-value parameter.
  paper_specific_etas <- c("etaruv")
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight (time-varying)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying: weight was recorded at screening and repeatedly during",
        "the phase 3 and long-term extension studies, carried forward when",
        "missing (Methods, Study design). Enters CL/F and V/F as power",
        "functions normalized to 46.3 kg, the cohort median (Methods; supplement",
        "control stream `(WT/46.3)**COEFCL`). Cohort range 11.1-121.8 kg",
        "(Table 2); the dosing simulations extend down to 10 kg."
      ),
      source_name = "WT"
    ),
    FORM_SOLUTION = list(
      description = "Oral solution (1 mg/mL) formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "tablet (FORM_SOLUTION = 0)",
      notes = paste(
        "Source column FORM: 2 = tablet (the reference, 'Most common' in the",
        "control stream), 4 = oral solution, so FORM_SOLUTION = (FORM == 4).",
        "The solution multiplies ka by (1 + 1.64) = 2.64; it has no effect on",
        "bioavailability, CL/F or V/F. In the studies every dose below 5 mg was",
        "given as the oral solution and 5 mg was given as a tablet to patients",
        ">= 40 kg able to swallow it (Methods). Initially 100 patients received",
        "solution and 146 tablets; 11 switched formulation during the studies",
        "(Table 2 footnote). Set per dosing record when the formulation changes."
      ),
      source_name = "FORM"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age (time-varying and at baseline)",
      units = "years",
      type = "continuous",
      notes = "Screened on CL/F and V/F in the stepwise covariate search but not retained (Results, PK modeling; Table S2)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL/F and V/F but not retained (Results; Table S2)."
    ),
    CRCL_BASE = list(
      description = "Baseline creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Screened but not retained (Results, PK modeling)."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "tofacitinib",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "tofacitinib",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 246,
    n_studies = 3,
    n_observations = 1392,
    age_mean = "11.6 years",
    age_range = "2-17 years",
    weight_median = "46.3 kg",
    weight_range = "11.1-121.8 kg",
    sex_female_pct = 74.0,
    race_ethnicity = c(White = 87.8, Black = 2.0, Other = 10.2, Asian = 0),
    disease_state = paste(
      "juvenile idiopathic arthritis: extended oligoarthritis 12.6%, RF+",
      "polyarthritis 15.4%, RF- polyarthritis 48.8%, systemic JIA without",
      "active systemic features 4.9%, juvenile psoriatic arthritis 8.9%,",
      "enthesitis-related arthritis 9.3%"
    ),
    dose_range = "1-5 mg twice daily by body-weight band (Table S1); oral solution 1 mg/mL below 5 mg, tablet for 5 mg in patients >= 40 kg",
    co_medication = "NSAIDs 16.3%, glucocorticoids 5.3%, methotrexate 15.4% (Table 2)",
    regions = "multinational phase 1 (NCT01513902), phase 3 (NCT02592434) and long-term extension (NCT01500551, PK cutoff 21 May 2019) studies",
    notes = paste(
      "Demographics from Table 2 of Chang 2024. 65 of 1392 records were below",
      "the 0.100 ng/mL LLOQ and were ignored. Sampling was rich in the phase 1",
      "study (predose and 0.5, 1, 4, 8 h on day 5) and sparse thereafter",
      "(Table 1). 16.3% of patients were overweight and 11.4% obese by CDC",
      "criteria (Results)."
    )
  )

  ini({
    # Typical values are for the reference patient: 46.3 kg taking the tablet
    # (Results, PK modeling). Table 3 rows follow the supplement control
    # stream's THETA order (CL, V, KA, COEFCL, COEFV, SCAL, PROP, ALAG1, PROP2,
    # KAFORM).
    lka <- log(2.78); label("Absorption rate constant for the tablet (1/h)")              # Table 3, 'k a (h-1)' = 2.78 (RSE 12.04%), SIR 95% CI 2.29-3.44
    lcl <- log(26.05); label("Apparent oral clearance CL/F at 46.3 kg (L/h)")             # Table 3, 'CL/F (L/h)' = 26.05 (RSE 2.30%), SIR 95% CI 25.20-27.13
    lvc <- log(89.22); label("Apparent volume of distribution V/F at 46.3 kg (L)")        # Table 3, 'V/F (L)' = 89.22 (RSE 2.55%), SIR 95% CI 86.12-92.91
    ltlag <- log(0.19); label("Absorption lag time (h)")                                  # Table 3, 'Lag time (h)' = 0.19 (RSE 5.02%), SIR 95% CI 0.17-0.20

    e_wt_cl <- 0.31; label("Allometric exponent on (WT/46.3) for CL/F (unitless)")        # Table 3, 'Exponent on CL/F' = 0.31 (RSE 13.40%), SIR 95% CI 0.25-0.37
    e_wt_vc <- 0.54; label("Allometric exponent on (WT/46.3) for V/F (unitless)")         # Table 3, 'Exponent on V/F' = 0.54 (RSE 9.57%), SIR 95% CI 0.48-0.60

    # Control stream: KAFORM = 1 + THETA(11) for FORM 4 (solution). Table 3
    # prints THETA(11) itself; see the vignette for why 1.64 is the fractional
    # increase and not the ka ratio the Results prose states.
    e_form_solution_ka <- 1.64; label("Fractional increase in ka for the oral solution vs tablet (unitless)")  # Table 3, 'Formulation effect' = 1.64 (RSE 38.37%), SIR 95% CI 0.95-2.68

    # V = TVV * EXP(ETA(1) * SCAL): V/F has no eta of its own. Table 3 labels
    # this row 'Covariance between CL/F and V/F'; it is THETA(6), SCAL.
    vc_eta_scale <- 0.27; label("Scaling factor relating the V/F deviation to etalcl (unitless)")  # Table 3, 'Covariance between CL/F and V/F' = 0.27 (RSE 38.91%), SIR 95% CI 0.14-0.41

    # OMEGA BLOCK(3) on ETA(1) CL/F, ETA(2) ka, ETA(3) residual SD. The two
    # structural diagonals use the Results-text CVs (Table 3 prints them to two
    # decimals) under CV = sqrt(omega); off-diagonals and the residual-SD
    # variance are Table 3 as printed.
    etalcl + etalka + etaruv ~ c(
      0.055225,             # Results 'IIV ... for CL/F ... 23.5%': 0.235^2; Table 3 'IIVCL/F' = 0.06
      -0.10, 1.461681,      # Table 3 'IIVCL/F-ka' = -0.10; Results 'IIV ... ka ... 120.9%': 1.209^2, Table 3 'IIVka' = 1.46
      0.07, -0.57, 0.25     # Table 3 'IIVCL/F-residual error' = 0.07, 'IIVka-residual error' = -0.57, 'IIVresidual error' = 0.25
    )

    # Y = LOG(IPRED) + W*EPS(1) with SIGMA 1 FIX, W = PROP or PROP2 times
    # EXP(ETA(3)): the thetas are log-scale residual SDs.
    expSd_early <- 0.6796; label("Log-scale residual SD, time after dose <= 1.08 h")     # Table 3, 'Residual error before TAD = 1.08 h (%)' = 67.96 (RSE 6.17%); control stream PROP2
    expSd_late <- 0.4098; label("Log-scale residual SD, time after dose > 1.08 h")        # Table 3, 'Residual error (%)' = 40.98 (RSE 5.95%); control stream PROP
  })

  model({
    # Individual parameters (supplement control stream $PK)
    ka <- exp(lka + etalka) * (1 + e_form_solution_ka * FORM_SOLUTION)
    cl <- exp(lcl + etalcl) * (WT / 46.3)^e_wt_cl
    vc <- exp(lvc + vc_eta_scale * etalcl) * (WT / 46.3)^e_wt_vc
    tlag <- exp(ltlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    # Dose in mg and vc in L give mg/L; x1000 for ng/mL (control stream: AMT
    # in ug, S2 = V).
    Cc <- 1000 * central / vc

    # Residual SD switches on time after dose at 1.08 h and is scaled by the
    # per-subject etaruv. tad() is evaluated once on its own line.
    early_flag <- tad() <= 1.08
    expSdTad <- expSd_early * early_flag + expSd_late * (1 - early_flag)
    expSdCc <- expSdTad * exp(etaruv)
    Cc ~ lnorm(expSdCc)
  })
}
