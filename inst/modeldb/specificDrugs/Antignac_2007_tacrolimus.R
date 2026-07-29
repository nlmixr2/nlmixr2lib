Antignac_2007_tacrolimus <- function() {
  description <- paste(
    "One-compartment population PK model for oral and intravenous",
    "tacrolimus in adult kidney transplant recipients (Antignac 2007).",
    "First-order absorption (ka fixed at 4.5 1/h from Jusko 1995),",
    "linear elimination, simultaneous fit of IV and oral data to",
    "estimate bioavailability. Clearance increases sigmoidally with",
    "days postoperation from a baseline CLmin (at POD = 0) to 2 *",
    "CLmin asymptotically, with half-recovery at TCL50 = 3.81 days",
    "and Hill exponent 2.54; clearance is multiplied by",
    "(1 + theta_PRD) when concomitant prednisone dose exceeds 25",
    "mg/day. No covariates retained on V or F."
  )
  reference <- paste(
    "Antignac M, Barrou B, Farinotti R, Lechat P, Urien S.",
    "Population pharmacokinetics and bioavailability of tacrolimus",
    "in kidney transplant patients.",
    "Br J Clin Pharmacol. 2007;64(6):750-757.",
    "doi:10.1111/j.1365-2125.2007.02895.x"
  )
  vignette <- "Antignac_2007_tacrolimus"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    POD = list(
      description        = "Days elapsed since kidney transplantation (post-operative day)",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying within subject. Day 1 corresponds to the day",
        "of transplantation per Antignac 2007 Methods 'Patients and",
        "data collection' (item 3). Enters CL through a Hill-type",
        "recovery curve: CL = CLmin * (1 + POD^gamma_cl /",
        "(POD^gamma_cl + tcl50^gamma_cl)). At POD -> infinity the",
        "multiplier saturates to 2 (so CLmax = 2 * CLmin); at POD =",
        "TCL50 = 3.81 days the multiplier is 1.5. Cohort range 1 -",
        "158 days (Table 1, median 33 days, mean 38 days).",
        "Use POD >= 0 to avoid 0^gamma_cl evaluation; the paper's",
        "POD = 1 anchor sidesteps the issue in practice."
      ),
      source_name        = "POD"
    ),
    PRED_DOSE = list(
      description        = paste(
        "Concomitant oral prednisone daily dose. Used in a",
        "threshold-form (> 25 mg/day) binary indicator on CL."
      ),
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Antignac 2007 Methods records prednisone as",
        "'PRD (mg)' in the daily dose register (item 4); the daily",
        "reading is the standard transplant-immunosuppressant",
        "tapering schedule. Table 1 mean 22 mg, SD 4 mg, range 0 -",
        "30 mg. Final-model effect (Antignac 2007 Eq. p. 754):",
        "E_PRD = 1 + theta_PRD when PRED_DOSE > 25 mg/day, else 1.",
        "theta_PRD = 0.575 (Table 4 Original dataset Mean), giving",
        "a ~1.6-fold rise in CL above the threshold. The binary",
        "indicator is derived inside model() from the continuous",
        "PRED_DOSE column, following the canonical-register pattern",
        "established in TerHeine_2018_everolimus.R."
      ),
      source_name        = "PRD"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 83L,
    n_studies      = 1L,
    n_observations = 1589L,
    age_range      = "16 - 67 years",
    age_median     = "46 years",
    weight_range   = "29 - 100 kg",
    weight_median  = "69 kg",
    sex_female_pct = 34.9,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Adult kidney transplant recipients on a triple",
      "immunosuppressive regimen (tacrolimus + mycophenolate mofetil",
      "+ corticosteroids). Surgery between January 2003 and",
      "September 2004 at a single French centre (Pitie-Salpetriere",
      "Hospital, Paris); 14 grafts from living donors, 69 from",
      "deceased (brain-death) donors."
    ),
    dose_range     = paste(
      "Initial tacrolimus dose ~0.015 mg/kg twice daily.",
      "Subsequent oral doses titrated by TDM (target trough 5 - 10",
      "ng/mL during the first 3 months posttransplant). Table 2:",
      "per-record dose mean 2792 mg, range 500 - 7000 mg (note:",
      "Table 2 'Dose (mg)' is the cumulative total recorded per",
      "subject, not a single dose). Initial dose mean 0.02 mg/kg",
      "(0.013 mg/kg IV, 0.02 mg/kg oral). 25 patients received an",
      "initial IV bolus and were switched to oral 24 - 48 h after",
      "treatment initiation; the remaining 58 patients received",
      "oral only."
    ),
    regions        = "France",
    pod_range      = "1 - 158 days (Table 1)",
    pred_dose_range = "0 - 30 mg/day (Table 1; mean 22, SD 4)",
    notes          = paste(
      "Routine therapeutic drug monitoring (TDM) data, mostly trough",
      "samples (predose). Mean 19.1 +/- 4.7 samples per subject",
      "(range 1 - 33). Whole-blood tacrolimus measured by",
      "microparticle enzyme immunoassay (MEIA) on the IMx platform",
      "(LLOQ 1.5 ng/mL, linear 1.5 - 30 ng/mL; cross-reactive with",
      "metabolites M-II, M-III, MV). NONMEM V level 1.1 double",
      "precision, FOCE, subroutines ADVAN2 TRANS2. Bootstrap",
      "validation 910 of 1000 successful runs (Table 4)."
    )
  )

  ini({
    # Structural PK parameters -- Antignac 2007 Table 4 Original dataset
    # Mean column. CLmin is the baseline (POD = 0) clearance prior to
    # the Hill-type recovery; V is the central volume; F is oral
    # bioavailability; ka was fixed at 4.5 1/h to the literature value
    # of Jusko 1995 (ref 10) because only trough samples were available
    # in the dataset and ka could not be identified from the data.
    lcl <- log(1.81)
    label("Baseline CL at POD = 0 (CLmin; L/h)")
    # Antignac 2007 Table 4 TV(CLmin) Mean = 1.81 L/h
    lvc <- log(98.4)
    label("Central volume V (L)")
    # Antignac 2007 Table 4 V Mean = 98.4 L
    lfdepot <- log(0.137)
    label("Oral bioavailability F (fraction)")
    # Antignac 2007 Table 4 F Mean = 13.7% = 0.137
    lka <- fixed(log(4.5))
    label("Absorption rate constant ka (1/h; FIXED from Jusko 1995)")
    # Antignac 2007 Table 4 'ka, fixed (1/h)' = 4.5; held fixed per
    # Results 'Population pharmacokinetics' (ka could not be estimated
    # from trough-only data; refitting with ka in {1.5, 2, 3} produced
    # no significant OFV change)

    # POD covariate effect on CL (Antignac 2007 Eq. p. 754,
    # CL = CLmin * [1 + POD^gamma_CL / (POD^gamma_CL + TCL50^gamma_CL)]
    # * E_PRD). The half-recovery time and Hill exponent are reported
    # on the linear scale with no IIV (Table 4 Original dataset Mean).
    tcl50 <- 3.81
    label("Days postoperation at which CL has recovered to 1.5 * CLmin (days)")
    # Antignac 2007 Table 4 TCL50 Mean = 3.81 days
    gamma_cl <- 2.54
    label("Hill exponent on the POD-driven CL recovery curve (unitless)")
    # Antignac 2007 Table 4 gCl Mean = 2.54

    # Prednisone-dose covariate effect on CL (Antignac 2007 Eq. p. 754,
    # E_PRD = 1 + theta_PRD if PRD > 25 mg, else 1). The threshold
    # indicator is derived inside model() from the continuous PRED_DOSE
    # column.
    e_pred_dose_high_cl <- 0.575
    label("Fractional increase in CL when PRED_DOSE > 25 mg/day (unitless)")
    # Antignac 2007 Table 4 qPRD Mean = 0.575

    # Inter-individual variability -- Antignac 2007 Table 4 reports IIV
    # as %CV (rows 'wClmin (%)', 'wV (%)', 'wF (%)') under the
    # exponential IIV model stated in Methods 'Population
    # pharmacokinetic modelling' (CL_j = TV(CL) * exp(eta_CL)). Convert
    # to log-scale variance via omega^2 = log(CV^2 + 1):
    #   CLmin  CV 31% -> log(0.31^2 + 1) = 0.09178
    #   V      CV 79% -> log(0.79^2 + 1) = 0.48490
    #   F      CV 32% -> log(0.32^2 + 1) = 0.09753
    etalcl ~ log(0.31^2 + 1)
    # Antignac 2007 Table 4 wClmin (%) = 31
    etalvc ~ log(0.79^2 + 1)
    # Antignac 2007 Table 4 wV (%) = 79
    etalfdepot ~ log(0.32^2 + 1)
    # Antignac 2007 Table 4 wF (%) = 32

    # Residual error -- Antignac 2007 Table 4 reports the proportional
    # + additive combined model. The 'Residual additive variance, s
    # (ng/mL)' and 'Residual proportional variance, s (%)' labels in
    # Table 4 use the symbol s (sigma) with units of ng/mL and percent,
    # consistent with standard deviations (not variances) under the
    # 'proportional plus additive' error model stated in Results
    # 'Population pharmacokinetics'.
    propSd <- 0.186
    label("Proportional residual error (fraction)")
    # Antignac 2007 Table 4 Residual proportional variance, s (%) = 18.6
    addSd <- 0.96
    label("Additive residual error (ng/mL)")
    # Antignac 2007 Table 4 Residual additive variance, s (ng/mL) = 0.96
  })

  model({
    # Individual structural PK parameters with exponential (log-normal)
    # IIV. ka is held fixed at the population value (no IIV in the
    # source).
    cl_baseline <- exp(lcl + etalcl)
    vc          <- exp(lvc + etalvc)
    fdepot      <- exp(lfdepot + etalfdepot)
    ka          <- exp(lka)

    # POD-driven Hill recovery of CL (Antignac 2007 Eq. p. 754, first
    # bracketed term). At POD = 0 the multiplier is 1 (CL = CLmin); as
    # POD -> infinity it saturates to 2 (CLmax = 2 * CLmin); at POD =
    # TCL50 it equals 1.5.
    pod_effect <- 1 + POD^gamma_cl / (POD^gamma_cl + tcl50^gamma_cl)

    # High-prednisone-dose threshold effect on CL (Antignac 2007 Eq. p.
    # 754, second term). Derived as a binary indicator inside model()
    # so the underlying continuous PRED_DOSE column remains available;
    # same pattern as TerHeine_2018_everolimus.R.
    pred_high   <- (PRED_DOSE > 25)
    pred_effect <- 1 + e_pred_dose_high_cl * pred_high

    cl <- cl_baseline * pod_effect * pred_effect

    # One-compartment disposition with first-order absorption.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - (cl / vc) * central

    # Oral bioavailability (applied to the depot input only). IV doses
    # entering the central compartment directly are unaffected.
    f(depot) <- fdepot

    # Observation: tacrolimus whole-blood concentration in ng/mL.
    # Dose in mg, vc in L => central / vc has units of mg/L; multiply
    # by 1000 to convert to ng/mL (the source paper's reporting unit).
    Cc <- central / vc * 1000
    Cc ~ prop(propSd) + add(addSd)
  })
}
