Leegwater_2025_trimethoprim <- function() {
  description <- paste(
    "One-compartment population PK model for trimethoprim with first-order",
    "absorption in hospitalized adults treated with oral or intravenous",
    "cotrimoxazole (trimethoprim/sulfamethoxazole), built from routine",
    "therapeutic-drug-monitoring data in three Dutch university medical",
    "centers (Leegwater 2025). Renal function enters as a power function of",
    "BSA-normalized eGFR on apparent clearance, but only in patients NOT",
    "receiving continuous renal replacement therapy: CRRT replaces the eGFR",
    "term outright with a single multiplicative factor, and CRRT patients",
    "also carry their own between-subject variability on clearance. Oral",
    "bioavailability was estimated at ~100% and then fixed to 1, so the",
    "clearance and volume are apparent values that apply to both routes.",
    "Dose amounts are the TRIMETHOPRIM component of the combination product",
    "(cotrimoxazole 1,920 mg delivers 320 mg trimethoprim).",
    sep = " "
  )
  reference <- paste(
    "Leegwater E, Baidjoe L, Wilms EB, Visser LG, Touw DJT, de Winter BCM,",
    "de Boer MGJ, van Paassen J, van den Berg CHSB, van Prehn J,",
    "van Gelder T, Moes DJAR. Population Pharmacokinetics of",
    "Trimethoprim/Sulfamethoxazole: Dosage Optimization for Patients with",
    "Renal Insufficiency or Receiving Continuous Renal Replacement Therapy.",
    "Clin Pharmacol Ther. 2025;117(1):184-192. doi:10.1002/cpt.3421.",
    "Fixed effects, between-subject variability and residual error from",
    "Table 2 and its clearance footnote; the piecewise CRRT-vs-eGFR",
    "clearance structure, the CRRT-specific eta and the OMEGA variances",
    "from the trimethoprim NONMEM control stream reproduced in the",
    "Supplement ('Supplement NONMEM code', ADVAN2 TRANS2).",
    sep = " "
  )
  vignette <- "Leegwater_2025_cotrimoxazole"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Estimated glomerular filtration rate calculated with the CKD-EPI equation and reported BSA-normalized",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters apparent clearance as the power function (CRCL/68)^e_crcl_cl",
        "(Table 2 footnote a: 'in case of no CRRT = 4.21 x (EGFR/68)^0.317').",
        "The normalizing value is 68 mL/min/1.73 m^2 as printed in the",
        "footnote and in the control stream ('(EGFR/68)**THETA(5)'); note",
        "this differs slightly from the 70 mL/min/1.73 m^2 median and the",
        "70.8 mean reported for the cohort in Table 1 and the Results, and",
        "68 is the value the model actually uses.",
        "The eGFR term is switched OFF entirely for patients on CRRT -- the",
        "control stream evaluates TVCL = THETA(1)*(THETA(6)**rel) with no",
        "eGFR term in that branch -- so CRCL is ignored whenever",
        "RRT_CRRT_STATUS = 1. A CRCL value must nevertheless be supplied for",
        "every subject because the arithmetic switch in model() evaluates",
        "both branches; any positive placeholder works for CRRT subjects.",
        "The Discussion notes a known upward bias in serum creatinine (and",
        "therefore a downward bias in eGFR) of roughly 10-30% during",
        "trimethoprim treatment because trimethoprim inhibits tubular",
        "creatinine secretion; the model was fitted to CKD-EPI eGFR as",
        "measured, so this bias is inside the estimate rather than",
        "corrected for.",
        sep = " "
      ),
      source_name        = "EGFR"
    ),
    RRT_CRRT_STATUS = list(
      description        = "Continuous renal replacement therapy during cotrimoxazole treatment",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CRRT)",
      notes              = paste(
        "1 = concomitant CRRT, 0 = no CRRT. Called 'rel' in the NONMEM",
        "control stream and 'CRRT' in Table 2. Subject-level and time-fixed:",
        "patients treated with intermittent hemodialysis or ECMO were",
        "excluded from the study, so the flag never changes within a",
        "subject. 14 of the 52 trimethoprim subjects (27%) were on CRRT",
        "(Table 1).",
        "This covariate does two things at once, both taken from the control",
        "stream: it multiplies typical clearance by e_rrt_crrt_status_cl in",
        "place of (not on top of) the eGFR power term, and it selects which",
        "of the two estimated clearance etas applies to the subject.",
        sep = " "
      ),
      source_name        = "rel"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "trimethoprim", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "trimethoprim", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 52,
    n_studies      = 1,
    age_mean       = "60.3 years (SD 11.9)",
    weight_mean    = "80.6 kg (SD 18.5)",
    sex_female_pct = 36.6,
    race_ethnicity = "Not reported.",
    disease_state  = paste(
      "Hospitalized adults (>= 18 years) treated with therapeutic doses of",
      "oral or intravenous cotrimoxazole, sampled as part of routine",
      "therapeutic drug monitoring. Indications include Pneumocystis",
      "jirovecii pneumonia and other infections requiring high-dose",
      "cotrimoxazole. In the 168-patient full cohort, 20.8% were solid-organ",
      "transplant recipients, 13.7% had a malignancy, 11.3% were living with",
      "HIV, 7.7% had a stem-cell transplant and 36.9% received",
      "corticosteroids.",
      sep = " "
    ),
    renal_function = paste(
      "eGFR (CKD-EPI) mean 49.2 mL/min/1.73 m^2 (SD 34.7) in the 52",
      "trimethoprim subjects, versus 70.8 (SD 33.2) in the full 168-patient",
      "cohort; serum creatinine mean 170.1 umol/L (SD 101.5). 14 of 52 (27%)",
      "received concomitant CRRT. Patients on intermittent hemodialysis or",
      "ECMO were excluded.",
      sep = " "
    ),
    dose_range     = paste(
      "Routine care rather than protocol-assigned. Daily cotrimoxazole",
      "starting doses in the trimethoprim subset: <= 960 mg 11.5%,",
      "1,920-2,400 mg 36.6%, 2,880 mg 11.5%, 3,840-4,800 mg 13.5%,",
      "5,760 mg 27.0%. 42.3% started oral, 57.7% intravenous.",
      sep = " "
    ),
    regions        = "The Netherlands (Leiden University Medical Center and University Medical Center Groningen; trimethoprim was assayed in only two of the three participating centers).",
    n_observations = "137 trimethoprim plasma concentrations from 52 patients, peaks and troughs, ranging 0.2-15.6 mg/L.",
    notes          = paste(
      "Retrospective multicenter observational cohort, January 2016 to",
      "December 2021 (Methods, 'Study design' and 'Participants and data",
      "collection'); demographics from Table 1, right-hand column. The",
      "sulfamethoxazole / N-acetyl sulfamethoxazole model of the same paper",
      "was fitted to the larger 168-patient cohort -- see",
      "modellib('Leegwater_2025_sulfamethoxazole').",
      sep = " "
    )
  )

  ini({
    # Structural parameters - Leegwater 2025 Table 2, cross-checked against
    # the trimethoprim $THETA block of the supplement control stream (all
    # values there are the final estimates entered as FIX for simulation).
    lka <- log(0.337); label("First-order absorption rate constant (ka, 1/h)")                                       # Table 2: 0.337 1/h, RSE 54%, bootstrap 95% CI 0.095-0.877; control stream $THETA(3) "0.337 FIX". The Discussion flags this as the one imprecise parameter, attributed to routinely-collected dosing times.
    lcl <- log(4.21); label("Apparent clearance at eGFR 68 mL/min/1.73 m^2 without CRRT (CL/F, L/h)")                # Table 2: 4.21 L/h, RSE 13%, bootstrap 95% CI 3.21-5.42; control stream $THETA(1) "4.21 FIX"
    lvc <- log(134); label("Apparent central volume of distribution (Vc/F, L)")                                      # Table 2: 134 L, RSE 10%, bootstrap 95% CI 105-161; control stream $THETA(2) "134 FIX"
    lfdepot <- fixed(log(1)); label("Oral bioavailability (F, unitless)")                                            # Table 2: "Biological availability 1 Fixed". Results: estimated first, found to be ~100%, then fixed to 100% with similar model performance. The control stream carries no F1 term, i.e. F = 1.

    # Covariate effects on apparent clearance. These are mutually exclusive
    # branches, not multiplicative layers: control stream $PK evaluates
    # IF(rel==0) TVCL = THETA(1)*((EGFR/68)**THETA(5)) and
    # IF(rel==1) TVCL = THETA(1)*(THETA(6)**rel), so a CRRT patient's
    # clearance carries no eGFR term at all.
    e_crcl_cl <- 0.317; label("Power exponent on (CRCL/68) for CL in patients not receiving CRRT (unitless)")        # Table 2 "eGFR on CL": 0.317, RSE 36%, bootstrap 95% CI 0.049-0.558; control stream $THETA(5) "0.317 FIX"; footnote a
    e_rrt_crrt_status_cl <- 1.12; label("Multiplicative factor on CL for patients receiving CRRT (unitless)")        # Table 2 "CRRT on CL": 1.12, RSE 16%, bootstrap 95% CI 0.77-1.57; control stream $THETA(6) "1.12 FIX"; footnote a gives "in case of CRRT: 4.21 x 1.12"

    # Between-subject variability. The supplement's $OMEGA block carries
    # variances on the log scale (CL and V both enter as TV*EXP(ETA)), and
    # Table 2 reports each as a percentage that is exactly sqrt(omega^2):
    # sqrt(0.161) = 40.1%, sqrt(0.108) = 32.9%, sqrt(0.102) = 31.9%. The
    # percentages are therefore omega itself, and ini() takes omega^2
    # directly from the control stream.
    #
    # Clearance carries TWO etas because the Results state "As the variation
    # in clearance was small for patients on CRRT, we estimated the IIV for
    # the clearance for patients treated with CRRT and patients without
    # CRRT"; the control stream implements this as ETA(1) in the rel==0
    # branch and ETA(3) in the rel==1 branch. Both strata therefore carry an
    # explicit suffix (references/parameter-names.md "Stratum-suffixed
    # parameters"): a bare etalcl silently meaning "the non-CRRT value" is
    # exactly the ambiguity the suffix removes.
    etalcl_nocrrt ~ 0.161                                                                                           # control stream $OMEGA(1) '0.161 FIX'; Table 2 'IIV CL (%)' 40.1, RSE 12%, bootstrap 95% CI 26.9-48.5
    etalcl_crrt   ~ 0.108                                                                                           # control stream $OMEGA(3) '0.108 FIX'; Table 2 'IIV CL patients on CRRT (%)' 32.9, RSE 32%, bootstrap 95% CI 7.54-48.7
    etalvc        ~ 0.102                                                                                           # control stream $OMEGA(2) '0.102 FIX'; Table 2 'IIV Vd (%)' 31.9, RSE 32%, bootstrap 95% CI 9.77-49.2

    # Residual error. Control stream $ERROR: W = SQRT(THETA(4)**2*IPRED**2),
    # Y = IPRED + W*EPS(1) with $SIGMA 1 FIX, i.e. a pure proportional error
    # whose SD is THETA(4) = 0.169.
    propSd <- 0.169; label("Proportional residual error SD (fraction)")                                             # Table 2 "Proportional error / Trimethoprim": 0.169, RSE 10%, bootstrap 95% CI 0.131-0.198; control stream $THETA(4) "0.169 FIX"
  })

  model({
    # Normalizing eGFR for the power term (Table 2 footnote a and the control
    # stream both use 68, not the cohort median of 70).
    ref_crcl <- 68

    ka <- exp(lka)

    # Apparent clearance. The two NONMEM IF branches are written here as an
    # arithmetic switch on the binary covariate so the expression stays
    # differentiable and avoids branch-dependent parameter binding:
    #   RRT_CRRT_STATUS = 0 -> exp(lcl + etalcl_nocrrt) * (CRCL/68)^e_crcl_cl
    #   RRT_CRRT_STATUS = 1 -> exp(lcl + etalcl_crrt)   * e_rrt_crrt_status_cl
    # The eGFR power term is multiplied by (1 - RRT_CRRT_STATUS) INSIDE the
    # sum rather than raised to it, so a CRRT subject never evaluates a
    # fractional power of a possibly-zero eGFR.
    cl <-
      exp(lcl + etalcl_nocrrt * (1 - RRT_CRRT_STATUS) + etalcl_crrt * RRT_CRRT_STATUS) *
      ((1 - RRT_CRRT_STATUS) * (CRCL / ref_crcl)^e_crcl_cl +
         RRT_CRRT_STATUS * e_rrt_crrt_status_cl)

    vc  <- exp(lvc + etalvc)
    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    f(depot) <- exp(lfdepot)

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
