AlQurain_2022_tramadol <- function() {
  description <- "Joint parent-and-metabolite population PK model for extended-release oral tramadol and its active CYP2D6-derived metabolite O-desmethyltramadol (ODT, M1) in older hospital inpatients (Al-Qurain 2022). Two-compartment tramadol disposition with first-order absorption and linear elimination; CL/F is the total apparent tramadol clearance, and ODT is formed at a first-order rate (Kt) from the tramadol central amount into a one-compartment apparent ODT pool that shares the tramadol central volume and is cleared linearly (CLm/F). Because the fraction of tramadol metabolised to ODT is not identifiable, the ODT pool is scaled by it and the formation term does not deplete tramadol. The Identification of Seniors At Risk (ISAR) frailty score raises the tramadol inter-compartmental clearance, and Cockcroft-Gault creatinine clearance raises both the tramadol apparent clearance and the peripheral volume, each through an uncentred exponential term."
  reference <- "Al-Qurain AA, Upton RN, Tadros R, Roberts MS, Wiese MD. Population Pharmacokinetic Model for Tramadol and O-desmethyltramadol in Older Patients. Eur J Drug Metab Pharmacokinet. 2022;47(3):387-402. doi:10.1007/s13318-022-00756-x"
  vignette <- "AlQurain_2022_tramadol"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(
      analyte = "tramadol",
      units = "mg",
      specimen = "administration site",
      verified = TRUE,
      notes = "Oral extended-release tramadol dose site (Figure 1 'Ka'); absorbed first-order into central. Bioavailability is set to 1 (Methods 2.4.2), so every volume and clearance is apparent (/F)."
    ),
    central = list(
      analyte = "tramadol",
      units = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes = "Tramadol central compartment (Figure 1 'V1'). Loses drug only by linear elimination (CL/F, the total apparent tramadol clearance, which includes the unidentified CYP2D6 route to ODT) and by inter-compartmental exchange with peripheral1 (Q). ODT formation (Kt * central) feeds central_m1 without a matching loss term here; see central_m1 notes. Concentrations were measured in finger-prick blood by volumetric absorptive microsampling and converted to plasma-equivalent values (Methods 2.3)."
    ),
    peripheral1 = list(
      analyte = "tramadol",
      units = "mg",
      specimen = "tissue",
      verified = TRUE,
      notes = "Tramadol peripheral compartment (Figure 1 'V2')."
    ),
    central_m1 = list(
      analyte = "O-desmethyltramadol (M1)",
      units = "mg (tramadol-equivalent; see notes)",
      specimen = "plasma",
      verified = TRUE,
      notes = paste(
        "ODT compartment (Figure 1, metabolite box). Fed by Kt * central and",
        "eliminated by CLm/F. The fraction of tramadol metabolised to ODT (Fm) is",
        "not identifiable from parent-only oral dosing; Methods 2.4.2 states that",
        "only Vm/(F km) and CLm/(F km) can be determined. The state therefore holds",
        "an APPARENT amount scaled by 1/Fm, and the formation term does not",
        "deplete the tramadol central compartment (tramadol is cleared by CL/F",
        "alone). Kt * V1 (18.4 L/h) exceeding CL/F (6-8 L/h) is a consequence of",
        "that scaling, not a mass-balance violation. This no-depletion reading",
        "reproduces the paper's steady-state VPC (Figure 6) and the CrCL dependence",
        "of Figure 11. No molecular-weight correction is applied on the",
        "transfer (tramadol 263.4 g/mol, ODT 249.3 g/mol), and no ODT volume is",
        "reported anywhere in the paper or supplement: the Monolix parent-metabolite",
        "parameterisation the authors used (Bertrand et al.) gives the metabolite",
        "the parent central volume because the metabolite volume is not",
        "identifiable from parent-only oral dosing (Methods 2.4.2). The ODT state",
        "therefore shares vc, and CLm/F is an apparent clearance that absorbs the",
        "unidentifiable fraction metabolised. See the vignette Assumptions and",
        "deviations.",
        sep = " "
      )
    )
  )

  covariateData <- list(
    SCORE_ISAR = list(
      description = "Identification of Seniors At Risk (ISAR) frailty screening score (integer 0-6; higher = frailer)",
      units = "(SCORE_ISAR units, 0-6 score)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Six yes/no screening questions (McCusker 1999); the paper classifies",
        "ISAR 0-2 as fit and ISAR >= 3 as frail. Cohort mean 3 (SD 1.2, range",
        "1-5; Table 1). Enters the tramadol inter-compartmental clearance as the",
        "uncentred Monolix-default exponential Q = Q_pop * exp(0.255 * ISAR), so",
        "the Table 3 Q of 42.6 L/h is the value at ISAR = 0. The paper prints no",
        "covariate equation; the uncentred form is supported by the Figure 10",
        "simulated profiles at ISAR 0 ('fit'), 3 and 6.",
        sep = " "
      ),
      source_name = "ISAR"
    ),
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation (raw, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Raw Cockcroft-Gault creatinine clearance in mL/min, NOT BSA-normalized",
        "(Methods 2.2). Cohort mean 65.1 (SD 69.1, range 16-351; Table 1).",
        "Enters the tramadol apparent clearance and peripheral volume as the",
        "uncentred Monolix-default exponentials CL/F = CL_pop * exp(0.00498 * CrCL)",
        "and V2/F = V2_pop * exp(0.0119 * CrCL), so the Table 3 typical values",
        "are the values at CrCL = 0. At the cohort mean (65.1 mL/min) CL/F is",
        "8.3 L/h and V2/F is 820 L.",
        sep = " "
      ),
      source_name = "CrCL"
    )
  )

  # Screened in the covariate model (Supplementary Table S5) but not retained
  # after backward deletion (Supplementary Table S7), so no coefficient exists.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on ka, CL/F, V1/F, Q, V2/F, Kt and CLm (Supplementary Table S5); not retained."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened as 'Gender' on every parameter (Supplementary Table S5); not retained."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Weight effects on CL/F, Q and V2/F entered the full model but were removed at backward deletion (Supplementary Table S7)."
    ),
    SCORE_CCI = list(
      description = "Charlson Comorbidity Index total score",
      units = "(SCORE_CCI units, weighted comorbidity count)",
      type = "continuous",
      notes = "The CCI effect on Q entered the full model but was removed at backward deletion (Supplementary Table S7)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 21L,
    n_studies = 1L,
    age_range = "69-93 years (mean 83.2, SD 7.5)",
    weight_range = "33-112 kg (mean 68.3, SD 20.6)",
    sex_female_pct = 57,
    disease_state = "Older hospital inpatients (aged 65 years or more) prescribed extended-release oral tramadol for pain; admission reasons included infection, falls, musculoskeletal pain, fractures, chest, abdominal and cancer pain. All had normal liver function tests.",
    dose_range = "Extended-release oral tramadol 25 mg (n = 3), 50 mg (n = 6), 75 mg (n = 1) or 100 mg (n = 11)",
    regions = "Australia (single centre: Royal Adelaide Hospital)",
    renal_function = "Cockcroft-Gault creatinine clearance 16-351 mL/min (mean 65.1, SD 69.1)",
    frailty = "ISAR score 1-5 (mean 3, SD 1.2)",
    co_medication = "CYP3A4 inhibitors 8 (38%); CYP2D6 inhibitors 8 (38%); 3-10 prescribed medications (mean 7.2)",
    notes = paste(
      "Table 1 baseline characteristics. BMI 15.5-39.4 kg/m2 (mean 25.8);",
      "Charlson Comorbidity Index 4-8 (mean 5.9). Up to five finger-prick",
      "volumetric absorptive microsamples per patient within 12 h of a dose (99",
      "samples, all above the 10 ng/mL limit of quantitation). Observed",
      "steady-state trough ODT/tramadol ratio 0.26 (SD 0.2, range 0.08-0.82;",
      "n = 17). Estimation: Monolix 2020, SAEM.",
      sep = " "
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # Units. Table 3 prints the volumes in 'l' and the clearances in 'l/h',
    # but the magnitudes (V1/F 0.373, CL/F 0.00604) are about 1000-fold below
    # any tramadol value in litres. Monolix applies no unit scaling, so with
    # dose in mg and concentrations in ng/mL the volume unit of the fit is
    # mg / (ng/mL) = 1000 L and the clearance unit is 1000 L/h. Every volume
    # and clearance below is the printed value x 1000. The Figure 10 typical
    # 100 mg profile (Cmax about 230 ng/mL) confirms the scaling.
    # -----------------------------------------------------------------------
    lka <- log(2.96); label("First-order absorption rate constant ka (1/h)") # Table 3 final model 'K a, /h' = 2.96 (RSE 63.7)
    lvc <- log(373); label("Apparent central volume of tramadol V1/F (L)") # Table 3 final model 'V 1 / F , l' = 0.373 printed; x 1000 = 373 L (RSE 19.6)
    lq <- log(42.6); label("Apparent inter-compartmental clearance of tramadol Q at ISAR = 0 (L/h)") # Table 3 final model 'Q , l/h' = 0.0426 printed; x 1000 = 42.6 L/h (RSE 6.09)
    lvp <- log(379); label("Apparent peripheral volume of tramadol V2/F at CrCL = 0 (L)") # Table 3 final model 'V 2 / F , l' = 0.379 printed; x 1000 = 379 L (RSE 35.5)
    lcl <- log(6.04); label("Apparent total clearance of tramadol CL/F at CrCL = 0 (L/h)") # Table 3 final model 'CL/ F , l/h' = 0.00604 printed; x 1000 = 6.04 L/h (RSE 6.71)

    # Kt is a first-order formation rate constant acting on the tramadol
    # central amount (Table 3 footnote 'Kt first-order rate constant for
    # tramadol metabolism to O-desmethyltramadol'; Results 3.3 'a first-order
    # metabolism rate constant'); the 'l/h' unit label in Table 3 is a
    # misprint, so the value is used as printed in 1/h with no x 1000. It feeds
    # the Fm-scaled apparent ODT pool and does not deplete tramadol.
    lkmet <- log(0.0492); label("First-order formation rate constant of ODT from the tramadol central amount Kt (1/h)") # Table 3 final model 'K t , l/h' = 0.0492 (RSE 21.8)
    lcl_m1 <- log(143); label("Apparent clearance of O-desmethyltramadol CLm/F (L/h)") # Table 3 final model 'CL m / F , l/h' = 0.143 printed; x 1000 = 143 L/h (RSE 11.9)

    # Covariate effects (uncentred exponential; see covariateData notes).
    e_score_isar_q <- 0.255; label("ISAR effect on Q, exponential coefficient (1/score unit)") # Table 3 final model 'ISAR effect on Q' = 0.255 (RSE 7.49)
    e_crcl_vp <- 0.0119; label("CrCL effect on V2/F, exponential coefficient (1/(mL/min))") # Table 3 final model 'CrCL effect on V 2 / F' = 0.0119 (RSE 33.6)
    e_crcl_cl <- 0.00498; label("CrCL effect on CL/F, exponential coefficient (1/(mL/min))") # Table 3 final model 'CrCL effect on CL/F' = 0.00498 (RSE 14.3)

    # Between-subject variability: log-normal (Methods Eq. 1). Table 3
    # footnote: 'Omega between-subject variability presented as standard
    # deviation', so each variance below is the printed omega squared.
    etalka ~ 0.142884 # Table 3 final model 'Omega K a' SD 0.378, squared (RSE 87.2)
    etalvc ~ 0.544644 # Table 3 final model 'Omega V 1 / F' SD 0.738, squared (RSE 20.5)
    etalq ~ 0.00047961 # Table 3 final model 'Omega Q' SD 0.0219, squared (RSE 52.3)
    etalvp ~ 0.948676 # Table 3 final model 'Omega V 2 / F' SD 0.974, squared (RSE 22.9)
    etalcl ~ 0.00436921 # Table 3 final model 'Omega CL/ F' SD 0.0661, squared (RSE 46.9)
    etalkmet ~ 0.616225 # Table 3 final model 'Omega K t' SD 0.785, squared (RSE 20.7)
    etalcl_m1 ~ 0.046656 # Table 3 final model 'Omega CL m / F' SD 0.216, squared (RSE 47.2)

    # Residual error (Results 3.3: proportional for tramadol, combined for
    # ODT). The '(%)' labels on the Table 3 residual rows are misprints: b is
    # a proportional fraction and the ODT a is in ng/mL.
    propSd <- 0.157; label("Proportional residual error for tramadol (fraction)") # Table 3 final model 'Tramadol, b (%)' = 0.157 (RSE 8.82)
    addSd_m1 <- 5.19; label("Additive residual error for O-desmethyltramadol (ng/mL)") # Table 3 final model 'O -desmethyltramadol, a (%)' = 5.19 (RSE 20.5)
    propSd_m1 <- 0.189; label("Proportional residual error for O-desmethyltramadol (fraction)") # Table 3 final model 'O -desmethyltramadol, b (%)' = 0.189 (RSE 12.8)
  })

  model({
    # --- Individual parameters (Methods Eq. 1, log-normal BSV) -------------
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)
    # Uncentred Monolix-default continuous covariate model:
    # log(theta_i) = log(theta_pop) + beta * COV + eta_i
    q <- exp(lq + e_score_isar_q * SCORE_ISAR + etalq)
    vp <- exp(lvp + e_crcl_vp * CRCL + etalvp)
    cl <- exp(lcl + e_crcl_cl * CRCL + etalcl)
    kmet <- exp(lkmet + etalkmet)
    cl_m1 <- exp(lcl_m1 + etalcl_m1)

    # --- Micro-constants ----------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    # ODT shares the tramadol central volume (see compartmentData notes).
    kel_m1 <- cl_m1 / vc

    # --- ODEs (Figure 1) ----------------------------------------------------
    # CL/F is the total tramadol clearance, so ODT formation carries no loss
    # term in d/dt(central): central_m1 is an apparent pool scaled by the
    # unidentifiable fraction metabolised (Methods 2.4.2, Vm/(F km) and
    # CLm/(F km)).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_m1) <- kmet * central - kel_m1 * central_m1

    # --- Observation --------------------------------------------------------
    # Amount in mg over volume in L is mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc_m1 <- 1000 * central_m1 / vc
    Cc ~ prop(propSd)
    # Monolix 'combined1' error, y = f + (a + b * f) * e, the Monolix default
    # combined form (the paper names only a 'combined residual model').
    Cc_m1 ~ add(addSd_m1) + prop(propSd_m1) + combined1()
  })
}
