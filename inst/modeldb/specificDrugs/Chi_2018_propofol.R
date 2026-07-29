Chi_2018_propofol <- function() {
  description <- paste(
    "Two-compartment population PK model for propofol target-controlled IV",
    "infusion in Chinese adults with hepatic insufficiency undergoing",
    "elective liver transplantation, with additive body-weight effect on",
    "clearance and power Child-Turcotte-Pugh score effect on peripheral",
    "volume (Chi 2018 final regression model). Typical-value-only model:",
    "the source paper reports the six final-model THETAs but provides no",
    "OMEGA (IIV), no SIGMA (residual error), and no GOF / VPC, so all",
    "etas are fixed at zero and no residual error term is included. See",
    "vignette Assumptions and deviations for the resulting limitations on",
    "VPC-style validation and the recommendation to consult the",
    "modellib('Ye_2012_propofol') companion (when extracted) for a",
    "fully-reported 3-compartment propofol popPK fit in a larger",
    "Chinese-multicenter cohort that shares the Chi 2018 first author."
  )
  reference <- paste(
    "Chi X, Pan J, Cai J, Luo G, Li S, Yuan D, Rui J, Chen W, Hei Z.",
    "Pharmacokinetic Analysis of Propofol Target-Controlled Infusion",
    "Models in Chinese Patients with Hepatic Insufficiency.",
    "Med Sci Monit 2018; 24:6925-6933.",
    "doi:10.12659/MSM.910103.",
    "Modeling methodology follows Ye HB, Li JH, Rui JZ et al,",
    "Propofol pharmacokinetics in China: A multicentric study,",
    "Indian J Pharmacol 2012; 44(3):393-7;",
    "pmid:22701253; this Chi 2018 paper re-estimates structural",
    "parameters in a 32-patient hepatic-insufficiency sub-cohort with",
    "Child-Turcotte-Pugh stratification.",
    sep = " "
  )
  vignette <- "Chi_2018_propofol"
  units <- list(
    time          = "minute",
    dosing        = "mg",
    concentration = "ug/mL"
  )
  # Dose units mg; central / vc has units mg/L = ug/mL, matching the
  # Chi 2018 propofol plasma concentration units (Methods, Anesthesia:
  # propofol plasma concentrations were set at 3 ug/mL).

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. Enters CL via the additive linear form",
        "CL (L/min) = 0.737 + 0.0163 * WT per Chi 2018 page 6930 final",
        "regression model. The intercept is non-zero at WT = 0, so this is",
        "not a multiplicative or allometric covariate effect: the slope",
        "0.0163 L/min/kg is applied directly on the linear scale and added",
        "to exp(lcl). Cohort group-mean body weights 58.6 +/- 11.2 kg (CTP",
        "A), 66.1 +/- 16.4 kg (CTP B), 63.1 +/- 8.4 kg (CTP C) per Table 1;",
        "pooled range approximately 40-100 kg."
      ),
      source_name        = "BW"
    ),
    CTP_SCORE = list(
      description        = "Composite Child-Turcotte-Pugh score (integer 5-15 points)",
      units              = "points",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. Sum of five subscores (total bilirubin,",
        "serum albumin, INR / prothrombin time, ascites severity, hepatic",
        "encephalopathy severity). Class A = 5-6, Class B = 7-9, Class C",
        "= 10-15. Enters V2 via the power form",
        "V2 (L) = 52.9 * (CTP_SCORE / 5)^(-0.848) per Chi 2018 page 6930",
        "final regression model; the reference value 5 is the Class A",
        "lower bound. Cohort stratified A / B / C with n = 11 / 10 / 11",
        "(Methods, Study population). The negative exponent means worse",
        "hepatic function (higher CTP_SCORE) is associated with a smaller",
        "peripheral volume; the discussion attributes the joint covariate",
        "pattern to shunted portal blood flow and altered tissue",
        "redistribution in cirrhosis."
      ),
      source_name        = "CTP"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 32L,
    n_studies        = 1L,
    age_range        = "18-65 years (inclusion); group means 43.7 +/- 7.4 (CTP A), 48.3 +/- 8.5 (CTP B), 46.4 +/- 8.7 (CTP C)",
    weight_range     = "Group means 58.6 +/- 11.2 kg (CTP A), 66.1 +/- 16.4 kg (CTP B), 63.1 +/- 8.4 kg (CTP C)",
    sex_female_pct   = 100 * 5 / 32,
    race_ethnicity   = "Chinese (Han majority; not separately reported)",
    disease_state    = paste(
      "Chinese adults with hepatic insufficiency scheduled for elective",
      "liver transplantation, ASA II-IV. Underlying liver disease was",
      "cirrhosis in 18 of 32 (56%) and hepatic carcinoma in 14 of 32",
      "(44%). Stratified by Child-Turcotte-Pugh (CTP) class: A n = 11",
      "(34.4%), B n = 10 (31.2%), C n = 11 (34.4%)."
    ),
    dose_range       = paste(
      "Propofol target-controlled IV infusion via Diprifusor TCI pump",
      "(P6003, Alaris, USA) with Marsh parameters; 3 ug/mL plasma target",
      "for the 30-minute induction window, then discontinued for 30",
      "minutes of washout sampling."
    ),
    regions          = "China (Third Affiliated Hospital, Sun Yat-sen University, Guangzhou; trial registration ChiCTR-OCH-12002255)",
    notes            = paste(
      "Single-center prospective observational study at the Third",
      "Affiliated Hospital of Sun Yat-sen University, conducted between",
      "May 2014 and March 2016 in patients undergoing elective liver",
      "transplantation. Demographics in Table 1. 32 enrolled; 416",
      "pre-administration blank-control samples (13 per patient) and",
      "matched PK sampling at 1, 2, 5, 10, 20, 30 min during the propofol",
      "TCI infusion plus 1, 2, 5, 10, 20, 30 min after discontinuation",
      "(13 samples per subject for the PK fit). Propofol assayed by",
      "reversed-phase HPLC with fluorescence detection (Fan 2006, Knibbe",
      "1998 references). NONMEM Version V, Level 1.1 (double precision);",
      "covariate selection by forward inclusion / backward elimination",
      "with delta-OFV threshold 6.64 (alpha = 0.01, df = 1); bootstrap",
      "validation with Wings for NONMEM. The paper reports only the six",
      "typical-value THETAs from the final regression model; OMEGA,",
      "SIGMA, GOF / VPC, and individual parameter estimates are not",
      "provided in the published text. The structural-model form follows",
      "the methodology of Ye et al. 2012 (Indian J Pharmacol 44:393-7)",
      "but Chi 2018 fits a 2-compartment model with covariates on CL and",
      "V2 rather than Ye 2012's 3-compartment model with covariates on",
      "V1 and Q3, so OMEGA / SIGMA values are not transferable between",
      "the two papers. See vignette Assumptions and deviations."
    )
  )

  ini({
    # Structural parameters from Chi 2018 page 6930 final regression
    # model. The full equation set written out in the paper is:
    #   CL = theta_CL + BW * theta_BW              (additive linear in BW)
    #   V1 = theta_V1
    #   Q  = theta_Q
    #   V2 = theta_V2 * (CTP / 5)^theta_CTP        (power, normalised to CTP = 5)
    # All six THETAs are wrapped in fixed() because the paper provides no
    # standard errors / RSE / confidence intervals for the final-model
    # estimates and the structural form does not separate "estimated" from
    # "structural-anchor" parameters; treating the values as fixed point
    # estimates is the conservative encoding (see vignette Assumptions and
    # deviations).

    # CL: additive linear in body weight. The intercept value 0.737 L/min
    # is the CL at WT = 0 (extrapolation; the cohort minimum is roughly
    # 40 kg). Stored as log(0.737) on the canonical lcl scale; the WT
    # slope 0.0163 L/min/kg is applied on the linear scale in model().
    lcl       <- fixed(log(0.737))   ; label("CL intercept at WT = 0 (L/min)")             # Chi 2018 p.6930 theta_CL = 0.737 L/min (additive linear form)
    e_wt_cl   <- fixed(0.0163)       ; label("Additive WT slope on CL (L/min/kg)")         # Chi 2018 p.6930 theta_BW = 0.0163 L/min/kg
    lvc       <- fixed(log(9.94))    ; label("Central volume of distribution V1 (L)")      # Chi 2018 p.6930 theta_V1 = 9.94 L
    lq        <- fixed(log(1.2))     ; label("Intercompartmental clearance Q (L/min)")     # Chi 2018 p.6930 theta_Q  = 1.2  L/min
    lvp       <- fixed(log(52.9))    ; label("Peripheral volume V2 at CTP_SCORE = 5 (L)")  # Chi 2018 p.6930 theta_V2 = 52.9 L
    e_ctp_score_vp <- fixed(-0.848)  ; label("Power exponent on (CTP_SCORE / 5) for V2")   # Chi 2018 p.6930 theta_CTP = -0.848

    # Inter-individual variability. Chi 2018 ran NONMEM final-model
    # estimation with Wings for NONMEM bootstrap validation in a 32-patient
    # cohort but does not report any OMEGA values in the published text;
    # operator decision (sidecar request-002 q1 = C) is to fix every eta
    # to zero rather than borrow from Ye 2012 or insert placeholders. The
    # etas are declared structurally so the typical-value-only intent is
    # encoded in the ini() block (Oniki 2018 nafld_risk precedent for the
    # `~ fixed(0)` pattern).
    etalcl ~ fixed(0)
    etalvc ~ fixed(0)
    etalq  ~ fixed(0)
    etalvp ~ fixed(0)

    # Residual error. Chi 2018 does not report any SIGMA values; per
    # operator decision (sidecar request-002 q1 = C, with the explicit
    # modification "do NOT add a propSd placeholder") no residual error
    # term is included. The model is therefore a deterministic
    # typical-value predictor; VPC-style validation is not possible. See
    # vignette Assumptions and deviations.
  })

  model({
    # Individual clearance: additive linear in body weight per Chi 2018
    # page 6930 (CL = theta_CL + BW * theta_BW). The intercept is stored on
    # the log scale (lcl) per nlmixr2lib convention; etalcl is fixed to 0
    # so exp(etalcl) = 1 and the typical-value structural form
    # CL = 0.737 + 0.0163 * WT is reproduced exactly. Precedent for the
    # `(exp(lcl) + e * (cov - ref)) * exp(etalcl)` additive-on-linear
    # encoding is Royer 2010 HuHMFG1 (AST on CL) and Weatherley 2018
    # fosdagrocorat (AGE on CL).
    cl <- (exp(lcl) + e_wt_cl * WT) * exp(etalcl)

    # Central volume - no covariate in Chi 2018 final model.
    vc <- exp(lvc + etalvc)

    # Intercompartmental clearance - no covariate in Chi 2018 final model.
    q  <- exp(lq + etalq)

    # Peripheral volume: power effect of CTP score normalised to the
    # Class A border value 5. The negative exponent -0.848 means V2
    # decreases monotonically as hepatic dysfunction worsens
    # (Chi 2018 Discussion attributes this to portal-systemic shunting
    # and impaired peripheral redistribution in cirrhosis).
    vp <- exp(lvp + etalvp) * (CTP_SCORE / 5)^e_ctp_score_vp

    # Micro-constants for the two-compartment central-disposition.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV with linear elimination from the central
    # compartment. Propofol TCI infusion enters `central` directly via
    # the event table's rate column; no depot / first-order absorption
    # is modelled because propofol is administered IV in Chi 2018
    # (Methods, Anesthesia).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Plasma propofol concentration. Dose units mg, vc units L; central /
    # vc has units mg/L = ug/mL, matching the Chi 2018 propofol plasma
    # concentration units (3 ug/mL TCI target; Methods, Anesthesia).
    # No residual-error declaration (Cc ~ ...) is included by operator
    # decision; the model is a deterministic typical-value predictor.
    Cc <- central / vc
  })
}
