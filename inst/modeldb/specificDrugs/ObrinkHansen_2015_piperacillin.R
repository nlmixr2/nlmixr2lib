ObrinkHansen_2015_piperacillin <- function() {
  description <- "Two-compartment population PK model for piperacillin in critically ill adults with septic shock (Obrink-Hansen 2015); linear first-order elimination with an additive linear effect of plasma creatinine on clearance, IIV on CL and central volume, and a proportional residual error."
  reference <- paste(
    "Obrink-Hansen K, Juul RV, Storgaard M, Thomsen MK, Hardlei TF,",
    "Brock B, Kreilgaard M, Gjedsted J.",
    "Population pharmacokinetics of piperacillin in the early phase of",
    "septic shock: does standard dosing result in therapeutic plasma",
    "concentrations?",
    "Antimicrob Agents Chemother. 2015;59(11):7018-7026.",
    "doi:10.1128/AAC.01347-15.",
    sep = " "
  )
  vignette <- "ObrinkHansen_2015_piperacillin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CREAT = list(
      description        = "Plasma creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at the level measured on the day of piperacillin sampling.",
        "Enters CL additively (linear-deviation form):",
        "CL_i = (TVCL + beta_Pcrea * (CREAT - 170)) * exp(eta_CL),",
        "with beta_Pcrea = -0.011 (L/h)/(umol/L) (Obrink-Hansen 2015 Table 2",
        "and Table 2 footnote). Reference 170 umol/L is the cohort median",
        "(Table 1, IQR 119-282 umol/L; observed range 53-446 umol/L).",
        "Source paper column 'p-creatinine'."
      ),
      source_name        = "p-creatinine"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 15,
    n_studies      = 1,
    age_range      = "59-79 years (IQR); cohort min/max not reported",
    age_median     = "66 years",
    weight_range   = "70.2-95 kg (IQR); cohort min/max not reported",
    weight_median  = "80 kg",
    sex_female_pct = 27,
    race_ethnicity = "Not reported",
    disease_state  = "Critically ill adults with known or suspected septic shock requiring noradrenaline infusion, treated empirically with piperacillin-tazobactam; patients on renal replacement therapy and those under 18 years were excluded",
    dose_range     = "4 g piperacillin (with 0.5 g tazobactam) as a 3-min IV infusion every 8 h; samples drawn during the third consecutive dosing interval",
    regions        = "Denmark (single-centre intensive care unit, Aarhus University Hospital, Skejby; ClinicalTrials.gov NCT02306928)",
    apache_score_median = "19 (IQR 14-23)",
    sofa_score_median   = "9 (IQR 7-10)",
    aki_pct             = 67,
    creat_range         = "53-446 umol/L observed; cohort median 170 umol/L (IQR 119-282)",
    albumin_median      = "30 g/L (IQR 27-32)",
    notes          = paste(
      "15 critically ill adults sampled prospectively September 2014 -",
      "January 2015 during the third consecutive piperacillin-tazobactam",
      "(4 g/0.5 g) dose. Eight free-piperacillin plasma concentrations per",
      "patient (pre-dose; 10, 20, 30 min; 1, 2, 4, 8 h post-dose) measured by",
      "UHPLC after 30 kDa ultrafiltration. Baseline characteristics in",
      "Obrink-Hansen 2015 Table 1."
    )
  )

  ini({
    # ===== Structural disposition (Obrink-Hansen 2015 Table 2 final 2-cmt model) =====
    lcl <- log(3.6);  label("Typical clearance at CREAT = 170 umol/L (L/h)")        # Table 2: CL 3.6 L/h (RSE 15.7%)
    lvc <- log(7.3);  label("Typical central volume of distribution (L)")           # Table 2: V1 7.3 L (RSE 11.8%)
    lq  <- log(6.58); label("Typical inter-compartmental clearance (L/h)")          # Table 2: Q 6.58 L/h (RSE 16.4%)
    lvp <- log(3.9);  label("Typical peripheral volume of distribution (L)")        # Table 2: V2 3.9 L (RSE 9.7%)

    # ===== Covariate effect (Obrink-Hansen 2015 Table 2 + Table 2 footnote) =====
    # Additive linear form on CL in absolute units (L/h)/(umol/L), NOT the
    # multiplicative [1 + beta*(X - Xmed)] form described generically in
    # Methods. The Table 2 footnote is explicit: CL_i = CL + beta_Pcrea *
    # (CREAT - 170). Higher plasma creatinine (worse renal function) lowers
    # CL; the negative coefficient reproduces this.
    e_creat_cl <- -0.011; label("Additive effect of plasma creatinine on CL ((L/h)/(umol/L))")  # Table 2: beta_Pcrea = -0.011 (RSE 11.9%)

    # ===== Reference covariate values (Obrink-Hansen 2015 Table 1 cohort median) =====
    creat_ref <- 170; label("Reference plasma creatinine (umol/L, cohort median)")  # Table 1 + Table 2 footnote: median 170 umol/L

    # ===== IIV (Obrink-Hansen 2015 Table 2) =====
    # Exponential (log-normal) IIV; CV%-to-variance: omega^2 = log(1 + CV^2).
    # 71.2% -> log(1 + 0.712^2) = log(1.50694) = 0.41001
    # 57.8% -> log(1 + 0.578^2) = log(1.33408) = 0.28818
    etalcl ~ 0.41001  # Table 2: IIV on CL 71.2% CV
    etalvc ~ 0.28818  # Table 2: IIV on V1 57.8% CV
    # Q and V2 had no IIV reported (paper: "interindividual variability (IIV)
    # adequately described as variance in clearance and central volume of
    # distribution").

    # ===== Residual error (Obrink-Hansen 2015 Table 2) =====
    propSd <- 0.147; label("Proportional residual error (fraction)")  # Table 2: proportional error 14.7% (RSE 14.4%)
  })

  model({
    # ----- Individual PK parameters -----
    # Additive linear creatinine effect on CL (Table 2 footnote):
    #   CL_i = (TVCL + beta_Pcrea * (CREAT - 170)) * exp(eta_CL).
    # eta_CL is applied multiplicatively on the linear scale (log-normal IIV).
    cl <- (exp(lcl) + e_creat_cl * (CREAT - creat_ref)) * exp(etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- ODE system -----
    # IV piperacillin-tazobactam dosed directly into the central compartment
    # (3-min infusion in the source study; no depot or absorption phase for
    # the systemic-disposition model). Two patients in the dataset showed an
    # extended absorption phase that was handled individually with a
    # first-order ka + lag time to avoid biasing systemic-parameter estimates
    # (Results, "PK/PD analysis"); the systemic parameters in Table 2 are
    # the population values reported for all 15 patients and apply unchanged
    # under direct IV dosing.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # ----- Output -----
    # Dose in mg, vc in L -> mg/L plasma piperacillin concentration. The
    # UHPLC assay measured the free (unbound) drug after 30 kDa
    # ultrafiltration, so Cc is interpreted as free piperacillin
    # concentration in plasma (Methods, "UHPLC analysis").
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
