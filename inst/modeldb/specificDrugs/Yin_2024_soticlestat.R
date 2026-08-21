Yin_2024_soticlestat <- function() {
  description <- paste(
    "Joint population PK / CH24H enzyme-occupancy (EO) / 24S-hydroxycholesterol (24HC)",
    "pharmacodynamic model for soticlestat (TAK-935), a cholesterol 24-hydroxylase inhibitor,",
    "in healthy volunteers and patients with developmental and epileptic encephalopathies (DEE)",
    "including Dravet syndrome and Lennox-Gastaut syndrome (Yin 2024).",
    "Two-compartment PK with first-order absorption and an absorption lag time; a",
    "plasma-to-brain effect compartment whose concentration drives both a fixed sigmoid Emax",
    "CH24H enzyme-occupancy read-out and a semimechanistic sigmoid Imax inhibitory",
    "indirect-response turnover model for 24HC.",
    "Covariates: formulation on lag time and ka; dose on ka, Q, Vp and relative",
    "bioavailability; BMI on ka (participants aged <= 18 y only); strong CYP3A-inducing",
    "antiseizure comedication and Chinese descent on ka; patient (vs healthy-volunteer)",
    "status on CL; Japanese and Chinese descent on Q; eGFR on Vp; body weight and",
    "alpha-1-acid glycoprotein on relative bioavailability; and age (hockey stick below",
    "17.5 y), alpha-1-acid glycoprotein and body weight on baseline 24HC."
  )
  reference <- paste(
    "Yin W, Facius A, Asgharnejad M, Lahu G, Vakilynejad M.",
    "Population pharmacokinetics, enzyme occupancy, and pharmacodynamic modeling of",
    "soticlestat in patients with developmental and epileptic encephalopathies.",
    "Clin Transl Sci. 2024;17(3):e13722. doi:10.1111/cts.13722.",
    "The CH24H enzyme-occupancy sub-model (ke0, Emax, EC50, gamma) is carried forward",
    "unchanged and FIXED from the healthy-volunteer model of Yin W, Facius A, Asgharnejad M,",
    "Wang S, Rosen L, Bhattacharya A, Lahu G, Vakilynejad M.",
    "Modeling and simulation of soticlestat pharmacokinetics, brain enzyme occupancy, and",
    "pharmacodynamics in healthy volunteers. Clin Transl Sci. 2023;16(8):1422-1434.",
    "doi:10.1111/cts.13517.",
    sep = " "
  )
  vignette <- "Yin_2024_soticlestat"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    DOSE = list(
      description        = "Soticlestat dose level for the current regimen",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Empirical dose-nonlinearity covariate entering as (DOSE/300)^exponent on ka, Q, Vp",
        "and relative bioavailability. All four exponents were FIXED during estimation",
        "(Table 1a 'Dose effect, exponent' rows are flagged 'Fixed'; NONMEM $THETA 13-16 carry",
        "FIX). Reference 300 mg. Studied dose range 15-1350 mg (Table S1)."
      ),
      source_name        = "DOSE"
    ),
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on relative bioavailability, reference 65.9 kg, and on baseline 24HC,",
        "reference 55.25 kg. The authors estimated the weight effect on Frel rather than",
        "fixing an allometric exponent on CL/V (paper section 'PopPK model')."
      ),
      source_name        = "WEIGHT"
    ),
    BMI = list(
      description        = "Baseline body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on ka, reference 23.6 kg/m^2, applied ONLY to participants aged",
        "<= 18 years (NONMEM 'IF( AGE.LE.18 )'). For older participants the factor is 1."
      ),
      source_name        = "BMI"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Two distinct roles. (1) Gate for the BMI-on-ka effect, which applies only when",
        "AGE <= 18 y. (2) Hockey-stick power effect on baseline 24HC applied only when",
        "AGE < 17.5 y (the breakpoint itself was estimated); above the breakpoint the",
        "factor is 1."
      ),
      source_name        = "AGE"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate (BSA-normalized renal function)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on peripheral volume, reference 162.4 mL/min/1.73 m^2.",
        "The source reports eGFR; the estimating equation is not stated in the paper.",
        "Cohort mean (SD) 166 (54.1) mL/min/1.73 m^2 in the popPK analysis set (Table S2)."
      ),
      source_name        = "EGFR"
    ),
    AAG = list(
      description        = "Baseline alpha-1-acid glycoprotein concentration",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reported in mg/dL by the source, NOT in the register's default g/L; the reference",
        "value is 20 mg/dL (= 0.2 g/L). Power effect on relative bioavailability and on",
        "baseline 24HC. Missing values were imputed at the observed median of 20 mg/dL,",
        "flagged by AAG_MISSING."
      ),
      source_name        = "A1AGLP"
    ),
    AAG_MISSING = list(
      description        = "Indicator that baseline alpha-1-acid glycoprotein was not measured and was imputed",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (AAG measured)",
      notes              = paste(
        "1 = AAG was not collected in the participant's trial and was imputed at the median",
        "of 20 mg/dL; 0 = AAG measured. Does not shift the typical value: it inflates the",
        "between-subject SD of relative bioavailability by the factor",
        "(1 + e_aag_missing_fdepot), so imputed participants carry extra unexplained",
        "variability. Studies TAK-935-1003, TAK-935-2001, TAK-935-2002 and TAK-935-18-002",
        "have AAG imputed for all participants (Table S2 footnote a)."
      ),
      source_name        = "A1AGLP_MISS"
    ),
    FORM_TABLET = list(
      description        = "Tablet (vs oral solution) formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral solution, the typical-value reference)",
      notes              = paste(
        "1 = any non-oral-solution formulation (tablet, including tablets administered via",
        "G-tube / PEG tube); 0 = oral solution. NONMEM encodes this as FORM != 1, with",
        "FORM = 1 the oral solution. Percent effects on both the absorption lag time",
        "(+23.8%) and ka (-43.7%). The original model represented this delay with a",
        "separate delay compartment; the present model replaced it with the",
        "formulation-on-lag-time covariate."
      ),
      source_name        = "FORM"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant (vs DEE patient) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = paste(
        "1 (healthy volunteer) is the TYPICAL-VALUE reference in this paper, which is the",
        "complement of the register default; the -22.8% CL effect is carried by patients",
        "(DIS_HEALTHY = 0)"
      ),
      notes              = paste(
        "The source column is PATIENT (1 = patient with a DEE, 0 = healthy volunteer), so",
        "DIS_HEALTHY = 1 - PATIENT. The published typical CL of 4.2 L/h is the",
        "HEALTHY-VOLUNTEER value; the model applies the patient effect on (1 - DIS_HEALTHY)",
        "so that the published anchor is preserved. Patient cohort comprises DEE, Dravet",
        "syndrome, Lennox-Gastaut syndrome, 15q duplication syndrome and CDKL5 deficiency",
        "disorder (Table S2)."
      ),
      source_name        = "PATIENT"
    ),
    RACE_JAPANESE = list(
      description        = "Japanese ethnic-background indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Japanese)",
      notes              = paste(
        "Percent effect on inter-compartmental clearance. NONMEM encodes this as",
        "COUNTRY == 2 within the Asian subpopulation. 24/218 (11%) of the popPK analysis",
        "set (Table S2)."
      ),
      source_name        = "COUNTRY (= 2)"
    ),
    RACE_CHINESE = list(
      description        = "Chinese ethnic-background indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Chinese)",
      notes              = paste(
        "Percent effects on ka and on inter-compartmental clearance. NONMEM encodes this as",
        "COUNTRY == 4 within the Asian subpopulation. 20/218 (9%) of the popPK analysis set",
        "(Table S2)."
      ),
      source_name        = "COUNTRY (= 4)"
    ),
    CONMED_EIAED = list(
      description        = "Concomitant enzyme-inducing antiseizure medication indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no strong CYP3A-inducing antiseizure comedication)",
      notes              = paste(
        "1 = concomitant treatment with an antiseizure medication grouped by the authors as",
        "a strong CYP3A enzyme inducer (carbamazepine, phenobarbital or phenytoin);",
        "0 otherwise. Percent effect on ka."
      ),
      source_name        = "AED_STR_IND"
    )
  )

  covariatesDataExcluded <- list(
    RACE_ASIAN = list(
      description = "Asian ethnic-background indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Present in the final NONMEM control stream as a covariate on peripheral volume",
        "(V3~ASIAN), but its coefficient is '0 FIX', making the term an exact no-op; the",
        "row does not appear in Table 1a. The retained ethnic-background effects are the",
        "Japanese and Chinese indicators. Omitted from model() rather than encoded as a",
        "fixed-zero multiplier."
      )
    ),
    FORM_CRUSHED_TABLET = list(
      description = "Crushed-tablet administration indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Carried in the NONMEM $INPUT and $TABLE records (CRUSHED) but never referenced in",
        "$PK or $ERROR of either final model, and not reported in Table 1. Screened but not",
        "retained."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as an individual characteristic in the covariate analysis; not retained in either final model."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened as a protein-binding factor; not retained in either final model."
    ),
    BILI = list(
      description = "Total bilirubin",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened as a liver-function factor (popPK model only); not retained."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte  = "soticlestat",
      units    = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte  = "soticlestat",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte  = "soticlestat",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    effect = list(
      analyte  = "soticlestat",
      units    = "ng/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    hc24 = list(
      analyte  = "24S-hydroxycholesterol",
      units    = "ng/mL",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 306,
    n_studies      = 8,
    age_range      = "paediatric through adult; mean (SD) 24.1 (14.5) years in the popPK set and 21.8 (14.6) years in the PK/EO/PD set",
    weight_range   = "mean (SD) 56.6 (24.4) kg in the popPK set and 52.3 (24.6) kg in the PK/EO/PD set",
    race_ethnicity = c(White = 67, Black = 11, `Asian (Japan)` = 11, `Asian (China)` = 9, Other = 2),
    disease_state  = "healthy volunteers and patients with developmental and epileptic encephalopathies (DEE, Dravet syndrome, Lennox-Gastaut syndrome, 15q duplication syndrome, CDKL5 deficiency disorder)",
    dose_range     = "15-1350 mg oral soticlestat, single dose or 100-300 mg b.i.d.",
    notes          = paste(
      "Two overlapping analysis sets from eight phase I/II trials (Table S1, Table S2).",
      "The population PK model used 218 individuals (110 healthy, 108 patients),",
      "3288 soticlestat concentrations and 8732 dosing events. The PK/EO/PD model used",
      "306 individuals (132 healthy, 174 patients), 2621 plasma 24HC concentrations and",
      "8703 dosing events; 17 of 2625 observations with |CWRES| > 5 were excluded during",
      "model development. n_subjects records the larger (PK/EO/PD) analysis set.",
      "Race percentages are from the PK/EO/PD column of Table S2 and are rounded, so they",
      "do not sum to exactly 100."
    )
  )

  ini({
    # ---- Structural PK, Table 1a (a) PopPK model, 'Estimate' column ----
    lka     <- log(8.39)   ; label("Absorption rate constant (ka, 1/h)")                     # Table 1a, Absorption rate (ka), TV
    lcl     <- log(4.2)    ; label("Linear elimination clearance in a healthy volunteer (CL, L/h)")  # Table 1a, Elimination clearance (CL), TV
    lvc     <- log(3.01)   ; label("Central volume of distribution (Vc, L)")                 # Table 1a, Central volume (Vc), TV
    lq      <- log(1.15)   ; label("Inter-compartmental clearance (Q, L/h)")                 # Table 1a, Distribution clearance (Q), TV
    lvp     <- log(7.8)    ; label("Peripheral volume of distribution (Vp, L)")              # Table 1a, Peripheral volume (Vp), TV
    ltlag   <- log(0.133)  ; label("Absorption lag time for the oral solution (h)")          # Table 1a, Lag time of the first compartment (ALAG1), TV
    lfdepot <- fixed(log(0.0216)) ; label("Relative bioavailability of the depot (unitless)")  # Table 1a, Bioavailability (F1), TV, Fixed; NONMEM $THETA 7 FIX

    # ---- PK covariate effects, Table 1a ----
    # Percent effects enter as (1 + theta/100); power effects as (cov/ref)^theta.
    e_form_tablet_tlag   <- 23.8            ; label("Non-oral-solution formulation effect on absorption lag time (%)")  # Table 1a, ALAG1 Non-OS formulation effect
    e_form_tablet_ka     <- -43.7           ; label("Non-oral-solution formulation effect on ka (%)")                   # Table 1a, ka Non-OS formulation effect
    e_dose_ka            <- fixed(-0.753)   ; label("Exponent on (DOSE/300) for ka (unitless)")                         # Table 1a, ka Dose effect, exponent, Fixed
    e_bmi_ka             <- 2.24            ; label("Exponent on (BMI/23.6) for ka, participants aged <= 18 y (unitless)")  # Table 1a, ka BMI effect, exponent
    e_conmed_eiaed_ka    <- -66.2           ; label("Strong CYP3A enzyme-inducing antiseizure comedication effect on ka (%)")  # Table 1a, ka Strong CYP3A enzyme inducer effect
    e_race_chinese_ka    <- -63.7           ; label("Chinese-descent effect on ka (%)")                                 # Table 1a, ka Chinese descent effect
    e_dis_healthy_cl     <- -22.8           ; label("Patient (DIS_HEALTHY = 0) effect on CL relative to a healthy volunteer (%)")  # Table 1a, CL Patient effect
    e_dose_q             <- fixed(-0.218)   ; label("Exponent on (DOSE/300) for Q (unitless)")                          # Table 1a, Q Dose effect, exponent, Fixed
    e_race_japanese_q    <- -42.7           ; label("Japanese-descent effect on Q (%)")                                 # Table 1a, Q Japanese descent effect
    e_race_chinese_q     <- -75.7           ; label("Chinese-descent effect on Q (%)")                                  # Table 1a, Q Chinese descent effect
    e_dose_vp            <- fixed(-0.214)   ; label("Exponent on (DOSE/300) for Vp (unitless)")                         # Table 1a, Vp Dose effect, exponent, Fixed
    e_crcl_vp            <- -0.406          ; label("Exponent on (eGFR/162.4) for Vp (unitless)")                       # Table 1a, Vp eGFR effect, exponent
    e_dose_fdepot        <- fixed(0.204)    ; label("Exponent on (DOSE/300) for relative bioavailability (unitless)")   # Table 1a, F1 Dose effect, exponent, Fixed
    e_wt_fdepot          <- -0.593          ; label("Exponent on (WT/65.9) for relative bioavailability (unitless)")    # Table 1a, F1 Body weight effect, exponent
    e_aag_fdepot         <- 0.544           ; label("Exponent on (AAG/20) for relative bioavailability (unitless)")     # Table 1a, F1 AGP effect, exponent
    e_aag_missing_fdepot <- 0.42            ; label("Inflation of the relative-bioavailability BSV scale when AAG was imputed (unitless)")  # Table 1a, F1 'BSV explained by AGP'; NONMEM $THETA 20

    # ---- PK between-subject variability, Table 1a ----
    # Table 1 footnote a states 'Reported on variance scale', but the panel (a)
    # values are standard deviations: (0.527)^2 = 0.2777 reproduces the control
    # stream's $OMEGA for F1 exactly, and the paper's own simulated steady-state
    # AUC 90% prediction intervals (Table 2 footnote: 2470/432 = 5.72) imply
    # omega_F1 = log(5.72)/(2*1.645) = 0.530, not 0.727 = sqrt(0.527).
    # Panel (b) IS on the variance scale. See the vignette Errata.
    etalka      ~ 1.0404    # Table 1a, BSV ka   = 1.02  (SD) -> variance 1.02^2
    etalq       ~ 0.190096  # Table 1a, BSV Q    = 0.436 (SD) -> variance 0.436^2
    etalvp      ~ 0.390625  # Table 1a, BSV Vp   = 0.625 (SD) -> variance 0.625^2
    etalfdepot  ~ 0.277729  # Table 1a, BSV F1   = 0.527 (SD) -> variance 0.527^2

    # ---- Soticlestat plasma residual error, Table 1a ----
    propSd <- 0.483            ; label("Proportional residual error for soticlestat plasma concentration (fraction)")  # Table 1a, Residual variability Proportional 48.3%
    addSd  <- fixed(0.001)     ; label("Additive residual error for soticlestat plasma concentration (ng/mL)")         # Table 1a, Residual variability Additive, Fixed

    # ---- CH24H enzyme-occupancy sub-model ----
    # Not re-estimated in this paper ('The PK/EO model was not updated because no
    # additional EO data were available'); all four values carry FIX in the
    # NONMEM $THETA of the final PK/EO/PD run and are inherited from Yin 2023.
    lke0     <- fixed(log(0.254)) ; label("Plasma-to-brain effect-site equilibration rate constant (ke0, 1/h)")  # Appendix S1 PK/EO/PD $THETA 3 KPLBR, FIX
    lemax    <- fixed(log(100))   ; label("Maximum CH24H enzyme occupancy (Emax, %)")                            # Appendix S1 PK/EO/PD $THETA 4 EMAX, FIX
    lec50    <- fixed(log(5.86))  ; label("Effect-site concentration for 50% CH24H enzyme occupancy (EC50, ng/mL)")  # Appendix S1 PK/EO/PD $THETA 5 EC50, FIX
    lhill_eo <- fixed(log(0.769)) ; label("Sigmoidicity of the CH24H enzyme-occupancy relationship (unitless)")  # Appendix S1 PK/EO/PD $THETA 6 EGAM, FIX

    # ---- 24HC turnover, Table 1b 'Original estimate' column ----
    lrbase          <- log(50.5)   ; label("Typical baseline 24HC concentration (ng/mL)")                # Table 1b, Baseline 24HC (BL24HC), TV
    e_age_bp_rbase  <- 17.5        ; label("Age breakpoint below which the age effect on baseline 24HC applies (years)")  # Table 1b, Age effect cutoff
    e_age_rbase     <- -0.511      ; label("Exponent on (AGE/17.5) for baseline 24HC below the breakpoint (unitless)")    # Table 1b, Age effect, exponent
    e_aag_rbase     <- 0.215       ; label("Exponent on (AAG/20) for baseline 24HC (unitless)")           # Table 1b, AGP effect, exponent
    e_wt_rbase      <- -0.256      ; label("Exponent on (WT/55.25) for baseline 24HC (unitless)")         # Table 1b, Body weight effect, exponent
    lkout           <- log(0.0199) ; label("First-order 24HC degradation rate constant (kout, 1/h)")      # Table 1b, 24HC degradation rate (kout), TV
    limax           <- log(92)     ; label("Maximum inhibition of 24HC production (Imax, %)")             # Table 1b, Maximum inhibition of 24HC production (Imax), TV
    lic50           <- log(9.85)   ; label("Effect-site concentration for 50% of maximum 24HC inhibition (IC50, ng/mL)")  # Table 1b, Effect-site concentration for 50% maximum effect (IC50), TV
    lhill_hc24      <- log(0.881)  ; label("Sigmoidicity of the 24HC inhibition relationship (unitless)")  # Table 1b, Shape parameter (gamma), TV
    # NOTE THE SCALE: I0 in the control stream is a PERCENT and enters as
    # (100 - I0)/100, whereas the canonical iplac is the FRACTION I0/100.
    # Both are zero here, but a future non-zero value must be divided by 100.
    iplac           <- fixed(0)    ; label("Constant fractional placebo inhibition of 24HC production (unitless)")  # Appendix S1 PK/EO/PD $THETA 9 I0 = 0 FIX; iplac = I0/100

    # ---- 24HC between-subject variability, Table 1b (variance scale) ----
    etalrbase ~ 0.0811   # Table 1b, BL24HC BSV
    etalic50  ~ 0.636    # Table 1b, IC50 BSV

    # ---- 24HC residual error, Table 1b ----
    # The proportional component was fixed to zero (NONMEM $THETA 1 '0 FIX'), so the
    # 24HC residual is purely additive.
    addSd_hc24 <- 3.91 ; label("Additive residual error for plasma 24HC (ng/mL)")  # Table 1b, Residual variability Additive
  })

  model({
    # ---------------------------------------------------------------------
    # 1. Individual PK parameters (final PopPK model, Appendix S1 $PK)
    # ---------------------------------------------------------------------
    # Absorption lag time; the FORM effect is a percent change (ALAG1 in NONMEM).
    tlag <- exp(ltlag) * (1 + e_form_tablet_tlag / 100 * FORM_TABLET)

    # The BMI effect on ka applies only to participants aged <= 18 y; the logical
    # returns 1/0 so the exponent collapses to 0 (factor 1) for older participants.
    ka <- exp(lka + etalka) *
      (1 + e_form_tablet_ka / 100 * FORM_TABLET) *
      (DOSE / 300)^e_dose_ka *
      (BMI / 23.6)^(e_bmi_ka * (AGE <= 18)) *
      (1 + e_conmed_eiaed_ka / 100 * CONMED_EIAED) *
      (1 + e_race_chinese_ka / 100 * RACE_CHINESE)

    # The typical CL is the healthy-volunteer value; the patient effect is carried
    # on the complement of DIS_HEALTHY (source column PATIENT = 1 - DIS_HEALTHY).
    cl <- exp(lcl) * (1 + e_dis_healthy_cl / 100 * (1 - DIS_HEALTHY))

    vc <- exp(lvc)

    q <- exp(lq + etalq) *
      (DOSE / 300)^e_dose_q *
      (1 + e_race_japanese_q / 100 * RACE_JAPANESE) *
      (1 + e_race_chinese_q / 100 * RACE_CHINESE)

    vp <- exp(lvp + etalvp) *
      (DOSE / 300)^e_dose_vp *
      (CRCL / 162.4)^e_crcl_vp

    # Relative bioavailability. The between-subject SD is inflated by
    # (1 + e_aag_missing_fdepot) for participants whose AAG was imputed.
    fdepot <- exp(lfdepot + (1 + e_aag_missing_fdepot * AAG_MISSING) * etalfdepot) *
      (DOSE / 300)^e_dose_fdepot *
      (WT / 65.9)^e_wt_fdepot *
      (AAG / 20)^e_aag_fdepot

    # ---------------------------------------------------------------------
    # 2. Effect-site and 24HC turnover parameters (Appendix S1 PK/EO/PD $PK)
    # ---------------------------------------------------------------------
    ke0     <- exp(lke0)
    emax    <- exp(lemax)
    ec50    <- exp(lec50)
    hill_eo <- exp(lhill_eo)

    # Hockey-stick age effect on baseline 24HC: applied only below the breakpoint.
    rbase <- exp(lrbase + etalrbase) *
      (AGE / e_age_bp_rbase)^(e_age_rbase * (AGE < e_age_bp_rbase)) *
      (AAG / 20)^e_aag_rbase *
      (WT / 55.25)^e_wt_rbase

    kout      <- exp(lkout)
    kin       <- rbase * kout
    imax      <- exp(limax)
    ic50      <- exp(lic50 + etalic50)
    hill_hc24 <- exp(lhill_hc24)

    # ---------------------------------------------------------------------
    # 3. ODE system (Appendix S1 PK/EO/PD $DES)
    # ---------------------------------------------------------------------
    # Dose in mg / volume in L gives mg/L; the factor 1000 converts to ng/mL and
    # reproduces the NONMEM scale factor S2 = V2/1000.
    Cc <- 1000 * central / vc

    # Effect-site (brain) concentration, driven by the plasma concentration.
    # NONMEM: DADT(4) = KPLBR * (A(2)/S2 - A(4)).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - cl / vc * central -
                          q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <-  q / vc * central - q / vp * peripheral1
    d/dt(effect)      <-  ke0 * (Cc - effect)

    # Fractional inhibition of 24HC synthesis. NONMEM guards the sigmoid with
    # 'IF( A(4) > 0 )'; max() reproduces that guard and keeps the fractional power
    # defined if the solver returns a small negative effect-site concentration.
    ce      <- max(effect, 0)
    drugEff <- imax * ce^hill_hc24 / (ce^hill_hc24 + ic50^hill_hc24) / 100
    eff     <- 1 - iplac - drugEff

    d/dt(hc24) <- kin * eff - kout * hc24
    hc24(0)    <- rbase

    f(depot)    <- fdepot
    alag(depot) <- tlag

    # ---------------------------------------------------------------------
    # 4. Derived read-outs and observations (Appendix S1 PK/EO/PD $ERROR)
    # ---------------------------------------------------------------------
    # CH24H enzyme occupancy in brain (%), a read-out of the same effect site.
    occ <- emax * ce^hill_eo / (ec50^hill_eo + ce^hill_eo)
    # Percent change from the individual's own baseline 24HC.
    hc24Chg <- 100 * (hc24 / rbase - 1)

    Cc   ~ prop(propSd) + add(addSd)
    hc24 ~ add(addSd_hc24)
  })
}
