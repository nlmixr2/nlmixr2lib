Jeong_2024_levocetirizine <- function() {
  description <- "Co-linked population PK-PD model for a single oral 5 mg levocetirizine immediate-release tablet in 40 healthy Korean adults (24 men, 16 women) from a two-way crossover bioequivalence study (Jeong 2024). PK is a two-compartment disposition model (central + peripheral1) fed by a chain of five sequential first-order absorption compartments (depot -> transit1 -> transit2 -> transit3 -> transit4 -> central, with distinct rate constants Ka1 to Ka5) preceded by an absorption lag time Tlag; volumes and clearances are apparent (Vc/F, CLc/F, Vp/F, CLp/F). Retained covariates: body surface area (power), serum albumin (power) and female sex (linear fractional decrease) on Vp/F, body surface area (power) on CLp/F, and body weight (power) on the fifth absorption rate constant Ka5, each continuous covariate centred on the cohort value. Residual error is log-additive, i.e. log-normal on the linear concentration scale. PD is two parallel indirect-response (input-inhibition) models for the histamine-induced skin wheal and flare sizes, expressed as a percentage of the pre-dose response (baseline kin / kout = 100), in which plasma levocetirizine inhibits the zero-order production rate through a sigmoid Imax function. The PD parameters were estimated from published mean data in children aged 6 to 11 years (Simons 2005) and carry no inter-individual variability."
  reference <- paste(
    "Jeong SH, Jang JH, Lee YB.",
    "Is Gender an Important Factor in the Precision Medicine Approach to Levocetirizine?",
    "Pharmaceutics. 2024;16(1):146.",
    "doi:10.3390/pharmaceutics16010146.",
    "The wheal and flare observations used to fit the pharmacodynamic sub-model were taken from",
    "Simons FER, Simons KJ.",
    "Levocetirizine: pharmacokinetics and pharmacodynamics in children age 6 to 11 years.",
    "J Allergy Clin Immunol. 2005;116(2):355-361. doi:10.1016/j.jaci.2005.04.010.",
    sep = " "
  )
  vignette <- "Jeong_2024_levocetirizine"

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
    # units$concentration documents the plasma levocetirizine output Cc. The
    # two PD outputs `wheal` and `flare` are skin-response SIZES expressed as a
    # percentage of the pre-dose histamine-induced response (100 = no
    # inhibition, 0 = complete suppression), not concentrations.
  )

  paper_specific_compartments <- c("wheal", "flare")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot = list(analyte = "levocetirizine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "levocetirizine", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "levocetirizine", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "levocetirizine", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "levocetirizine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "levocetirizine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "levocetirizine", units = "mg", specimen = "plasma", verified = TRUE),
    wheal = list(
      analyte = "histamine-induced skin wheal size",
      units = "percent of pre-dose response",
      specimen = "not applicable",
      verified = TRUE
    ),
    flare = list(
      analyte = "histamine-induced skin flare size",
      units = "percent of pre-dose response",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed. Supplementary Information 7 Eq. 5 prints",
        "Ka5 = tvKa5 * (Weight / median Weight)^dKa5dWeight * exp(eta_Ka5). The cohort",
        "median weight is not reported (Table S1 gives only mean (SD) 66.97 (12.59) kg), so",
        "the cohort mean 66.97 kg is used as the centring value (see vignette Assumptions).",
        "Table S7 footnote c attributes the weight effect to Ka4, but Eq. 5, Table S6",
        "model 5, Results 3.2 and the Conclusions all place it on Ka5.",
        sep = " "
      ),
      source_name = "Weight (Jeong 2024 Table S7 'dKa5dWeight'; Supplementary Information 7)"
    ),
    BSA = list(
      description = "Body surface area by the Mosteller formula, sqrt(height (cm) * body weight (kg) / 3600)",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed at baseline. Supplementary Information 7 Eqs. 9 and 10 apply",
        "(BSA / median BSA)^exponent to Vp/F and CLp/F. The cohort median is not reported",
        "(Table S1 gives only mean (SD) 1.76 (0.20) m^2), so the cohort mean 1.76 m^2 is used",
        "as the centring value (see vignette Assumptions). Supplementary Information 1 states",
        "the Mosteller derivation.",
        sep = " "
      ),
      source_name = "BSA (Jeong 2024 Table S7 'dVp/FdBSA', 'dCLp/FdBSA'; Table S1; Supplementary Information 1 and 7)"
    ),
    ALB = list(
      description = "Serum albumin concentration at baseline",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed at baseline. Supplementary Information 7 Eq. 9 applies",
        "(Albumin / median Albumin)^dVp/FdAlbumin to Vp/F. The paper reports albumin in g/dL",
        "(Table S1 mean (SD) 4.70 (0.24) g/dL); the canonical column is g/L, and the ratio",
        "form is unit-invariant once the centre is converted too, so the centre is",
        "4.70 g/dL = 47.0 g/L. The cohort median is not reported, so the cohort mean is used",
        "as the centring value (see vignette Assumptions).",
        sep = " "
      ),
      source_name = "Albumin (Jeong 2024 Table S7 'dVp/FdAlbumin'; Table S1)"
    ),
    SEXF = list(
      description = "Biological sex indicator; 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Time-fixed. Supplementary Information 7 Eq. 9 prints the sex term as",
        "[1 + dVp/FdGender * (If Male = 0 or Female = 1)], i.e. LINEAR in",
        "(1 + e_sexf_vp * SEXF) with e_sexf_vp = -0.18, so the typical female Vp/F is 82",
        "percent of the male value at equal BSA and albumin.",
        sep = " "
      ),
      source_name = "Gender (Jeong 2024 Table S7 'dVp/FdGender'; Supplementary Information 7)"
    )
  )

  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation, not body-surface-area normalised",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Renal-function indices (Cockcroft-Gault CrCL on actual, ideal and BMI-selected",
        "weight; MDRD and CKD-EPI eGFR) were screened and not retained (Discussion: the",
        "healthy cohort's CrCL of 80-190 mL/min gave no detectable effect on V/F or CL/F).",
        "Cohort mean (SD) 125.91 (23.70) mL/min (Table S1).",
        sep = " "
      ),
      source_name = "Creatinine clearance based on weight (Jeong 2024 Table S1; Discussion)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 40L,
    n_studies = 1L,
    age_range = "19 years or older by inclusion criterion; mean (SD) 28.38 (8.30) years (men 28.58 (8.08), women 28.06 (8.88))",
    weight_range = "mean (SD) 66.97 (12.59) kg (men 71.58 (12.81), women 60.05 (8.68)); inclusion weight > 50 kg (men) and > 45 kg (women)",
    height_range = "mean (SD) 168.33 (9.10) cm (men 173.41 (7.50), women 160.70 (5.08))",
    bsa_range = "mean (SD) 1.76 (0.20) m^2 (men 1.85 (0.19), women 1.63 (0.13))",
    albumin_range = "mean (SD) 4.70 (0.24) g/dL (men 4.79 (0.23), women 4.57 (0.17))",
    bmi_range = "18-30 kg/m^2 by inclusion criterion; mean (SD) 23.53 (3.36)",
    sex_female_pct = 40,
    race_ethnicity = "Korean",
    disease_state = "Healthy adult volunteers",
    dose_range = "Single oral 5 mg levocetirizine immediate-release tablet (reference formulation, UCB Pharma)",
    regions = "Republic of Korea (single centre; Ministry of Food and Drug Safety approval MB22-013)",
    renal_function = "Normal; Cockcroft-Gault creatinine clearance mean (SD) 125.91 (23.70) mL/min, MDRD eGFR 108.25 (15.85) mL/min/1.73 m^2 (Table S1)",
    notes = paste(
      "Randomised, single-dose, open-label, two-way crossover bioequivalence study with a",
      "7-day washout (Supplementary Information 3); only the reference-formulation period",
      "was modelled. Subjects fasted for more than 10 h and took the tablet with 150 mL",
      "of water. Plasma sampling at 0, 0.17, 0.33, 0.5, 0.75, 1, 1.25, 1.5, 2, 3, 4, 6, 8,",
      "12, 24 and 36 h post-dose; LLOQ 0.5 ng/mL. Demographics from Table S1. The",
      "pharmacodynamic sub-model was NOT fitted to data from this cohort: wheal and flare",
      "observations after 5 mg levocetirizine in children aged 6 to 11 years were taken",
      "from Simons 2005 (Methods 2.4) and fitted with the adult PK structure.",
      sep = " "
    )
  )

  ini({
    # ========================================================================
    # PK. Two-compartment disposition fed by five sequential first-order
    # absorption compartments with an absorption lag time on the dosing depot.
    # Structure from Jeong 2024 Results 3.2, Table S5 (model 02-05-02) and
    # Figure S2; estimates from Table S7; equations from Supplementary
    # Information 7. All disposition parameters are apparent (divided by the
    # unknown oral bioavailability F), so no bioavailability term is estimated.
    #
    # The chain is dose -> depot -(Ka1)-> transit1 -(Ka2)-> transit2 -(Ka3)->
    # transit3 -(Ka4)-> transit4 -(Ka5)-> central, matching Supplementary
    # Information 7 "Ka1 ... dosing depot-to-depot 1 ... Ka5 ... depot 4-to-
    # central compartment". The five rate constants are estimated separately,
    # so this is a sequential absorption chain rather than a Savic equal-ktr
    # transit chain.
    # ========================================================================
    lka1 <- log(9.87); label("First sequential absorption rate constant Ka1, depot to transit1 (1/h)") # Jeong 2024 Table S7: tv Ka1 = 9.87 1/h (SE 1.53, CV 15.53%; bootstrap median 9.86, 95% CI 7.07-12.74)
    lka2 <- log(28.41); label("Second sequential absorption rate constant Ka2, transit1 to transit2 (1/h)") # Jeong 2024 Table S7: tv Ka2 = 28.41 1/h (SE 5.39, CV 18.96%; bootstrap median 28.60, 95% CI 20.87-40.82)
    lka3 <- log(11.53); label("Third sequential absorption rate constant Ka3, transit2 to transit3 (1/h)") # Jeong 2024 Table S7: tv Ka3 = 11.53 1/h (SE 1.87, CV 16.20%; bootstrap median 11.45, 95% CI 7.91-15.35)
    lka4 <- log(11.29); label("Fourth sequential absorption rate constant Ka4, transit3 to transit4 (1/h)") # Jeong 2024 Table S7: tv Ka4 = 11.29 1/h (SE 1.95, CV 17.26%; bootstrap median 10.71, 95% CI 8.29-15.67)
    lka5 <- log(0.90); label("Fifth sequential absorption rate constant Ka5, transit4 to central, at the centring weight (1/h)") # Jeong 2024 Table S7: tv Ka5 = 0.90 1/h (SE 0.15, CV 16.85%; bootstrap median 0.90, 95% CI 0.62-1.18)
    ltlag <- log(0.08); label("Absorption lag time Tlag (h)") # Jeong 2024 Table S7: tv Tlag = 0.08 h (SE 0.01, CV 10.15%; bootstrap median 0.08, 95% CI 0.07-0.10)
    lvc <- log(1.07); label("Apparent central volume of distribution Vc/F (L)") # Jeong 2024 Table S7: tv Vc/F = 1.07 L (SE 0.28, CV 26.00%; bootstrap median 1.03, 95% CI 0.63-1.44)
    lcl <- log(2.94); label("Apparent central clearance CLc/F (L/h)") # Jeong 2024 Table S7: tv CLc/F = 2.94 L/h (SE 0.07, CV 2.49%; bootstrap median 2.95, 95% CI 2.79-3.07)
    lvp <- log(28.12); label("Apparent peripheral volume of distribution Vp/F in men at the centring BSA and albumin (L)") # Jeong 2024 Table S7: tv Vp/F = 28.12 L (SE 1.07, CV 3.82%; bootstrap median 28.02, 95% CI 26.44-30.70)
    lq <- log(13.80); label("Apparent inter-compartmental clearance CLp/F at the centring BSA (L/h)") # Jeong 2024 Table S7: tv CLp/F = 13.80 L/h (SE 2.24, CV 16.21%; bootstrap median 13.88, 95% CI 9.57-18.16)

    # ---- Retained covariate effects ----------------------------------------
    # Supplementary Information 7:
    #   Ka5  = tvKa5  * (WT / median WT)^dKa5dWeight * exp(eta_Ka5)
    #   Vp/F = tvVp/F * (BSA / median BSA)^dVp/FdBSA
    #                 * (ALB / median ALB)^dVp/FdAlbumin
    #                 * [1 + dVp/FdGender * SEXF] * exp(eta_Vp/F)
    #   CLp/F = tvCLp/F * (BSA / median BSA)^dCLp/FdBSA * exp(eta_CLp/F)
    # The continuous effects are power models on the centred covariate; the
    # sex effect is a LINEAR fractional change.
    e_wt_ka5 <- 0.69; label("Power exponent on (WT / 66.97 kg) for Ka5 (unitless)") # Jeong 2024 Table S7: dKa5dWeight = 0.69 (SE 0.28, CV 40.62%; bootstrap median 0.56, 95% CI 0.51-1.41)
    e_bsa_vp <- 1.20; label("Power exponent on (BSA / 1.76 m^2) for Vp/F (unitless)") # Jeong 2024 Table S7: dVp/FdBSA = 1.20 (SE 0.14, CV 11.59%; bootstrap median 1.18, 95% CI 0.94-1.54)
    e_alb_vp <- -0.79; label("Power exponent on (ALB / 47.0 g/L) for Vp/F (unitless)") # Jeong 2024 Table S7: dVp/FdAlbumin = -0.79 (SE 0.19, CV 23.90%; bootstrap median -0.79, 95% CI -1.23 to -0.40)
    e_sexf_vp <- -0.18; label("Fractional change in Vp/F for women, applied as 1 + e * SEXF (unitless)") # Jeong 2024 Table S7: dVp/FdGender = -0.18 (SE 0.02, CV 11.35%; bootstrap median -0.18, 95% CI -0.24 to -0.15)
    e_bsa_q <- 1.94; label("Power exponent on (BSA / 1.76 m^2) for CLp/F (unitless)") # Jeong 2024 Table S7: dCLp/FdBSA = 1.94 (SE 0.49, CV 25.29%; bootstrap median 1.80, 95% CI 1.21-3.29)

    # ========================================================================
    # PD. Two parallel indirect-response (input-inhibition) models of identical
    # structure, one for the histamine-induced wheal size and one for the
    # flare size, each expressed as a percentage of the pre-dose response.
    # Supplementary Information 8 Eq. 11:
    #   dR/dt = Kin * (1 - Imax * C^gamma / (IC50^gamma + C^gamma)) - Kout * R
    # with C the plasma levocetirizine concentration. Estimates from Table S8.
    # The shape factor gamma is a Hill exponent in a sigmoid Imax function, so
    # it takes the canonical name hill. Imax is a fraction bounded by 1, so it
    # is carried untransformed. No IIV, covariate or residual error was
    # reported for any PD parameter (Discussion).
    # ========================================================================
    lkin_wheal <- log(87.9); label("Zero-order production rate of the wheal response Kin (%/h)") # Jeong 2024 Table S8, wheal: Kin = 87.9 %/h
    lkout_wheal <- log(0.879); label("First-order loss rate of the wheal response Kout (1/h)") # Jeong 2024 Table S8, wheal: Kout = 0.879 1/h
    imax_wheal <- 0.92; label("Maximum fractional inhibition of wheal production Imax (fraction)") # Jeong 2024 Table S8, wheal: Imax = 0.92
    lic50_wheal <- log(45.5); label("Plasma levocetirizine concentration giving half-maximal inhibition of wheal production IC50 (ng/mL)") # Jeong 2024 Table S8, wheal: IC50 = 45.5 ng/mL
    lhill_wheal <- log(1.11); label("Hill exponent gamma of the wheal inhibition (unitless)") # Jeong 2024 Table S8, wheal: gamma = 1.11

    lkin_flare <- log(99.9); label("Zero-order production rate of the flare response Kin (%/h)") # Jeong 2024 Table S8, flare: Kin = 99.9 %/h
    lkout_flare <- log(0.999); label("First-order loss rate of the flare response Kout (1/h)") # Jeong 2024 Table S8, flare: Kout = 0.999 1/h
    imax_flare <- 0.999; label("Maximum fractional inhibition of flare production Imax (fraction)") # Jeong 2024 Table S8, flare: Imax = 0.999
    lic50_flare <- log(10.5); label("Plasma levocetirizine concentration giving half-maximal inhibition of flare production IC50 (ng/mL)") # Jeong 2024 Table S8, flare: IC50 = 10.5 ng/mL
    lhill_flare <- log(2.51); label("Hill exponent gamma of the flare inhibition (unitless)") # Jeong 2024 Table S8, flare: gamma = 2.51

    # ========================================================================
    # Inter-individual variability. Supplementary Information 7 specifies an
    # exponential model P_i = P_tv * exp(eta_i) with eta_i ~ N(0, omega^2), and
    # Table S7 reports omega^2 directly, so the variances are transcribed
    # without conversion. Table S5 (models 02-05-02-01 to -10) retains IIV on
    # all ten PK parameters.
    # ========================================================================
    etalvc ~ 0.89 # Jeong 2024 Table S7: omega^2 Vc/F = 0.89 (SE 0.29, CV 33.13%; bootstrap median 0.87, 95% CI 0.05-1.70)
    etalcl ~ 0.03 # Jeong 2024 Table S7: omega^2 CLc/F = 0.03 (SE 0.01, CV 23.45%; bootstrap median 0.03, 95% CI 0.02-0.05)
    etaltlag ~ 0.37 # Jeong 2024 Table S7: omega^2 Tlag = 0.37 (SE 0.11, CV 29.46%; bootstrap median 0.37, 95% CI 0.10-0.65)
    etalvp ~ fixed(0) # Jeong 2024 Table S7: omega^2 Vp/F reported as 0.00 (SE 0.00; bootstrap 95% CI 0.00-0.00). Retained (Table S5 model 02-05-02-03: removing it worsened -2LL by 512.90) but below the two-decimal reporting precision, so carried as an explicit zero rather than invented.
    etalq ~ fixed(0) # Jeong 2024 Table S7: omega^2 CLp/F reported as 0.00 (SE 0.00; bootstrap 95% CI 0.00-0.00). Retained (Table S5 model 02-05-02-04: removing it worsened -2LL by 420.05) but below the two-decimal reporting precision.
    etalka1 ~ 2.09 # Jeong 2024 Table S7: omega^2 Ka1 = 2.09 (SE 0.29, CV 13.72%; bootstrap median 2.06, 95% CI 1.33-2.79)
    etalka2 ~ 4.14 # Jeong 2024 Table S7: omega^2 Ka2 = 4.14 (SE 0.98, CV 23.65%; bootstrap median 3.94, 95% CI 2.08-5.80)
    etalka3 ~ 2.36 # Jeong 2024 Table S7: omega^2 Ka3 = 2.36 (SE 0.43, CV 18.05%; bootstrap median 2.31, 95% CI 1.49-3.12)
    etalka4 ~ 2.11 # Jeong 2024 Table S7: omega^2 Ka4 = 2.11 (SE 0.41, CV 19.54%; bootstrap median 2.13, 95% CI 1.20-3.07)
    etalka5 ~ 0.02 # Jeong 2024 Table S7: omega^2 Ka5 = 0.02 (SE 0.01, CV 34.76%; bootstrap median 0.02, 95% CI 0.01-0.03)

    # ========================================================================
    # Residual error. Table S5 selects the "log additive" residual model
    # (02-05-02), which in Phoenix NLME is an additive error on the natural-log
    # scale, log(Cobs) = log(Cpred) + eps -- a log-normal residual on the
    # linear concentration scale, so it maps to `~ lnorm(...)`. Phoenix reports
    # the residual parameter as a standard deviation.
    #
    # The wheal and flare endpoints have NO reported residual error (Table S8
    # lists only Kin, Kout, Imax, IC50 and gamma), so they are carried as
    # fixed(0). See the vignette Assumptions and deviations section.
    # ========================================================================
    expSd <- 0.11; label("Log-scale residual SD on levocetirizine Cc (log-additive error)") # Jeong 2024 Table S7: epsilon = 0.11 (SE 0.01, CV 11.04%; bootstrap median 0.11, 95% CI 0.09-0.13)
    addSd_wheal <- fixed(0); label("Additive residual SD on the wheal size (% of pre-dose response)") # Jeong 2024 Table S8 reports no residual-error term for the pharmacodynamic model
    addSd_flare <- fixed(0); label("Additive residual SD on the flare size (% of pre-dose response)") # Jeong 2024 Table S8 reports no residual-error term for the pharmacodynamic model
  })

  model({
    # ---- 1. Individual PK parameters ---------------------------------------
    # Supplementary Information 7 covariate forms. The centring values are the
    # Table S1 cohort means (66.97 kg, 1.76 m^2, 4.70 g/dL = 47.0 g/L) because
    # the cohort medians named in the equations are not reported.
    ka1 <- exp(lka1 + etalka1)
    ka2 <- exp(lka2 + etalka2)
    ka3 <- exp(lka3 + etalka3)
    ka4 <- exp(lka4 + etalka4)
    ka5 <- exp(lka5 + etalka5) * (WT / 66.97)^e_wt_ka5
    tlag <- exp(ltlag + etaltlag)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)
    vp <- exp(lvp + etalvp) * (BSA / 1.76)^e_bsa_vp * (ALB / 47.0)^e_alb_vp * (1 + e_sexf_vp * SEXF)
    q <- exp(lq + etalq) * (BSA / 1.76)^e_bsa_q

    # ---- 2. Micro-constants ------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 3. ODE system -----------------------------------------------------
    d/dt(depot) <- -ka1 * depot
    d/dt(transit1) <- ka1 * depot - ka2 * transit1
    d/dt(transit2) <- ka2 * transit1 - ka3 * transit2
    d/dt(transit3) <- ka3 * transit2 - ka4 * transit3
    d/dt(transit4) <- ka4 * transit3 - ka5 * transit4
    d/dt(central) <- ka5 * transit4 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- 4. Lag time on the dosing depot (Figure S2) -----------------------
    alag(depot) <- tlag

    # ---- 5. Plasma concentration -------------------------------------------
    # Central amount is in mg and vc in L, so central/vc is mg/L; multiply by
    # 1000 to report ng/mL, the unit in which IC50 was estimated.
    Cc <- (central / vc) * 1000

    # ---- 6. Indirect-response PD (Supplementary Information 8 Eq. 11) ------
    kin_wheal <- exp(lkin_wheal)
    kout_wheal <- exp(lkout_wheal)
    ic50_wheal <- exp(lic50_wheal)
    hill_wheal <- exp(lhill_wheal)
    kin_flare <- exp(lkin_flare)
    kout_flare <- exp(lkout_flare)
    ic50_flare <- exp(lic50_flare)
    hill_flare <- exp(lhill_flare)

    inh_wheal <- imax_wheal * Cc^hill_wheal / (ic50_wheal^hill_wheal + Cc^hill_wheal)
    inh_flare <- imax_flare * Cc^hill_flare / (ic50_flare^hill_flare + Cc^hill_flare)

    d/dt(wheal) <- kin_wheal * (1 - inh_wheal) - kout_wheal * wheal
    d/dt(flare) <- kin_flare * (1 - inh_flare) - kout_flare * flare
    # Drug-free turnover steady state kin / kout = 100, the pre-dose response
    # to which Table 1 and Figures 3-5 normalise the wheal and flare sizes.
    wheal(0) <- kin_wheal / kout_wheal
    flare(0) <- kin_flare / kout_flare

    # ---- 7. Residual error -------------------------------------------------
    Cc ~ lnorm(expSd)
    wheal ~ add(addSd_wheal)
    flare ~ add(addSd_flare)
  })
}
