Willmann_2024_elinzanetant <- function() {
  description <- paste(
    "Joint population PK/PD model for oral elinzanetant, its intravenous",
    "[13C5]-labelled tracer, and its three principal metabolites M30/34, M27 and",
    "M18/21, built on integrated phase I and II data from healthy volunteers and",
    "women with vasomotor symptoms associated with menopause (Willmann 2024).",
    "Elinzanetant and the primary metabolites M30/34 and M27 are two-compartment;",
    "the secondary metabolite M18/21 is one-compartment. Absorption is first order",
    "with a lag time; a high-fat breakfast and the hard-gel capsule formulation",
    "modify the absorption rate constant, the lag time and (high-fat breakfast",
    "only) the bioavailability. Elinzanetant clearance is split by fixed fractions",
    "into formation of M30/34 (0.6) and M27 (0.3), the remaining 0.1 going to",
    "unmodelled metabolites; the linear clearances of M30/34 and M27 both form",
    "M18/21. Each metabolite carries an additional saturable Michaelis-Menten",
    "elimination. All clearances and Vmax terms follow a 24 h cosine circadian",
    "modulation anchored on the clock time of the first dose, with a 49.2%",
    "relative amplitude for elinzanetant and 22.5% for the metabolites. NK-1 and",
    "NK-3 receptor occupancies are derived from the elinzanetant plasma",
    "concentration through Emax relationships with fixed EC50 values.",
    sep = " "
  )
  reference <- paste(
    "Willmann S, Lloyd A, Austin R, Joseph S, Solms A, Zhang Y,",
    "Schneider ARP, Frechen S, Schultze-Mosgau MH.",
    "Population pharmacokinetic-pharmacodynamic model of elinzanetant based on",
    "integrated clinical phase I and II data.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(12):2137-2149.",
    "doi:10.1002/psp4.13226",
    sep = " "
  )
  vignette <- "Willmann_2024_elinzanetant"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    FED_HIGHFAT = list(
      description = paste(
        "High-fat breakfast at the time of the oral elinzanetant dose",
        "(1 = high-fat breakfast, 0 = fasted or any other meal)"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = paste(
        "0 = source FED level 0 (fasted) or 1 (fed with a meal other than a",
        "high-fat breakfast)"
      ),
      notes = paste(
        "Per-dose-record indicator. The source data file codes feeding status as",
        "FED = 0 (fasted), 1 (fed, not a high-fat breakfast) or 2 (high-fat",
        "breakfast); Willmann 2024 Data S1 defines FED2 = 1 when FED = 2, else 0,",
        "so FED_HIGHFAT = FED2. Modifies Ka, the absorption lag time and the",
        "oral bioavailability."
      ),
      source_name = "FED2"
    ),
    FORM_ELZ_HARDGEL = list(
      description = paste(
        "Hard-gelatin capsule elinzanetant formulation indicator",
        "(1 = hard gel capsule, 0 = soft-gel capsule, aqueous suspension or",
        "intravenous solution)"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = paste(
        "0 = source FORM level 1, 2 or 4 (soft-gel capsule, aqueous suspension",
        "and the intravenous tracer)"
      ),
      notes = paste(
        "Per-dose-record indicator. Willmann 2024 Data S1 defines FORM3 = 1 when",
        "FORM = 3, else 0, so FORM_ELZ_HARDGEL = FORM3. FORM = 3 is the hard gel",
        "capsule used in studies 21677 (100 mg) and RELENT-1 (50-400 mg).",
        "Modifies Ka and the absorption lag time but not the extent of",
        "absorption (Willmann 2024 Discussion)."
      ),
      source_name = "FORM3"
    ),
    T_FIRSTDOSE = list(
      description = paste(
        "Wall-clock time of the first elinzanetant dose, as a decimal between 0",
        "and 24 h (e.g. 9:45 a.m. is 9.75)"
      ),
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Per-subject, time-fixed. Converts the model integration axis (time since",
        "the first dose) to wall-clock time so the 24 h circadian cosine on the",
        "clearance and Vmax terms can be evaluated: clock time = time +",
        "T_FIRSTDOSE. Willmann 2024 Equation 1 calls this TFD. Set to 9 for the",
        "paper's 9:00 a.m. morning-dosing simulation and 21 for the 9:00 p.m.",
        "evening-dosing simulation (Willmann 2024 Table 3)."
      ),
      source_name = "TFD"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "elinzanetant", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "elinzanetant", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "elinzanetant", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    central_c13 = list(
      analyte = "[13C5]-elinzanetant", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1_c13 = list(
      analyte = "[13C5]-elinzanetant", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    central_m3034 = list(
      analyte = "elinzanetant metabolite M30/34", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1_m3034 = list(
      analyte = "elinzanetant metabolite M30/34", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    central_m27 = list(
      analyte = "elinzanetant metabolite M27", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1_m27 = list(
      analyte = "elinzanetant metabolite M27", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    central_m1821 = list(
      analyte = "elinzanetant metabolite M18/21", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 195,
    n_studies = 7,
    age_range = "22-65 years",
    age_median = "55 years",
    weight_range = "51-103.6 kg",
    weight_median = "70 kg",
    height_median = "162 cm",
    bmi_median = "26.8 kg/m2",
    sex_female_pct = 92.8,
    race_ethnicity = c(White = 78.5, Black = 19.5, Other = 2),
    disease_state = paste(
      "70.3% healthy volunteers, 29.7% women with moderate-to-severe vasomotor",
      "symptoms associated with menopause"
    ),
    dose_range = paste(
      "25-400 mg orally (soft-gel capsule, hard gel capsule and aqueous",
      "suspension) plus a 100 ug [13C5]-elinzanetant intravenous microdose"
    ),
    notes = paste(
      "Model-building dataset: 195 subjects in seven phase I studies contributing",
      "8408 elinzanetant PK observations (5.7% below the lower limit of",
      "quantification, handled with the M3 method) and 2020 observations for each",
      "of the three metabolites from 52 of those subjects; demographics from",
      "Willmann 2024 Table S1. A separate validation dataset of 171 subjects in",
      "three further phase I/II studies (including 139 women with vasomotor",
      "symptoms in SWITCH-1) was used for external visual predictive checks and",
      "did not contribute to the parameter estimates; across both datasets 366",
      "subjects were studied, 336 female and 197 of them women with vasomotor",
      "symptoms. Baseline laboratory medians in the model-building dataset were",
      "aspartate transaminase 19 U/L, alanine transaminase 18 U/L, alkaline",
      "phosphatase 76 U/L, total bilirubin 0.5 mg/dL and albumin 4.3 g/dL.",
      "Vasomotor-symptom patient status versus healthy-volunteer status was one of",
      "the intrinsic factors screened post hoc on model-derived AUC(0-24)ss",
      "(Willmann 2024 Table S2, Figure 3) and was not retained in the model."
    )
  )

  # All entries below were screened post hoc by Willmann 2024 as ratios of the
  # geometric mean model-derived AUC(0-24)ss between dichotomised subgroups of
  # the 336 female subjects, and none was carried into the model. The ratios and
  # 90% CIs quoted are read from the Figure 3 forest plot, which prints them
  # numerically; the continuous covariates are split at the female-population
  # medians given in the Figure 3 caption (height 1.62 m, AST 21 U/L, ALT
  # 20 U/L, ALP 79.5 U/L, bilirubin 0.491 mg/dL, albumin 4.4 g/dL).
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units = "year", type = "continuous",
      notes = paste(
        "Screened post hoc at two cut-points (Willmann 2024 Table S2, Figure 3):",
        ">= 55 vs < 55 years gave a geometric mean AUC(0-24)ss ratio of 0.939",
        "[90% CI 0.869-1.017] (174 vs 162 subjects) and >= 60 vs < 60 years gave",
        "0.956 [0.860-1.061] (66 vs 270). Both CIs include unity, so age was not",
        "retained in the model."
      )
    ),
    WT = list(
      description = "Body weight at baseline",
      units = "kg", type = "continuous",
      notes = paste(
        "Screened post hoc at two cut-points (Willmann 2024 Figure 3):",
        ">= 60 vs < 60 kg gave a ratio of 0.945 [90% CI 0.851-1.055] (278 vs 58",
        "subjects) and >= 80 vs < 80 kg gave 0.934 [0.830-1.041] (74 vs 262).",
        "Both CIs include unity; not retained."
      )
    ),
    HT = list(
      description = "Body height at baseline",
      units = "cm", type = "continuous",
      notes = paste(
        "Screened post hoc above vs below the female-population median of 1.62 m",
        "(Willmann 2024 Figure 3): ratio 0.884 [90% CI 0.815-0.961] (168 vs 168",
        "subjects), i.e. taller women had approximately 12% lower elinzanetant",
        "AUC(0-24)ss and the CI excludes unity. The authors judged the effect to",
        "lie within the bioequivalence range and to be confounded with race",
        "(black women were significantly taller, median 1.65 vs 1.60 m,",
        "p < 0.001), so no height term was added to the model. Only the subgroup",
        "ratio is reported; no structural covariate coefficient is available."
      )
    ),
    BMI = list(
      description = "Body mass index at baseline",
      units = "kg/m2", type = "continuous",
      notes = paste(
        "Screened post hoc at two cut-points (Willmann 2024 Figure 3):",
        ">= 25 vs < 25 kg/m2 gave a ratio of 0.977 [90% CI 0.907-1.067] (231 vs",
        "105 subjects) and >= 30 vs < 30 kg/m2 gave 1.035 [0.916-1.153] (77 vs",
        "259). Both CIs include unity; not retained."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)", type = "binary",
      reference_category = "0 = male",
      notes = paste(
        "Metabolite AUC(0-24)ss was significantly lower in the 30 male than in",
        "the 336 female subjects (one-way ANOVA p between 1.2e-6 and 9.1e-7);",
        "the same trend for elinzanetant itself was not significant (p = 0.238).",
        "Because the target population is women with vasomotor symptoms, the",
        "remaining covariate screen was restricted to female subjects and no sex",
        "term was carried into the model (Willmann 2024 Results)."
      )
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units = "(binary)", type = "binary",
      reference_category = "0 = White or Other",
      notes = paste(
        "Screened post hoc within the female subpopulation (Willmann 2024",
        "Figure 3): black vs white-or-other gave a geometric mean AUC(0-24)ss",
        "ratio of 0.876 [90% CI 0.783-0.974] (65 vs 271 subjects), i.e. black",
        "women had approximately 12% lower elinzanetant exposure and the CI",
        "excludes unity. Black women were also significantly taller (median 1.65",
        "vs 1.60 m, p < 0.001), so the authors attributed the race and height",
        "signals to a common cause and retained neither. Only the subgroup ratio",
        "is reported; no structural covariate coefficient is available."
      )
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units = "(binary)", type = "binary",
      reference_category = "0 = Black or Other",
      notes = paste(
        "Screened post hoc alongside RACE_BLACK (Willmann 2024 Figure 3): white",
        "vs black-or-other gave a ratio of 1.147 [90% CI 1.026-1.269] (262 vs 74",
        "subjects), the mirror image of the RACE_BLACK contrast. Not retained."
      )
    ),
    AST = list(
      description = "Aspartate transaminase at baseline",
      units = "U/L", type = "continuous",
      notes = paste(
        "Screened post hoc above vs below the female-population median of",
        "21 U/L (Willmann 2024 Figure 3): ratio 1.045 [90% CI 0.966-1.140]",
        "(178 vs 158 subjects). CI includes unity; not retained."
      )
    ),
    ALT = list(
      description = "Alanine transaminase at baseline",
      units = "U/L", type = "continuous",
      notes = paste(
        "Screened post hoc above vs below the female-population median of",
        "20 U/L (Willmann 2024 Figure 3): ratio 1.022 [90% CI 0.942-1.104]",
        "(180 vs 156 subjects). CI includes unity; not retained."
      )
    ),
    ALP = list(
      description = "Alkaline phosphatase at baseline",
      units = "U/L", type = "continuous",
      notes = paste(
        "Screened post hoc above vs below the female-population median of",
        "79.5 U/L (Willmann 2024 Figure 3): ratio 0.946 [90% CI 0.876-1.021]",
        "(168 vs 168 subjects). CI includes unity; not retained."
      )
    ),
    TBILI_BASE = list(
      description = "Total bilirubin at baseline",
      units = "mg/dL", type = "continuous",
      notes = paste(
        "Screened post hoc above vs below the female-population median of",
        "0.491 mg/dL (Willmann 2024 Figure 3): ratio 1.147 [90% CI 1.063-1.239]",
        "(169 vs 167 subjects), i.e. women with bilirubin ABOVE the median had",
        "approximately 15% HIGHER elinzanetant exposure and the CI excludes",
        "unity. Note that the Willmann 2024 Results and Discussion prose states",
        "the opposite direction ('women with bilirubin levels below the median",
        "were found to have approximately 15% higher elinzanetant exposure' and",
        "'Low bilirubin levels were associated with a slightly (approximately",
        "15%) higher elinzanetant exposure'). The Figure 3 direction is adopted",
        "here because it is the numerically printed result and because it is the",
        "one consistent with the authors' own mechanistic reading, that the",
        "association 'might reflect the direct causality between liver function",
        "and the hepatic metabolism of elinzanetant' -- higher bilirubin implies",
        "poorer hepatic function and therefore higher exposure. Not retained in",
        "the model, so this discrepancy does not affect model predictions."
      )
    ),
    ALB = list(
      description = "Serum albumin at baseline",
      units = "g/dL", type = "continuous",
      notes = paste(
        "Screened post hoc above vs below the female-population median of",
        "4.4 g/dL (Willmann 2024 Figure 3): ratio 1.060 [90% CI 0.976-1.153]",
        "(187 vs 149 subjects). CI includes unity; not retained."
      )
    )
  )

  # The two routes Willmann 2024 actually administered: oral elinzanetant into
  # `depot`, and the 100 ug [13C5]-elinzanetant intravenous microdose of study
  # 21772 into `central_c13`. Declared explicitly because the automatic
  # detection in buildModelDb() would otherwise report `depot, central` --
  # advertising a direct intravenous elinzanetant route that the paper never
  # used, and omitting the tracer route that identifies absolute
  # bioavailability.
  dosing <- c("depot", "central_c13")

  ini({
    # ---- Elinzanetant disposition (Willmann 2024 Table 2, fixed effects) ----
    lcl <- log(7.26)    ; label("Elinzanetant clearance at the mid-point of the circadian variation (L/h)")  # Table 2 CLpop = 7.26 (RSE 4.94%, 95% CI 6.56-7.97)
    lvc <- log(23.7)    ; label("Elinzanetant central volume of distribution (L)")                            # Table 2 Vcpop = 23.7 (RSE 4.64%, 95% CI 21.6-25.9)
    lvp <- log(168)     ; label("Elinzanetant peripheral volume of distribution (L)")                         # Table 2 Vp = 168 (RSE 2.00%, 95% CI 161-175)
    lq  <- log(6.57)    ; label("Elinzanetant inter-compartmental clearance (L/h)")                           # Table 2 Q = 6.57 (RSE 2.56%, 95% CI 6.24-6.90)

    # ---- Elinzanetant absorption ----
    lka     <- log(2.37)  ; label("Elinzanetant absorption rate constant with reference feeding status and formulation (1/h)")  # Table 2 KaREF = 2.37 (RSE 7.22%, 95% CI 2.04-2.71); the unit "L/h" in the Results text is a typographical error, Table 2 gives h^-1
    ltlag   <- log(0.292) ; label("Elinzanetant absorption lag time with reference feeding status and formulation (h)")          # Table 2 ALAGREF = 0.292 (RSE 0.163%, 95% CI 0.291-0.293)
    lfdepot <- log(0.367) ; label("Absolute oral bioavailability of elinzanetant with reference feeding status (fraction)")      # Table 2 FREF = 0.367 (RSE 1.78%, 95% CI 0.354-0.380)

    # ---- Covariate effects on absorption (Willmann 2024 Data S1) ----
    # Ka, ALAG and F are each multiplied by (1 + <indicator> * <coefficient>).
    e_fed_highfat_ka         <- -0.911 ; label("Proportional change in Ka with a high-fat breakfast (fraction)")                 # Table 2 KaFED2 = -0.911 (RSE 0.409%, 95% CI -0.919 to -0.904); a 91% decrease
    e_form_elz_hardgel_ka    <- -0.264 ; label("Proportional change in Ka with the hard gel capsule (fraction)")                 # Table 2 KaFORM3 = -0.264 (RSE 14.6%, 95% CI -0.340 to -0.189); a 26% decrease
    e_fed_highfat_tlag       <-  0.545 ; label("Proportional change in the absorption lag time with a high-fat breakfast (fraction)")  # Table 2 ALFED2 = 0.545 (RSE 2.91%, 95% CI 0.514-0.576); a 55% increase
    e_form_elz_hardgel_tlag  <-  0.548 ; label("Proportional change in the absorption lag time with the hard gel capsule (fraction)")  # Table 2 ALFORM3 = 0.548 (RSE 0.880%, 95% CI 0.539-0.558); a 55% increase
    e_fed_highfat_fdepot     <- -0.227 ; label("Proportional change in oral bioavailability with a high-fat breakfast (fraction)")     # Table 2 FFED2 = -0.227 (RSE 8.59%, 95% CI -0.265 to -0.189); a 23% decrease

    # ---- M30/34 (primary metabolite, two-compartment) ----
    lcl_m3034   <- log(4.14) ; label("M30/34 linear clearance to M18/21 at the mid-point of the circadian variation (L/h)")  # Table 2 CL3034pop = 4.14 (RSE 5.44%, 95% CI 3.70-4.58)
    lvc_m3034   <- log(3.58) ; label("M30/34 central volume of distribution (L)")                                            # Table 2 Vc3034 = 3.58 (RSE 7.76%, 95% CI 3.04-4.13)
    lvp_m3034   <- log(15.6) ; label("M30/34 peripheral volume of distribution (L)")                                         # Table 2 Vp3034 = 15.6 (RSE 4.56%, 95% CI 14.2-17.0)
    lq_m3034    <- log(3.65) ; label("M30/34 inter-compartmental clearance (L/h)")                                           # Table 2 Q3034 = 3.65 (RSE 9.05%, 95% CI 3.00-4.29)
    lvmax_m3034 <- log(24.1) ; label("M30/34 saturable-elimination Vmax at the mid-point of the circadian variation (ug/L/h)")  # Table 2 Vmax3034 = 24.1 (RSE 10.8%, 95% CI 19.0-29.2)
    lkm_m3034   <- log(7.87) ; label("M30/34 saturable-elimination Michaelis-Menten constant (ug/L)")                        # Table 2 Km3034 = 7.87 (RSE 8.17%, 95% CI 6.61-9.13)

    # ---- M27 (primary metabolite, two-compartment) ----
    lcl_m27   <- log(4.11) ; label("M27 linear clearance to M18/21 at the mid-point of the circadian variation (L/h)")  # Table 2 CL27pop = 4.11 (RSE 5.85%, 95% CI 3.64-4.58)
    lvc_m27   <- log(1.85) ; label("M27 central volume of distribution (L)")                                            # Table 2 Vc27 = 1.85 (RSE 7.84%, 95% CI 1.56-2.13)
    lvp_m27   <- log(40.4) ; label("M27 peripheral volume of distribution (L)")                                         # Table 2 Vp27 = 40.4 (RSE 5.83%, 95% CI 35.7-45.0)
    lq_m27    <- log(4.45) ; label("M27 inter-compartmental clearance (L/h)")                                           # Table 2 Q27 = 4.45 (RSE 5.77%, 95% CI 3.94-4.95)
    lvmax_m27 <- log(42.6) ; label("M27 saturable-elimination Vmax at the mid-point of the circadian variation (ug/L/h)")  # Table 2 Vmax27 = 42.6 (RSE 12.6%, 95% CI 32.1-53.1)
    lkm_m27   <- log(17.7) ; label("M27 saturable-elimination Michaelis-Menten constant (ug/L)")                        # Table 2 Km27 = 17.7 (RSE 11.0%, 95% CI 13.9-21.5)

    # ---- M18/21 (secondary metabolite, one-compartment) ----
    lcl_m1821   <- log(9.10)  ; label("M18/21 linear clearance at the mid-point of the circadian variation (L/h)")  # Table 2 CL1821pop = 9.10 (RSE 6.14%, 95% CI 8.01-10.2)
    lvc_m1821   <- log(0.636) ; label("M18/21 central volume of distribution (L)")                                  # Table 2 Vc1821 = 0.636 (RSE 9.86%, 95% CI 0.513-0.759)
    lvmax_m1821 <- log(335)   ; label("M18/21 saturable-elimination Vmax at the mid-point of the circadian variation (ug/L/h)")  # Table 2 Vmax18/21 = 335 (RSE 13.0%, 95% CI 250-421)
    lkm_m1821   <- log(45.8)  ; label("M18/21 saturable-elimination Michaelis-Menten constant (ug/L)")              # Table 2 Km18/21 = 45.8 (RSE 13.0%, 95% CI 34.1-57.4)

    # ---- Fractions of elinzanetant clearance forming each primary metabolite ----
    # Willmann 2024 Data S1 "Information about PopPK model parametrization"
    # gives these as constants, not as estimated THETAs; they do not appear in
    # Table 2 and carry no standard error.
    fm_m3034 <- fixed(0.6) ; label("Fraction of elinzanetant clearance forming M30/34 (unitless)")  # Data S1: FMET3034 = 0.6
    fm_m27   <- fixed(0.3) ; label("Fraction of elinzanetant clearance forming M27 (unitless)")     # Data S1: FMET27 = 0.3

    # ---- Circadian modulation of the clearance and Vmax terms ----
    amp_diurnal       <- 0.492 ; label("Relative amplitude of the circadian cosine on elinzanetant clearance (fraction)")            # Table 2 AMP = 0.492 (RSE 3.59%, 95% CI 0.457-0.526)
    amp_diurnal_met   <- 0.225 ; label("Relative amplitude of the circadian cosine on the metabolite clearances and Vmax (fraction)")  # Table 2 AMPMET = 0.225 (RSE 6.39%, 95% CI 0.197-0.253)
    shift_diurnal     <- 8.05  ; label("Phase-shift parameter positioning the circadian peaks and troughs relative to clock time (h)")  # Table 2 theta_SHIFT = 8.05 (RSE 1.34%, 95% CI 7.84-8.26); places the clearance maximum near 16:00 and the minimum near 04:00

    # ---- Box-Cox shape parameter for the IIV on Vc ----
    # Petersson 2009 form; Willmann 2024 Data S1 states the transformation was
    # used because the untransformed eta3 distribution was skewed.
    boxcox_lvc <- 2.07 ; label("Box-Cox shape parameter lambda for the IIV transformation on elinzanetant Vc (unitless)")  # Table 2 lambda = 2.07 (RSE 9.51%, 95% CI 1.69-2.46)

    # ---- NK-1 / NK-3 receptor-occupancy Emax model ----
    # Neither EC50 was estimated in this analysis: the NK-1 value is taken from
    # the elinzanetant PET study (Willmann 2024 reference 14, repeated-dosing
    # arm) and the NK-3 value is that value multiplied by 10 on the assumption
    # that the in vitro 10-fold affinity difference is maintained in vivo.
    lec50_nk1 <- fixed(log(0.97)) ; label("Elinzanetant plasma EC50 for NK-1 receptor occupancy, from the PET study (ug/L)")  # Willmann 2024 PK/PD modeling section: in vivo EC50 for NK-1 = 0.97 ug/L
    lec50_nk3 <- fixed(log(9.7))  ; label("Elinzanetant plasma EC50 for NK-3 receptor occupancy, assumed 10-fold the NK-1 value (ug/L)")  # Willmann 2024 PK/PD modeling section: EC50 for NK-3 set to 9.7 ug/L

    # ---- Inter-individual variability (Willmann 2024 Table 2, OMEGA) ----
    # eta1 (elinzanetant CL), eta2 (shared across the three metabolite
    # clearances) and eta3 (elinzanetant Vc) form a full 3x3 block; eta4 (Ka)
    # and eta5 (ALAG) are diagonal. This matches the Omega matrix drawn in
    # Willmann 2024 Data S1, which shows a full 3x3 block on (CL, CLM, VC) and
    # zeros elsewhere. The lower-triangle order in the c() below is
    # (1,1); (2,1), (2,2); (3,1), (3,2), (3,3), and the six values are, from
    # Willmann 2024 Table 2 OMEGA:
    #   0.215 = ETA(1,1) omega^2, elinzanetant CL (RSE 12.1%), CV 49.0%
    #   0.224 = ETA(2,1) omega, covariance CL vs metabolite CL (RSE 19.3%)
    #   0.318 = ETA(2,2) omega^2, metabolite clearances (RSE 23.6%), CV 61.2%
    #   0.164 = ETA(3,1) omega, covariance Vc vs CL (RSE 19.2%)
    #   0.214 = ETA(3,2) omega, covariance Vc vs metabolite CL (RSE 24.2%)
    #   0.221 = ETA(3,3) omega^2, elinzanetant Vc (RSE 21.7%), CV 49.7%
    # The off-diagonal entries are covariances, not correlations: Data S1
    # states "the off-diagonal elements are the covariances". Comments cannot
    # be placed inside the c() below because rxode2 rewrites each in-block
    # comment to ";" while parsing labels, which is a syntax error.
    etalcl + etalcl_met + etalvc ~ c(
      0.215,
      0.224, 0.318,
      0.164, 0.214, 0.221
    )
    etalka   ~ 0.539          # ETA(4,4) omega^2 = 0.539 (RSE 14.5%); CV 84.5%
    etaltlag ~ 0.0989         # ETA(5,5) omega^2 = 0.0989 (RSE 10.5%); CV 32.2%

    # ---- Residual error (Willmann 2024 Table 2, SIGMA) ----
    # NONMEM fitted an additive error on log-transformed concentrations, which
    # is a log-normal residual in nlmixr2; each SD below is sqrt(sigma^2) on the
    # natural-log scale. The reported CV% equals 100*sqrt(exp(sigma^2) - 1).
    expSd       <- sqrt(0.378)  ; label("Elinzanetant log-scale residual SD (log ug/L)")          # EPS(1,1) sigma^2 = 0.378 (RSE 0.872%); CV 67.8%
    expSd_c13   <- sqrt(0.122)  ; label("[13C5]-elinzanetant log-scale residual SD (log ug/L)")   # EPS(2,2) sigma^2 = 0.122 (RSE 10.3%); CV 36.0%
    expSd_m3034 <- sqrt(0.106)  ; label("M30/34 log-scale residual SD (log ug/L)")                # EPS(3,3) sigma^2 = 0.106 (RSE 4.45%); CV 33.5%
    expSd_m27   <- sqrt(0.143)  ; label("M27 log-scale residual SD (log ug/L)")                   # EPS(4,4) sigma^2 = 0.143 (RSE 5.20%); CV 39.3%
    expSd_m1821 <- sqrt(0.141)  ; label("M18/21 log-scale residual SD (log ug/L)")                # EPS(5,5) sigma^2 = 0.141 (RSE 4.69%); CV 38.9%
  })

  model({
    # ---------------------------------------------------------------------
    # Unit bridge. Doses and all compartment amounts are in mg, while the
    # paper reports concentrations, Km and Vmax in ug/L.
    # ---------------------------------------------------------------------
    ugPerMg <- 1000

    # ---------------------------------------------------------------------
    # 1. Circadian modulation (Willmann 2024 Equation 1 and Data S1).
    #    P(t) = Pavg * (1 + AMP * cos((2*pi/24) * (t + TFD + theta_SHIFT)))
    #    with t the time since the first dose and TFD the clock time of that
    #    first dose, so (time + T_FIRSTDOSE) is the wall-clock time.
    # ---------------------------------------------------------------------
    clockTime <- time + T_FIRSTDOSE
    circ_parent <- 1 + amp_diurnal *
      cos(2 * 3.141592653589793 / 24 * (clockTime + shift_diurnal))
    circ_met <- 1 + amp_diurnal_met *
      cos(2 * 3.141592653589793 / 24 * (clockTime + shift_diurnal))

    # ---------------------------------------------------------------------
    # 2. Individual parameters.
    # ---------------------------------------------------------------------
    # Box-Cox transformed IIV on Vc (Willmann 2024 Data S1):
    # VC = VCpop * exp((exp(eta3)^lambda - 1) / lambda), and
    # exp(eta3)^lambda == exp(lambda * eta3), i.e. the Petersson 2009 form.
    eta_lvc_bc <- (exp(boxcox_lvc * etalvc) - 1) / boxcox_lvc

    cl <- exp(lcl + etalcl) * circ_parent
    vc <- exp(lvc + eta_lvc_bc)
    vp <- exp(lvp)
    q  <- exp(lq)

    ka <- exp(lka + etalka) *
      (1 + e_fed_highfat_ka * FED_HIGHFAT) *
      (1 + e_form_elz_hardgel_ka * FORM_ELZ_HARDGEL)
    tlag <- exp(ltlag + etaltlag) *
      (1 + e_fed_highfat_tlag * FED_HIGHFAT) *
      (1 + e_form_elz_hardgel_tlag * FORM_ELZ_HARDGEL)
    fdepot <- exp(lfdepot) * (1 + e_fed_highfat_fdepot * FED_HIGHFAT)

    # The three metabolite clearances share eta2; their Vmax terms carry no IIV.
    cl_m3034   <- exp(lcl_m3034 + etalcl_met) * circ_met
    vc_m3034   <- exp(lvc_m3034)
    vp_m3034   <- exp(lvp_m3034)
    q_m3034    <- exp(lq_m3034)
    vmax_m3034 <- exp(lvmax_m3034) * circ_met
    km_m3034   <- exp(lkm_m3034)

    cl_m27   <- exp(lcl_m27 + etalcl_met) * circ_met
    vc_m27   <- exp(lvc_m27)
    vp_m27   <- exp(lvp_m27)
    q_m27    <- exp(lq_m27)
    vmax_m27 <- exp(lvmax_m27) * circ_met
    km_m27   <- exp(lkm_m27)

    cl_m1821   <- exp(lcl_m1821 + etalcl_met) * circ_met
    vc_m1821   <- exp(lvc_m1821)
    vmax_m1821 <- exp(lvmax_m1821) * circ_met
    km_m1821   <- exp(lkm_m1821)

    ec50_nk1 <- exp(lec50_nk1)
    ec50_nk3 <- exp(lec50_nk3)

    # ---------------------------------------------------------------------
    # 3. Micro-constants (Willmann 2024 Data S1 K-notation in parentheses).
    # ---------------------------------------------------------------------
    kel   <- cl / vc          # total elinzanetant elimination, = K20 + K26 + K28
    k12   <- q / vc           # K24 and K35
    k21   <- q / vp           # K42 and K53
    k12_m3034 <- q_m3034 / vc_m3034   # K67
    k21_m3034 <- q_m3034 / vp_m3034   # K76
    k12_m27   <- q_m27 / vc_m27       # K89
    k21_m27   <- q_m27 / vp_m27       # K98
    kel_m3034 <- cl_m3034 / vc_m3034  # K610, M30/34 -> M18/21
    kel_m27   <- cl_m27 / vc_m27      # K810, M27 -> M18/21
    kel_m1821 <- cl_m1821 / vc_m1821  # K100

    # Concentrations driving the saturable elimination terms (ug/L).
    Cc       <- ugPerMg * central       / vc
    Cc_c13   <- ugPerMg * central_c13   / vc
    Cc_m3034 <- ugPerMg * central_m3034 / vc_m3034
    Cc_m27   <- ugPerMg * central_m27   / vc_m27
    Cc_m1821 <- ugPerMg * central_m1821 / vc_m1821

    # Saturable (Michaelis-Menten) elimination of each metabolite, converted
    # from ug/h to mg/h. Vmax [ug/L/h] * Vc [L] * C/(Km + C) gives ug/h.
    mm_m3034 <- vmax_m3034 * vc_m3034 * Cc_m3034 / (km_m3034 + Cc_m3034) / ugPerMg
    mm_m27   <- vmax_m27   * vc_m27   * Cc_m27   / (km_m27   + Cc_m27)   / ugPerMg
    mm_m1821 <- vmax_m1821 * vc_m1821 * Cc_m1821 / (km_m1821 + Cc_m1821) / ugPerMg

    # ---------------------------------------------------------------------
    # 4. ODE system (Willmann 2024 Figure 1). Compartment numbers from the
    #    source are given in comments.
    # ---------------------------------------------------------------------
    # Oral elinzanetant (1) -> elinzanetant central (2) <-> peripheral (4)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - k12 * central + k21 * peripheral1 -
      kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Intravenous [13C5]-elinzanetant (3) <-> peripheral (5). Shares CL, Q, Vc
    # and Vp with the oral model; its metabolites were not measured, so the
    # whole of CL leaves the system (K30).
    d/dt(central_c13) <- -k12 * central_c13 + k21 * peripheral1_c13 -
      kel * central_c13
    d/dt(peripheral1_c13) <- k12 * central_c13 - k21 * peripheral1_c13

    # M30/34 central (6) <-> peripheral (7); formed by K26 = FMET3034 * CL/Vc.
    d/dt(central_m3034) <- fm_m3034 * kel * central -
      k12_m3034 * central_m3034 + k21_m3034 * peripheral1_m3034 -
      kel_m3034 * central_m3034 - mm_m3034
    d/dt(peripheral1_m3034) <- k12_m3034 * central_m3034 -
      k21_m3034 * peripheral1_m3034

    # M27 central (8) <-> peripheral (9); formed by K28 = FMET27 * CL/Vc.
    d/dt(central_m27) <- fm_m27 * kel * central -
      k12_m27 * central_m27 + k21_m27 * peripheral1_m27 -
      kel_m27 * central_m27 - mm_m27
    d/dt(peripheral1_m27) <- k12_m27 * central_m27 - k21_m27 * peripheral1_m27

    # M18/21 central (10); formed from the linear clearances of both primary
    # metabolites (K610 and K810).
    d/dt(central_m1821) <- kel_m3034 * central_m3034 + kel_m27 * central_m27 -
      kel_m1821 * central_m1821 - mm_m1821

    # ---------------------------------------------------------------------
    # 5. Bioavailability and lag time on the oral depot.
    # ---------------------------------------------------------------------
    f(depot) <- fdepot
    alag(depot) <- tlag

    # ---------------------------------------------------------------------
    # 6. Receptor occupancy (Willmann 2024 Equation 2), driven by the
    #    elinzanetant plasma concentration only.
    # ---------------------------------------------------------------------
    roNK1 <- Cc / (Cc + ec50_nk1)
    roNK3 <- Cc / (Cc + ec50_nk3)

    # ---------------------------------------------------------------------
    # 7. Residual error: additive on the log scale in NONMEM == lnorm here.
    # ---------------------------------------------------------------------
    Cc       ~ lnorm(expSd)
    Cc_c13   ~ lnorm(expSd_c13)
    Cc_m3034 ~ lnorm(expSd_m3034)
    Cc_m27   ~ lnorm(expSd_m27)
    Cc_m1821 ~ lnorm(expSd_m1821)
  })
}
