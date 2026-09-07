Youssef_2024_bepirovirsen <- function() {
  description <- paste(
    "Three-compartment population PK model for bepirovirsen (a novel",
    "antisense oligonucleotide in development for chronic",
    "hepatitis B virus (HBV) infection) in healthy participants and",
    "participants with chronic HBV infection (Youssef 2024; N = 479",
    "subjects, 11021 plasma concentrations pooled from three studies -- a",
    "phase 1 single- and multiple-ascending-dose study in healthy",
    "volunteers, a phase 2a study and the phase 2b B-Clear study in chronic",
    "HBV infection). Disposition is linear and apparent throughout (CL/F,",
    "V/F): first-order subcutaneous absorption from a depot with an",
    "absorption lag time, distribution to a shallow and a deep peripheral",
    "compartment, and first-order elimination from the central compartment.",
    "Typical values for a 73 kg healthy participant are CL/F = 3.10 L/h,",
    "V2/F = 11.7 L, ka = 0.237 1/h, ALAG1 = 0.232 h, Q3/F = 0.107 L/h,",
    "V3/F = 33.2 L (shallow), Q4/F = 0.0428 L/h and V4/F = 63.7 L (deep).",
    "Two covariates were retained. Body weight acts on CL/F and V2/F as a",
    "power model referenced to 73 kg with ESTIMATED exponents of 0.494 and",
    "1.01 -- the clearance exponent is well below the theoretical 0.75.",
    "Chronic HBV infection acts as a fractional shift on the absorption lag",
    "time (-46.6%) and on the shallow peripheral volume (-33.7%), so the",
    "tabulated lag time and V3/F are the HEALTHY-participant values. Age,",
    "albumin, race, nucleos(t)ide-analogue treatment status and baseline",
    "HBsAg category were screened and not retained; despite the abstract's",
    "claim to the contrary, Asian versus non-Asian race is NOT in the final",
    "model. Inter-individual variability is estimated on all eight",
    "disposition parameters and is very large on the deep-compartment terms",
    "(146% on Q4/F, 198% on V4/F); residual error is proportional (21.2%).",
    sep = " "
  )
  reference <- paste(
    "Youssef AS, Ismail M, Han K, Magee M, Nader A. (2024).",
    "Population pharmacokinetics of bepirovirsen in healthy participants",
    "and participants with chronic hepatitis B virus infection: results",
    "from phase 1, 2a, and 2b studies.",
    "Infect Dis Ther 13:1515-1530.",
    "doi:10.1007/s40121-024-00980-9",
    sep = " "
  )
  vignette <- "Youssef_2024_bepirovirsen"
  units    <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Youssef 2024 Fig. 1 ("Schematic of the
  # starting base three-compartment PK model after SC administration") and the
  # Methods 'Base Structural Model Development' paragraph.
  #
  # The two peripheral compartments are recorded as `tissue` rather than
  # `plasma` because Youssef 2024 attributes the multi-exponential disposition
  # explicitly to tissue uptake: "bepirovirsen is expected to rapidly
  # distribute into tissues and to be eliminated via nuclease-mediated
  # metabolism in tissues" and "a rapid decline in plasma concentration and a
  # low level of renal elimination of full-length bepirovirsen over the first
  # 24-h post dose, reflecting the distribution in tissues" (Discussion). The
  # tissue is not identified -- the paper's stated limitation is "the absence
  # of measured liver concentrations" -- so these remain empirical mammillary
  # compartments and no organ is named here.
  compartmentData <- list(
    depot       = list(analyte = "bepirovirsen", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "bepirovirsen", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bepirovirsen", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "bepirovirsen", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline, time-fixed; one participant with a missing baseline",
        "weight was imputed with the population median (Youssef 2024",
        "Table 1 footnote a). Cohort median 70.2 kg, range 43.3-140.0 kg",
        "(Table 1). Applied as a power model to CL/F and V2/F only --",
        "weight on Q3/F and on V3/F was tested in forward selection and",
        "rejected (Supplementary Table 4, steps 1 and 3), and no weight",
        "effect was tested on the deep-compartment parameters",
        "(Supplementary Table 3).",
        "THE REFERENCE WEIGHT IS 73 kg, NOT the 70.2 kg cohort median.",
        "The Table 2 footnote prints the covariate model literally as",
        "'CL/F (l/h) = 3.10 x (WT/73)^0.494; V2/F (L) = 11.7 x",
        "(WT/73)^1.01', while Supplementary Equation 5 states only the",
        "general rule that continuous covariates were 'centered around",
        "their median values'. The printed 73 is used, and it is not a",
        "typo: for this linear model steady-state AUCtau = Dose/(CL/F)",
        "exactly, and against the three fixed-weight simulations of",
        "Table 4 (300 mg weekly at 40 / 70 / 100 kg, published median",
        "AUCtau 131.0 / 99.2 / 83.2 ug*h/mL) a 73 kg reference reproduces",
        "130.3 / 98.8 / 82.8 (-0.6% / -0.4% / -0.4%) whereas a 70 kg",
        "reference gives 127.6 / 96.8 / 81.1 (-2.6% / -2.4% / -2.5%) -- a",
        "systematic same-direction bias far outside the Monte Carlo error",
        "of a 1000-subject median. See the vignette Errata.",
        "Both exponents are ESTIMATED, not fixed (Table 2 rows 'WT on",
        "CL/F' and 'WT on V2/F', each with an RSE). Supplementary Table 2",
        "records that estimating them (run Poppk-base-v1.1) beat both no",
        "weight effect (v1) and the theoretical 0.75 / 1 pair (v1.2) on",
        "MVOF."
      ),
      source_name        = "WT"
    ),
    DIS_CHB = list(
      description        = "Chronic hepatitis B virus infection indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy participant; the complement group is the 21 healthy volunteers of phase 1 study 213725, 4.2 percent of the analysis population)",
      notes              = paste(
        "Per-subject, time-fixed. 96% of the analysis population had",
        "chronic HBV infection (Youssef 2024 Results 'Participants'); the",
        "healthy complement is study 213725 alone (N = 21 of 499,",
        "Table 1).",
        "Youssef 2024's structural reference subject is the HEALTHY",
        "participant -- Table 2 labels its rows 'ALAG1 (h), healthy' and",
        "'V3/F (L), healthy' and gives the chronic-HBV group a separate",
        "'Proportional shift' row -- so the indicator flags the patient",
        "group and the tabulated 0.232 h and 33.2 L are quoted unchanged.",
        "Both retained effects are fractional multipliers of the form",
        "(1 + theta * DIS_CHB): the Table 2 footnote prints 'ALAG1 (h) =",
        "0.232 * (1-0.466)' and 'V3/F (l/h) = 33.2 * (1-0.337)' for",
        "participants with chronic HBV infection.",
        "Disease status was tested on all eight disposition parameters",
        "(Supplementary Table 3) and survived on exactly two: V3/F entered",
        "at forward-selection step 1 (dMVOF -29.08, p = 6.93e-08) and",
        "ALAG1 at step 2 (dMVOF -15.20, p = 9.70e-05); both survived",
        "backward elimination (Supplementary Tables 4 and 5). Note the",
        "step-ordering artefact worth knowing when reusing these numbers:",
        "in step 1, before V3/F carried the effect, disease status on",
        "ALAG1 was the SECOND-WORST candidate in the table (dMVOF +0.74);",
        "it only became significant once the V3/F effect was in the model.",
        "Despite being statistically significant, neither effect is",
        "clinically relevant -- simulated steady-state exposures are",
        "essentially identical between healthy and chronic-HBV",
        "participants (Table 4: AUCtau 98.3 vs 98.2 ug*h/mL, Cmax 8.4 vs",
        "8.4 ug/mL)."
      ),
      source_name        = "Disease status (chronic HBV infection vs healthy)"
    )
  )

  # Covariates that Youssef 2024 screened but did NOT retain in the final
  # model. Documented here rather than in covariateData because no point
  # estimate exists for any of them and none is referenced in model().
  # Two tiers are pooled below: those carried into the formal stepwise search
  # (Supplementary Table 3) and those screened only graphically by
  # ETA-covariate correlation plots (Methods 'Assessment of Covariate
  # Effects'). Results: "After accounting for weight and disease status,
  # other demographic and baseline characteristics of clinical interest,
  # baseline HBsAg level, age, and albumin, were not identified as PK
  # covariates during formal covariate testing."
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Median 45.0 years, range 18-77 (Youssef 2024 Table 1); mean (SD)",
        "45.4 (11.6). Tested as a power model on CL/F and rejected at",
        "forward-selection step 3 (dMVOF -1.67, p = 0.196; Supplementary",
        "Tables 3 and 4)."
      )
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "271/499 (54.3%) Asian, 196 (39.3%) Caucasian, 31 (6.2%) Black,",
        "1 (0.2%) Native Indian or Alaska Native (Youssef 2024 Table 1).",
        "Tested as a proportional shift on CL/F and on V2/F and rejected",
        "at forward-selection step 3 (dMVOF +8.15 and +5.32, both",
        "p = 1; Supplementary Table 4). NOT in the final model, despite",
        "the abstract listing 'Asian versus non-Asian race' among the",
        "'key covariates included in the final model' -- see the vignette",
        "Errata. The paper nonetheless simulated Asian versus non-Asian",
        "exposures post hoc (Table 4) to check the weight-mediated",
        "difference, and found it below 20% and not clinically relevant."
      )
    ),
    ALBUMIN = list(
      description = "Baseline serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = paste(
        "Median 4.7 g/dL, range 3.8-5.5 (Youssef 2024 Table 1). Tested as",
        "a power model on CL/F and rejected at forward-selection step 3",
        "(dMVOF +4.76, p = 1; Supplementary Tables 3 and 4)."
      )
    ),
    CONMED_NA_HBV = list(
      description = "Concomitant nucleos(t)ide analogue (NA) treatment indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "237/499 (47.5%) on stable NA therapy, 241 (48.3%) not, 21 (4.2%)",
        "not applicable (the healthy cohort) -- Youssef 2024 Table 1.",
        "Tested as a proportional shift on CL/F and on V2/F and rejected",
        "at forward-selection step 3 (dMVOF +7.77 and -0.48, p = 1 and",
        "0.489; Supplementary Table 4). The Discussion concludes 'NA",
        "treatment does not affect bepirovirsen PK'. This name is",
        "provisional and is NOT a registered canonical -- it is recorded",
        "only as documentation of the covariate screen and is never",
        "referenced in model()."
      )
    ),
    HBSAG_BL_LOG10 = list(
      description = "Baseline hepatitis B surface antigen",
      units       = "log10 IU/mL",
      type        = "continuous",
      notes       = paste(
        "Median 3.48 log10 IU/mL, range 1.64-5.72; missing for the 21",
        "healthy participants (Youssef 2024 Table 1). Entered the formal",
        "search DICHOTOMISED, not as the continuous value: Supplementary",
        "Table 3 lists 'Baseline HBsAg category (>1000 IU/mL vs <=1000",
        "IU/mL)' as a categorical covariate on CL/F, while the",
        "Supplementary Table 4 footnote c instead describes the split as",
        "'>3 log10 IU/mL vs <=3 log10 IU/mL' (numerically the same",
        "threshold, 1000 IU/mL = 3 log10 IU/mL). Rejected at",
        "forward-selection step 3 as the single worst candidate in the",
        "table (dMVOF +15.16, p = 1)."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "318/499 (63.7%) male, hence 36.3% female (Youssef 2024 Table 1).",
        "Screened only graphically -- Methods 'Assessment of Covariate",
        "Effects' lists sex among the categorical covariates 'evaluated by",
        "visual inspection for differences between groups' -- and it does",
        "not appear in the formal stepwise search of Supplementary",
        "Table 3."
      )
    ),
    CREAT = list(
      description = "Baseline serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste(
        "Median 0.8 mg/dL, range 0.3-1.3 (Youssef 2024 Table 1). Screened",
        "graphically by ETA-covariate linear regression and scatter plots",
        "(Methods) and not carried into the formal search."
      )
    ),
    CRCL = list(
      description = "Baseline creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Median 120 mL/min, range 51.7-270 (Youssef 2024 Table 1);",
        "estimated GFR median 106 mL/min/1.73m2, range 64.1-220.",
        "Screened graphically and not carried into the formal search. The",
        "Discussion notes that participants with mild renal impairment",
        "were included and that renal function was not identified as a",
        "covariate, but that moderate-to-severe renal disease was an",
        "exclusion criterion."
      )
    ),
    BILI = list(
      description = "Baseline total bilirubin",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste(
        "Median 0.5 mg/dL, range 0.2-2.0 (Youssef 2024 Table 1, the",
        "unlabelled row following albumin). Screened graphically and not",
        "carried into the formal search."
      )
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Median 24 U/L, range 7-349 (Youssef 2024 Table 1, the unlabelled",
        "row preceding HBsAg). Screened graphically and not carried into",
        "the formal search."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 479,
    n_studies      = 3,
    age_range      = "18-77 years",
    age_median     = "45 years",
    weight_range   = "43.3-140.0 kg",
    weight_median  = "70.2 kg",
    sex_female_pct = 36.3,
    race_ethnicity = c(Asian = 54.3, Caucasian = 39.3, Black = 6.2, NativeIndianOrAlaskaNative = 0.2),
    disease_state  = "chronic hepatitis B virus infection (96%, both on and off stable nucleos(t)ide-analogue therapy) pooled with healthy volunteers (4%); participants with cirrhosis or liver failure, or with moderate-to-severe renal disease, were excluded",
    dose_range     = "75-450 mg subcutaneous: single ascending dose, six doses over 3 weeks, or 150-300 mg weekly for up to 24 weeks",
    regions        = "China Mainland, East Asia, Japan and Other (the four strata of Youssef 2024 Supplementary Table 7)",
    notes          = paste(
      "Youssef 2024 Table 1 (demographics by study and overall).",
      "Three study cohorts: 213725 / NCT03020745 (phase 1, healthy,",
      "N = 21 analysed), 205695 / NCT02981602 (phase 2a, chronic HBV,",
      "N = 23) and 209668 B-Clear / NCT04449029 (phase 2b, chronic HBV,",
      "N = 455). Note the two subject counts that appear in the paper and",
      "why n_subjects is 479 rather than 499: 516 were randomised and 499",
      "were eligible for the analysis (the N of Table 1), but 11021 of",
      "12140 concentrations (90.8%) from 479 of 499 participants (96.0%)",
      "were actually used to develop the model (Results 'Participants').",
      "525 of the used samples (4.8%) were below the 1 ng/mL LLOQ and",
      "were handled by the Beal M3 method; outliers identified",
      "graphically and by CWRES were excluded from the fit and",
      "reintroduced only as a sensitivity analysis, which left the fixed",
      "effects essentially unchanged but raised every IIV except that on",
      "ka by >= 42%."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # Structural parameters. Youssef 2024 Table 2, 'Estimate' column. Every
    # clearance and volume is APPARENT (CL/F, V/F): bepirovirsen was given
    # only subcutaneously in these three studies, so no absolute
    # bioavailability is identifiable and none is estimated. Accordingly no
    # f(depot) term appears in model() -- F is folded into the tabulated
    # values, exactly as the paper's own CL/F and V2/F notation states.
    #
    # Compartment mapping. The paper numbers compartments in NONMEM style
    # with the depot as 1, so its V2/V3/V4 and Q3/Q4 map onto the nlmixr2lib
    # canonical names as:
    #   V2/F = 11.7 L                                  -> vc   (central)
    #   V3/F = 33.2 L  with Q3/F = 0.107  L/h  (shallow) -> vp,  q   (peripheral1)
    #   V4/F = 63.7 L  with Q4/F = 0.0428 L/h  (deep)    -> vp2, q2  (peripheral2)
    # peripheral1 is the SHALLOW (fast) compartment and peripheral2 the DEEP
    # (slow) one, following the paper's own numbering.
    #
    # The lag-time and shallow-volume values below are the HEALTHY-participant
    # values; Table 2 labels those two rows explicitly 'ALAG1 (h), healthy'
    # and 'V3/F (L), healthy'. The chronic-HBV shift is applied in model().
    # -----------------------------------------------------------------------
    lcl  <- log(3.10)   ; label("Apparent systemic clearance CL/F, 73 kg reference (L/h)")                                  # Youssef 2024 Table 2 'CL/F (l/h)' = 3.10 (RSE 2.89%)
    lvc  <- log(11.7)   ; label("Apparent central volume of distribution V2/F, 73 kg reference (L)")                        # Youssef 2024 Table 2 'V2/F (L)' = 11.7 (RSE 2.47%)
    lka  <- log(0.237)  ; label("First-order subcutaneous absorption rate constant (1/h)")                                  # Youssef 2024 Table 2 'KA (h-1)' = 0.237 (RSE 3.48%)
    ltlag <- log(0.232) ; label("Absorption lag time ALAG1, healthy participant (h)")                                       # Youssef 2024 Table 2 'ALAG1 (h), healthy' = 0.232 (RSE 4.68%)
    lq   <- log(0.107)  ; label("Apparent intercompartmental clearance to the shallow peripheral compartment Q3/F (L/h)")   # Youssef 2024 Table 2 'Q3/F (l/h)' = 0.107 (RSE 2.91%)
    lvp  <- log(33.2)   ; label("Apparent shallow peripheral volume of distribution V3/F, healthy participant (L)")         # Youssef 2024 Table 2 'V3/F (L), healthy' = 33.2 (RSE 1.98%)
    lq2  <- log(0.0428) ; label("Apparent intercompartmental clearance to the deep peripheral compartment Q4/F (L/h)")      # Youssef 2024 Table 2 'Q4/F (l/h)' = 0.0428 (RSE 2.69%)
    lvp2 <- log(63.7)   ; label("Apparent deep peripheral volume of distribution V4/F (L)")                                 # Youssef 2024 Table 2 'V4/F (L)' = 63.7 (RSE 2.25%)

    # -----------------------------------------------------------------------
    # Allometry. Both exponents are ESTIMATED, not fixed at theoretical
    # values: Table 2 gives them their own rows with RSEs, and Supplementary
    # Table 2 records that estimating them (MVOF 50023.25) beat both omitting
    # the weight effect (50172.45) and fixing them at 0.75 / 1 (50062.22).
    # The reference weight of 73 kg comes from the Table 2 footnote's printed
    # equations and is verified against the paper's own Table 4 simulations --
    # see covariateData$WT$notes and the vignette Errata.
    #
    # Do not confuse these with the 0.477 / 0.92 pair in Supplementary
    # Table 2 or the "0.48 and 0.92" quoted in the main-text Methods: those
    # belong to the BASE model, before the two disease-status effects entered.
    # -----------------------------------------------------------------------
    e_wt_cl <- 0.494 ; label("Power exponent on (WT/73) for CL/F (unitless)")   # Youssef 2024 Table 2 'WT on CL/F' = 0.494 (RSE 0.0237); Table 2 footnote 'CL/F (l/h) = 3.10 x (WT/73)^0.494'
    e_wt_vc <- 1.01  ; label("Power exponent on (WT/73) for V2/F (unitless)")   # Youssef 2024 Table 2 'WT on V2/F' = 1.01 (RSE 0.0247); Table 2 footnote 'V2/F (L) = 11.7 x (WT/73)^1.01'

    # -----------------------------------------------------------------------
    # Covariate effects of chronic HBV infection. Table 2 calls both rows a
    # 'Proportional shift', and the Table 2 footnote prints the two resulting
    # equations literally as 'ALAG1 (h) = 0.232 * (1-0.466)' and 'V3/F (l/h)
    # = 33.2 * (1-0.337)' for participants with chronic HBV infection. Each
    # coefficient is therefore used as a (1 + theta * DIS_CHB) multiplier, so
    # a chronic-HBV participant gets a 0.534x lag time and a 0.663x shallow
    # peripheral volume relative to a healthy participant.
    # -----------------------------------------------------------------------
    e_dis_chb_tlag <- -0.466 ; label("Effect of chronic HBV infection on the absorption lag time (fraction)")     # Youssef 2024 Table 2 'Chronic HBV infection on ALAG1' = -0.466 (RSE 8.40%)
    e_dis_chb_vp   <- -0.337 ; label("Effect of chronic HBV infection on the shallow peripheral volume (fraction)") # Youssef 2024 Table 2 'Chronic HBV infection on V3/F' = -0.337 (RSE 2.56%)

    # -----------------------------------------------------------------------
    # Inter-individual variability. Youssef 2024 Table 2 column
    # 'IIV (CV%)', which prints each cell as "<estimate> (<CV%>)".
    #
    # The printed estimate is the VARIANCE omega^2, not an SD, and is entered
    # here unchanged. This is not an assumption -- it is forced by the
    # paper's own two columns. Supplementary Methods Equation 2 defines the
    # reported CV% for a log-normal parameter, i.e. CV = sqrt(exp(w2) - 1),
    # and back-transforming all eight estimates through it reproduces the
    # printed CV% column to within 0.4 percentage points:
    #                        KA   ALAG1   CL/F   V2/F   Q3/F   V3/F    Q4/F    V4/F
    #   sqrt(exp(w2)-1)    46.0    38.2   31.7   40.8   51.2   75.2   145.8   197.6
    #   printed CV%        45.9    38.2   31.7   40.8   51.2   75.2   146.0   198.0
    # Reading the same numbers as SDs instead would give 19.4 / 13.7 / 9.6 /
    # 15.5 / 23.6 / 47.1 / 163.3 / 339.6, which matches nothing.
    #
    # Table 2 reports no off-diagonal covariance terms, so all eight etas are
    # independent. Shrinkage is substantial on several (54.6% on V2/F, 54.5%
    # on ALAG1, 54.1% on Q3/F, 47.9% on KA) and the sensitivity analysis in
    # Results shows every IIV except that on KA rises by >= 42% when the
    # excluded outliers are reintroduced -- so the deep-compartment CVs of
    # 146% and 198% in particular should be treated as fitted descriptions of
    # a sparse tail, not as transferable population variability.
    # -----------------------------------------------------------------------
    etalcl   ~ 0.096  # Youssef 2024 Table 2 'CL/F (l/h)' IIV = 0.096 (CV 31.7%, shrinkage 31.9%)
    etalvc   ~ 0.154  # Youssef 2024 Table 2 'V2/F (L)' IIV = 0.154 (CV 40.8%, shrinkage 54.6%)
    etalka   ~ 0.192  # Youssef 2024 Table 2 'KA (h-1)' IIV = 0.192 (CV 45.9%, shrinkage 47.9%)
    etaltlag ~ 0.136  # Youssef 2024 Table 2 'ALAG1 (h), healthy' IIV = 0.136 (CV 38.2%, shrinkage 54.5%)
    etalq    ~ 0.233  # Youssef 2024 Table 2 'Q3/F (l/h)' IIV = 0.233 (CV 51.2%, shrinkage 54.1%)
    etalvp   ~ 0.448  # Youssef 2024 Table 2 'V3/F (L), healthy' IIV = 0.448 (CV 75.2%, shrinkage 41.8%)
    etalq2   ~ 1.140  # Youssef 2024 Table 2 'Q4/F (l/h)' IIV = 1.140 (CV 146.0%, shrinkage 25.8%)
    etalvp2  ~ 1.590  # Youssef 2024 Table 2 'V4/F (L)' IIV = 1.590 (CV 198.0%, shrinkage 15.4%)

    # -----------------------------------------------------------------------
    # Residual error. PROPORTIONAL ONLY. Methods: "Residual variability (RV)
    # for plasma bepirovirsen concentration-time data was estimated using a
    # proportional error model". A combined additive-plus-CCV model was
    # tested and rejected -- Supplementary Table 2 run Poppk-base-v2 reached
    # an MVOF of 6098546.10 with the comment "Model was not stable, resulting
    # in unreliable model fits (all concentrations close to 0)".
    #
    # Like the IIV column, the printed 0.0449 is a variance: Table 2 prints
    # "0.0449 (21.2)" and sqrt(0.0449) = 0.21190 -> 21.2% CV, so propSd is
    # the square root.
    # -----------------------------------------------------------------------
    propSd <- 0.21190 ; label("Proportional residual error (fraction)")   # Youssef 2024 Table 2 'RV (CV%)' = 0.0449 (21.2) (RSE 1.59%); sqrt(0.0449) = 0.21190
  })

  model({
    # ---------------------------------------------------------------------
    # 1. Individual parameters.
    #
    # Weight enters CL/F and V2/F as the Table 2 footnote's power model
    # referenced to 73 kg; chronic HBV infection enters the lag time and the
    # shallow peripheral volume as that footnote's fractional shifts. At
    # WT = 73 kg and DIS_CHB = 0 every multiplier collapses to 1 and all
    # eight parameters reduce exactly to their Table 2 typical values.
    #
    # No weight effect is applied to Q3/F, V3/F, Q4/F or V4/F: weight on
    # Q3/F and V3/F was tested and rejected (Supplementary Table 4) and
    # weight on the deep-compartment parameters was never tested
    # (Supplementary Table 3). Inventing an allometric term for them would
    # not be in the paper.
    # ---------------------------------------------------------------------
    cl   <- exp(lcl + etalcl) * (WT / 73)^e_wt_cl
    vc   <- exp(lvc + etalvc) * (WT / 73)^e_wt_vc
    ka   <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag) * (1 + e_dis_chb_tlag * DIS_CHB)
    q    <- exp(lq + etalq)
    vp   <- exp(lvp + etalvp) * (1 + e_dis_chb_vp * DIS_CHB)
    q2   <- exp(lq2 + etalq2)
    vp2  <- exp(lvp2 + etalvp2)

    # 2. Micro-constants for the three-compartment system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ---------------------------------------------------------------------
    # 3. ODE system (Youssef 2024 Fig. 1). A subcutaneous dose enters the
    #    depot, is held for the absorption lag, then transfers first-order
    #    at ka into the central compartment, which exchanges with a shallow
    #    and a deep peripheral compartment and eliminates first-order.
    #
    #    Dose records should be placed in `depot`. Because the lag is
    #    implemented with alag(depot), a dose record's time is the
    #    ADMINISTRATION time -- do not pre-offset it by tlag or the lag is
    #    applied twice.
    # ---------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    alag(depot) <- tlag

    # ---------------------------------------------------------------------
    # 4. Observation. Dose is in mg and vc in L, so central / vc is mg/L,
    #    i.e. ug/mL -- the unit Youssef 2024 uses for Cmax and Ctrough in
    #    Tables 3 and 4 (AUCtau in ug*h/mL). No scaling factor is needed.
    #    Note that the assay LLOQ of 1 ng/mL and Supplementary Tables 6-7
    #    are stated in ng/mL, i.e. 1000x this scale.
    # ---------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
