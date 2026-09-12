FernandezTeruel_2025_capivasertib <- function() {
  description <- "Three-compartment population PK model for capivasertib (oral pan-AKT inhibitor) with parallel first-order and zero-order absorption, absorption lag time, and sigmoidal time- and dose-dependent auto-inhibition of apparent clearance, with power effects of body weight and age on CL0/F, in patients with advanced solid tumours and HR-positive/HER2-negative advanced breast cancer receiving capivasertib plus fulvestrant (Fernandez Teruel 2025)"
  reference <- "Fernandez Teruel C, Cullberg M, Gonzalez-Garcia I, Schiavon G, Zhang L, Zhou D. Population pharmacokinetics and exposure-response analyses for capivasertib in combination with fulvestrant in patients with breast cancer. Clin Transl Sci. 2025;18:e70286. doi:10.1111/cts.70286"
  vignette <- "FernandezTeruel_2025_capivasertib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(
      analyte = "capivasertib", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "capivasertib", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "capivasertib", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral2 = list(
      analyte = "capivasertib", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column BBW (baseline body weight). Normalised to 67 kg in both",
        "of the paper's printed covariate equations. Note that 67 kg is the",
        "NORMALISATION CONSTANT carried over from the parent model",
        "(FernandezTeruel_2024_capivasertib.R), not this analysis's cohort",
        "median, which is 65 kg (Table S1); the equations printed on p. 7 read",
        "(BBW/67) in both places, and the equation governs. Enters twice:",
        "a POWER effect on the initial apparent clearance,",
        "CL0/F * (WT/67)^e_wt_cl, and a POWER effect on the logit-scale",
        "first-order absorbed fraction, LogitF1 * (WT/67)^e_wt_logitffo.",
        "Both forms CHANGED from the 2024 parent model, where weight entered",
        "CL0/F as a linear deviation (1 + (WT - 67) * 0.00585); reading the",
        "2025 coefficient 0.302 as a per-kg linear term would be absurd.",
        "Baseline-only (time-fixed); the analysis carried a single per-patient",
        "weight."
      ),
      source_name        = "BBW"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column AGE. New in this analysis: age was NOT a covariate in",
        "the 2024 parent model. Enters CL0/F as a power term normalised to",
        "57 years, (AGE/57)^e_age_cl, with e_age_cl negative, so apparent",
        "clearance FALLS with increasing age. 57 years is both the printed",
        "normalisation constant in the p. 7 equation and the cohort median",
        "(Table S1, range 26-87). The paper reports the resulting exposure",
        "difference between the 5th and 95th percentiles of the covariate",
        "distribution as < 20%, and concludes it is not clinically relevant."
      ),
      source_name        = "AGE"
    ),
    FASTED_STRICT = list(
      description        = "Overnight-fasted dosing indicator (1 = dose taken after an overnight fast; 0 = semi-fasted or fed)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (semi-fasted or fed)",
      notes              = paste(
        "Source column FASTED. Gates the absorption lag time: the whole ALAG1",
        "term is multiplied by (1 - FASTED_STRICT), so there is no lag under an",
        "overnight fast. Carried unchanged from the 2024 parent model, whose",
        "OAK arm supplied the food-effect data; the equation on p. 7 retains",
        "the (1 - FASTED) factor. The reference population used for this",
        "paper's covariate forest plot (Figure 2) is semi-fasted, i.e.",
        "FASTED_STRICT = 0, which is the state in which the lag is active.",
        "Per dose record, not per subject."
      ),
      source_name        = "FASTED"
    ),
    FORM_CAPSULE = list(
      description        = "Capsule formulation indicator (1 = capsule; 0 = tablet)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tablet)",
      notes              = paste(
        "Source column CAP. The only retained formulation effect is on the",
        "absorption lag time (0.468 h capsule vs 0.287 h tablet, Table 2), both",
        "under non-overnight-fasted conditions; the formulation acts on neither",
        "bioavailability nor Ka. CAPItello-291, the Phase III trial that",
        "motivates this update, dosed the tablet exclusively",
        "(FORM_CAPSULE = 0). Per dose record: the OAK study compared the two",
        "formulations within patient."
      ),
      source_name        = "CAP"
    ),
    DOSE_CAPIVASERTIB_MG = list(
      description        = "Planned capivasertib dose level per administration",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column DOSE. The planned twice-daily dose level, 80-800 mg",
        "(Table S1), not the daily total. Enters the maximal auto-inhibition of",
        "CL/F as a linear-deviation term centred on 480 mg,",
        "(1 + (DOSE_CAPIVASERTIB_MG - 480) * e_dose_cl_time_max), applied to",
        "the LOG-scale magnitude. Because e_dose_cl_time_max is negative and",
        "the log-scale magnitude is itself negative, the extent of",
        "auto-inhibition GROWS with dose. This encoding reproduces both of the",
        "paper's printed reductions exactly: 11.2% at 400 mg and 29.1% at",
        "640 mg (Sect. 3.2 states 11.2% and 29.2%). Carried as a drug-specific",
        "column rather than a bare `DOSE` because rxode2's etTrans consumes a",
        "column literally named DOSE and never exposes it to model()."
      ),
      source_name        = "DOSE"
    ),
    CONMED_PACLITAXEL = list(
      description        = "Concomitant paclitaxel indicator (1 = capivasertib given with paclitaxel)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (capivasertib monotherapy, or with fulvestrant)",
      notes              = paste(
        "Source column PACL. Paclitaxel was co-administered in 90 of 851",
        "patients (10.6%, Table S1), in the BEECH study and part B of the",
        "China PK study. The term (1 + PACL * Imax_pacl) is RETAINED in the",
        "paper's printed Imax equation (p. 7) and Imax_pacl is defined in the",
        "equation legend, but ITS ESTIMATE IS REPORTED NOWHERE IN THIS PAPER",
        "-- not in Table 2, not in the supplement, not in any figure panel.",
        "The coefficient used here is therefore BORROWED from the 2024 parent",
        "model (FernandezTeruel_2024_capivasertib.R, Imax_pacl = 1.15), held",
        "constant, and is the only parameter in this file not taken from",
        "Fernandez Teruel 2025. Setting CONMED_PACLITAXEL = 1 therefore",
        "simulates a cross-fit approximation, NOT a 2025 estimate; every other",
        "parameter shared by the two fits was re-estimated. See the",
        "`e_pacl_cl_time_max` entry in ini() and the vignette Errata for the",
        "full provenance. The reference level pools capivasertib monotherapy",
        "and capivasertib + fulvestrant, since fulvestrant was screened and",
        "not retained."
      ),
      source_name        = "PACL"
    )
  )

  # Covariates that Fernandez Teruel 2025 screened in the PsN stepwise
  # covariate-model build (forward p = 0.005, backward p = 0.001, Sect. 2.3) and
  # did NOT retain in the final model. Documentation only -- the paper reports
  # no point estimate for any of them, so none can be encoded, and none is
  # referenced in model(). Sources: Sect. 2.3 (the screened list), Sect. 3.2 and
  # Figure 2b (the not-significant verdict), Table S1 (the frequencies).
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator", units = "(binary)", type = "binary",
      notes = "Screened; not significant (Figure 2b, 'gender'). 88.8% female, Table S1."
    ),
    CONMED_FULVESTRANT = list(
      description = "Concomitant fulvestrant indicator", units = "(binary)", type = "binary",
      notes = paste(
        "Screened explicitly as the lead covariate of interest for this update",
        "(Sect. 2.3) and NOT retained (Figure 2b). 468 of 851 patients (55.0%,",
        "Table S1). This is a notable null result: the paper's whole purpose is",
        "the capivasertib + fulvestrant combination, and it finds fulvestrant",
        "does not alter capivasertib PK."
      )
    ),
    HEPATIC_IMPAIR = list(
      description = "NCI-ODWG hepatic function category", units = "(category)", type = "categorical",
      notes = "Screened; not significant. Normal 64.0%, mild 35.0%, moderate 0.8%, no severe, Table S1."
    ),
    RENAL_IMPAIR = list(
      description = "Renal function category", units = "(category)", type = "categorical",
      notes = "Screened; not significant. Normal 54.4%, mild 34.3%, moderate 10.6%, no severe, Table S1."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator", units = "(binary)", type = "binary",
      notes = paste(
        "Screened as part of the race/region covariate; not significant.",
        "27.8% of patients, Table S1. The paper additionally compared Chinese",
        "vs non-Chinese and Asian vs non-Asian populations explicitly",
        "(Sect. 2.3) and found steady-state exposure in Chinese/Asian patients",
        "only < 15% higher than the rest of the world (Sect. 3.2)."
      )
    ),
    RACE_BLACK = list(
      description = "Black race indicator", units = "(binary)", type = "binary",
      notes = paste(
        "Screened as part of the race covariate; not significant. 12 of 851",
        "patients (1.4%), Table S1. The Discussion cautions that the number of",
        "Black patients was limited."
      )
    ),
    REGION_CHINA = list(
      description = "Mainland China / Taiwan region indicator", units = "(binary)", type = "binary",
      notes = "Screened; not significant. 11.0% of patients, Table S1."
    ),
    SMOKING = list(
      description = "Current smoking status indicator", units = "(binary)", type = "binary",
      notes = "Screened; not significant (Figure 2b). Never 42.5%, current 4.5%, former 19.7%, missing 33.3%, Table S1."
    ),
    CONMED_ARA = list(
      description = "Concomitant acid-reducing agent indicator", units = "(binary)", type = "binary",
      notes = "Screened; not significant (Figure 2b, 'ARA'). 216 of 851 patients (25.4%), Table S1."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 851L,
    n_studies      = 6L,
    n_observations = 5960L,
    age_range      = "26-87 years",
    age_median     = "57 years",
    weight_range   = "32-150 kg",
    weight_median  = "65 kg",
    sex_female_pct = 88.8,
    race_ethnicity = c(
      White = 60.9, Black = 1.4, Asian = 27.8,
      `American Indian or Alaska Native` = 1.8,
      `Native Hawaiian or Other Pacific Islander` = 0.1,
      Other = 7.5, Missing = 0.5
    ),
    disease_state  = "Advanced solid malignancies, and HR-positive / HER2-negative locally advanced or metastatic breast cancer resistant to aromatase inhibitors",
    dose_range     = paste(
      "80-800 mg orally twice daily, most commonly 400 mg (63.8%) or 480 mg",
      "(22.4%); continuous dosing (8.0%) or one of two intermittent schedules,",
      "4 days on / 3 days off (4/3, 86.4%) or 2 days on / 5 days off",
      "(2/5, 5.6%). The Phase III regimen is 400 mg twice daily [4/3] with",
      "fulvestrant 500 mg."
    ),
    renal_function = "Normal 54.4%, mild impairment 34.3%, moderate impairment 10.6%, no severe impairment",
    hepatic_function = "Normal 64.0%, mild impairment 35.0%, moderate impairment 0.8%, no severe impairment",
    co_medication  = "Fulvestrant 55.0%, paclitaxel 10.6%, acid-reducing agent 25.4%",
    regions        = "Global; China 11.0%, Asia excluding China 15.3%, rest of world 73.7%",
    notes          = paste(
      "Pooled from six Phase I-III studies: Study 1, BEECH, Study 4",
      "(all-Japanese), OAK, a China PK study, and the Phase III CAPItello-291",
      "trial. This analysis UPDATES the 441-patient, four-study model of",
      "FernandezTeruel_2024_capivasertib.R by adding CAPItello-291 and the",
      "China PK study (Sect. 2.3). 6630 concentrations were collected and 5960",
      "(89.9%) were analysed; 670 were excluded, mainly 436 (6.6%) drawn before",
      "treatment start, plus 220 (3.3%) below the 1.00 ng/mL limit of",
      "quantification and 9 (0.1%) unexpectedly high. Demographics are",
      "Table S1. The exposure-response analyses used subsets of this cohort:",
      "798 patients for efficacy (Table S2) and 468 for safety (Table S3)."
    )
  )

  ini({
    # --- Structural disposition -------------------------------------------
    # All parameters are apparent (/F): capivasertib was dosed orally only in
    # these six studies, so F is not separately identifiable and is absorbed
    # into CL, V and Q. (Absolute bioavailability is 29% from a separate
    # clinical pharmacology study cited in Sect. 1, but that study is not in
    # this dataset and F is NOT estimated here.)
    #
    # Cross-check on the disposition block: V2/F + V3/F + V4/F
    # = 49 + 113 + 103 = 265 L, exactly the "apparent volume of distribution at
    # steady state was 265 L" reported in Sect. 3.2.
    lcl  <- log(59.6); label("Initial apparent clearance CL0/F, before auto-inhibition (L/h)")            # Fernandez Teruel 2025 Table 2: CL0/F = 59.6 L/h (bootstrap median 59.3, 95% CI 56.8-62)
    lvc  <- log(49);   label("Apparent central volume V2/F (L)")                                          # Fernandez Teruel 2025 Table 2: V2/F = 49 L (bootstrap median 49.7, 95% CI 36.7-68.9)
    lvp  <- log(113);  label("Apparent first peripheral volume V3/F (L)")                                 # Fernandez Teruel 2025 Table 2: V3/F = 113 L (bootstrap median 116, 95% CI 95.7-137)
    lq   <- log(2.52); label("Apparent intercompartmental clearance to peripheral1, Q3/F (L/h)")          # Fernandez Teruel 2025 Table 2: Q3/F = 2.52 L/h (bootstrap median 2.36, 95% CI 1.43-3.38)
    lvp2 <- log(103);  label("Apparent second peripheral volume V4/F (L)")                                # Fernandez Teruel 2025 Table 2: V4/F = 103 L (bootstrap median 98.7, 95% CI 74-159)
    lq2  <- log(23.4); label("Apparent intercompartmental clearance to peripheral2, Q4/F (L/h)")          # Fernandez Teruel 2025 Table 2: Q4/F = 23.4 L/h (bootstrap median 23, 95% CI 16.5-31.6)

    # --- Absorption --------------------------------------------------------
    # Two parallel input routes from the same dose: a fraction F1 enters a
    # first-order depot with rate Ka and lag time ALAG1, and the remaining
    # F2 = 1 - F1 is delivered as a zero-order input straight into the central
    # compartment over a duration D2 (Fernandez Teruel 2025 p. 7 equations for
    # F1i and F2i, and Sect. 3.2).
    lka <- log(0.441); label("First-order absorption rate constant Ka (1/h)")                             # Fernandez Teruel 2025 Table 2: Ka = 0.441 1/h (bootstrap median 0.441, 95% CI 0.358-0.538)
    ld1 <- log(45);    label("Duration of the zero-order input into central, D2 (h)")                     # Fernandez Teruel 2025 Table 2: D2 = 45 h (bootstrap median 45, 95% CI 43.6-45.6)
    # The paper estimates the first-order fraction on the LOGIT scale (p. 7
    # equation for F1i), so `logitffo` -- not `lffo` -- is the correct
    # encoding, and the between-subject random effect is additive on the logit
    # scale as printed. expit(1.6) = 0.832, i.e. about 83% of the dose is
    # absorbed by the first-order route in a 67 kg patient.
    logitffo <- 1.6;   label("Logit of the fraction of the dose absorbed by the first-order route, F1 (unitless)")  # Fernandez Teruel 2025 Table 2: LogitF1 = 1.6 (bootstrap median 1.47, 95% CI 1.16-2.18)

    # Two absorption lag times, selected by formulation. Both apply only when
    # the dose was NOT taken after an overnight fast (p. 7 equation for
    # ALAG1i, which carries the (1 - FASTED) factor).
    ltlag_tab <- log(0.287); label("Absorption lag time, tablet, non-overnight-fasted (h)")               # Fernandez Teruel 2025 Table 2: Lag1_tab = 0.287 h (bootstrap median 0.287, 95% CI 0.223-0.334)
    ltlag_cap <- log(0.468); label("Absorption lag time, capsule, non-overnight-fasted (h)")              # Fernandez Teruel 2025 Table 2: Lag1_cap = 0.468 h (bootstrap median 0.461, 95% CI 0.405-0.505)

    # --- Time- and dose-dependent auto-inhibition of CL/F ------------------
    # Fernandez Teruel 2025 p. 7:
    #   CL/Fi = CL0/Fi * (1 - exp(Imax_i) * Time^5 / (T50^5 + Time^5))
    # The paper's `Imax` is therefore the LOG of the maximal fractional
    # inhibition, not the fraction itself. That reading is confirmed
    # arithmetically by the paper's own two printed reductions (Sect. 3.2):
    # with Imax = -1.87 and Imax_dose = -0.00213,
    #   400 mg: exp(-1.87 * (1 + (400 - 480) * -0.00213)) = 0.1120 -> 11.2%
    #   640 mg: exp(-1.87 * (1 + (640 - 480) * -0.00213)) = 0.2914 -> 29.1%
    # against the reported 11.2% and 29.2%. The value is therefore stored here
    # on the log scale rather than back-transformed.
    lcl_time_max  <- -1.87;      label("log of the maximal fractional inhibition of CL/F at t >> T50 (unitless); sign applied in model()")  # Fernandez Teruel 2025 Table 2: Imax = -1.87 (bootstrap median -1.81, 95% CI -2.53 to -1.44)
    lcl_t50       <- log(126);   label("T50, time at which half of the maximal inhibition of CL/F is reached (h)")                          # Fernandez Teruel 2025 Table 2: T50 = 126 h (bootstrap median 127, 95% CI 63.3-144)
    lcl_time_hill <- fixed(log(5)); label("Hill shape coefficient of the time-on-CL/F function (unitless)")                                 # Fernandez Teruel 2025 p. 7: the exponent 5 appears literally in the CL/Fi equation (Time^5, T50^5); fixed, as in the 2024 parent model

    # --- Covariate effects -------------------------------------------------
    # Fernandez Teruel 2025 p. 7:
    #   CL0/Fi = CL0/F * (BBW/67)^CL0_BBW * (AGE/57)^CL0_AGE * exp(eta)
    # BOTH are POWER exponents. This is a change of functional form from the
    # 2024 parent model, where weight entered CL0/F as a linear deviation.
    e_wt_cl            <- 0.302;    label("Power exponent of (WT / 67 kg) acting on CL0/F (unitless)")                              # Fernandez Teruel 2025 Table 2: CL0_BBW = 0.302 (bootstrap median 0.334, 95% CI 0.189-0.492)
    e_age_cl           <- -0.314;   label("Power exponent of (AGE / 57 years) acting on CL0/F (unitless)")                          # Fernandez Teruel 2025 Table 2: CL0_AGE = -0.314 (bootstrap median -0.319, 95% CI -0.476 to -0.151)
    e_wt_logitffo      <- -1.27;    label("Power exponent of (WT / 67 kg) acting on the logit-scale first-order fraction (unitless)")  # Fernandez Teruel 2025 Table 2: F1_BBW = -1.27 (bootstrap median -1.2, 95% CI -1.62 to -0.811)
    e_dose_cl_time_max <- -0.00213; label("Linear-deviation coefficient of planned dose on log-Imax, per mg above 480 mg (1/mg)")   # Fernandez Teruel 2025 Table 2: Imax_dose = -0.00213 (bootstrap median -0.00209, 95% CI -0.00253 to -0.00133)

    # NOT REPORTED IN THIS PAPER -- the value below is BORROWED from the 2024
    # parent model, held constant, and is the ONE parameter in this file that
    # does not come from Fernandez Teruel 2025. Read the provenance before
    # using the paclitaxel arm.
    #
    # The paper's printed Imax equation (p. 7) retains the factor
    # (1 + PACL * Imax_pacl), and the equation legend states "Imax_pacl
    # represents the relationship between concomitant paclitaxel and Imax",
    # but no estimate for it appears anywhere on disk: not in Table 2 (which
    # tabulates every OTHER coefficient in every equation), not in the
    # supplement (Tables S1-S5, Figures S1-S6), and not inside any figure
    # panel (the publisher's native-resolution figure files were checked).
    #
    # Carrying the 2024 estimate across is a cross-fit borrow, not a 2025
    # estimate: every other shared parameter WAS re-estimated between the two
    # fits (Imax -1.54 -> -1.87, Imax_dose -0.00183 -> -0.00213), so the true
    # 2025 value is unknown and is not necessarily 1.15. What the borrow does
    # buy is a usable paclitaxel arm instead of an inert one, and the equation
    # form is identical between the two papers, so the value is at least
    # dimensionally and structurally transportable.
    #
    # The borrowed value reproduces the 2024 paper's own printed covariate
    # result, which is the check that the ENCODING (as opposed to the
    # transported value) is right. Fernandez-Teruel 2024 Sect. 4 reports
    # patients on concomitant paclitaxel had 20% higher CL_ss/F (median ratio
    # 1.20). On the 2024 fit at 400 mg:
    #   no paclitaxel: exp(-1.54 * 1.1464)          = 0.1711 -> 17.1% inhibition
    #   + paclitaxel : exp(-1.54 * 1.1464 * 2.15)   = 0.0225 ->  2.3% inhibition
    #   CL_ss/F ratio = (1 - 0.0225) / (1 - 0.1711) = 1.179   vs printed 1.20
    # (1.179 is the typical-subject value; the printed 1.20 is a median over a
    # cohort carrying IIV on Imax.)
    e_pacl_cl_time_max <- fixed(1.15); label("Fractional change in log-Imax with concomitant paclitaxel (unitless)")  # NOT from Fernandez Teruel 2025, which reports no estimate: borrowed from the 2024 parent, Fernandez-Teruel 2024 Table 3 Imax_pacl = 1.15 (RSE 12.5%, bootstrap 95% CI 1-6.91), per operator decision 2026-09-11. See vignette Errata.

    # --- Between-subject variability ---------------------------------------
    # Table 2 reports every IIV in a single column headed "CV (%)", and one of
    # those five etas -- LogitF1 -- is ADDITIVE on the logit scale (it appears
    # as `+ eta_LogitF1` inside the p. 7 logit expression). A log-normal
    # 100 * sqrt(exp(omega^2) - 1) back-transform is undefined for an additive
    # random effect, so the column can only be reporting 100 * omega. That one
    # convention is therefore applied uniformly: omega = CV / 100 for every eta
    # below. This is the same reading, on the same table layout by the same
    # first author, as FernandezTeruel_2024_capivasertib.R.
    #
    # Only five etas are reported. The 2024 parent model additionally held the
    # Ka and D2 IIVs at a fixed 15% CV; the 2025 Table 2 lists no IIV for
    # either, so none is encoded here.
    etalcl          ~ 0.382^2;    # Fernandez Teruel 2025 Table 2: CL0/F IIV CV = 38.2% (bootstrap median 37.5, 95% CI 32.4-42.4)
    etalvc          ~ 1.03^2;     # Fernandez Teruel 2025 Table 2: V2/F IIV CV = 103% (bootstrap median 103, 95% CI 75.7-143)
    etalcl_time_max ~ 0.953^2;    # Fernandez Teruel 2025 Table 2: Imax IIV CV = 95.3% (bootstrap median 89.1, 95% CI 54.5-132)
    etaltlag        ~ 0.837^2;    # Fernandez Teruel 2025 Table 2: Lag1 IIV CV = 83.7% (bootstrap median 73.2, 95% CI 60.2-95.5)
    etalogitffo     ~ 0.713^2;    # Fernandez Teruel 2025 Table 2: LogitF1 IIV CV = 71.3% (bootstrap median 66.4, 95% CI 51.2-78.3); additive on the logit scale per the p. 7 equation

    # --- Residual unexplained variability ----------------------------------
    # Combined proportional + additive, rendered here on the linear
    # concentration scale, which is nlmixr2's `prop() + add()` form.
    propSd <- 0.43;   label("Proportional residual error SD (fraction)")        # Fernandez Teruel 2025 Table 2: RUV proportional CV = 43% (bootstrap median 44.9, 95% CI 40.4-49.2)
    addSd  <- 1.12;   label("Additive residual error SD (ng/mL)")               # Fernandez Teruel 2025 Table 2: RUV additive = 1.12 ug/L (bootstrap median 0.928, 95% CI 0.516-1.58); 1 ug/L = 1 ng/mL
  })

  model({
    # --- 1. Absorption -----------------------------------------------------
    ka <- exp(lka)
    d1 <- exp(ld1)

    # Fraction absorbed by the first-order route, on the logit scale, with a
    # power effect of body weight and an additive between-subject random
    # effect. Fernandez Teruel 2025 p. 7:
    #   F1i = exp(LogitF1 * (BBW/67)^F1_BBW + eta) /
    #         (1 + exp(LogitF1 * (BBW/67)^F1_BBW + eta))
    # Note the exponent acts on the NORMALISED weight and multiplies the whole
    # logit, exactly as printed.
    logitffo_ind <- logitffo * (WT / 67)^e_wt_logitffo + etalogitffo
    ffo <- expit(logitffo_ind)

    # Absorption lag time. Formulation selects the magnitude and the overnight
    # fast switches it off entirely. Fernandez Teruel 2025 p. 7:
    #   ALAG1i = (Lag1_tab*(1-CAP) + Lag1_cap*CAP) * (1 - FASTED) * exp(eta)
    tlag <- (exp(ltlag_tab) * (1 - FORM_CAPSULE) + exp(ltlag_cap) * FORM_CAPSULE) *
      (1 - FASTED_STRICT) * exp(etaltlag)

    # --- 2. Disposition ----------------------------------------------------
    # Fernandez Teruel 2025 p. 7:
    #   CL0/Fi = CL0/F * (BBW/67)^CL0_BBW * (AGE/57)^CL0_AGE * exp(eta)
    cl_base <- exp(lcl + etalcl) * (WT / 67)^e_wt_cl * (AGE / 57)^e_age_cl
    vc      <- exp(lvc + etalvc)
    vp      <- exp(lvp)
    q       <- exp(lq)
    vp2     <- exp(lvp2)
    q2      <- exp(lq2)

    # --- 3. Time-dependent auto-inhibition of CL/F -------------------------
    # Fernandez Teruel 2025 p. 7: the dose and paclitaxel effects are
    # MULTIPLICATIVE on the log-scale magnitude, and so is the random effect:
    #   Imax_i = Imax * (1 + (DOSE - 480)*Imax_dose) * (1 + PACL*Imax_pacl) * exp(eta)
    # Imax is the LOG of the fractional inhibition and is negative, so the
    # multiplier acts in the COUNTER-INTUITIVE direction: a multiplier > 1
    # drives Imax_i further negative, shrinking exp(Imax_i), which SHALLOWS
    # the inhibition and raises CL/F; a multiplier < 1 DEEPENS it. The paper's
    # own two printed reductions confirm the direction -- 640 mg carries the
    # SMALLER multiplier (0.659) and gives the LARGER inhibition (29.1%),
    # while 400 mg carries 1.170 and gives 11.2%.
    lcl_time_max_ind <- lcl_time_max *
      (1 + (DOSE_CAPIVASERTIB_MG - 480) * e_dose_cl_time_max) *
      (1 + CONMED_PACLITAXEL * e_pacl_cl_time_max) *
      exp(etalcl_time_max)
    cl_time_max_i <- -exp(lcl_time_max_ind)
    cl_t50        <- exp(lcl_t50)
    cl_time_hill  <- exp(lcl_time_hill)

    # Fernandez Teruel 2025 p. 7. `t` is time since the first dose, matching
    # the paper's `Time`.
    cl <- cl_base *
      (1 + cl_time_max_i * t^cl_time_hill / (cl_t50^cl_time_hill + t^cl_time_hill))

    # --- 4. Micro-constants and ODEs ---------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <- ka * depot - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # --- 5. Parallel input routing -----------------------------------------
    # The same dose is recorded twice: a fraction `ffo` enters `depot` as a
    # lagged first-order input, and the complementary fraction 1 - ffo
    # (Fernandez Teruel 2025 p. 7, F2i = 1 - F1i) enters `central` as a
    # zero-order input of duration d1. The central dose record must request a
    # modelled duration (rate = -2 in the event table).
    f(depot)     <- ffo
    alag(depot)  <- tlag
    f(central)   <- 1 - ffo
    dur(central) <- d1

    # --- 6. Observation ----------------------------------------------------
    # Dose is in mg and vc in L, so central/vc is mg/L; x 1000 gives ng/mL
    # (= ug/L), the unit Fernandez Teruel 2025 reports concentrations in.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
