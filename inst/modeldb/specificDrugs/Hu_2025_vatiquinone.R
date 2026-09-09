Hu_2025_vatiquinone <- function() {
  description <- paste(
    "Two-compartment population pharmacokinetic model for vatiquinone",
    "(PTC743), a first-in-class 15-lipoxygenase inhibitor developed for",
    "Friedreich's ataxia and other mitochondrial diseases, with parallel",
    "zero-order and first-order oral absorption and linear elimination",
    "(Hu 2025; 343 participants and 4,608 quantifiable plasma samples",
    "pooled from eight phase I/II/III studies in adult healthy volunteers",
    "and adult and pediatric patients with Friedreich's ataxia or other",
    "mitochondrial diseases). 74.4% of the absorbed dose enters through a",
    "first-order arm delayed by a 2.79 h lag time; the remaining 25.6%",
    "enters the central compartment as a 6.03 h zero-order input, a dual",
    "pathway the authors attribute to the compound's extreme lipophilicity",
    "(cLogP 7.8) and partial lymphatic uptake via chylomicrons. Vatiquinone",
    "exposure is dominated by prandial state: relative to the reference",
    "medium-fat meal, a liquid PediaSure supplement gives 6.9% and the",
    "fasted state 3.6% of the reference exposure, i.e. a medium-fat meal",
    "raises exposure roughly 14-fold and 28-fold respectively. Strong",
    "CYP3A4 modulation moves apparent clearance in both directions",
    "(itraconazole to 23.5%, rifampicin to 202% of the monotherapy value).",
    "Patients with Friedreich's ataxia carry a 50.1% lower relative",
    "bioavailability and a 40.6% lower apparent clearance, which combine to",
    "a net 19% lower steady-state AUC. Apparent clearance also scales with",
    "body weight (power exponent 0.915, reference 65 kg) and inversely with",
    "body mass index (exponent -0.975, reference 21.6 kg/m^2), and the",
    "central volume scales linearly with body weight."
  )
  reference <- paste(
    "Hu Y, Gao L, Lee L, Cherry JJ, Kong R. Characterizing Population",
    "Pharmacokinetics of Vatiquinone in Healthy Volunteers and Patients",
    "with Friedreich's Ataxia. Pharmaceuticals. 2025;18(9):1339.",
    "doi:10.3390/ph18091339"
  )
  vignette <- "Hu_2025_vatiquinone"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: checked against Hu 2025 Figure 5
  # (model schematic) and the Figure 2 VPC concentration axis (ng/mL).
  compartmentData <- list(
    depot       = list(analyte = "vatiquinone", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "vatiquinone", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vatiquinone", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Two separate power effects, both referenced to 65 kg, which Hu 2025",
        "Section 2.3 states explicitly ('the reference body weight value is",
        "set to 65 kg in the population') -- a rounded standard value, not the",
        "cohort median, which is 58.6 kg (Supplementary Table S2). On apparent",
        "clearance the exponent is theta14 = 0.915 (Table 2, estimated). On the",
        "apparent central volume the exponent is 1 and is written as a literal",
        "in the Section 2.3 V equation with no corresponding theta in Table 2,",
        "so it is fixed rather than estimated. Cohort range 6.30-119 kg",
        "(Supplementary Table S2), a 19-fold span driven by the pediatric",
        "Friedreich's-ataxia and mitochondrial-disease studies."
      ),
      source_name        = "BWT"
    ),
    BMI = list(
      description        = "Baseline body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on apparent clearance, exponent theta17 = -0.975 (Hu 2025",
        "Table 2), referenced to 21.6 kg/m^2, which Section 2.3 identifies as",
        "the population median and Supplementary Table S2 confirms as the",
        "overall median. The negative exponent means clearance FALLS as BMI",
        "rises, so exposure rises with BMI -- consistent with a highly",
        "lipophilic compound partitioning into adipose tissue. BMI and WT act",
        "on clearance simultaneously and in opposite directions; they are",
        "correlated in the cohort, so the two exponents are not independently",
        "interpretable as marginal effects. Cohort range 11.9-41.4 kg/m^2."
      ),
      source_name        = "BMI"
    ),
    DIS_FRDA = list(
      description        = "Friedreich's ataxia disease indicator: 1 = patient with Friedreich's ataxia, 0 = healthy volunteer or patient with another mitochondrial disease.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Time-fixed per subject. Hu 2025 Section 2.3 defines the indicator",
        "as 'FAi is set to 1 for participants with FA and is set to 0 for",
        "healthy volunteers'. The pooled dataset also contains 54 patients",
        "with other mitochondrial diseases (Supplementary Table S1); the",
        "Discussion reports that 'the effect of other mitochondrial diseases",
        "on vatiquinone exposures was examined, no significant differences",
        "were observed in PK exposures compared to healthy volunteers', so",
        "that stratum shares the healthy-volunteer reference level and is",
        "coded DIS_FRDA = 0. Two linear-deviation effects, both from Table 2:",
        "theta16 = -0.501 on relative bioavailability and theta15 = -0.406 on",
        "apparent clearance. Cohort composition 173/343 (50.4%) FA, 116",
        "(33.8%) healthy volunteers, 54 (15.7%) other mitochondrial disease."
      ),
      source_name        = "FA"
    ),
    FED = list(
      description        = "Fed-vs-fasted dose-record indicator: 1 = dose taken with food, 0 = fasted",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Used COMPLEMENTED in model(): the paper's own indicator is FST",
        "('FSTi is set to 1 in the fasted state', Hu 2025 Section 2.3), so",
        "the canonical column maps as FST = 1 - FED and the effect enters as",
        "exp(e_fasted_fdepot * (1 - FED)). Only the 6 participants in the",
        "EPI743-12-001 food-effect crossover contributed fasted records",
        "(Supplementary Table S1: 6/343 fasted, 6 liquid, 331 medium fat), so",
        "the fasted arm is a deliberate crossover challenge rather than a",
        "population stratum. Both FED = 1 levels (the reference medium-fat",
        "meal and the liquid PediaSure supplement) are distinguished by the",
        "companion FED_LIQUIDSUPP indicator; the model's prandial reference",
        "is the medium-fat meal (FED = 1, FED_LIQUIDSUPP = 0), NOT the fasted",
        "state, which is the usual orientation for this family and is the",
        "reason the fasted effect is carried on the complement."
      ),
      source_name        = "FST"
    ),
    FED_LIQUIDSUPP = list(
      description        = "Liquid nutritional-supplement meal at dosing indicator: 1 = dose taken with a liquid nutritional supplement (PediaSure), 0 = any other prandial state.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Hu 2025 Section 2.3: 'LQDi assumes a value of 1 when a liquid",
        "PediaSure meal is ingested'. Multiplicative log-scale effect on",
        "relative bioavailability, theta12 = -2.671 (Table 2), i.e. exposure",
        "falls to exp(-2.671) = 6.9% of the medium-fat-meal reference. The",
        "protocol reference meal is defined in the Table 5 footnote as a",
        "medium-fat meal, 'approximately 25% to 50% fat'; PediaSure is a",
        "ready-to-drink pediatric nutritional supplement supplying far less",
        "lipid than that solid meal. Because vatiquinone is extremely",
        "lipophilic (cLogP 7.8, Discussion) and is thought to reach the",
        "systemic circulation partly by chylomicron-mediated lymphatic",
        "uptake, the available lipid load rather than the fat PERCENTAGE",
        "drives bioavailability -- which is why this indicator is not a",
        "member of the fat-fraction axis (FED_LOWFAT / FED_HIGHFAT). Only 6",
        "of 343 participants contributed liquid-supplement records, all from",
        "the EPI743-12-001 three-way food-effect crossover."
      ),
      source_name        = "LQD"
    ),
    CONMED_ITRACONAZOLE = list(
      description        = "Concomitant itraconazole (strong CYP3A4 inhibitor) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Hu 2025 Section 2.3: 'ITRi is set to 1 when there is a concomitant",
        "administration of itraconazole'. Multiplicative log-scale effect on",
        "apparent clearance, theta10 = -1.446 (Table 2): CL/F falls to",
        "exp(-1.446) = 23.5% of the monotherapy value, a 4.25-fold increase",
        "in steady-state AUC. Sourced from the dedicated drug-drug",
        "interaction study EPI743-18-002 (49 healthy adults, 400 mg single",
        "dose, medium-fat meals, Table 5). Time-varying by dose record in",
        "that crossover design."
      ),
      source_name        = "ITR"
    ),
    CONMED_RIFAMPICIN = list(
      description        = "Concomitant rifampicin (rifampin; strong CYP3A4 inducer) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Hu 2025 Section 2.3: 'RFMi equals 1 when rifampin is",
        "co-administered'. Multiplicative log-scale effect on apparent",
        "clearance, theta11 = 0.704 (Table 2): CL/F rises to exp(0.704) =",
        "202% of the monotherapy value, halving steady-state AUC. The source",
        "paper uses the USAN name 'rifampin'; the canonical column uses the",
        "INN 'rifampicin'. Same crossover study as CONMED_ITRACONAZOLE",
        "(EPI743-18-002), so the two indicators are mutually exclusive per",
        "dose record and both are 0 under vatiquinone monotherapy."
      ),
      source_name        = "RFM"
    )
  )

  # Screened by Hu 2025 but NOT retained in the final model. The Conclusions
  # state that "race, sex, age, pediatric maturation, and various hepatic and
  # renal function biomarkers had no apparent impact on vatiquinone
  # pharmacokinetic exposures"; the Discussion adds that maturation and age on
  # systemic clearance were assessed specifically because 29 participants were
  # under 7 years old, and that the other-mitochondrial-disease stratum showed
  # no exposure difference from healthy volunteers. Baseline distributions for
  # every entry below are in Supplementary Tables S1 (categorical) and S2
  # (continuous). Documented here rather than in covariateData because a
  # covariateData entry that model() never references is a convention warning.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age", units = "years", type = "continuous",
      notes = "Median 20 years, range 1-67 (Supplementary Table S2). Screened on CL/F together with a pediatric maturation function; neither reached the forward-addition criterion."
    ),
    SEXF = list(
      description = "Female sex indicator", units = "(binary)", type = "binary",
      notes = "178/343 (51.9%) female (Supplementary Table S1). Screened, not retained."
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator", units = "(binary)", type = "binary",
      notes = "35/343 (10.2%); White 279 (81.3%), Asian 12 (3.5%), missing 17 (5.0%) (Supplementary Table S1). Screened, not retained."
    ),
    ALB = list(
      description = "Baseline serum albumin", units = "g/L", type = "continuous",
      notes = "Median 46 g/L, range 37-54 (Supplementary Table S2; reported there under a 'g/dL' column header whose values are unambiguously g/L). Screened as a hepatic marker, not retained."
    ),
    BILI = list(
      description = "Baseline total bilirubin", units = "umol/L", type = "continuous",
      notes = "Median 6.84, range 2.00-25.7 (Supplementary Table S2; reported there under a 'mg/dL' column header whose values are unambiguously umol/L). Screened as a hepatic marker, not retained."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase", units = "IU/L", type = "continuous",
      notes = "Median 20, range 9-89 (Supplementary Table S2). Screened, not retained."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase", units = "IU/L", type = "continuous",
      notes = "Median 17, range 5-127 (Supplementary Table S2). Screened, not retained."
    ),
    ALP = list(
      description = "Baseline alkaline phosphatase", units = "IU/L", type = "continuous",
      notes = "Median 85, range 26-352 (Supplementary Table S2). Screened, not retained."
    ),
    CREAT = list(
      description = "Baseline serum creatinine", units = "umol/L", type = "continuous",
      notes = "Median 55, range 18-115 (Supplementary Table S2; reported there under a 'mg/dL' column header whose values are unambiguously umol/L). Screened as a renal marker, not retained."
    ),
    CRCL = list(
      description = "Baseline creatinine clearance", units = "mL/min", type = "continuous",
      notes = "Median 123, range 47-472 (Supplementary Table S2). Computed by Hu 2025 as (140-AGE)*WEIGHT*1.23/CREAT, multiplied by 0.85 for females (Supplementary Table S2 footnote). Screened as a renal marker, not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 343L,
    n_studies      = 8L,
    n_samples      = 4608L,
    age_range      = "1-67 years (median 20); 144/343 (42.0%) pediatric (<18 years), including 29 participants under 7 years of age",
    weight_range   = "6.30-119 kg (median 58.6)",
    bmi_range      = "11.9-41.4 kg/m^2 (median 21.6)",
    sex_female_pct = 51.9,
    race_ethnicity = "White 279 (81.3%), Black or African American 35 (10.2%), Asian 12 (3.5%), not reported 17 (5.0%); Hispanic or Latino ethnicity 76 (22.2%).",
    disease_state  = paste(
      "116 (33.8%) adult healthy volunteers, 173 (50.4%) patients with",
      "Friedreich's ataxia (adult and pediatric), and 54 (15.7%) pediatric",
      "patients with other mitochondrial diseases (epilepsy indication)."
    ),
    dose_range     = paste(
      "120-1,400 mg (median 400 mg), given as single doses or three times",
      "daily. Capsule 284 (82.8%) and oral solution 59 (17.2%). Pediatric",
      "studies dosed 15 mg/kg below 13 kg body weight and 200 mg at or",
      "above 13 kg."
    ),
    regions        = "Not reported in the article.",
    notes          = paste(
      "Pooled analysis of eight phase I/II/III studies (NONMEM 7.5):",
      "EPI743-12-001 (three-way food-effect crossover: fasted, liquid",
      "PediaSure, medium-fat meal), EPI743-18-002 (itraconazole /",
      "rifampin drug-drug interaction crossover), PTC743-NEU-004-FA part 1,",
      "PTC743-CNS-006-HV (thorough-QT run-in, 400 and 1,400 mg),",
      "EPI-2010-006, PTC743-NEU-003-FA (NCT04577352),",
      "PTC743-NEU-005-FA (NCT05485987, children under 7) and",
      "PTC743-MIT-001-EP (NCT04378075). Every study other than",
      "EPI743-12-001 dosed exclusively with a medium-fat meal, so the",
      "prandial covariates are identified almost entirely by the 18-subject",
      "food-effect crossover. Samples below the limit of quantification were",
      "discarded (M1 method). No observation had |CWRES| > 5."
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # STRUCTURAL DISPOSITION  --  Hu 2025 Table 2 (final model), apparent
    # (oral, /F) parameters throughout. The paper never fitted intravenous
    # data, so absolute bioavailability is not identifiable and every volume
    # and clearance below is an apparent value.
    # ------------------------------------------------------------------------
    lcl <- log(162.721) ;   label("Apparent clearance CL/F at the reference covariate values (L/h)")
    # Hu 2025 Table 2 theta 6: CL/F = 162.721 L/h (ASE 10.14, 90% CI 145.990; 179.452)
    lvc <- log(180.748) ;   label("Apparent central volume of distribution V/F at 65 kg (L)")
    # Hu 2025 Table 2 theta 5: V/F = 180.748 L (ASE 21.38, 90% CI 145.471; 216.025)
    lq  <- log(67.896) ;    label("Apparent intercompartmental clearance Q/F (L/h)")
    # Hu 2025 Table 2 theta 8: Q/F = 67.896 L/h (ASE 6.136, 90% CI 57.772; 78.019)
    lvp <- log(4852.69) ;   label("Apparent peripheral volume of distribution V2/F (L)")
    # Hu 2025 Table 2 theta 7: V2/F = 4852.69 L (ASE 773.404, 90% CI 3576.573; 6128.807)

    # ------------------------------------------------------------------------
    # PARALLEL ZERO-ORDER + FIRST-ORDER ORAL ABSORPTION  --  Hu 2025 Figure 5
    # and Table 2. The dose is split between a first-order arm (fraction FK0,
    # rate Ka, lag TLAG1) and a zero-order arm delivered straight into the
    # central compartment over TK0 hours. In the paper's NONMEM ADVAN
    # numbering the first-order arm is compartment 1 (hence TLAG1) and the
    # central compartment is compartment 2 (hence the Table 2 footnote
    # "TK0: absorption duration, D2").
    #
    # RE-PARAMETERISATION -- IMPORTANT, see the vignette Assumptions and
    # deviations section. Hu 2025 gives a SINGLE parameter FK0 two distinct
    # jobs in the Section 2.3 equation
    #     FK0i = TVFK0 * exp(LQDi*theta12) * exp(FSTi*theta13) * (1 + FAi*theta16)
    # namely (a) the 74.4% / 25.6% split between the two absorption arms and
    # (b) the carrier of the prandial and disease relative-bioavailability
    # effects. Those two jobs are separated here into `logitffo` (the split,
    # estimated) and `lfdepot` (relative bioavailability, anchored at 1 and
    # carrying the covariate effects), because a covariate that only moved the
    # SPLIT could not change total exposure at all -- both arms deliver to the
    # same central compartment, so the fractions would still sum to 1.
    #
    # The paper's own Table 3 settles which reading is operative. It reports
    # the FK0-mediated exposure ratios as Cmax,ss / Cmin,ss / AUC0-24h,ss =
    # 0.50 / 0.50 / 0.50 for Friedreich's ataxia, 0.07 / 0.07 / 0.07 for the
    # liquid PediaSure meal and 0.04 / 0.04 / 0.04 for the fasted state. Three
    # different exposure metrics moving by an IDENTICAL factor is the
    # signature of a pure scalar multiplier on the whole concentration-time
    # curve, which happens only when the multiplier scales BOTH absorption
    # arms and leaves their ratio untouched. Those three published ratios are
    # reproduced to two decimals by the encoding below and are re-derived as
    # regression tests in the validation vignette.
    # ------------------------------------------------------------------------
    lka       <- log(0.200) ;    label("First-order absorption rate constant Ka (1/h)")
    # Hu 2025 Table 2 theta 2: Ka = 0.200 1/h (ASE 0.016, 90% CI 0.173; 0.228)
    ltlag1    <- log(2.787) ;    label("First-order absorption lag time TLAG1 (h)")
    # Hu 2025 Table 2 theta 4: TLAG1 = 2.787 h (ASE 0.044, 90% CI 2.715; 2.859)
    ld0       <- log(6.034) ;    label("Zero-order absorption duration TK0 / D2 (h)")
    # Hu 2025 Table 2 theta 3: TK0 = 6.034 h (ASE 0.102, 90% CI 5.866; 6.203)
    logitffo  <- logit(0.744) ;  label("Logit fraction of the absorbed dose entering the first-order arm (FK0)")
    # Hu 2025 Table 2 theta 1: FK0 = 0.744 (ASE 0.021, 90% CI 0.710; 0.779);
    # the Discussion states the complementary 25.6% goes through the
    # zero-order arm, so the two arms sum to the full relative dose.
    lfdepot   <- fixed(log(1)) ; label("Relative bioavailability at the reference prandial and disease state (unitless, log-scale)")
    # Anchor, not a paper-reported theta. Hu 2025 has no absolute-F data, so
    # relative bioavailability is defined as 1 at the model's reference
    # condition (medium-fat meal, no Friedreich's ataxia); every prandial and
    # disease effect below is expressed relative to that anchor.

    # ------------------------------------------------------------------------
    # COVARIATE EFFECTS  --  Hu 2025 Table 2 thetas 10-17, in the exact
    # functional forms printed in the Section 2.3 equations:
    #
    #   FK0i = TVFK0 * e^(LQDi*t12) * e^(FSTi*t13) * (1 + FAi*t16)
    #   CLi  = TVCL * (1 + FAi*t15) * e^(ITRi*t10) * e^(RFMi*t11)
    #                 * (BWTi/65)^t14 * (BMIi/21.6)^t17 * exp(eta_CL)
    #   Vi   = TVV * (BWTi/65)^1 * exp(eta_V)
    #
    # Note the deliberate mix of forms: the two disease effects are LINEAR
    # deviations (1 + FA*theta), while the prandial and comedication effects
    # are LOG-LINEAR (exp(x*theta)). That is what the paper prints, and it is
    # what makes theta15 = -0.406 read directly as "40.6% lower clearance"
    # and theta16 = -0.501 as "50.1% lower relative bioavailability", both of
    # which the Discussion quotes back as "a 50% reduction in relative
    # bioavailability (FK0) and a 40% decrease in clearance".
    # ------------------------------------------------------------------------
    e_wt_cl                   <- 0.915 ;    label("Power exponent of (WT / 65 kg) on CL/F (unitless)")
    # Hu 2025 Table 2 theta 14: BWT on CL/F = 0.915 (ASE 0.123, 90% CI 0.712; 1.118)
    e_bmi_cl                  <- -0.975 ;   label("Power exponent of (BMI / 21.6 kg/m^2) on CL/F (unitless)")
    # Hu 2025 Table 2 theta 17: BMI on CL/F = -0.975 (ASE 0.249, 90% CI -1.386; -0.564)
    e_wt_vc                   <- fixed(1) ; label("Power exponent of (WT / 65 kg) on V/F (unitless)")
    # Hu 2025 Section 2.3 V equation: the exponent is printed as the literal 1
    # and has no corresponding row in Table 2, so it is fixed, not estimated.
    e_dis_frda_cl             <- -0.406 ;   label("Linear-deviation effect of Friedreich's ataxia on CL/F (unitless)")
    # Hu 2025 Table 2 theta 15: Disease FA on CL/F = -0.406 (ASE 0.082, 90% CI -0.541; -0.272)
    e_dis_frda_fdepot         <- -0.501 ;   label("Linear-deviation effect of Friedreich's ataxia on relative bioavailability (unitless)")
    # Hu 2025 Table 2 theta 16: Disease FA on FK0 = -0.501 (ASE 0.057, 90% CI -0.595; -0.408)
    e_conmed_itraconazole_cl  <- -1.446 ;   label("Log-scale effect of concomitant itraconazole on CL/F (unitless)")
    # Hu 2025 Table 2 theta 10: Itraconazole on CL/F = -1.446 (ASE 0.175, 90% CI -1.735; -1.158)
    e_conmed_rifampicin_cl    <- 0.704 ;    label("Log-scale effect of concomitant rifampicin on CL/F (unitless)")
    # Hu 2025 Table 2 theta 11: Rifampin on CL/F = 0.704 (ASE 0.099, 90% CI 0.541; 0.867)
    e_fed_liquidsupp_fdepot   <- -2.671 ;   label("Log-scale effect of a liquid PediaSure meal on relative bioavailability (unitless)")
    # Hu 2025 Table 2 theta 12: Liquid PediaSure on FK0 = -2.671 (ASE 0.130, 90% CI -2.886; -2.457)
    e_fasted_fdepot           <- -3.324 ;   label("Log-scale effect of the fasted state on relative bioavailability (unitless)")
    # Hu 2025 Table 2 theta 13: Fasted statuses on FK0 = -3.324 (ASE 0.167, 90% CI -3.599; -3.050).
    # Applied to (1 - FED) because the canonical column is oriented fed = 1
    # while the paper's FST indicator is oriented fasted = 1.

    # ------------------------------------------------------------------------
    # INTER-INDIVIDUAL VARIABILITY  --  Hu 2025 Table 2 IIV block.
    #
    # Only ONE random effect is reported: IIV-CL/F = 0.191, with the
    # accompanying "(%CV)" column reading 45.880 and shrinkage 24.370. That
    # column confirms the estimate is an OMEGA VARIANCE on the log scale
    # rather than an SD or a CV, because sqrt(exp(0.191) - 1) = 0.4588 = the
    # printed 45.880%. (Had 0.191 been an SD the implied CV would be 19.5%;
    # had it been a CV the variance would be 0.0365. Neither reproduces the
    # printed column.)
    #
    # The Section 2.3 equations additionally carry eta terms on Ka and on V,
    # but Table 2 reports no estimate for either and the text never mentions
    # them again, so their magnitudes are simply unreported. They are OMITTED
    # here rather than written as `~ fixed(0)`, because a zero-variance
    # diagonal makes OMEGA singular and breaks the Cholesky sampler used by
    # rxSolve (same reasoning and precedent as Thoueille_2026_salmeterol.R).
    # The gap is recorded in the vignette Assumptions and deviations section.
    # No IIV is reported on FK0, TK0, TLAG1, Q/F or V2/F either.
    # ------------------------------------------------------------------------
    etalcl ~ 0.191
    # Hu 2025 Table 2 IIV-CL/F: omega^2 = 0.191 (ASE 0.025); printed %CV 45.880 = sqrt(exp(0.191) - 1) * 100

    # ------------------------------------------------------------------------
    # RESIDUAL UNEXPLAINED VARIABILITY  --  Hu 2025 Table 2, "9 Additive
    # residual" = 1.062 (ASE 0.011, %RSE 1.036, 90% CI 1.043; 1.081).
    #
    # The paper's Methods never write out an $ERROR block, so the scale has to
    # be read off the reported number. Three facts settle it as an additive
    # residual on LOG-TRANSFORMED concentration, i.e. log-normal residual
    # error on the linear scale:
    #   1. The row is numbered 9, inside the THETA sequence (between Q/F = 8
    #      and Itraconazole on CL/F = 10) and is reported with a %RSE and a
    #      90% CI in the same style as the thetas -- the NONMEM idiom
    #      W = THETA(9); Y = LOG(IPRED) + W*EPS(1) with EPS(1) fixed to 1,
    #      which makes THETA(9) the residual SD directly.
    #   2. An additive residual of 1.062 ng/mL on the LINEAR scale would be
    #      ~0.05% of the observed median concentration and would leave the
    #      model with effectively no residual error at all, which cannot be
    #      reconciled with the Figure 2 VPC.
    #   3. The Figure 2 log-scale VPC spans roughly 10-50 ng/mL at the 5th
    #      percentile against 2,000-5,000 ng/mL at the 95th. A ~100-fold
    #      5th-to-95th spread implies a total log-scale SD near
    #      ln(100)/(2*1.645) = 1.40; removing the CL random effect
    #      (sqrt(0.191) = 0.437) leaves ~1.33 for the residual, which brackets
    #      the reported 1.062 once the covariate- and dose-driven spread that
    #      the model explains is also accounted for.
    # Note that whether 1.062 is read as the SD or as the variance barely
    # matters numerically here (SD 1.062 vs sqrt(1.062) = 1.031); the
    # decision that matters is log scale versus linear scale.
    #
    # A residual this large is the direct consequence of the model carrying a
    # random effect on clearance ONLY: all unexplained absorption variability
    # in a compound whose exposure moves ~28-fold with the meal has nowhere
    # to go but the residual.
    # ------------------------------------------------------------------------
    expSd <- 1.062 ;  label("Log-scale additive (log-normal) residual SD (unitless)")
    # Hu 2025 Table 2 row 9 "Additive residual" = 1.062 (ASE 0.011, 90% CI 1.043; 1.081)
  })

  model({
    # ---------------------------------------------------------------------
    # Covariate reference values, both stated verbatim in Hu 2025
    # Section 2.3: "The median BMI value of the population is 21.6 kg/m2,
    # while the reference body weight value is set to 65 kg in the
    # population."
    # ---------------------------------------------------------------------
    wt_ref  <- 65      # kg
    bmi_ref <- 21.6    # kg/m^2

    # ---------------------------------------------------------------------
    # Relative bioavailability. Reference condition = medium-fat meal
    # (FED = 1, FED_LIQUIDSUPP = 0) in a participant without Friedreich's
    # ataxia (DIS_FRDA = 0), at which frel = 1 exactly.
    #
    # The fasted term uses (1 - FED) because the canonical FED column is
    # oriented fed = 1 whereas the paper's FST indicator is oriented
    # fasted = 1; see covariateData[["FED"]]$notes.
    # ---------------------------------------------------------------------
    frel <- exp(lfdepot +
                  e_fed_liquidsupp_fdepot * FED_LIQUIDSUPP +
                  e_fasted_fdepot * (1 - FED)) *
      (1 + e_dis_frda_fdepot * DIS_FRDA)

    # ---------------------------------------------------------------------
    # Individual disposition parameters.
    # ---------------------------------------------------------------------
    cl <- exp(lcl + etalcl +
                e_conmed_itraconazole_cl * CONMED_ITRACONAZOLE +
                e_conmed_rifampicin_cl * CONMED_RIFAMPICIN) *
      (1 + e_dis_frda_cl * DIS_FRDA) *
      (WT / wt_ref) ^ e_wt_cl *
      (BMI / bmi_ref) ^ e_bmi_cl
    vc <- exp(lvc) * (WT / wt_ref) ^ e_wt_vc
    q  <- exp(lq)
    vp <- exp(lvp)

    # Absorption parameters.
    ka    <- exp(lka)
    tlag1 <- exp(ltlag1)
    d0    <- exp(ld0)
    ffo   <- expit(logitffo)

    # Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---------------------------------------------------------------------
    # ODE system. Two parallel absorption arms feed `central`:
    #   - `depot`   first-order, rate Ka, lag TLAG1, dose fraction frel * ffo
    #   - `central` zero-order over TK0 hours,       dose fraction frel * (1 - ffo)
    #
    # The zero-order arm is dosed DIRECTLY into central, matching the paper's
    # D2 (duration on compartment 2). Simulations must therefore supply TWO
    # dose records per administration: a bolus to `depot` and a modelled-
    # duration record (rate = -2) of the same amount to `central`; f() splits
    # the two arms and no dose is double-counted. Hu 2025 fitted no
    # intravenous data, so declaring f(central) cannot contaminate an IV
    # route here.
    #
    # Figure 5 labels a TLAG2 (lag time for the zero-order process) in the
    # schematic, but Table 2 has no TLAG2 row and the theta numbering 1-17 is
    # gapless, so it was not estimated; the zero-order input starts at the
    # dose time.
    # ---------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-               k12 * central - k21 * peripheral1

    f(depot)     <- frel * ffo
    alag(depot)  <- tlag1

    f(central)   <- frel * (1 - ffo)
    dur(central) <- d0

    # ---------------------------------------------------------------------
    # Observation. Doses are in mg and vc is in L, so central/vc is mg/L;
    # the factor 1000 converts to the ng/mL used on every concentration axis
    # in Hu 2025 (Figure 2 VPC, Supplementary Figures S1 and S4).
    # ---------------------------------------------------------------------
    Cc <- central / vc * 1000

    Cc ~ lnorm(expSd)
  })
}
