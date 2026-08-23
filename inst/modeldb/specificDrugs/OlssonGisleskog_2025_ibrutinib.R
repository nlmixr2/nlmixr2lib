OlssonGisleskog_2025_ibrutinib <- function() {
  description <- paste(
    "Two-compartment population PK model for oral ibrutinib (a covalent",
    "Bruton's tyrosine kinase inhibitor) in patients aged 65 years and",
    "older with previously untreated mantle cell lymphoma, treated with",
    "ibrutinib 560 mg once daily on top of bendamustine plus rituximab in",
    "the phase 3 SHINE study (Olsson Gisleskog 2025; N = 259 patients,",
    "2070 plasma concentrations). Absorption is sequential zero-order",
    "input into the gut compartment over a duration D1 followed by",
    "first-order transfer at ka into the central compartment, all after an",
    "absorption lag time ALAG1. Typical values are CL/F = 1123 L/h,",
    "V2/F = 613 L, Q/F = 1128 L/h, V3/F = 6673 L, ka = 0.719 1/h,",
    "ALAG1 = 0.167 h and D1 = 2.45 h under fasting or modified-fasting",
    "conditions (3.29 h fed). Relative bioavailability F1 is 1 in the",
    "modified-fasting and fed states and 0.666 under strict fasting, and",
    "carries two covariate effects that this analysis could not re-estimate",
    "and therefore held at the values of the previously developed model: a",
    "1.59-fold increase with concomitant CYP3A inhibitors and an",
    "(AGE/65)^0.699 power effect of age. Between-subject variability is",
    "carried on V2/F (133.4% CV), Q/F (55.02%), V3/F (27.74%), ALAG1",
    "(68.04%) and F1 (65.64%); CL/F, ka and D1 carry none in the final",
    "model. Residual unexplained variability is exponential and was fixed",
    "at 81.3% from the previous analysis rather than estimated. The paper's",
    "paired exposure-response analysis is a pair of static logistic",
    "regressions on AUCtau,ss (atrial fibrillation and any hemorrhage);",
    "being non-ODE statistical models they are reproduced in the paired",
    "vignette narrative rather than encoded here as separate model files.",
    sep = " "
  )
  reference <- paste(
    "Olsson Gisleskog P, Valenzuela B, Treijtel N, Deshpande S,",
    "Henninger T, Perez-Ruixo JJ. (2025). Population Pharmacokinetic and",
    "Exposure-Response Analyses of Ibrutinib Combined With Bendamustine",
    "and Rituximab in Patients With Mantle Cell Lymphoma.",
    "CPT Pharmacometrics Syst Pharmacol.",
    "doi:10.1002/psp4.70061",
    sep = " "
  )
  vignette <- "OlssonGisleskog_2025_ibrutinib"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Confirmed against Olsson Gisleskog 2025 Methods 2.2
  # ("sequential zero-order input into the gut compartment ... followed by
  # first-order absorption") and Methods 2.1 (ibrutinib PLASMA concentrations
  # by LC-MS/MS, LLOQ 0.5 ng/mL).
  compartmentData <- list(
    depot       = list(analyte = "ibrutinib", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "ibrutinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ibrutinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FASTED_STRICT = list(
      description        = "Strict-fasting dose-record indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (modified fasting or fed)",
      notes              = paste(
        "Per-dose-record indicator distinguishing the STRICTER of the two",
        "non-fed prandial states that Olsson Gisleskog 2025 Table 2",
        "contrasts. 1 = strict fasting, defined by Methods 2.1 as the",
        "PK-sampling-day instruction to 'fast from midnight prior (or at a",
        "minimum, 2 h prior) to dosing and continue fasting until",
        "approximately 30 min after capsule intake'. 0 = the paper's",
        "'modified fasting' routine-dosing state (Methods 2.1: ibrutinib",
        "'taken at least 30 min before or at least 2 h after a meal') or",
        "the fed state. Retained on F1 only, as the multiplicative factor",
        "0.666: Table 2 gives 'F1 mod fast/fed' = 1 FIX and 'F1 fast' =",
        "0.666 FIX, so strict fasting reduces relative bioavailability by",
        "33.4% while leaving the zero-order duration unchanged (Table 2",
        "pools 'D1 fast/mod fast' at a single 2.45 h). Two indicators are",
        "therefore needed to span the paper's three prandial states, and",
        "FASTED_STRICT is the one that separates strict from modified",
        "fasting while FED separates fed from both: strict fasting is",
        "FASTED_STRICT = 1, FED = 0; modified fasting is FASTED_STRICT = 0,",
        "FED = 0; fed is FASTED_STRICT = 0, FED = 1. Results 3.1 selected",
        "the strict-fasting condition for the analysis on objective-function",
        "grounds (-221.216 vs fed, -100.337 vs modified fasting), and",
        "Discussion paragraph 2 confirms 'it was assumed that ibrutinib was",
        "administered in the fasting condition', so FASTED_STRICT = 1 is the",
        "state the published typical exposures correspond to even though it",
        "is not this column's reference category. IMPORTANT: do NOT identify",
        "this indicator with the bare phrase 'modified fasting'. The phrase",
        "is not portable between papers -- Mauro 2025 nilotinib uses it for",
        "the opposite protocol (food first, dose about 2 h later, encoded as",
        "MEAL_PREDOSE_2H = 1), whereas here it means dose first, food later.",
        "Always read each paper's own definition."
      ),
      source_name        = "prandial state (Table 2 row labels 'fast' / 'mod fast' / 'fed')"
    ),
    FED = list(
      description        = "Fed-state dose-record indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasting or modified fasting)",
      notes              = paste(
        "Per-dose-record indicator. 1 = the dose was taken with food; 0 =",
        "strict fasting or modified fasting. Retained on the zero-order",
        "input duration D1 only: Olsson Gisleskog 2025 Table 2 gives",
        "'D1 fast/mod fast' = 2.45 h (%SEM 6.11) and 'D1 fed' = 3.29 h FIX,",
        "i.e. food lengthens the zero-order absorption window by 34% with no",
        "effect on relative bioavailability ('F1 mod fast/fed' is a single",
        "1 FIX row covering both states). The paper carries no meal-",
        "composition covariate, so FED rather than FED_HIGHFAT /",
        "FED_LOWFAT is the right column. In SHINE itself no dose was fed",
        "by protocol (Methods 2.1), so FED = 0 for every record in the",
        "analysis dataset; the fed D1 is inherited from the previously",
        "developed model (Marostica 2015, reference 15) and is retained here",
        "for completeness so the model can simulate the fed state. Pairs",
        "with FASTED_STRICT to span the paper's three prandial states."
      ),
      source_name        = "prandial state (Table 2 row labels 'fast' / 'mod fast' / 'fed')"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline, time-fixed. SHINE enrolled patients 65 years of age or",
        "older by protocol (Methods 2.1); the ibrutinib arm median is 71",
        "years, range 65-86 (Table 1). Enters F1 as the power effect",
        "(AGE/65)^0.699, held FIX at the previously developed model's",
        "estimate (Table 2 'Effect of age (power)' = 0.699 FIX; Results 3.1",
        "explains that re-estimating it 'led to minimization terminated due",
        "to rounding errors, indicating that insufficient data was",
        "available'). The 65-year centering value is NOT stated anywhere in",
        "the paper or reachable supplement and was back-solved from the",
        "paper's own reported typical exposure -- see the vignette Errata",
        "for the derivation and for the rejected alternatives. An uncentered",
        "power is arithmetically excluded (71^0.699 = 19.7 would give F1 =",
        "13 and AUCtau,ss of roughly 6500 ng.h/mL against a reported median",
        "of 349). Because the effect is on F1 rather than on CL/F it scales",
        "exposure directly: over the observed 65-86 year range the factor",
        "spans 1.00 to 1.20."
      ),
      source_name        = "AGE"
    ),
    CONMED_CYP3A4_INH = list(
      description        = "Concomitant CYP3A inhibitor coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant CYP3A inhibitor)",
      notes              = paste(
        "1 = the dose was taken with a concomitant CYP3A inhibitor of any",
        "strength; 0 = no concomitant CYP3A inhibitor. Table 1 splits the",
        "ibrutinib arm into none 79 (31%), weak or moderate 128 (49%) and",
        "strong 52 (20%); the model pools weak, moderate and strong into a",
        "single indicator, so no CONMED_CYP3A4_INH_STRONG /",
        "CONMED_CYP3A4_INH_MOD stratification is used. Enters F1 as the",
        "multiplicative ratio 1.59, held FIX at the previously developed",
        "model's estimate (Table 2 'Effect of CYP3A inhibitors (ratio)' =",
        "1.59 FIX); Discussion paragraph 8 restates it as 'a 59% increase in",
        "bioavailability'. TIME-VARYING, not subject-level: Methods 2.3",
        "states that the analysis accounted for comedication by 'having a",
        "weighted average dose based on the incidence of the doses given",
        "with and without the comedications; the same approach was applied",
        "to the apparent CL in deriving the AUC'. Reading it as a",
        "subject-level flag is arithmetically falsified -- applying 1.59 to",
        "the 69% of the arm that Table 1 records as inhibitor-exposed would",
        "put the typical AUCtau,ss near 528 ng.h/mL against the reported",
        "median of 349. Set the column per dose record from the",
        "comedication history."
      ),
      source_name        = "coadministered CYP3A inhibitor (Table 1)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tabulated in Olsson Gisleskog 2025 Table 1 (ibrutinib arm median",
        "71.9 kg, range 35-126) and screened as a prognostic factor in the",
        "exposure-response analysis, but carries no effect in the final",
        "population PK model. Note that this is a difference from the",
        "upstream Marostica 2015 model (reference 15), whose abstract",
        "reports a body-weight effect on volume; Results 3.1 states that",
        "only the CYP3A4-inhibitor and age effects on F1 were carried",
        "forward, and Discussion paragraph 2 that 'covariate effects on the",
        "PK parameters could not be formally quantified and were assumed to",
        "be the same as in the previous analysis'. Whether a weight effect",
        "was dropped or was never in the version of the previous model that",
        "this analysis started from cannot be settled from the sources on",
        "disk (Marostica 2015 is not open access)."
      ),
      source_name        = "Weight"
    ),
    BMI = list(
      description        = "Baseline body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Table 1 ibrutinib arm median 25.7 kg/m^2 (range 16.7-49.2).",
        "Reported as a baseline characteristic only; not a covariate in the",
        "population PK model and not among the prognostic factors carried",
        "into the exposure-response regressions."
      ),
      source_name        = "BMI"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Table 1 ibrutinib arm 82 / 259 female (32%). Screened as an",
        "exposure-response prognostic factor (Table 1 'E, S' annotation;",
        "subgroups in Table S3) and found not significant for either",
        "atrial fibrillation (Table S4) or any hemorrhage (Table S5). No",
        "effect on any PK parameter."
      ),
      source_name        = "Sex"
    ),
    CRCL = list(
      description        = "Baseline creatinine clearance",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Table 1 dichotomises the ibrutinib arm as CRCL >= 60 mL/min 177",
        "(68%) and CRCL < 60 mL/min 82 (32%); the continuous values are not",
        "tabulated. Screened as an exposure-response safety prognostic",
        "factor and not retained (Tables S4, S5). No effect on any PK",
        "parameter -- consistent with ibrutinib's negligible renal",
        "elimination."
      ),
      source_name        = "Renal function"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 259L,
    n_studies       = 1L,
    n_observations  = 2070L,
    age_range       = "65-86 years",
    age_median      = "71 years",
    weight_range    = "35-126 kg",
    weight_median   = "71.9 kg",
    bmi_median      = "25.7 kg/m^2 (range 16.7-49.2)",
    sex_female_pct  = 31.7,
    race_ethnicity  = c(
      `White, not Hispanic or Latino` = 71.0,
      `White, Hispanic or Latino`     = 4.6,
      Asian                           = 17.8,
      `Black, of African heritage or African American` = 0.8,
      Other                           = 5.8
    ),
    disease_state   = paste(
      "Previously untreated mantle cell lymphoma in patients 65 years of",
      "age or older. Simplified MCL international prognostic index low",
      "risk 52 (20%), intermediate 119 (46%), high 88 (34%); ECOG 0 133",
      "(51%) and 1-2 126 (49%); blastoid or pleomorphic histology 19 (7%);",
      "tumour bulk >= 5 cm 95 (37%); baseline TP53 mutated 26, unmutated",
      "114, unknown 119 (Olsson Gisleskog 2025 Table 1, ibrutinib column)."
    ),
    dose_range      = paste(
      "Ibrutinib 560 mg orally once daily, taken at least 30 min before or",
      "at least 2 h after a meal, continued until progressive disease or",
      "unacceptable toxicity. Reduction to 420 mg or 280 mg once daily was",
      "permitted to manage toxicity, so the individual average daily dose",
      "used for the exposure-response analysis spans 280-560 mg. All",
      "patients also received open-label bendamustine 90 mg/m2 IV on days 1",
      "and 2 of each 28-day cycle plus rituximab 375 mg/m2 IV on day 1, for",
      "up to 6 cycles, and responders continued rituximab maintenance",
      "every 8 weeks for up to 12 further doses."
    ),
    regions         = "Global (phase 3 SHINE, NCT01776840)",
    hepatic_function = paste(
      "By NCI ODWG liver-function classification (Table 1): normal 231",
      "(89%), mild or moderate impairment 28 (11%), severe 0."
    ),
    renal_function  = paste(
      "Creatinine clearance >= 60 mL/min 177 (68%), < 60 mL/min 82 (32%)",
      "(Table 1)."
    ),
    co_medication   = paste(
      "Concomitant CYP3A inhibitor at any point: none 79 (31%), weak or",
      "moderate 128 (49%), strong 52 (20%) (Table 1). All patients received",
      "bendamustine and rituximab background therapy."
    ),
    notes           = paste(
      "The PK dataset is the 259 patients of the SHINE ibrutinib arm; only",
      "ibrutinib-arm samples were assayed (Methods 2.1). Sampling was",
      "sparse: pre-dose on day 2 of cycles 1, 2 and 3, plus 1, 2 and 4 h",
      "post-dose on day 2 of cycles 1 and 2. Plasma ibrutinib was measured",
      "by LC-MS/MS with an LLOQ of 0.5 ng/mL. The exposure-response dataset",
      "adds the 260 placebo patients for a total of 519. Parameter values",
      "implemented here are the 'Final (updated) population PK model'",
      "columns of Table 2, fit in NONMEM 7.3.0 -- a re-estimation of the",
      "previously developed model of reference 15 (Marostica 2015, patients",
      "with other B-cell malignancies) after external evaluation of that",
      "model on the SHINE data was judged unsuccessful. The final model",
      "improved the objective function by 110.626 points relative to the",
      "previous model. Reported individual model-predicted exposures are",
      "AUCtau,ss 425 (267) ng.h/mL, Cmax 74.5 (48.3) ng/mL and Ctrough",
      "3.90 (2.64) ng/mL, mean (SD); Table 3 gives AUCtau,ss quartile",
      "boundaries of 57.7 / 234 / 349 / 581 / 1864 ng.h/mL."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # Structural PK -- Olsson Gisleskog 2025 Table 2, "Final (updated)
    # population PK model / Population mean estimate (%SEM)" column. All
    # disposition parameters are apparent (i.e. /F): the analysis used oral
    # data only, and F1 is a RELATIVE bioavailability whose reference state
    # (modified fasting or fed) is fixed to 1.
    #
    # Units are h / L / L/h; concentration is converted to ng/mL inside
    # model() (see the observation comment there).
    # -----------------------------------------------------------------------
    lcl   <- log(1123) ; label("Apparent clearance CL/F (L/h)")                             # Table 2 final CL/F = 1123 L/h (%SEM 4.83; previous model 1100, i.e. +2%, matching Results 3.1 "no large difference ... for CL/F (2%)")
    lvc   <- log(613)  ; label("Apparent central volume of distribution V2/F (L)")          # Table 2 final V2/F = 613 L (%SEM 19.7; previous 436)
    lq    <- log(1128) ; label("Apparent intercompartmental clearance Q/F (L/h)")           # Table 2 final Q/F = 1128 L/h (%SEM 10.9; previous 709)
    lvp   <- log(6673) ; label("Apparent peripheral volume of distribution V3/F (L)")       # Table 2 final V3/F = 6673 L (%SEM 8.30; previous 8660). V2/F + V3/F = 7286 L, the Vss quoted in Results 3.1, and 7286/9096 = 0.80, matching the stated 20% reduction
    lka   <- log(0.719); label("First-order absorption rate constant ka (1/h)")             # Table 2 final ka = 0.719 1/h (%SEM 14.1; previous 0.451)
    lalag <- log(0.167); label("Absorption lag time ALAG1 (h)")                             # Table 2 final ALAG1 = 0.167 h (%SEM 21.38; previous 0.226, and 0.167/0.226 = 0.74, matching the stated 26% shorter Tlag in Results 3.1)

    # D1, the duration of the zero-order input into the gut compartment, is
    # tabulated as two absolute values rather than as a covariate ratio, so
    # both printed values are entered as their own parameters and mixed by
    # the binary FED indicator in model(). Table 2 pools strict fasting and
    # modified fasting into a single "D1 fast/mod fast" row, which is why
    # FED (and not FASTED_STRICT) is the indicator that acts on D1.
    ld1    <- log(2.45)        ; label("Zero-order input duration D1 when not fed (h)")     # Table 2 final "D1 fast/mod fast" = 2.45 h (%SEM 6.11; previous 1.29, and 2.45/1.29 = 1.90, matching the stated 90% longer zero-order duration in Results 3.1)
    ld1fed <- fixed(log(3.29)) ; label("Zero-order input duration D1 when fed (h)")         # Table 2 "D1 fed" = 3.29 FIX in BOTH the previous and the final model column (no %SEM reported)

    # -----------------------------------------------------------------------
    # Relative bioavailability F1 and its covariate effects.
    #
    # Table 2 reports F1 as two rows: "F1 mod fast/fed (fixed)" = 1 FIX and
    # "F1 fast" = 0.666 FIX. The reference state of the model is therefore
    # the modified-fasting / fed state, where F1 = 1, and strict fasting is
    # the deviation. lfdepot is entered as fixed(log(1)) accordingly -- a
    # log-scale zero, kept explicit so the FIX status is visible rather than
    # implied by the parameter's absence.
    #
    # Both covariate effects were held at the previously developed model's
    # estimates. Results 3.1: "the previously identified covariate effects,
    # CYP3A4 inhibitors and age effects on F1, were fixed to previous
    # estimates as the attempts to re-estimate these effects led to
    # minimization terminated due to rounding errors".
    # -----------------------------------------------------------------------
    lfdepot                    <- fixed(log(1))     ; label("Relative bioavailability F1 in the reference (modified-fasting or fed) state")  # Table 2 "F1 mod fast/fed (fixed)" = 1 FIX
    e_fasted_strict_fdepot     <- fixed(0.666)      ; label("Strict-fasting multiplicative effect on F1 (ratio)")                            # Table 2 "F1 fast" = 0.666 FIX -- a 33.4% lower relative bioavailability than the modified-fasting/fed reference
    e_conmed_cyp3a4_inh_fdepot <- fixed(1.59)       ; label("Concomitant CYP3A inhibitor multiplicative effect on F1 (ratio)")                # Table 2 "Effect of CYP3A inhibitors (ratio)" = 1.59 FIX in the final model (1.59, %SEM 4.82, when estimated in the previous model); Discussion paragraph 8 restates it as "a 59% increase in bioavailability"
    e_age_fdepot               <- fixed(0.699)      ; label("Power effect of age on F1 (unitless, reference 65 years)")                       # Table 2 "Effect of age (power)" = 0.699 FIX in the final model (0.699, %SEM 19.2, when estimated in the previous model). The 65-year centering is NOT reported and was back-solved -- see the vignette Errata

    # -----------------------------------------------------------------------
    # Between-subject variability. Table 2 reports BSV as "%CV" for exactly
    # five parameters in the final model -- V2/F 133.4, Q/F 55.02, V3/F
    # 27.74, ALAG1 68.04 and F1 65.64 -- and a dash for CL/F, ka and D1. The
    # Table 2 shrinkage note lists shrinkage for the same five and only
    # those five ("for V2, Q, V3, ALAG1, and F1 was 39%, 65%, 66%, 67%, and
    # 97%"), independently confirming that the final model carries five etas
    # and that CL/F, ka and D1 carry none. No eta is invented for them here.
    #
    # Methods 2.2 specifies a LOGNORMAL between-subject error model, so the
    # tabulated %CV is converted with the exact lognormal identity
    #     omega^2 = log(1 + CV^2)
    # rather than the approximate omega = CV/100. The choice matters a great
    # deal at the V2/F magnitude (omega 1.011 exact vs 1.334 approximate)
    # and is settled empirically against the paper's own reported exposure
    # distribution: because there is no eta on CL/F, AUCtau,ss = dose * F1 /
    # (CL/F) inherits the F1 lognormal exactly, so the F1 row predicts the
    # AUCtau,ss coefficient of variation directly. The exact reading predicts
    # CV 65.6% and a mean/median ratio of 1.196; the approximate reading
    # predicts 73.4% and 1.240. Results 3.1 and Table 3 report mean (SD)
    # 425 (267) ng.h/mL, i.e. CV 62.8%, with median 349, i.e. mean/median
    # 1.218, and a lower quartile of 234 against the exact reading's 236.
    # The exact reading is the one consistent with the published numbers.
    # See the vignette Errata for the full comparison.
    #
    #   V2/F   CV 1.3340 -> omega^2 = log(1 + 1.3340^2) = 1.0222912
    #   Q/F    CV 0.5502 -> omega^2 = log(1 + 0.5502^2) = 0.2644544
    #   V3/F   CV 0.2774 -> omega^2 = log(1 + 0.2774^2) = 0.0741337
    #   ALAG1  CV 0.6804 -> omega^2 = log(1 + 0.6804^2) = 0.3804510
    #   F1     CV 0.6564 -> omega^2 = log(1 + 0.6564^2) = 0.3582763
    #
    # Table 2 reports no off-diagonal covariance terms and no OMEGA block
    # structure, so the five etas are entered as independent variances.
    # -----------------------------------------------------------------------
    etalvc     ~ 1.0222912  ; label("BSV on V2/F")   # Table 2 final V2/F BSV = 133.4 %CV (%SEM 10.41; previous 154.6); shrinkage 39%
    etalq      ~ 0.2644544  ; label("BSV on Q/F")    # Table 2 final Q/F BSV = 55.02 %CV (%SEM 14.45; previous 64.9); shrinkage 65%
    etalvp     ~ 0.0741337  ; label("BSV on V3/F")   # Table 2 final V3/F BSV = 27.74 %CV (%SEM 14.07; previous 56.0); shrinkage 66%
    etalalag   ~ 0.3804510  ; label("BSV on ALAG1")  # Table 2 final ALAG1 BSV = 68.04 %CV (%SEM 12.62; previous 72.5); shrinkage 67%
    etalfdepot ~ 0.3582763  ; label("BSV on F1")     # Table 2 final F1 BSV = 65.64 %CV (%SEM 5.37; previous 67.5); shrinkage 97%

    # -----------------------------------------------------------------------
    # Residual unexplained variability. Methods 2.2: "An exponential error
    # model was used to quantify the residual unexplained variability of the
    # ibrutinib plasma concentrations", which is nlmixr2's lnorm() error.
    #
    # Table 2 tabulates NO residual-error row for either model, so the only
    # RUV magnitude anywhere on disk is Results 3.1: "given the sparseness
    # of the PK data, a residual variability (RUV) of 81.3%, determined from
    # previous studies [15], was selected". It is entered as fixed()
    # accordingly -- it was imported from the previous analysis, not
    # estimated here. Methods 2.2 does list "the residual variability" among
    # the things that were candidates for re-estimation in the model
    # refinement step, so it is possible that the final model's RUV differs
    # from 81.3%; no such value is reported. Recorded in the vignette Errata.
    # -----------------------------------------------------------------------
    expSd <- fixed(0.813) ; label("Exponential residual error, SD on the log-transformed concentration scale")  # Results 3.1 "a residual variability (RUV) of 81.3%, determined from previous studies [15], was selected"
  })

  model({
    # ---------------------------------------------------------------------
    # 1. Individual PK parameters.
    #
    # CL/F, ka, Q/F and V3/F are covariate-free; V2/F carries only its eta.
    # ---------------------------------------------------------------------
    cl  <- exp(lcl)
    vc  <- exp(lvc + etalvc)
    q   <- exp(lq + etalq)
    vp  <- exp(lvp + etalvp)
    ka  <- exp(lka)
    tlag <- exp(lalag + etalalag)

    # Zero-order input duration: 2.45 h fasting or modified fasting, 3.29 h
    # fed. FED is binary, so mixing the two absolute durations linearly
    # reproduces each tabulated value exactly. No eta -- Table 2 reports no
    # BSV on D1 in the final model.
    d1 <- exp(ld1) * (1 - FED) + exp(ld1fed) * FED

    # Relative bioavailability. The reference state is modified fasting or
    # fed (F1 = 1 FIX); strict fasting scales it by 0.666, a concomitant
    # CYP3A inhibitor by 1.59, and age by the power (AGE/65)^0.699. The two
    # ratio effects are raised to the power of their binary indicator so
    # each collapses to 1 in its reference category, which is the
    # algebraically identical multiplicative-ratio form the paper prints.
    f1 <- exp(lfdepot + etalfdepot) *
      e_fasted_strict_fdepot^FASTED_STRICT *
      e_conmed_cyp3a4_inh_fdepot^CONMED_CYP3A4_INH *
      (AGE / 65)^e_age_fdepot

    # 2. Micro-constants for the two-compartment system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---------------------------------------------------------------------
    # 3. Two-compartment ODE system with sequential zero-order then
    #    first-order oral absorption after a lag time (Methods 2.2). The
    #    dose enters `depot` at a constant rate over the zero-order window
    #    of duration D1, starting ALAG1 after the dose record, and `depot`
    #    then drains into `central` first-order at ka. Dose records must
    #    carry rate = -2 so that rxode2 uses the model's dur(depot); a plain
    #    bolus would collapse the zero-order phase and bias Cmax upward.
    # ---------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    dur(depot)  <- d1
    alag(depot) <- tlag
    f(depot)    <- f1

    # ---------------------------------------------------------------------
    # 4. Observation. Dose is in mg and vc in L, so central / vc is mg/L
    #    (= ug/mL); multiplying by 1000 gives ng/mL, the assay unit used
    #    throughout Olsson Gisleskog 2025 (Methods 2.1 LLOQ 0.5 ng/mL;
    #    Results 3.1 mean Cmax 74.5 ng/mL and Ctrough 3.90 ng/mL).
    # ---------------------------------------------------------------------
    Cc <- central / vc * 1000
    Cc ~ lnorm(expSd)
  })
}
