Riggs_2013_empagliflozin <- function() {
  description <- paste(
    "Population pharmacokinetic model for empagliflozin in patients with",
    "type 2 diabetes mellitus, pooled from five randomised placebo-controlled",
    "Phase I/II multiple-oral-dose studies (N = 974; 1-100 mg once daily for",
    "up to 12 weeks). Two-compartment disposition with lagged first-order oral",
    "absorption (absorption lag fixed at 0.5 h) and first-order elimination.",
    "Encodes the prespecified FULL covariate model: an imposed allometric body",
    "weight effect on CL/F, V2/F, Q/F and V3/F, plus age, sex, Asian race,",
    "total protein, serum creatinine, smoking history and alcohol history on",
    "the apparent disposition parameters and an Asian-race effect on ka.",
    "Residual error is study-group dependent (Studies A/B/C vs D/E)."
  )
  reference <- paste(
    "Riggs MM, Staab A, Seman L, MacGregor TR, Bergsma TT, Gastonguay MR,",
    "Macha S. Population Pharmacokinetics of Empagliflozin, a Sodium Glucose",
    "Cotransporter 2 Inhibitor, in Patients With Type 2 Diabetes.",
    "J Clin Pharmacol. 2013;53(10):1028-1038. doi:10.1002/jcph.147."
  )
  vignette <- "Riggs_2013_empagliflozin"
  units <- list(
    time          = "h",
    dosing        = "mg empagliflozin (oral, once daily)",
    concentration = "Cc in nmol/L (converted from mg/L via MW 450.91 g/mol)"
  )

  compartmentData <- list(
    depot       = list(analyte = "empagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "empagliflozin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "empagliflozin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric effect imposed (not estimated) on CL/F, V2/F, Q/F and V3/F,",
        "normalised to a 70 kg reference individual: exponent FIXED at 0.75 for the",
        "clearances (CL/F, Q/F) and at 1 for the volumes (V2/F, V3/F) -- Riggs 2013",
        "Table 2 rows theta_20, theta_21, theta_22, theta_23 (all marked FIXED) and",
        "Methods equation 1. Cohort mean (min-max) weight was 85 (44-152) kg."
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-form effect (AGE/50)^theta on CL/F, V2/F and V3/F; reference 50 years",
        "(Riggs 2013 Figure 3 caption reference individual). Cohort mean (min-max)",
        "age was 58 (28-80) years."
      ),
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = male (the reference individual is a man)",
      notes              = paste(
        "Multiplicative effect on CL/F, V2/F and V3/F (Riggs 2013 Table 2 rows",
        "theta_24, theta_25, theta_26). The paper's dummy is named SEX with the",
        "reference individual defined as a man (Figure 3 caption), and Figure 3",
        "labels the contrast 'Female' -- so SEX = 1 denotes female and maps directly",
        "onto the canonical SEXF with no value inversion. Cohort was 574 male /",
        "400 female (41.1% female)."
      ),
      source_name        = "SEX"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = non-Asian",
      notes              = paste(
        "Multiplicative effect on CL/F, V2/F, V3/F and ka (Riggs 2013 Table 2 rows",
        "theta_11, theta_15, theta_18, theta_27). The paper defines the group as",
        "Japanese, Taiwanese and Korean patients versus non-Asian (Methods,",
        "'Race was specifically investigated for exposure differences between the",
        "Asian (Japanese, Taiwanese, and Koreans) and non-Asian patients'), and the",
        "reference individual is Caucasian. 213/974 (22%) of the cohort were Asian."
      ),
      source_name        = "ASIAN"
    ),
    TPRO = list(
      description        = "Total serum protein at baseline",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-form effect (TPRO/68 g/L)^theta on CL/F, V2/F and V3/F. The paper",
        "reports total protein in g/dL and normalises to 6.8 g/dL (Riggs 2013",
        "Table 2 rows theta_12, theta_16, theta_19; Figure 3 x-axis labels read",
        "'Total protein = 5.0 / 6.5 / 8.0 g/dL'). 6.8 g/dL = 68 g/L in the canonical",
        "SI unit, matching the reference used in Johnston_2021_empagliflozin_popPK.",
        "Note the Figure 3 caption's 'total protein of 6.8 mg/dL' is a typographical",
        "error -- Table 1 reports total protein in g/dL (study means 6.62-7.42)."
      ),
      source_name        = "TPRO"
    ),
    CREAT = list(
      description        = "Serum creatinine at baseline",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters CL/F as the INVERTED ratio (0.8 / CREAT)^theta_13 -- the paper uses",
        "'a simple inverse relationship with serum creatinine ... as a surrogate for",
        "renal filtration-associated clearance' (Methods) because estimated creatinine",
        "clearance was collinear with weight and age. Reference 0.8 mg/dL (Figure 3",
        "caption reference individual; the caption's 'serum creatinine of 0.8 g/dL' is",
        "a typographical error -- Table 1 and the Figure 3 x-axis labels both use",
        "mg/dL). Encoded as printed rather than algebraically flipped, so a rising",
        "serum creatinine lowers CL/F.",
        "UNITS ARE mg/dL FOR THIS MODEL, not the umol/L alternative the canonical",
        "register also permits."
      ),
      source_name        = "SCR"
    ),
    SMOKE_NEVER = list(
      description        = "Never-smoker indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = former or current smoker",
      notes              = paste(
        "Riggs 2013 codes smoking history as two dummies SMK1 and SMK2 against a",
        "NEVER-smoker reference (the reference individual is 'non-smoking', Figure 3",
        "caption), whereas the canonical SMOKE_NEVER + SMOKE_CURRENT pair is defined",
        "against a FORMER-smoker reference. Rather than rebase the published",
        "coefficients, model() reconstructs the paper's former-smoker dummy as",
        "(1 - SMOKE_NEVER - SMOKE_CURRENT). This keeps the published base CL/F",
        "(9.87 L/h) and both published multipliers exactly as printed -- the",
        "derivation the SMOKE_NEVER register entry itself recommends ('keep the",
        "published coefficient and intercept unchanged, rather than re-parameterising').",
        "Supply SMOKE_NEVER = 1, SMOKE_CURRENT = 0 for a never-smoker (the paper's",
        "reference); both 0 for a former smoker; SMOKE_CURRENT = 1 for a current smoker."
      ),
      source_name        = "SMK1 (as 1 - SMOKE_NEVER - SMOKE_CURRENT)"
    ),
    SMOKE_CURRENT = list(
      description        = "Current-smoker indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = former or never smoker",
      notes              = paste(
        "The paper's SMK2 dummy (Riggs 2013 Table 2 theta_9 = 1.02). Table 2 lists the",
        "CL/F covariate rows in the order SMK1, SMK2, ALC, and Figure 3A lists the",
        "matching categorical effects in the order 'Ex-smoker', 'Current smoker',",
        "'Consumes alcohol' -- which pins SMK1 = ex-smoker and SMK2 = current smoker.",
        "The pairing is corroborated by the bootstrap percentages printed in Figure 3A:",
        "the Ex-smoker row reads '66 34' (majority below the null, consistent with",
        "theta_8 = 0.991 < 1) and the Current smoker row reads '12 88' (majority above",
        "the null, consistent with theta_9 = 1.02 > 1).",
        "See SMOKE_NEVER notes for the reference-category reconstruction."
      ),
      source_name        = "SMK2"
    ),
    ALCOHOL_USE = list(
      description        = "Any history of alcohol use",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = no history of alcohol use",
      notes              = paste(
        "The paper's ALC dummy (Riggs 2013 Table 2 theta_10 = 1.02), labelled",
        "'Consumes alcohol' in Figure 3A, against a reference individual with 'no",
        "history of alcohol use' (Figure 3 caption). This is ANY alcohol consumption,",
        "NOT the chronic-abuse concept encoded by the separate ALCOHOL_ABUSE canonical",
        "(>= 6 units/day, Swart 2004 / VUMC criteria) -- the two columns must not be",
        "pooled. The paper does not state an ascertainment threshold beyond",
        "'alcohol histories' (Methods) and concludes the effect is null",
        "(95% CI 0.984-1.06, wholly inside the 75-125% null band)."
      ),
      source_name        = "ALC"
    ),
    STUDY_DE = list(
      description        = "Riggs 2013 Study D / Study E cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = Studies A, B or C",
      notes              = paste(
        "Selects which of the two estimated residual-error pairs applies to an",
        "observation. Riggs 2013 estimated separate residual variance terms",
        "'to account for the greater variability in the empagliflozin concentrations",
        "observed from Studies D and E compared with Studies A, B, and C' (Methods).",
        "1 = Study D (NCT00789035, 12-week Phase II) or Study E (NCT00749190,",
        "12-week metformin add-on Phase II); 0 = Study A (EudraCT 2007-000654-32,",
        "8-day), Study B (NCT00558571, 4-week) or Study C (NCT00885118, 4-week",
        "Japanese). Affects the residual error only -- no structural or disposition",
        "parameter depends on it. Supply per observation record."
      ),
      source_name        = "STUDY"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 974L,
    n_studies      = 5L,
    age_range      = "28-80 years",
    age_mean       = "58 years",
    weight_range   = "44-152 kg",
    weight_mean    = "85 kg",
    sex_female_pct = 41.1,
    race_ethnicity = c(White = 77.1, Asian = 21.9, Black = 0.8, `Hawaiian/Pacific` = 0.2),
    disease_state  = "type 2 diabetes mellitus",
    dose_range     = "1-100 mg empagliflozin once daily for 8 days to 12 weeks",
    n_observations = 8289L,
    renal_function = paste(
      "Mostly preserved: >67% of patients had estimated creatinine clearance",
      ">90 mL/min and only 14 patients (<1.5%) had estimated creatinine clearance",
      "<50 mL/min (Riggs 2013 Results)."
    ),
    notes          = paste(
      "Baseline demographics and laboratory values are in Riggs 2013 Table 1,",
      "reported per study. Study A (EudraCT 2007-000654-32, n = 48, 2.5/10/25/100 mg",
      "q.d. x 8 days), Study B (NCT00558571, n = 78, 10/25/100 mg q.d. x 4 weeks),",
      "Study C (NCT00885118, n = 100, Japanese patients, 1/5/10/25 mg q.d. x 4 weeks),",
      "Study D (NCT00789035, n = 324, 5/10/25 mg q.d. x 12 weeks) and Study E",
      "(NCT00749190, n = 424, 1/5/10/25/50 mg q.d. added to metformin x 12 weeks).",
      "All five were double-blind and placebo-controlled. Patients contributed 8,289",
      "quantifiable plasma concentrations (1,116 / 2,240 / 1,029 / 1,390 / 2,514 from",
      "Studies A-E). The Asian patients came from Study C (n = 100, Japanese) and",
      "Study D (n = 112). Mean weight in Study C was lower (68 kg) than in the other",
      "studies. Note the model was fit to observations collected AFTER 1 hour",
      "post-dose only, with the absorption lag fixed at 0.5 h (Methods)."
    )
  )

  ini({
    # ---- Structural PK parameters (Riggs 2013 Table 2, bold rows = reference ----
    # ---- population estimates). The reference individual (Figure 3 caption) is a ----
    # ---- 50-year-old, non-smoking, Caucasian man of 70 kg with no history of ----
    # ---- alcohol use, total protein 6.8 g/dL and serum creatinine 0.8 mg/dL. ----
    lcl <- log(9.87);  label("Apparent oral clearance CL/F (L/h) for the reference individual")   # Table 2 theta_1 = 9.87 L/h (bootstrap median 9.86; 95% CI 9.33-10.4)
    lvc <- log(3.02);  label("Apparent central volume V2/F (L) for the reference individual (V2 -> vc per canonical convention)")      # Table 2 theta_2 = 3.02 L (bootstrap median 3.06; 95% CI 2.40-3.76)
    lq  <- log(5.16);  label("Apparent inter-compartmental clearance Q/F (L/h) for the reference individual")  # Table 2 theta_3 = 5.16 L/h (bootstrap median 5.16; 95% CI 4.82-5.56)
    lvp <- log(60.4);  label("Apparent peripheral volume V3/F (L) for the reference individual (V3 -> vp per canonical convention)")   # Table 2 theta_4 = 60.4 L (bootstrap median 60.2; 95% CI 56.4-64.3)
    lka <- log(0.224); label("First-order absorption rate constant ka (1/h) for the reference individual")     # Table 2 theta_5 = 0.224 1/h (bootstrap median 0.227; 95% CI 0.176-0.279)

    # Absorption lag time. FIXED, not estimated: "the absorption lag in the current
    # analysis was fixed at 0.5 hours" (Methods) and Table 2 footnote "An oral
    # absorption lag time (ALAG1) was included and fixed at 0.5 hours."
    lalag <- fixed(log(0.5)); label("Absorption lag time ALAG1 (h)")  # Table 2 footnote; Methods

    # ---- Imposed allometric body-weight exponents. All four are FIXED, not ----
    # ---- estimated: Methods equation 1 states theta_allo "is a fixed allometric ----
    # ---- power parameter, which was assigned a value of 0.75 for physiologic ----
    # ---- processes, such as clearances, and a value of 1 for anatomical volumes", ----
    # ---- and Table 2 marks each of these four rows FIXED. ----
    e_wt_cl <- fixed(0.75); label("Power exponent of (WT/70) on CL/F (unitless)")  # Table 2 theta_20 = 0.75 FIXED
    e_wt_vc <- fixed(1);    label("Power exponent of (WT/70) on V2/F (unitless)")  # Table 2 theta_21 = 1 FIXED
    e_wt_q  <- fixed(0.75); label("Power exponent of (WT/70) on Q/F (unitless)")   # Table 2 theta_22 = 0.75 FIXED
    e_wt_vp <- fixed(1);    label("Power exponent of (WT/70) on V3/F (unitless)")  # Table 2 theta_23 = 1 FIXED

    # ---- Continuous (power-form) covariate effects on CL/F (Riggs 2013 Table 2) ----
    # NOTE ON SIGNS: this journal's PDF renders the minus sign as a C0 control byte
    # that most text extractors drop. Every negative exponent below was recovered
    # from the raw byte stream (0x05 = minus, 0x04 = multiplication sign) and
    # independently corroborated by the direction of the bootstrap densities printed
    # in Figure 3 -- e.g. a 30-year-old has HIGHER CL/F than the 50-year reference
    # (Figure 3A "AGE = 30 Yr" density sits right of the null), which is only
    # possible if theta_7 is negative.
    e_age_cl   <- -0.192; label("Power exponent of (AGE/50) on CL/F")            # Table 2 theta_7 = -0.192 (95% CI -0.338, -0.0313)
    e_tpro_cl  <- -0.345; label("Power exponent of (TPRO/68 g/L) on CL/F")       # Table 2 theta_12 = -0.345 (95% CI -0.490, -0.190)
    e_creat_cl <-  0.249; label("Power exponent of the INVERTED ratio (0.8 mg/dL / CREAT) on CL/F")  # Table 2 theta_13 = 0.249 (95% CI 0.168, 0.334)

    # ---- Categorical (multiplicative) covariate effects on CL/F (Table 2) ----
    cl_smoke_former  <- 0.991; label("Multiplicative effect of former (ex-) smoker on CL/F; ref = never smoker")  # Table 2 theta_8 (SMK1) = 0.991 (95% CI 0.953, 1.04)
    cl_smoke_current <- 1.02;  label("Multiplicative effect of current smoker on CL/F; ref = never smoker")       # Table 2 theta_9 (SMK2) = 1.02 (95% CI 0.981, 1.07)
    cl_alcohol       <- 1.02;  label("Multiplicative effect of any alcohol-use history on CL/F; ref = no history") # Table 2 theta_10 (ALC) = 1.02 (95% CI 0.984, 1.06)
    cl_asian         <- 0.98;  label("Multiplicative effect of Asian race on CL/F; ref = non-Asian")              # Table 2 theta_11 = 0.98 (95% CI 0.938, 1.03)
    cl_female        <- 0.988; label("Multiplicative effect of female sex on CL/F; ref = male")                   # Table 2 theta_24 = 0.988 (95% CI 0.942, 1.04)

    # ---- Covariate effects on V2/F (Table 2) ----
    e_age_vc  <-  0.0807; label("Power exponent of (AGE/50) on V2/F")       # Table 2 theta_14 = 0.0807 (95% CI -0.0985, 0.256)
    e_tpro_vc <- -0.0331; label("Power exponent of (TPRO/68 g/L) on V2/F")  # Table 2 theta_16 = -0.0331 (95% CI -0.292, 0.234)
    vc_asian  <-  1.20;   label("Multiplicative effect of Asian race on V2/F; ref = non-Asian")  # Table 2 theta_15 = 1.20 (95% CI 0.925, 1.58)
    vc_female <-  0.948;  label("Multiplicative effect of female sex on V2/F; ref = male")       # Table 2 theta_25 = 0.948 (95% CI 0.901, 0.996)

    # ---- Covariate effects on V3/F (Table 2) ----
    e_age_vp  <-  0.246; label("Power exponent of (AGE/50) on V3/F")       # Table 2 theta_17 = 0.246 (95% CI -0.0146, 0.528)
    e_tpro_vp <- -0.336; label("Power exponent of (TPRO/68 g/L) on V3/F")  # Table 2 theta_19 = -0.336 (95% CI -0.678, -0.0217)
    vp_asian  <-  1.15;  label("Multiplicative effect of Asian race on V3/F; ref = non-Asian")  # Table 2 theta_18 = 1.15 (95% CI 1.06, 1.28)
    vp_female <-  1.02;  label("Multiplicative effect of female sex on V3/F; ref = male")       # Table 2 theta_26 = 1.02 (95% CI 0.916, 1.12)

    # ---- Covariate effect on ka (Table 2) ----
    ka_asian <- 1.24; label("Multiplicative effect of Asian race on ka; ref = non-Asian")  # Table 2 theta_27 = 1.24 (95% CI 0.939, 1.65)

    # ---- Interindividual variability ----
    # Riggs 2013 models IIV with "exponential variance models with covariance terms
    # between the apparent oral PK parameters (CL/F, V3/F, KA)" (Methods) and reports
    # the diagonal as CV% in Table 2 (CL/F 26.9, V3/F 30.8, KA 15.2) with the
    # correlations in Results: "0.30, 0.27 and -0.13 between eta_CL/F (ETA1) -
    # eta_V3/F (ETA2), eta_CL/F - eta_KA (ETA3), and eta_V3/F - eta_KA, respectively."
    # Back-transformed for a log-normal eta via omega^2 = log(1 + CV^2):
    #   CL/F  omega^2 = log(1 + 0.269^2) = 0.0698628
    #   V3/F  omega^2 = log(1 + 0.308^2) = 0.0906302
    #   KA    omega^2 = log(1 + 0.152^2) = 0.0228411
    # and the off-diagonals as corr * sqrt(omega^2_i * omega^2_j):
    #   cov(CL,V3) =  0.30 * sqrt(0.0698628 * 0.0906302) =  0.0238715
    #   cov(CL,KA) =  0.27 * sqrt(0.0698628 * 0.0228411) =  0.0107856
    #   cov(V3,KA) = -0.13 * sqrt(0.0906302 * 0.0228411) = -0.0059148
    # The resulting 3x3 block is positive definite (eigenvalues 0.1063, 0.0584, 0.0186).
    # V2/F and Q/F carry no IIV -- Table 2 leaves their IIV cells blank.
    etalcl + etalvp + etalka ~ c(
       0.0698628,
       0.0238715,  0.0906302,
       0.0107856, -0.0059148,  0.0228411
    )  # Table 2 IIV column (CV%) + Results correlations

    # ---- Residual error, stratified by study group ----
    # "Separate residual variance terms were included during model development to
    # account for the greater variability in the empagliflozin concentrations
    # observed from Studies D and E compared with Studies A, B, and C" (Methods).
    # Table 2 reports each group as a proportional CV% plus an additive SD in nmol/L;
    # the error model is "a combined additive (eps_a,ij) and proportional (eps_p,ij)
    # model" (Methods). CV% is taken as the proportional SD directly (19.6% -> 0.196).
    propSd_studyabc <- 0.196; label("Proportional residual error, Studies A/B/C (fraction)")  # Table 2 residual variance, Studies A/B/C: CV% 19.6 (95% CI 18.5, 20.8)
    addSd_studyabc  <- 0.179; label("Additive residual error, Studies A/B/C (nmol/L)")        # Table 2 residual variance, Studies A/B/C: SD 0.179 nmol/L (95% CI 0.010, 0.334)
    propSd_studyde  <- 0.357; label("Proportional residual error, Studies D/E (fraction)")    # Table 2 residual variance, Studies D/E: CV% 35.7 (95% CI 33.7, 38.0)
    addSd_studyde   <- 0.010; label("Additive residual error, Studies D/E (nmol/L)")          # Table 2 residual variance, Studies D/E: SD 0.010 nmol/L (95% CI 0.010, 0.010; at the estimation lower bound)
  })

  model({
    # ---- Constants ----
    # Empagliflozin molecular weight; small molecule (BI 10773; SGLT2 inhibitor;
    # CAS 864070-44-0, molecular formula C23H27ClO7) -> 450.91 g/mol. NOT reported
    # in Riggs 2013 -- an external chemical constant, carried here so that Cc is in
    # the paper's nmol/L unit. It is the same value used by the sibling models
    # Baron_2016_empagliflozin and Johnston_2021_empagliflozin_popPK, and it is
    # corroborated by the paper's own numbers: Results reports steady-state AUC for
    # 25 mg q.d. of 6,830 / 4,860 / 3,700 nmol*h/L at 54 / 85 / 123 kg, and
    # Dose / (CL/F) with this MW reproduces 6,824 / 4,856 / 3,681 (within 0.5%).
    mw_empa     <- 450.91
    nmol_per_mg <- 1e6 / mw_empa

    # ---- Covariate reference values (Riggs 2013 Table 2 normalisers and the ----
    # ---- Figure 3 caption reference individual) ----
    ref_wt        <- 70    # kg
    ref_age       <- 50    # years
    ref_tpro_gL   <- 68    # g/L (the paper's 6.8 g/dL in canonical SI units)
    ref_creat     <- 0.8   # mg/dL

    # Reconstruct the paper's former-smoker dummy (SMK1) from the canonical
    # never/current pair. Riggs 2013 uses a NEVER-smoker reference; the canonical
    # SMOKE_NEVER + SMOKE_CURRENT pair uses a FORMER-smoker reference. Deriving the
    # missing indicator keeps the published base CL/F and both published multipliers
    # exactly as printed, instead of rebasing them.
    #   never smoker   -> SMOKE_NEVER = 1, SMOKE_CURRENT = 0 -> smokeFormer = 0
    #   former smoker  -> SMOKE_NEVER = 0, SMOKE_CURRENT = 0 -> smokeFormer = 1
    #   current smoker -> SMOKE_NEVER = 0, SMOKE_CURRENT = 1 -> smokeFormer = 0
    smokeFormer <- 1 - SMOKE_NEVER - SMOKE_CURRENT

    # ---- Individual PK parameters ----
    # Riggs 2013 Methods equation 1: a typical value is the reference estimate times
    # an imposed allometric weight term, times a normalised power term per continuous
    # covariate, times a multiplicative factor raised to each 0/1 categorical dummy.
    cl <- exp(lcl + etalcl) *
          (WT / ref_wt)^e_wt_cl *
          (AGE / ref_age)^e_age_cl *
          (TPRO / ref_tpro_gL)^e_tpro_cl *
          (ref_creat / CREAT)^e_creat_cl *
          cl_smoke_former^smokeFormer *
          cl_smoke_current^SMOKE_CURRENT *
          cl_alcohol^ALCOHOL_USE *
          cl_asian^RACE_ASIAN *
          cl_female^SEXF

    vc <- exp(lvc) *
          (WT / ref_wt)^e_wt_vc *
          (AGE / ref_age)^e_age_vc *
          (TPRO / ref_tpro_gL)^e_tpro_vc *
          vc_asian^RACE_ASIAN *
          vc_female^SEXF

    vp <- exp(lvp + etalvp) *
          (WT / ref_wt)^e_wt_vp *
          (AGE / ref_age)^e_age_vp *
          (TPRO / ref_tpro_gL)^e_tpro_vp *
          vp_asian^RACE_ASIAN *
          vp_female^SEXF

    q  <- exp(lq) * (WT / ref_wt)^e_wt_q

    ka <- exp(lka + etalka) * ka_asian^RACE_ASIAN

    # ---- Micro-constants (NONMEM ADVAN4 TRANS4 parameterisation) ----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- ODE system: two-compartment disposition with first-order oral input ----
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Lagged absorption; ALAG1 fixed at 0.5 h.
    alag(depot) <- exp(lalag)

    # ---- Observation ----
    # The central state holds mg; convert to the paper's nmol/L reporting unit.
    Cc <- central / vc * nmol_per_mg

    # Study-group-dependent combined residual error. STUDY_DE selects between the
    # two estimated pairs; it enters nowhere else in the model.
    ruvProp <- propSd_studyabc * (1 - STUDY_DE) + propSd_studyde * STUDY_DE
    ruvAdd  <- addSd_studyabc  * (1 - STUDY_DE) + addSd_studyde  * STUDY_DE
    Cc ~ prop(ruvProp) + add(ruvAdd)
  })
}
