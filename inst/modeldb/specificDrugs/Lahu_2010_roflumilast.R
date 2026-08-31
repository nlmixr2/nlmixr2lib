Lahu_2010_roflumilast <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral roflumilast and",
    "its primary active metabolite roflumilast N-oxide in adult healthy",
    "volunteers and patients with moderate-to-severe COPD (Lahu 2010).",
    "Roflumilast is described by a two-compartment model with first-order",
    "absorption and a lag time; the absolute parent bioavailability is",
    "not identifiable and is fixed at F1 = 1. Roflumilast N-oxide is",
    "described by a one-compartment model with zero-order absorption",
    "(duration D1) and a lag time, with relative bioavailability Frel",
    "fixed at 1 for the null-covariate reference (also non-identifiable).",
    "Retained covariates on roflumilast parameters are food on tlag and",
    "ka, sex / smoking / race-Black / race-Hispanic / COPD on CL, and",
    "COPD on V1. Retained covariates on roflumilast N-oxide parameters",
    "are food on D1; age / sex / smoking / COPD on CL; body weight and",
    "COPD on Vd; and age / sex / race-Black / race-Hispanic on Frel.",
    "Inter-individual variability is reported on parent tlag, ka, CL,",
    "V1, Q, V2 (with a Q-V2 covariance) and on N-oxide D1, CL, Vd",
    "(with a full 3x3 covariance block); no IIV is reported on N-oxide",
    "tlag or on Frel. Residual error is proportional on the linear-",
    "concentration scale (additive on the log-transformed observation)",
    "for both observed analytes, fitted on the phase I dataset (the more",
    "data-rich layer).")
  reference <- "Lahu G, Hunnemeyer A, Diletti E, Elmlinger M, Ruth P, Zech K, McCracken N, Facius A. Population pharmacokinetic modelling of roflumilast and roflumilast N-oxide by total phosphodiesterase-4 inhibitory activity and development of a population pharmacodynamic-adverse event model. Clin Pharmacokinet. 2010;49(9):589-606. doi:10.2165/11536600-000000000-00000"
  vignette <- "Lahu_2010_roflumilast"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot          = list(analyte = "roflumilast", units = "ug", specimen = "administration site", verified = FALSE),
    central        = list(analyte = "roflumilast", units = "ug", specimen = "plasma", verified = FALSE),
    peripheral1    = list(analyte = "roflumilast", units = "ug", specimen = "plasma", verified = FALSE),
    central_noxide = list(analyte = "roflumilast N-oxide", units = "ug", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    AGE = list(
      description        = "Subject age (baseline; constant within an individual in the source dataset).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling (AGE/40)^-0.471 on N-oxide CL and (AGE/40)^-0.269 on N-oxide Frel; reference age 40 years per Lahu 2010 Methods (page 591). The two age effects on N-oxide CL and Frel partially cancel so that AUC_N-oxide increases approximately as (AGE/40)^0.202 with age (Lahu 2010 Results page 598). Age on the PARENT CL was tested in the full model and dropped by the WAM selection (theta_12 = -0.140, SE 0.0550; FDA NDA 22-522 review Table 6), which is why no age effect appears on parent CL in the shipped model. Cohort age was 38.2 +/- 14.7 years, so the 40-year reference sits essentially at the cohort mean.",
      source_name        = "Age"
    ),
    WT = list(
      description        = "Total body weight (baseline; constant within an individual in the source dataset).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling (WT/70)^1.00 on N-oxide Vd; reference weight 70 kg per Lahu 2010 equation 7. Body weight on parent V1 was tested in the full model (theta_17 = 0.497, SE 0.161 per FDA NDA 22-522 review Table 6) but dropped from the WAM-selected final model. Body weight on N-oxide CL was likewise tested and dropped (theta_9 = -0.497, SE 0.185; FDA review Table 7). Body weight on parent CL competed with sex in the full model and was removed to stabilise the full-model fit (Lahu 2010 Discussion page 602). The N-oxide Vd exponent itself moved from 0.824 (SE 0.142) in the full model to exactly 1.00 (SE 0.117) in the final-with-race model. Cohort weight was 75.4 +/- 11.2 kg, so the 70 kg reference is a round anchor rather than the cohort mean.",
      source_name        = "Weight"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male).",
      notes              = "Lahu 2010 uses the inverted convention Sex = 1 for male in equations 6 and 7. The canonical SEXF column inverts this; the model() block applies the source coefficient via the male indicator (1 - SEXF). Effects: +19.1% on parent CL when male (theta_13); +46.7% on N-oxide CL when male (theta_8); +23.1% on N-oxide Frel when male (theta_14). The paper's model-equation intercept (when all binary covariates equal 0) is therefore the FEMALE, non-smoking, White, healthy reference; the paper's narrative text instead describes the reference as 'male, non-smoking, White, healthy, 40-year-old subjects'. That contradiction is RESOLVED in favour of the equation reading: the FDA NDA 22-522 Clinical Pharmacology review (printed pages 77-80), reproducing the sponsor's own study report 114/2005, states the female reference three separate times -- 'the clearance (theta_3) is 10.5 +/- 0.490 L/hr for nonsmoking, non-black, non-Hispanic females', 'the metabolite CL (theta_3) is 0.883 +/- 0.0472 L/hr for a 40-year-old, 70 kg, nonsmoking female', and 'The bioavailability (FRel) ... was fixed to 1 for a 40-year-old non-black, non-Hispanic female'. An independent regulatory reviewer working from the sponsor's source report read the intercept the same way this extraction does, so the paper's prose is the erroneous side. Female sex distribution of the model-building dataset is 87 of 338 (25.7%). See vignette Errata.",
      source_name        = "Sex"
    ),
    SMOKE = list(
      description        = "Current-smoker indicator, 1 = current smoker at baseline, 0 = non-smoker.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-smoker).",
      notes              = "Linear additive effects: +30.7% on parent CL (theta_15) and +23.5% on N-oxide CL (theta_10). Lahu 2010 Methods page 591: 'Smoking status was defined as current smoking or non-smoking, irrespective of previous smoking status.' Mechanism is increased CYP1A2 activity in smokers (~79% higher per Funck-Brentano 2006; Lahu 2010 Discussion page 602). Smoking on N-oxide Frel was tested and dropped (theta_16 = 0.0129, SE 0.0335; FDA NDA 22-522 review Table 7). Smokers were 109 of 338 (32.2%) of the model-building dataset.",
      source_name        = "Smoking"
    ),
    FED = list(
      description        = "Fed-state indicator at the time of dose, 1 = fed, 0 = fasted.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted).",
      notes              = "Linear additive effects: -30.8% on parent tlag (theta_9; surprising negative sign, see vignette Errata), -69.9% on parent ka (theta_11, food slows absorption), and +236% on N-oxide D1 (theta_6, food prolongs metabolite zero-order input duration). The food effect on tlag had robustness issues in the source bootstrap (only 61% of bootstrap replicates negative; Lahu 2010 page 599), but the paper notes that tlag does not influence steady-state exposure so the finding was not considered critical. The large SE the FDA NDA 22-522 review reports for this coefficient (theta_9 = -0.308, SE 0.296, i.e. the interval spans zero) independently corroborates that it is data-poor. Food on parent Frel (theta_19 = 0.0214, SE 0.0720) and on N-oxide Frel (theta_18 = 0.0133, SE 0.0347) were also tested and dropped.",
      source_name        = "Food"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator, 1 = Black, 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any non-Black race; in Lahu 2010 the reference is 'White (any non-Black or non-Hispanic race)' per Methods page 591).",
      notes              = "Linear additive effects: -14.0% on parent CL (theta_20) and +43.1% on N-oxide Frel (theta_21). Race was not initially in the final parsimonious model selected by WAM; it was added back to test as a covariate on parent and metabolite clearance and on N-oxide Frel after the final model was selected (Lahu 2010 page 592). Race on parent CL had minor bootstrap robustness issues (Lahu 2010 page 599).",
      source_name        = "RaceBlack"
    ),
    RACE_HISPANIC = list(
      description        = "Hispanic / Latino race indicator, 1 = Hispanic, 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any non-Hispanic race; in Lahu 2010 the reference is 'White (any non-Black or non-Hispanic race)' per Methods page 591).",
      notes              = "Linear additive effects: -29.7% on parent CL (theta_21) and +26.7% on N-oxide Frel (theta_22). Hispanic effect on parent CL was estimated with greater uncertainty than the Black effect because of limited Hispanic representation (Lahu 2010 Discussion page 602).",
      source_name        = "RaceHispanic"
    ),
    DIS_COPD = list(
      description        = "Chronic obstructive pulmonary disease patient indicator, 1 = moderate-to-severe COPD patient, 0 = healthy volunteer.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer; reference is the pooled phase I healthy-volunteer cohort).",
      notes              = "Linear additive effects estimated in the phase II / III extension fit: parent CL -39.4% (theta_22) and parent V1 +184% (theta_25); N-oxide CL -7.85% (theta_25 in N-oxide model) and N-oxide Vd -21.4% (theta_26). The COPD covariate effects were fitted holding the phase I-derived covariates and structural parameters fixed; only the COPD coefficients and the residual error were re-estimated (Lahu 2010 Methods page 595, 'Extension to COPD Patients'). Inflammatory-cytokine downregulation of CYP3A4 and CYP1A2 in COPD is the proposed mechanism (Lahu 2010 Discussion page 603).",
      source_name        = "COPD"
    )
  )

  population <- list(
    n_subjects     = 338L,
    n_studies      = 28L,
    age_range      = "Mean 38.2 years, SD 14.7 years in the 338-subject phase I healthy-volunteer model-building dataset; individual minimum and maximum are not reported. The reference age of 40 years used in the covariate equations therefore sits essentially at the cohort mean. Ages of the phase II / III COPD patients contributing the COPD extension are not tabulated in either source (FDA NDA 22-522 Clinical Pharmacology review Table 5, printed page 75).",
    weight_range   = "Mean 75.4 kg, SD 11.2 kg in the 338-subject phase I healthy-volunteer model-building dataset; individual minimum and maximum are not reported. Note that the reference weight of 70 kg used in the N-oxide Vd equation is a round allometric anchor rather than the cohort mean (FDA NDA 22-522 Clinical Pharmacology review Table 5, printed page 75).",
    sex_female_pct = 25.7,
    race_ethnicity = "Pooled White, Black, and Hispanic with White (any non-Black or non-Hispanic race) as the reference category per Lahu 2010 Methods page 591. Phase I model-building dataset composition (n = 338): non-Black / non-Hispanic 296 (87.6%), Black 27 (7.99%), Hispanic 15 (4.44%) (FDA NDA 22-522 Clinical Pharmacology review Table 5, printed page 75).",
    disease_state  = "Pooled adult healthy volunteers (21 phase I index studies + 5 phase I validation studies for roflumilast) and patients with moderate-to-severe COPD (1 phase II [IN-108] + 1 phase III [BY217/M2-110; ClinicalTrials.gov NCT00062582] study used for the COPD extension). The phase III COPD trial enrolled patients with moderate-to-severe COPD treated with roflumilast 500 mcg or placebo orally once daily for 24 weeks.",
    dose_range     = "Oral roflumilast tablets 250-1000 ug in dose-proportionality and dose-escalation phase I studies; standard once-daily dose 500 ug (= 0.5 mg) in the majority of phase I studies and in the phase II / III COPD studies. Note that the Lahu 2010 paper Methods text reports the standard dose as '500 mg', which is a typographical error in the source -- the marketed roflumilast (Daxas / Daliresp) is 500 micrograms once daily and the trial protocol (NCT00062582) and dose-proportionality reference (Bethke 2007) both confirm 250-1000 micrograms. The model file uses ug as the dosing unit and ug/L as the concentration unit so that user-supplied amt values are read directly in micrograms.",
    regions        = "Multinational; specific country composition not reported in the on-disk source.",
    n_observations = "Roflumilast parent: 7705 observations from 338 subjects in the phase I index dataset; an additional 771 observations from 228 subjects contributed to the COPD extension. Roflumilast N-oxide: 7112 observations from 298 subjects in the phase I index dataset; an additional 703 observations from 208 subjects contributed to the COPD extension. Concentrations were quantified by HPLC with tandem mass-spectrometry detection (Lahu 2010 Methods page 591).",
    notes          = "Demographic table S-2 of the Lahu 2010 Supplemental Digital Content is not available on disk for this extraction, but its phase I equivalent was recovered from the FDA NDA 22-522 Clinical Pharmacology review (Table 5, printed page 75), which reproduces the sponsor's own population PK study report 114/2005 -- the same analysis Lahu published. That table gives, for the n = 338 phase I model-building dataset: female 87 (25.7%) / male 251 (74.3%); non-smokers 229 (67.8%) / smokers 109 (32.2%); alcohol non-users 104 (30.8%) / users 234 (69.2%); non-Black non-Hispanic 296 (87.6%) / Black 27 (7.99%) / Hispanic 15 (4.44%); age 38.2 +/- 14.7 years; weight 75.4 +/- 11.2 kg. Every stratum sums to 338. ALCOHOL USE was a tested covariate that the published paper never mentions: the FDA review's full-model columns carry Alcohol on N-oxide CL (theta_11 = 0.00354, SE 0.0626) and on N-oxide Frel (theta_17 = -0.0239, SE 0.0338), both indistinguishable from zero and both dropped by the WAM selection. Alcohol is therefore deliberately absent from covariateData -- it is not a covariate of the shipped model. Twenty-eight studies (21 + 5 + 1 + 1) contributed observations across the joint dataset; the phase II (IN-108) study identifier and the phase III study (BY217/M2-110, NCT00062582) are named in Lahu 2010 Methods. NONMEM version V level 1.1 was used for fitting; the HYBRID method was used (FOCE for all parameters except tlag for the parent and D1 for the metabolite, both estimated by FO)."
  )

  ini({
    # ------------------------------------------------------------------
    # PARENT (roflumilast) STRUCTURAL PARAMETERS
    # ------------------------------------------------------------------
    # Lahu 2010 Table I, 'Final for COPD patients' column: theta values
    # carry over from the phase I 'Final with race' column except for
    # the COPD-specific coefficients (theta_22 on CL, theta_25 on V1).
    # The model-equation intercept (when all binary covariates = 0 and
    # age = 40) is FEMALE, non-smoker, White, healthy, fasted.
    #
    # SECONDARY CORROBORATION. Every non-COPD value below was verified
    # against an independent reproduction of the same analysis: the FDA
    # NDA 22-522 Clinical Pharmacology review, Tables 6 (parent, printed
    # page 77) and 7 (N-oxide, printed page 79), which reprint the
    # sponsor's population PK study report 114/2005 -- Lahu was the
    # sponsor's modeller, so the NDA and the journal article report the
    # same fit. All 40 parent + N-oxide values below reproduce the FDA
    # "Final With Race" columns exactly, with ZERO disagreements. The
    # FDA tables additionally carry standard errors that the paper's
    # Table I omits; they are not transcribed here parameter-by-
    # parameter, but the full digitisation is recorded in the vignette
    # "Source trace" section. The COPD coefficients are the exception:
    # the FDA review's COPD section has no parameter table, so those
    # remain sourced solely from the paper's Table I.

    ltlag <- log(0.158)
    label("Roflumilast absorption lag time at fasted reference (h)")        # Lahu 2010 Table I final-for-COPD: theta_1 = 0.158 h

    lka <- log(0.533)
    label("Roflumilast first-order absorption rate constant at fasted reference (1/h)")  # Lahu 2010 Table I final-for-COPD: theta_2 = 0.533 1/h

    lcl <- log(10.5)
    label("Roflumilast apparent clearance CL/F at female / non-smoker / White / healthy reference (L/h)")  # Lahu 2010 Table I final-for-COPD: theta_3 = 10.5 L/h

    lvc <- log(14.3)
    label("Roflumilast apparent central volume V1/F at healthy reference (L)")  # Lahu 2010 Table I final-for-COPD: theta_4 = 14.3 L

    lq <- log(20.3)
    label("Roflumilast apparent inter-compartmental clearance Q/F (L/h)")    # Lahu 2010 Table I final-for-COPD: theta_5 = 20.3 L/h

    lvp <- log(201)
    label("Roflumilast apparent peripheral volume V2/F (L)")                 # Lahu 2010 Table I final-for-COPD: theta_6 = 201 L

    lfdepot <- fixed(log(1))
    label("Roflumilast bioavailability into depot (F1, 1 -- absolute F not identifiable without IV data)")  # Lahu 2010 Methods page 591: "apparent fraction absorbed for roflumilast ... were unidentifiable"

    # ------------------------------------------------------------------
    # PARENT (roflumilast) COVARIATE EFFECTS
    # ------------------------------------------------------------------

    e_fed_tlag <- -0.308
    label("Linear coefficient for fed-vs-fasted on roflumilast tlag (unitless)")  # Lahu 2010 Table I final-for-COPD: theta_9 = -0.308 (food effect with poor bootstrap robustness; see vignette Errata)

    e_fed_ka <- -0.699
    label("Linear coefficient for fed-vs-fasted on roflumilast ka (unitless)")    # Lahu 2010 Table I final-for-COPD: theta_11 = -0.699 (food slows parent absorption)

    e_sex_m_cl <- 0.191
    label("Linear coefficient for male-vs-female on roflumilast CL (unitless)")   # Lahu 2010 Table I final-for-COPD: theta_13 = 0.191 (male sex; canonical SEXF inverts via (1 - SEXF))

    e_smoke_cl <- 0.307
    label("Linear coefficient for current-smoker-vs-non-smoker on roflumilast CL (unitless)")  # Lahu 2010 Table I final-for-COPD: theta_15 = 0.307

    e_race_black_cl <- -0.140
    label("Linear coefficient for Black-vs-non-Black on roflumilast CL (unitless)")   # Lahu 2010 Table I final-for-COPD: theta_20 = -0.140

    e_race_hisp_cl <- -0.297
    label("Linear coefficient for Hispanic-vs-non-Hispanic on roflumilast CL (unitless)")  # Lahu 2010 Table I final-for-COPD: theta_21 = -0.297

    e_copd_cl <- -0.394
    label("Linear coefficient for COPD-vs-healthy on roflumilast CL (unitless)")  # Lahu 2010 Table I final-for-COPD: theta_22 = -0.394 (CL reduced by 39.4% in COPD)

    e_copd_vc <- 1.84
    label("Linear coefficient for COPD-vs-healthy on roflumilast V1 (unitless)")  # Lahu 2010 Table I final-for-COPD: theta_25 = 1.84 (V1 increased by 184% in COPD)

    # ------------------------------------------------------------------
    # PARENT (roflumilast) INTER-INDIVIDUAL VARIABILITY
    # ------------------------------------------------------------------
    # Lahu 2010 reports IIV as variance (omega^2) directly in Table I
    # footnote ("xi = variance"). Parent has a single off-diagonal
    # covariance between Q and V2 (no other off-diagonals are listed in
    # the table). IIV(tlag) for the parent is 1.73, which corresponds
    # to SD = 1.32 (132% lognormal CV); the source notes that tlag
    # could only be fitted with the FO estimation method and had poor
    # bootstrap robustness (see vignette Errata).

    etaltlag ~ 1.73                                                          # Lahu 2010 Table I final-for-COPD: omega^2(eta_tlag) = 1.73
    etalka   ~ 0.154                                                         # Lahu 2010 Table I final-for-COPD: omega^2(eta_ka) = 0.154
    etalcl   ~ 0.136                                                         # Lahu 2010 Table I final-for-COPD: omega^2(eta_CL) = 0.136
    etalvc   ~ 0.734                                                         # Lahu 2010 Table I final-for-COPD: omega^2(eta_V1) = 0.734

    etalq + etalvp ~ c(0.0726,
                       0.0703, 0.117)                                        # Lahu 2010 Table I final-for-COPD: omega^2(eta_Q) = 0.0726, omega(eta_Q,eta_V2) = 0.0703, omega^2(eta_V2) = 0.117

    # ------------------------------------------------------------------
    # PARENT (roflumilast) RESIDUAL ERROR
    # ------------------------------------------------------------------
    # Phase I 'Final with race' residual = 25.1% CV. The 'Final for
    # COPD patients' column re-fits the residual on the smaller phase II /
    # III COPD-only dataset and reports a much larger 54.5% CV, which
    # primarily reflects sparse-sampling design rather than a different
    # intrinsic analytical noise mechanism. The model file carries the
    # phase I (richer-data) residual; the COPD-specific value is
    # discussed in the vignette Errata.

    propSd <- 0.251
    label("Roflumilast proportional residual SD on linear concentration (fraction)")  # Lahu 2010 Table I 'Final with race': theta_7 CV{sigma} = 25.1%; equivalent to log-additive SD on log-Cc

    # ------------------------------------------------------------------
    # METABOLITE (roflumilast N-oxide) STRUCTURAL PARAMETERS
    # ------------------------------------------------------------------

    ltlag_noxide <- log(0.156)
    label("Roflumilast N-oxide absorption lag time (h)")                    # Lahu 2010 Table I final-for-COPD: theta_1 = 0.156 h (N-oxide model)

    ld1_noxide <- log(2.21)
    label("Roflumilast N-oxide zero-order input duration D1 at fasted reference (h)")  # Lahu 2010 Table I final-for-COPD: theta_2 = 2.21 h (N-oxide model)

    lcl_noxide <- log(0.883)
    label("Roflumilast N-oxide apparent clearance CL/F at age 40 / female / non-smoker / healthy reference (L/h)")  # Lahu 2010 Table I final-for-COPD: theta_3 = 0.883 L/h (N-oxide model)

    lvd_noxide <- log(65.8)
    label("Roflumilast N-oxide apparent volume of distribution Vd/F at WT 70 kg / healthy reference (L)")  # Lahu 2010 Table I final-for-COPD: theta_4 = 65.8 L (N-oxide model)

    # ------------------------------------------------------------------
    # METABOLITE (roflumilast N-oxide) COVARIATE EFFECTS
    # ------------------------------------------------------------------

    e_fed_d1_noxide <- 2.36
    label("Linear coefficient for fed-vs-fasted on N-oxide D1 (unitless)")  # Lahu 2010 Table I final-for-COPD: theta_6 = 2.36 (food prolongs metabolite zero-order input)

    e_age_cl_noxide <- -0.471
    label("Power exponent of (AGE/40) on N-oxide CL (unitless)")            # Lahu 2010 Table I final-for-COPD: theta_7 = -0.471

    e_sex_m_cl_noxide <- 0.467
    label("Linear coefficient for male-vs-female on N-oxide CL (unitless)")  # Lahu 2010 Table I final-for-COPD: theta_8 = 0.467 (canonical SEXF inverts via (1 - SEXF))

    e_smoke_cl_noxide <- 0.235
    label("Linear coefficient for current-smoker-vs-non-smoker on N-oxide CL (unitless)")  # Lahu 2010 Table I final-for-COPD: theta_10 = 0.235

    e_wt_vd_noxide <- 1.00
    label("Power exponent of (WT/70) on N-oxide Vd (unitless)")             # Lahu 2010 Table I final-for-COPD: theta_12 = 1.00

    e_copd_cl_noxide <- -0.0785
    label("Linear coefficient for COPD-vs-healthy on N-oxide CL (unitless)")  # Lahu 2010 Table I final-for-COPD: theta_25 = -0.0785 (CL reduced by 7.85% in COPD; N-oxide model uses theta_25 distinct from parent's theta_25 on V1)

    e_copd_vd_noxide <- -0.214
    label("Linear coefficient for COPD-vs-healthy on N-oxide Vd (unitless)")  # Lahu 2010 Table I final-for-COPD: theta_26 = -0.214 (Vd reduced by 21.4% in COPD)

    # ------------------------------------------------------------------
    # METABOLITE Frel COVARIATE EFFECTS
    # ------------------------------------------------------------------
    # Frel baseline = 1 (fixed because absolute F is non-identifiable
    # without IV roflumilast or directly-administered N-oxide data;
    # Lahu 2010 Methods page 591). Covariate effects scale Frel
    # multiplicatively but the baseline itself is not a free parameter.

    e_age_frel <- -0.269
    label("Power exponent of (AGE/40) on N-oxide Frel (unitless)")          # Lahu 2010 Table I final-for-COPD: theta_13 = -0.269

    e_sex_m_frel <- 0.231
    label("Linear coefficient for male-vs-female on N-oxide Frel (unitless)")  # Lahu 2010 Table I final-for-COPD: theta_14 = 0.231

    e_race_black_frel <- 0.431
    label("Linear coefficient for Black-vs-non-Black on N-oxide Frel (unitless)")  # Lahu 2010 Table I final-for-COPD: theta_21 = 0.431 (N-oxide model)

    e_race_hisp_frel <- 0.267
    label("Linear coefficient for Hispanic-vs-non-Hispanic on N-oxide Frel (unitless)")  # Lahu 2010 Table I final-for-COPD: theta_22 = 0.267 (N-oxide model)

    # ------------------------------------------------------------------
    # METABOLITE (roflumilast N-oxide) INTER-INDIVIDUAL VARIABILITY
    # ------------------------------------------------------------------
    # Lahu 2010 Table I reports a 3x3 covariance block on D1, CL, Vd
    # with three off-diagonals. No IIV is reported on N-oxide tlag
    # (the paper's eq 7 shows tlag = theta * exp(eta_tlag) but Table I
    # does not list an omega^2 for the N-oxide tlag; see vignette
    # Errata). The omission is REAL, not a typesetting loss in the
    # journal: the FDA NDA 22-522 review's Table 7 (printed page 79),
    # reproducing the sponsor's own study report 114/2005, likewise
    # lists an omega block over D1, CL and V only, across all four of
    # its Base / Full / Initial Final / Final With Race columns.

    etald1_noxide + etalcl_noxide + etalvd_noxide ~
      c(0.268,
        0.0221,  0.150,
        0.0536, -0.0110, 0.0449)                                             # Lahu 2010 Table I final-for-COPD: omega^2(D1)=0.268, omega(D1,CL)=0.0221, omega^2(CL)=0.150, omega(D1,Vd)=0.0536, omega(CL,Vd)=-0.0110, omega^2(Vd)=0.0449

    # ------------------------------------------------------------------
    # METABOLITE (roflumilast N-oxide) RESIDUAL ERROR
    # ------------------------------------------------------------------
    # Phase I 'Final with race' residual = 24.1% CV. The 'Final for
    # COPD patients' column reports a refit value of 20.9% CV from
    # the COPD-only data (lower than the phase I value, in contrast
    # to the parent which increased to 54.5%). The model file carries
    # the phase I residual.

    propSd_noxide <- 0.241
    label("Roflumilast N-oxide proportional residual SD on linear concentration (fraction)")  # Lahu 2010 Table I 'Final with race': theta_5 CV{sigma} = 24.1%; equivalent to log-additive SD on log-Cc_noxide
  })

  model({
    # Reference covariate values for the power / linear effect equations
    # (Lahu 2010 equations 6 and 7).
    ref_age <- 40
    ref_wt  <- 70

    # ------------------------------------------------------------------
    # PARENT (roflumilast) individual parameters -- Lahu 2010 equation 6
    # ------------------------------------------------------------------
    # SEXF inverts the paper's Sex = 1 = male convention; the male-
    # indicator effect is applied via (1 - SEXF) so the model-equation
    # intercept (when SEXF = 1) is the female reference.

    tlag_parent <- exp(ltlag + etaltlag) * (1 + e_fed_tlag * FED)
    ka          <- exp(lka + etalka) * (1 + e_fed_ka * FED)
    cl          <- exp(lcl + etalcl) *
                   (1 + e_sex_m_cl * (1 - SEXF)) *
                   (1 + e_smoke_cl * SMOKE) *
                   (1 + e_race_black_cl * RACE_BLACK) *
                   (1 + e_race_hisp_cl  * RACE_HISPANIC) *
                   (1 + e_copd_cl * DIS_COPD)
    vc          <- exp(lvc + etalvc) * (1 + e_copd_vc * DIS_COPD)
    q           <- exp(lq  + etalq)
    vp          <- exp(lvp + etalvp)

    # ------------------------------------------------------------------
    # METABOLITE (roflumilast N-oxide) individual parameters -- Lahu 2010 equation 7
    # ------------------------------------------------------------------

    tlag_noxide <- exp(ltlag_noxide)
    d1_noxide   <- exp(ld1_noxide + etald1_noxide) *
                   (1 + e_fed_d1_noxide * FED)
    cl_noxide   <- exp(lcl_noxide + etalcl_noxide) *
                   (AGE / ref_age)^e_age_cl_noxide *
                   (1 + e_sex_m_cl_noxide * (1 - SEXF)) *
                   (1 + e_smoke_cl_noxide * SMOKE) *
                   (1 + e_copd_cl_noxide * DIS_COPD)
    vd_noxide   <- exp(lvd_noxide + etalvd_noxide) *
                   (WT / ref_wt)^e_wt_vd_noxide *
                   (1 + e_copd_vd_noxide * DIS_COPD)
    frel        <- (AGE / ref_age)^e_age_frel *
                   (1 + e_sex_m_frel * (1 - SEXF)) *
                   (1 + e_race_black_frel * RACE_BLACK) *
                   (1 + e_race_hisp_frel  * RACE_HISPANIC)

    # Micro-constants (parent two-compartment + metabolite one-compartment).
    kel_parent <- cl        / vc
    k12        <- q         / vc
    k21        <- q         / vp
    kel_noxide <- cl_noxide / vd_noxide

    # ODE system. The parent has standard first-order absorption from a
    # depot compartment; the N-oxide receives zero-order input directly
    # into its central compartment (rxode2 dur() / f() / lag() set the
    # zero-order duration, relative bioavailability, and lag time when
    # the user dosing event targets cmt = "central_noxide"). Both depot
    # and central_noxide receive an event of amt = dose-of-roflumilast
    # at the same time; F1 = 1 (parent) and F2 = Frel (N-oxide) handle
    # the relative scaling.
    d/dt(depot)          <- -ka * depot
    d/dt(central)        <-  ka * depot - kel_parent * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1)    <-  k12 * central - k21 * peripheral1
    d/dt(central_noxide) <- -kel_noxide * central_noxide

    # Absorption modifiers.
    lag(depot)          <- tlag_parent
    f(depot)            <- exp(lfdepot)

    lag(central_noxide) <- tlag_noxide
    dur(central_noxide) <- d1_noxide
    f(central_noxide)   <- frel

    # Observations. Dose is in ug, volumes in L, so central / vc gives
    # ug/L directly -- the unit Lahu 2010 uses for both analytes (see
    # Figure 1 axis labels).
    Cc        <- central        / vc
    Cc_noxide <- central_noxide / vd_noxide

    Cc        ~ prop(propSd)
    Cc_noxide ~ prop(propSd_noxide)
  })
}
