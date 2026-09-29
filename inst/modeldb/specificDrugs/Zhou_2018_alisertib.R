Zhou_2018_alisertib <- function() {
  description <- paste0(
    "Global two-compartment population PK model for the investigational ",
    "Aurora A kinase inhibitor alisertib (MLN8237) given orally to adult ",
    "cancer patients, with a FOUR-TRANSIT-COMPARTMENT absorption chain and ",
    "linear elimination (Zhou 2018, n = 422 in the model-development set ",
    "pooled from 10 phase I / I-II / II studies in Western countries, Japan ",
    "and other East Asian countries; a further 249 patients formed an ",
    "external validation set). The dose enters transit1 and passes through ",
    "four sequential first-order transfers at the common rate ktr, so the ",
    "net mean oral transit time is 4/ktr = 0.96 h. The one clinically ",
    "important covariate is REGION: patients enrolled in East Asia were ",
    "estimated to have 52 percent higher relative bioavailability than ",
    "Western patients, encoded exactly as Zhou 2018 Figure 1 writes it -- a ",
    "common multiplier (1 + e_region_eastasia_f) applied to ALL FOUR ",
    "apparent disposition parameters CL/F, V1/F, Q/F and V2/F, which is what ",
    "a change in F alone does to an oral model. Body surface area acts on ",
    "V1/F only and therefore changes peak-to-trough fluctuation without ",
    "changing AUC; body weight, BSA on CL/F, UGT1A1 *28 / *6 genotype, ",
    "creatinine clearance, sex, age, race, ALT, AST and bilirubin were all ",
    "screened and none explained the apparent oral clearance. Between-subject ",
    "variability is a full 4x4 correlated block on CL/F, V1/F, V2/F and KTR; ",
    "Q/F deliberately carries none. Concentrations are MOLAR (nmol/L), so ",
    "doses must be supplied in nmol -- see the vignette for the mg-to-nmol ",
    "conversion, and the dose must be given into transit1. PASS ",
    "useLinCmt = FALSE TO EVERY rxSolve() CALL: the disposition is written ",
    "with explicit k12 / k21 micro-constants, and rxSolve's default ",
    "ODE-to-linCmt auto-conversion silently discards peripheral1, which ",
    "leaves AUC unchanged and so is invisible to any exposure-only check -- ",
    "only the terminal half-life moves. Companion exposure-safety logistic models for the three ",
    "mechanism-related antiproliferative toxicities are packaged as ",
    "Zhou_2018_alisertib_neutropenia, Zhou_2018_alisertib_stomatitis and ",
    "Zhou_2018_alisertib_diarrhea."
  )
  reference <- paste(
    "Zhou X, Mould DR, Takubo T, Sheldon-Waniga E, Huebner D, Milton A,",
    "Venkatakrishnan K. Global population pharmacokinetics of the",
    "investigational Aurora A kinase inhibitor alisertib in cancer patients:",
    "rationale for lower dosage in Asia.",
    "Br J Clin Pharmacol. 2018;84(1):35-51. doi:10.1111/bcp.13430.",
    "Structural and covariate equations from Figure 1; final parameter",
    "estimates from Table 5; the updated final model refitted on the combined",
    "analysis plus validation data is Supplementary Table S1 and is NOT the",
    "model packaged here.",
    sep = " "
  )
  vignette <- "Zhou_2018_alisertib"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")
  # The oral dose enters transit1, not a state named depot or central, so
  # buildModelDb()'s two-name heuristic cannot infer it and the registry would
  # otherwise report the dosing compartment as `central`.
  dosing <- "transit1"

  covariateData <- list(
    BSA = list(
      description = "Baseline body surface area.",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power-scaled on the apparent central volume V1/F only, referenced to",
        "1.84 m^2 -- the median BSA of both the combined (n = 671) and the",
        "analysis (n = 422) data sets (Table 2). The BSA computation formula is",
        "not stated in the source. Zhou 2018 found BSA significant on V1/F but",
        "NOT on CL/F, so BSA changes Cmax and peak-to-trough fluctuation while",
        "leaving AUC untouched; the paper leans on this to argue that the",
        "region effect on bioavailability is not an artefact of the smaller",
        "body size of the East Asian cohort. Figure 4 simulates the 2.5th, 50th",
        "and 97.5th BSA percentiles as 1.44, 1.84 and 2.43 m^2. The two",
        "simulated regional cohorts have geometric mean BSA 1.88 m^2 (West, log",
        "SD 0.135) and 1.63 m^2 (East Asia, log SD 0.0862).",
        "Analysis-set range 1.36 to 2.97 m^2; combined-set range 1.34 to 3.28."
      ),
      source_name = "BSA"
    ),
    REGION_EASTASIA = list(
      description = "East Asian enrolling-region indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Western countries)",
      notes = paste(
        "1 = enrolled in Japan or in the other East Asian study (Singapore,",
        "Taiwan, Hong Kong, South Korea); 0 = enrolled in a Western country.",
        "59 of 422 patients (14 percent) in the analysis set and 59 of 671 (9",
        "percent) in the combined set; the external validation set is 100",
        "percent Western.",
        "This is REGION and NOT race. Zhou 2018 separates the two deliberately:",
        "the 8 patients of Asian race enrolled in the West had dose-normalized",
        "AUC of 475 nM.h/mg, matching the 482 of Western non-Asian patients and",
        "not the 797 of Asian patients in East Asia (Table 6, Figure 5), so the",
        "authors attribute the effect to extrinsic regional factors rather than",
        "to ancestry. With only 8 such patients the paper explicitly declines to",
        "call the race-versus-region attribution conclusive. Race was screened",
        "separately and not retained; see covariatesDataExcluded.",
        "The coefficient is applied as a shared multiplier on CL/F, V1/F, Q/F",
        "and V2/F, which is algebraically identical to raising F."
      ),
      source_name = "RGN"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight.",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened on CL/F, V1/F and V2/F (Table 4) and not retained. Median 73.3",
        "kg (range 42.6 to 175.0) in the analysis set. The Discussion states",
        "body weight 'was not important in explaining the variability in the",
        "apparent oral clearance of alisertib', supporting fixed rather than",
        "body-size-adjusted dosing."
      ),
      source_name = "WT"
    ),
    CRCL = list(
      description = "Baseline creatinine clearance.",
      units = "mL/min",
      type = "continuous",
      notes = paste(
        "Screened on CL/F (Table 4) and not retained. Median 83.2 mL/min (range",
        "27.1 to 241.0) in the analysis set. Zhou 2018 reads the null result as",
        "evidence that mild or moderate renal impairment (down to 27 mL/min) has",
        "no clinically meaningful effect on alisertib exposure."
      ),
      source_name = "CCL"
    ),
    ALB = list(
      description = "Baseline plasma albumin.",
      units = "g/L",
      type = "continuous",
      notes = "Screened on CL/F, V1/F and V2/F (Table 4) and not retained. Median 40 g/L (range 20.0 to 51.4).",
      source_name = "ALB"
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL/F (Table 4) and not retained. Median 22 U/L (range 5.0 to 229.0).",
      source_name = "ALT"
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL/F (Table 4) and not retained. Median 26 U/L (range 9.0 to 341.0).",
      source_name = "AST"
    ),
    BILI = list(
      description = "Baseline total bilirubin.",
      units = "umol/L",
      type = "continuous",
      notes = "Screened on CL/F (Table 4) and not retained. Median 7.0 umol/L (range 1.71 to 38.0 in the analysis set).",
      source_name = "BILI"
    ),
    AGE = list(
      description = "Baseline age.",
      units = "year",
      type = "continuous",
      notes = "Screened on CL/F (Table 4) and not retained. Median 62 years (range 21 to 85).",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL/F, V1/F and V2/F (Table 4) and not retained. 196 of 422 (46 percent) female in the analysis set.",
      source_name = "Gender"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened on CL/F and on F (Table 4) and NOT retained -- the retained",
        "covariate is the enrolling region REGION_EASTASIA, which is a distinct",
        "concept. See the REGION_EASTASIA notes for the Table 6 / Figure 5",
        "evidence the authors used to separate the two."
      ),
      source_name = "Race"
    ),
    UGT1A1_28_N = list(
      description = "Number of UGT1A1*28 variant alleles carried (0, 1 or 2).",
      units = "alleles",
      type = "count",
      notes = paste(
        "Screened on CL/F (Table 4) and not retained; 203 / 138 / 37 patients in",
        "the analysis set carried 0 / 1 / 2 copies and 44 were not genotyped.",
        "Alisertib is glucuronidated by several UGT isoforms and the *28 variant",
        "only partially reduces UGT1A1 expression, which Zhou 2018 offers as the",
        "explanation for the null effect. Reported here for provenance only; no",
        "canonical count column is introduced, because the register requires a",
        "count covariate to be decomposed into binary indicators and this one is",
        "not used by the model."
      ),
      source_name = "UGT1A1*28 alleles, n"
    ),
    UGT1A1_6_N = list(
      description = "Number of UGT1A1*6 variant alleles carried (0, 1 or 2).",
      units = "alleles",
      type = "count",
      notes = paste(
        "Screened on CL/F (Table 4) and not retained; genotyped only in the",
        "Japanese and other East Asian studies, where 36 / 11 / 3 patients",
        "carried 0 / 1 / 2 copies. Combined *6-plus-*28 allele-count models were",
        "also tested and were likewise not significant."
      ),
      source_name = "UGT1A1*6 alleles, n"
    ),
    FORM_ALISERTIB_ECT = list(
      description = "Enteric-coated-tablet formulation indicator (versus powder-in-capsule).",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened on KTR / lag and on F (Table 4) and not retained. 123 of 422",
        "(29 percent) of the analysis set received the enteric-coated tablet and",
        "299 (71 percent) the powder-in-capsule; the external validation set is",
        "100 percent enteric-coated tablet."
      ),
      source_name = "Alisertib formulation"
    )
  )

  compartmentData <- list(
    transit1 = list(analyte = "alisertib", units = "nmol", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "alisertib", units = "nmol", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "alisertib", units = "nmol", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "alisertib", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "alisertib", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "alisertib", units = "nmol", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 422L,
    n_studies = 10L,
    age_range = "21 to 85 years (analysis set); 21 to 88 years across the combined 671-patient data set",
    age_median = "62 years",
    weight_range = "42.6 to 175.0 kg (analysis set)",
    weight_median = "73.3 kg",
    bsa_range = "1.36 to 2.97 m^2 (analysis set)",
    bsa_median = "1.84 m^2",
    sex_female_pct = 46,
    race_ethnicity = c(White = 77, `Asian (all)` = 14, Black = 5, Other = 1, Missing = 2),
    disease_state = "Advanced haematological (43 percent) and nonhaematological (57 percent) malignancies",
    dose_range = "5 to 150 mg alisertib orally once or twice daily, under nil-per-os conditions (no food from 2 h before until 1 h after each dose), in 7-, 14- or 21-day dosing schedules within 21-, 28- or 35-day cycles; the recommended single-agent schedule is 7 days of twice-daily dosing in a 21-day cycle",
    regions = "Western countries (7 studies, 363 of 422 analysis-set patients); Japan (2 studies) and Singapore / Taiwan / Hong Kong / South Korea (1 study), together 59 of 422",
    renal_function = "Creatinine clearance 27.1 to 241.0 mL/min (analysis set); no effect on CL/F was detected over this range",
    notes = paste0(
      "Baseline demographics are Table 2 (continuous) and Table 3 ",
      "(categorical) of Zhou 2018; the study list is Table 1. The 671 pooled ",
      "patients were split into a 422-patient model-development set (rich plus ",
      "sparse sampling) and a 249-patient external validation set (the phase II ",
      "portion of study C14007, entirely Western and entirely enteric-coated ",
      "tablet). The parameters packaged here are the FINAL MODEL fitted to the ",
      "422-patient analysis set (Table 5). Post-hoc CL/F in the validation set ",
      "ran 10.9 percent low with 23.7 percent precision, which the authors ",
      "judged unimportant against the 52 percent between-subject variability. ",
      "The maximum tolerated dose is 50 mg twice daily in the West and 30 mg ",
      "twice daily in East Asia."
    )
  )

  ini({
    # ================================================================
    # Zhou 2018 Table 5, 'Final model parameters' -- the model fitted to
    # the 422-patient analysis data set. Estimated by SAEM in NONMEM
    # 7.2 with a log-transform-both-sides approach; 95 percent CIs come
    # from a 1000-replicate nonparametric bootstrap, all of which
    # converged.
    #
    # Do NOT substitute Supplementary Table S1: that is the 'updated
    # final model' refitted after the external validation set was
    # folded in, and it differs on every row (CL/F 3.91, V1/F 51.8,
    # Q/F 7.53, V2/F 27.4, KTR 3.93, RGNF -0.325, BSAV1 0.79,
    # CCV 0.477, ADD 4.8 nmol/L). The two are discussed in the vignette.
    #
    # Scale of the 'IIV, ratio' column: Table 5 prints 0.518 for CL/F
    # and the Results text reads it back as 'interpatient coefficient
    # of variation (CV): 51.8%', i.e. CV percent = 100 x the tabulated
    # ratio. That is the usual CV ~= omega approximation, so the column
    # is the LOG-SCALE SD and the variances below are its square. The
    # reading is confirmed against the paper's own simulation: Figure
    # 4's inset reports geometric-mean Cmax,ss / Cmin,ss and their CVs
    # for 200 Western patients at three BSA percentiles, and the
    # log-SD reading reproduces all twelve numbers substantially better
    # than treating the column as an exact lognormal CV. The vignette
    # re-runs that comparison.
    # ================================================================

    # ----- Structural disposition, apparent (oral) parameters -----
    lcl <- log(4.11); label("Apparent oral clearance CL/F (L/h)")  # Table 5, CL/F 4.11, SE 3.1 percent
    lvc <- log(54.3); label("Apparent central volume V1/F (L)")  # Table 5, V1/F 54.3, SE 3.9 percent
    lq <- log(7.07); label("Apparent intercompartmental clearance Q/F (L/h)")  # Table 5, Q/F 7.07, SE 10.9 percent
    lvp <- log(28.7); label("Apparent peripheral volume V2/F (L)")  # Table 5, V2/F 28.7, SE 10 percent

    # ----- Transit absorption -----
    # Four transit compartments in series, all sharing the single rate
    # ktr (Figure 1). Net mean oral transit time = 4/ktr = 0.959 h,
    # matching the 0.96 h quoted in Results.
    lktr <- log(4.17); label("Transit-compartment transfer rate constant KTR (1/h)")  # Table 5, KTR 4.17, SE 3 percent

    # ----- Covariate effects -----
    # Figure 1: every apparent parameter is multiplied by (1 + f), with
    # f = 0 in the West and f = RGNF in East Asia. Table 5 prints the
    # magnitude 0.341 with a footnote reading 'the affected parameters
    # were multiplied by (1-0.341) for the East Asian region'; the
    # supplement prints the same coefficient WITH its minus sign
    # (RGNF = -0.325 in the updated final model), which settles the
    # sign that main-text typesetting dropped. The multiplier is
    # therefore 1 - 0.341 = 0.659, and 1/0.659 = 1.517 reproduces the
    # '51.7 percent (95 percent CI 36.5, 70.1) higher F' quoted in
    # Results.
    e_region_eastasia_f <- -0.341; label("Fractional change in CL/F, V1/F, Q/F and V2/F for enrolment in East Asia versus the West (unitless); a negative value is higher relative bioavailability")  # Table 5, RGNF 0.341, SE 11 percent, sign from the Table 5 footnote and Supplementary Table S1
    e_bsa_vc <- 0.899; label("Power exponent for body surface area on the apparent central volume V1/F (unitless)")  # Table 5, BSAV1 0.899, SE 26.3 percent

    # ----- Between-subject variability -----
    # Full 4x4 block on CL/F, V1/F, V2/F and KTR (Results). Q/F carries
    # none: Figure 1 writes Q/F without an exponential term and Table 5
    # leaves its IIV row blank, because Methods says the effect of
    # removing BSV from Q/F was examined and it was dropped.
    # Diagonals are the squares of the Table 5 'IIV, ratio' column;
    # off-diagonals are r * sd_i * sd_j from the Table 5 correlation
    # matrix, which is lower-triangular in the order CL/F, V1/F, V2/F,
    # KTR.
    etalcl + etalvc + etalvp + etalktr ~ c(
      0.268324,
      0.124208, 0.169744,
      0.275263, 0.036561, 1.089936,
      0.032448, 0.056732, 0.014094, 0.291600
    )  # Table 5: log-SD 0.518, 0.412, 1.044, 0.54 (SE 10.3, 15.8, 11.6, 17.5 percent); correlations 0.582, 0.509, 0.085, 0.116, 0.255, 0.025

    # ----- Residual error -----
    # Methods declare a combined proportional-plus-additive model. The
    # additive term went to zero, so the fitted residual is purely
    # proportional; it is kept here at its reported value of zero so
    # the declared structure is not silently lost.
    propSd <- 0.491; label("Proportional residual error (fraction)")  # Table 5, CCV 0.491, SE 3.5 percent
    addSd <- fixed(0); label("Additive residual error (nmol/L); driven to zero by the SAEM algorithm, per the Table 5 footnote")  # Table 5, ADD 0
  })

  model({
    # ---- Regional multiplier on the apparent parameters --------------
    # Zhou 2018 Figure 1. Raising relative bioavailability by a factor
    # 1/(1 + e_region_eastasia_f) divides every APPARENT parameter by
    # that same factor, which is why one coefficient moves CL/F, V1/F,
    # Q/F and V2/F together. Concentrations are unchanged by any other
    # reading of the same effect.
    regmult <- 1 + e_region_eastasia_f * REGION_EASTASIA

    # ---- Individual parameters ---------------------------------------
    cl <- exp(lcl + etalcl) * regmult
    vc <- exp(lvc + etalvc) * regmult * (BSA / 1.84)^e_bsa_vc
    q <- exp(lq) * regmult
    vp <- exp(lvp + etalvp) * regmult
    ktr <- exp(lktr + etalktr)

    # ---- Micro-constants ---------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- ODE system ---------------------------------------------------
    # Dose enters transit1; four sequential transfers at rate ktr
    # deliver drug to central (Figure 1 schematic).
    d/dt(transit1) <- -ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(central) <- ktr * transit4 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- Observation ---------------------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
