Comisar_2025_zavegepant <- function() {
  description <- paste(
    "Three-compartment population PK model for zavegepant (ZAVZPRET, a",
    "small-molecule calcitonin gene-related peptide (CGRP) receptor",
    "antagonist approved for acute treatment of migraine) in healthy",
    "adults and patients with migraine (Comisar 2025; N = 277 subjects,",
    "10819 concentrations pooled from nine phase I studies, with a tenth",
    "study used for external validation). Disposition is a",
    "three-compartment model with first-order elimination from the",
    "central compartment; absorption after intranasal or oral dosing is",
    "sequential zero-order then first-order (the dose enters the depot",
    "over a zero-order window D1 and the depot then drains first-order at",
    "ka), while intravenous doses enter the central compartment directly.",
    "Typical values for a 70 kg subject are CL = 13.3 L/h, Vc = 12.1 L,",
    "Q = 2.6 L/h and Vp = 66.0 L for the first (deep) peripheral",
    "compartment, and Q2 = 5.2 L/h and Vp2 = 10.7 L for the second",
    "(shallow) peripheral compartment, giving Vss = 88.8 L. Absorption is",
    "route-specific: intranasal F = 5.1%, ka = 5.8 1/h, D1 = 8.6 min and",
    "no lag; oral F = 0.65%, ka = 0.81 1/h, D1 = 57.2 min and a 12.2 min",
    "lag. Body-weight allometry uses the standard fixed exponents",
    "(0.75 on all clearances, 1 on all volumes) referenced to 70 kg.",
    "Three covariates act on CL as fractional shifts: moderate hepatic",
    "impairment (Child-Pugh 7-9) -43.9%, rifampin co-administration",
    "-41.1%, and itraconazole co-administration -27.9% but only for",
    "orally administered zavegepant. Fed state lowers oral bioavailability",
    "by 50.9%. Age, sex, race, ethnicity, migraine status, creatinine",
    "clearance, oral contraceptives and sumatriptan were screened and not",
    "retained. Inter-individual variability is estimated on F (66.1%),",
    "ka (44.0%), D1 (115%), Vc (32.2%), Vp (26.1%) and CL (26.8%), each",
    "shared across routes; residual error is proportional (37.8%) plus a",
    "negligible fixed additive term.",
    sep = " "
  )
  reference <- paste(
    "Comisar CM, Francis J, Hughes JH, Bhardwaj R, Bertz R, Liu J. (2025).",
    "Population pharmacokinetic modeling of zavegepant, a calcitonin",
    "gene-related peptide receptor antagonist, in healthy adults and",
    "patients with migraine.",
    "CPT Pharmacometrics Syst Pharmacol 14(1):179-191.",
    "doi:10.1002/psp4.13257",
    sep = " "
  )
  vignette <- "Comisar_2025_zavegepant"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Comisar 2025 Figure 1 (schematic of
  # the selected model) and the Results 'Population pharmacokinetic model
  # development' paragraph.
  compartmentData <- list(
    depot       = list(analyte = "zavegepant", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "zavegepant", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "zavegepant", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "zavegepant", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline, time-fixed. Population median 74 kg (range 49-131;",
        "Comisar 2025 Table 2 'Body weight (kg)'). The allometric",
        "reference weight is 70 kg, NOT the cohort median -- Comisar 2025",
        "Methods 'Population pharmacokinetic modeling' states 'Body",
        "weight-based empirical allometric scaling was applied using",
        "standard exponents (0.75 for clearance and 1 for volume of",
        "distribution) and a reference body weight of 70 kg', and Table 3",
        "Note confirms the tabulated values are 'for a typical body weight",
        "of 70 kg'. Applied to every disposition parameter: CL, Q, Q2 with",
        "exponent 0.75 and Vc, Vp, Vp2 with exponent 1. Both exponents are",
        "fixed, not estimated (Table 3 Note: 'body weight exponents fixed",
        "to 1.0 and 0.75'). Comisar 2025 Discussion reports the resulting",
        "population mean ratios versus the 70 kg typical subject: CL 85.7%",
        "at 57 kg and 126% at 95 kg; volume 81.4% at 57 kg and 136% at",
        "95 kg -- all four reproduce exactly from (WT/70)^0.75 and",
        "(WT/70)^1."
      ),
      source_name        = "WT"
    ),
    ROUTE_ORAL = list(
      description        = "Oral (soft gelatin capsule) administration indicator on the dose record",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-oral: intranasal nasal spray, or intravenous infusion)",
      notes              = paste(
        "Per-dose-record indicator selecting the oral absorption parameter",
        "set. Comisar 2025 Table 3 reports a separate bioavailability,",
        "first-order absorption rate constant and zero-order absorption",
        "duration for intranasal and for oral zavegepant, and an",
        "absorption lag time for oral only. ROUTE_ORAL = 1 selects",
        "F = 0.65%, ka = 0.81 1/h, D1 = 57.2 min and tlag = 12.2 min;",
        "ROUTE_ORAL = 0 selects the intranasal set F = 5.1%,",
        "ka = 5.8 1/h, D1 = 8.6 min and no lag. It additionally gates the",
        "itraconazole clearance effect, which Table 3 restricts to",
        "'clearance when zavegepant administered orally'.",
        "Intravenous doses are distinguished from intranasal ones not by",
        "this column but by the dose record's cmt: IV doses are placed",
        "directly into `central` and therefore never read the depot",
        "parameters at all. Set ROUTE_ORAL = 0 on IV records. In the",
        "model-development cohort 62/277 subjects (22.4%) received oral",
        "zavegepant, 209 (75.5%) intranasal and 6 (2.2%) intravenous",
        "(Comisar 2025 Results 'Baseline population characteristics')."
      ),
      source_name        = "Administration route"
    ),
    FED = list(
      description        = "Fed-versus-fasted state at the time of the oral dose",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = paste(
        "Per-dose-record indicator; 22/277 subjects (7.9%) contributed fed",
        "records (Comisar 2025 Table 2 'Fed, n (%)'). Acts only on oral",
        "bioavailability: Table 3 'Food effect on bioavailability of oral",
        "zavegepant (fed state as compared with fasted)' = -50.9%.",
        "Comisar 2025 Discussion states 'The effect of fed vs. fasted",
        "state was not evaluated on intranasal zavegepant administration",
        "as the intranasal absorption is not expected to be impacted', so",
        "the effect is applied only when ROUTE_ORAL = 1 and FED should be",
        "left at 0 on intranasal and intravenous records."
      ),
      source_name        = "FED"
    ),
    HEPIMP_MOD = list(
      description        = "Moderate hepatic impairment indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function; no mild or severe strata were represented in the analysis dataset)",
      notes              = paste(
        "Classification scheme is Child-Pugh, NOT NCI ODWG: Comisar 2025",
        "Methods lists 'hepatic function (hepatic impairment categorized",
        "by Child-Pugh scores)' and the Abstract specifies 'Moderate",
        "hepatic impairment (Child-Pugh score 7-9)', i.e. Child-Pugh",
        "Class B. Per-subject, time-fixed. 8/277 subjects (2.9%) were",
        "moderately impaired (Table 2 'Hepatic function, moderate",
        "impairment'), all from study BHV3500-108 which enrolled 8",
        "impaired and 8 matched normal-function subjects (Table 1).",
        "Effect on CL is a fractional reduction of -43.9% (Table 3).",
        "This is the only covariate the paper judged clinically relevant",
        "for the approved intranasal formulation, with predicted mean",
        "increases in AUC0-24 and Cmax of 63% and 12% respectively",
        "(Results 'Covariate analysis', Figure 2)."
      ),
      source_name        = "Hepatic function, moderate impairment"
    ),
    CONMED_RIF = list(
      description        = "Concomitant rifampin (rifampicin) administration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no rifampin)",
      notes              = paste(
        "Per-subject / per-occasion indicator from the BHV3500-111",
        "drug-drug interaction study; 15/277 subjects (5.4%) were",
        "co-treated (Comisar 2025 Table 2, which annotates the row 'oral",
        "zavegepant' because that DDI cohort received zavegepant orally).",
        "The retained effect is a -41.1% fractional reduction in CL",
        "(Table 3 'Co-administration of rifampin on CL'). Unlike the",
        "itraconazole row, the rifampin row carries NO route restriction",
        "in Table 3 and none in the Discussion ('Co-administration of",
        "rifampin ... reduced zavegepant clearance by 41.1%'), so it is",
        "encoded here as a route-independent effect on systemic CL.",
        "Comisar 2025 Table 1 explains why the two rows differ and why",
        "the distinction is moot in practice: BHV3500-111 evaluated",
        "rifampin on the PK of ORAL zavegepant only, but itraconazole on",
        "zavegepant 'administered orally or as a nasal spray'. There are",
        "therefore no intranasal-plus-rifampin observations in the",
        "dataset, so the authors had no route contrast to test for",
        "rifampin (and no need to restrict the row), whereas for",
        "itraconazole they did have one, tested it, and found the effect",
        "only for oral. On any dataset resembling the paper's own -- one",
        "where CONMED_RIF is set only on oral records -- the",
        "route-independent and oral-restricted encodings are numerically",
        "identical. See the vignette Errata for the several places where",
        "the paper's prose nonetheless attaches 'oral' or 'nasal spray'",
        "to this effect.",
        "Mechanistically the paper attributes the direction (a DECREASE in",
        "clearance despite rifampin being a CYP3A4 inducer) to a composite",
        "of OATP1B3 and NTCP inhibition dominating CYP3A induction."
      ),
      source_name        = "Co-treated with rifampin"
    ),
    CONMED_ITRACONAZOLE = list(
      description        = "Concomitant itraconazole administration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no itraconazole)",
      notes              = paste(
        "Per-subject / per-occasion indicator from the BHV3500-111 study;",
        "18/277 subjects (6.5%) received itraconazole with intranasal",
        "zavegepant and 18/277 (6.5%) with oral zavegepant (Comisar 2025",
        "Table 2). The effect is a -27.9% fractional reduction in CL that",
        "applies ONLY to orally administered zavegepant -- Table 3 labels",
        "the row 'Co-administration of itraconazole on clearance when",
        "zavegepant administered orally' and the Discussion states",
        "itraconazole 'was not a significant covariate for intranasal",
        "zavegepant clearance, but decreased oral zavegepant clearance by",
        "27.9%'. It is therefore encoded as the product",
        "CONMED_ITRACONAZOLE * ROUTE_ORAL, so setting",
        "CONMED_ITRACONAZOLE = 1 on an intranasal or intravenous record",
        "correctly produces no effect."
      ),
      source_name        = "Co-treated with itraconazole"
    )
  )

  # Covariates that Comisar 2025 screened by stepwise forward inclusion /
  # backward elimination but did NOT retain in the final model. Documented
  # here (rather than in covariateData) because no point estimate exists for
  # them and they are never referenced in model(). Comisar 2025 Results
  # 'Covariate analysis': 'Covariates including age, race, ethnicity, sex,
  # renal function (CrCL), oral contraceptives, and sumatriptan were not
  # significant covariates in affecting zavegepant exposure.'
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Median 40 years (range 18-71; Comisar 2025 Table 2). Screened and",
        "not retained. The Discussion cautions that only 6 subjects (1.0%)",
        "were >= 65 years, limiting the evaluation in older individuals."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "126/277 (45.5%) female per Comisar 2025 Table 2 'Sex, female'.",
        "Screened and not retained. Note that the Results narrative",
        "instead states 'female (49.3%)', which disagrees with Table 2;",
        "see the vignette Errata."
      )
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "32/277 (11.6%) Black, 245/277 (88.4%) Caucasian, no Asian or",
        "Other subjects (Comisar 2025 Table 2). Screened and not retained."
      )
    ),
    CRCL = list(
      description = "Creatinine clearance (renal function)",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Median 115 mL/min (range 53-211; Comisar 2025 Table 2). Screened",
        "and not retained, consistent with the limited role of renal",
        "clearance in zavegepant elimination. The Discussion notes only 17",
        "subjects had mild and 1 had moderate renal impairment, so the",
        "renal-function screen was underpowered at the impaired end."
      )
    ),
    CONMED_BIRTHCONTROL = list(
      description = "Concomitant oral contraceptive use indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "25/277 (9.0%) co-treated (Comisar 2025 Table 2). Screened in the",
        "BHV3500-109 ethinyl estradiol / levonorgestrel interaction study",
        "and not retained as a covariate on zavegepant PK."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 277,
    n_studies      = 9,
    age_range      = "18-71 years",
    age_median     = "40 years",
    weight_range   = "49-131 kg",
    weight_median  = "74 kg",
    sex_female_pct = 45.5,
    race_ethnicity = c(Caucasian = 88.4, Black = 11.6, Asian = 0, Other = 0),
    disease_state  = "healthy adults (85.9%) and patients with migraine (14.1%); 2.9% with moderate (Child-Pugh B) hepatic impairment",
    dose_range     = "intranasal 10-100 mg, oral 10-100 mg (soft gelatin capsule), and intravenous infusion; single and multiple dose",
    regions        = "North America",
    notes          = paste(
      "Comisar 2025 Table 2 (baseline demographic and laboratory data for",
      "the model-development dataset). 594 individuals were screened",
      "across all 10 phase I studies contributing 12253 concentration",
      "records; 1434 (11.7%) were excluded (226 pre-dose, 1208 below the",
      "0.4 ng/mL LLOQ), leaving 10819 samples from 277 participants for",
      "model development (Results 'Data summary'). Route split of the 277:",
      "209 intranasal (75.5%), 62 oral (22.4%), 6 intravenous (2.2%).",
      "Ethnicity was 159/277 (57.4%) Hispanic/Latino. n_studies = 9 counts",
      "only the model-development studies (BHV3500-101, -102, -104, -105,",
      "-107, -108, -109, -110, -111); the tenth study BHV3500-106 (317",
      "subjects, sparse sampling, 100 mg oral daily for 8 weeks) was held",
      "out for external validation by VPC and is not part of N = 277."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # Disposition. Comisar 2025 Table 3, 'Distribution' and 'Clearance'
    # blocks. The Table 3 Note states these are 'for intravenous dosing
    # (bioavailability of 100%) for a typical body weight of 70 kg', so they
    # are true (not apparent) CL and V values.
    #
    # Compartment mapping. The paper numbers compartments in NONMEM style
    # with the depot as 1, so its V2/V3/V4 and Q3/Q4 map onto the nlmixr2lib
    # canonical names as:
    #   V2 = 12.1 L  -> vc          (central)
    #   V3 = 66.0 L  -> vp   with Q3 = 2.6 L/h -> q    (peripheral1)
    #   V4 = 10.7 L  -> vp2  with Q4 = 5.2 L/h -> q2   (peripheral2)
    # Note that peripheral1 here is the LARGE, SLOWLY equilibrating
    # compartment (66.0 L at 2.6 L/h) and peripheral2 the small fast one
    # (10.7 L at 5.2 L/h); the ordering follows the paper's own numbering,
    # not a fast-then-slow convention. Figure 1's caption labels the same two
    # intercompartmental clearances Q2 and Q3 rather than Q3 and Q4 -- see
    # the vignette Errata; Table 3 and the Results text agree with each other
    # and are used here.
    # -----------------------------------------------------------------------
    lcl <- log(13.3) ; label("Elimination clearance from the central compartment (L/h)")                  # Comisar 2025 Table 3 'Elimination clearance, CL' = 13.3 L/h (RSE 9.9%)
    lvc <- log(12.1) ; label("Central volume of distribution (L)")                                        # Comisar 2025 Table 3 'Central volume of distribution, V2' = 12.1 L (RSE 8.8%)
    lq  <- log(2.6)  ; label("Intercompartmental clearance to the first peripheral compartment (L/h)")    # Comisar 2025 Table 3 'Intercompartmental clearance between central compartment and first peripheral compartment, Q3' = 2.6 L/h (RSE 10.5%)
    lvp <- log(66.0) ; label("First peripheral volume of distribution (L)")                               # Comisar 2025 Table 3 'First peripheral volume of distribution, V3' = 66.0 L (RSE 10.7%)
    lq2 <- log(5.2)  ; label("Intercompartmental clearance to the second peripheral compartment (L/h)")   # Comisar 2025 Table 3 'Intercompartmental clearance between central compartment and second peripheral compartment, Q4' = 5.2 L/h (RSE 16.3%)
    lvp2 <- log(10.7); label("Second peripheral volume of distribution (L)")                              # Comisar 2025 Table 3 'Second peripheral volume of distribution, V4' = 10.7 L (RSE 10.1%)

    # -----------------------------------------------------------------------
    # Absorption. Comisar 2025 Table 3 'Absorption' block. Every quantity is
    # reported once for intranasal and once for oral, so each carries an
    # explicit stratum suffix (parameter-names.md 'Stratum-suffixed
    # parameters'); none keeps the bare canonical name.
    #
    # D1 and the absorption lag are tabulated in MINUTES while the rest of
    # the model is on an hours time base (ka in 1/h, CL in L/h). They are
    # converted here by dividing by 60, which is exactly what the paper does
    # itself: Results 'Population pharmacokinetic model development' quotes
    # 'a D1 of 0.14 h' for intranasal (8.6/60 = 0.14333) and 'a D1 of 0.95 h'
    # for oral (57.2/60 = 0.95333). The minute values are used rather than
    # the paper's own 2-significant-figure hour roundings.
    # -----------------------------------------------------------------------
    lfdepot_intranasal <- log(0.051)     ; label("Bioavailability of intranasal zavegepant (fraction)")             # Comisar 2025 Table 3 'Bioavailability of intranasal zavegepant' = 5.1% (RSE 11.0%)
    lfdepot_oral       <- log(0.0065)    ; label("Bioavailability of oral zavegepant, fasted (fraction)")           # Comisar 2025 Table 3 'Bioavailability of oral zavegepant' = 0.65% (RSE 17.5%)
    lka_intranasal     <- log(5.8)       ; label("First-order absorption rate constant, intranasal (1/h)")          # Comisar 2025 Table 3 'First-order absorption rate constant (ka) for intranasal zavegepant' = 5.8 1/h (RSE 5.1%)
    lka_oral           <- log(0.81)      ; label("First-order absorption rate constant, oral (1/h)")                # Comisar 2025 Table 3 'First-order absorption rate constant (ka) for oral zavegepant' = 0.81 1/h (RSE 6.1%)
    ld1_intranasal     <- log(8.6 / 60)  ; label("Duration of zero-order absorption, intranasal (h)")               # Comisar 2025 Table 3 'Duration of zero-order absorption for intranasal zavegepant' = 8.6 min (RSE 8.7%); = 0.14 h in Results
    ld1_oral           <- log(57.2 / 60) ; label("Duration of zero-order absorption, oral (h)")                     # Comisar 2025 Table 3 'Duration of zero-order absorption for oral zavegepant' = 57.2 min (RSE 14.4%); = 0.95 h in Results
    ltlag_oral         <- log(12.2 / 60) ; label("Absorption lag time, oral (h)")                                   # Comisar 2025 Table 3 'Absorption lag for oral zavegepant' = 12.2 min (RSE 4.3%); no lag is reported for intranasal

    # -----------------------------------------------------------------------
    # Allometry. Both exponents are FIXED at the empirical standard values,
    # not estimated: Comisar 2025 Methods states 'standard exponents (0.75
    # for clearance and 1 for volume of distribution) and a reference body
    # weight of 70 kg', and the Table 3 Note repeats 'Volumes of distribution
    # and clearance rates have body weight exponents fixed to 1.0 and 0.75'.
    # Neither appears as an estimated row in Table 3 and neither carries an
    # RSE, which is the other fixed-parameter signal. A single exponent is
    # applied to all three clearances and a single one to all three volumes
    # ('tested and applied to all disposition parameters', Results
    # 'Covariate analysis').
    # -----------------------------------------------------------------------
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent on all clearances, reference 70 kg (unitless)")   # Comisar 2025 Methods 'Population pharmacokinetic modeling'; Table 3 Note
    e_wt_vc <- fixed(1)    ; label("Allometric exponent on all volumes, reference 70 kg (unitless)")      # Comisar 2025 Methods 'Population pharmacokinetic modeling'; Table 3 Note

    # -----------------------------------------------------------------------
    # Covariate effects. Comisar 2025 Table 3 'Covariates' block reports all
    # four as signed PERCENTAGES, which identifies the fractional form of the
    # paper's categorical covariate equation (Methods: 'theta_x is the
    # proportional constant of the covariate effect model', applied to a
    # 0/1 indicator X_i). Each coefficient is therefore entered as the
    # decimal fraction and used as a (1 + theta * X) multiplier, so that
    # e.g. moderate hepatic impairment gives CL * (1 - 0.439) = 0.561 * CL.
    # -----------------------------------------------------------------------
    e_hepimp_mod_cl          <- -0.439 ; label("Effect of moderate hepatic impairment on CL (fraction)")                    # Comisar 2025 Table 3 'Moderate hepatic impairment on CL' = -43.9% (RSE 13.5%)
    e_conmed_rif_cl          <- -0.411 ; label("Effect of concomitant rifampin on CL (fraction)")                           # Comisar 2025 Table 3 'Co-administration of rifampin on CL' = -41.1% (RSE 12.1%)
    e_conmed_itraconazole_cl <- -0.279 ; label("Effect of concomitant itraconazole on CL, oral dosing only (fraction)")     # Comisar 2025 Table 3 'Co-administration of itraconazole on clearance when zavegepant administered orally' = -27.9% (RSE 26.3%)
    e_fed_fdepot_oral        <- -0.509 ; label("Effect of fed state on oral bioavailability (fraction)")                    # Comisar 2025 Table 3 'Food effect on bioavailability of oral zavegepant (fed state as compared with fasted)' = -50.9% (RSE 9.8%)

    # -----------------------------------------------------------------------
    # Inter-individual variability. Comisar 2025 Table 3 column
    # '%IIV (%RSE) [%Shrinkage]'. Each %IIV cell in that column is merged
    # across the intranasal and oral rows of the quantity it belongs to, so a
    # SINGLE eta is shared between the two routes for bioavailability, ka and
    # D1 -- consistent with Methods, which lists the IIV terms evaluated as
    # 'bioavailability of intranasal and soft gelatin capsules, absorption
    # rate constant (ka) for intranasal and soft gelatin capsules, and
    # duration of absorption (D1) for intranasal and soft gelatin capsules'
    # (one term per quantity, not one per route x quantity).
    #
    # The %IIV values are read as log-normal CV%, so
    #   omega^2 = log(CV^2 + 1):
    #     F   66.1% -> log(0.661^2 + 1) = 0.3625026
    #     ka  44.0% -> log(0.440^2 + 1) = 0.1769740
    #     D1  115%  -> log(1.150^2 + 1) = 0.8426442
    #     Vc  32.2% -> log(0.322^2 + 1) = 0.0986537
    #     Vp  26.1% -> log(0.261^2 + 1) = 0.0659010
    #     CL  26.8% -> log(0.268^2 + 1) = 0.0693619
    # Table 3 gives no footnote naming the convention; see the vignette
    # Errata for the alternative (omega^2 = CV^2) reading and its numerical
    # consequence. No IIV is reported for Vp2, Q, Q2, the oral absorption lag
    # or any covariate coefficient, and none is invented here. Table 3
    # reports no off-diagonal covariance terms, so all etas are independent.
    # -----------------------------------------------------------------------
    etalfdepot ~ 0.3625026   # Comisar 2025 Table 3 bioavailability rows, %IIV = 66.1% (RSE 4.8%, shrinkage 6.8%); shared by the intranasal and oral bioavailabilities
    etalka     ~ 0.1769740   # Comisar 2025 Table 3 ka rows, %IIV = 44.0% (RSE 11.1%, shrinkage 41.1%); shared by the intranasal and oral ka
    etald1     ~ 0.8426442   # Comisar 2025 Table 3 D1 rows, %IIV = 115% (RSE 7.7%, shrinkage 11.6%); shared by the intranasal and oral D1
    etalvc     ~ 0.0986537   # Comisar 2025 Table 3 'Central volume of distribution, V2', %IIV = 32.2% (RSE 10.9%, shrinkage 28%)
    etalvp     ~ 0.0659010   # Comisar 2025 Table 3 'First peripheral volume of distribution, V3', %IIV = 26.1% (RSE 23.1%, shrinkage 49.9%)
    etalcl     ~ 0.0693619   # Comisar 2025 Table 3 'Elimination clearance, CL', %IIV = 26.8% (RSE 6.1%, shrinkage 17.8%)

    # -----------------------------------------------------------------------
    # Residual error. Comisar 2025 Methods: 'A combined additive and
    # proportional residual error model was used.' Table 3 'Residual errors'
    # reports the proportional term as 37.8% (RSE 2.9%) and the additive term
    # as '0.001 ng/L (FIXED)' with no RSE, hence the fixed() wrapper.
    #
    # The additive term's printed unit is almost certainly a typo for ng/mL
    # (the assay unit; the LLOQ is 0.4 ng/mL per Results 'Data summary').
    # Under either reading it is a numerically negligible stabilising
    # placeholder five orders of magnitude below the LLOQ, so it is entered
    # as the printed 0.001 on the model's ng/mL concentration scale. See the
    # vignette Errata.
    # -----------------------------------------------------------------------
    propSd <- 0.378         ; label("Proportional residual error (fraction)")   # Comisar 2025 Table 3 'Proportional error' = 37.8% (RSE 2.9%)
    addSd  <- fixed(0.001)  ; label("Additive residual error (ng/mL as printed)")  # Comisar 2025 Table 3 'Additive error (FIXED)' = 0.001, printed as ng/L
  })

  model({
    # ---------------------------------------------------------------------
    # 1. Route-specific absorption parameters.
    #
    # The route switch is applied on the LOG scale inside a single exp() so
    # that each parameter keeps a clean mu-referenced form. With
    # ROUTE_ORAL = 0 the intranasal value is selected and with ROUTE_ORAL = 1
    # the oral value; the shared eta rides on both, matching the merged %IIV
    # cells of Comisar 2025 Table 3.
    #
    # The fed-state effect is multiplied by ROUTE_ORAL as well, because
    # Comisar 2025 evaluated food only on the oral formulation (Discussion:
    # 'The effect of fed vs. fasted state was not evaluated on intranasal
    # zavegepant administration'). tlag likewise collapses to 0 for the
    # intranasal route, which reports no absorption lag.
    # ---------------------------------------------------------------------
    fdepot <- exp((1 - ROUTE_ORAL) * lfdepot_intranasal +
                  ROUTE_ORAL * lfdepot_oral + etalfdepot) *
      (1 + e_fed_fdepot_oral * FED * ROUTE_ORAL)
    ka <- exp((1 - ROUTE_ORAL) * lka_intranasal +
              ROUTE_ORAL * lka_oral + etalka)
    d1 <- exp((1 - ROUTE_ORAL) * ld1_intranasal +
              ROUTE_ORAL * ld1_oral + etald1)
    tlag <- ROUTE_ORAL * exp(ltlag_oral)

    # ---------------------------------------------------------------------
    # 2. Disposition parameters. Allometry on every clearance and every
    #    volume, referenced to 70 kg. The three CL covariates enter as the
    #    paper's fractional (1 + theta * X) multipliers; the itraconazole
    #    term is gated on ROUTE_ORAL so that it contributes nothing on
    #    intranasal or intravenous records.
    #
    #    At WT = 70 kg with all indicators 0 every multiplier collapses to 1
    #    and the parameters reduce exactly to the Table 3 typical values.
    # ---------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl *
      (1 + e_hepimp_mod_cl * HEPIMP_MOD) *
      (1 + e_conmed_rif_cl * CONMED_RIF) *
      (1 + e_conmed_itraconazole_cl * CONMED_ITRACONAZOLE * ROUTE_ORAL)
    vc  <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q   <- exp(lq)   * (WT / 70)^e_wt_cl
    vp  <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc
    q2  <- exp(lq2)  * (WT / 70)^e_wt_cl
    vp2 <- exp(lvp2) * (WT / 70)^e_wt_vc

    # 3. Micro-constants for the three-compartment system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ---------------------------------------------------------------------
    # 4. ODE system. Sequential zero-order then first-order absorption into
    #    the central compartment (Comisar 2025 Figure 1): the dose is
    #    delivered into `depot` at a constant rate over the zero-order window
    #    of duration D1, and `depot` then drains into `central` first-order
    #    at ka. Extravascular dose records must therefore carry rate = -2 so
    #    rxode2 uses the model's dur(depot); a plain bolus would collapse the
    #    zero-order phase and overstate Cmax.
    #
    #    Intravenous doses bypass all of this and are placed directly into
    #    `central` with cmt = 'central' (with ROUTE_ORAL = 0), which is why
    #    no bioavailability term is applied to the central compartment --
    #    Table 3's disposition values are already the IV (F = 100%) values.
    # ---------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    f(depot)    <- fdepot
    dur(depot)  <- d1
    alag(depot) <- tlag

    # ---------------------------------------------------------------------
    # 5. Observation. Dose is in mg and vc in L, so central / vc is mg/L
    #    (= ug/mL); multiplying by 1000 gives ng/mL, the assay unit used
    #    throughout Comisar 2025 (LLOQ 0.4 ng/mL, Results 'Data summary').
    # ---------------------------------------------------------------------
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
