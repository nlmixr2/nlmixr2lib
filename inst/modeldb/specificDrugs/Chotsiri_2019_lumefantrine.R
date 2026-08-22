# Population pharmacokinetic-pharmacodynamic model for oral lumefantrine in
# children with and without severe acute malnutrition treated with
# artemether-lumefantrine for uncomplicated Plasmodium falciparum malaria
# (Chotsiri 2019, Clin Pharmacol Ther 106(6):1299-1309;
# doi:10.1002/cpt.1531; MAL-NUT trial NCT01958905).

Chotsiri_2019_lumefantrine <- function() {
  description <- paste(
    "Population PK-PD model for oral lumefantrine in 131 severely acutely",
    "malnourished (SAM) and 160 non-SAM children aged 6-59 months with",
    "uncomplicated Plasmodium falciparum malaria in Mali and Niger",
    "(Chotsiri 2019 Clin Pharmacol Ther; MAL-NUT trial NCT01958905). Two",
    "transit-absorption compartments followed by two-compartment",
    "disposition, with the absorption rate constant set equal to the",
    "transit rate constant (Savic 2007) so ktr = ka = 3 / MTT. Relative",
    "bioavailability F is fixed at 1 with Box-Cox-transformed IIV",
    "(lambda = -0.373, CV 64.0%; Petersson 2009 form). Fixed allometric",
    "scaling of CL/F and Q/F (power 3/4) and of Vc/F and Vp/F (power 1) on",
    "body weight centered at the model-building median 9.62 kg, plus a",
    "first-order enzyme-maturation effect on CL/F (TM50 = 2.91 months,",
    "shape alpha fixed at 1) driven by postnatal age. Mid-upper arm",
    "circumference enters as an exponential covariate on F (25.4% lower",
    "bioavailability per 1 cm reduction in MUAC) and is the only",
    "malnutrition covariate retained in the final model. A dose-saturable",
    "absorption term (Dose50 = 3.86 mg/kg, fixed from Kloprogge 2018) is",
    "carried for dose-optimisation simulations, normalised to each",
    "child's standard 120 mg per-dose exposure so that it is inert at the",
    "standard regimen and reduces bioavailability only at higher dose",
    "levels (see vignette Errata). Coupled interval-censored time-to-event model",
    "for P. falciparum reinfection over 42 days: constant baseline hazard",
    "(5.25 reinfections/year) multiplied by a sigmoid Emax protective",
    "effect of lumefantrine (IC50 = 156 ng/mL, gamma = 4.77, Emax fixed",
    "at 1), exposed as the survivor function `sur`. Companion",
    "lumefantrine models: modellib('Kloprogge_2018_lumefantrine').",
    sep = " "
  )
  reference <- paste(
    "Chotsiri P, Denoeud-Ndam L, Baudin E, Guindo O, Diawara H, Attaher O,",
    "et al. (2019). Severe acute malnutrition results in lower",
    "lumefantrine exposure in children treated with artemether-lumefantrine",
    "for uncomplicated malaria. Clinical Pharmacology and Therapeutics",
    "106(6):1299-1309. doi:10.1002/cpt.1531. The packaged vignette also",
    "reproduces the dosing-regimen simulations of Simeon S, Hughes E,",
    "Wallender E, Solans BP, Savic R (2024). Optimizing lumefantrine",
    "dosing for young children in high-malaria-burden countries using",
    "pharmacokinetic-pharmacodynamic simulations. Open Forum Infectious",
    "Diseases 11(11):ofae627. doi:10.1093/ofid/ofae627, which applies this",
    "model without re-estimating it.",
    sep = " "
  )
  vignette <- "Chotsiri_2019_lumefantrine"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "lumefantrine", units = "mg", specimen = "administration site", verified = FALSE),
    transit1    = list(analyte = "lumefantrine", units = "mg", specimen = "administration site", verified = FALSE),
    transit2    = list(analyte = "lumefantrine", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "lumefantrine", units = "mg", specimen = "whole blood", verified = FALSE),
    peripheral1 = list(analyte = "lumefantrine", units = "mg", specimen = "whole blood", verified = FALSE),
    cumhaz      = list(analyte = "Cumulative hazard of P. falciparum reinfection", units = NA_character_, specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject at admission. Chotsiri 2019 Table 1",
        "reports median (95% CI) body weight 6.78 (5.07, 10.5) kg in the",
        "131 SAM children and 10.9 (6.84, 15.9) kg in the 160 non-SAM",
        "children of the PK-PD arm. Applied as a fixed allometric scaler",
        "with exponent 3/4 on CL/F and Q/F and exponent 1 on Vc/F and",
        "Vp/F, centered at the pooled model-building median 9.62 kg",
        "(supplementary methods, 'Population pharmacokinetic analysis':",
        "'Individual body-weight (BWi) was introduced into the",
        "pharmacokinetic model as a fixed allometric function on all",
        "volume (exponent of n = 1.00) and clearance (exponent of",
        "n = 0.75) parameters, scaled to the median body weight (9.62 kg)",
        "of the study population', theta_i = theta * exp(eta) *",
        "(BW_i / 9.62)^n). Body weight also determines the weight-band",
        "tablet count (1 tablet below 15 kg, 2 tablets 15-25 kg; Methods,",
        "'Study design') and therefore the per-dose mg/kg exposure that",
        "drives the dose-saturable bioavailability term.",
        sep = " "
      ),
      source_name        = "BW"
    ),
    PNA = list(
      description        = "Postnatal age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Postnatal (chronological) age in months. Chotsiri 2019 Table 1",
        "reports median (95% CI) age 15 (6.3, 39) months in SAM children",
        "and 27 (8.0, 53) months in non-SAM children; the trial enrolled",
        "children aged 6-59 months. Drives the first-order enzyme-",
        "maturation effect on CL/F (supplementary methods: MF = AGE^alpha",
        "/ (TM50^alpha + AGE^alpha), CL_i = theta_CL * exp(eta_CL) *",
        "(BW_i / 9.62)^0.75 * MF), with TM50 = 2.91 months and the shape",
        "factor alpha fixed at 1 (Table 2). The canonical column PNA",
        "already carries months, so the source paper's AGE (months) maps",
        "to it with no unit conversion. Note the maturation factor is NOT",
        "normalised to a reference age: theta_CL = 2.34 L/h is the fully",
        "mature apparent clearance of a 9.62 kg child, and MF < 1 at every",
        "age in the studied range (MF = 0.838 at the SAM median 15 months",
        "and 0.903 at the non-SAM median 27 months).",
        sep = " "
      ),
      source_name        = "AGE"
    ),
    MUAC = list(
      description        = "Mid-upper arm circumference",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject at admission. Chotsiri 2019 Table 1",
        "reports median (95% CI) MUAC 116 (104, 131) mm in SAM children",
        "and 140 (122, 163) mm in non-SAM children. MUAC was the most",
        "significant of the three malnutrition indicators screened (MUAC,",
        "weight-for-height z-score, weight-for-age z-score; largest drop",
        "in objective function, dOFV = -64.4) and is the only malnutrition",
        "covariate retained in the final model (Results, 'Covariate",
        "model'). Applied as a continuous exponential effect on relative",
        "bioavailability: F_muac = exp(e_muac_f * (MUAC / 10 - muac_ref)),",
        "with e_muac_f = -log(1 - 0.254) per cm so that each 1 cm",
        "reduction in MUAC multiplies F by (1 - 0.254), i.e. the '25.4%",
        "decreased absorption per 1 cm reduction' of the Abstract and",
        "Table 2. MUAC is stored here in millimetres (the unit Table 1",
        "tabulates and the unit of the WHO SAM threshold, MUAC < 115 mm)",
        "and converted to centimetres inside model() because the",
        "coefficient is published per centimetre. The centering value",
        "muac_ref is NOT reported by the source paper and is back-solved",
        "from the paper's own Table 3 secondary parameters -- see the",
        "ini() comment on muac_ref and the vignette Errata. Distinct from",
        "the binary MAL_NOURISH canonical: MAL_NOURISH is a paper-defined",
        "malnourished/not indicator, whereas MUAC is the continuous",
        "anthropometric measurement itself (and is one of the criteria",
        "papers use to SET MAL_NOURISH).",
        sep = " "
      ),
      source_name        = "MUAC"
    ),
    DOSE = list(
      description        = "Per-dose lumefantrine amount administered (mg)",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-dose lumefantrine amount in milligrams, supplied as a",
        "per-dose-record covariate aligned with the corresponding event-",
        "table dosing row (use case (b) of the canonical DOSE entry,",
        "applied per dose record). Required to compute the per-dose",
        "milligram-per-kilogram exposure dose_mgkg = DOSE / WT, which",
        "drives the dose-saturable bioavailability term. Fixed-dose",
        "combination tablets of 20 mg artemether + 120 mg lumefantrine",
        "(Coartem) were given twice daily for 3 days, 1 tablet below",
        "15 kg and 2 tablets for 15-25 kg (Methods, 'Study design'), so",
        "DOSE = 120 for essentially every child in this cohort and",
        "DOSE = 240 in the 'increased' dose-optimisation regimen",
        "(supplementary methods, 'In silico lumefantrine dose",
        "optimisation'). Set DOSE to the milligram amount administered at",
        "each dose event in the rxode2 event table, alongside amt (mg).",
        "Matches the usage in Kloprogge_2018_lumefantrine.R, the source of",
        "the Dose50 = 3.86 mg/kg saturation constant.",
        sep = " "
      ),
      source_name        = "Dosage"
    )
  )

  # Screened but not retained in the final model. Chotsiri 2019 evaluated
  # these in a stepwise covariate search and, for the three malnutrition
  # indicators, additionally in a full covariate approach (Results,
  # 'Covariate model'; supplementary methods). None is referenced in
  # model(), so they are documented here rather than in covariateData.
  covariatesDataExcluded <- list(
    WHZ = list(
      description = "Weight-for-height z-score",
      units       = "(z-score)",
      type        = "continuous",
      notes       = paste(
        "Screened as a malnutrition indicator on relative bioavailability",
        "and found significant, but MUAC gave the larger drop in objective",
        "function (dOFV = -64.4) and was retained instead; adding any",
        "further malnutrition indicator alongside MUAC gave no additional",
        "improvement. The full covariate approach reports a 15.5% (95% CI",
        "3.07-33.2%) reduction in median bioavailability per 1 unit",
        "weight-for-age z-score reduction. Chotsiri 2019 Table 1 medians:",
        "-3.52 (SAM), -1.07 (non-SAM). Not retained -> not encoded.",
        sep = " "
      )
    ),
    WAZ = list(
      description = "Weight-for-age z-score",
      units       = "(z-score)",
      type        = "continuous",
      notes       = "Screened alongside MUAC and WHZ; not retained in the final model (see WHZ note). Chotsiri 2019 Table 1 medians: -3.40 (SAM), -1.30 (non-SAM)."
    ),
    HAZ = list(
      description = "Height-for-age z-score",
      units       = "(z-score)",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate search; not significant and not retained. Chotsiri 2019 Table 1 medians: -1.67 (SAM), -1.27 (non-SAM). Children with severe stunting (HAZ < -3) were excluded from the trial."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in the stepwise covariate search (supplementary methods lists sex among the covariates investigated); not retained. Chotsiri 2019 Table 1 reports 50.4% male in SAM and 46.9% male in non-SAM children."
    ),
    TEMP = list(
      description = "Axillary body temperature at admission",
      units       = "degC",
      type        = "continuous",
      notes       = "Screened as 'body temperature at admission' in the stepwise covariate search; not retained. Chotsiri 2019 Table 1 medians: 38.0 degC (SAM), 38.2 degC (non-SAM)."
    ),
    PARA = list(
      description = "Plasmodium falciparum parasitaemia at admission (asexual parasites/uL)",
      units       = "parasites/uL",
      type        = "continuous",
      notes       = paste(
        "Admission parasitaemia; Chotsiri 2019 Table 1 medians 11,040",
        "(SAM) and 10,780 (non-SAM) parasites/uL. Not a covariate in this",
        "model. Recurrent parasite density is used to back-extrapolate the",
        "interval-censoring window of the time-to-event model, and the",
        "sibling Kloprogge_2018_lumefantrine.R does carry an exponential",
        "parasitaemia effect on F, but Chotsiri 2019 did not retain one.",
        sep = " "
      )
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 291L,
    n_subjects_pk   = 291L,
    n_subjects_pd   = 380L,
    n_studies       = 1L,
    age_range       = "6-59 months (inclusion criterion); median 15 months (95% CI 6.3-39) in SAM and 27 months (95% CI 8.0-53) in non-SAM children (Table 1)",
    weight_range    = "median 6.78 kg (95% CI 5.07-10.5) in SAM and 10.9 kg (95% CI 6.84-15.9) in non-SAM children (Table 1); allometric centering value 9.62 kg",
    sex_female_pct  = 51.2,
    disease_state   = paste(
      "Uncomplicated Plasmodium falciparum malaria, with and without",
      "severe acute malnutrition (SAM). SAM was defined per WHO criteria",
      "as weight-for-height z-score < -3 and/or mid-upper arm",
      "circumference < 115 mm. Children with kwashiorkor, severe stunting",
      "(height-for-age z-score < -3), severe anaemia, known underlying or",
      "chronic disease, and complications requiring hospitalisation were",
      "excluded.",
      sep = " "
    ),
    dose_range      = paste(
      "Fixed-dose combination tablets of 20 mg artemether + 120 mg",
      "lumefantrine (Coartem, non-dispersible), weight-band dosed",
      "(1 tablet below 15 kg, 2 tablets 15-25 kg) twice daily for 3 days",
      "(6 doses). Total lumefantrine dose 105 mg/kg (95% CI 64.8-142) in",
      "SAM and 68.2 mg/kg (95% CI 48.3-106) in non-SAM children (Table 1).",
      "Doses were given with fat: one glass of milk (~250 mL, 2.5 g fat)",
      "in non-SAM children and one 92 g bag of ready-to-use therapeutic",
      "food (Plumpy'Nut, 32.9 g fat) in SAM children.",
      sep = " "
    ),
    regions         = "Mali (Oulessebougou District Hospital, Koulikoro region) and Niger (primary healthcare centre of Andoume, Maradi City)",
    notes           = paste(
      "MAL-NUT trial, ClinicalTrials.gov NCT01958905. The PK analysis used",
      "1,342 capillary dried-blood-spot lumefantrine concentrations (642",
      "from 131 SAM and 700 from 160 non-SAM children); 84 of 1,341",
      "samples (6.26%) were below the 39.1 ng/mL LLOQ and were handled by",
      "the M6 method (first LLOQ record per patient imputed at half LLOQ).",
      "Five samples per child were drawn at a randomly allocated time",
      "among 6, 12, 24, 36 or 48 h, then at 60 h, 72 h, day 7, and either",
      "day 14 or day 21. The PD (time-to-reinfection) analysis pooled the",
      "PK-PD arm with a separate PD-only arm of 108 non-SAM children for",
      "380 analysable children and 95 reinfections over 42 days of",
      "follow-up. Concentrations are capillary blood measured by LC-MS/MS",
      "after solid-phase extraction. Between-subject variability on CL/F,",
      "Vc/F and Vp/F was estimated close to zero and removed from the",
      "final model; eta shrinkage was 53.2% on MTT, 13.2% on F and 24.4%",
      "on Q/F, and epsilon shrinkage was 18.0%. Percentage of female",
      "children is derived from the Table 1 male percentages weighted by",
      "the SAM (n = 131) and non-SAM (n = 160) PK-arm group sizes.",
      sep = " "
    )
  )

  ini({
    # ---- Absorption -------------------------------------------------
    # Two transit compartments with the absorption rate constant set
    # equal to the transit rate constant (supplementary methods: 'For the
    # transit absorption model, the absorption constant (Ka) and transit
    # rate constant (Ktr) were assumed to be equal'; Savic 2007). With
    # n = 2 transit compartments the mean transit time from the dosing
    # compartment to the central compartment spans n + 1 = 3 first-order
    # transfers, so ktr = ka = 3 / MTT = 3 / 3.48 = 0.862069 /h. MTT is
    # the estimated parameter and carries the (large) IIV, so it is the
    # quantity parameterised here.
    lmtt <- log(3.48)  ; label("Mean transit time MTT (h)")                                   # Chotsiri 2019 Table 2: MTT = 3.48 h (RSE 11.7%; bootstrap median 3.48, 95% CI 2.52-4.35)

    # ---- Disposition ------------------------------------------------
    # Two-compartment disposition. Table 2 population estimates are for a
    # typical individual with a body weight of 9.62 kg (Table 2 footnote
    # a). CL/F is the FULLY MATURE apparent clearance: the maturation
    # factor MF multiplies it un-normalised (see model()).
    lcl  <- log(2.34)  ; label("Apparent elimination clearance CL/F (L/h, fully mature, typical 9.62-kg child)")  # Chotsiri 2019 Table 2: CL/F = 2.34 L/h (RSE 6.25%; bootstrap median 2.34, 95% CI 2.20-2.81)
    lvc  <- log(110)   ; label("Apparent central volume of distribution Vc/F (L, typical 9.62-kg child)")          # Chotsiri 2019 Table 2: VC/F = 110 L (RSE 4.46%; bootstrap median 110, 95% CI 101-123)
    lq   <- log(1.10)  ; label("Apparent intercompartmental clearance Q/F (L/h, typical 9.62-kg child)")           # Chotsiri 2019 Table 2: QP/F = 1.10 L/h (RSE 8.10%; bootstrap median 1.10, 95% CI 0.942-1.33)
    lvp  <- log(872)   ; label("Apparent peripheral volume of distribution Vp/F (L, typical 9.62-kg child)")       # Chotsiri 2019 Table 2: VP/F = 872 L (RSE 13.9%; bootstrap median 872, 95% CI 635-1,200)

    # Relative bioavailability, fixed to unity for the population
    # (supplementary methods: 'Relative bioavailability was fixed to unity
    # for the population, but allowed for estimation of inter-individual
    # variability in the same parameter'). Absolute bioavailability is not
    # identifiable from oral data alone; all of the F signal is carried by
    # the Box-Cox IIV and the MUAC covariate below.
    lfdepot <- fixed(log(1)) ; label("Relative bioavailability F (unitless)")                 # Chotsiri 2019 Table 2: F = 1 (fixed)

    # Box-Cox shape parameter for the IIV on F (Petersson 2009 form; see
    # model()). The negative estimate indicates a left-skewed distribution
    # of relative bioavailability, which the authors attribute to the
    # dose-limited absorption of lumefantrine (Discussion).
    boxcox_fdepot <- -0.373 ; label("Box-Cox shape parameter lambda for the IIV on relative bioavailability (unitless)") # Chotsiri 2019 Table 2: Box-Cox on F = -0.373 (RSE 54.5%; bootstrap median -0.373, 95% CI -0.742 to 0.205)

    # ---- Allometric scaling -----------------------------------------
    # Fixed exponents, centered on the model-building median body weight
    # 9.62 kg (supplementary methods; theta_i = theta * exp(eta) *
    # (BW_i / 9.62)^n with n = 0.75 for clearances and n = 1.00 for
    # volumes). Introduced a priori on biological grounds and retained in
    # the final model despite not improving the fit (Results, 'Covariate
    # model': dOFV = 66.1).
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent on CL/F and Q/F (unitless)")          # Chotsiri 2019 supplementary methods: clearance exponent n = 0.75 (fixed)
    e_wt_vc <- fixed(1)    ; label("Allometric exponent on Vc/F and Vp/F (unitless)")         # Chotsiri 2019 supplementary methods: volume exponent n = 1.00 (fixed)

    # ---- Enzyme maturation on CL/F ----------------------------------
    # MF = AGE^alpha / (TM50^alpha + AGE^alpha) with AGE in months
    # (supplementary methods; Anderson & Holford 2009 form). alpha is
    # fixed at 1 by the source paper, reducing MF to the hyperbolic form
    # PNA / (TM50 + PNA).
    tm50_cl <- 2.91     ; label("Maturation half-life TM50 for CL/F (months); postnatal age at 50% mature CL/F") # Chotsiri 2019 Table 2: TM50 = 2.91 months (RSE 31.5%; bootstrap median 2.91, 95% CI 2.86-5.98)
    hill_cl <- fixed(1) ; label("Shape factor alpha for the enzyme-maturation effect on CL/F (unitless)")         # Chotsiri 2019 Table 2: alpha = 1 (fixed)

    # ---- MUAC effect on relative bioavailability --------------------
    # The final model's only malnutrition covariate. Table 2 reports the
    # effect as '25.4% (per 1 cm)' and the Abstract as '25.4% decreased
    # absorption per 1 cm reduction' in MUAC. Encoded as a constant
    # proportional (exponential) effect so that the published '% per 1 cm'
    # holds at every MUAC rather than only near the reference:
    #   F_muac = exp(e_muac_f * (MUAC_cm - muac_ref))
    # with e_muac_f = -log(1 - 0.254) = 0.29303 per cm, i.e. F falls by a
    # factor of (1 - 0.254) for each 1 cm reduction. The source paper does
    # not print the covariate equation, so the functional form is inferred
    # (see vignette Errata). The exponential reading is preferred because
    # (i) the paper states the effect as a single '% per 1 cm' figure
    # applied across the whole MUAC range, in both the final model (25.4%)
    # and the full covariate model (21.0%), which is only self-consistent
    # for a constant proportional effect; (ii) it keeps F strictly
    # positive across the observed 104-163 mm MUAC range; and (iii) it
    # matches the multiplicative exp(...) form used for every other
    # covariate in this model family. A linear form gives numerically
    # similar predictions over the observed range and is not excluded by
    # any published quantity.
    e_muac_f <- -log(1 - 0.254) ; label("Exponential effect of MUAC on relative bioavailability (per cm; F falls 25.4% per 1 cm reduction)") # Chotsiri 2019 Table 2: MUAC on F = 25.4% per 1 cm (RSE 5.24%; bootstrap median 25.4%, 95% CI 21.3-27.1%)

    # Centering value for the MUAC effect. NOT REPORTED ANYWHERE in
    # Chotsiri 2019 or its supplement, and load-bearing: because theta_F
    # is fixed to 1, muac_ref alone sets the absolute scale of F and hence
    # of every predicted concentration (a 1 cm shift moves all
    # concentrations by 29%). Two independent lines of evidence agree on
    # 130 mm:
    #  (a) the pooled central MUAC of the PK population -- 129.1 mm as the
    #      group-size-weighted mean of the Denoeud-Ndam 2016 Table 1 group
    #      means, and ~131 mm from the pooled-median position implied by
    #      the Chotsiri Table 1 group medians (116 and 140 mm) and group
    #      sizes (131 and 160). This mirrors how the paper centres the
    #      allometry on the pooled median body weight, 9.62 kg.
    #  (b) back-solving from the paper's own Table 3 secondary parameters:
    #      reproducing each group's published median AUC0-28d, and
    #      independently its published median day-7 concentration, at that
    #      group's median body weight, age and MUAC gives 128 mm from the
    #      SAM rows and 130 mm from the non-SAM rows.
    # 13.0 cm is adopted as the rounded standard per the undefined-
    # centering-value convention. At this value the model reproduces every
    # published Table 3 secondary parameter except the terminal half-life
    # to within 8% -- see the vignette source-trace and Errata.
    muac_ref <- fixed(13.0) ; label("Centering value for the MUAC effect on F (cm); not reported by the source -- derived, see vignette Errata")

    # ---- Dose-saturable absorption ----------------------------------
    # Lumefantrine absorption saturates with increasing dose. Chotsiri
    # 2019 did NOT estimate this term: it was added to the final model
    # only for the dose-optimisation simulations (supplementary methods,
    # 'In silico lumefantrine dose optimisation': 'Lumefantrine exhibit
    # dose-dependent absorption and a saturation model was implemented on
    # the relative bioavailability, according to previously published
    # results, in order to avoid bias in dose optimisation simulations'),
    # with theta_Dosage50 fixed to 3.86 mg/kg from Kloprogge 2018.
    # Because the term was bolted on AFTER estimation, CL/F and Vc/F
    # already absorb the standard-dose bioavailability; applying the
    # printed form F = theta_F * (1 - Dosage / (Dose50 + Dosage))
    # un-normalised would multiply every concentration by 0.236 at the
    # standard dose and is contradicted by the paper's own Table 3 and
    # Figure 4. It is therefore normalised to the per-subject standard
    # per-dose exposure dose_std / WT, which makes the factor exactly 1
    # for every child on the standard regimen -- so the base model
    # reproduces the published Table 3 secondary parameters unchanged --
    # while giving the intended saturation for alternative dose levels.
    # The normalisation cancels in every dose RATIO, so the published dose
    # comparison is unaffected: for the typical 9.62-kg child, doubling
    # the dose gives F x 0.567, i.e. the '43% decrease in bioavailability'
    # quoted by Simeon 2024. See vignette Errata.
    dose50   <- fixed(3.86) ; label("Per-dose exposure at 50% saturation of absorption (mg/kg)")             # Chotsiri 2019 supplementary methods: theta_Dosage50 fixed to 3.86 mg/kg, from Kloprogge 2018 Table 2
    dose_std <- fixed(120)  ; label("Standard per-dose lumefantrine amount at which the saturation term is inert (mg)") # Chotsiri 2019 Methods, 'Study design': 1 tablet (120 mg lumefantrine) per dose below 15 kg, the band covering essentially the whole cohort

    # ---- Pharmacodynamics: time to P. falciparum reinfection --------
    # Constant baseline hazard, reported per year and converted to the
    # model's hourly time base: 5.25 / (365.25 * 24) = 5.9890e-4 /h.
    lbase <- log(5.25 / (365.25 * 24)) ; label("Log baseline hazard of P. falciparum reinfection (1/h)")          # Chotsiri 2019 Table 2: BASE = 5.25 reinfections/year (RSE 11.8%; bootstrap median 5.58, 95% CI 4.19-6.85)

    # Sigmoid Emax protective effect of lumefantrine on the hazard
    # (supplementary methods): LFEFF = 1 - Emax * Cp^gamma /
    # (IC50^gamma + Cp^gamma). IC50 is published as 156 ng/mL and is
    # expressed here in ug/mL to match the model's concentration units.
    lic50 <- log(0.156) ; label("Log lumefantrine concentration giving 50% reduction of the reinfection hazard (ug/mL)") # Chotsiri 2019 Table 2: IC50 = 156 ng/mL (RSE 12.1%; bootstrap median 156, 95% CI 141-214)
    lhill <- log(4.77)  ; label("Log shape factor gamma of the sigmoid drug effect on the reinfection hazard (unitless)") # Chotsiri 2019 Table 2: gamma = 4.77 (RSE 38.8%; bootstrap median 4.90, 95% CI 2.30-9.78)

    # Emax is not tabulated in Chotsiri 2019 Table 2 and no value appears
    # anywhere in the paper or supplement. It is fixed to 1 here, which is
    # the only value consistent with the published parameterisation: with
    # Emax = 1 the protective effect LFEFF spans 1 (no drug) down to 0
    # (saturating concentrations) and IC50 is, as the supplement defines
    # it, 'the lumefantrine capillary blood concentration needed in order
    # to reduce the hazard of malaria reinfection by 50%'. Recorded in
    # vignette Errata as an unreported-parameter assumption.
    emax <- fixed(1) ; label("Maximum fractional reduction of the reinfection hazard (unitless)")                  # Not reported by Chotsiri 2019 -- fixed to 1; see vignette Errata

    # ---- Between-subject variability --------------------------------
    # Table 2 reports BSV as %CV with footnote a giving the log-normal
    # conversion 100 * sqrt(exp(omega^2) - 1), so the internal variances
    # are omega^2 = log((CV/100)^2 + 1). BSV on CL/F, Vc/F and Vp/F was
    # estimated close to zero and removed from the final model (Results,
    # 'Covariate model'), so those parameters carry no eta here.
    #   F    CV 64.0% -> log(0.640^2 + 1) = 0.343306
    #   MTT  CV 192%  -> log(1.92^2  + 1) = 1.544665
    #   Q/F  CV 67.8% -> log(0.678^2 + 1) = 0.378220
    etalfdepot ~ 0.343306  # Chotsiri 2019 Table 2: BSV on F = 64.0% CV (RSE 7.08%; bootstrap median 63.9%, 95% CI 52.2-76.9%). Enters F through the Box-Cox transformation, not as a plain log-normal eta -- see model().
    etalmtt    ~ 1.544665  # Chotsiri 2019 Table 2: BSV on MTT = 192% CV (RSE 8.88%; bootstrap median 192%, 95% CI 104-259%); eta shrinkage 53.2%
    etalq      ~ 0.378220  # Chotsiri 2019 Table 2: BSV on QP/F = 67.8% CV (RSE 7.69%; bootstrap median 67.8%, 95% CI 54.0-82.2%); eta shrinkage 24.4%

    # ---- Residual error ---------------------------------------------
    # Chotsiri 2019 modelled the natural logarithm of the capillary blood
    # concentration with an additive error on the log scale
    # (supplementary methods: 'Residual unexplained variability was
    # modelled as an additive error on the log-transformed observed
    # concentrations (equivalent to an exponential error on an arithmetic
    # scale)'), which maps to proportional error in nlmixr2's linear
    # concentration space. Table 2's footnote defines the tabulated
    # quantity as 'sigma, residual error VARIANCE of lumefantrine
    # concentrations', so the reported 0.339 is squared units and the
    # log-scale SD is sqrt(0.339) = 0.582237. See vignette Errata: the
    # sibling Kloprogge 2018 reports a numerically similar sigma = 0.323
    # that its own footnote calls an 'additive residual error on log
    # scale' (i.e. an SD), so the two papers' conventions differ.
    propSd <- 0.582237 ; label("Proportional residual SD on the linear concentration scale (= additive SD on the log scale)") # Chotsiri 2019 Table 2: sigma = 0.339 (RSE 5.20; bootstrap median 0.339, 95% CI 0.265-0.426), reported as a VARIANCE -> SD = sqrt(0.339)

    # Placeholder additive residual on the survivor-function output. The
    # source model is an interval-censored time-to-event likelihood, not
    # an additively-observed continuous endpoint; `sur` is exposed as a
    # derived output so the protective-efficacy simulations of the
    # vignette can be run, and this residual exists only so the nlmixr2
    # likelihood machinery accepts a second endpoint. Same convention as
    # Lindauer_2017_lacosamide_seizure.R.
    addSd <- 0.001 ; label("Placeholder additive residual error on the reinfection-free survival output `sur` (unitless); not from the source -- see vignette Assumptions")
  })

  model({
    # ---- Individual parameters --------------------------------------
    # Allometric body-weight scaling on all disposition parameters,
    # centered at the model-building median 9.62 kg, plus first-order
    # enzyme maturation on CL/F driven by postnatal age in months. The
    # maturation factor is NOT normalised to a reference age (see the
    # covariateData note on PNA): exp(lcl) is the fully mature value.
    mf  <- PNA^hill_cl / (tm50_cl^hill_cl + PNA^hill_cl)
    mtt <- exp(lmtt + etalmtt)
    cl  <- exp(lcl)         * (WT / 9.62)^e_wt_cl * mf
    vc  <- exp(lvc)         * (WT / 9.62)^e_wt_vc
    q   <- exp(lq + etalq)  * (WT / 9.62)^e_wt_cl
    vp  <- exp(lvp)         * (WT / 9.62)^e_wt_vc

    # Transit-absorption rate constant. Two transit compartments plus the
    # dosing compartment give three first-order transfers between the
    # dosing site and the central compartment, so the mean transit time
    # MTT = 3 / ktr and ka = ktr (Savic 2007).
    ktr <- 3 / mtt
    ka  <- ktr

    # Two-compartment disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- Relative bioavailability -----------------------------------
    # Three multiplicative factors on F (which is fixed to 1 for the
    # population):
    #
    # 1. Box-Cox-transformed IIV (Petersson 2009, as printed in the
    #    supplementary methods):
    #      eta_transformed = ((exp(eta_F))^lambda - 1) / lambda
    #      F_i             = theta_F * exp(eta_transformed)
    #    (exp(eta))^lambda is written here as exp(lambda * eta), which is
    #    the same quantity and is numerically better behaved.
    #
    # 2. Exponential MUAC effect, 25.4% lower F per 1 cm reduction. MUAC
    #    is carried in mm and converted to cm because the published
    #    coefficient is per cm.
    #
    # 3. Dose-saturable absorption, normalised so the factor is exactly 1
    #    for a child on the standard per-dose amount dose_std, whatever
    #    that child weighs. From the printed form
    #    g(D) = 1 - D / (dose50 + D) = dose50 / (dose50 + D), the
    #    normalised factor is g(D) / g(D_std) =
    #    (dose50 + D_std) / (dose50 + D), with D = DOSE / WT and
    #    D_std = dose_std / WT both in mg/kg. See the ini() comment on
    #    dose_std and the vignette Errata for why the normalisation is
    #    required.
    eta_fdepot_bc <- (exp(boxcox_fdepot * etalfdepot) - 1) / boxcox_fdepot
    fmuac         <- exp(e_muac_f * (MUAC / 10 - muac_ref))
    dose_mgkg     <- DOSE / WT
    fdose         <- (dose50 + dose_std / WT) / (dose50 + dose_mgkg)
    f(depot)      <- exp(lfdepot + eta_fdepot_bc) * fmuac * fdose

    # ---- ODE system -------------------------------------------------
    # Dose enters `depot`, the first compartment of the transit chain;
    # `transit1` and `transit2` are the two transit compartments of the
    # source model; absorption into `central` proceeds at ka = ktr.
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(central)     <-  ka  * transit2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                   k12 * central - k21 * peripheral1

    # Lumefantrine capillary blood concentration in ug/mL (dose in mg and
    # Vc in L give mg/L = ug/mL). The paper tabulates concentrations in
    # ng/mL; multiply by 1,000 to compare against the published day-7
    # concentrations, Cmax values and the 164-182 ng/mL clinical MIC.
    Cc <- central / vc

    # ---- Time to P. falciparum reinfection --------------------------
    # Sigmoid Emax protective effect and constant baseline hazard
    # (supplementary methods):
    #   LFEFF = 1 - Emax * Cp^gamma / (IC50^gamma + Cp^gamma)
    #   Hz(t) = theta_BASE * LFEFF
    #   S(t)  = exp(-integral_0^t Hz(t) dt)
    # The supplement prints the survivor function without its minus sign,
    # S(t) = exp(integral Hz dt), which would exceed 1 for all t > 0; the
    # standard (and only coherent) form is used here. Flagged in vignette
    # Errata.
    hill  <- exp(lhill)
    ic50  <- exp(lic50)
    lfeff <- 1 - emax * Cc^hill / (ic50^hill + Cc^hill)
    hazard <- exp(lbase) * lfeff

    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    sur          <- exp(-cumhaz)

    # Proportional residual error on the linear concentration scale
    # (NONMEM additive-on-log-scale maps to proportional here) plus the
    # placeholder residual that exposes the survivor function as a second
    # endpoint (see ini()).
    Cc  ~ prop(propSd)
    sur ~ add(addSd)
  })
}
