Marier_2014_dutogliptin <- function() {
  description <- "Two-compartment population PK model with a lag time and formulation-specific first-order absorption for dutogliptin, a selective dipeptidyl peptidase-4 inhibitor, in healthy subjects and patients with type 2 diabetes mellitus (Marier 2014)"
  reference <- "Marier JF, Mouksassi MS, Gosselin NH, Li J. Population Pharmacokinetic Analysis of Dutogliptin, a Selective Dipeptidyl Peptidase-4 Inhibitor. Clin Pharmacol Drug Dev. 2014;3(4):297-304. doi:10.1002/cpdd.87"
  vignette <- "Marier_2014_dutogliptin"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "dutogliptin", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "dutogliptin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dutogliptin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance calculated with the Cockcroft-Gault formula; NOT BSA-normalized",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Raw Cockcroft-Gault CrCL in mL/min (Marier 2014 Methods, 'Population PK Modeling'), not the",
        "BSA-normalized mL/min/1.73 m^2 form that this canonical column more commonly carries.",
        "Power effect on CL/F with reference 115.7 mL/min. Marier 2014 Results, 'Subject Demographics':",
        "'The upper range of CrCL was truncated to 150 mL/min for the covariate analysis in order to",
        "better represent upper physiological range of CrCL' -- the model reproduces that truncation,",
        "so CRCL values above 150 mL/min have no further effect on CL/F.",
        "Cohort mean 117.5 mL/min, median 111.0, range 25.4-329.3 (Table 1)."
      ),
      source_name        = "CrCL"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect with reference 82.5 kg (the cohort median, Table 1) and exponent 1.0 on the",
        "apparent steady-state volume Vss/F (Table 2), applied here to Vc/F and Vp/F together so that",
        "the printed Table 2 equation Vss/F = 2041 * (WT/82.5)^1.0 is reproduced exactly.",
        "No weight effect on CL/F or Q/F was retained."
      ),
      source_name        = "Body Weight"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Caucasian or Asian; the two groups share one typical Vss/F in Marier 2014)",
      notes              = paste(
        "Multiplicative effect on Vss/F. Marier 2014 Results p.302-303 reports Vss/F = 2041 L for Caucasian",
        "and Asian subjects, 2432 L for Black subjects and 1998 L for other races, so the reference",
        "category pools Caucasian and Asian rather than being Caucasian alone.",
        "Paired with RACE_OTHER; both = 0 selects the Caucasian/Asian reference.",
        "Cohort composition (Table 1): Caucasian 63.99%, Asian 12.66%, Black 3.74%, Other 19.61%."
      ),
      source_name        = "Race"
    ),
    RACE_OTHER = list(
      description        = "Race category 'Other' indicator (American Indian or Alaskan Native, Native Hawaiian or other Pacific Islander, American Hispanic, Mixed Race)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Caucasian or Asian)",
      notes              = paste(
        "Multiplicative effect on Vss/F (1998 L vs the 2041 L Caucasian/Asian reference).",
        "The constituent groups are enumerated in Marier 2014 Methods, 'Population PK Modeling'.",
        "Paired with RACE_BLACK; both = 0 selects the Caucasian/Asian reference."
      ),
      source_name        = "Race"
    ),
    FORM_TABLET = list(
      description        = "Dutogliptin tartrate-salt tablet formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pooled uncoated free-base capsule and aqueous solution)",
      notes              = paste(
        "Per-dose covariate selecting the tablet first-order absorption rate constant",
        "(Ka = 0.558 1/h) instead of the pooled capsule / aqueous-solution reference (1.55 1/h),",
        "per Marier 2014 Table 2. Formulation acts on Ka only: Marier 2014 Results p.300 states that",
        "'no formulation or fed/fasting effects were observed on ALAG and Frel parameters'.",
        "The non-tablet comparator here is a solid capsule pooled with an aqueous solution, not a",
        "liquid alone. Paired with FORM_DUTOGLIPTIN_ECCAPSULE."
      ),
      source_name        = "Formulation"
    ),
    FORM_DUTOGLIPTIN_ECCAPSULE = list(
      description        = "Enteric-coated dutogliptin capsule formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pooled uncoated free-base capsule and aqueous solution)",
      notes              = paste(
        "Per-dose covariate selecting the enteric-coated-capsule first-order absorption rate constant",
        "(Ka = 0.0204 1/h), per Marier 2014 Results p.302. Given in period 4 of the PROT103 crossover",
        "to 6 subjects only (Supplemental Table S1); the paper's bootstrap median for this Ka was",
        "16.3% above the point estimate because of that small n.",
        "Paired with FORM_TABLET; both = 0 selects the capsule / aqueous-solution reference."
      ),
      source_name        = "Formulation"
    )
  )

  # Screened in the Marier 2014 covariate analysis but NOT retained in the
  # final model, so they are documented rather than referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Evaluated as a power-function covariate on CL/F and Vc/F (Marier 2014 Results, 'Covariates Analysis'); not retained by forward selection at P = .01. Cohort mean 52.7 years, range 18-77 (Table 1); 52 subjects were older than 65 years."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Correlated with CL/F (r = 0.22) in the exploratory screen but explicitly NOT carried into the formal covariate analysis: 'The effect of BMI on PK parameters of dutogliptin was not evaluated due to its strong correlation with body weight.' Cohort mean 31.1 kg/m^2 (Table 1)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Evaluated on CL/F and Vc/F; 'PK parameters of dutogliptin were not influenced by sex'. Cohort 49.73% female (Table 1)."
    ),
    DIS_DIAB = list(
      description = "Type 2 diabetes mellitus disease-status indicator (1 = T2DM patient, 0 = healthy subject)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Disease status was evaluated as a covariate; 'No differences were observed between healthy subjects and patients with T2DM' (Abstract). 144 Phase I subjects were healthy except PROT107 (renal impairment) and PROT109 (T2DM); all 417 Phase II subjects had T2DM."
    ),
    CONMED_METFORMIN = list(
      description = "Concomitant metformin co-administration indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Evaluated on CL/F and Vc/F in the drug-drug-interaction screen (PROT109 crossover); 'The coadministration of metformin did not affect any PK parameters of dutogliptin'. Not retained."
    ),
    FASTED_STRICT = list(
      description = "Fed / fasting status at the time of dosing",
      units       = "(binary)",
      type        = "binary",
      notes       = "Three fed/fasting conditions were tested on Ka, ALAG and Frel (administration with breakfast, breakfast/meal 1 hour after administration, meal 3-4 hours after administration); 'no formulation or fed/fasting effects were observed on ALAG and Frel parameters' and no fed/fasting effect on Ka was retained in the final model."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 561,
    n_studies      = 9,
    age_range      = "18-77 years",
    age_median     = "54 years",
    weight_range   = "44.6-158.6 kg",
    weight_median  = "82.5 kg",
    sex_female_pct = 49.73,
    race_ethnicity = c(Caucasian = 63.99, Asian = 12.66, Black = 3.74, Other = 19.61),
    disease_state  = "Healthy subjects and patients with type 2 diabetes mellitus; includes a dedicated renal-impairment study (PROT107, normal to severe impairment) and a metformin drug-drug-interaction study (PROT109).",
    dose_range     = "50-500 mg orally, single and once-daily repeated doses",
    renal_function = "Creatinine clearance (Cockcroft-Gault) mean 117.5 mL/min, median 111.0, range 25.4-329.3 mL/min; truncated at 150 mL/min for the covariate analysis.",
    co_medication  = "Phase II subjects were on stable metformin, or metformin plus a thiazolidinedione.",
    formulations   = "Aqueous solution, free-base hard-shell capsule, enteric-coated capsule, and tartrate-salt tablet.",
    notes          = "Demographics from Marier 2014 Table 1 (pooled Phase I and Phase II). 144 subjects contributed 2,933 measurable concentrations in 7 Phase I studies (PROT101, PROT102, PROT103, PROT104, PROT107, PROT109, PROT110; Supplemental Table S1); 417 subjects contributed 2,153 measurable concentrations in 2 Phase II studies (PROT201, PROT202; Supplemental Table S2). Assay LOQ 1 ng/mL. Model fitted with NONMEM VI."
  )

  ini({
    # ---- Absorption -------------------------------------------------------
    # Formulation-specific first-order absorption rate constants. Marier 2014
    # Results, 'Covariates Analysis': "Statistically significant formulations
    # effects were observed on the Ka of dutogliptin for all formulation tested."
    lka_cap <- log(1.55);   label("Absorption rate constant, free-base capsule and aqueous solution (Ka, 1/h)")  # Table 2 "1.55 for capsules"; Results p.302 "The Ka of dutogliptin following oral administrations of aqueous solutions was comparable to that of the capsules" -- no separate solution value is printed, so the solution shares the capsule Ka
    lka_tab <- log(0.558);  label("Absorption rate constant, tartrate-salt tablet (Ka, 1/h)")                    # Table 2 "0.558 for tablets"
    lka_ec  <- log(0.0204); label("Absorption rate constant, enteric-coated capsule (Ka, 1/h)")                  # Results p.302 "the enteric-coated capsules displayed much slower absorption (0.0204 hour^-1)"

    # Absorption lag time. Marier 2014 names ALAG three times as a component of
    # the final model (Abstract "a two-compartment model with a first-order rate
    # constant of absorption (Ka) and a lag time"; Methods p.299; Results p.302
    # "Frel and ALAG of dutogliptin for the capsule and tablet formulations were
    # similar") and reports that no formulation or fed/fasting effect acted on
    # it, but NEVER prints its value -- not in Table 2, not in the text, and not
    # in the supplement, which contains no parameter table. No printed quantity
    # depends on it, so it is not back-solvable, and Figures 1-2 have 6-12 h
    # x-axis ticks, far too coarse to read a sub-hour lag from. Encoded as
    # fixed(0) per the standing rule that an unreported value is never replaced
    # by a class-typical or extractor-chosen placeholder; the parameter is kept
    # named and present so the published structure stays visible and a user can
    # perturb it. Consequence: simulated Tmax is earlier than the published
    # model's by the (small) unreported lag. See the vignette Errata.
    tlag <- fixed(0);       label("Absorption lag time (ALAG, h) at 0; value not reported in Marier 2014")  # NOT REPORTED anywhere in the paper or supplement; operator ruling 2026-08-26 (sidecar request-001 q1 = A)

    # ---- Disposition ------------------------------------------------------
    lcl <- log(176);        label("Apparent clearance at CrCL 115.7 mL/min (CL/F, L/h)")  # Table 2: CL/F = 176 * (CrCL/115.7)^0.848
    e_crcl_cl <- 0.848;     label("Power exponent of creatinine clearance on CL/F (unitless)")  # Table 2 exponent; reproduces the printed CL/F of 121 L/h at CrCL 75 (121.9) and 79 L/h at CrCL 45 (79.0), Results p.302

    # Vss/F = 2041 L is the only distribution parameter Marier 2014 prints
    # (Table 2 and Results p.302-303, which defines it as the "sum of central and
    # peripheral volume of distribution"). Table S3 nonetheless carries SEPARATE
    # etas on Vc/F and Vp/F, so the fitted model had distinct typical values.
    # The split below is BACK-SOLVED, not paper-printed: the three published
    # (CL/F, terminal t1/2) pairs -- (176, 12.2 h), (121, 15.4 h), (79, 21.3 h)
    # from Results p.302-303 -- together with Vc + Vp = 2041 L over-determine
    # the two unknowns (Vc, Q) and reproduce all three half-lives to four
    # significant figures (12.199 / 15.430 / 21.290 h). Sweeping the full
    # 3-significant-figure rounding envelope of the printed half-lives admits
    # only Vc/F 1167-1306 L, Vp/F 735-874 L, Q/F 70.2-91.6 L/h, so the solve is
    # tightly bounded. Operator ruling 2026-08-26 (sidecar request-001 q2 = A).
    lvc <- log(1250);       label("Apparent central volume at 82.5 kg, Caucasian/Asian (Vc/F, L)")            # BACK-SOLVED from Table 2 Vss/F = 2041 L plus the three printed terminal half-lives (Results p.302-303); not printed in Marier 2014
    lvp <- log(791);        label("Apparent peripheral volume at 82.5 kg, Caucasian/Asian (Vp/F, L)")         # BACK-SOLVED; 2041 - 1250 = 791 L preserves the printed Vss/F exactly
    lq  <- log(78.6);       label("Apparent intercompartmental clearance (Q/F, L/h)")                         # BACK-SOLVED from the same three (CL/F, t1/2) pairs
    e_wt_vc_vp <- 1.0;          label("Power exponent of body weight on Vss/F (unitless)")                        # Table 2: Vss/F = 2041 * (Body Weight/82.5)^1.0

    # Race effects on Vss/F, expressed as fractional changes from the pooled
    # Caucasian/Asian reference. Results p.302-303: "The Vss/F of dutogliptin (sum
    # of central and peripheral volume of distribution) in Caucasian and Asian
    # (2,041 L) were consistent with those derived in Black (2,432 L) and other
    # races (1,998 L)."
    e_race_black_vc_vp <- 0.1916;   label("Fractional effect of Black race on Vss/F (unitless)")   # Results p.302-303: 2432 / 2041 - 1 = 0.19157
    e_race_other_vc_vp <- -0.02107; label("Fractional effect of Other race on Vss/F (unitless)")   # Results p.302-303: 1998 / 2041 - 1 = -0.021068

    # ---- Between-subject variability --------------------------------------
    # Supplemental Table S3 is the authoritative NONMEM OMEGA variance-
    # covariance matrix. Its diagonal is a log-scale VARIANCE: Table 2's BSV
    # column is exactly its square root for CL/F (sqrt(0.238) = 48.8%) and for
    # Ka (sqrt(0.920) = 95.9%). Results p.302: "An Omega block was used for
    # CL/F and Vc/F (refer to Table S3)."
    etalcl + etalvc ~ c(0.238,
                        0.253, 0.414)  # Table S3: var(CL/F) 0.238, cov(CL/F, Vc/F) 0.253, var(Vc/F) 0.414; correlation 0.806, determinant 0.0345 (positive definite)
    etalvp ~ 0.463   # Table S3: var(Vp/F), off-diagonals with CL/F, Vc/F and Ka all zero
    etalka ~ 0.920   # Table S3: var(Ka); matches Table 2 BSV 95.9% = sqrt(0.920)

    # ---- Residual error ---------------------------------------------------
    addSd  <- 1.13;  label("Additive residual error (ng/mL)")        # Table 2 "Additive error 1.13 ng/mL"
    propSd <- 0.277; label("Proportional residual error (fraction)") # Table 2 "Proportional error 27.7%"
  })
  model({
    # Creatinine clearance truncated at 150 mL/min, per Marier 2014 Results,
    # 'Subject Demographics'. Written as a weighted sum rather than min() so it
    # parses identically across rxode2 back-ends.
    crcl_cov <- CRCL * (CRCL <= 150) + 150 * (CRCL > 150)

    # Formulation selector. The two indicators are mutually exclusive; both
    # zero selects the pooled free-base-capsule / aqueous-solution reference.
    form_ref <- 1 - FORM_TABLET - FORM_DUTOGLIPTIN_ECCAPSULE

    # Race multiplier on Vss/F (Caucasian/Asian reference = 1).
    race_v <- 1 + e_race_black_vc_vp * RACE_BLACK + e_race_other_vc_vp * RACE_OTHER

    # Allometric-style weight scaling of the total volume of distribution.
    wt_v <- (WT / 82.5)^e_wt_vc_vp

    # Individual parameters. Ka is selected on the linear scale so the single
    # published Ka IIV (Table S3) applies to whichever formulation was given.
    ka <- (exp(lka_cap) * form_ref +
           exp(lka_tab) * FORM_TABLET +
           exp(lka_ec)  * FORM_DUTOGLIPTIN_ECCAPSULE) * exp(etalka)

    cl <- exp(lcl + etalcl) * (crcl_cov / 115.7)^e_crcl_cl
    vc <- exp(lvc + etalvc) * wt_v * race_v
    vp <- exp(lvp + etalvp) * wt_v * race_v
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Absorption lag time (see the ini() note: value never reported, fixed at 0).
    alag_d <- tlag
    alag(depot) <- alag_d

    # Dose in mg, volumes in L -> mg/L = ug/mL; x 1000 gives the ng/mL scale
    # that Marier 2014 reports throughout (assay LOQ 1 ng/mL, additive residual
    # error 1.13 ng/mL, average steady-state concentrations 47.3-210 ng/mL).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
