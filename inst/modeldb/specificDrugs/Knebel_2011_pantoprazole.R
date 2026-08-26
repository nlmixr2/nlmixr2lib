Knebel_2011_pantoprazole <- function() {
  description <- "Two-compartment population PK model with first-order absorption for pantoprazole in 202 pediatric patients from birth to 16 years, pooled across six clinical trials of an intravenous formulation, a delayed-release tablet, and a delayed-release granule (spheroid) formulation (Knebel 2011). Body-weight allometric scaling is fixed (0.75 on CL and Q, 1 on Vc and Vp, reference 10 kg). Allometrically scaled clearance carries a sigmoid Emax maturation function of postnatal age with the maximum effect fixed at 1, a Hill coefficient of 1.48, and an age at 50% of mature CL of 0.153 years in full-term infants and 1.38-fold higher (0.211 years) in preterm infants; three multiplicative categorical effects also act on CL (male sex 1.06, CYP2C19 poor metabolizer 0.0716, African American race 1.29). The tablet is the oral bioavailability anchor (F1 = 1) and the granule has F1 = 0.295 with 56.7% inter-occasion variability and a slower absorption rate constant (0.613 vs 1.32 per hour); both oral forms share a 0.444 hour absorption lag. Residual error is proportional-only for intravenous data and proportional plus a shared additive term for the two oral formulations. The reference subject is a female, full-term, 10 kg, extensive or unknown CYP2C19 metabolizer, non-African American patient receiving the intravenous or tablet formulation."
  reference   <- paste(
    "Knebel W, Tammara B, Udata C, Comer G, Gastonguay MR, Meng X.",
    "Population pharmacokinetic modeling of pantoprazole in pediatric",
    "patients from birth to 16 years.",
    "J Clin Pharmacol. 2011;51(3):333-345.",
    "doi:10.1177/0091270010366146.",
    sep = " "
  )
  vignette    <- "Knebel_2011_pantoprazole"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot       = list(analyte = "pantoprazole", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "pantoprazole", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "pantoprazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling fixed a priori at 0.75 on clearances (CL, Q) and",
        "1 on volumes (V2, V3), normalised to a reference weight of 10 kg",
        "(Knebel 2011 Methods equation 3 and Results equation 6). Cohort",
        "range 1.57-127 kg, median 7.93 kg, mean 19.2 kg (Table I).",
        "Body surface area was also investigated as the size descriptor but",
        "body weight was carried into the final model; see",
        "covariatesDataExcluded."
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Postnatal age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives the sigmoid Emax maturation function on allometrically",
        "scaled CL in Knebel 2011 equation 6:",
        "AGE^HILL / (AGE^HILL + AG50^HILL), with the maximum effect fixed at",
        "1 so the factor asymptotes to the mature (adult) allometrically",
        "scaled CL. AG50 is 0.153 years in full-term and 0.211 years in",
        "preterm infants. Because the ratio AGE / AG50 is what enters, AGE",
        "must be supplied in YEARS to match the published AG50. Cohort range",
        "0.025-16 years (0.3-192 weeks), median 0.65 years (Table I). This is",
        "postnatal, not postmenstrual, age -- the preterm adjustment is",
        "carried by TERM_BIRTH rather than by shifting the age axis."
      ),
      source_name        = "AGE"
    ),
    TERM_BIRTH = list(
      description        = "Term-vs-preterm birth indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (full-term birth) is the typical-value reference in Knebel 2011; the preterm stratum (0) multiplies AG50 by 1.38",
      notes              = paste(
        "Knebel 2011 defines preterm as gestational age less than 38 weeks",
        "(Results, Analysis Population). NOTE this is the 38-week cutoff used",
        "by the paper, not the 37-week cutoff named in the TERM_BIRTH",
        "canonical register entry -- assemble the column with the paper's",
        "38-week rule when driving this model. 77 of 202 patients were",
        "preterm, from 1 to 15 weeks preterm. Knebel 2011 estimated a",
        "separate AG50 for preterm infants because the base-model eta on CL",
        "showed a more prominent age effect in preterm infants (Figure 4).",
        "Enters as a multiplicative factor on AG50 raised to (1 -",
        "TERM_BIRTH), i.e. the factor applies only to preterm subjects."
      ),
      source_name        = "gestational age < 38 weeks (preterm indicator; canonical orientation is inverted)"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female) is the Knebel 2011 typical-value reference; the estimated 1.06 factor applies to males",
      notes              = paste(
        "Knebel 2011 Table II codes sex as 0 = female / 1 = male and the",
        "abstract names female as a reference covariate, so the published",
        "THETA15 = 1.06 is the MALE multiplier. The canonical SEXF column is",
        "1 = female, so the effect is applied as e_sexf_cl^(1 - SEXF) to",
        "preserve the verbatim published coefficient with the canonical",
        "column orientation (same construction as Bajaj_2017_nivolumab.R and",
        "Wada_2023_sparsentan.R). 81 females / 121 males (40% / 60%).",
        "Knebel 2011 classes this effect as clinically unimportant: the 95%",
        "CI (0.832, 1.34) sits inside the prespecified +/- 25% window."
      ),
      source_name        = "SEX (0 = female, 1 = male; SEXF = 1 - SEX)"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-African American: White, Asian, Hispanic, and Other pooled)",
      notes              = paste(
        "Knebel 2011 Table II codes ethnic origin as 1 = White, 2 = African",
        "American, 3 = Asian, 4 = Hispanic, 5 = Other, and equation 6 carries",
        "only the RACE2 (African American) contrast, so all four non-African",
        "American categories form the reference. 38 of 202 patients (19%)",
        "were African American. Knebel 2011 calls the effect indeterminate",
        "with respect to clinical relevance: the 95% CI (0.995, 1.63)",
        "straddles the +25% boundary."
      ),
      source_name        = "RACE (RACE2 = African American; RACE_BLACK = as.integer(RACE == 2))"
    ),
    CYP2C19_PM = list(
      description        = "CYP2C19 poor-metabolizer phenotype indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (extensive metabolizer OR unknown / not-determined phenotype, pooled)",
      notes              = paste(
        "Knebel 2011 Table III codes CYP2C19 phenotype as 0 = unknown / not",
        "determined, 1 = poor metabolizer, 2 = extensive metabolizer, and",
        "equation 6 carries only the CPH = 1 (poor metabolizer) contrast, so",
        "the reference pools extensive metabolizers with unknown phenotype",
        "(the abstract names 'extensive/unknown CYP2C19 metabolizer status'",
        "as a reference covariate). 4 poor / 165 extensive / 33 unknown.",
        "Three of the four poor metabolizers were preterm infants and all",
        "four were under 7 months old, so this estimate is confounded with",
        "the maturation term and rests on very few subjects."
      ),
      source_name        = "CPH (CYP2C19 phenotype; CYP2C19_PM = as.integer(CPH == 1))"
    ),
    FORM_TABLET = list(
      description        = "Delayed-release tablet formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not the tablet). Paired with FORM_GRANULE so that both = 0 selects the INTRAVENOUS route; see notes",
      notes              = paste(
        "Knebel 2011 pooled three administration modes: intravenous,",
        "delayed-release tablet, and delayed-release granule (called",
        "'spheroid' in the paper's equations). They are encoded here as the",
        "two-indicator decomposition FORM_TABLET / FORM_GRANULE with",
        "(0, 0) = intravenous, the same three-level construction used by",
        "Kleideiter_2017_cebranopadol.R and Wada_2023_sparsentan.R. NOTE the",
        "comparator here is the INTRAVENOUS route, not the non-tablet oral",
        "liquid named as the default reference in the FORM_TABLET register",
        "entry. Exactly one of {intravenous, FORM_TABLET, FORM_GRANULE} is",
        "selected per dose record. The tablet is the oral bioavailability",
        "anchor (F1 fixed to 1) and supplies the default absorption rate",
        "constant that the granule branch overrides. The formulation also",
        "selects the residual-error model: proportional-only for",
        "intravenous, proportional plus a shared additive term for both oral",
        "formulations."
      ),
      source_name        = "formulation ('tablet' branch of Knebel 2011 equation 6)"
    ),
    FORM_GRANULE = list(
      description        = "Delayed-release granule (spheroid) formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not the granule). Paired with FORM_TABLET so that both = 0 selects the INTRAVENOUS route; see FORM_TABLET notes",
      notes              = paste(
        "Knebel 2011's 'spheroid' formulation: delayed-release granules to be",
        "sprinkled on applesauce, mixed with apple juice, or blended into a",
        "suspension, developed for pediatric patients unable to swallow a",
        "tablet. It overrides the tablet defaults for both bioavailability",
        "(F1 = 0.295 vs 1) and absorption rate (Ka = 0.613 vs 1.32 per hour),",
        "and carries the only inter-occasion variability in the model (56.7%",
        "CV on F1). The median age of granule recipients was 0.3 years with",
        "116 of 137 under 1 year; the authors attribute the low F1 estimate",
        "partly to incomplete intake rather than to true formulation",
        "bioavailability, noting that an adult bioequivalence study found the",
        "granule bioequivalent to the tablet (Discussion). The comparator",
        "here is the intravenous route (see FORM_TABLET notes), not the",
        "tablet named as the default reference in the FORM_GRANULE register",
        "entry."
      ),
      source_name        = "formulation ('spheroid' branch of Knebel 2011 equation 6)"
    ),
    OCC = list(
      description        = "Occasion index for inter-occasion variability on granule bioavailability",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Knebel 2011 Table IV footnote: 'IOV-spheroid, interoccasion",
        "variability for spheroid (occasion defined as first dose vs all",
        "others)'. Two occasions only: OCC = 1 for the first granule dose,",
        "OCC = 2 for every later granule dose. Equation 6 writes the IOV as",
        "F1_spheroid = THETA8 * exp(eta_OCC1 + eta_OCC2), the NONMEM idiom in",
        "which exactly one of the two etas is active on any given occasion",
        "and both share one variance (NONMEM $OMEGA BLOCK(1) SAME).",
        "Decomposed inside model() into oc1 / oc2 indicators multiplexing",
        "etaiov_fdepot_1 / etaiov_fdepot_2, per the pattern in",
        "Chen_2023_nemonoxacin.R and Jonsson_2011_ethambutol.R. The IOV",
        "applies to the granule only; intravenous and tablet records may pass",
        "any OCC value because the term is multiplied by FORM_GRANULE."
      ),
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Knebel 2011 Covariate Model: 'a measure of body size (weight or body",
        "surface area [BSA]) on CL, central volume of distribution (V2),",
        "intercompartmental clearance (Q), and peripheral volume of",
        "distribution (V3)'. Body weight was the size descriptor carried into",
        "the final model, so BSA does not appear in equation 6 and no BSA",
        "coefficient is reported. Cohort range 0.13-2.5 m^2, median 0.39,",
        "mean 0.636 (Table I). Recorded here for provenance only."
      ),
      source_name        = "Body surface area"
    ),
    GA = list(
      description        = "Gestational age at birth",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reported for the 119 patients aged birth to 11 months (Table I:",
        "range 23-41 weeks, median 34, mean 33.6). Knebel 2011 does NOT use",
        "gestational age as a continuous covariate; it enters the final model",
        "only through the dichotomy 'preterm = gestational age < 38 weeks',",
        "which is carried by TERM_BIRTH. Recorded here so the derivation rule",
        "for TERM_BIRTH is auditable."
      ),
      source_name        = "Gestational age"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 202L,
    n_studies        = 6L,
    n_observations   = 922L,
    age_range        = "0.025-16 years (Table I row 'Age, y'); birth to 16 years",
    age_median       = "0.65 years (Table I row 'Age, y'; mean 3.69 years)",
    age_notes        = paste(
      "Table I reports two age rows for the same 202 patients and they are",
      "NOT unit conversions of each other: 'Age, y' = 0.025-16, median 0.65,",
      "mean 3.69, while 'Age, wks' = 0.3-192, median 7.8, mean 44.3. Sixteen",
      "years is 835 weeks, not 192, and 0.65 years is 33.9 weeks, not 7.8.",
      "Only the years row is consistent with the paper's own 'birth to 16",
      "years' framing and with the 1-to-16-year intravenous and tablet arms,",
      "so the years row is the one quoted here and is the row that",
      "corresponds to the AGE covariate driving the maturation function."
    ),
    weight_range     = "1.57-127 kg",
    weight_median    = "7.93 kg (mean 19.2 kg)",
    bsa_range        = "0.13-2.5 m^2 (median 0.39, mean 0.636)",
    ga_range         = "23-41 weeks in the 119 patients with gestational age recorded (median 34, mean 33.6)",
    sex_female_pct   = 40,
    race_ethnicity   = c(White = 72, Black = 19, Asian = 2, Hispanic = 2, Other = 5),
    disease_state    = paste(
      "Gastroesophageal reflux disease (GERD) in the four oral single- and",
      "multiple-dose trials; the two intravenous single-dose trials enrolled",
      "pediatric patients aged 1 to 16 years. 77 of 202 patients were preterm",
      "(gestational age < 38 weeks), from 1 to 15 weeks preterm. CYP2C19",
      "phenotype: 4 poor metabolizers, 165 extensive metabolizers, 33 unknown",
      "or not determined; three of the four poor metabolizers were preterm",
      "infants and all four were under 7 months old."
    ),
    dose_range       = paste(
      "Intravenous 0.8 or 1.6 mg/kg single dose (22 patients aged 1-16",
      "years); oral single and multiple doses of 1.25, 2.5, 20, or 40 mg",
      "fixed or 0.6 or 1.2 mg/kg (119 patients from birth to 11 months and",
      "61 patients aged 1-16 years)."
    ),
    regions          = "Not reported",
    nonmem_method    = "NONMEM VI (ICON Development Solutions), FOCE with interaction, installed and patch-tracked via NMQual 6.3",
    sampling_schema  = paste(
      "Plasma sampled at various times across the 24-hour dosing interval;",
      "922 quantifiable concentrations from 202 patients. Assay LC-MS/MS",
      "(AAI Pharma, Shawnee KS) with a lower limit of quantification of 10",
      "ng/mL. Observations below the quantification limit, missing values,",
      "and observations without a corresponding dosing time were excluded",
      "from the analysis (no BLQ likelihood method was applied)."
    ),
    notes            = paste(
      "Demographics from Knebel 2011 Tables I, II, and III; six pooled",
      "clinical trials (2 single-dose intravenous, 4 single- and",
      "multiple-dose oral). Model evaluated by a 500-replicate predictive",
      "check on the within-individual median concentration and a",
      "1000-replicate nonparametric bootstrap stratified by sex, age, CYP2C19",
      "phenotype, and race. The predictive check showed some overprediction",
      "of the tablet and underprediction of the granule, so the authors ran",
      "their dosing simulations from the post-hoc individual CL estimates",
      "rather than from a fresh Monte Carlo draw. Covariate effects were",
      "judged by whether the bootstrap 95% CI implied a change of more than",
      "+/- 25% of the typical value: only CYP2C19 poor-metabolizer status and",
      "age reached that bar."
    )
  )

  ini({
    # ===== Structural PK (Knebel 2011 Table IV, 'Final (Full) Model' column) =====
    # Reference subject: female, full term, extensive/unknown CYP2C19
    # metabolizer, non-African American, 10 kg, intravenous or tablet
    # (Knebel 2011 abstract). The typical CL below is the value at the
    # maturation asymptote; the age factor multiplies it down for infants.
    lcl <- log(1.93);  label("Typical mature allometrically scaled clearance CL at 10 kg (L/h)")      # Knebel 2011 Table IV: CL = 1.93 L/h, %SE 13, bootstrap 95% CI 1.53-2.61
    lvc <- log(1.3);   label("Typical central volume V2 at 10 kg (L)")                                # Knebel 2011 Table IV: V2 = 1.3 L, %SE 9, bootstrap 95% CI 0.925-1.56
    lq  <- log(0.23);  label("Typical intercompartmental clearance Q at 10 kg (L/h)")                 # Knebel 2011 Table IV: Q = 0.23 L/h, %SE 23, bootstrap 95% CI 0.155-0.953
    lvp <- log(0.596); label("Typical peripheral volume V3 at 10 kg (L)")                             # Knebel 2011 Table IV: V3 = 0.596 L, %SE 31, bootstrap 95% CI 0.297-0.974

    # ===== Absorption, per oral formulation (Knebel 2011 Table IV) =====
    # Equation 6 is a NONMEM $PK default-then-override block: Ka, ALAG1 and
    # F1 are set to the tablet values first, and the If(spheroid) branch
    # overrides Ka and F1 only. ALAG1 is NOT overridden, so the lag applies
    # to both oral formulations. Both formulations share one eta on Ka.
    lka_tablet  <- log(1.32);       label("Typical absorption rate constant Ka for the delayed-release tablet (1/h)")   # Knebel 2011 Table IV: Ka tablet = 1.32 1/h, %SE 9, bootstrap 95% CI 1.05-1.92
    lka_granule <- log(0.613);      label("Typical absorption rate constant Ka for the delayed-release granule (1/h)")  # Knebel 2011 Table IV: Ka spheroid = 0.613 1/h, %SE 18, bootstrap 95% CI 0.428-1.40
    ltlag       <- log(0.444);      label("Absorption lag time for oral dosing (h)")                                    # Knebel 2011 Table IV: Lag Time = 0.444 h, %SE 3, bootstrap 95% CI 0.400-0.491

    # F1 is relative to intravenous. The tablet is the anchor and was held
    # at 1; the granule bioavailability was estimated against it.
    lfdepot_tablet  <- fixed(log(1)); label("Relative bioavailability of the delayed-release tablet (unitless)")   # Knebel 2011 Table IV: F1 tablet = 1 (fixed), NA bootstrap CI
    lfdepot_granule <- log(0.295);    label("Relative bioavailability of the delayed-release granule (unitless)")  # Knebel 2011 Table IV: F1 spheroid = 0.295, %SE 17, bootstrap 95% CI 0.175-0.405

    # ===== Allometric exponents (Knebel 2011 equation 3, fixed) =====
    # "theta_allo is a fixed allometric power parameter, which was assigned a
    # value of 0.75 for physiologic parameters, such as clearances, and a
    # value of 1 for anatomical volumes." Reference weight WTref = 10 kg.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT/10) for CL (unitless)") # Knebel 2011 Methods equation 3; Table IV row '(WT/10)^THETA10' = 0.75 (fixed)
    e_wt_q  <- fixed(0.75); label("Allometric exponent on (WT/10) for Q (unitless)")  # Knebel 2011 Methods equation 3; Table IV row '(WT/10)^THETA10' = 0.75 (fixed)
    e_wt_vc <- fixed(1);    label("Allometric exponent on (WT/10) for V2 (unitless)") # Knebel 2011 Methods equation 3; Table IV row '(WT/10)^THETA11' = 1 (fixed)
    e_wt_vp <- fixed(1);    label("Allometric exponent on (WT/10) for V3 (unitless)") # Knebel 2011 Methods equation 3; Table IV row '(WT/10)^THETA11' = 1 (fixed)

    # ===== Age maturation of allometrically scaled CL (Knebel 2011 equation 6) =====
    # Sigmoid Emax in postnatal age with Emax fixed at 1:
    #   AGE^HILL / (AGE^HILL + AG50^HILL)
    # AG50 is the full-term value; preterm infants get AG50 * 1.38.
    # The Emax = 1 fixing is a structural choice, not an estimated value, so
    # it is written as the literal 1 in the sigmoid rather than as a
    # parameter -- Table IV reports no Emax row.
    mat_hill              <- 1.48;  label("Hill coefficient of the age maturation function on CL (unitless)")                          # Knebel 2011 Table IV: Hill_CL = 1.48, %SE 13, bootstrap 95% CI 0.979-1.90
    mat_age50             <- 0.153; label("Postnatal age at 50% of mature allometrically scaled CL, full-term birth (years)")          # Knebel 2011 Table IV: AG50 = 0.153 year, %SE 32, bootstrap 95% CI 0.0896-0.554
    mat_age50_preterm_ratio <- 1.38; label("Multiplicative factor applied to the age at 50% mature CL for preterm birth (unitless)")   # Knebel 2011 Table IV: AG50P_preterm = 1.38, %SE 24, bootstrap 95% CI 0.805-2.08; Discussion: 0.153 * 1.38 = 0.211 year for preterm

    # ===== Categorical covariate effects on CL (Knebel 2011 equation 6) =====
    # Equation 4 form: a multiplicative factor raised to the 0/1 indicator,
    # i.e. the factor when the indicator is 1 and 1 when it is 0.
    e_sexf_cl       <- 1.06;   label("Multiplicative factor on CL for male sex, relative to the female reference (unitless)")               # Knebel 2011 Table IV: THETA15_SEX = 1.06, %SE 12, bootstrap 95% CI 0.832-1.34; applied to (1 - SEXF) because the paper's reference is female
    e_cyp2c19_pm_cl <- 0.0716; label("Multiplicative factor on CL for CYP2C19 poor metabolizers (unitless)")                                # Knebel 2011 Table IV: THETA16_CPH1 = 0.0716, %SE 41, bootstrap 95% CI 0.0274-0.199
    e_race_black_cl <- 1.29;   label("Multiplicative factor on CL for African American race (unitless)")                                    # Knebel 2011 Table IV: THETA17_RACE2 = 1.29, %SE 12, bootstrap 95% CI 0.995-1.63

    # ===== Inter-individual variability (Knebel 2011 Table IV, final model) =====
    # Table IV reports OMEGA on the VARIANCE scale: the printed 'CV%' equals
    # the square root of the printed variance (e.g. sqrt(0.412) = 0.642 =
    # 64.2%), and the CL-V2 covariance reproduces the printed correlation
    # (0.0898 / sqrt(0.412 * 0.25) = 0.28). So the numbers below are used
    # verbatim, with no CV-to-variance conversion.
    etalcl + etalvc ~ c(0.412,
                        0.0898, 0.25)  # Knebel 2011 Table IV: OMEGA1.1 CL = 0.412 (%SE 18, CV% 64.2, 95% CI 0.242-0.573); OMEGA1.2 COV CL-V2 = 0.0898 (%SE 115, r = 0.28, 95% CI -0.234-0.292); OMEGA2.2 V2 = 0.25 (%SE 31, CV% 50, 95% CI 0.113-0.897)
    etalka          ~ 0.586            # Knebel 2011 Table IV: OMEGA3.3 Ka = 0.586 (%SE 32, CV% 76.5, 95% CI 1.3e-11-1.42); one eta shared by the tablet and granule Ka per equation 6

    # ===== Inter-occasion variability on granule F1 (Knebel 2011 equation 6) =====
    # F1_spheroid = THETA8 * exp(eta_OCC1 + eta_OCC2) with two occasions
    # (first dose vs all others) sharing one variance. Expanded into two
    # occasion-multiplexed etas because rxode2 cannot simulate the
    # `eta ~ var | OCC` form; the second is fixed to the first's variance,
    # which is what NONMEM `$OMEGA BLOCK(1) SAME` encodes.
    etaiov_fdepot_1 ~ 0.321         # Knebel 2011 Table IV: OMEGA4.4 F1 IOV-spheroid = 0.321 (%SE 27, CV% 56.7, 95% CI 0.142-0.519)
    etaiov_fdepot_2 ~ fixed(0.321)  # SAME-equivalent: equal to the occasion-1 IOV variance

    # ===== Residual error (Knebel 2011 equation 6 and Table IV) =====
    # Equation 6 gives a different residual model per formulation:
    #   If (IV)       C = Chat * (1 + eps1_p)
    #   If (tablet)   C = Chat * (1 + eps2_p) + eps3_a
    #   If (spheroid) C = Chat * (1 + eps4_p) + eps3_a
    # so the two oral formulations SHARE one additive term (Table IV names
    # the row 'addTAB-SPH') and each formulation has its own proportional
    # term. Table IV reports SIGMA on the variance scale with the SD / CV%
    # alongside; the SDs below are the printed SD / CV% values, which equal
    # the square root of the printed variance in every row.
    propSd_iv      <- 0.260; label("Proportional residual SD for intravenous observations (fraction)")   # Knebel 2011 Table IV: SIGMA1.1 proIV = 0.0678 (%SE 23), CV% = 26.0 = sqrt(0.0678); 95% CI 0.0275-0.101
    propSd_tablet  <- 0.586; label("Proportional residual SD for tablet observations (fraction)")        # Knebel 2011 Table IV: SIGMA2.2 proTAB = 0.344 (%SE 17), CV% = 58.6 = sqrt(0.344); 95% CI 0.241-0.476
    propSd_granule <- 0.560; label("Proportional residual SD for granule observations (fraction)")       # Knebel 2011 Table IV: SIGMA4.4 proSPH = 0.314 (%SE 10), CV% = 56.0 = sqrt(0.314); 95% CI 0.244-0.377
    addSd_oral     <- 6.08;  label("Additive residual SD shared by the tablet and granule observations (ng/mL)") # Knebel 2011 Table IV: SIGMA3.3 addTAB-SPH = 37 (%SE 58), SD = 6.08 = sqrt(37); 95% CI 8.27-1940
  })

  model({
    # ----- 1. Derived covariate terms -----
    # Three-level formulation decomposition: (0, 0) selects intravenous.
    is_iv <- (1 - FORM_TABLET) * (1 - FORM_GRANULE)

    # Occasion indicators multiplexing the granule-F1 IOV etas (OCC = 1 is
    # the first granule dose, OCC = 2 every later granule dose).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2

    # Sigmoid Emax maturation of allometrically scaled CL, Emax fixed at 1.
    # AG50 is inflated 1.38-fold in preterm infants (TERM_BIRTH == 0).
    mat_age50_i   <- mat_age50 * mat_age50_preterm_ratio^(1 - TERM_BIRTH)
    maturation_cl <- AGE^mat_hill / (AGE^mat_hill + mat_age50_i^mat_hill)

    # Categorical CL multipliers. The published sex factor is the MALE
    # effect, so it is raised to (1 - SEXF) against the canonical
    # 1 = female column.
    cov_cl <- e_sexf_cl^(1 - SEXF) *
              e_cyp2c19_pm_cl^CYP2C19_PM *
              e_race_black_cl^RACE_BLACK

    # ----- 2. Individual PK parameters -----
    cl <- exp(lcl + etalcl) * (WT / 10)^e_wt_cl * maturation_cl * cov_cl
    vc <- exp(lvc + etalvc) * (WT / 10)^e_wt_vc
    q  <- exp(lq)           * (WT / 10)^e_wt_q
    vp <- exp(lvp)          * (WT / 10)^e_wt_vp

    # Tablet supplies the default Ka; the granule branch overrides it.
    ka   <- exp(lka_tablet * (1 - FORM_GRANULE) + lka_granule * FORM_GRANULE + etalka)
    tlag <- exp(ltlag)

    # ----- 3. Micro-constants -----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ----- 4. ODE system -----
    # Oral doses (tablet, granule) enter `depot`; intravenous doses are
    # given directly into `central`.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # ----- 5. Bioavailability and lag time -----
    # Tablet F1 is the fixed anchor of 1; the granule overrides it with
    # 0.295 and carries the inter-occasion variability.
    f(depot)    <- exp(lfdepot_tablet * (1 - FORM_GRANULE) +
                       (lfdepot_granule + iov_fdepot) * FORM_GRANULE)
    alag(depot) <- tlag

    # ----- 6. Observation and error -----
    # Amounts are in mg and volumes in L, so central/vc is mg/L; the
    # factor of 1000 converts to the paper's ng/mL assay scale.
    Cc <- 1000 * central / vc

    # Formulation-specific residual error per equation 6: proportional-only
    # for intravenous, proportional plus the shared additive term for both
    # oral formulations.
    propSd <- propSd_iv * is_iv +
              propSd_tablet * FORM_TABLET +
              propSd_granule * FORM_GRANULE
    addSd  <- addSd_oral * (1 - is_iv)
    Cc ~ prop(propSd) + add(addSd)
  })
}
