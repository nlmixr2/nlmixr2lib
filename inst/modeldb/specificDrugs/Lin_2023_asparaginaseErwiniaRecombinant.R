Lin_2023_asparaginaseErwiniaRecombinant <- function() {
  description <- "One-compartment population PK model for intramuscular recombinant Erwinia chrysanthemi asparaginase (JZP458, marketed as Rylaze) in pediatric and young-adult patients with acute lymphoblastic leukemia or lymphoblastic lymphoma who developed hypersensitivity to E. coli-derived asparaginases (Lin 2023, phase II/III study AALL1931). The measured quantity is serum asparaginase activity (SAA), so all amounts are activity units (IU) rather than mass. Absorption is mixed-order: the dose enters the depot as a zero-order input at a constant rate R1, running simultaneously with first-order absorption Ka out of the depot, which makes the terminal phase absorption rate limited (flip-flop). Ka, R1 and the relative bioavailability F were fixed to values carried over from the phase I intensive-sampling PopPK model (Lin 2021) because the sparse AALL1931 sampling could not characterize the absorption phase. Body surface area is an allometric covariate on both clearance and volume; Black or African American race and T-cell ALL disease subtype are multiplicative fractional-change covariates on clearance. Interindividual variability is exponential on clearance and volume, with a combined proportional and additive residual error."
  reference   <- "Lin T, Whigham T, Fernando I, Choi MR, Wang Q, Silverman JA. Population pharmacokinetics of intramuscular recombinant Erwinia chrysanthemi asparaginase (JZP458) in patients with acute lymphoblastic leukemia. Clin Transl Sci. 2023;16(5):898-909. doi:10.1111/cts.13499"
  vignette    <- "Lin_2023_asparaginaseErwiniaRecombinant"
  units       <- list(time = "h", dosing = "IU", concentration = "IU/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Amounts are asparaginase ACTIVITY units (IU), not mass:
  # the assay reads serum asparaginase activity and the model's zero-order
  # absorption rate R1 is reported in IU/h (Lin 2023 Table 2; Lin 2021 Table 2).
  compartmentData <- list(
    depot   = list(analyte = "recombinant Erwinia chrysanthemi asparaginase (JZP458)", units = "IU", specimen = "administration site", verified = TRUE),
    central = list(analyte = "recombinant Erwinia chrysanthemi asparaginase (JZP458)", units = "IU", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric (power) covariate on BOTH clearance (exponent 1.48) and central volume (exponent 1.61), normalized to a 1.2 m^2 reference. Lin 2023 Results, Base model: 'The BSA standard used for scaling was 1.2 m^2 to reflect an average pediatric patient.' The printed final-model equations (Lin 2023 Results, Covariate analysis and final PopPK model selection) are CL [mL/h] = 146 * (BSA/1.2)^1.48 * 0.674^(Black/African American) * 0.771^(T-ALL) and V [mL] = 445 * (BSA/1.2)^1.61. Analysis-set BSA 0.44-2.53 m^2 (mean 1.23, median 1.17; Lin 2023 Table 1). Note both exponents exceed 1, so clearance and volume rise faster than proportionally with BSA; because dose is also BSA-proportional (mg/m^2), trough activity DECREASES with increasing BSA. The dose-regimen simulations spanned BSA 0.25-3.03 m^2 (NHANES virtual population).",
      source_name        = "BSA"
    ),
    RACE_BLACK = list(
      description        = "Black or African American race indicator (1 = Black/African American, 0 = other)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Black/African American; in the AALL1931 analysis set this pools White/Caucasian 68.7%, Declined to state 11.4%, Asian 4.2%, American Indian/Alaska Native 1.8% and Other 0.6%)",
      notes              = "Time-fixed per subject. Multiplies clearance by 0.674, i.e. a 32.6% LOWER SAA clearance in Black/African American patients (Lin 2023 Table 2, RSE 10.6%, 95% CI 0.534-0.814; Discussion: 'modeled as fractional changes of 0.674'). 22 of 166 patients (13.3%; Lin 2023 Table 1). Race was carried into the model-based simulations deliberately: Lin 2023 Methods, Model-based simulations, notes the NHANES virtual population was constructed so its Black/African American proportion matched real-world ALL epidemiology precisely because race is a significant covariate on JZP458 PK. Despite the covariate being retained, the paper's subgroup simulations concluded no dosage modification is recommended based on race.",
      source_name        = "Race (Black/African American)"
    ),
    DIS_TALL = list(
      description        = "T-cell acute lymphoblastic leukemia disease-subtype indicator (1 = T-ALL, 0 = other ALL/LBL subtype)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-T-ALL; in the AALL1931 analysis set this pools B-ALL 74.1%, T-LBL 9.6% and B-LBL 0.6%)",
      notes              = "Time-fixed per subject. Multiplies clearance by 0.771, i.e. a 22.9% LOWER SAA clearance in T-ALL patients (Lin 2023 Table 2, RSE 8.89%, 95% CI 0.637-0.905; Discussion: 'modeled as fractional changes of ... 0.771'). 26 of 166 patients (15.7%; Lin 2023 Table 1). NOTE the indicator is T-ALL specifically and NOT T-lineage: the 16 T-LBL patients (9.6%) sit in the reference group, because Lin 2023 screened primary disease (ALL vs LBL) and disease subtype (B-cell vs T-cell) as separate covariates and retained only the T-ALL cell of that cross-classification (Lin 2023 Methods, Population PK modeling; Results, Covariate analysis).",
      source_name        = "Disease subtype (T-ALL)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as an intrinsic covariate on JZP458 PK but not retained in the final model (Lin 2023 Methods, Population PK modeling; Discussion). Analysis set 1.7-25 years."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained (Lin 2023 Discussion). 63 of 166 patients female (38.0%; Table 1)."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened but not retained (Lin 2023 Discussion). Analysis set 74.4-195.0 cm."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened but not retained; BSA was the retained body-size metric (Lin 2023 Discussion). Analysis set 9.33-131.0 kg."
    ),
    RACE_HISPANIC = list(
      description = "Hispanic or Latino ethnicity indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Ethnicity screened but not retained (Lin 2023 Discussion). 52 of 166 patients Hispanic or Latino (31.3%; Table 1)."
    ),
    DIS_ALL = list(
      description = "Primary disease indicator (ALL vs LBL)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Primary disease (ALL or LBL) screened but not retained; only the disease-SUBTYPE T-ALL effect survived backward elimination (Lin 2023 Discussion). 149 ALL / 17 LBL (Table 1)."
    ),
    ADA_POS = list(
      description = "Antidrug antibody positive indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained (Lin 2023 Discussion). Exactly 83 of 166 patients (50.0%) were ADA positive (Table 1); neutralizing-antibody positive in only 4 of 166 (2.4%)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 166L,
    n_studies      = 1L,
    n_observations = 2687L,
    age_range      = "1.7-25 years",
    age_median     = "10.0 years",
    weight_range   = "9.33-131.0 kg",
    weight_median  = "36.5 kg",
    bsa_range      = "0.44-2.53 m^2",
    bsa_median     = "1.17 m^2",
    sex_female_pct = 38.0,
    race_ethnicity = c(
      `White/Caucasian`               = 68.7,
      `Black/African American`        = 13.3,
      `Asian`                         = 4.2,
      `American Indian/Alaska Native` = 1.8,
      `Other`                         = 0.6,
      `Declined to state`             = 11.4
    ),
    disease_state  = "Newly diagnosed acute lymphoblastic leukemia (89.8%) or lymphoblastic lymphoma (10.2%) with a grade 3 or greater allergic reaction, or silent inactivation, to a long-acting E. coli-derived asparaginase. Disease subtype B-ALL 74.1%, T-ALL 15.7%, T-LBL 9.6%, B-LBL 0.6%. All 166 patients had received pegaspargase before study entry.",
    dose_range     = "Intramuscular JZP458 in three Monday/Wednesday/Friday cohorts: 25 mg/m^2 (n = 32), 37.5 mg/m^2 (n = 83) and 25/25/50 mg/m^2 (n = 51), six doses per two-week course. Injection volume at a single site capped at 2 mL.",
    regions        = "Children's Oncology Group sites (AALL1931; NCT04145531)",
    notes          = "Lin 2023 Table 1 (baseline demographics of the PopPK analysis set) and Results, Patient demographics. 2687 SAA observations, of which 2145 quantifiable and 542 below the limit of quantitation; BLOQ handled by the Beal M3 method. Assay calibration range 0.0349-0.2096 IU/mL with dilution linearity to 467.72-fold. NOTE Table 1 prints the age range as '(1.7-2.5)', which contradicts the same table's mean of 10.2 years and the Results text 'ranged from 1.7 to 25 years of age'; the range recorded here is the text value 1.7-25 years."
  )

  ini({
    # -----------------------------------------------------------------
    # Structural PK -- Lin 2023 Table 2 (Final JZP458 PopPK model
    # parameter estimates). Reference subject: BSA 1.2 m^2, non-Black/
    # African American, non-T-ALL.
    #
    # AMOUNT UNITS ARE ACTIVITY UNITS (IU), NOT MASS. The zero-order
    # absorption rate is reported in IU/h in both Lin 2023 Table 2 and the
    # upstream Lin 2021 Table 2, and the observation is serum asparaginase
    # ACTIVITY in IU/mL, so central/vc is IU/mL only if the dose is entered
    # in IU. Neither publication states the mg-to-IU specific activity of
    # JZP458, so a clinical mg/m^2 dose CANNOT be converted to model units
    # from any on-disk source -- see the vignette Errata section.
    # -----------------------------------------------------------------
    lcl <- log(146)                ; label("Clearance CL (mL/h)")                              # Lin 2023 Table 2 (CL 146 mL/h, SE 4.93, RSE 3.38%, 95% CI 136-156)
    lvc <- log(445)                ; label("Central volume of distribution V (mL)")            # Lin 2023 Table 2 (Central compartment V 445 mL, SE 34.9, RSE 7.85%, 95% CI 377-513)

    # Absorption parameters were NOT estimated here. Lin 2023 Methods,
    # Population PK modeling: "Due to sparse sampling in the AALL1931 study,
    # the absorption phase could not be fully characterized with data from
    # the AALL1931 study alone; as a result, three absorption-related
    # parameters (absorption rate constant [Ka], zero-order rate [R1], and
    # bioavailability [F]) were fixed to the values determined from a PopPK
    # model developed with intensive phase I data." Accordingly Lin 2023
    # Table 2 reports no SE / RSE / CI for these three rows.
    lka     <- fixed(log(0.0369))  ; label("First-order absorption rate constant Ka (1/h)")    # Lin 2023 Table 2 (Absorption rate constant 0.0369 1/h; no SE/RSE/CI reported); carried from the phase I model of Lin 2021 (doi:10.1002/cpdd.1002), which estimated 0.0348 1/h
    lr1     <- fixed(log(1810))    ; label("Zero-order input rate R1 into depot (IU/h)")       # Lin 2023 Table 2 (Zero-order rate 1810 IU/h; no SE/RSE/CI reported); Lin 2021 Table 2 estimated 4000 IU/h
    # NOTE Lin 2023 Table 2 labels this row "Relative bioavailability (%)" but
    # the value 0.359 is a FRACTION, not a percentage. Two checks in the
    # upstream Lin 2021 settle it: its Table 2 gives F = 0.365 while its
    # Results text calls the same quantity "bioavailability at 36.5%", and its
    # Table 2 footnote reconciles Vd = 3030 mL (i.v.) with Vd/F = 8.30 L (i.m.),
    # which requires 3030 / 0.365 = 8301 mL. Reading 0.359 as 0.359% would
    # inflate the apparent volume a hundredfold.
    lfdepot <- fixed(log(0.359))   ; label("Relative bioavailability F into depot (unitless fraction)") # Lin 2023 Table 2 (Relative bioavailability 0.359; no SE/RSE/CI reported); Lin 2021 Table 2 estimated 0.365

    # -----------------------------------------------------------------
    # Covariate effects -- Lin 2023 Table 2 and the printed final-model
    # equations in Results, Covariate analysis and final PopPK model
    # selection:
    #   CL [mL/h] = 146 * (BSA[m^2] / 1.2[m^2])^1.48
    #                   * 0.674 [if Black or African American]
    #                   * 0.771 [if T-ALL disease subtype]
    #   V  [mL]   = 445 * (BSA[m^2] / 1.2[m^2])^1.61
    # The two race / disease terms are multiplicative FRACTIONAL CHANGES,
    # not log-additive shifts, so they are written <estimate>^<indicator>
    # (following Kloos_2021_pegasparaginase.R).
    # -----------------------------------------------------------------
    e_bsa_cl   <- 1.48             ; label("Allometric exponent of BSA on CL (unitless)")      # Lin 2023 Table 2 (BSA on CL 1.48, SE 0.0686, RSE 4.63%, 95% CI 1.35-1.61)
    e_bsa_vc   <- 1.61             ; label("Allometric exponent of BSA on V (unitless)")       # Lin 2023 Table 2 (BSA on V 1.61, SE 0.151, RSE 9.4%, 95% CI 1.31-1.91)
    e_black_cl <- 0.674            ; label("Fractional change in CL for Black/African American race") # Lin 2023 Table 2 (Black/African American race on CL 0.674, SE 0.0715, RSE 10.6%, 95% CI 0.534-0.814)
    e_tall_cl  <- 0.771            ; label("Fractional change in CL for T-ALL disease subtype")       # Lin 2023 Table 2 (Disease subtype (T-ALL) on CL 0.771, SE 0.0686, RSE 8.89%, 95% CI 0.637-0.905)

    # -----------------------------------------------------------------
    # Interindividual variability -- Lin 2023 Table 2, "Exponential" error
    # type. The Estimate column carries the VARIANCE and the CV column the
    # corresponding coefficient of variation, which reconciles exactly via
    # CV = sqrt(exp(omega^2) - 1):
    #   sqrt(exp(0.125) - 1) = 0.365 -> 36.5% (printed)
    #   sqrt(exp(0.479) - 1) = 0.784 -> 78.4% (printed)
    # so 0.125 and 0.479 are variances, not SDs. The paper reports no
    # CL-V covariance, so the two etas are left uncorrelated.
    # -----------------------------------------------------------------
    etalcl ~ 0.125                 ; label("IIV variance on lcl")                             # Lin 2023 Table 2 (IIV on CL 0.125, CV 36.5%, SE 0.0173, RSE 13.8%, 95% CI 0.0911-0.159)
    etalvc ~ 0.479                 ; label("IIV variance on lvc")                             # Lin 2023 Table 2 (IIV on V 0.479, CV 78.4%, SE 0.0867, RSE 18.1%, 95% CI 0.309-0.649)

    # -----------------------------------------------------------------
    # Residual error -- combined proportional and additive, both reported
    # directly as standard deviations in Lin 2023 Table 2. The table's last
    # row, "Within-subject variability | Proportional | 1 | 100%", is the
    # NONMEM idiom of fixing $SIGMA to 1 so the residual magnitudes are
    # carried by the two $THETA rows above; it is not a third error term
    # and needs no separate parameter here.
    # -----------------------------------------------------------------
    propSd <- 0.46                 ; label("Proportional residual error SD (fraction)")       # Lin 2023 Table 2 (Proportional residual error (SD) 0.46, SE 0.00665, RSE 1.45%, 95% CI 0.447-0.473)
    addSd  <- 0.0166               ; label("Additive residual error SD (IU/mL)")              # Lin 2023 Table 2 (Additive residual error (SD) 0.0166, SE 0.000539, RSE 3.25%, 95% CI 0.0155-0.0177)
  })

  model({
    # 1. Covariate factors. Reference subject is BSA 1.2 m^2, non-Black/
    #    African American, non-T-ALL, for which every factor collapses to 1.
    bsa_cl  <- (BSA / 1.2)^e_bsa_cl
    bsa_vc  <- (BSA / 1.2)^e_bsa_vc
    race_cl <- e_black_cl^RACE_BLACK
    tall_cl <- e_tall_cl^DIS_TALL

    # 2. Individual parameters.
    cl     <- exp(lcl + etalcl) * bsa_cl * race_cl * tall_cl
    vc     <- exp(lvc + etalvc) * bsa_vc
    ka     <- exp(lka)
    r1     <- exp(lr1)
    fdepot <- exp(lfdepot)

    # 3. One-compartment disposition with mixed-order absorption. The dose
    #    is delivered into the depot as a zero-order input at rate r1 while
    #    first-order absorption ka drains the depot into the central
    #    compartment at the same time -- Lin 2023 Results, Base model:
    #    "simultaneous mixed-order absorption (zero-order infusion followed
    #    by a first-order absorption)"; Lin 2021 Results, Base model calls
    #    the same construct a "sequential mixed order absorption function".
    #    Because ka (0.0369 /h) is far smaller than kel = cl/vc (0.328 /h at
    #    the reference BSA), disposition is absorption rate limited and the
    #    terminal half-life is ln(2)/ka = 18.8 h, matching the post-hoc
    #    t1/2 of 19.1 h reported in Lin 2023 Results.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - (cl / vc) * central

    # 4. Bioavailability and the zero-order depot input, giving the NONMEM R1
    #    semantics the source model used: the input duration is
    #    fdepot * amt / r1.
    #
    #    IMPORTANT -- dose records MUST carry `rate = -1`. rxode2 only consults
    #    a modelled `rate()` when the dose record requests it; with the default
    #    `rate = 0` the modelled rate is SILENTLY IGNORED and the dose is given
    #    as an instantaneous bolus into the depot. That failure is quiet and
    #    changes the profile materially -- verified on this model at BSA 1.2
    #    with a 20000 IU dose: `rate = -1` puts the depot peak at 4.0 h
    #    (= fdepot * amt / r1 = 3.97 h) and Tmax at 9.7 h, whereas `rate = 0`
    #    peaks the depot at t = 0 and pulls Tmax down to 7.5 h
    #    (= ln(kel/ka) / (kel - ka), the bolus value). Cmax differs by under 1%
    #    between the two, so Cmax alone will NOT reveal the mistake.
    f(depot)    <- fdepot
    rate(depot) <- r1

    # 5. Observation: serum asparaginase activity (SAA) in IU/mL.
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
