Verrest_2023_miltefosine <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order oral absorption",
    "for miltefosine in 265 Eastern African children and adults with",
    "visceral leishmaniasis (Verrest 2023), enrolled in the phase III",
    "randomized controlled trial NCT03129646 across six sites in Kenya,",
    "Sudan, Ethiopia and Uganda and treated with allometric oral",
    "miltefosine for 14 or 28 days combined with intramuscular paromomycin",
    "for 14 days. CL/F, Vc/F and Vp/F are allometrically scaled on",
    "fat-free mass (exponents fixed at 0.75 and 1.0; reference 18 kg);",
    "Q/F is not scaled. Relative bioavailability is structurally fixed at",
    "1 and modified by two non-linearities: it is 65% lower during the",
    "first week of treatment, with large between-subject variability",
    "(74.8% CV), and it falls as a power function of the cumulative",
    "miltefosine dose once that dose exceeds 60 mg/kg, which reproduces",
    "the observed stagnation of miltefosine accumulation in the third",
    "week of treatment. Between-subject variability is log-normal on CL/F",
    "(16.3% CV) and on the first-week bioavailability reduction; residual",
    "error is proportional (31.5% CV). The companion paromomycin model",
    "from the same trial is Verrest_2023_paromomycin."
  )
  reference <- paste(
    "Verrest L, Roseboom IC, Wasunna M, Mbui J, Njenga S, Musa AM, Olobo J,",
    "Mohammed R, Ritmeijer K, Chu WY, Huitema ADR, Solomos A, Alves F,",
    "Dorlo TPC. Population pharmacokinetics of a combination of miltefosine",
    "and paromomycin in Eastern African children and adults with visceral",
    "leishmaniasis. J Antimicrob Chemother. 2023;78(11):2702-2714.",
    "doi:10.1093/jac/dkad286. ClinicalTrials.gov NCT03129646.",
    sep = " "
  )
  vignette <- "Verrest_2023_miltefosine_paromomycin"
  units    <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Amounts are mg of miltefosine, matching the trial
  # dataset's AMT column ("MF: mg", supplementary NONMEM control stream
  # $INPUT dataset description).
  compartmentData <- list(
    depot       = list(analyte = "miltefosine", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "miltefosine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "miltefosine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric power scaling on CL/F (exponent fixed 0.75) and on",
        "Vc/F and Vp/F (exponent fixed 1.00), normalized to a reference",
        "of 18 kg. Q/F is NOT scaled on FFM. The 18 kg reference is not",
        "printed in Verrest 2023 Table 4, whose equations write only the",
        "symbol FFM_med; it is read from the supplementary NONMEM",
        "control stream ('ALLOCL = (FFM/18)**0.75', 'ALLOV =",
        "(FFM/18)**1'). Computed per subject from body weight, height",
        "and sex: the control stream's dataset description records 'FFM:",
        "Fat free mass (kg) according Maturation Model Al-Sallami 2015',",
        "i.e. the Janmahasatian adult formula with the Al-Sallami",
        "paediatric age-correction multiplier, which is the appropriate",
        "choice for a cohort that is 59% children. The same FFM",
        "descriptor also determined the allometric miltefosine dose",
        "given to children under 30 kg (Verrest 2023 Figure S1)."
      ),
      source_name        = "FFM"
    ),
    DOSE_MF_CUM_MGKG = list(
      description        = "Cumulative administered miltefosine dose per kg body weight",
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TIME-VARYING: the running total of miltefosine administered to",
        "the subject up to the current record, divided by body weight",
        "(control stream dataset description: 'AMTT: total MF dose',",
        "'DDOS: AMTT/body weight'). Drives the second bioavailability",
        "non-linearity, which the paper attributes to slow, saturable",
        "transport of miltefosine across the gastrointestinal membrane",
        "(Discussion). The effect is switched on only once the",
        "cumulative dose reaches 60 mg/kg and is then a power function",
        "normalized to 70 mg/kg. Neither constant is printed in Table 4:",
        "the 60 mg/kg threshold appears only in Table 4 footnote c",
        "('Applied after a cumulative dose of 60 mg/kg is reached'), and",
        "the 70 mg/kg normalizing value appears only in the",
        "supplementary NONMEM control stream ('IF(DDOS.GE.60) COVF2=",
        "((DDOS/70)**THETA(8))'). Both are confirmed by the paper's own",
        "worked example -- a 35 kg patient on 100 mg/day reaches 77.1",
        "mg/kg after 27 days of dosing, and (77.1/70)^-2.40 = 0.792,",
        "i.e. the 21% lower bioavailability on Day 28 stated in Results."
      ),
      source_name        = "DDOS"
    )
  )

  covariatesDataExcluded <- list(
    MICR = list(
      description        = "Leishmania parasite load by microscopy score at treatment start (0-6)",
      units              = "(ordinal score)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Screened as a marker of disease severity on the first-week",
        "bioavailability reduction and not retained: 'No other",
        "covariates could be identified to explain the non-linear",
        "effects on bioavailability' (Verrest 2023 Results). Scored in",
        "spleen, bone marrow or lymph node aspirate."
      ),
      source_name        = "MICR"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened as a disease-severity / malnutrition marker on the",
        "first-week bioavailability reduction and not retained (Verrest",
        "2023 Methods 'Covariate analysis'; Results). Cohort value 28.3",
        "g/L mean (range 1.5-54.6)."
      ),
      source_name        = "ALB"
    ),
    NEUT = list(
      description        = "Absolute neutrophil count",
      units              = "cells/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened as a disease-severity marker on the miltefosine",
        "first-week bioavailability reduction and not retained. It IS",
        "retained in the companion paromomycin model from the same",
        "trial, where it drives the decrease in clearance over time",
        "(Verrest_2023_paromomycin)."
      ),
      source_name        = "NEUTR"
    ),
    WBC = list(
      description        = "Total white blood cell count",
      units              = "cells/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened as a disease-severity marker on the first-week",
        "bioavailability reduction and not retained (Verrest 2023",
        "Methods 'Covariate analysis')."
      ),
      source_name        = "WBC"
    ),
    HAZ = list(
      description        = "Height-for-age Z-score",
      units              = "(Z-score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened as a nutritional-status marker on the first-week",
        "bioavailability reduction for patients up to 19 years old and",
        "not retained: 'variables associated with malnutrition, such as",
        "height-for-age or BMI-for-age, were not identified as",
        "explanatory covariates on the initially decreased",
        "bioavailability' (Verrest 2023 Discussion). Derived with the R",
        "package 'zscorer' against the WHO Child Growth Standards."
      ),
      source_name        = "HAZ"
    ),
    BAZ = list(
      description        = "BMI-for-age Z-score",
      units              = "(Z-score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened as a nutritional-status marker on the first-week",
        "bioavailability reduction for patients up to 19 years old and",
        "not retained (Verrest 2023 Discussion). Derived with the R",
        "package 'zscorer' against the WHO Child Growth Standards."
      ),
      source_name        = "BAZ"
    ),
    BMI = list(
      description        = "Body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened as a nutritional-status marker on the first-week",
        "bioavailability reduction for patients from 19 years old and",
        "not retained (Verrest 2023 Methods 'Covariate analysis').",
        "Cohort value 17.7 kg/m^2 mean (range 14.6-21.3)."
      ),
      source_name        = "BMI"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 265L,
    n_studies      = 1L,
    age_range      = "4-45 years",
    age_mean       = "13.6 years (sparse sampling group); 13.8 years (intensive sampling group)",
    weight_range   = "11.0-71.0 kg",
    weight_mean    = "32.9 kg (sparse sampling group); 33.6 kg (intensive sampling group)",
    sex_female_pct = 19.2,
    race_ethnicity = "Eastern African (Kenyan, Sudanese, Ethiopian and Ugandan cohorts; no further breakdown reported)",
    disease_state  = paste(
      "Eastern African children and adults with symptomatic,",
      "parasitologically confirmed visceral leishmaniasis; 59% were",
      "paediatric (<=12 years). Patients with relapse, severe",
      "malnutrition, severe VL, HIV co-infection or concomitant severe",
      "infection were excluded. Patients were severely ill at treatment",
      "start, frequently with anaemia, leucopenia and malnourishment,",
      "which the authors propose as the explanation for the reduced",
      "first-week bioavailability."
    ),
    dose_range     = paste(
      "Oral miltefosine twice daily for 14 days (PM+MF14D arm) or 28",
      "days (PM+MF28D arm). Children under 30 kg received an allometric",
      "dose based on sex, weight and height, i.e. on fat-free mass",
      "(Figure S1); patients 30 to under 45 kg received 100 mg/day and",
      "patients 45 kg and over received 150 mg/day. A re-dose was given",
      "if the patient vomited within 30 min of administration. All",
      "patients also received intramuscular paromomycin sulphate 20",
      "mg/kg/day for 14 days."
    ),
    regions        = "Eastern Africa (Kenya: Kacheliba; Uganda: Amudat; Sudan: Doka, Um El Kher; Ethiopia: Gondar, Abdurafi)",
    co_medication  = paste(
      "Intramuscular paromomycin given simultaneously in every patient",
      "for 14 days. The exposures achieved for both drugs matched",
      "previous monotherapy studies, which the authors read as evidence",
      "against a drug-drug interaction; no interaction term is present",
      "in this model. The companion paromomycin model is",
      "Verrest_2023_paromomycin."
    ),
    samples        = paste(
      "910 miltefosine plasma concentrations after exclusions (927 per",
      "the Results text; 309 from the 26-patient intensive-sampling",
      "cohort plus 601 sparse per Table 2). Sparse sampling was on Day",
      "14 (PM+MF14D only), Day 28 and Day 56; the intensive cohort added",
      "0, 1, 2, 4 and 24 h or 0, 1, 2, 8 and 24 h on Day 1 and Day 14.",
      "Seventeen observations were excluded: three physiologically",
      "implausible and 14 (2.8%) below the 4 ng/mL limit of",
      "quantification. Where a patient vomited, the pre-vomiting dose",
      "was excluded and the re-dose retained (n = 20)."
    ),
    notes          = paste(
      "NONMEM 7.5, ADVAN13, FOCE-I with interaction; parameter precision",
      "by sampling importance resampling. The starting point was the",
      "earlier Eastern African miltefosine model of Dorlo 2017 (Verrest",
      "2023 reference 15), which is in this library as",
      "Dorlo_2017_miltefosine and shares the structure of a two-",
      "compartment model with an FFM-allometric disposition and a",
      "reduced first-week bioavailability. Efficacy target: the in vitro",
      "susceptibility value EC90 of 10.6 ug/mL, against which the paper",
      "reports Time > EC90. Demographics are reported in Table 1 as mean",
      "(range) at baseline; the sex percentage recorded here is derived",
      "from the 193 of 239 sparsely sampled patients reported as male."
    )
  )

  ini({
    # ============================================================
    # Structural PK parameters -- Verrest 2023 Table 4, 'Population
    # parameters' block, expressed per DAY to match the paper's own
    # exposure units (AUC in ug*day/mL, Time > EC90 in days).
    #
    # UNITS NOTE, load-bearing: the supplementary NONMEM control stream
    # ('NONMEM control stream miltefosine.docx') runs on the dataset's
    # TIME column, which is in HOURS ('TIME: time after first dose
    # (hr)'; the first-week bioavailability switch is written
    # 'IF(TIME.LE.168)', i.e. 7 days expressed in hours). Every rate
    # constant in that stream is therefore per hour. Table 4 converted
    # the clearances to per-day -- its CL of 1.85 L/day and Q of 0.17
    # L/day are exactly 24x the stream's $THETA values of 0.077 and
    # 0.007 L/h -- but printed ka as the raw hourly value 0.037 while
    # labelling the row 'ka (day-1)'. That label is a units error: an
    # absorption rate constant of 0.037/day is an absorption half-life
    # of 18.7 DAYS, which cannot produce the miltefosine measured 1-4 h
    # after the first dose in the intensive-sampling cohort and would
    # put AUC_D0-7 near zero against the observed median of ~20
    # ug*day/mL. ka is encoded here as 0.037/h = 0.888/day. See the
    # vignette's Assumptions and deviations section.
    # ============================================================
    lcl <- log(1.85)
    label("Apparent clearance CL/F at FFM = 18 kg (L/day)")                  # Table 4 'CL (L/day)': 1.85 (95% CI 1.75-1.94); control stream $THETA(1) 0.077 L/h x 24 = 1.848
    lvc <- log(13.6)
    label("Apparent central volume Vc/F at FFM = 18 kg (L)")                 # Table 4 'Vc (L)': 13.6 (95% CI 12.8-14.4); control stream $THETA(2)
    lq  <- log(0.17)
    label("Apparent inter-compartmental clearance Q/F (L/day)")              # Table 4 'Q (L/day)': 0.17 (95% CI 0.13-0.21); control stream $THETA(4) 0.007 L/h x 24 = 0.168
    lvp <- log(2.22)
    label("Apparent peripheral volume Vp/F at FFM = 18 kg (L)")              # Table 4 'Vp (L)': 2.22 (95% CI 1.96-2.59); control stream $THETA(5)
    lka <- log(0.888)
    label("First-order oral absorption rate constant ka (1/day)")            # Table 4 'ka' estimate 0.037 with the row labelled day-1; the control stream runs on an hours TIME axis so 0.037 is per hour -> 0.037 * 24 = 0.888 /day. See the units note above.

    # ============================================================
    # Bioavailability. F is structurally fixed at 1 (unidentifiable
    # from oral-only data) and is then multiplied by two independent
    # non-linear factors, exactly as the control stream composes them:
    # F1 = TVF1 * COVF * COVF2.
    # ============================================================
    lfdepot <- fixed(log(1))
    label("Reference relative bioavailability F1 (unitless)")         # Table 4 'F1': 1 (fixed); control stream $THETA(6) '(1) FIX ; F fix'

    # Factor 1: the first-week reduction. Table 4 reports COV_F,W1 as a
    # 'fractional change' of -0.65, i.e. a 65% DECREASE, which the
    # Results state directly ('bioavailability was 65% (95% CI 57%-73%)
    # lower during the first week of treatment'). The Table 4 equation
    # prints the factor as (1 - COV_F,W1), which with the tabulated
    # sign of -0.65 would read as a 65% INCREASE; the control stream
    # settles it -- THETA(7) is bounded (0, 0.65, 1) and the factor is
    # '(1 - THETA(7))*EXP(ETA(4))', giving a multiplier of 0.35.
    # Re-parameterised here on the log scale so that the individual
    # multiplier stays positive across the log-normal BSV draw, which
    # is also how the sibling Dorlo_2017_miltefosine encodes it.
    lfred_mult <- log(1 - 0.65)
    label("Log of the F multiplier during the first week of treatment (unitless)") # Table 4 'COV F,W1 (fractional change)': -0.65 (95% CI -0.73 to -0.57) -> multiplier 1 - 0.65 = 0.35; control stream 'IF(TIME.LE.168) COVF = (1 - THETA(7))*EXP(ETA(4))' with $THETA(7) = 0.65

    # Factor 2: the cumulative-dose reduction, a power function applied
    # only once the cumulative dose reaches 60 mg/kg.
    e_dose_mf_cum_mgkg_fdepot <- -2.40
    label("Power exponent of cumulative miltefosine dose (mg/kg) on relative bioavailability") # Table 4 'COV F,CD (exponent of power relationship)': -2.40 (95% CI -3.79 to -1.21); footnote c 'Applied after a cumulative dose of 60 mg/kg is reached'; control stream $THETA(8) with 'IF(DDOS.GE.60) COVF2 = ((DDOS/70)**THETA(8))'

    # ============================================================
    # Covariate effects -- allometric exponents, not estimated but
    # fixed at the standard theory-based values (Methods, 'Covariate
    # analysis'). Note that the control stream applies allometry to
    # CL, V2 and V3 but NOT to Q, so no shared CL/Q exponent exists
    # here (unlike the companion paromomycin model).
    # ============================================================
    e_ffm_cl    <- fixed(0.75)
    label("Allometric exponent on fat-free mass for CL/F (unitless)")        # Table 4 equation exponent 0.75; control stream 'ALLOCL = (FFM/18)**0.75'
    e_ffm_vc_vp <- fixed(1.00)
    label("Allometric exponent on fat-free mass shared by Vc/F and Vp/F (unitless)") # Table 4 equation exponent 1.00; control stream 'ALLOV = (FFM/18)**1'

    # ============================================================
    # Inter-individual variability -- Verrest 2023 Table 4
    # 'Between-subject variability'. BSV is carried by CL/F and by the
    # first-week bioavailability reduction only; the control stream
    # fixes the V2 and ka etas to zero ($OMEGA '0 FIX'). As in the
    # paromomycin model, the paper's CV% is sqrt(omega^2): 16.3%
    # squares to 0.0266 and 74.8% squares to 0.5595, reproducing the
    # control stream's $OMEGA values of 0.0265 and 0.56.
    # ============================================================
    etalcl        ~ 0.163^2
    # Table 4 BSV 'CL (CV%)': 16.3 (95% CI 14.3-18.5), shrinkage 14% -> omega^2 = 0.163^2 = 0.02657; control stream $OMEGA '0.0265 ; IIV CL'
    etalfred_mult ~ 0.748^2
    # Table 4 BSV 'COV F,W1 (CV%)': 74.8 (95% CI 62.0-90.3), shrinkage 48% -> omega^2 = 0.748^2 = 0.5595; control stream $OMEGA '0.56 ; IIV F'. Applied on the log-scale F multiplier so F stays positive on every draw.

    # ============================================================
    # Residual error -- Verrest 2023 Table 4 'Residual variability'.
    # Control stream: Y = IPRED * (1 + EPS(1)), a pure proportional
    # error; $SIGMA 0.0992 has sqrt 0.315, matching Table 4.
    # ============================================================
    propSd <- 0.315
    label("Proportional residual error on Cc (fraction)")                    # Table 4 'Proportional error (CV%)': 31.5 (95% CI 29.7-33.6), shrinkage 14%
  })

  model({
    # ------------------------------------------------------------
    # Individual PK parameters. Allometry on fat-free mass normalized
    # to 18 kg (control stream 'ALLOCL'/'ALLOV'). Q/F is deliberately
    # NOT scaled -- the control stream assigns it as a bare 'Q = TVQ'.
    # ------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (FFM / 18)^e_ffm_cl
    vc <- exp(lvc) * (FFM / 18)^e_ffm_vc_vp
    q  <- exp(lq)
    vp <- exp(lvp) * (FFM / 18)^e_ffm_vc_vp
    ka <- exp(lka)

    # ------------------------------------------------------------
    # Bioavailability: F1 = 1 (fixed) * covf * covf2.
    #
    # covf -- the first-week reduction. The control stream's predicate
    # is TIME <= 168 h, i.e. t <= 7 days on this model's day axis, and
    # it INCLUDES t = 0 (unlike Dorlo 2017, which used t > 0).
    #
    # covf2 -- the cumulative-dose reduction, held at exactly 1 below
    # 60 mg/kg and switched to (CD/70)^-2.40 at and above it. NOTE that
    # the published parameterisation is DISCONTINUOUS at the threshold,
    # because the switch point (60) and the normalizing value (70) are
    # different numbers: (60/70)^-2.40 = 1.45, so F jumps UP by 45% at
    # 60 mg/kg, decays back through 1 at 70 mg/kg, and only then falls
    # below 1 (0.79 at 77 mg/kg, the value behind the paper's own
    # 'bioavailability was 21% lower on Day 28' example). This is the
    # model exactly as published and fitted -- both constants are
    # confirmed against the paper's worked example -- and it is
    # reproduced faithfully here rather than smoothed. See the
    # vignette's Assumptions and deviations section.
    # ------------------------------------------------------------
    fred_active <- (t <= 7)
    covf  <- exp(lfred_mult + etalfred_mult) * fred_active + 1 * (1 - fred_active)

    cd_active <- (DOSE_MF_CUM_MGKG >= 60)
    covf2 <- (DOSE_MF_CUM_MGKG / 70)^e_dose_mf_cum_mgkg_fdepot * cd_active +
      1 * (1 - cd_active)

    # ------------------------------------------------------------
    # Two-compartment oral PK with first-order absorption and linear
    # elimination from the central compartment (control stream $DES).
    # ------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka  * depot
    d/dt(central)     <-  ka  * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot) * covf * covf2

    # ------------------------------------------------------------
    # Observation. Amounts are mg and vc is in L, so central / vc is
    # mg/L = ug/mL, the unit in which the paper reports the EC90
    # target of 10.6 ug/mL and exposures in ug*day/mL.
    # ------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
