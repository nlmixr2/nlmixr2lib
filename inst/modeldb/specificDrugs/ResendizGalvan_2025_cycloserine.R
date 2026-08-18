ResendizGalvan_2025_cycloserine <- function() {
  description <- paste(
    "One-compartment population PK model for oral cycloserine in Indian",
    "adolescents and adults treated for multidrug-resistant tuberculosis",
    "(Resendiz-Galvan 2025). Savic transit-compartment absorption (N = 4",
    "fixed, MTT = 0.610 h) feeds first-order absorption into a",
    "one-compartment disposition model. Total clearance is split into a",
    "renal arm scaled linearly by a size- and sex-neutral renal-function",
    "ratio (Cockcroft-Gault creatinine clearance recomputed for a 56 kg",
    "male, normalised to the cohort median of 122 mL/min) and a non-renal",
    "arm; the sum is then scaled allometrically by fat-free mass",
    "(reference 38.6 kg). Between-occasion variability on all absorption",
    "parameters is inflated 2.26-fold on occasions following unobserved,",
    "self-reported doses taken at home.",
    sep = " "
  )
  reference <- paste(
    "Resendiz-Galvan JE, Arora PR, Lokhande RV, Udwadia ZF, Rodrigues C,",
    "Gupta A, Tornheim JA, Denti P, Ashavaid TF. Evaluation of cycloserine",
    "dose regimens in an Indian cohort with multidrug-resistant",
    "tuberculosis: a population pharmacokinetic analysis. Antimicrob Agents",
    "Chemother. 2025 Oct;69(10):e00101-25. doi:10.1128/aac.00101-25.",
    "PMCID: PMC12486832. (Version of record posted 1 October 2025 correcting",
    "the author-contributions statement of the 2 September 2025 original; no",
    "model parameter was revised by that correction.)",
    sep = " "
  )
  vignette <- "ResendizGalvan_2025_cycloserine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot   = list(analyte = "cycloserine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cycloserine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass (Janmahasatian formula)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Resendiz-Galvan 2025 Methods ('Model development and PK analysis'):",
        "FFM estimated from sex, height, and body weight per the paper's",
        "reference 25 (Janmahasatian et al. Clin Pharmacokinet",
        "2005;44:1051-1065). Allometric scaling on FFM described body size",
        "better than total body weight (dOFV 20.9 vs 12.8). Cohort median",
        "FFM = 38.7 kg (Table 1, IQR 32.3-47.1); the reference value used",
        "for the typical individual is 38.6 kg (Table 2 footnote a, and",
        "Results: 'For a typical male with an FFM of 38.6 kg'), which is",
        "the value the model normalises by. Applied to the summed renal +",
        "non-renal clearance (exponent 0.75) and to central volume",
        "(exponent 1); see e_ffm_cl / e_ffm_vc in ini() for the",
        "exponent-provenance note.",
        sep = " "
      ),
      source_name        = "FFM"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Cohort median 27 years (Table 1, IQR 21-35); adolescents and",
        "adults aged >= 15 years were eligible (Methods, 'Study",
        "population'). Age enters the model only through the",
        "Cockcroft-Gault numerator (140 - AGE) of the size- and",
        "sex-neutral creatinine clearance CLcr,56M used to build the renal",
        "function ratio RF; it carries no separate effect on any PK",
        "parameter.",
        sep = " "
      ),
      source_name        = "AGE"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Cohort median 0.7 (Table 1, IQR 0.6-0.9). Table 1 labels the units",
        "'mg/L', but the printed Cockcroft-Gault equation in Methods reads",
        "'72 . SCr (mg/dL)' and the arithmetic confirms mg/dL: the median",
        "individual (AGE 27, SCr 0.7 mg/dL) gives CLcr,56M =",
        "(140 - 27) * 56 / (72 * 0.7) = 125.6 mL/min, consistent with the",
        "stated median CLcr,56M of 122 mL/min (Table 2 footnote b), whereas",
        "mg/L would give ~1,256 mL/min. The Table 1 unit label is therefore",
        "a typographical error; mg/dL is used here. Enters only through",
        "CLcr,56M / RF (see AGE).",
        sep = " "
      ),
      source_name        = "SCr"
    ),
    OCC = list(
      description        = "Dose-observation occasion: 1 = occasion following unobserved self-reported doses taken at home, 2 = occasion following a dose directly observed at the PK visit",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Resendiz-Galvan 2025 Methods: 'For between-occasion variability",
        "estimation, non-observed doses before the PK visits and pre-dose",
        "concentrations were considered as an independent occasion,",
        "separate from the observed doses during the visits and the",
        "subsequent concentrations.' The two occasion types are decomposed",
        "inside model() into binary indicators oc1 (unobserved / home) and",
        "oc2 (observed / clinic) that (a) select the per-occasion BOV etas",
        "on ka, MTT, and bioavailability and (b) apply the 2.26-fold",
        "inflation of those BOV standard deviations on the unobserved-dose",
        "occasion (Table 2, 'Variability factor for unobserved doses').",
        "For simulating a prescribed, fully-observed regimen set OCC = 2",
        "throughout; that is the setting under which the paper's Table 3",
        "simulated exposures are reproduced.",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 180L,
    n_studies       = 1L,
    n_observations  = 1281L,
    age_range       = "Adolescents and adults >= 15 years; Table 1 median 27 years (IQR 21-35)",
    weight_range    = "Table 1: median 55.5 kg (IQR 46.0-65.9)",
    height_range    = "Table 1: median 1.59 m (IQR 1.53-1.69)",
    ffm_range       = "Table 1: median 38.7 kg (IQR 32.3-47.1); model reference 38.6 kg",
    creat_range     = "Table 1: median 0.7 mg/dL (IQR 0.6-0.9); Table 1 mislabels the unit as mg/L",
    crcl_range      = "Table 1: median 108 mL/min (IQR 91.7-133.0), Cockcroft-Gault with actual weight and sex",
    sex_female_pct  = 65.0,
    n_hiv_positive  = 4L,
    n_current_smokers = 8L,
    renal_function  = paste(
      "Largely normal; the size- and sex-neutral CLcr,56M used to build the",
      "renal-function ratio had a cohort median of 122 mL/min",
      "(Table 2 footnote b)."
    ),
    disease_state   = paste(
      "Treatment-naive multidrug-resistant pulmonary tuberculosis on",
      "individualised, susceptibility-guided longer MDR-TB regimens.",
      "Pre-treatment sputum-isolate cycloserine MICs (n = 171) had a median",
      "of 16 mg/L (range 1-64); 76.5% of isolates had MIC >= 16 mg/L."
    ),
    dose_range      = paste(
      "Oral cycloserine (Lupin 250 mg and 500 mg capsules) started at 250 mg",
      "once daily and escalated as tolerated to a total daily dose of 500 or",
      "750 mg following Indian national weight-band guidance. Regimens",
      "received were 250 or 500 mg QD, 250 mg BID or TID, 250/500 mg",
      "a.m./p.m., and 500 mg BID; at the first PK visit 85% of participants",
      "were on 250 mg BID and 9% on 250/500 mg a.m./p.m."
    ),
    regions         = "India (single tertiary hospital, Mumbai)",
    co_medication   = paste(
      "Individualised MDR-TB regimens: linezolid (91% of participants),",
      "moxifloxacin (90%), clofazimine (79%), pyrazinamide (42%),",
      "ethambutol (34%), para-aminosalicylic acid (25%), bedaquiline (28%),",
      "ethionamide (22%), kanamycin (20%). No drug-drug interaction with",
      "cycloserine was identified."
    ),
    notes           = paste(
      "Prospective observational cohort enrolled October 2017 - February",
      "2022 (MDR-TB MUKT / Indo-South Africa study teams). 1,281 cycloserine",
      "observations: 312 intensive (pre-dose, 1, 2, 4, 6, 8 h) and 969",
      "sparse (pre-dose and 2 h post-dose) from visits at 1, 2, 6, and 12",
      "months of treatment. Assay linear over 0.782-50 mg/L; three",
      "below-limit-of-detection values (0.2%) were imputed at LOD/2 =",
      "0.195 mg/L under an adaptation of the M6 method. More than 80% of",
      "participants took cycloserine with food at each PK visit, with no",
      "detectable food effect."
    )
  )

  ini({
    # ---- Structural disposition parameters ------------------------------
    # All typical values are FINAL estimates from Resendiz-Galvan 2025
    # Table 2, reported for a typical male with FFM 38.6 kg and the cohort
    # median size-/sex-neutral creatinine clearance CLcr,56M of 122 mL/min
    # (Table 2 footnotes a and b; Results, 'Pharmacokinetic model').
    lcl_renal  <- log(0.589); label("Renal clearance CLr at the reference individual (L/h)")      # Table 2: renal clearance 0.589 L/h (95% CI 0.449-0.745)
    lcl_nonren <- log(0.901); label("Non-renal clearance CLnr at the reference individual (L/h)") # Table 2: non-renal clearance 0.901 L/h (95% CI 0.767-1.04)
    lvc        <- log(37.0);  label("Central volume of distribution at the reference individual (L)") # Table 2: central volume 37.0 L (95% CI 35.1-38.9)

    # ---- Absorption ------------------------------------------------------
    lka  <- log(2.15);  label("Absorption rate constant ka out of the depot (1/h)")   # Table 2: absorption rate constant 2.15 1/h (95% CI 1.51-3.55)
    lmtt <- log(0.610); label("Mean transit time of the Savic transit chain (h)")     # Table 2: mean transit time 0.610 h (95% CI 0.362-0.833)
    lnn  <- fixed(log(4)); label("Number of transit compartments N (unitless)")       # Table 2: NN = 4 fixed (footnote c: fixed to 4 to improve parameter-estimate stability; sensitivity analysis confirmed no significant effect on fit)

    # Bioavailability was fixed to 1; it is not identifiable from oral-only
    # data, so all disposition parameters above are apparent (CL/F, V/F).
    lfdepot <- fixed(log(1)); label("Bioavailability (unitless)")                     # Table 2: bioavailability = 1 fixed

    # ---- Allometric exponents -------------------------------------------
    # The paper states that allometric scaling on FFM was applied to the
    # disposition parameters and cites Holford et al. (reference 27) for the
    # renal / non-renal separation approach, but it does not print the
    # exponents and Table 2 contains no exponent row (i.e. they were not
    # estimated). They are therefore the standard theory-based allometric
    # values, fixed here as such and flagged in the vignette Errata.
    e_ffm_cl <- fixed(0.75); label("FFM allometric exponent on total clearance (unitless)")   # Methods: 'Allometric scaling ... on disposition parameters'; standard theory-based value, not printed in the source
    e_ffm_vc <- fixed(1.0);  label("FFM allometric exponent on central volume (unitless)")    # Methods: 'Allometric scaling ... on disposition parameters'; standard theory-based value, not printed in the source

    # ---- Effect of an unobserved dose on the BOV magnitude ---------------
    # Multiplies the BOV standard deviation (not the variance) of every
    # absorption parameter on occasions that follow a self-reported dose
    # taken at home. Same construction as the covariate-on-variability
    # thetas in Wojciechowski_2023_ritlecitinib_*.R.
    e_occunobs_iov <- 2.26; label("Fold change in the BOV SD of the absorption parameters on unobserved-dose occasions (unitless)") # Table 2: variability factor for unobserved doses = 2.26 (95% CI 1.97-2.54); footnote d: 'applies to the BOV of absorption rate constant, mean transit time, and bioavailability for self-reported doses taken at home on days prior to blood sampling'

    # ---- Between-subject variability ------------------------------------
    # Table 2 footnote f states the reported %CV is computed as
    # %CV = sqrt(omega^2) * 100 -- i.e. the tabulated number is 100 * omega
    # on the log scale, NOT 100 * sqrt(exp(omega^2) - 1). The variances
    # below are therefore (CV/100)^2 and not log((CV/100)^2 + 1).
    # BSV in CL was estimated as a single random effect shared by both
    # elimination routes (Results: 'Between-subject variability in CL was
    # included on both elimination routes as a single random effect'), so it
    # multiplies the summed renal + non-renal clearance in model().
    etalcl ~ 0.116281  # Table 2: BSV on clearance 34.1% CV (95% CI 30.2-38.7) -> omega^2 = 0.341^2 = 0.116281
    # Central volume carries no BSV (Table 2 reports 'NA' in the variability
    # column for the central-volume row).

    # ---- Between-occasion variability -----------------------------------
    # BOV was placed on all absorption parameters (Results). Two occasion
    # types are declared -- unobserved (home) and observed (clinic) -- so
    # each gets an independent draw; the two variances are equal, encoded
    # with fix() on the second in the NONMEM '$OMEGA BLOCK(1) SAME' idiom
    # used by Jonsson_2011_ethambutol.R and Smythe_2013_gatifloxacin.R.
    # The 2.26-fold inflation on the unobserved occasion is applied to the
    # SD in model() rather than being folded into these variances.
    etaiov_ka_1  ~ 0.894916       # Table 2: BOV on ka 94.6% CV (95% CI 79.4-121) -> omega^2 = 0.946^2 = 0.894916 (unobserved-dose occasion)
    etaiov_ka_2  ~ fix(0.894916)  # same BOV variance on the observed-dose occasion
    etaiov_mtt_1 ~ 0.697225       # Table 2: BOV on mean transit time 83.5% CV (95% CI 67.3-111) -> omega^2 = 0.835^2 = 0.697225 (unobserved-dose occasion)
    etaiov_mtt_2 ~ fix(0.697225)  # same BOV variance on the observed-dose occasion
    etaiov_fdepot_1 ~ 0.042849       # Table 2: BOV on bioavailability 20.7% CV (95% CI 18.9-22.7) -> omega^2 = 0.207^2 = 0.042849 (unobserved-dose occasion)
    etaiov_fdepot_2 ~ fix(0.042849)  # same BOV variance on the observed-dose occasion

    # ---- Residual unexplained variability --------------------------------
    # Combined proportional + additive on the linear concentration scale.
    # The source additionally inflated the additive term by LOD/2 for the
    # three below-limit-of-detection records it imputed (M6 adaptation,
    # Methods); that record-specific inflation affects 0.2% of observations
    # and is not carried in the packaged model (see vignette Errata).
    propSd <- 0.0458; label("Proportional residual error (fraction)")  # Table 2: proportional error 4.58% (95% CI 3.23-5.66)
    addSd  <- 0.856;  label("Additive residual error (mg/L)")          # Table 2: additive error 0.856 mg/L (95% CI 0.716-1.03)
  })

  model({
    # 1. Occasion indicators (binary decomposition of the OCC column).
    oc1 <- (OCC == 1)  # occasion following an unobserved, self-reported home dose
    oc2 <- (OCC == 2)  # occasion following a dose observed at the PK visit

    # 2. Per-occasion BOV etas on the absorption parameters. The
    #    unobserved-dose occasion carries the 2.26-fold inflated SD; the
    #    observed-dose occasion carries the Table 2 SD unmodified.
    iov_ka  <- oc1 * e_occunobs_iov * etaiov_ka_1  + oc2 * etaiov_ka_2
    iov_mtt <- oc1 * e_occunobs_iov * etaiov_mtt_1 + oc2 * etaiov_mtt_2
    iov_fdepot <- oc1 * e_occunobs_iov * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2

    # 3. Renal function. The Cockcroft-Gault creatinine clearance is
    #    recomputed with the cohort median weight (56 kg) and the male
    #    formula for every individual, which strips the weight and sex
    #    dependence out of the renal term so that it does not compete with
    #    the FFM allometry (Methods, adapting Holford et al.). Normalising
    #    by the cohort median CLcr,56M of 122 mL/min gives RF = 1 for the
    #    typical individual.
    clcr56m <- (140 - AGE) * 56 / (72 * CREAT)
    rf      <- clcr56m / 122

    # 4. Body-size scaling on the disposition parameters, referenced to the
    #    typical FFM of 38.6 kg.
    fsize_cl <- (FFM / 38.6)^e_ffm_cl
    fsize_vc <- (FFM / 38.6)^e_ffm_vc

    # 5. Individual disposition parameters. Total clearance is the renal arm
    #    (scaled by RF) plus the non-renal arm, with the allometric size
    #    factor applied to the sum:
    #        CL = (CLnr + RF * CLr) * Fsize
    #    (Methods equation; the parenthesisation is confirmed by the
    #    typeset equation image and by Table 2 footnote a, which states
    #    that ALL disposition parameters are allometrically scaled). The
    #    single BSV eta on CL multiplies the whole expression, matching the
    #    paper's single shared random effect across both routes.
    cl <- (exp(lcl_nonren) + rf * exp(lcl_renal)) * fsize_cl * exp(etalcl)
    vc <- exp(lvc) * fsize_vc

    # 6. Absorption. Savic transit-compartment chain feeding a first-order
    #    depot; ktr = (N + 1) / MTT is derived internally by rxode2's
    #    analytical transit() function.
    ka   <- exp(lka  + iov_ka)
    mtt  <- exp(lmtt + iov_mtt)
    nn   <- exp(lnn)
    fbio <- exp(lfdepot + iov_fdepot)

    # 7. ODE system. The dose lands on depot but the bolus is suppressed via
    #    f(depot) <- 0 so the whole dose enters through the transit-chain
    #    input rate -- the standard nlmixr2lib transit-absorption idiom
    #    (Smythe_2013_gatifloxacin.R, Vinnard_2017_rifampicin.R).
    kel <- cl / vc
    d/dt(depot)   <- transit(nn, mtt, fbio) - ka * depot
    d/dt(central) <-                          ka * depot - kel * central
    f(depot) <- 0

    # 8. Observation and error. Dose in mg and vc in L give Cc in mg/L, the
    #    unit the source reports throughout (Table 2, Table 3, Figure 2).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
