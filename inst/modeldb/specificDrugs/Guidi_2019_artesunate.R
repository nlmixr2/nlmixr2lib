Guidi_2019_artesunate <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model of artesunate (AS) and its",
    "active metabolite dihydroartemisinin (DHA) in 48 Kenyan children aged",
    "6-59 months with uncomplicated Plasmodium falciparum malaria, treated",
    "with the fixed-dose artesunate-mefloquine (ASMQ) dispersible tablet",
    "(Guidi 2019). AS and DHA are each described by a one-compartment",
    "disposition with first-order absorption into AS; AS elimination is",
    "assumed to occur exclusively by irreversible mole-for-mole conversion",
    "to DHA, so the AS elimination clearance is the metabolic formation",
    "clearance. Relative bioavailability F1 is fixed at 100% and carries all",
    "the between-subject variability, and is reduced with increasing age",
    "and on the second and third treatment days. Body weight scales DHA",
    "clearance and volume allometrically. The companion mefloquine model",
    "from the same trial is modellib('Guidi_2019_mefloquine')."
  )
  reference <- paste(
    "Guidi M, Mercier T, Aouri M, Decosterd LA, Csajka C, Ogutu B, Carn G,",
    "Kiechel JR. Population pharmacokinetics and pharmacodynamics of the",
    "artesunate-mefloquine fixed dose combination for the treatment of",
    "uncomplicated falciparum malaria in African children.",
    "Malaria Journal 2019;18:139. doi:10.1186/s12936-019-2754-6"
  )
  vignette <- "Guidi_2019_artesunate_mefloquine"

  # The source analysis was performed in MOLAR units throughout ("Molar
  # units were used for AS/DHA pharmacokinetic analyses", Methods
  # "Pharmacokinetics analysis"), with concentrations reported as nmol/mL
  # (Table 2 sigma_add, Table 4 Cmax, Fig. 1 LOQ). This file keeps amounts
  # in umol and volumes in L, so `central / vc` is umol/L, which is
  # NUMERICALLY IDENTICAL to the paper's nmol/mL (1 umol/L = 1 nmol/mL).
  # Every published concentration therefore applies directly with no
  # rescaling. Because both species are tracked on a molar basis, the
  # mole-for-mole AS -> DHA conversion needs no molecular-weight factor.
  units <- list(
    time = "h",
    dosing = "umol",
    concentration = "umol/L"
  )

  compartmentData <- list(
    depot = list(analyte = "artesunate", units = "umol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "artesunate", units = "umol", specimen = "plasma", verified = TRUE),
    central_dihydroart = list(analyte = "dihydroartemisinin", units = "umol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of the DHA (metabolite) disposition only, on the",
        "median population body weight MBW = 12.2 kg given in the Table 2",
        "footnote d: CLM_ind = CLM * (BW/MBW)^0.75 and VM_ind = VM * BW/MBW.",
        "Both exponents were FIXED, not estimated (Methods 'Covariate",
        "analysis': 'PWR the function power fixed to 0.75 for clearances and",
        "1 for volumes of distribution'). Note that AS clearance and AS",
        "central volume are NOT weight-scaled in the final model: only CLM",
        "and VM carry the footnote-d marker in Table 2, matching the",
        "Abstract statement 'Body weight affected DHA PK parameters'.",
        "Table 1 rounds the model-building cohort median weight to 12 kg;",
        "the allometric reference of 12.2 kg is the unrounded value used by",
        "the authors and is the one encoded here."
      ),
      source_name = "BW"
    ),
    AGE = list(
      description = "Chronological age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear effect on relative bioavailability F1, centred on the median",
        "population age MAGE = 2.6 years, in the fractional-deviation form",
        "given verbatim in the Table 2 footnote: '(1 + theta_age_F1",
        "(AGE-MAGE)/MAGE) with MAGE = 2.6 years'. With theta = -0.68 a child",
        "of twice the median age (5.2 y) has covariate value (5.2-2.6)/2.6 =",
        "1 and so F1 = 1 - 0.68 = 0.32, reproducing the Results statement",
        "that 'F1 is reduced by 68% upon doubling child age with respect to",
        "the population median (2.6 years)'. This fractional-deviation form",
        "is NOT a log2 or power parameterisation; both would also give a 68%",
        "drop at exactly double the median age but diverge everywhere else,",
        "and only the table footnote distinguishes them."
      ),
      source_name = "AGE"
    ),
    DAY2 = list(
      description = "Day-after-first-dose landmark indicator for the 3-day ASMQ course",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the first treatment day, the paper's day 0)",
      notes = paste(
        "1 = the record falls on the second or third treatment day (the",
        "paper's days 1 and 2), 0 = the first treatment day (day 0).",
        "Time-varying within subject. Encodes the Table 2 footnote",
        "'theta_day F1, day effect on F1 expressed as (1 + theta_day F1 Q1)",
        "with Q1 = 0 for the first treatment day; 1 for subsequent therapy",
        "days'. The authors read the treatment day as a surrogate for the",
        "rapid improvement in health produced by the first artesunate dose",
        "(Methods 'Covariate analysis'). See the model file's ini() comment",
        "for the reconciliation of theta = -0.29 with the Results sentence",
        "'29% higher in the first day of therapy'."
      ),
      source_name = "DAY"
    )
  )

  # Screened in the covariate analysis but NOT retained in the final model.
  # Documented rather than dropped so the paper's covariate screen is not
  # lost; none of these is referenced in model().
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Significantly affected the DHA volume of distribution VM in",
        "univariate analysis (dOFV = -8.8, p < 0.01) but was discarded in",
        "the complete multivariate analysis (Results 'Covariate analysis')."
      )
    ),
    PARA = list(
      description = "Plasmodium falciparum parasitaemia (asexual parasites/uL)",
      units = "parasites/uL",
      type = "continuous",
      notes = paste(
        "F1 increased with log10 baseline parasite count (dOFV = -13.2,",
        "p < 0.01) in univariate analysis, but parasitaemia and treatment",
        "day are correlated and combining them gave no further improvement",
        "(dOFV < 3.8, p > 0.05), so only treatment day (DAY2) was retained."
      )
    ),
    BILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = paste(
        "Significant on F1 in univariate analysis but dropped because the",
        "effect was poorly estimated (RSE = 155%; Results 'Covariate",
        "analysis')."
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = paste(
        "Significant on F1 in univariate analysis, but a sensitivity",
        "analysis showed the effect was driven entirely by a single patient",
        "with the highest ALT and AST, so it was not retained."
      )
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = paste(
        "Significant on F1 in univariate analysis, but a sensitivity",
        "analysis showed the effect was driven entirely by a single patient",
        "with the highest ALT and AST, so it was not retained."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 48L,
    n_studies = 1L,
    n_observations = list(AS = 117L, DHA = 134L),
    age_range = "0.6-5.0 years (median 2.6)",
    weight_range = "7-17 kg (median 12; allometric reference 12.2)",
    weight_median = "12.2 kg (MBW, the allometric reference in Table 2 footnote d)",
    age_median = "2.6 years (MAGE, the centring value for the age effect on F1)",
    sex_female_pct = 60,
    race_ethnicity = "Black African (Kenyan)",
    disease_state = "Uncomplicated Plasmodium falciparum malaria (parasite density 2000-200,000 asexual parasites/uL, fever >= 37.5 C)",
    dose_range = "One or two dispersible tablets of 25 mg artesunate / 55 mg mefloquine once daily for 3 days: one tablet for ages 6-11 months, two tablets for ages 12-59 months",
    regions = "Kenya (the intensively sampled sub-study centre)",
    notes = paste(
      "Demographics from Guidi 2019 Table 1, 'Model-building dataset",
      "(n = 48)' column. The AS/DHA model was developed only on the first",
      "50 Kenyan children enrolled in the ASMQ arm who underwent intensive",
      "sampling, of whom 48 remained after exclusions. Unlike the companion",
      "mefloquine model, the AS/DHA model received NO external validation:",
      "'Because of the very fast rate of AS and DHA elimination and the",
      "selection of the trial sampling times, an external model validation",
      "could only be performed for MQ'. AS/DHA data were sparse and heavily",
      "censored -- 71% of AS and 57% of DHA samples were below the limit of",
      "quantification (LOQ 0.005 nmol/mL for both), with a median of 1",
      "quantifiable AS and 2 quantifiable DHA samples per child. The authors",
      "note that the prediction-corrected VPCs 'evidence model",
      "misspecification' but judged the model acceptable given the paucity",
      "of AS/DHA data (Results 'Model evaluation and assessment')."
    )
  )

  ini({
    # ---- Relative bioavailability -------------------------------------
    # F1 was fixed to 100% with an estimated BSV (Table 2, row 'F1 (%)':
    # '100 fixed'). Methods: relative bioavailability was introduced 'to
    # account for dose variation with respect to the nominal value due to
    # the administration of water dispersible tablets'. Its BSV absorbed
    # all the between-subject variability previously carried by AS and DHA
    # clearance (Results: 'Inclusion of relative F1 ... explained all the
    # BSV on AS and DHA clearance'), which is why no eta sits on cl or
    # cl_dihydroart below.
    lfdepot <- fixed(log(1))
    label("Relative bioavailability of artesunate, F1 (fraction)") # Guidi 2019 Table 2: F1 = 100% fixed

    # Age effect on F1. Table 2 footnote form:
    #   F1 = F1_pop * (1 + theta_age_F1 * (AGE - MAGE) / MAGE), MAGE = 2.6 y
    e_age_fdepot <- -0.68
    label("Effect of age on artesunate F1 (fractional deviation from 2.6 years)") # Guidi 2019 Table 2: theta_age F1 = -0.68 (RSE 19%)

    # Treatment-day effect on F1. Table 2 footnote form:
    #   F1 = F1_pop * (1 + theta_day_F1 * Q1), Q1 = 0 on the first
    #   treatment day and 1 on subsequent therapy days.
    # The Results sentence reads 'F1 ... is 29% higher in the first day of
    # therapy than in the subsequent treatment days'. Read strictly against
    # the footnote, theta = -0.29 makes days 1-2 29% LOWER than day 0, so
    # day 0 is 1/0.71 = 41% higher than days 1-2. The footnote equation is
    # the authoritative encoding (the prose quotes the coefficient
    # magnitude rather than the back-transformed ratio); it is reproduced
    # exactly here.
    e_day2_fdepot <- -0.29
    label("Effect of treatment day 2-3 on artesunate F1 (multiplicative deviation)") # Guidi 2019 Table 2: theta_day F1 = -0.29 (RSE 54%)

    # ---- Artesunate disposition ---------------------------------------
    # Ka could not be estimated from these data because at most one sample
    # per child was drawn shortly after a dose, so it was fixed to 3.2 /h,
    # 'the mean of previously published estimates retrieved from papers
    # using a first-order process to depict AS absorption' (Methods,
    # references Tan 2009 [17] and Morris 2013 [18]). Tan 2009 is itself in
    # this library as modellib('Tan_2009_artesunate') (Ka = 3.85 /h fasted).
    lka <- fixed(log(3.2))
    label("First-order absorption rate constant of artesunate, Ka (1/h)") # Guidi 2019 Table 2: Ka = 3.2 /h fixed; Methods 'Structural and statistical model'

    # AS is assumed to be eliminated ONLY by irreversible conversion to
    # DHA (Methods: 'since AS is rapidly and almost completely hydrolysed
    # in DHA, its elimination was assumed to occur exclusively via
    # irreversible conversion to DHA'), so this apparent clearance is the
    # DHA formation clearance. Not weight-scaled in the final model.
    lcl <- log(146)
    label("Apparent artesunate clearance by conversion to DHA, CL/F (L/h)") # Guidi 2019 Table 2: CL = 146 L/h (RSE 20%)
    lvc <- log(139)
    label("Apparent artesunate central volume of distribution, VC/F (L)") # Guidi 2019 Table 2: VC = 139 L (RSE 23%)

    # ---- Dihydroartemisinin disposition -------------------------------
    # Typical values at the allometric reference weight MBW = 12.2 kg
    # (Table 2 footnote d).
    lcl_dihydroart <- log(11)
    label("Apparent DHA clearance at 12.2 kg, CLM/F (L/h)") # Guidi 2019 Table 2: CLM = 11 L/h (RSE 15%)
    lvc_dihydroart <- log(11)
    label("Apparent DHA volume of distribution at 12.2 kg, VM/F (L)") # Guidi 2019 Table 2: VM = 11 L (RSE 20%)

    # Allometric exponents, both FIXED by the authors rather than
    # estimated (Methods 'Covariate analysis': 'PWR the function power
    # fixed to 0.75 for clearances and 1 for volumes of distribution',
    # citing Holford 1996 [19]).
    e_wt_cl_dihydroart <- fixed(0.75)
    label("Allometric exponent on DHA clearance (unitless)") # Guidi 2019 Table 2 footnote d; Methods 'Covariate analysis'
    e_wt_vc_dihydroart <- fixed(1)
    label("Allometric exponent on DHA volume of distribution (unitless)") # Guidi 2019 Table 2 footnote d; Methods 'Covariate analysis'

    # ---- Between-subject variability ----------------------------------
    # Methods: 'Exponential errors were assumed to capture BSV in all the
    # pharmacokinetic parameters', and Table 2 reports BSV as a percentage.
    # Converted to the variance scale with omega^2 = log(CV^2 + 1):
    #   F1  CV 56% -> log(1 + 0.56^2) = 0.2727715
    #   VM  CV 60% -> log(1 + 0.60^2) = 0.3074847
    # The CV (rather than the naive sqrt(omega^2)) reading is confirmed by
    # the companion mefloquine model, where it reproduces the published
    # AUC prediction interval to within 0.1%; see the vignette.
    # Table 2 shows BSV on F1 and VM only: BSV on VC 'did not improve data
    # description' and the F1 BSV absorbed the CL and CLM variability.
    etalfdepot ~ 0.2727715 # Guidi 2019 Table 2: BSV on F1 = 56% (RSE 16%), bootstrap 53% (CI95 29 to 70)
    etalvc_dihydroart ~ 0.3074847 # Guidi 2019 Table 2: BSV on VM = 60% (RSE 29%), bootstrap 49% (CI95 16 to 75)

    # ---- Residual variability -----------------------------------------
    # Combined proportional + additive on each analyte. Table 2 reports
    # sigma_prop as a CV% and sigma_add in nmol/mL, which equals the
    # umol/L concentration unit of this file exactly.
    propSd <- 0.79
    label("Proportional residual SD for artesunate (fraction)") # Guidi 2019 Table 2: sigma_prop,AS = 79% CV (RSE 8%)
    addSd <- 0.0023
    label("Additive residual SD for artesunate (umol/L = nmol/mL)") # Guidi 2019 Table 2: sigma_add,AS = 0.0023 nmol/mL (RSE 3%)
    propSd_dihydroart <- 0.60
    label("Proportional residual SD for dihydroartemisinin (fraction)") # Guidi 2019 Table 2: sigma_prop,DHA = 60% CV (RSE 9%)
    addSd_dihydroart <- 0.0042
    label("Additive residual SD for dihydroartemisinin (umol/L = nmol/mL)") # Guidi 2019 Table 2: sigma_add,DHA = 0.0042 nmol/mL (RSE 13%)
  })

  model({
    # 1. Relative bioavailability, carrying the age and treatment-day
    #    effects and all of the between-subject variability.
    fdepot <- exp(lfdepot + etalfdepot) *
      (1 + e_age_fdepot * (AGE - 2.6) / 2.6) *
      (1 + e_day2_fdepot * DAY2)

    # 2. Individual parameters. AS disposition is not weight-scaled; DHA
    #    disposition is, on MBW = 12.2 kg.
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    cl_dihydroart <- exp(lcl_dihydroart) * (WT / 12.2)^e_wt_cl_dihydroart
    vc_dihydroart <- exp(lvc_dihydroart + etalvc_dihydroart) *
      (WT / 12.2)^e_wt_vc_dihydroart

    # 3. Micro-constants. AS has no elimination pathway other than
    #    conversion to DHA, so kel is entirely a formation rate constant.
    kel <- cl / vc
    kel_dihydroart <- cl_dihydroart / vc_dihydroart

    # 4. ODE system. Both states are in umol, so the mole-for-mole
    #    AS -> DHA conversion transfers amount one-for-one with no
    #    molecular-weight factor.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(central_dihydroart) <- kel * central - kel_dihydroart * central_dihydroart

    # 5. Bioavailability
    f(depot) <- fdepot

    # 6. Observations. umol/L is numerically identical to the paper's
    #    reported nmol/mL.
    Cc <- central / vc
    Cc_dihydroart <- central_dihydroart / vc_dihydroart

    Cc ~ add(addSd) + prop(propSd)
    Cc_dihydroart ~ add(addSd_dihydroart) + prop(propSd_dihydroart)
  })
}
