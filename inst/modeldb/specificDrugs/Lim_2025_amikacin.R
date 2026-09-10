Lim_2025_amikacin <- function() {
  description <- "One-compartment population PK model of intravenous and intramuscular amikacin in (pre)term neonates and NICU infants (Lim 2025), with an estimated power body-weight effect on clearance and volume (reference 1 kg) and an estimated power postmenstrual-age effect on clearance (reference 30 weeks). Externally validated in an independent cohort of 91 neonates."
  reference <- paste(
    "Lim CP, Candra SR, Tseng SH, Edison EP, Tan MGS, Yong MHA, Chen Y, Yeo CL.",
    "Externally validated population pharmacokinetics of amikacin and evaluation",
    "of dosage regimen based on achieved serum concentrations in neonates.",
    "Antimicrob Agents Chemother. 2025;69(8):e00818-25.",
    "doi:10.1128/aac.00818-25.",
    sep = " "
  )
  vignette <- "Lim_2025_amikacin"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "amikacin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at amikacin initiation",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Lim 2025 Methods ('Population PK analysis') states the weight effect was",
        "explored with WT 'expressed in grams', and both final equations normalize",
        "it as (WT/1,000). The reference is therefore 1,000 g = 1 kg, and the",
        "canonical kilogram-valued WT column enters the model directly as (WT / 1).",
        "Cohort weight at initiation: median 1.34 kg (IQR 0.90-1.90), range",
        "0.58-3.77 kg (Table 1, PK modeling column). Birth weight median 1.00 kg",
        "(IQR 0.87-1.44), range 0.54-3.77 kg. Time-varying in the source data set:",
        "the Methods record that a missing weight was imputed by last observation",
        "carried forward."
      ),
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age plus postnatal age) at amikacin initiation",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "WEEKS, not the register-default months. Lim 2025 Methods list the age",
        "covariates as 'gestational age (GA [weeks]), PMA (weeks), and PNA (days)',",
        "and the final CL equation is written as the bare power (PMA/30)^1.53 on a",
        "30-week reference, which is only meaningful on the week scale. The register",
        "explicitly permits a weeks-scaled PAGE for models whose source equations are",
        "written that way (see the PAGE entry in inst/references/covariate-columns.md,",
        "and Germovsek_2018_meropenem.R / Alsultan_2023_vancomycin.R).",
        "",
        "CAUTION on the cohort PMA distribution. Table 1 prints 'PMA at initiation'",
        "as median 45+6 weeks (IQR 39+5, 53+5; range 25+6 to 113+6), which is",
        "internally inconsistent with the rest of the paper and appears to carry a",
        "GA_weeks + PNA_days unit error rather than GA_weeks + PNA_days/7. Four",
        "independent checks agree that the true median PMA at initiation is near 31",
        "weeks; the reconstruction and its arithmetic are set out in full in the",
        "vignette's 'Assumptions and deviations' section. This does NOT affect the",
        "model equations, which are unambiguous -- it affects only what PMA values a",
        "user should regard as representative of the source cohort. Supply PMA on the",
        "true postmenstrual scale (GA_weeks + PNA_days/7)."
      ),
      source_name        = "PMA"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 90L,
    n_studies       = 1L,
    n_observations  = 254L,
    ga_range        = "24+0 to 41+0 weeks (median 29+0, IQR 26+0 to 31+0)",
    page_range      = "Table 1 prints 25+6 to 113+6 weeks (median 45+6); see covariateData$PAGE notes -- the printed column appears to carry a GA_weeks + PNA_days unit error, and the reconstructed true median is approximately 31 weeks",
    weight_range    = "0.58-3.77 kg at amikacin initiation (median 1.34, IQR 0.90-1.90)",
    birth_weight_range = "0.54-3.77 kg (median 1.00, IQR 0.87-1.44)",
    sex_female_pct  = 37.8,
    race_ethnicity  = c(Chinese = 62.2, Malay = 13.4, Indian = 12.2, Other = 12.2),
    disease_state   = paste(
      "Neonates and NICU infants treated with amikacin for suspected or proven",
      "septicaemia, predominantly extremely and very preterm (median gestational",
      "age 29 weeks). 25.6% small for gestational age, 11.1% with intrauterine",
      "growth restriction. None had renal impairment; neonates with congenital",
      "kidney disease, major congenital heart disease, acute kidney injury,",
      "unstable renal function, or ECMO were excluded."
    ),
    dose_range      = paste(
      "11 mg/kg every 36 h (postmenstrual age < 29 weeks) or 11 mg/kg every 24 h",
      "(postmenstrual age >= 29 weeks), given either as a 30-min intravenous",
      "infusion or as an intramuscular injection."
    ),
    creatinine_summary = "Serum creatinine at initiation 15-91 umol/L (median 43, IQR 29-57). Screened as a power covariate on CL but not retained in the final model.",
    regions         = "Single centre (Department of Neonatal and Developmental Medicine, Singapore General Hospital, Singapore). Retrospective chart review, November 2012 to October 2017.",
    notes           = paste(
      "Baseline demographics from Lim 2025 Table 1, 'PK modeling' column. The 181",
      "eligible neonates were split into a model-building set (90 neonates, 254",
      "serum concentrations) and an external-validation set (91 neonates, 280",
      "concentrations); the demographics above describe the model-building set. The",
      "external-validation set was comparable (median GA 29+0 weeks, median weight",
      "at initiation 1.39 kg, 54.9% male). External validation met the predefined",
      "criteria for MPE% (-4.99%), MAPE% (22.46%), F20 (44.64%) and F30 (64.29%)",
      "but not RMSE% (48.82%); the NPDE global test P value was < 0.01. Fit in",
      "NONMEM 7.5 with FOCE-I; bootstrap (2,000 replicates) and pcVPC performed in",
      "PsN 5.3.0. The data set was reused from Lim CP et al. 2022 (reference 35 of",
      "the paper)."
    )
  )

  ini({
    # Structural parameters. Lim 2025 reports the entire final model as two
    # printed equations in Results, 'Population PK analysis' (reproduced as
    # image aac.00818-25.m001 in the PMC supplementary file set); the paper
    # contains NO parameter-estimate table, so these two equations are the
    # complete published parameter set.
    #
    #   CL (L/H) = 0.0487(WT/1,000)^0.919 (PMA/30)^1.53 exp[eta ~ N(0, sigma^2 = 0.0188)]
    #   V (L/kg) = 0.497(WT/1,000)^0.851
    #
    # Reference subject: WT = 1,000 g = 1 kg, PMA = 30 weeks.
    lcl <- log(0.0487); label("Clearance at 1 kg and 30 weeks postmenstrual age (CL, L/h)")   # Lim 2025 Results, final model equation 1: 0.0487 L/h
    lvc <- log(0.497);  label("Volume of distribution at 1 kg (V, L)")                        # Lim 2025 Results, final model equation 2: 0.497 L

    # Power covariate exponents, both estimated (Lim 2025 Results, final model
    # equations). The paper reports no standard errors or bootstrap intervals
    # for any parameter, so none are recorded here.
    e_wt_cl   <- 0.919; label("Power exponent of (WT/1 kg) on CL (unitless)")   # Lim 2025 Results, final model equation 1
    e_page_cl <- 1.53;  label("Power exponent of (PMA/30 weeks) on CL (unitless)") # Lim 2025 Results, final model equation 1
    e_wt_vc   <- 0.851; label("Power exponent of (WT/1 kg) on V (unitless)")    # Lim 2025 Results, final model equation 2

    # Between-subject variability. The CL equation carries the term
    # exp[eta ~ N(0, sigma^2 = 0.0188)]; the variance is already on the log
    # scale, so it is used directly (0.0188 corresponds to an SD of 0.137 and
    # an approximate CV of 13.8%). The paper writes this variance with the
    # symbol sigma^2, but it is the between-subject (eta) variance -- it sits
    # inside the exponentiated random-effect term of the structural equation,
    # and the Methods separately describe BSV as "an exponential error model".
    #
    # No BSV on V. The printed V equation has no exp[eta] term, while the CL
    # equation printed immediately above it does; the omission is therefore
    # informative rather than a reporting gap. It is corroborated by the
    # Discussion ("there was a good correlation between body WT and volume of
    # distribution (V), which matches our expectation"), i.e. weight alone
    # accounted for the between-subject spread in V. No eta is declared on lvc.
    etalcl ~ 0.0188  # Lim 2025 Results, final model equation 1: eta ~ N(0, sigma^2 = 0.0188)

    # Residual unexplained variability. Lim 2025 Methods and Results state that
    # "a combined proportional-additive error model for residual unexplained
    # variability (RUV) resulted in the largest decrease in the objective
    # function value", so the STRUCTURE is reported, but neither coefficient is
    # given anywhere in the paper, and the paper has no parameter table and no
    # supplementary control stream (the PMC supplementary file set contains only
    # figure and equation images). Both coefficients are therefore encoded at
    # fixed(0) rather than invented; the combined structure is retained so a
    # user who obtains the estimates can substitute them directly. See the
    # vignette's "Assumptions and deviations" section.
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- structure reported but magnitude not published)") # Lim 2025 Results/Methods: combined proportional-additive RUV, no value reported
    addSd  <- fixed(0); label("Additive residual SD (mg/L; 0 -- structure reported but magnitude not published)")         # Lim 2025 Results/Methods: combined proportional-additive RUV, no value reported
  })

  model({
    # Individual PK parameters (Lim 2025 Results, final model equations).
    # WT is the canonical kilogram-valued column; the paper's (WT/1,000) with
    # WT in grams is numerically identical to (WT_kg / 1 kg).
    cl <- exp(lcl + etalcl) * (WT / 1)^e_wt_cl * (PAGE / 30)^e_page_cl
    vc <- exp(lvc)          * (WT / 1)^e_wt_vc

    kel <- cl / vc

    # One-compartment disposition with first-order elimination. Amikacin was
    # given either as a 30-min intravenous infusion or as an intramuscular
    # injection; the paper models both routes with this same one-compartment
    # structure and reports no absorption-rate parameter, so doses of either
    # route enter the central compartment directly. Infusion duration is not
    # hard-coded -- users specify rate or dur per dose in their event table.
    d/dt(central) <- -kel * central

    # Serum concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
