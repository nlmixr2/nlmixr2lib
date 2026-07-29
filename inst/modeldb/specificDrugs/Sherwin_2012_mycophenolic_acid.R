Sherwin_2012_mycophenolic_acid <- function() {
  description <- paste(
    "Pediatric / adolescent enterohepatic-recirculation population PK",
    "model for mycophenolic acid (MPA) and its main inactive metabolite",
    "7-O-MPA-glucuronide (MPAG) in patients with childhood-onset systemic",
    "lupus erythematosus (cSLE) on oral mycophenolate mofetil (MMF)",
    "(Sherwin 2012; n = 19, age 10-28 years). Absorption uses a Savic",
    "2007 transit-compartment chain (number of compartments NN = 8.2",
    "estimated, mean transit time MTT = 1.1 h estimated) feeding the gut",
    "compartment, which delivers MPA into the central compartment at a",
    "fixed first-order absorption rate constant Ka = 1.5 1/h. MPA",
    "follows two-compartment disposition (CL1/F, V3/F central; CL2/F,",
    "V4/F peripheral). The fraction FM of total MPA elimination is",
    "converted to MPAG (FM fixed at 0.85; the remaining 0.15 is",
    "metabolism to AcMPAG, not contained in the model). MPAG follows",
    "one-compartment disposition with V_MPAG fixed equal to V3 MPA",
    "and apparent renal clearance CLM/F. Total MPAG elimination is",
    "partitioned into renal ((1 - FMPAG) fraction) and biliary (FMPAG",
    "= 0.65 fraction, fixed). Biliary MPAG enters a gallbladder",
    "compartment that empties to the gut during fixed meal-time",
    "windows (1-2 h and 4-6 h post-dose); only a fraction EHC = 0.35",
    "(fixed) of the emptied gallbladder content is reabsorbed (the",
    "rest is excreted in feces). Upon return to the gut compartment",
    "MPAG is reconverted to MPA and re-enters the absorption pathway",
    "via Ka, generating the characteristic secondary peak. The",
    "model() block hardcodes the BID dose interval (tau = 12 h) so",
    "the meal-time windows recur each interval; CL1/F, V3/F, CL2/F,",
    "V4/F, and CLM/F carry exponential IIV (paper-reported CV%).",
    "Inter-individual variability and residual error follow the source",
    "Table 3. No covariates entered the final model: bodyweight,",
    "age, sex, race, ethnicity, and disease duration were screened",
    "and rejected (see covariatesDataExcluded). Dose unit is",
    "MPA-mass-equivalent mg (MMF mg x 0.739 by molecular-weight",
    "ratio MPA/MMF = 320.3 / 433.5)."
  )
  reference <- paste(
    "Sherwin CMT, Sagcal-Gironella ACP, Fukuda T, Brunner HI, Vinks AA.",
    "Development of population PK model with enterohepatic circulation",
    "for mycophenolic acid in patients with childhood-onset systemic",
    "lupus erythematosus.",
    "British Journal of Clinical Pharmacology. 2012;73(5):727-740.",
    "doi:10.1111/j.1365-2125.2011.04140.x.",
    sep = " "
  )
  vignette <- "Sherwin_2012_mycophenolic_acid"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened in the covariate analysis (Sherwin 2012 Results,",
        "Covariate analysis). Allometric-scaled bodyweight (reference",
        "70 kg, fixed exponents 0.75 / 1) significantly improved a",
        "two-compartment first-order-absorption submodel (CL1/F, V3/F",
        "MPA and CL2/F MPAG) but did NOT improve the full EHC model",
        "and generally caused minimization failures, attributed by the",
        "authors to over-parameterisation. Not retained in the final",
        "model."
      )
    ),
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened as a continuous covariate and as a categorical",
        "covariate (child 2-12 / adolescent 12-21 / adult >21).",
        "Sherwin 2012 Results, Subjects and pharmacokinetics: 'there",
        "were no statistically significant differences related to age",
        "identified in the model.' Not retained in the final model."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened in the covariate analysis (Sherwin 2012 Results,",
        "Covariate analysis). 'No covariate provided a significant",
        "improvement to the EHC model.' Not retained in the final",
        "model. Note the cohort was 95% female (18 / 19 subjects)",
        "so a sex effect was effectively unidentifiable."
      )
    ),
    RACE_BLACK = list(
      description = "African American race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened in the covariate analysis (Sherwin 2012 Results,",
        "Covariate analysis; cohort distribution 58% African American /",
        "42% Caucasian per Table 1). Not retained in the final model."
      )
    ),
    ETHNICITY_HISPANIC = list(
      description = "Hispanic ethnicity indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened in the covariate analysis (Sherwin 2012 Results,",
        "Covariate analysis; cohort distribution 21% Hispanic /",
        "79% non-Hispanic per Table 1). Not retained in the final",
        "model."
      )
    ),
    DISEASE_DURATION = list(
      description = "Time since cSLE diagnosis",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened in the covariate analysis (Sherwin 2012 Results,",
        "Covariate analysis; cohort mean 3.3 +/- 3 years, range",
        "0.2-12.8 per Table 1). Not retained in the final model."
      )
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 19L,
    n_studies       = 1L,
    n_observations  = "186 MPA + 186 MPAG plasma concentrations (372 total)",
    age_range       = "10.6-28.2 years (Sherwin 2012 Table 1; children 10-12 y n=2, adolescents 12-21 y n=14, adults >21 y n=3)",
    age_mean        = "16.9 +/- 4 years (mean +/- SD)",
    age_median      = "16.5 years",
    weight_range    = "43.4-103 kg (Sherwin 2012 Table 1)",
    weight_mean     = "66.6 +/- 15 kg (mean +/- SD)",
    sex_female_pct  = 95.0,
    race_ethnicity  = c(
      `African American` = 58.0,
      Caucasian          = 42.0,
      Hispanic           = 21.0,
      `Non-Hispanic`     = 79.0
    ),
    disease_state   = "Childhood-onset systemic lupus erythematosus (cSLE), all subjects meeting the American College of Rheumatology classification criteria for SLE prior to age 16 years.",
    disease_duration = "Mean 3.3 +/- 3 years (range 0.2-12.8 years)",
    treatment_duration = "Mean MMF duration 1.5 +/- 1.3 years (range 0.14-6.4 years) at study entry; stable oral regimen for at least 3 weeks; therapy for at least 2 months.",
    dose_range      = "MMF 1973 +/- 634 mg/day (mean +/- SD; range 1000-3000 mg/day), administered orally as Cellcept; typical prescribed regimen is 600 mg/m^2 BID. Sampling on a single fasting outpatient visit with standardized meals at 1 h and 4 h post-dose.",
    sampling_window = "Pre-dose (1 h prior to morning MMF dose), 20 min, 40 min, 1, 1.5, 2, 3, 4, 6, 9 h post-dose.",
    co_medications  = "Oral prednisone n=18 (17.2 +/- 10.4 mg/day); high-dose i.v. methylprednisolone n=3; hydroxychloroquine n=17; NSAIDs n=11; antihypertensives n=8 (Sherwin 2012 Table 1).",
    regions         = "USA; recruitment at Cincinnati Children's Hospital Medical Center (Ohio) and Children's Memorial Hospital (Chicago, Illinois).",
    notes           = "Baseline serum chemistries assessed at study screening (visit 1) only, not at the PK sampling visit, and therefore were NOT included in the model: serum albumin 3.5 +/- 0.26 g/dL (n=17), AST 47 +/- 59 U/L (n=13), ALT 28 +/- 29 U/L (n=13), serum creatinine 0.69 +/- 0.17 mg/dL (n=17), urine protein:creatinine 0.33 +/- 0.34 (n=13). Renal function was a stable-renal-function inclusion criterion."
  )

  ini({
    # Final parameter estimates from Sherwin 2012 Table 3 ('Summary of
    # parameter estimates for the final enterohepatic recycling
    # pharmacokinetic model of MPA and MPAG'). Inter-individual
    # variability is exponential (Eq. 1); CV% values are translated to
    # log-normal log-scale variance via omega^2 = log(1 + CV^2). Residual
    # error is described in the Methods as a combined additive +
    # proportional model (Eq. 2) but Table 3 reports a single 'error CV%'
    # per analyte, which is encoded here as proportional residual SD;
    # the additive component is not reported and is set to 0 (see vignette
    # Assumptions and deviations).

    # ---- Structural parameters (point estimates from Table 3) ----
    lka      <- fixed(log(1.5));  label("Absorption rate constant Ka (1/h), gut -> central MPA, fixed")          # Table 3 K_A = 1.5 1/h (fixed)
    lcl      <- log(25.3);        label("Apparent total clearance of MPA CL1/F (L/h)")                            # Table 3 CL1 MPA = 25.3 L/h
    lvc      <- log(20.9);        label("Apparent central volume of MPA V3/F (L); also used as V_MPAG = V3 MPA") # Table 3 V3 MPA = 20.9 L
    lq       <- log(19.8);        label("Apparent intercompartmental clearance of MPA CL2/F (L/h)")               # Table 3 CL2 MPA = 19.8 L/h
    lvp      <- log(234);         label("Apparent peripheral volume of MPA V4/F (L)")                             # Table 3 V4 MPA = 234 L
    lcl_mpag <- log(2.5);         label("Apparent renal clearance of MPAG CLM/F (L/h); equals (1 - FMPAG) * total MPAG elimination") # Table 3 CLM MPAG = 2.5 L/h
    lmtt     <- log(1.1);         label("Mean transit time into the gut compartment MTT (h)")                     # Table 3 MTT = 1.1 h
    lnn      <- log(8.2);         label("Number of transit compartments NN (unitless, non-integer in Savic 2007)")# Table 3 n = 8.2

    # ---- Fixed metabolite-fraction and EHC parameters ----
    # FM (0.85), FMPAG (0.65), and EHC (0.35) were all reported as
    # 'fixed' in Table 3, with the values pre-set from prior literature
    # because the rate constants (k_30, k_35, k_50, k_56) and the
    # biliary-recirculation fraction were unidentifiable from the cohort
    # data. See Sherwin 2012 Results: 'These elimination rates were
    # unidentifiable by the model due to a lack of data and estimations
    # were based on previously reported values in the literature [45].'
    e_fm           <- fixed(0.85); label("Fraction of MPA elimination converted to MPAG (FM, unitless, fixed; remainder 0.15 is AcMPAG, not modelled)")      # Table 3 FM     fixed at 85%
    e_fmpag        <- fixed(0.65); label("Fraction of MPAG total elimination via biliary pathway (FMPAG, unitless, fixed; complement is renal)")             # Table 3 FMPAG  fixed at 65%
    e_ehc          <- fixed(0.35); label("Fraction of gallbladder content reaching the gut during meal-time emptying (EHC, unitless, fixed; complement is feces)") # Table 3 EHC    fixed at 35%

    # ---- Inter-individual variability (Eq. 1; CV% per Table 3) ----
    # omega^2 = log(1 + CV^2):
    #   CV = 48.6% -> 0.21201    CV = 59.2% -> 0.30048
    #   CV = 42.9% -> 0.16895    CV = 60.0% -> 0.30748
    #   CV = 55.9% -> 0.27193
    etalcl       ~ 0.21201    # Table 3 IIV CL1 MPA  = 48.6% CV -> log(1 + 0.486^2)
    etalvc       ~ 0.30048    # Table 3 IIV V3 MPA   = 59.2% CV -> log(1 + 0.592^2)
    etalq        ~ 0.16895    # Table 3 IIV CL2 MPA  = 42.9% CV -> log(1 + 0.429^2)
    etalvp       ~ 0.30748    # Table 3 IIV V4 MPA   = 60.0% CV -> log(1 + 0.600^2)
    etalcl_mpag  ~ 0.27193    # Table 3 IIV CLM MPAG = 55.9% CV -> log(1 + 0.559^2)

    # ---- Residual error (Eq. 2; Table 3 reports a single CV% per analyte) ----
    # The Methods state a combined additive + proportional model but
    # Table 3 only reports 'MPA error CV%' = 41.2% and 'MPAG error CV%'
    # = 45.4%. The single-CV reporting is consistent with the
    # proportional component dominating; the additive component is not
    # reported and is omitted here. See vignette Assumptions and
    # deviations.
    propSd       <- 0.412;        label("Proportional residual error on MPA plasma concentration (fraction)")    # Table 3 MPA  error CV = 41.2%
    propSd_mpag  <- 0.454;        label("Proportional residual error on MPAG plasma concentration (fraction)")   # Table 3 MPAG error CV = 45.4%
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual PK parameters. Inter-individual variability is
    #    exponential (log-normal) on the log-transformed parameter.
    # ------------------------------------------------------------------
    ka          <- exp(lka)
    cl          <- exp(lcl       + etalcl)
    vc          <- exp(lvc       + etalvc)
    q           <- exp(lq        + etalq)
    vp           <- exp(lvp      + etalvp)
    cl_mpag     <- exp(lcl_mpag  + etalcl_mpag)
    mtt         <- exp(lmtt)
    nn          <- exp(lnn)

    # Bare names for the fixed metabolite fractions and EHC factor.
    fm          <- e_fm
    fmpag       <- e_fmpag
    ehc         <- e_ehc

    # MPAG central volume is fixed equal to V3 MPA per Table 3:
    # 'V_M = V3 MPA'. The same etalvc IIV therefore propagates to V_MPAG.
    v_mpag      <- vc

    # ------------------------------------------------------------------
    # 2. Micro-constants.
    #    kel        = total elimination rate constant of MPA from central.
    #    kel_mpag   = total elimination rate constant of MPAG from central
    #                 (renal + biliary). CLM is the apparent renal CL of
    #                 MPAG = (1 - FMPAG) * kel_mpag_total * V_MPAG, so
    #                 kel_mpag = CLM / (V_MPAG * (1 - FMPAG)).
    # ------------------------------------------------------------------
    kel         <- cl       / vc
    k12         <- q        / vc
    k21         <- q        / vp
    kel_mpag    <- cl_mpag  / (v_mpag * (1.0 - fmpag))

    # ------------------------------------------------------------------
    # 3. Gallbladder-emptying windows. The source paper times two
    #    bile-release windows to the standardized meal times during the
    #    PK visit: 1-2 h and 4-6 h post-dose (Sherwin 2012 Results,
    #    'Gallbladder emptying was modelled to simulate two release
    #    times, 1 and 4 h post dose ... Gallbladder emptying was
    #    "turned off" at 2 h and 6 h post dose'). The packaged model
    #    hardcodes the BID dosing interval tau = 12 h so the meal-time
    #    windows recur each dosing interval, matching the source paper's
    #    typical regimen (600 mg/m^2 MMF BID). For other dosing intervals
    #    edit the `tau` constant.
    # ------------------------------------------------------------------
    tau         <- 12.0
    tpost       <- t - floor(t / tau) * tau
    in_meal_1   <- (tpost >= 1.0) * (tpost <= 2.0)
    in_meal_2   <- (tpost >= 4.0) * (tpost <= 6.0)
    k_empty     <- 5.0
    empty_rate  <- k_empty * (in_meal_1 + in_meal_2)

    # ------------------------------------------------------------------
    # 4. ODE system. State variables (mg of MPA- or MPAG-equivalent):
    #      depot         - gut compartment (Sherwin 2012 cmt 2) -- the
    #                      Savic transit chain feeds this compartment;
    #                      ka transports its contents into central MPA.
    #                      EHC-recycled mass also enters here as MPA.
    #      central       - MPA central (cmt 3)
    #      peripheral1   - MPA peripheral (cmt 4)
    #      central_mpag  - MPAG central (cmt 5)
    #      gallbladder   - MPAG gallbladder (cmt 6)
    #
    # The dose record places a bolus on `depot`; f(depot) <- 0 suppresses
    # the bolus so all input arrives via the analytical Savic transit
    # macro `transit(nn, mtt, 1)` (bioavailability = 1; dose units are
    # MPA-mass-equivalent mg, the user externally converts MMF mass by
    # MMF * 0.739 = MMF * MW_MPA / MW_MMF). The transit() input rate is
    # added to the d/dt(depot) equation and feeds the gut at chain rate
    # ktr = (nn + 1) / mtt.
    # ------------------------------------------------------------------
    d/dt(depot)        <- transit(nn, mtt, 1.0) - ka * depot + ehc * empty_rate * gallbladder
    d/dt(central)      <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1)  <-                                k12 * central - k21 * peripheral1
    # MPA -> MPAG metabolic formation: FM fraction of total MPA
    # elimination becomes MPAG (1:1 molar transfer; the model is mass-
    # based with MPA and MPAG units interchanged at this junction,
    # consistent with the source paper's choice to track both species
    # in mg/L).
    d/dt(central_mpag) <-  fm * kel * central - kel_mpag * central_mpag
    # MPAG biliary excretion to gallbladder: fmpag fraction of total
    # MPAG elimination enters the gallbladder. During meal windows,
    # empty_rate * gallbladder leaves the gallbladder at a fast rate
    # (k_empty = 5/h); only an EHC = 0.35 fraction of that emptied
    # content reaches the gut, the rest is excreted to feces.
    d/dt(gallbladder)  <-  fmpag * kel_mpag * central_mpag - empty_rate * gallbladder

    # Suppress the bolus into depot so all dose input is routed through
    # the transit chain.
    f(depot) <- 0

    # ------------------------------------------------------------------
    # 5. Observed plasma concentrations (mg/L). The same V3 MPA is used
    #    for both the MPA and MPAG central compartments per the source
    #    paper's V_M = V3 assumption.
    # ------------------------------------------------------------------
    Cc      <- central      / vc
    Cc_mpag <- central_mpag / v_mpag

    Cc      ~ prop(propSd)
    Cc_mpag ~ prop(propSd_mpag)
  })
}
