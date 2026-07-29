ArabAlameddine_2012_raltegravir <- function() {
  description <- "Two-compartment first-order-absorption population PK model for oral raltegravir (RAL) in 145 HIV-positive adults and 19 healthy volunteers, with two HIV-status-specific absorption rate constants (ka HIV+ slower than HIV-), HIV-status-specific proportional residual error, a fixed reference bioavailability F=1 for healthy volunteers, and an estimated relative bioavailability for HIV+ subjects modified linearly by sex (female +55%), atazanavir coadministration (+39%), and total bilirubin centered at 30 umol/L (+36% per doubling), plus a -59% race effect on the central volume of distribution for Caucasian relative to non-Caucasian subjects (Arab-Alameddine 2012)."
  reference <- "Arab-Alameddine M, Fayet-Mello A, Lubomirov R, Neely M, di Iulio J, Owen A, Boffito M, Cavassini M, Gunthard HF, Rentsch K, Buclin T, Aouri M, Telenti A, Decosterd LA, Rotger M, Csajka C, and the Swiss HIV Cohort Study Group. Population Pharmacokinetic Analysis and Pharmacogenetics of Raltegravir in HIV-Positive and Healthy Individuals. Antimicrob Agents Chemother. 2012;56(6):2959-2966. doi:10.1128/AAC.05424-11"
  vignette <- "ArabAlameddine_2012_raltegravir"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    HIV_POS = list(
      description        = "HIV-1 antibody-positive cohort indicator at study entry",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HIV-negative healthy volunteers)",
      notes              = "Time-fixed per subject. Drives the HIV-status-specific absorption rate (ka HIV+ = 0.21 1/h vs ka HIV- = 0.65 1/h), the relative bioavailability anchor (F_HIV- = 1 fixed, F_HIV+ = 0.75 typical with covariate modifiers), and the proportional residual error magnitude (60% CV in HIV+ vs 83.3% CV in HIV-). Source column name in the paper Table 2 is the implicit cohort-membership indicator distinguishing the 19 healthy volunteers (HIV-) from the 145 HIV-infected patients (HIV+).",
      source_name        = "HIV"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Time-fixed per subject. Linear additive effect on the HIV+ relative bioavailability (theta_female = 0.55, +55% relative to male reference); 65% higher RAL exposure in females per Discussion paragraph 4. Source-paper Table 2 covariate symbol theta_female with values 1 = female, 0 = male; same orientation as the canonical.",
      source_name        = "female"
    ),
    CONMED_ATAZANAVIR = list(
      description        = "Concomitant atazanavir (HIV protease inhibitor) coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant atazanavir)",
      notes              = "Time-fixed per subject within the analysis window. Linear additive effect on the HIV+ relative bioavailability (theta_ATV = 0.39, +39% relative to non-ATV reference) consistent with ATV-mediated UGT1A1 inhibition of raltegravir glucuronidation (Discussion paragraph 2). 11 of 145 HIV+ subjects in Table 1 were on concomitant atazanavir.",
      source_name        = "ATV"
    ),
    TBILI = list(
      description        = "Total serum bilirubin concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (single baseline assay per the paper). Reference centring value 30 umol/L (the cohort mean per Methods 'Covariate model'; median 12, range 5-91 umol/L per Table 1). Linear additive effect on the HIV+ relative bioavailability in normalized form (1 + theta_bilirubin * (TBILI / 30 - 1)) with theta_bilirubin = 0.36; a doubling of bilirubin from 30 to 60 umol/L yields a 36% increase in F, consistent with the paper's 'approximately 30% increase in drug levels in case of grade 1 hyperbilirubinaemia (total bilirubin > 30 umol/L)' (Results paragraph 4). The exact functional form is documented in the vignette Assumptions and deviations -- the paper text ('linear, centered on the mean') is consistent with a normalized centered linear coding given the magnitudes; a raw-additive form (1 + theta * (BIL - 30)) would yield non-physical predictions at low bilirubin and is ruled out.",
      source_name        = "BIL"
    ),
    RACE_WHITE = list(
      description        = "Caucasian race indicator (paper-dichotomized Caucasian vs non-Caucasian)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Caucasian; predominantly Black per Table 1 and Results paragraph 2)",
      notes              = "Time-fixed per subject. Linear additive effect on the central volume of distribution (theta_race = -0.59, -59% relative to non-Caucasian reference) per Results paragraph 2: 60% lower V1 in Caucasians than in other ethnicities. The Caucasian-as-effect-group / non-Caucasian-as-reference parameterisation matches the Hu 2014 bapineuzumab precedent in the canonical RACE_WHITE register entry. Source column name in the paper is 'theta_race'; the paper's RACE = 1 indicator means Caucasian, identical orientation to the canonical RACE_WHITE = 1 (White).",
      source_name        = "RACE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 164L,
    n_subjects_hiv_pos = 145L,
    n_subjects_hiv_neg = 19L,
    n_studies      = 3L,
    n_observations = 544L,
    age_range      = "18-72 years (HIV+ cohort)",
    age_median     = "48.5 years (HIV+ cohort)",
    weight_range   = "45-114 kg (HIV+ cohort)",
    weight_median  = "70 kg (HIV+ cohort)",
    height_range   = "152-194 cm (HIV+ cohort)",
    height_median  = "175 cm (HIV+ cohort)",
    sex_female_pct = 21.4,
    race_ethnicity = c(White = 91.0, Black = 6.9, Hispanic = 0.7,
                       Asian = 0.0, Other = 0.7, Unknown = 0.7),
    disease_state  = "Pooled HIV-positive adults on stable raltegravir-containing antiretroviral therapy in routine therapeutic drug monitoring (n=145 Swiss HIV Cohort Study) plus healthy adult volunteers from two pharmacokinetic-interaction trials (n=19, of which 19 in an open-label crossover with and without atazanavir, and 10 HIV-infected subjects in a cellular-disposition study).",
    dose_range     = "Oral raltegravir 400 mg twice daily (n=137 HIV+) or 800 mg once daily (n=8 HIV+); 400 mg twice daily +/- atazanavir 400 mg once daily for the healthy-volunteer crossover.",
    regions        = "Switzerland (Swiss HIV Cohort Study and University Hospital Lausanne) and the United States (USC Laboratory of Applied Pharmacokinetics, healthy volunteer study).",
    notes          = "Demographics summarised from Arab-Alameddine 2012 Table 1 (HIV+ cohort) and the Methods 'Study design and population' paragraph. Concomitant antiretroviral medications in the HIV+ cohort included ritonavir (49.7%), darunavir (50.3%), tenofovir (51.7%), etravirine (36.5%), and atazanavir (7.6%); median CD4 count 335 cells/mm^3 (29-968); median HIV RNA 130 copies/mL (0-277,000). Median total bilirubin 12 umol/L (5-91); the covariate model was centred at the mean 30 umol/L. Model fit with NONMEM v7.1 / NM-TRAN II (FOCE-I) per Methods 'Parameter estimation and selection'."
  )

  ini({
    # Structural parameters (Arab-Alameddine 2012 Table 2 'Final population pharmacokinetic parameters' column).
    # Values are apparent (oral) quantities at the HIV+ subject reference (SEXF = 0, CONMED_ATAZANAVIR = 0, TBILI = 30 umol/L) and the HIV- reference (RACE_WHITE = 0).
    lcl       <- log(60.2); label("Apparent oral clearance (CL/F, L/h)")                                   # Table 2: CL/F = 60.2 L/h, RSE 44.3%
    lvc       <- log(223);  label("Apparent central volume of distribution at non-Caucasian reference (V1/F, L)") # Table 2: V1/F = 223 L, RSE 52.6%
    lvp       <- log(113);  label("Apparent peripheral volume of distribution (V2/F, L)")                  # Table 2: V2/F = 113 L, RSE 41.0%
    lq        <- log(8.5);  label("Apparent inter-compartmental clearance (Q/F, L/h)")                     # Table 2: Q/F = 8.5 L/h, RSE 54.1%
    lka_neg   <- log(0.65); label("First-order absorption rate constant in healthy volunteers (ka HIV-, 1/h)") # Table 2: ka HIV- = 0.65 1/h, RSE 33.0%
    lka_pos   <- log(0.21); label("First-order absorption rate constant in HIV+ patients (ka HIV+, 1/h)")  # Table 2: ka HIV+ = 0.21 1/h, RSE 15.6%
    lfdepot   <- log(0.75); label("Relative bioavailability for HIV+ at reference covariates (F_HIV+, fraction; F_HIV- = 1 fixed)") # Table 2: F_HIV+ = 0.75, RSE 40.4%; F_HIV- fixed = 1

    # Covariate effects on the HIV+ relative bioavailability (Results paragraph 4 / Table 2).
    # Each enters via the linear additive form `(1 + theta * X)` for the binary covariates and
    # via the normalized centered linear form `(1 + theta * (TBILI / 30 - 1))` for total bilirubin.
    e_sexf_fdepot        <-  0.55; label("Linear effect of female sex on F_HIV+ (fraction; +55% relative to male)")                 # Table 2: theta_female = 0.55, RSE 52.4%
    e_atazanavir_fdepot  <-  0.39; label("Linear effect of concomitant atazanavir on F_HIV+ (fraction; +39% relative to non-ATV)") # Table 2: theta_ATV = 0.39, RSE 17.5%
    e_tbili_fdepot       <-  0.36; label("Linear effect of total bilirubin on F_HIV+ at normalized centring (fraction per unit fold-change above 30 umol/L)") # Table 2: theta_bilirubin = 0.36, RSE 47.9%

    # Covariate effect on the central volume of distribution (Results paragraph 2 / Table 2).
    e_white_vc  <- -0.59; label("Linear effect of Caucasian race on V1/F (fraction; -59% relative to non-Caucasian)") # Table 2: theta_race = -0.59, RSE 30.7%

    # IIV (Table 2 'omega' block reports CV%; log-normal variance = log(1 + CV^2)).
    # The paper reports diagonal omegas on V1/F, ka, and F_HIV+ with no off-diagonal covariances; the
    # final model removed IIV on CL/F (text: BSV on CL decreased from 72% to 19% after assigning IIV on
    # F and did not remain statistically relevant).
    etalvc      ~ 0.5675  # log(1 + 0.874^2); Table 2: omega_V1/F = 87.4% CV, RSE 93.4%
    etalka      ~ 0.6353  # log(1 + 0.942^2); Table 2: omega_ka = 94.2% CV, RSE 52.7% (shared between HIV+ and HIV-)
    etalfdepot  ~ 0.5605  # log(1 + 0.867^2); Table 2: omega_F = 86.7% CV, RSE 40.0% (shared between HIV+ and HIV-)

    # Residual error: HIV-status-conditional proportional model (Table 2 'sigma' block).
    propSd_neg <- 0.833; label("Proportional residual error in healthy volunteers (fraction)") # Table 2: sigma HIV- = 83.3% CV, RSE 40.9%
    propSd_pos <- 0.600; label("Proportional residual error in HIV+ patients (fraction)")     # Table 2: sigma HIV+ = 60.0% CV, RSE 36.2%
  })

  model({
    # 1. Derived covariate factors on the HIV+ relative bioavailability.
    # Linear additive form per Methods 'Covariate model'; bilirubin enters via the normalized centred
    # form (TBILI / 30 - 1) so the magnitude reported in the paper text (~36% increase per doubling of
    # bilirubin above the mean of 30 umol/L) reproduces directly. Each cov_* factor evaluates to 1 at
    # the reference covariate (SEXF = 0, CONMED_ATAZANAVIR = 0, TBILI = 30 umol/L).
    cov_sexf       <- 1 + e_sexf_fdepot       * SEXF
    cov_atazanavir <- 1 + e_atazanavir_fdepot * CONMED_ATAZANAVIR
    cov_tbili      <- 1 + e_tbili_fdepot      * (TBILI / 30 - 1)
    cov_fdepot     <- cov_sexf * cov_atazanavir * cov_tbili

    # Race effect on V1 (linear additive form; reference category is non-Caucasian).
    cov_white_vc <- 1 + e_white_vc * RACE_WHITE

    # 2. Individual PK parameters.
    # ka is HIV-status-specific (two distinct typical values) with a shared IIV magnitude.
    # The blend `(1 - HIV_POS) * lka_neg + HIV_POS * lka_pos` selects the correct typical
    # log-ka at each subject; etalka adds the shared between-subject deviation on the log scale.
    ka <- exp((1 - HIV_POS) * lka_neg + HIV_POS * lka_pos + etalka)
    cl <- exp(lcl)
    vc <- exp(lvc + etalvc) * cov_white_vc
    vp <- exp(lvp)
    q  <- exp(lq)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. ODE system (two-compartment with first-order absorption from depot).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                  k12 * central - k21 * peripheral1

    # 5. Bioavailability: F = 1 for HIV-, F = exp(lfdepot) * cov_fdepot for HIV+; shared multiplicative IIV.
    # The conditional form preserves the paper's fixed F_HIV- = 1 anchor while letting the covariate
    # modifiers act only on the HIV+ relative bioavailability term.
    f(depot) <- ((1 - HIV_POS) + HIV_POS * exp(lfdepot) * cov_fdepot) * exp(etalfdepot)

    # 6. Observation and HIV-status-conditional proportional residual error.
    # central is in mg, vc is in L -> central/vc is in mg/L; multiply by 1000 for ng/mL output consistent
    # with the paper Figure 1 / Figure 3 axes.
    Cc <- (central / vc) * 1000

    # Residual error: per-record proportional SD selected by HIV_POS (matches the Cirincione 2017
    # exenatide and Kyhl 2016 nalmefene precedents for indicator-switched per-study assay residuals).
    propSd <- propSd_neg * (1 - HIV_POS) + propSd_pos * HIV_POS
    Cc ~ prop(propSd)
  })
}
