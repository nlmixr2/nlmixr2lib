Morcos_2023b_copanlisib <- function() {
  description <- "Three-compartment population PK model for intravenous copanlisib adapted from the adult phase I-III model to children and adolescents (4-21 years) with relapsed/refractory solid tumors: every clearance and volume is scaled by body surface area relative to 1.88 m^2 with a fixed exponent of 1, with the adult typical values, categorical covariate effects (rifampicin and itraconazole comedication, sex, hepatic impairment, Japanese region, CHRONOS-3 study membership), IIV and infusion-time / study-phase stratified log-additive residual error carried over unchanged"
  reference <- paste(
    "Morcos PN, Schlender J, Burghaus R, Moss J, Lloyd A, Childs BH, Macy ME,",
    "Reid JM, Chung J, Garmann D. Model-informed approach to support",
    "pediatric dosing for the pan-PI3K inhibitor copanlisib in children and",
    "adolescents with relapsed/refractory solid tumors.",
    "Clin Transl Sci. 2023;16(7):1197-1209.",
    "doi:10.1111/cts.13523.",
    "Adult backbone parameters: Morcos PN, Moss J, Austin R, Hiemeyer F,",
    "Zinzani PL, Beckert V, Mongay Soler L, Childs BH, Garmann D. Copanlisib",
    "population pharmacokinetics from phase I-III studies and",
    "exposure-response relationships in combination with rituximab.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(11):1666-1686.",
    "doi:10.1002/psp4.13000.",
    sep = " "
  )
  vignette <- "Morcos_2023b_copanlisib"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = "n/a (power scaling relative to the typical adult value 1.88 m^2)",
      notes = paste(
        "The pediatric adaptation of Morcos 2023 CTS (Methods, 'Population PK",
        "modeling'): 'The adult popPK model was adapted through allometric",
        "scaling using BSA (centered around a typical adult BSA value of 1.88",
        "m2) with an exponent of 1 on all rate (including clearance) and volume",
        "parameters.' The exponent was fixed, not estimated; exponents such as",
        "0.75 were tested as sensitivity analyses and gave no improvement",
        "(Results, 'PopPK modeling'). The BSA formula is unspecified. Pediatric",
        "cohort median 1.37 m^2 (range 0.64-2.44). The adult backbone itself",
        "carries no body-size covariate (body weight was screened and not",
        "retained in Morcos 2023 CPT:PSP), so this term is the only size",
        "effect.",
        sep = " "
      ),
      source_name = "BSA"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Adult backbone covariate, carried unchanged into the pediatric",
        "adaptation. Retained on BOTH clearance and central volume in Morcos",
        "2023 CPT:PSP ('females had 42.9% lower V1 than males' and 'females had",
        "16.7% lower CL than males'; Table 2 footnotes d and e name a male",
        "reference patient). The pediatric Figure 5 simulation used 250 male",
        "and 250 female virtual patients per age.",
        sep = " "
      ),
      source_name = "Sex"
    ),
    HEPIMP = list(
      description = "Any hepatic impairment (NCI ODWG group >= 2) versus normal hepatic function",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal hepatic function, NCI ODWG group 1)",
      notes = paste(
        "Adult backbone covariate (Morcos 2023 CPT:PSP Table 2 theta_NCICL:",
        "'influence of NCI for any hepatic impairment category relative to",
        "normal hepatic function on clearance'). The pediatric paper does not",
        "report hepatic function for its cohort.",
        sep = " "
      ),
      source_name = "NCI"
    ),
    REGION_JAPAN = list(
      description = "Enrolled at a Japanese study site",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all non-Japan regions)",
      notes = paste(
        "Adult backbone covariate (Morcos 2023 CPT:PSP: 'patients from Japan",
        "had 20.4% lower CL than patients from other regions'). The pediatric",
        "paper does not report enrollment regions for its cohort.",
        sep = " "
      ),
      source_name = "Region = Japan"
    ),
    CONMED_RIFAMPICIN = list(
      description = "Concomitant rifampicin (strong CYP3A4 inducer) coadministration",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no rifampicin coadministration)",
      notes = paste(
        "Adult backbone covariate, time-varying; set to 1 only during",
        "rifampicin comedication in the adult DDI study 16270 (Morcos 2023",
        "CPT:PSP Table S2). Retained on both CL (+191%) and V1 (+108%). 0 for",
        "the pediatric study.",
        sep = " "
      ),
      source_name = "Comedication with rifampin in study 16270"
    ),
    CONMED_ITRACONAZOLE = list(
      description = "Concomitant itraconazole (strong CYP3A4 inhibitor) coadministration",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no itraconazole coadministration)",
      notes = paste(
        "Adult backbone covariate, time-varying; set to 1 only during",
        "itraconazole comedication in the adult DDI study 16270 (Morcos 2023",
        "CPT:PSP Table S2). Retained on CL only (-36.1%). 0 for the pediatric",
        "study.",
        sep = " "
      ),
      source_name = "Comedication with itraconazole in study 16270"
    ),
    STUDY_CHRONOS3 = list(
      description = "Enrolled in Bayer study 17067 (CHRONOS-3, NCT02367040)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (any other study, including the pediatric phase I study 19176)",
      notes = paste(
        "Adult backbone covariate with two roles: a -18.4% effect on clearance",
        "(Morcos 2023 CPT:PSP Table 2, theta_17067CL) and selection of the",
        "phase III residual magnitude beyond the first 20 min of an infusion.",
        "0 for the pediatric phase I study 19176, which therefore uses the",
        "phase I/II residual stratum.",
        sep = " "
      ),
      source_name = "Study 17067"
    )
  )

  compartmentData <- list(
    central = list(analyte = "copanlisib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "copanlisib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "copanlisib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 30,
    n_studies = 1,
    n_observations = 160,
    age_range = "4-21 years",
    age_median = "13 years",
    bsa_range = "0.64-2.44 m^2",
    bsa_median = "1.37 m^2",
    disease_state = "children and adolescents with relapsed/refractory solid tumors or lymphoma for which there is no standard curative treatment",
    dose_range = "28 mg/m2 (n = 24) or 35 mg/m2 (n = 7) as an intravenous infusion on days 1, 8 and 15 of a 28-day cycle",
    notes = paste(
      "Morcos 2023 CTS Results, 'PopPK modeling': 31 pediatric patients in",
      "the open-label phase I dose-escalation study 19176 (Rolling-6 design)",
      "gave 191 samples; 30 pre-dose BLQ samples and 1 kinetically",
      "implausible sample were removed, leaving 160 observations in 30",
      "patients. Samples were drawn pre-infusion, at 1-1.25 h, 1.5-3 h (age",
      ">= 6 years only) and 22-24 h after the start of the infusion on cycle 1",
      "days 1 and 15; LC/MS LLOQ 0.5 ng/mL. The pediatric data were used only",
      "to evaluate the adapted model: no parameter was re-estimated. Every",
      "typical value, covariate effect, IIV and residual magnitude below is",
      "the adult estimate of Morcos 2023 CPT:PSP (712 adults, nine phase I-III",
      "studies). The sex split of the pediatric cohort is not reported.",
      sep = " "
    )
  )

  ini({
    # --------------------------------------------------------------------
    # Structural parameters for the adult reference patient. Morcos 2023 CTS
    # (Methods, 'Population PK modeling') used 'the comprehensive copanlisib
    # adult popPK model which has been previously established based on data
    # collected from phase I-III studies' (its reference 9, the 2022 AACR
    # poster whose figures Morcos 2023 CPT:PSP reproduces) 'without further
    # modification or any parameter re-estimation'. The values are therefore
    # those of Morcos 2023 CPT:PSP Table 2, 'Fixed effects (THETA)'. Reference
    # patient: male, non-Japanese phase I / II study, no rifampicin or
    # itraconazole, normal hepatic function, BSA 1.88 m^2.
    #
    # NONMEM V1/V2/V3 and Q2/Q3 map onto V1 -> vc, V2 -> vp, V3 -> vp2,
    # Q2 -> q, Q3 -> q2.
    # --------------------------------------------------------------------
    lcl  <- log(22.2) ; label("Clearance for the adult reference patient (L/h)")                    # Morcos 2023 CPT:PSP Table 2: CL_pop = 22.2 L/h (RSE 3.18%)
    lvc  <- log(92.1) ; label("Central volume of distribution for the adult reference patient (L)") # Morcos 2023 CPT:PSP Table 2: V1_pop = 92.1 L (RSE 7.47%)
    lq   <- log(79.3) ; label("Intercompartmental clearance, central to peripheral1 (L/h)")         # Morcos 2023 CPT:PSP Table 2: Q2 = 79.3 L/h (RSE 1.58%)
    lvp  <- log(508)  ; label("First peripheral volume of distribution (L)")                        # Morcos 2023 CPT:PSP Table 2: V2 = 508 L (RSE 2.54%)
    lq2  <- log(7.34) ; label("Intercompartmental clearance, central to peripheral2 (L/h)")         # Morcos 2023 CPT:PSP Table 2: Q3 = 7.34 L/h (RSE 6.96%)
    lvp2 <- log(522)  ; label("Second peripheral volume of distribution (L)")                       # Morcos 2023 CPT:PSP Table 2: V3 = 522 L (RSE 4.26%)

    # --------------------------------------------------------------------
    # Body-surface-area scaling, the pediatric adaptation. Morcos 2023 CTS
    # Methods: 'allometric scaling using BSA (centered around a typical adult
    # BSA value of 1.88 m2) with an exponent of 1 on all rate (including
    # clearance) and volume parameters', held constant rather than estimated.
    # --------------------------------------------------------------------
    e_bsa_cl  <- fixed(1) ; label("Power exponent on (BSA / 1.88 m^2) for CL (unitless)") # Morcos 2023 CTS Methods: exponent of 1 on all rate and volume parameters
    e_bsa_vc  <- fixed(1) ; label("Power exponent on (BSA / 1.88 m^2) for V1 (unitless)") # Morcos 2023 CTS Methods: exponent of 1 on all rate and volume parameters
    e_bsa_q   <- fixed(1) ; label("Power exponent on (BSA / 1.88 m^2) for Q2 (unitless)") # Morcos 2023 CTS Methods: exponent of 1 on all rate and volume parameters
    e_bsa_vp  <- fixed(1) ; label("Power exponent on (BSA / 1.88 m^2) for V2 (unitless)") # Morcos 2023 CTS Methods: exponent of 1 on all rate and volume parameters
    e_bsa_q2  <- fixed(1) ; label("Power exponent on (BSA / 1.88 m^2) for Q3 (unitless)") # Morcos 2023 CTS Methods: exponent of 1 on all rate and volume parameters
    e_bsa_vp2 <- fixed(1) ; label("Power exponent on (BSA / 1.88 m^2) for V3 (unitless)") # Morcos 2023 CTS Methods: exponent of 1 on all rate and volume parameters

    # --------------------------------------------------------------------
    # Adult categorical covariate effects, Morcos 2023 CPT:PSP Table 2, in
    # that paper's form P = P_TV * (1 + N1 * theta_1 + ...) * exp(eta): each
    # theta is a fractional change.
    # --------------------------------------------------------------------
    e_conmed_rifampicin_cl   <-  1.91   ; label("Fractional change in CL with rifampicin coadministration")   # Morcos 2023 CPT:PSP Table 2: theta_RIFCL = 1.91 (RSE 3.56%)
    e_conmed_itraconazole_cl <- -0.361  ; label("Fractional change in CL with itraconazole coadministration") # Morcos 2023 CPT:PSP Table 2: theta_ITRACL = -0.361 (RSE 5.07%)
    e_study_chronos3_cl      <- -0.184  ; label("Fractional change in CL for CHRONOS-3 (study 17067)")        # Morcos 2023 CPT:PSP Table 2: theta_17067CL = -0.184 (RSE 17.6%)
    e_sexf_cl                <- -0.167  ; label("Fractional change in CL for female sex")                     # Morcos 2023 CPT:PSP Table 2: theta_SEXCL = -0.167 (RSE 17.2%)
    e_hepimp_cl              <- -0.192  ; label("Fractional change in CL with any hepatic impairment")        # Morcos 2023 CPT:PSP Table 2: theta_NCICL = -0.192 (RSE 17.4%)
    e_region_japan_cl        <- -0.204  ; label("Fractional change in CL for patients from Japan")            # Morcos 2023 CPT:PSP Table 2: theta_JAPCL = -0.204 (RSE 29.3%)

    e_sexf_vc                <- -0.429  ; label("Fractional change in V1 for female sex")                     # Morcos 2023 CPT:PSP Table 2: theta_SEXV1 = -0.429 (RSE 14.6%)
    e_conmed_rifampicin_vc   <-  1.08   ; label("Fractional change in V1 with rifampicin coadministration")   # Morcos 2023 CPT:PSP Table 2: theta_RIFV1 = 1.08 (RSE 39.0%)

    # --------------------------------------------------------------------
    # Inter-individual variability, Morcos 2023 CPT:PSP Table 2 (omega^2
    # variances, diagonal, on CL and V1 only). Morcos 2023 CTS Figure 5
    # legend: 'CL and V1 were derived by sampling random effects from the
    # adjusted popPK model'.
    # --------------------------------------------------------------------
    etalcl ~ 0.124  # Morcos 2023 CPT:PSP Table 2: omega^2 on CL_pop = 0.124, i.e. CV 36.3%
    etalvc ~ 0.846  # Morcos 2023 CPT:PSP Table 2: omega^2 on V1_pop = 0.846, i.e. CV 115%

    # --------------------------------------------------------------------
    # Residual error, Morcos 2023 CPT:PSP Table 2 and footnote g:
    # log-transformed observations and predictions with an additive error,
    # i.e. lnorm() with a log-scale SD = sqrt(sigma^2). Three strata: the
    # first 20 min of an infusion, then phase I/II versus phase III.
    # --------------------------------------------------------------------
    expSd_early   <- 2.2583 ; label("Log-scale residual SD, first 20 min of an infusion")         # Morcos 2023 CPT:PSP Table 2: sigma^2 = 5.10 -> sqrt(5.10) = 2.2583
    expSd_phase12 <- 0.4195 ; label("Log-scale residual SD, phase I/II studies beyond 20 min")    # Morcos 2023 CPT:PSP Table 2: sigma^2 = 0.176 -> sqrt(0.176) = 0.4195
    expSd_phase3  <- 0.7950 ; label("Log-scale residual SD, phase III study 17067 beyond 20 min") # Morcos 2023 CPT:PSP Table 2: sigma^2 = 0.632 -> sqrt(0.632) = 0.7950
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters: adult categorical covariate factors times
    #    the BSA power term. BSA = 1.88 m^2 with every indicator 0 recovers
    #    the adult reference patient of Morcos 2023 CPT:PSP exactly.
    # ------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (BSA / 1.88)^e_bsa_cl *
      (1 + e_conmed_rifampicin_cl   * CONMED_RIFAMPICIN)   *
      (1 + e_conmed_itraconazole_cl * CONMED_ITRACONAZOLE) *
      (1 + e_study_chronos3_cl      * STUDY_CHRONOS3)      *
      (1 + e_sexf_cl                * SEXF)                *
      (1 + e_hepimp_cl              * HEPIMP)              *
      (1 + e_region_japan_cl        * REGION_JAPAN)

    vc <- exp(lvc + etalvc) * (BSA / 1.88)^e_bsa_vc *
      (1 + e_sexf_vc              * SEXF) *
      (1 + e_conmed_rifampicin_vc * CONMED_RIFAMPICIN)

    q   <- exp(lq)   * (BSA / 1.88)^e_bsa_q
    vp  <- exp(lvp)  * (BSA / 1.88)^e_bsa_vp
    q2  <- exp(lq2)  * (BSA / 1.88)^e_bsa_q2
    vp2 <- exp(lvp2) * (BSA / 1.88)^e_bsa_vp2

    # ------------------------------------------------------------------
    # 2. Micro-constants.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ------------------------------------------------------------------
    # 3. Three-compartment IV disposition with first-order elimination from
    #    the central compartment. Copanlisib is given only as an intravenous
    #    infusion, so dose records target `central` with a rate or duration.
    # ------------------------------------------------------------------
    d/dt(central)     <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <-   k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-   k13 * central - k31 * peripheral2

    # ------------------------------------------------------------------
    # 4. Observation and stratified log-additive residual error.
    #
    #    Doses in mg and volumes in L give Cc in mg/L; the paper reports
    #    concentrations in ug/L and AUC in ug*h/L (1 mg/L = 1000 ug/L).
    #    tad() < 20/60 h selects the 'first 20 min of the infusion' stratum
    #    because a copanlisib dose record starts the infusion.
    # ------------------------------------------------------------------
    Cc <- central / vc

    earlyInfusion <- tad() < (20 / 60)

    expSd <- expSd_early * earlyInfusion +
      (1 - earlyInfusion) * (expSd_phase12 * (1 - STUDY_CHRONOS3) +
                               expSd_phase3 * STUDY_CHRONOS3)

    Cc ~ lnorm(expSd)
  })
}
