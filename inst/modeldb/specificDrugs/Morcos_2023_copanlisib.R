Morcos_2023_copanlisib <- function() {
  description <- "Three-compartment population PK model for intravenous copanlisib in adults with advanced solid tumors or non-Hodgkin lymphoma, pooled across nine phase I-III studies (n = 712), with categorical covariate effects of rifampicin and itraconazole comedication, sex, hepatic impairment, Japanese region and CHRONOS-3 study membership on clearance and central volume, and an infusion-time / study-phase stratified log-additive residual error"
  reference   <- paste(
    "Morcos PN, Moss J, Austin R, Hiemeyer F, Zinzani PL, Beckert V,",
    "Mongay Soler L, Childs BH, Garmann D.",
    "Copanlisib population pharmacokinetics from phase I-III studies and",
    "exposure-response relationships in combination with rituximab.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(11):1666-1686.",
    "doi:10.1002/psp4.13000.",
    sep = " "
  )
  vignette    <- "Morcos_2023_copanlisib"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Retained on BOTH clearance and central volume. Morcos 2023 Results",
        "'PopPK meta-analyses': 'females had 42.9% lower V1 than males' and",
        "'females had 16.7% lower CL than males'. The Table 2 footnotes d and e",
        "both name a male reference patient, so SEXF = 1 for female and the",
        "fractional effects are negative. Note the paper's Table S2 codes sex",
        "as 1 = male for the separate exposure-response covariate screen; that",
        "opposite coding applies only to the ER analysis, which is not encoded",
        "here (see the vignette Errata).",
        sep = " "
      ),
      source_name        = "Sex"
    ),
    HEPIMP = list(
      description        = "Any hepatic impairment (NCI ODWG group >= 2) versus normal hepatic function",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function, NCI ODWG group 1)",
      notes              = paste(
        "Morcos 2023 Table 2 describes theta_NCICL as the 'influence of NCI for",
        "any hepatic impairment category relative to normal hepatic function on",
        "clearance', i.e. mild, moderate and severe are pooled into a single",
        "indicator against the normal-function reference. Table 1 reports the",
        "pooled cohort as 84.1% normal, 14.9% mild, 0.4% moderate, 0.6% severe,",
        "so HEPIMP = 1 in 15.9% of the 712 patients.",
        sep = " "
      ),
      source_name        = "NCI"
    ),
    REGION_JAPAN = list(
      description        = "Enrolled at a Japanese study site",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (all non-Japan regions: North America, Europe, mainland China, other Asia, other)",
      notes              = paste(
        "Morcos 2023 Results: 'patients from Japan had 20.4% lower CL than",
        "patients from other regions'. Table 1 gives 61 of 712 (8.6%) Japanese",
        "patients. Table S2 lists Europe, North America, mainland China, Japan",
        "and Others as five separate dichotomous region covariates screened;",
        "only Japan survived backward elimination, so the remaining four are",
        "recorded in covariatesDataExcluded.",
        sep = " "
      ),
      source_name        = "Region = Japan"
    ),
    CONMED_RIFAMPICIN = list(
      description        = "Concomitant rifampicin (strong CYP3A4 inducer) coadministration",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no rifampicin coadministration)",
      notes              = paste(
        "Time-varying. Morcos 2023 Table S2: 'Covariate is set to 1 during",
        "periods of rifampin comedication in study 16270 and set to 0 at all",
        "other times for patients in study 16270, and set to 0 at all times for",
        "patients not in study 16270.' Study 16270 (NCT02253420) was the",
        "dedicated DDI study; arm B received rifampicin 600 mg once daily on",
        "cycle 1 days 10-21 (Table S1). Retained on both CL (+191%) and V1",
        "(+108%). The paper uses the US name 'rifampin'; the canonical column",
        "uses the INN 'rifampicin'.",
        sep = " "
      ),
      source_name        = "Comedication with rifampin in study 16270"
    ),
    CONMED_ITRACONAZOLE = list(
      description        = "Concomitant itraconazole (strong CYP3A4 inhibitor) coadministration",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no itraconazole coadministration)",
      notes              = paste(
        "Time-varying, defined exactly as CONMED_RIFAMPICIN but for arm A of",
        "study 16270, which received itraconazole 200 mg twice daily on cycle 1",
        "day 12 and 200 mg once daily on cycle 1 days 13-21 (Morcos 2023",
        "Tables S1 and S2). Retained on CL only (-36.1%).",
        sep = " "
      ),
      source_name        = "Comedication with itraconazole in study 16270"
    ),
    STUDY_CHRONOS3 = list(
      description        = "Enrolled in Bayer study 17067 (CHRONOS-3, NCT02367040)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any of the eight phase I / phase II studies in the pooled analysis)",
      notes              = paste(
        "Carries two distinct roles in this model. (1) Structural: a -18.4%",
        "effect on clearance (Morcos 2023 Table 2, theta_17067CL). The authors",
        "note this 'could not be robustly assigned to co-administration with",
        "rituximab due to differences in study design and patient population',",
        "and Table S2 records that rituximab comedication is entirely",
        "confounded with study 17067 because rituximab was given only in that",
        "study. (2) Residual error: selects the phase III residual magnitude",
        "(CV 93.8%) instead of the pooled phase I/II magnitude (CV 43.9%) for",
        "samples drawn more than 20 min after the start of an infusion",
        "(Table 2, SIGMA rows). 289 of the 712 pooled PK patients (40.6%) are",
        "from study 17067 (Table 1).",
        sep = " "
      ),
      source_name        = "Study 17067"
    )
  )

  # Every covariate below was screened in the Morcos 2023 full stepwise
  # forward-inclusion / backward-elimination procedure (Table S2) but did NOT
  # survive into the final model of Table 2, so none is referenced in model().
  # They are recorded here so the provenance of the paper's covariate search is
  # not lost. Body weight is the notable omission: allometric scaling was
  # investigated (Methods, "Allometric scaling used for investigation of body
  # weight") but not retained, and the Discussion instead attributes the lower
  # clearance in Japanese patients to "general differences in body weight".
  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened as an allometric term P = Ptv * (BW/BWmed)^theta (Morcos 2023 Methods); not retained. Pooled median 70.0 kg (range 41.1-165), Table 1."
    ),
    AGE = list(
      description = "Baseline age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened (Table S2); not retained. Pooled median 63 years (range 20-91), Table 1."
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened (Table S2, reported there in g/dL); not retained. Pooled median 4.14 g/dL = 41.4 g/L (range 1.6-6.5 g/dL), Table 1."
    ),
    CRCL = list(
      description = "Baseline estimated glomerular filtration rate",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened as eGFR (Table S2); not retained. Pooled median 87.9 mL/min (range 13.64-155.91), Table 1."
    ),
    REGION_EUROPE = list(
      description = "Enrolled at a European study site",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (Table S2); not retained. 351 of 712 patients (49.3%), Table 1."
    ),
    REGION_USA = list(
      description = "Enrolled at a North American study site",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as 'Region = North America' (Table S2); not retained. 147 of 712 patients (20.6%), Table 1."
    ),
    REGION_CHINA = list(
      description = "Enrolled at a mainland China study site",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (Table S2); not retained. 70 of 712 patients (9.8%), Table 1."
    ),
    REGION_ROW = list(
      description = "Enrolled outside Europe, North America, mainland China and Japan",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as 'Region = Others' (Table S2); not retained. Table 1 splits this into other Asia (28, 3.9%) and other (55, 7.7%)."
    )
  )

  compartmentData <- list(
    central     = list(analyte = "copanlisib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "copanlisib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "copanlisib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 712,
    n_studies      = 9,
    n_observations = 5958,
    age_range      = "20-91 years",
    age_median     = "63 years",
    weight_range   = "41.1-165 kg",
    weight_median  = "70.0 kg",
    sex_female_pct = 52.4,
    disease_state  = "advanced solid tumors, aggressive non-Hodgkin lymphoma, or indolent non-Hodgkin lymphoma (chiefly relapsed follicular lymphoma); one phase I study also enrolled healthy participants alongside hepatic- and renal-impairment cohorts",
    dose_range     = "0.1-1.2 mg/kg or 12-60 mg flat, given as a 1-h intravenous infusion on days 1, 8 and 15 of a 28-day cycle (3 weeks on / 1 week off); the approved and phase III regimen is 60 mg flat",
    regions        = "Europe (49.3%), North America (20.6%), mainland China (9.8%), Japan (8.6%), other Asia (3.9%), other (7.7%)",
    hepatic_function = "84.1% normal, 14.9% mild, 0.4% moderate, 0.6% severe (NCI ODWG)",
    renal_function   = "46.8% normal, 40.2% mild, 11.8% moderate, 1.3% severe (NCI criteria)",
    albumin_median   = "4.14 g/dL (range 1.6-6.5)",
    egfr_median      = "87.9 mL/min (range 13.64-155.91)",
    co_medication    = "rituximab 375 mg/m2 in study 17067 (CHRONOS-3) only; rifampicin or itraconazole in the dedicated DDI study 16270 only",
    notes          = paste(
      "Baseline demographics are Morcos 2023 Table 1, column 'Pooled PopPK",
      "analyses'; the per-study designs, dosing regimens and PK sampling",
      "schedules are Table S1. The nine studies are 12871, 15205, 16270,",
      "16349 (parts A and B, CHRONOS-1), 16790, 16866, 17067 (CHRONOS-3),",
      "17792 and 18041. Concentrations were measured by validated LC/MS with",
      "an LLOQ of 2 ng/mL; 276 of 5958 observations were below that limit and",
      "were handled with the Beal M3 method, which is a likelihood-estimation",
      "device and has no counterpart in a forward simulation.",
      sep = " "
    )
  )

  ini({
    # --------------------------------------------------------------------
    # Structural parameters. Morcos 2023 Table 2, "Fixed effects (THETA)".
    # These are the typical values for the paper's reference patient: male,
    # in a non-Japanese phase I or phase II study, without rifampicin or
    # itraconazole comedication, and with normal hepatic function
    # (NCI ODWG group 1) -- Table 2 footnotes d and e.
    #
    # NONMEM V1/V2/V3 and Q2/Q3 map onto the canonical nlmixr2lib names as
    # V1 -> vc, V2 -> vp, V3 -> vp2, Q2 -> q, Q3 -> q2.
    # --------------------------------------------------------------------
    lcl  <- log(22.2) ; label("Clearance for the reference patient (L/h)")                    # Table 2: CL_pop = 22.2 L/h (RSE 3.18%, 95% CI 20.8-23.5)
    lvc  <- log(92.1) ; label("Central volume of distribution for the reference patient (L)") # Table 2: V1_pop = 92.1 L (RSE 7.47%, 95% CI 78.6-106)
    lq   <- log(79.3) ; label("Intercompartmental clearance, central to peripheral1 (L/h)")   # Table 2: Q2 = 79.3 L/h (RSE 1.58%, 95% CI 76.8-81.7)
    lvp  <- log(508)  ; label("First peripheral volume of distribution (L)")                  # Table 2: V2 = 508 L (RSE 2.54%, 95% CI 483-534)
    lq2  <- log(7.34) ; label("Intercompartmental clearance, central to peripheral2 (L/h)")   # Table 2: Q3 = 7.34 L/h (RSE 6.96%, 95% CI 6.34-8.34)
    lvp2 <- log(522)  ; label("Second peripheral volume of distribution (L)")                 # Table 2: V3 = 522 L (RSE 4.26%, 95% CI 478-565)

    # --------------------------------------------------------------------
    # Covariate effects. All eight covariates retained by the paper's full
    # stepwise procedure are categorical, so every one uses the Morcos 2023
    # Methods categorical form
    #
    #   P_i = theta_TV * (1 + N1 * theta_cov1 + N2 * theta_cov2 + ...) * exp(eta_i)
    #
    # i.e. each theta is a FRACTIONAL change, not a multiplier and not a
    # log-scale offset. Each of the eight indicators below is its own
    # covariate with a single non-reference level, so each contributes an
    # independent (1 + theta * N) factor. Every value is checked against the
    # percentage the Results section quotes for it.
    # --------------------------------------------------------------------
    e_conmed_rifampicin_cl   <-  1.91   ; label("Fractional change in CL with rifampicin coadministration")   # Table 2: theta_RIFCL = 1.91 (RSE 3.56%); Results "increased CL by 191%"
    e_conmed_itraconazole_cl <- -0.361  ; label("Fractional change in CL with itraconazole coadministration") # Table 2: theta_ITRACL = -0.361 (RSE 5.07%); Results "decreased CL by 36.1%"
    e_study_chronos3_cl      <- -0.184  ; label("Fractional change in CL for CHRONOS-3 (study 17067)")        # Table 2: theta_17067CL = -0.184 (RSE 17.6%); Results "18.4% lower CL"
    e_sexf_cl                <- -0.167  ; label("Fractional change in CL for female sex")                     # Table 2: theta_SEXCL = -0.167 (RSE 17.2%); Results "females had 16.7% lower CL than males"
    e_hepimp_cl              <- -0.192  ; label("Fractional change in CL with any hepatic impairment")        # Table 2: theta_NCICL = -0.192 (RSE 17.4%); Results "19.2% lower CL than individuals with normal hepatic function"
    e_region_japan_cl        <- -0.204  ; label("Fractional change in CL for patients from Japan")            # Table 2: theta_JAPCL = -0.204 (RSE 29.3%); Results "patients from Japan had 20.4% lower CL"

    e_sexf_vc                <- -0.429  ; label("Fractional change in V1 for female sex")                     # Table 2: theta_SEXV1 = -0.429 (RSE 14.6%); Results "females had 42.9% lower V1 than males"
    e_conmed_rifampicin_vc   <-  1.08   ; label("Fractional change in V1 with rifampicin coadministration")   # Table 2: theta_RIFV1 = 1.08 (RSE 39.0%); Results "comedication with rifampin increased V1 by 108%"

    # --------------------------------------------------------------------
    # Inter-individual variability. Morcos 2023 Table 2, "Random effects:
    # inter-individual variability (OMEGA)". Exponential IIV on CL and V1
    # only; the paper reports no IIV on Q2, V2, Q3 or V3 and no OMEGA block,
    # so the two etas are encoded as independent diagonal variances.
    # The values below are the omega^2 variances, which is what nlmixr2's
    # `~` form expects; Table 2's CV% rows are the derived
    # 100*sqrt(exp(omega^2) - 1) transformation of the same numbers.
    # --------------------------------------------------------------------
    etalcl ~ 0.124  # Table 2: omega^2 on CL_pop = 0.124 (RSE 6.07%), i.e. CV 36.3%; shrinkage 22.5%
    etalvc ~ 0.846  # Table 2: omega^2 on V1_pop = 0.846 (RSE 8.94%), i.e. CV 115%;  shrinkage 27.3%

    # --------------------------------------------------------------------
    # Residual error. Morcos 2023 Table 2, "Residual error (SIGMA)", and
    # footnote g: "Both the observations and the model predictions were
    # log-transformed, and an additive residual error model was used. This is
    # equivalent to an exponential residual error model on untransformed
    # data". Log-additive maps onto nlmixr2's `lnorm()` and the canonical
    # `expSd` name, whose SD is on the log scale -- hence sqrt(sigma^2).
    #
    # Three mutually exclusive strata (Discussion: "a residual error model
    # with stratifications related to infusion time (for the first 20 min of
    # the infusion and thereafter)"). The early-infusion stratum is not
    # qualified by study, so it applies to phase I, II and III alike.
    # --------------------------------------------------------------------
    expSd_early   <- 2.2583 ; label("Log-scale residual SD, first 20 min of an infusion")                # Table 2: sigma^2 = 5.10 (RSE 4.87%) -> sqrt(5.10) = 2.2583; CV 1280%
    expSd_phase12 <- 0.4195 ; label("Log-scale residual SD, phase I/II studies beyond 20 min")           # Table 2: sigma^2 = 0.176 (RSE 1.22%) -> sqrt(0.176) = 0.4195; CV 43.9%
    expSd_phase3  <- 0.7950 ; label("Log-scale residual SD, phase III study 17067 beyond 20 min")        # Table 2: sigma^2 = 0.632 (RSE 2.44%) -> sqrt(0.632) = 0.7950; CV 93.8%
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters.
    #
    # Each retained covariate contributes its own multiplicative
    # (1 + theta * indicator) factor, per the Morcos 2023 Methods
    # categorical covariate equation. Setting every indicator to 0
    # recovers the Table 2 reference patient exactly.
    # ------------------------------------------------------------------
    cl <- exp(lcl + etalcl) *
      (1 + e_conmed_rifampicin_cl   * CONMED_RIFAMPICIN)   *
      (1 + e_conmed_itraconazole_cl * CONMED_ITRACONAZOLE) *
      (1 + e_study_chronos3_cl      * STUDY_CHRONOS3)      *
      (1 + e_sexf_cl                * SEXF)                *
      (1 + e_hepimp_cl              * HEPIMP)              *
      (1 + e_region_japan_cl        * REGION_JAPAN)

    vc <- exp(lvc + etalvc) *
      (1 + e_sexf_vc              * SEXF) *
      (1 + e_conmed_rifampicin_vc * CONMED_RIFAMPICIN)

    q   <- exp(lq)
    vp  <- exp(lvp)
    q2  <- exp(lq2)
    vp2 <- exp(lvp2)

    # ------------------------------------------------------------------
    # 2. Micro-constants.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ------------------------------------------------------------------
    # 3. Three-compartment IV disposition with first-order elimination
    #    from the central compartment (Results: "Copanlisib PK were best
    #    described by a linear three-compartment model with first-order
    #    elimination from the central compartment following intravenous
    #    infusion"). Copanlisib is given only as an intravenous infusion,
    #    so there is no depot state and no bioavailability term; dose
    #    records target `central` with a rate or duration.
    # ------------------------------------------------------------------
    d/dt(central)     <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <-   k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-   k13 * central - k31 * peripheral2

    # ------------------------------------------------------------------
    # 4. Observation and stratified log-additive residual error.
    #
    #    Doses are in mg and volumes in L, so Cc is in mg/L. The paper
    #    reports concentrations in ng/mL and exposures in ug*h/L (Figure 3b)
    #    or ng*h/mL (Discussion); 1 mg/L = 1000 ng/mL = 1000 ug/L, so
    #    published values divide by 1000 to compare against Cc.
    #
    #    tad() is rxode2's time since the most recent dose. Because a
    #    copanlisib dose record starts the 1-h infusion, tad() < 20/60 h
    #    selects exactly the paper's "first 20 min of the infusion"
    #    stratum. Beyond 20 min the stratum is chosen by study phase.
    #    This follows the per-record residual-switch precedent in
    #    Bukkems_2021_raltegravir.R and the study-stratified residual
    #    precedent in Xie_2025_aztreonam_avibactam.R.
    # ------------------------------------------------------------------
    Cc <- central / vc

    earlyInfusion <- tad() < (20 / 60)

    expSd <- expSd_early * earlyInfusion +
      (1 - earlyInfusion) * (expSd_phase12 * (1 - STUDY_CHRONOS3) +
                               expSd_phase3 * STUDY_CHRONOS3)

    Cc ~ lnorm(expSd)
  })
}
