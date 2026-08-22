Kong_2025_sudapyridine <- function() {
  description <- "Three-compartment population PK model for the anti-tuberculosis diarylquinoline sudapyridine (WX-081) with a three-transit-compartment absorption chain, jointly fitted with a coupled two-compartment model for its QT-liability metabolite WX-081-M3, in Chinese healthy volunteers and drug-susceptible / multidrug-resistant tuberculosis patients receiving 150-450 mg once-daily oral doses, with a power effect of baseline alkaline phosphatase on the apparent metabolite clearance."
  reference <- paste(
    "Kong W., Liang H., Zhang Y., Li L., Li Y., Yan X., Liu D. (2025).",
    "Population pharmacokinetic and exposure-response study of a novel",
    "anti-tuberculosis drug to inform its dosage design in phase III",
    "clinical trial.",
    "European Journal of Pharmaceutical Sciences 212:107160.",
    "doi:10.1016/j.ejps.2025.107160.",
    sep = " "
  )
  vignette <- "Kong_2025_sudapyridine"

  # Kong 2025 Methods 2.3 states that plasma concentrations "were converted to
  # molar units" before fitting, so the reported clearances and volumes are
  # amount-unit agnostic. Dosing here is in mg of WX-081 and the parent
  # observation Cc is therefore a true WX-081 mass concentration in mg/L
  # (the paper tabulates the same quantity in ug/L; 1 mg/L = 1000 ug/L).
  # See compartmentData for the WX-081-M3 amount convention.
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    ALP = list(
      description        = "Baseline serum alkaline phosphatase activity",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained by Kong 2025 after stepwise forward",
        "inclusion (dOFV > 3.84) and backward elimination (dOFV > 6.63):",
        "baseline ALP on the apparent metabolite clearance CL_M3 (dOFV 8.4;",
        "Kong 2025 Results 3.3). Entered as the paper's continuous power",
        "function (Kong 2025 eqs. 3), P_i = theta1 * (Cov / Cov_median)^theta2,",
        "with theta2 = 0.44, so a subject with above-median ALP has a HIGHER",
        "apparent metabolite clearance. The authors note ALP is unlikely to be",
        "directly involved in WX-081 metabolism and interpret it as an indirect",
        "marker of hepatic clearance (Kong 2025 Discussion).",
        "REFERENCE VALUE: Kong 2025 does not report the pooled-cohort median",
        "ALP. Table 1 reports per-cohort medians of 64.5 U/L (healthy",
        "volunteers, n = 24), 90.0 U/L (drug-susceptible TB, n = 28) and",
        "80.5 U/L (multidrug-resistant TB, n = 20); the subject-count-weighted",
        "mean of those three medians is 78.9 U/L. 80 U/L is used here as the",
        "rounded normalisation constant. The choice is close to immaterial:",
        "because the exponent is only 0.44, a 10% error in the reference shifts",
        "the typical CL_M3 by 4.2%, and setting ALP = 80 reproduces the",
        "published typical value exactly. See vignette Assumptions and",
        "deviations.",
        sep = " "
      ),
      source_name        = "ALP"
    )
  )

  # Covariates Kong 2025 screened but did not retain in the final model
  # (Kong 2025 Methods 2.4 candidate list; Results 3.3 reports that only ALP
  # on CL_M3 survived). Documented here so the paper's covariate screen is
  # preserved without declaring covariates the model never references.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age", units = "years", type = "continuous",
      notes = "Screened as a demographic candidate; not retained (Kong 2025 Methods 2.4, Results 3.3). Cohort medians 28.0 / 44.0 / 36.5 years (Table 1)."
    ),
    WT = list(
      description = "Body weight", units = "kg", type = "continuous",
      notes = "Screened as a demographic candidate; not retained. Kong 2025 Discussion notes that TIME-VARYING weight (a significant covariate for bedaquiline in Svensson 2016) could not be evaluated because sequential weights were not collected."
    ),
    SEXF = list(
      description = "Female sex indicator", units = "(binary)", type = "binary",
      notes = "Screened as a demographic candidate; not retained. Sex on Vc/F was reported for bedaquiline by McLeay 2014 but was not identified here (Kong 2025 Discussion)."
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer versus tuberculosis-patient indicator (the paper's HP covariate)", units = "(binary)", type = "binary",
      notes = "Screened as a candidate; not retained. Kong 2025 Discussion additionally reports that disease state defined as DS-TB versus MDR-TB was negligible (Supplementary Figures S3-S4)."
    ),
    SMOKE_CURRENT = list(
      description = "Current-smoker indicator (the paper's SMOKE covariate)", units = "(binary)", type = "binary",
      notes = "Screened as a candidate; not retained (Kong 2025 Methods 2.4)."
    ),
    ALT = list(
      description = "Alanine aminotransferase", units = "U/L", type = "continuous",
      notes = "Screened as a liver-function candidate; not retained. Cohort medians 12.0 / 18.5 / 10.0 U/L (Kong 2025 Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase", units = "U/L", type = "continuous",
      notes = "Screened as a liver-function candidate; not retained. Cohort medians 17.0 / 16.5 / 15.5 U/L (Kong 2025 Table 1)."
    ),
    TBILI = list(
      description = "Total bilirubin", units = "umol/L", type = "continuous",
      notes = "Screened as a liver-function candidate; not retained. Cohort medians 11.1 / 9.45 / 9.46 umol/L (Kong 2025 Table 1)."
    ),
    TCHOL = list(
      description = "Total serum cholesterol", units = "mmol/L", type = "continuous",
      notes = "Screened as a liver-function candidate; not retained (Kong 2025 Methods 2.4). Baseline values are not tabulated in the paper."
    )
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  #
  # AMOUNT CONVENTION FOR THE METABOLITE STATES. Kong 2025 fitted the joint
  # model in molar units and the metabolite formation flux is FM * CL / Vc
  # (Kong 2025 Fig. 1), so with dosing in mg of WX-081 the two WX-081-M3
  # states carry WX-081 MOLAR EQUIVALENTS expressed on the parent's mass
  # scale. Cc_m3 is therefore an M3 concentration in parent-molar-equivalent
  # mg/L; the true WX-081-M3 mass concentration is Cc_m3 * MW_M3 / MW_WX-081.
  # Neither molecular weight is reported anywhere in Kong 2025 or its
  # supplement, so the conversion factor cannot be sourced and is deliberately
  # NOT applied here. Every unit-free metabolite quantity (the M3-to-parent
  # molar AUC ratio, and ratios of M3 Cmax between dosing regimens) is
  # reproduced exactly; see the vignette.
  compartmentData <- list(
    depot           = list(analyte = "sudapyridine (WX-081)",   units = "mg", specimen = "administration site", verified = TRUE),
    transit1        = list(analyte = "sudapyridine (WX-081)",   units = "mg", specimen = "administration site", verified = TRUE),
    transit2        = list(analyte = "sudapyridine (WX-081)",   units = "mg", specimen = "administration site", verified = TRUE),
    transit3        = list(analyte = "sudapyridine (WX-081)",   units = "mg", specimen = "administration site", verified = TRUE),
    central         = list(analyte = "sudapyridine (WX-081)",   units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1     = list(analyte = "sudapyridine (WX-081)",   units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2     = list(analyte = "sudapyridine (WX-081)",   units = "mg", specimen = "plasma", verified = TRUE),
    central_m3      = list(analyte = "WX-081-M3 metabolite",    units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_m3  = list(analyte = "WX-081-M3 metabolite",    units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 72L,
    n_studies      = 2L,
    age_range      = "20-60 years",
    age_median     = "28.0 years (healthy volunteers), 44.0 years (DS-TB), 36.5 years (MDR-TB)",
    weight_range   = "41.0-77.0 kg",
    weight_median  = "62.6 kg (healthy volunteers), 58.0 kg (DS-TB), 55.1 kg (MDR-TB)",
    sex_female_pct = 27.8,
    race_ethnicity = c(Chinese = 100),
    disease_state  = paste(
      "24 healthy volunteers (study A, phase I, NCT06117514), 28 patients with",
      "drug-susceptible pulmonary tuberculosis and 20 patients with",
      "multidrug-resistant tuberculosis (study B, phase II, NCT04608955).",
      "MDR-TB patients received a multi-drug background treatment of",
      "levofloxacin, isoniazid, cycloserine and pyrazinamide during the second",
      "treatment stage.",
      sep = " "
    ),
    dose_range     = paste(
      "Study A: single oral 200, 300 or 400 mg on D1 then once daily D4-D14.",
      "Study B: 150, 300 or 450 mg once daily for 14 days (DS-TB); 400 mg once",
      "daily for 14 days followed by 150 mg once daily for 6 weeks (MDR-TB).",
      sep = " "
    ),
    regions        = "China",
    notes          = paste(
      "Baseline demographics from Kong 2025 Table 1. 1610 WX-081 and 1580",
      "WX-081-M3 plasma concentrations entered the analysis; 0.86% of WX-081",
      "and 2.71% of WX-081-M3 observations were below the limit of",
      "quantification (LLOQ 10 ug/L for WX-081 and 0.1 ug/L for WX-081-M3) and",
      "were excluded. Race was not evaluated as a covariate because both",
      "studies enrolled only Chinese participants.",
      sep = " "
    )
  )

  ini({
    # Final parameter estimates from Kong 2025 Table 2 "Estimated PPK final
    # model parameters with 500 bootstraps"; the parenthesised pairs in that
    # table are the 500-replicate bootstrap intervals and the RSE column is
    # from the NONMEM covariance step.
    #
    # APPARENT-PARAMETER CONVENTION. Bioavailability F was not identifiable, so
    # every WX-081 parameter is an apparent (/F) value. The metabolite
    # parameters are likewise apparent, and their scaling follows Kong 2025
    # Fig. 1, which is the authoritative statement of the fitted structure:
    # the metabolite central compartment A8 is drawn carrying the PARENT
    # central volume V_C, it is fed at rate F_M * CL / V_c * A5, and it is
    # emptied at rate CL_M / V_C * A8. That is exactly the identifiability
    # device described in Methods 2.3 -- "either FM or Vc,M need to be fixed
    # ... we adopted the assumption that the volume of distribution of
    # WX-081-M3 is identical with that of WX-081" -- which fixes Vc,M and
    # thereby makes FM estimable. The metabolite rows of Table 2 are
    # consequently labelled "/(F*FM)" but are numerically CL_M/F, Q_M/F and
    # V9/F on a parent-molar-equivalent amount scale (Results 3.3 writes the
    # same parameter as "CLM/F"). Taking the "/(F*FM)" labels literally instead
    # would place the metabolite formation flux at the full CL * Cc and give a
    # steady-state M3-to-parent molar AUC ratio of CL/CL_M = 3.96/0.783 = 5.06,
    # 63-fold above the 0.08 the paper reports in its Discussion; the encoding
    # used here gives FM * CL / CL_M = 0.0235 * 3.96 / 0.783 = 0.119 at steady
    # state, consistent with an observed day-14 ratio of 0.08 for a metabolite
    # that is still accumulating. See vignette Assumptions and deviations.

    # Absorption: dose enters the depot and passes through three transit
    # compartments at rate ktr, with the final transit compartment absorbing
    # into central at rate ka (Kong 2025 Fig. 1 and Results 3.2).
    lktr    <- log(0.708)  ; label("Transit-chain rate constant ktr (1/h, on log scale)")                                   # Kong 2025 Table 2 'ktr = 0.708 1/hr' (RSE 3.5%; bootstrap 0.647-0.789)
    lka     <- log(1.49)   ; label("First-order absorption rate constant ka from the last transit compartment (1/h, on log scale)") # Kong 2025 Table 2 'ka = 1.49 1/hr' (RSE 16.1%; bootstrap 1.16-2.00)

    # WX-081 disposition: three compartments.
    lcl     <- log(3.96)   ; label("Apparent WX-081 clearance CL/F (L/h, on log scale)")                                    # Kong 2025 Table 2 'CL/F = 3.96 L/hr' (RSE 7.7%; bootstrap 3.11-4.68)
    lvc     <- log(13.1)   ; label("Apparent WX-081 central volume Vc/F (L, on log scale)")                                 # Kong 2025 Table 2 'Vc/F = 13.1 L' (RSE 13.7%; bootstrap 9.64-18.8)
    lq      <- log(14.8)   ; label("Apparent WX-081 first inter-compartmental clearance Q1/F (L/h, on log scale)")          # Kong 2025 Table 2 'Q1/F = 14.8 L/hr' (RSE 6.6%; bootstrap 13.2-16.4)
    lvp     <- log(140)    ; label("Apparent WX-081 first peripheral volume V6/F (L, on log scale)")                        # Kong 2025 Table 2 'V6/F = 140 L' (RSE 4.7%; bootstrap 111-161)
    lq2     <- log(4.66)   ; label("Apparent WX-081 second inter-compartmental clearance Q2/F (L/h, on log scale)")         # Kong 2025 Table 2 'Q2/F = 4.66 L/hr' (RSE 13.3%; bootstrap 4.04-5.59)
    lvp2    <- log(1040)   ; label("Apparent WX-081 second peripheral volume V7/F (L, on log scale)")                       # Kong 2025 Table 2 'V7/F = 1040 L' (RSE 16.5%; bootstrap 598-1948)

    # WX-081-M3 formation and disposition: two compartments sharing the parent
    # central volume (Kong 2025 Fig. 1, compartment A8 labelled 'V_C').
    lfm     <- log(0.0235) ; label("Fraction of WX-081 clearance forming WX-081-M3, FM (unitless, on log scale)")           # Kong 2025 Table 2 'FM = 2.35 %' (RSE 9%; bootstrap 1.62-4.02%)
    lcl_m3  <- log(0.783)  ; label("Apparent WX-081-M3 clearance CL_M3/F at the reference ALP (L/h, on log scale)")         # Kong 2025 Table 2 'CLM/(F*FM) = 0.783 L/hr' (RSE 9.1%; bootstrap 0.529-1.31)
    lq_m3   <- log(23.4)   ; label("Apparent WX-081-M3 inter-compartmental clearance Q_M3/F (L/h, on log scale)")           # Kong 2025 Table 2 'QM/(F*FM) = 23.4 L/hr' (RSE 12.8%; bootstrap 11.6-50.5)
    lvp_m3  <- log(92.9)   ; label("Apparent WX-081-M3 peripheral volume V9/F (L, on log scale)")                           # Kong 2025 Table 2 'V9/(F*FM) = 92.9 L' (RSE 7.7%; bootstrap 61.9-173)

    # Covariate effect: power function on the apparent metabolite clearance,
    # normalised to the population median ALP (Kong 2025 eqs. 3).
    e_alp_cl_m3 <- 0.44    ; label("Power exponent for baseline ALP on apparent WX-081-M3 clearance (unitless)")            # Kong 2025 Table 2 'ALP on CLM = 0.44' (RSE 33.2%; bootstrap 0.225-0.755)

    # Between-subject variability, exponential (Kong 2025 eqs. 1). Table 2
    # reports the "Omega / IIV" column as the eta STANDARD DEVIATION, not the
    # variance: with 72 subjects the asymptotic RSE of a variance component
    # cannot fall below sqrt(2/72) = 16.7%, yet every reported omega RSE is
    # between 8.6% and 13.8%, and the bootstrap intervals imply the same
    # (e.g. CL/F: (0.458-0.314)/(2*1.96*0.40) = 9.2%). On the SD scale the
    # corresponding bound is sqrt(1/(2*72)) = 8.3%, which every reported value
    # clears. nlmixr2's ini() takes VARIANCES, so each entry below is the
    # published omega squared. See vignette Assumptions and deviations.
    etalcl    ~ 0.160000   # omega = 0.40  -> var = 0.40^2  ; Kong 2025 Table 2 IIV CL/F = 0.40 (RSE 11%, shrinkage 5.6%)
    etalvc    ~ 1.144900   # omega = 1.07  -> var = 1.07^2  ; Kong 2025 Table 2 IIV Vc/F = 1.07 (RSE 10.9%, shrinkage 31.6%)
    etalq     ~ 0.148996   # omega = 0.386 -> var = 0.386^2 ; Kong 2025 Table 2 IIV Q1/F = 0.386 (RSE 12.4%, shrinkage 18.4%)
    etalq2    ~ 0.885481   # omega = 0.941 -> var = 0.941^2 ; Kong 2025 Table 2 IIV Q2/F = 0.941 (RSE 10.7%, shrinkage 7.4%)
    etalka    ~ 0.992016   # omega = 0.996 -> var = 0.996^2 ; Kong 2025 Table 2 IIV ka = 0.996 (RSE 13.8%, shrinkage 27.8%)
    etalcl_m3 ~ 0.159201   # omega = 0.399 -> var = 0.399^2 ; Kong 2025 Table 2 IIV CLM = 0.399 (RSE 10.9%, shrinkage 8.5%)
    etalfm    ~ 0.216225   # omega = 0.465 -> var = 0.465^2 ; Kong 2025 Table 2 IIV FM = 0.465 (RSE 8.6%, shrinkage 0.5%)
    # No IIV was reported on ktr, V6/F, V7/F, Q_M3 or V9/F (Kong 2025 Table 2
    # records 'NA' in the Omega columns for those rows).

    # Residual error. Kong 2025 Methods 2.3 eqs. 2 specifies a PROPORTIONAL
    # error model for both analytes, but Table 2 does not report the sigma
    # estimates and no supplement, control stream or figure in the on-disk
    # sources carries them (the RSE sentence in Results 3.4 covers only "the
    # fixed and random effects"). Rather than invent a residual magnitude the
    # proportional error is encoded with a magnitude of zero, so simulations
    # from this model are IPRED-only. Users refitting or simulating with
    # residual error should re-estimate propSd / propSd_m3 or supply their own.
    propSd    <- fixed(0)  ; label("WX-081 proportional residual error (fraction; magnitude not reported)")     # Kong 2025 eqs. 2 (proportional form stated; sigma not reported in Table 2)
    propSd_m3 <- fixed(0)  ; label("WX-081-M3 proportional residual error (fraction; magnitude not reported)")  # Kong 2025 eqs. 2 (proportional form stated; sigma not reported in Table 2)
  })

  model({
    # 1. Individual parameters. Between-subject variability is exponential
    #    (Kong 2025 eqs. 1: P_i = TVP * exp(eta_i)).
    ktr   <- exp(lktr)
    ka    <- exp(lka  + etalka)
    cl    <- exp(lcl  + etalcl)
    vc    <- exp(lvc  + etalvc)
    q     <- exp(lq   + etalq)
    vp    <- exp(lvp)
    q2    <- exp(lq2  + etalq2)
    vp2   <- exp(lvp2)
    fm    <- exp(lfm  + etalfm)

    # 2. Baseline ALP power covariate on the apparent metabolite clearance
    #    (Kong 2025 eqs. 3, Table 2 'ALP on CLM'). 80 U/L is the rounded
    #    population-median normalisation constant; Kong 2025 does not tabulate
    #    the pooled median, and the subject-count-weighted mean of the three
    #    per-cohort medians in Table 1 is 78.9 U/L. See covariateData$ALP.
    cl_m3 <- exp(lcl_m3 + etalcl_m3) * (ALP / 80)^e_alp_cl_m3
    q_m3  <- exp(lq_m3)
    vp_m3 <- exp(lvp_m3)

    # 3. The WX-081-M3 central compartment carries the PARENT central volume:
    #    Kong 2025 Fig. 1 draws compartment A8 as 'A8, V_C' and its outflows as
    #    CL_M/V_C and Q_M, which is the Methods 2.3 identifiability assumption
    #    "the volume of distribution of WX-081-M3 is identical with that of
    #    WX-081" made explicit. IIV on Vc therefore propagates to the
    #    metabolite central compartment, as in the published model.
    vc_m3 <- vc

    # 4. ODE system. Depot -> three transit compartments (rate ktr) -> central
    #    (rate ka); three-compartment WX-081 disposition; metabolite formed at
    #    FM * CL * Cc and disposed of in two compartments (Kong 2025 Fig. 1).
    d/dt(depot)          <- -ktr * depot
    d/dt(transit1)       <-  ktr * depot    - ktr * transit1
    d/dt(transit2)       <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)       <-  ktr * transit2 - ka  * transit3
    d/dt(central)        <-  ka  * transit3 -
                             cl  * central / vc -
                             q   * (central / vc - peripheral1 / vp) -
                             q2  * (central / vc - peripheral2 / vp2)
    d/dt(peripheral1)    <-  q   * (central / vc - peripheral1 / vp)
    d/dt(peripheral2)    <-  q2  * (central / vc - peripheral2 / vp2)
    d/dt(central_m3)     <-  fm  * cl * central / vc -
                             cl_m3 * central_m3 / vc_m3 -
                             q_m3  * (central_m3 / vc_m3 - peripheral1_m3 / vp_m3)
    d/dt(peripheral1_m3) <-  q_m3  * (central_m3 / vc_m3 - peripheral1_m3 / vp_m3)

    # 5. Observations. Cc is a true WX-081 plasma concentration; Cc_m3 is the
    #    WX-081-M3 concentration expressed in WX-081 molar equivalents (see
    #    compartmentData for why the molecular-weight conversion is not
    #    applied).
    Cc    <- central    / vc
    Cc_m3 <- central_m3 / vc_m3

    Cc    ~ prop(propSd)
    Cc_m3 ~ prop(propSd_m3)
  })
}
