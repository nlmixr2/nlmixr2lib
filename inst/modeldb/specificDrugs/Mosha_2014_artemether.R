# Joint parent-metabolite population PK model of oral artemether and its
# active metabolite dihydroartemisinin (DHA) in pregnant and non-pregnant
# Tanzanian women with uncomplicated Plasmodium falciparum malaria
# (Mosha 2014, Antimicrobial Agents and Chemotherapy 58(8):4583-4592;
# doi:10.1128/AAC.02595-14).

Mosha_2014_artemether <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral artemether (AM)",
    "and dihydroartemisinin (DHA) in 33 pregnant (2nd or 3rd trimester) and",
    "22 non-pregnant women with uncomplicated Plasmodium falciparum malaria",
    "in Rufiji, Tanzania after standard fixed-dose artemether-lumefantrine",
    "(Mosha 2014). One-compartment AM disposition with first-order",
    "absorption and linear metabolism to a one-compartment DHA disposition,",
    "including a presystemic AM-to-DHA conversion fraction (1 - F1) with",
    "F1 = expit(logit_F1) parameterised on the logit scale. Absorption ka",
    "fixed at 0.70 1/h and DHA volume fixed equal to AM Vc per the source",
    "paper. IIV is present only on AM CL (99% CV); the remaining structural",
    "parameters carry no IIV in the final model. None of the available",
    "covariates (pregnancy, body weight, BMI, age, gestational age,",
    "diarrhoea) reached statistical significance on AM or DHA PK in",
    "Mosha 2014 and so none are encoded. AM residual error is combined",
    "(proportional plus additive); DHA residual error is proportional.",
    sep = " "
  )
  reference <- paste(
    "Mosha D, Guidi M, Mwingira F, Abdulla S, Mercier T, Decosterd LA,",
    "Csajka C, Genton B (2014).",
    "Population pharmacokinetics and clinical response for",
    "artemether-lumefantrine in pregnant and nonpregnant women with",
    "uncomplicated Plasmodium falciparum malaria in Tanzania.",
    "Antimicrobial Agents and Chemotherapy 58(8):4583-4592.",
    "doi:10.1128/AAC.02595-14.",
    sep = " "
  )
  vignette <- "Mosha_2014_artemether_lumefantrine"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  population <- list(
    species         = "human",
    n_subjects      = 55L,
    n_studies       = 1L,
    n_pregnant      = 33L,
    n_nonpregnant   = 22L,
    n_observations_am  = 146L,
    n_observations_dha = 98L,
    age_range       = "18-41 years (pregnant median 25; non-pregnant median 21.5; Table 1)",
    weight_range    = "40-80 kg (pregnant median 52; non-pregnant median 48.5; Table 1)",
    sex_female_pct  = 100,
    disease_state   = paste(
      "Uncomplicated Plasmodium falciparum malaria; pregnant women in",
      "the second or third trimester (median gestational age 27 weeks,",
      "range 14-37) compared concurrently with non-pregnant women",
      "(controls) recruited from the same clinic during the same period."
    ),
    dose_range      = paste(
      "Coartem (Novartis): 20 mg artemether + 120 mg lumefantrine per",
      "tablet; four tablets per dose (80 mg artemether) given twice daily",
      "for 3 days (oral; dose times 0, 8, 24, 36, 48, 60 hours) with",
      "200 mL milk containing 4.5 g fat to optimise lumefantrine",
      "bioavailability."
    ),
    regions         = "Tanzania (Rufiji district, Kibiti health center; April-September 2012)",
    notes           = paste(
      "Demographics from Mosha 2014 Table 1. 25% (n = 37) of AM and 7%",
      "(n = 7) of DHA concentrations were below the limit of quantification",
      "and handled via the M3 method (Beal). Companion lumefantrine",
      "extraction from the same cohort:",
      "modellib('Mosha_2014_lumefantrine'). The desbutyl-lumefantrine",
      "metabolite arm of the published joint model is not encoded; see",
      "the validation vignette for the rationale."
    )
  )

  ini({
    # Structural population mean parameters from Mosha 2014 Table 3
    # "Artemether" section, "Population pharmacokinetics analysis" column.
    # The paper reports CL, V_c, CL_met on the linear scale (L/h, L) and
    # K_23 on a rate-constant scale (1/h); log() is applied here for the
    # nlmixr2 internal scale. The reported logit-F1 typical of 1.4
    # corresponds to expit(1.4) = 0.802, i.e. 80.2% of the AM dose enters
    # AM central directly and the remaining 19.8% is presystemically
    # converted to DHA -- consistent with the paper's Results sentence
    # "Our results show that 21% of the AM dose is converted presystemically
    # into DHA".
    lcl     <- log(98)                ; label("Apparent artemether elimination clearance, CL/F (L/h)")               # Mosha 2014 Table 3: CL = 98 (RSE 24%; bootstrap 102, 95% CI 69-140)
    lvc     <- log(373)               ; label("Apparent artemether central volume of distribution, Vc/F (L)")        # Mosha 2014 Table 3: Vc = 373 (RSE 16%; bootstrap 354, 95% CI 225-492)
    logitfdepot <- 1.4                ; label("Logit-transformed fraction of AM dose entering AM central, logit(F1); 1 - F1 = presystemic AM-to-DHA fraction")  # Mosha 2014 Table 3: Logit F1 = 1.4 (RSE 27%; bootstrap 1.5, 95% CI 0.7 to 2.6); F1 = expit(1.4) = 0.802
    lka     <- fixed(log(0.70))       ; label("Absorption rate constant, ka (1/h); fixed")                            # Mosha 2014 Table 3: Ka fixed at 0.70 1/h (Methods: "the absorption rate constants ... were thus fixed ... to achieve peak plasma AM and LF concentrations, respectively, 2 h and 6 to 8 h after drug intake")
    lk23    <- log(0.084)             ; label("AM-to-DHA formation rate constant, K23 (1/h)")                         # Mosha 2014 Table 3: K23 = 0.084 (bootstrap 0.088, 95% CI 0.05-0.16)
    lcl_dha <- log(71)                ; label("Apparent dihydroartemisinin elimination clearance, CL_met/F (L/h)")    # Mosha 2014 Table 3: CL_met = 71 (RSE 46%; bootstrap 69, 95% CI 38-136)

    # DHA central volume Vm is fixed equal to AM Vc in the source paper
    # ("Owing to identification problems, the volume of distributions of
    # DLF and DHA could not be estimated and were assumed to be equal to
    # those of LF and AM, respectively"). Encoded inside model() by
    # setting vc_dha <- vc; no separate ini() parameter is declared.

    # IIV. Mosha 2014 reports IIV as %CV (Table 3 column "IIV (%)").
    # The internal variances are recovered by omega^2 = log((CV/100)^2 + 1):
    #   CL_AM CV 99% -> log(0.99^2 + 1) = 0.68273
    # The paper retained IIV only on AM CL ("Inclusion of interpatient
    # variability of Vc, CL_met, K23, or F1 in addition to AM CL did not
    # improve description of the data") so a single eta is carried.
    etalcl  ~ 0.68273  # Mosha 2014 Table 3: BSV on AM CL/F = 99% CV (RSE 65%; bootstrap 93%, 95% CI 66-120)

    # Residual error. AM observation has a combined proportional + additive
    # error model on the linear-concentration scale (Results: "A mixed-error
    # model best described residual intrapatient variability for AM").
    # Table 3 reports the proportional component on AM as 72% CV and the
    # additive component as 0.13 umol/L; converting the additive component
    # to ng/mL using the AM molecular weight 298.4 g/mol gives
    # 0.13 umol/L * 298.4 g/mol = 38.792 ng/mL. DHA observation has a pure
    # proportional error model with 53% CV.
    propSd     <- 0.72    ; label("Proportional residual SD for AM plasma concentration (fraction)")    # Mosha 2014 Table 3: sigma_prop,AM = 72% CV (RSE 26%; bootstrap 69%, 95% CI 49-87)
    addSd      <- 38.792  ; label("Additive residual SD for AM plasma concentration (ng/mL)")           # Mosha 2014 Table 3: sigma_add,AM = 0.13 umol/L * 298.4 g/mol = 38.79 ng/mL (RSE 7%; bootstrap 0.13, 95% CI 0.03-0.20 umol/L)
    propSd_dha <- 0.53    ; label("Proportional residual SD for DHA plasma concentration (fraction)")   # Mosha 2014 Table 3: sigma_prop,DHA = 53% CV (RSE 14%; bootstrap 51%, 95% CI 44-59)
  })

  model({
    # Molecular weights (g/mol). Source paper Methods (Drug assay) does not
    # state MW values explicitly; canonical values from the IUPAC compound
    # database are used (AM C16H26O5 = 298.37 g/mol, DHA C15H24O5 = 284.35
    # g/mol). These are used at the AM-to-DHA presystemic split and at the
    # systemic AM-to-DHA formation step under the assumption of 1:1 molar
    # conversion.
    mw_arm <- 298.4
    mw_dha <- 284.36

    # Individual PK parameters. Only AM CL carries IIV; ka is fixed; Vc,
    # K23, CL_met, and logit(F1) are point estimates with no eta.
    ka      <- exp(lka)
    cl      <- exp(lcl + etalcl)
    vc      <- exp(lvc)
    vc_dha  <- vc                          # Vm = Vc, fixed structurally
    cl_dha  <- exp(lcl_dha)
    k23     <- exp(lk23)

    # Fraction of AM dose entering AM central via the logit-transformed F1.
    # f1 is unitless on (0, 1); (1 - f1) of the dose bypasses AM central
    # and feeds directly into central_dha with AM-to-DHA molar conversion.
    f1      <- expit(logitfdepot)

    # Elimination rate constants. AM is eliminated at kel = cl / vc
    # (apparent total clearance), and DHA is eliminated at kel_dha =
    # cl_dha / vc_dha. K23 is the formation-rate constant feeding DHA from
    # AM central; it is reported alongside CL as a separate parameter in
    # Mosha 2014 Table 3 and is interpreted here as the AM-to-DHA flux
    # constant rather than as the AM elimination rate (which is fully
    # captured by cl / vc).
    kel     <- cl     / vc
    kel_dha <- cl_dha / vc_dha

    # ODE system: first-order absorption from depot with F1-fraction
    # routing into AM central (mass units mg) and (1 - F1)-fraction
    # routing into DHA central (mass units mg, converted from AM mass via
    # the MW ratio). Systemic AM is eliminated at kel; DHA is fed by both
    # the presystemic-conversion arm and the systemic K23 arm and is
    # eliminated at kel_dha.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  f1       * ka * depot                            - kel * central
    d/dt(central_dha) <-  (1 - f1) * ka * depot * (mw_dha / mw_arm) +
                          k23       * central   * (mw_dha / mw_arm) -
                          kel_dha   * central_dha

    # Plasma concentrations in ng/mL. With doses in mg and volumes in L,
    # central / vc has units mg/L = ug/mL; multiplying by 1000 yields
    # ng/mL, the unit the paper reports (e.g. Mosha 2014 Figure 2A).
    Cc     <- 1000 * central     / vc
    Cc_dha <- 1000 * central_dha / vc_dha

    # Residual error. AM observation: combined proportional + additive on
    # linear-concentration scale (Results, "A mixed-error model best
    # described residual intrapatient variability for AM"). DHA
    # observation: pure proportional ("A ... proportional one for DHA").
    Cc     ~ add(addSd) + prop(propSd)
    Cc_dha ~ prop(propSd_dha)
  })
}
