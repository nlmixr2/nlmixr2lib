ParraGuillen_2025_v937 <- function() {
  description <- paste(
    "Semi-mechanistic population PK / anti-drug antibody (PK-ADA) model for",
    "the oncolytic virus V937 (bioselected Coxsackievirus A21) in adults with",
    "advanced solid tumours given V937 intravenously or intratumorally, alone",
    "or with pembrolizumab. Three-compartment viral kinetics (serum, a",
    "theoretical tumour compartment and other tissues sharing one",
    "intercompartmental clearance; a logit-scale fraction FIT of an",
    "intratumoral dose reaches the tumour compartment). Serum virus drives",
    "B-cell activation through a saturable Emax term; three B-cell states in",
    "series (first-order transit KTR) delay the appearance of circulating",
    "neutralising ADA, which is degraded first-order (KDEG). A subject's",
    "observed pre-existing ADA titre is added to the predicted ADA when it is",
    "above the 1:16 detection limit. ADA has no effect on viral clearance."
  )
  reference <- paste(
    "Parra-Guillen ZP, Troconiz IF, Freshwater T. Role of Antidrug Antibodies",
    "in Oncolytic Viral Therapy: A Dynamic Modelling Approach in Cancer",
    "Patients Treated with V937 Alone or in Combination. Clin Pharmacokinet.",
    "2025;64(10):1549-1559. doi:10.1007/s40262-025-01546-9"
  )
  vignette <- "ParraGuillen_2025_v937"

  # The theoretical tumour compartment is carried as `peripheral1` and the
  # other-tissue compartment as `peripheral2` (Fig. 1 TUM / PER). The B-cell
  # transit chain and the circulating-ADA pool are not drug-disposition
  # compartments; Fig. 1 notes that Bcel1-3 and ADA share the ADA unit
  # (dilution factor).
  paper_specific_compartments <- c("bcell1", "bcell2", "bcell3", "ada", "ada_total")

  units <- list(
    time = "h",
    dosing = "copies",
    concentration = "copies/mL"
  )
  dosing <- c("peripheral1", "central")

  compartmentData <- list(
    central = list(analyte = "V937", units = "copies", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "V937", units = "copies", specimen = "tumor", verified = TRUE),
    peripheral2 = list(analyte = "V937", units = "copies", specimen = "tissue", verified = TRUE),
    bcell1 = list(
      analyte = "anti-V937 neutralising antibody precursor (activated B cells)",
      units = "dilution factor",
      specimen = "not applicable",
      verified = TRUE
    ),
    bcell2 = list(
      analyte = "anti-V937 neutralising antibody precursor (activated B cells)",
      units = "dilution factor",
      specimen = "not applicable",
      verified = TRUE
    ),
    bcell3 = list(
      analyte = "anti-V937 neutralising antibody precursor (activated B cells)",
      units = "dilution factor",
      specimen = "not applicable",
      verified = TRUE
    ),
    ada = list(
      analyte = "anti-V937 neutralising antibody",
      units = "dilution factor",
      specimen = "serum",
      verified = TRUE
    )
  )

  covariateData <- list(
    ADA_TITER_BASE = list(
      description = paste(
        "Observed pre-existing (pre-treatment) neutralising anti-V937 antibody",
        "titre, as the reciprocal dilution (dilution factor). One value per",
        "subject, measured before the first V937 dose."
      ),
      units = "dilution factor (reciprocal titre, e.g. 16, 32, 64)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters only the ADA observation: when ADA_TITER_BASE >= 16 (the 1:16",
        "assay detection limit) the observed baseline is added to the",
        "predicted ADA (ESM 2 $PK 'IF(NAB0>=16) BASENAB = NAB0'; $ERROR",
        "'IPRED = LOG(NAB + 0.001 + BASENAB)'); values below 16 contribute",
        "nothing. Set ADA_TITER_BASE = 0 for ADA-naive subjects. About 10% of",
        "the 208 patients (n = 20) had pre-existing ADA (Results, PK-ADA",
        "Model). The term is a fixed per-subject offset, not a dynamic state:",
        "the authors found a mixture model unstable and fixing ADA0 to the",
        "observed value gave the lowest OFV (ESM 7, model 20)."
      ),
      source_name = "NAB0"
    )
  )

  covariatesDataExcluded <- list(
    ROUTE_IT = list(
      description = "Intratumoral (1) versus intravenous (0) V937 administration",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (intravenous)",
      notes = paste(
        "Screened univariately on KSYN, KTR and KDEG (ESM 3 Table 1, models",
        "5-7; dOFV -4.34 to -0.19) and not retained. Route enters the model",
        "only through the dosing compartment (central vs peripheral1) and the",
        "intratumoral bioavailability FIT."
      ),
      source_name = "RT"
    ),
    CONMED_PEMBROLIZUMAB = list(
      description = "V937 given with pembrolizumab (1) versus V937 monotherapy (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (V937 monotherapy)",
      notes = paste(
        "Screened univariately on KSYN, KTR and KDEG (ESM 3 Table 1, models",
        "2-4). TRT on KDEG was imprecise (RSE > 100%); TRT on KTR inflated the",
        "IIV of KSYN and KDEG. Not retained. The PK model built on",
        "monotherapy data predicted the combination data without adjustment",
        "(Results, PK Model)."
      ),
      source_name = "TRT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 208L,
    n_studies = 4L,
    age_range = "28-94 years (mean 65.3, SD 12.5)",
    weight_range = "45.7-164 kg (mean 84.6, SD 21.3)",
    sex_female_pct = 30.3,
    race_ethnicity = c(White = 92.8, Asian = 4.3, Hawaiian = 0.5, Black = 0.5, Other = 1.9),
    disease_state = paste(
      "Advanced solid tumours: melanoma 52.4%, non-small cell lung cancer",
      "24.5%, bladder (urothelial) cancer 21.2%, metastatic castration-",
      "resistant prostate cancer 1.9%. Baseline ECOG 0 in 55.8%, 1 in 40.4%."
    ),
    dose_range = paste(
      "Intratumoral V937 up to 3 x 10^8 TCID50 (hyper-fractionated over",
      "injected lesions) on days 1, 3, 5 and 8 of cycle 1 then day 1 of each",
      "21-day cycle (V937-006 monotherapy, n = 70; V937-007 with",
      "pembrolizumab 2 mg/kg, n = 35). Intravenous V937 1 x 10^8, 3 x 10^8 or",
      "1 x 10^9 TCID50 as a 30-min infusion on days 1, 3, 5 (and 8 in Part B)",
      "of cycle 1 then day 1 of later cycles (MK3475-200 Part A monotherapy,",
      "n = 18; Part B with pembrolizumab 200 mg, n = 85)."
    ),
    regions = "Multinational (four Merck Sharp & Dohme sponsored trials).",
    notes = paste(
      "Pooled V937-006 (NCT01227551) plus extension (NCT01636882),",
      "V937-007 (NCT02565992) and MK3475-200 Parts A/B (NCT02043665).",
      "3628 serum V937 RNA measurements (1613 below the 3000 copies/mL",
      "detection limit) and 2163 neutralising ADA titres (704 below the 1:16",
      "limit). Demographics from Table 1; study design from ESM 5 Table S1;",
      "measurement counts from ESM 6 Table S2. NONMEM 7.4, FOCE-I with the",
      "Laplacian method and M3 censoring."
    )
  )

  ini({
    # PK model of V937 -- Table 2. All PK values were fixed when the ADA
    # model was estimated (ESM 2 $THETA lines 2-8, each 'FIX' with the same
    # value), so ESM 2 reproduces the Table 2 PK block exactly.
    lvc <- log(87.9); label("Apparent volume of the central serum compartment VS (L)") # Table 2 VS = 87.9 L (RSE 31.2%); ESM 2 THETA(2) 87.9
    lcl <- log(111); label("Viral clearance from serum CL (L/h)") # Table 2 CL = 111 L/h (RSE 15.5%); ESM 2 THETA(3) 111
    lq <- log(35.3); label("Intercompartmental clearance QTUM = QPER, shared by both peripheral compartments (L/h)") # Table 2 QTUM = QPER = 35.3 L/h (RSE 33.9%); ESM 2 THETA(4) 35.3 and Q3 = Q
    lvp <- log(103); label("Apparent volume of the theoretical tumour compartment VTUM (L)") # Table 2 VTUM = 103 L (RSE 26.9%); ESM 2 THETA(5) 103
    lvp2 <- log(893); label("Apparent volume of the other-tissue compartment VPER (L)") # Table 2 VPER = 893 L (RSE 48.2%); ESM 2 THETA(7) 893
    logitfdepot <- logit(0.201); label("Fraction FIT of an intratumoral dose available to the tumour compartment (unitless, logit scale)") # Table 2 FIT = 0.201 (RSE 32.0%); ESM 2 THETA(8) 0.201 with LGF2 = LOG(TVF2/(1-TVF2))

    # ADA model -- Table 2.
    lksyn <- log(113); label("First-order B-cell synthesis rate constant KSYN (dilution factor/h)") # Table 2 KSYN = 113 DF/h (RSE 10.3%)
    lec50 <- log(0.754 * 1e8); label("Serum V937 concentration giving half-maximal B-cell synthesis EC50 (copies/L)") # Table 2 EC50 = 0.754 x 10^8 copies/L (RSE 13.8%); ESM 2 TVEC50 = THETA(12)*1e8 against CP = A(1)/V1 in copies/L
    hill <- fixed(1); label("Hill coefficient of the B-cell activation Emax term (unitless)") # ESM 2 THETA(13) GA = 1 FIX; ESM 7 model 17 (sigmoid) rejected
    lktr <- log(5.14e-3); label("First-order transit rate constant between B-cell states KTR (1/h)") # Table 2 KTR = 5.14 x 10^-3 1/h (RSE 7.98%)
    lkdeg <- log(0.250e-3); label("First-order degradation rate constant of circulating ADA KDEG (1/h)") # Table 2 KDEG = 0.250 x 10^-3 1/h (RSE 12.8%)

    # IIV. The paper reports IIV as CV% = sqrt(exp(omega^2) - 1) (Table 2
    # footnote b), which reproduces the ESM 2 value fixed for CL exactly:
    # log(1 + 0.675^2) = 0.375 = $OMEGA 'CL' 0.375 FIX.
    etalcl ~ log(1 + 0.675^2) # Table 2 IIV CL 67.5% CV (RSE 14.4%); ESM 2 $OMEGA 0.375 FIX for CL
    etalogitfdepot ~ 4.17 # Table 2 IIV FIT 4.17 on the logit scale (footnote b); ESM 2 $OMEGA 4.17 FIX for F2 is the variance
    etalksyn ~ log(1 + 3.95^2) # Table 2 IIV KSYN 395% CV (RSE 7.9%)
    etalktr ~ log(1 + 0.952^2) # Table 2 IIV KTR 95.2% CV (RSE 7.70%)
    etalkdeg ~ log(1 + 0.968^2) # Table 2 IIV KDEG 96.8% CV (RSE 29.3%)

    # Residual error: additive on the natural-log scale for each endpoint
    # (Data Analysis; ESM 2 Y = IPRED + W*EPS(1) with $SIGMA 1 FIX).
    expSd <- 1.33; label("Additive residual SD on log(V937 copies/mL)") # Table 2 PK 'Error (log copies/mL)' = 1.33 (RSE 3.88%); ESM 2 THETA(1) 1.33
    expSd_ada_total <- 1.01; label("Additive residual SD on log(ADA dilution factor)") # Table 2 ADA 'Error (log DF)' = 1.01 (RSE 6.76%)
  })
  model({
    # Individual PK parameters (ESM 2 $PK)
    vc <- exp(lvc)
    cl <- exp(lcl + etalcl)
    q <- exp(lq)
    vp <- exp(lvp)
    vp2 <- exp(lvp2)
    fit <- expit(logitfdepot + etalogitfdepot)

    # Individual ADA parameters
    ksyn <- exp(lksyn + etalksyn)
    ec50 <- exp(lec50)
    ktr <- exp(lktr + etalktr)
    kdeg <- exp(lkdeg + etalkdeg)

    # Observed pre-existing ADA is added to the prediction only above the
    # 1:16 detection limit (ESM 2 'IF(NAB0>=16) BASENAB = NAB0').
    ada_base <- ADA_TITER_BASE * (ADA_TITER_BASE >= 16)

    # Serum V937 concentration in copies/L drives B-cell activation
    # (ESM 2 CP_DES = A(1)/V1; EFF = CP^GA / (CP^GA + EC50^GA)).
    cs <- central / vc
    eff <- cs^hill / (cs^hill + ec50^hill)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q / vc
    k31 <- q / vp2

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(bcell1) <- ksyn * eff - ktr * bcell1
    d/dt(bcell2) <- ktr * bcell1 - ktr * bcell2
    d/dt(bcell3) <- ktr * bcell2 - ktr * bcell3
    d/dt(ada) <- ktr * bcell3 - kdeg * ada

    # Intratumoral doses enter the tumour compartment; only the fraction FIT
    # is available for distribution (ESM 2 F2 on compartment 2).
    f(peripheral1) <- fit

    # Observations. Both carry the ESM 2 0.001 offset inside the log
    # (IPRED = LOG(F + 0.001) and LOG(NAB + 0.001 + BASENAB)), which also
    # keeps the log-normal error defined at zero prediction.
    Cc <- central / vc / 1000 + 0.001
    ada_total <- ada + 0.001 + ada_base

    Cc ~ lnorm(expSd)
    ada_total ~ lnorm(expSd_ada_total)
  })
}
