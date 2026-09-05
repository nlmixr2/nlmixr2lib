Rietveld_2024_irinotecan <- function() {
  description <- paste(
    "Joint five-compartment population PK model for intraperitoneally (IP)",
    "administered irinotecan and its active metabolite SN-38 in patients with",
    "peritoneal metastases from colorectal cancer (INTERACT phase I trial,",
    "Rietveld 2024). Irinotecan is instilled into the peritoneal cavity",
    "(depot_ip) and leaves it by two parallel first-order routes: absorption",
    "into plasma (ka_ip) and in-situ conversion to SN-38 by peritoneal",
    "carboxylesterases (kmet_ip). Irinotecan plasma disposition is",
    "two-compartment (central + peripheral1, CL / Vc / Q / Vp). A fixed",
    "fraction fm = 3% of the irinotecan cleared from plasma is converted to",
    "SN-38 systemically; the remaining (1 - fm) is eliminated. SN-38 has a",
    "one-compartment plasma disposition (central_sn38, CL_M / V4) fed by both",
    "the systemic conversion and absorption from the peritoneal SN-38 pool",
    "(depot_ip_sn38, V5 fixed at 487 L), the latter measured directly in",
    "peritoneal fluid. Weight enters SN-38 central volume as a median-83.4 kg",
    "normalized power function and gamma-glutamyltransferase enters SN-38",
    "plasma clearance as a median-32 U/L normalized power function.",
    "Inter-individual variability on SN-38 clearance and SN-38 central volume",
    "only; irinotecan carries no IIV, so its between-subject spread is",
    "absorbed into the residual error. Residual error is additive on the",
    "natural-log scale (equivalent to proportional in linear space) and is",
    "estimated separately for the three observed streams. Concentrations are",
    "molar, so parent-to-metabolite transfer is equimolar; see the vignette",
    "Errata for the derived SN-38 AUC summary values that do not reproduce",
    "from the published parameters."
  )
  reference <- paste(
    "Rietveld PCS, Sassen SDT, Guchelaar NAD, van Eerden RAG, de Boer NL,",
    "van den Heuvel TBM, Burger JWA, Mathijssen RHJ, Koch BCP, Koolen SLW.",
    "Population pharmacokinetics of intraperitoneal irinotecan and SN-38 in",
    "patients with peritoneal metastases from colorectal origin.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(6):1006-1016.",
    "doi:10.1002/psp4.13136."
  )
  vignette <- "Rietveld_2024_irinotecan"

  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  dosing <- c("depot_ip")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE means checked against Rietveld 2024
  # Figure 1 (final structural model diagram) and the Table 2 footnote.
  #
  # Both peritoneal states are depot states and carry the vocabulary's
  # designated depot label "administration site" -- the peritoneal cavity is
  # where irinotecan is instilled, and depot_ip_sn38 holds SN-38 formed in
  # situ at that same site. Note that SN-38 (unlike irinotecan) WAS assayed
  # there: 263 peritoneal-fluid samples. The specimen vocabulary in
  # R/conventions.R has no "peritoneal fluid" matrix, so that detail lives
  # here and in the vignette rather than in the specimen field.
  compartmentData <- list(
    depot_ip = list(
      analyte = "irinotecan", units = "umol", specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "irinotecan", units = "umol", specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "irinotecan", units = "umol", specimen = "plasma", verified = TRUE
    ),
    depot_ip_sn38 = list(
      analyte = "SN-38", units = "umol", specimen = "administration site", verified = TRUE
    ),
    central_sn38 = list(
      analyte = "SN-38", units = "umol", specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      source_name = "WT",
      description = "Body weight, entering SN-38 central volume of distribution as a power function normalized to the 83.4 kg population median",
      units = "kg",
      type = "continuous",
      reference_value = 83.4,
      notes = paste(
        "Rietveld 2024 Equation 2 normalizes to 83.4 kg and Figure 2 plots",
        "Q1 / median / Q3 = 70 / 83.4 / 98.5 kg. Table 1 reports a baseline",
        "weight median of 79.7 kg (range 59-105); the 83.4 kg value printed",
        "in the equation is used here. The authors state the model is only",
        "applicable over 60-100 kg."
      )
    ),
    GGT = list(
      source_name = "GGT",
      description = "Serum gamma-glutamyltransferase, entering SN-38 plasma clearance as a power function normalized to the 32 U/L population median",
      units = "U/L",
      type = "continuous",
      reference_value = 32,
      notes = paste(
        "Rietveld 2024 Equation 3 normalizes to 32 U/L and Figure 3 plots",
        "Q1 / median / Q3 = 26 / 32 / 62 U/L. GGT is not tabulated in",
        "Table 1; the quantiles come from the Figure 3 legend."
      )
    )
  )

  # Screened in the univariate covariate analysis but not retained in the
  # final model (Rietveld 2024 Results, "Covariates"). Documented for
  # provenance only; none is referenced in model().
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex, tested on irinotecan plasma clearance",
      units = "unitless",
      type = "categorical",
      notes = "Significant univariately (lower CL in females) but removed during backward elimination (p > 0.01)."
    ),
    AGE = list(
      description = "Age, tested on irinotecan plasma clearance",
      units = "years",
      type = "continuous",
      notes = "Significant univariately but removed during backward elimination (p > 0.01)."
    ),
    BMI = list(
      description = "Body mass index, tested on SN-38 plasma volume of distribution",
      units = "kg/m^2",
      type = "continuous",
      notes = "Significant univariately; dropped for collinearity with WT, which gave the largest OFV drop."
    ),
    BSA = list(
      description = "Body surface area, tested on SN-38 plasma volume of distribution",
      units = "m^2",
      type = "continuous",
      notes = "Significant univariately; dropped for collinearity with WT, which gave the largest OFV drop."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 18L,
    n_studies = 1L,
    n_observations = 855L,
    age_range = "42-77 years (median 64)",
    weight_range = "59-105 kg (median 79.7)",
    sex_female_pct = 33.3,
    disease_state = paste(
      "Adults with extensive peritoneal metastases (Peritoneal Cancer Index",
      "> 20, median 29, range 17-39) from histologically proven colorectal",
      "cancer, ineligible for cytoreductive surgery with hyperthermic",
      "intraperitoneal chemotherapy because of extensive peritoneal disease",
      "or an unresectable primary tumour, with no extra-abdominal",
      "metastases. ECOG performance status 0 (n = 12) or 1 (n = 6).",
      "Ascites fluid was drained before each IP instillation."
    ),
    dose_range = paste(
      "Irinotecan 50 mg (n = 4), 75 mg (n = 9) or 100 mg (n = 5)",
      "administered intraperitoneally every 2 weeks over a 1.5 h instillation",
      "of fluid prewarmed to 37 C through a peritoneal access port, given",
      "concurrently with systemic FOLFOX plus bevacizumab. The maximum",
      "tolerated dose was 75 mg two-weekly. Data from the first two cycles",
      "were pooled and expressed as time after dose; no accumulation was",
      "observed."
    ),
    regions = "Netherlands (Erasmus MC Rotterdam, Catharina Hospital Eindhoven).",
    notes = paste(
      "INTERACT phase I trial, classic 3 + 3 dose escalation. Samples:",
      "irinotecan in plasma (334), SN-38 in plasma (258) and SN-38 in",
      "peritoneal fluid (263), 855 in total; the Abstract instead reports",
      "588 plasma and 267 peritoneal-fluid samples. Plasma and peritoneal",
      "samples were drawn predose and at 0.5, 1, 1.5, 2, 3, 4, 6, 22.5 and",
      "46.5 h after the start of the infusion, plus a 45 min plasma sample.",
      "LLOQ 1 ng/mL for irinotecan and SN-38 in plasma and 2 ng/mL for SN-38",
      "in peritoneal fluid. IP irinotecan itself was NOT measured, because",
      "peritoneal samples were withdrawn through the same port used for the",
      "instillation. NONMEM 7.5 with FOCE+I on natural-log-transformed",
      "concentrations; Pirana 2.9.9, R 4.2.1, Xpose4 4.7.2. Nonparametric",
      "bootstrap n = 1000. Inter-occasion variability was tested and gave no",
      "significant improvement; M3 BLOQ handling gave equal estimates and",
      "was not retained. Conditional number 400.41."
    )
  )

  ini({
    # ================================================================
    # Irinotecan: peritoneal instillation and plasma disposition.
    # Rietveld 2024 Table 2 (final parameter estimates) and Figure 1
    # (final structural model). CL is the TOTAL irinotecan plasma
    # clearance; Figure 1 splits it into MR * CL (conversion to SN-38)
    # and (1 - MR) * CL (elimination).
    # ================================================================
    lka_ip <- log(1.02)
    label("Irinotecan first-order absorption rate constant from the peritoneal cavity into plasma, Ka (1/h)")   # Table 2 row Ka = 1.02 1/h, RSE 9%, bootstrap 95% CI 0.85-1.2
    lcl <- log(33.2)
    label("Irinotecan total plasma clearance, CL (L/h)")                                                        # Table 2 row CL = 33.2 L/h, RSE 6%, bootstrap 95% CI 29.4-36.99
    lvc <- log(225)
    label("Irinotecan central volume of distribution, V2 (L)")                                                  # Table 2 row V2 = 225 L, RSE 11%, bootstrap 95% CI 175.44-269.42
    lq <- log(13.9)
    label("Irinotecan inter-compartmental clearance, Q (L/h)")                                                  # Table 2 row Q = 13.9 L/h, RSE 17%, bootstrap 95% CI 9.46-21.24
    lvp <- log(119)
    label("Irinotecan peripheral volume of distribution, V3 (L)")                                               # Table 2 row V3 = 119 L, RSE 14%, bootstrap 95% CI 87.06-159.97

    # Fraction of the irinotecan cleared from plasma that is converted to
    # SN-38 ("metabolic ratio"). Fixed at 3% from a pooled analysis of three
    # phase I IV-irinotecan studies (168 PK datasets in 107 patients) because
    # it could not be estimated from the INTERACT data.
    fm <- fixed(0.03)
    label("Fraction of irinotecan plasma clearance converted to SN-38, MR (unitless)")                          # Table 2 row MR = 0.03 FIX; Results: "fixed at 3%, based on the literature"

    # ================================================================
    # Peritoneal conversion of irinotecan to SN-38. Rietveld 2024
    # Figure 1: the CLPM arrow runs from IRI IP directly to SN-38 IP,
    # in parallel with the Ka arrow into plasma. Reported in 1/h, so it
    # is a first-order rate constant acting on the amount in depot_ip.
    # ================================================================
    lkmet_ip <- log(0.118)
    label("Peritoneal irinotecan-to-SN-38 conversion rate constant, CLPM (1/h)")                                # Table 2 row CLPM = 0.118 1/h, RSE 14%, bootstrap 95% CI 0.09-0.156

    # ================================================================
    # SN-38 disposition. Rietveld 2024 Table 2 and Figure 1.
    # ================================================================
    lcl_sn38 <- log(46)
    label("SN-38 plasma clearance, CL_M (L/h)")                                                                 # Table 2 row CL_M = 46 L/h, RSE 12%, bootstrap 95% CI 36.57-57.15
    lvc_sn38 <- log(15.9)
    label("SN-38 plasma volume of distribution, V4 (L)")                                                        # Table 2 row V4 = 15.9 L, RSE 26%, bootstrap 95% CI 5.187-27.50

    # Ka2 is tabulated in L/h, i.e. an absorption CLEARANCE out of the
    # peritoneal SN-38 pool, so the transfer rate constant is Ka2 / V5.
    # V5 was estimable but highly correlated with CLPM, so the authors
    # fixed it at its own estimated value for model stability.
    lq_ip_sn38 <- log(4.68)
    label("SN-38 absorption clearance from peritoneal fluid into plasma, Ka2 (L/h)")                            # Table 2 row Ka2 = 4.68 L/h, RSE 15%, bootstrap 95% CI 3.49-6.3
    lv_ip_sn38 <- fixed(log(487))
    label("SN-38 peritoneal volume of distribution, V5 (L)")                                                    # Table 2 row V5 = 487 L FIX; Results: fixed to its own estimated value because of high covariance with CLPM

    # ================================================================
    # Covariate effects. Rietveld 2024 Equations 2 and 3, both
    # median-normalized power functions on the SN-38 parameters.
    # ================================================================
    e_wt_vc_sn38 <- 5.31
    label("Exponent of (WT / 83.4 kg) on SN-38 plasma volume of distribution (unitless)")                       # Table 2 Covariates row WT = 5.31, RSE 30%, bootstrap 95% CI 1.63-11.31; Equation 2; dOFV -15
    e_ggt_cl_sn38 <- -0.26
    label("Exponent of (GGT / 32 U/L) on SN-38 plasma clearance (unitless)")                                    # Table 2 Covariates row GGT = -0.26, RSE 35%, bootstrap 95% CI -0.37 to -0.026; Equation 3; dOFV -20

    # ================================================================
    # Inter-individual variability. Table 2 reports these two IIVs as
    # CV%, so the log-scale variance is omega^2 = log(1 + CV^2).
    # The Table 2 footnote also defines IIV-CL and IIV-V2, but neither
    # appears in the table and the Results state IIV was retained only
    # on SN-38 plasma CL and the SN-38 central volume.
    # ================================================================
    etalcl_sn38 ~ 0.134217   # Table 2 IIV row CL_M = 37.9 CV% (RSE 14%, shrinkage 0.1%); log(1 + 0.379^2)
    etalvc_sn38 ~ 0.539830   # Table 2 IIV row V4 = 84.6 CV% (RSE 34%, shrinkage 20%); log(1 + 0.846^2)

    # ================================================================
    # Residual unexplained variability. Concentrations were natural-log
    # transformed and the error was additive on that log scale, which is
    # proportional error in nlmixr2's linear space. Values are taken as
    # log-scale standard deviations; see the vignette Errata for the
    # SD-versus-variance reading.
    # ================================================================
    propSd <- 0.427
    label("Proportional residual error for irinotecan in plasma (fraction)")                                    # Table 2 row Add ERR CMT 2 = 0.427 (RSE 8%, shrinkage 0.1%), additive on log-transformed concentration
    propSd_sn38 <- 0.247
    label("Proportional residual error for SN-38 in plasma (fraction)")                                         # Table 2 row Add ERR CMT 4 = 0.247 (RSE 12%, shrinkage 6%)
    propSd_Cip_sn38 <- 0.587
    label("Proportional residual error for SN-38 in peritoneal fluid (fraction)")                               # Table 2 row Add ERR CMT 5 = 0.587 (RSE 9%, shrinkage 0.1%)
  })

  model({
    # ---- Individual parameters -------------------------------------
    ka_ip <- exp(lka_ip)
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    kmet_ip <- exp(lkmet_ip)
    v_ip_sn38 <- exp(lv_ip_sn38)
    q_ip_sn38 <- exp(lq_ip_sn38)

    # Rietveld 2024 Equation 2: weight on the SN-38 central volume.
    vc_sn38 <- exp(lvc_sn38 + etalvc_sn38) * (WT / 83.4)^e_wt_vc_sn38
    # Rietveld 2024 Equation 3: GGT on SN-38 plasma clearance.
    cl_sn38 <- exp(lcl_sn38 + etalcl_sn38) * (GGT / 32)^e_ggt_cl_sn38

    # ---- ODEs, Rietveld 2024 Figure 1 -------------------------------
    # Irinotecan leaves the peritoneal cavity by absorption into plasma
    # (Ka) and by in-situ conversion to SN-38 (CLPM), in parallel.
    d/dt(depot_ip) <- -ka_ip * depot_ip - kmet_ip * depot_ip
    d/dt(central) <- ka_ip * depot_ip -
      (cl / vc) * central -
      (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1
    # Peritoneal SN-38 pool: formed in situ, absorbed into SN-38 plasma.
    d/dt(depot_ip_sn38) <- kmet_ip * depot_ip -
      (q_ip_sn38 / v_ip_sn38) * depot_ip_sn38
    # SN-38 plasma: fed by systemic conversion (MR * CL) and by
    # absorption from the peritoneal pool.
    d/dt(central_sn38) <- fm * (cl / vc) * central +
      (q_ip_sn38 / v_ip_sn38) * depot_ip_sn38 -
      (cl_sn38 / vc_sn38) * central_sn38

    # ---- Observations ------------------------------------------------
    Cc <- central / vc
    Cc_sn38 <- central_sn38 / vc_sn38
    Cip_sn38 <- depot_ip_sn38 / v_ip_sn38

    Cc ~ prop(propSd)
    Cc_sn38 ~ prop(propSd_sn38)
    Cip_sn38 ~ prop(propSd_Cip_sn38)
  })
}
