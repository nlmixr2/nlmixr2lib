Wojciechowski_2023_ritlecitinib_base <- function() {
  description <- "Base model (iteration 1 of 3) for oral ritlecitinib: two-compartment disposition with first-order absorption and a direct-response non-stationary (autoinhibitory) Imax effect of the peripheral-compartment concentration on both apparent clearance and bioavailability. Fitted to intensively sampled healthy participants plus sparsely sampled rheumatoid arthritis and alopecia areata patients (186 individuals, 2174 concentrations). Carries rheumatoid arthritis and alopecia areata effects on CL/F, high-fat-meal and 800 mg dose effects on ka, and inflammatory-disease scaling of the IIV and proportional residual error magnitudes. Allometric weight scaling with fixed exponents 0.75 on clearances and 1 on volumes referenced to 70 kg."
  reference <- paste(
    "Wojciechowski J, Purohit VS, Huh Y, Banfield C, Nicholas T.",
    "Evolution of Ritlecitinib Population Pharmacokinetic Models During",
    "Clinical Drug Development. Clin Pharmacokinet. 2023;62(12):1765-1779.",
    "doi:10.1007/s40262-023-01318-3. PMCID PMC10684409.",
    "Base-model column of Table 2; structural model-building steps in",
    "Table S2 of the Electronic Supplementary Material.",
    sep = " "
  )
  vignette <- "Wojciechowski_2023_ritlecitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "ritlecitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "ritlecitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ritlecitinib", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline total body weight; allometrically scales CL/F, Q/F, Vc/F and Vp/F",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline weight (NONMEM column BWT). Reference 70 kg with fixed exponents 0.75 on CL/F and Q/F and 1.00 on Vc/F and Vp/F; Wojciechowski 2023 Table 2 footnote and Methods Sect. 2.3. Base-model analysis population median 77.2 kg (range 46.0-164).",
      source_name        = "BWT"
    ),
    DIS_RA = list(
      description        = "Rheumatoid arthritis patient indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy participant in this base-model analysis population)",
      notes              = "Derived from the NONMEM PTST patient-type column (PTST = 1). Carries a multiplicative effect on CL/F and, jointly with DIS_ALOPECIA_AREATA, defines the inflammatory-disease group that scales the IIV and residual-error magnitudes. The base-model dataset contains only healthy participants (PTST = 0), RA (PTST = 1) and alopecia areata (PTST = 3), so the reference category here is healthy participants.",
      source_name        = "PTST"
    ),
    DIS_ALOPECIA_AREATA = list(
      description        = "Alopecia areata patient indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy participant in this base-model analysis population)",
      notes              = "Derived from the NONMEM PTST patient-type column (PTST = 3). Carries a multiplicative effect on CL/F and, jointly with DIS_RA, defines the inflammatory-disease group that scales the IIV and residual-error magnitudes.",
      source_name        = "PTST"
    ),
    FED_HIGHFAT = list(
      description        = "High-fat-meal administration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted, or food intake not controlled)",
      notes              = "Per-dose-record indicator derived from the NONMEM FOOD column (0 = fasted, 1 = high-fat meal, 2 = not controlled for food); the control stream applies the effect only when FOOD = 1, so both fasted and not-controlled records are the reference. High-fat-meal effect established in the phase I relative-bioavailability / food-effect study B7981003 (NCT02684760).",
      source_name        = "FOOD"
    ),
    DOSE = list(
      description        = "Administered ritlecitinib dose level for the current dose record",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used only as a step-function switch: the control stream applies a separate ka effect when DOSE = 800 exactly (IF (DOSE.EQ.800)), motivated by the prolonged absorption observed at the top single-ascending dose of the first-in-human study. Every dose level below 800 mg is the reference. Base-model dose range 5-800 mg.",
      source_name        = "DOSE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 186L,
    n_studies      = 4L,
    n_observations = 2174L,
    age_range      = "19.0-74.0 years",
    age_median     = "39 years",
    weight_range   = "46.0-164 kg",
    weight_median  = "77.2 kg",
    sex_female_pct = 43.5,
    race_ethnicity = c(White = 91.9, Asian = 3.8, Black = 2.2, Other = 2.2),
    disease_state  = "Healthy participants (74; 39.8%), rheumatoid arthritis patients (42; 22.6%) and alopecia areata patients (70; 37.6%).",
    dose_range     = "5-800 mg/day oral ritlecitinib (single and multiple dose)",
    regions        = "Global (12-trial ritlecitinib clinical development programme; the base-model subset draws on phase I studies in healthy participants plus phase II RA and AA studies)",
    renal_function = "Normal (severe renal impairment participants entered only at the final-model iteration)",
    hepatic_function = "Normal (moderate hepatic impairment participants entered only at the updated-model iteration)",
    notes          = "Demographics from Wojciechowski 2023 Table 1, base-model column. Below-limit-of-quantification observations (LLOQ 0.5 ng/mL in the phase I studies) were excluded during estimation. Studies contributing to the base model are listed in Table S1 of the Electronic Supplementary Material."
  )

  # Implementation notes (see the vignette 'Assumptions and deviations'
  # section for the full justification of each item):
  # * Structure. Wojciechowski 2023 Sect. 3.2 and Fig. 3: two-compartment
  #   disposition, first-order absorption from a depot, and a direct-
  #   response Emax (here Imax) relationship in which the PERIPHERAL
  #   compartment concentration simultaneously inhibits CL/F and raises F.
  #   The deposited NONMEM $DES block (Electronic Supplementary Material,
  #   'NONMEM Control stream (for Updated Model)') writes this as
  #     AUTOINH = 1 + IMAXP*PERIC/(IC50P + PERIC)
  #     INHF    = 1 + (1 - AUTOINH)
  #   with IMAXP negative, so AUTOINH < 1 scales clearance DOWN and
  #   INHF = 2 - AUTOINH > 1 scales the absorbed fraction UP by the same
  #   fractional amount. Sect. 3.2 states the same mechanism is assumed to
  #   act inversely on F (reduced first-pass extraction).
  # * Units. Doses are mg and volumes L, so amounts / volumes are mg/L.
  #   IC50 is reported in ng/mL, and the control stream converts it with
  #   POPIC50P = THETA(6)/1000 to mg/L before comparing it with PERIC.
  #   The same /1000 conversion is reproduced here, and the reported
  #   observation Cc is central/vc*1000 to return ng/mL.
  # * Inlined state expressions. The peripheral and central concentrations
  #   are written out inside each d/dt() rather than routed through named
  #   intermediates: in the nlmixr2 model-function form a named
  #   intermediate that reads an ODE state and feeds a nonlinear term
  #   inside d/dt() can silently evaluate to zero, which would delete the
  #   autoinhibition term without any error.
  # * Absorption covariates. The control stream computes
  #     KAD = KA * COVFOODKA * COVDOSEKA * COVFORMKA
  #   inside $DES purely so that the $TABLE output of KA stays
  #   time-invariant; all three covariates are step functions of
  #   record-level data, so computing KAD before the ODE block is exactly
  #   equivalent. The base model has no formulation term (COVFORMKA = 1).
  # * Inflammatory-disease scaling. THETA(9), THETA(10) and THETA(11)
  #   multiply, respectively, the proportional residual SD, ETA(1) and
  #   ETA(2) for every patient with PTST in {1, 2, 3, 5}, i.e. every
  #   inflammatory-disease patient. In the base-model dataset that set is
  #   exactly {RA, alopecia areata}, so `inflam` is derived here as
  #   DIS_RA + DIS_ALOPECIA_AREATA. Note that the effects act on the
  #   STANDARD DEVIATION (the control stream multiplies ETA and the
  #   residual SD, not their variances) even though Table 2 labels the
  #   rows as effects on omega^2 and sigma^2_pro.
  # * IIV scale. The Table 2 footnote defines the reported "% CV" as
  #   sqrt(omega^2) * 100, so the tabulated 20.1 and 11.3 are 100*omega
  #   and the internal variances are 0.201^2 and 0.113^2. This is
  #   confirmed by the final model's $OMEGAP prior block, which carries
  #   the updated model's 0.198^2 = 0.039204 and 0.125^2 = 0.015625.
  # * Residual error. $SIGMA is 1 FIX and the SD is carried by THETA(8),
  #   so propSd is the THETA(8) point estimate on the fraction scale
  #   (Table 2 reports it as 34.0 "% CV").
  # * No IIV covariance. Table S2 run 332 removed the CL/F-Vc/F
  #   covariance (correlation estimate -0.0675), so the omegas are
  #   diagonal.
  # * ka point estimate. Table 2 prints 7.1 (95% CI 5.70-9.52) for the
  #   base model. Unlike every other row in the table the interval is not
  #   symmetric about the estimate, so one of the three printed numbers
  #   contains a typographical error; the printed point estimate is used
  #   here because confidence bounds do not enter the model. See the
  #   vignette Errata.
  ini({
    # ----- Structural parameters (Wojciechowski 2023 Table 2, base-model column) -----
    lcl   <- log(129);   label("Apparent clearance CL/F at 70 kg, uninhibited (L/h)")                          # Table 2 base: 129 (95% CI 117-141)
    lvc   <- log(145);   label("Apparent central volume Vc/F at 70 kg (L)")                                    # Table 2 base: 145 (95% CI 138-152)
    lq    <- log(0.298); label("Apparent inter-compartmental clearance Q/F at 70 kg (L/h)")                    # Table 2 base: 0.298 (95% CI 0.256-0.340)
    lvp   <- log(4.43);  label("Apparent peripheral volume Vp/F at 70 kg (L)")                                 # Table 2 base: 4.43 (95% CI 4.07-4.79)
    lka   <- log(7.1);   label("First-order absorption rate constant ka (1/h)")                                # Table 2 base: 7.1 (95% CI 5.70-9.52); interval not symmetric about the estimate, see vignette Errata

    # ----- Non-stationary (autoinhibition) parameters -----
    # imax is a signed fractional effect and cannot be log-transformed.
    imax  <- -0.559;     label("Maximum non-stationary fractional effect Imax,P on CL/F and F (unitless)")     # Table 2 base: -0.559 (95% CI -0.598 to -0.520)
    lic50 <- log(11.8);  label("Peripheral concentration giving half the maximum non-stationary effect IC50,P (ng/mL)")  # Table 2 base: 11.8 (95% CI 9.02-14.6)

    # ----- Allometric exponents (fixed, Table 2 footnote) -----
    e_wt_cl_q  <- fixed(0.75); label("Allometric WT exponent shared by CL/F and Q/F (unitless)")               # Table 2 footnote: exponent 0.75 on clearance parameters, reference 70 kg
    e_wt_vc_vp <- fixed(1.00); label("Allometric WT exponent shared by Vc/F and Vp/F (unitless)")              # Table 2 footnote: exponent 1 on volume parameters, reference 70 kg

    # ----- Categorical covariate effects on CL/F -----
    e_ra_cl               <- -0.439; label("Fractional effect of rheumatoid arthritis on CL/F (unitless)")     # Table 2 base: -0.439 (95% CI -0.554 to -0.324); ESM updated-model $PK THETA(12) PTSTCL1
    e_alopecia_areata_cl  <- -0.258; label("Fractional effect of alopecia areata on CL/F (unitless)")          # Table 2 base: -0.258 (95% CI -0.398 to -0.118); ESM updated-model $PK THETA(14) PTSTCL3

    # ----- Categorical covariate effects on ka -----
    e_fed_highfat_ka <- -0.718; label("Fractional effect of a high-fat meal on ka (unitless)")                 # Table 2 base: -0.718 (95% CI -0.781 to -0.655); ESM updated-model $PK THETA(18) FOODKA1
    e_dose800_ka     <- -0.815; label("Fractional effect of the 800 mg dose level on ka (unitless)")           # Table 2 base: -0.815 (95% CI -0.868 to -0.762); ESM updated-model $PK THETA(19) DOSEKA800

    # ----- Inflammatory-disease scaling of the variability magnitudes -----
    # These multiply the SD (eta and residual SD), not the variance; see the
    # implementation notes above and the ESM $PK / $ERROR blocks.
    e_inflam_iiv_cl <- 1.53;  label("Fractional effect of inflammatory disease on the SD of IIV on CL/F (unitless)")  # Table 2 base: 1.53 (95% CI 0.801-2.26); ESM $PK THETA(10) VARCLPTST
    e_inflam_iiv_vc <- 4.24;  label("Fractional effect of inflammatory disease on the SD of IIV on Vc/F (unitless)")  # Table 2 base row 2 of the duplicated omega^2 label: 4.24 (95% CI 1.93-6.55); ESM $PK THETA(11) VARV2PTST
    e_inflam_propsd <- 0.522; label("Fractional effect of inflammatory disease on the proportional residual SD (unitless)")  # Table 2 base: 0.522 (95% CI 0.446-0.598); ESM $ERROR THETA(9) RUVPROPTST

    # ----- Between-subject variability (diagonal; Table S2 run 332 removed the covariance) -----
    etalcl ~ 0.040401  # Table 2 base: 20.1 percent-CV column (= 100*omega) per the Table 2 footnote, so omega^2 = 0.201^2
    etalvc ~ 0.012769  # Table 2 base: 11.3 percent-CV column (= 100*omega) per the Table 2 footnote, so omega^2 = 0.113^2

    # ----- Residual error -----
    propSd <- 0.340; label("Proportional residual SD for the reference (healthy participant) group (fraction)")  # Table 2 base: 34.0 "% CV"; ESM $ERROR THETA(8) RUVPRO with $SIGMA 1 FIX
  })
  model({
    ref_wt <- 70  # Table 2 footnote: allometric reference weight

    # ----- Derived covariate terms -----
    # Inflammatory-disease group = PTST in {1, 2, 3, 5}. Only RA and
    # alopecia areata are present in the base-model analysis dataset.
    inflam <- DIS_RA + DIS_ALOPECIA_AREATA

    covptstcl  <- 1 + e_ra_cl * DIS_RA + e_alopecia_areata_cl * DIS_ALOPECIA_AREATA
    covvarcl   <- 1 + e_inflam_iiv_cl * inflam
    covvarvc   <- 1 + e_inflam_iiv_vc * inflam
    covfoodka  <- 1 + e_fed_highfat_ka * FED_HIGHFAT
    covdoseka  <- 1 + e_dose800_ka * (DOSE == 800)

    # ----- Individual parameters -----
    cl <- exp(lcl + etalcl * covvarcl) * (WT / ref_wt)^e_wt_cl_q * covptstcl
    vc <- exp(lvc + etalvc * covvarvc) * (WT / ref_wt)^e_wt_vc_vp
    q  <- exp(lq)  * (WT / ref_wt)^e_wt_cl_q
    vp <- exp(lvp) * (WT / ref_wt)^e_wt_vc_vp

    # Absorption rate after the record-level covariate step functions.
    kad <- exp(lka) * covfoodka * covdoseka

    # IC50 is reported in ng/mL; the peripheral concentration below is in
    # mg/L, so convert exactly as the control stream does (THETA(6)/1000).
    ic50 <- exp(lic50) / 1000

    # ----- ODE system -----
    # AUTOINH = 1 + imax*Cperi/(ic50 + Cperi) multiplies CL (imax < 0, so
    # clearance falls) and INHF = 2 - AUTOINH multiplies the absorbed
    # amount entering the central compartment (so F rises by the same
    # fraction). Concentrations are written inline; see implementation
    # notes.
    d/dt(depot) <- -kad * depot
    d/dt(central) <-
      kad * depot *
        (2 - (1 + imax * (peripheral1 / vp) / (ic50 + peripheral1 / vp))) -
      q * (central / vc) + q * (peripheral1 / vp) -
      cl * (1 + imax * (peripheral1 / vp) / (ic50 + peripheral1 / vp)) *
        (central / vc)
    d/dt(peripheral1) <- q * (central / vc) - q * (peripheral1 / vp)

    # ----- Observation and residual error -----
    # central is mg and vc is L, so central/vc is mg/L = 1000 ng/mL.
    Cc <- central / vc * 1000
    propSdEff <- propSd * (1 + e_inflam_propsd * inflam)
    Cc ~ prop(propSdEff)
  })
}
