Bukkems_2022_factorviii_thrombinPotential <- function() {
  description <- paste(
    "Sequential PK-PD model relating factor VIII (FVIII) activity to",
    "normalized thrombin potential (% of normal pooled plasma) measured with",
    "the Nijmegen haemostasis assay, after an intravenous bolus of standard",
    "half-life FVIII concentrates in 24 children and adults with severe,",
    "moderate or mild hemophilia A (Bukkems 2022). The PK layer is the",
    "two-compartment allometric model of Bukkems 2022 (see",
    "modellib('Bukkems_2022_factorviii')); the PD layer is an",
    "additive-baseline sigmoid Emax model, E = Ebase + Emax * C^n / (EC50^n +",
    "C^n), driven by total one-stage FVIII activity (exogenous plus the",
    "patient's endogenous baseline), with baseline 37.5% of normal pooled",
    "plasma, EC50 13.9 IU/dL and Hill coefficient 1.62. Emax is 72.5% of",
    "normal pooled plasma in a 75 kg patient with severe or moderate",
    "hemophilia, scales with (WT / 75)^-0.28, and is 0.709x in mild",
    "hemophilia A. IIV on baseline (41.8%), EC50 (88.0%) and Emax (22.9%);",
    "additive residual error 8.62% of normal pooled plasma.",
    sep = " "
  )
  reference <- paste(
    "Bukkems LH, Valke LLFG, Barteling W, Laros-van Gorkom BAP, Blijlevens",
    "NMA, Cnossen MH, van Heerde WL, Schols SEM, Mathot RAA. Combining factor",
    "VIII levels and thrombin/plasmin generation: a population",
    "pharmacokinetic-pharmacodynamic model for patients with haemophilia A.",
    "Br J Clin Pharmacol. 2022;88(6):2757-2768. doi:10.1111/bcp.15185.",
    "PMCID:PMC9304184.",
    sep = " "
  )
  vignette <- "Bukkems_2022_factorviii_thrombin_plasmin"
  units <- list(
    time = "h",
    dosing = "IU",
    concentration = "IU/dL (observation Cc is FVIII activity; the PD output thrombinPotential is normalized thrombin potential as a percentage of normal pooled plasma)"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "PK layer: allometric power scaling at a 75 kg reference, exponent 0.75 on CL and 1 on V1; Q and V2 not scaled (Bukkems 2022 Table 2 footer). PD layer: Emax = theta_Emax * (Weight / 75)^-0.28 * 0.709^(mild hemophilia) (Table 2 footer; 'Coefficient bodyweight on Emax' -0.28, bootstrap 95% CI -0.42 to -0.11). Supplementary Figure S4 confirms the 75 kg typical patient. PD cohort median 85 kg (range 11-134 kg; Table 1).",
      source_name = "Weight"
    ),
    VWF = list(
      description = "Von Willebrand factor ristocetin cofactor activity (VWF:Act)",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "PK layer: power effect (VWF / 100)^-0.52 on CL (Table 2). PD cohort median 97% (range 55-142%; Table 1). Must be supplied for every simulated patient (the missing-VWF parameter of Supplementary Methods Eq. 6 is not reported).",
      source_name = "VWF:Act"
    ),
    FORM_FVIII_RECOMB_FULLLENGTH = list(
      description = "Full-length recombinant FVIII product indicator: 1 = full-length recombinant concentrate (Advate, Helixate, Kogenate), 0 = plasma-derived or B-domain-deleted recombinant concentrate",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (plasma-derived or B-domain-deleted recombinant product)",
      notes = "PK layer: multiplicative theta^flag effects 1.27 on CL and 1.17 on V1 (Table 2). Per-profile.",
      source_name = "recombinant product"
    ),
    ADA_POS = list(
      description = "FVIII inhibitor status by the Nijmegen low-titre inhibitor assay (NLTIA): 1 = positive, 0 = negative. Mapped onto the canonical ADA_POS column per the neutralising-antibody alias in inst/references/covariate-columns.md.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (NLTIA negative)",
      notes = "PK layer: multiplicative theta^flag effects 1.49 on CL and 1.14 on V1 (Table 2). 4 of 24 PD patients were NLTIA-positive. No inhibitor effect was retained on the thrombin potential PD parameters.",
      source_name = "Positive NLTIA"
    ),
    ASSAY_OSA = list(
      description = "One-stage activated partial thromboplastin time clotting assay (OSA) indicator (1 = OSA, 0 = chromogenic substrate assay, CSA)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CSA) per the canonical column definition; in THIS model the structural PK parameters are on the OSA scale, so the assay correction and the CSA residual-error magnitudes are applied to the CSA (ASSAY_OSA = 0) rows",
      notes = "Affects the FVIII activity observation Cc only (CSA activity = OSA activity * 1.20, with IIV; Table 2). The PD layer is always driven by the one-stage-scale FVIII activity: the PD parameters were estimated against the individual OSA PK predictions (Discussion). Per-observation (per-row).",
      source_name = "OSA / CSA"
    ),
    FVIII_BL = list(
      description = "Endogenous baseline FVIII activity: the lowest FVIII activity ever measured in the patient",
      units = "IU/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Added (x100, to IU/dL) to the predicted exogenous FVIII in the observation (Supplementary Methods Eq. 2) and in the PD driver: Equation 1 defines C as 'the FVIII activity level (including endogenous baseline)'. Time-fixed per patient. Severe < 0.01, moderate 0.01-0.05, mild > 0.05-0.40 IU/mL (Methods 2.1).",
      source_name = "baseline FVIII"
    ),
    DIS_HEMOPHILIA_A_MILD = list(
      description = "Mild hemophilia A indicator: 1 = mild (baseline FVIII activity > 5-40 IU/dL), 0 = moderate (1-5 IU/dL) or severe (< 1 IU/dL)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (moderate or severe hemophilia A)",
      notes = "Multiplicative 0.709^flag on the thrombin potential Emax (Table 2 footer; row 'Mild haemophilia on Emax (% of severe)' = 70.9, RSE 15.9%). The row label says '% of severe' but the equation indicator is mild versus not mild, so moderate patients share the severe reference; the 4 moderate PK patients (3 in the PD analysis) are a small group. Severity bands from Methods 2.1. 6 of 24 PD patients were mild (Table 1). Patient-level, time-fixed.",
      source_name = "mild hemophilia"
    )
  )

  compartmentData <- list(
    central = list(analyte = "factor VIII", units = "IU", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "factor VIII", units = "IU", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    age_range = "2-77 years",
    age_median = "51 years",
    weight_range = "11-134 kg",
    weight_median = "85 kg",
    sex_female_pct = NULL,
    race_ethnicity = NULL,
    disease_state = "Hemophilia A: 15 severe, 3 moderate and 6 mild patients (PD study column of Table 1). 4 (16%) NLTIA-positive.",
    dose_range = "Intravenous bolus of a standard half-life FVIII concentrate, median 32.1 IU/kg (range 12.8-90.9); 25 PD profiles.",
    regions = "The Netherlands (Radboud University Medical Center, Nijmegen; single centre).",
    notes = "Normalized thrombin potential from 24 of the 30 PK patients; 8 values (3.2%) could not be determined at very low FVIII activity and were excluded (Results 3.1). Each NHA value is the mean of two measurements, normalized to the normal pooled plasma run alongside the patient sample. Sequential analysis: the individual empirical Bayes PK parameters from the Bukkems 2022 PK model were fixed as input to the PD estimation; the PK layer here is that population PK model, so simulations propagate PK variability as well."
  )

  ini({
    # ---- PK layer: Bukkems 2022 Table 2 (see Bukkems_2022_factorviii) ----
    lcl <- log(1.69) ; label("Clearance CL at the reference 75 kg (dL/h)")                 # Table 2: CL = 1.69 dL/h (RSE 11.1%)
    lvc <- log(27.7) ; label("Central volume V1 at the reference 75 kg (dL)")              # Table 2: V1 = 27.7 dL (RSE 5.8%)
    lq  <- log(2.27) ; label("Inter-compartmental clearance Q (dL/h)")                     # Table 2: Q = 2.27 dL/h (RSE 44.5%)
    lvp <- log(5.63) ; label("Peripheral volume V2 (dL)")                                  # Table 2: V2 = 5.63 dL (RSE 19.2%)

    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of (WT / 75) on CL (unitless)")   # Table 2 footer: CL = theta_CL * (Weight / 75)^0.75 ...; Supplementary Eq. 3
    e_wt_vc <- fixed(1)    ; label("Allometric exponent of (WT / 75) on V1 (unitless)")   # Table 2 footer: V1 = theta_V1 * (Weight / 75) ...; Supplementary Eq. 3

    e_vwf_cl <- -0.52 ; label("Power exponent of (VWF / 100) on CL (unitless)")            # Table 2: VWF exponent on CL = -0.52 (RSE 26.6%)
    e_form_fviii_recomb_fulllength_cl <- 1.27 ; label("Multiplicative factor on CL for full-length recombinant products (unitless)")   # Table 2: Full-length recombinant product on CL (%) = 127
    e_form_fviii_recomb_fulllength_vc <- 1.17 ; label("Multiplicative factor on V1 for full-length recombinant products (unitless)")   # Table 2: Full-length recombinant product on V1 (%) = 117
    e_ada_pos_cl <- 1.49 ; label("Multiplicative factor on CL for NLTIA inhibitor-positive patients (unitless)")   # Table 2: Positive NLTIA on CL (%) = 149
    e_ada_pos_vc <- 1.14 ; label("Multiplicative factor on V1 for NLTIA inhibitor-positive patients (unitless)")   # Table 2: Positive NLTIA on V1 (%) = 114

    e_assay_osa_cc     <- 1.20   ; label("Typical multiplicative factor on predicted FVIII activity for CSA- vs OSA-assayed samples (unitless)")   # Table 2: Correction factor CSA = 1.20 (RSE 3.5%)
    e_assay_osa_propsd <- 0.9375 ; label("Multiplicative factor on the proportional residual SD for CSA- vs OSA-assayed samples (unitless)") # Table 2: proportional error CSA 10.5% / OSA 11.2% = 0.9375
    e_assay_osa_addsd  <- 1.0313 ; label("Multiplicative factor on the additive residual SD for CSA- vs OSA-assayed samples (unitless)")     # Table 2: additive error CSA 4.28 / OSA 4.15 IU/dL = 1.0313

    # ---- PD layer: normalized thrombin potential (Table 2) ----
    le0   <- log(37.5) ; label("Baseline normalized thrombin potential (% of normal pooled plasma)")                       # Table 2: Baseline effect = 37.5% of NPP (RSE 13.1%)
    lemax <- log(72.5) ; label("Maximal FVIII effect on normalized thrombin potential at 75 kg, severe or moderate (% of normal pooled plasma)")   # Table 2: Maximal effect (Emax) = 72.5% of NPP (RSE 9.5%)
    lec50 <- log(13.9) ; label("FVIII activity producing half the maximal thrombin potential effect (IU/dL)")             # Table 2: EC50 = 13.9 IU/dL (RSE 21.2%)
    lhill <- log(1.62) ; label("Hill coefficient of the thrombin potential Emax relationship (unitless)")                 # Table 2: Hill coefficient = 1.62 (RSE 20.8%)

    e_wt_emax <- -0.28 ; label("Power exponent of (WT / 75) on the thrombin potential Emax (unitless)")                  # Table 2: Coefficient bodyweight on Emax = -0.28 (RSE 21.0%); bootstrap 95% CI -0.42 to -0.11
    e_dis_hemophilia_a_mild_emax <- 0.709 ; label("Multiplicative factor on the thrombin potential Emax for mild hemophilia A (unitless)")   # Table 2: Mild haemophilia on Emax (% of severe) = 70.9 (RSE 15.9%)

    # ---- Random effects: omega^2 = log(CV^2 + 1) (Table 2 abbreviations) ----
    etalcl + etalvc ~ c(0.156785,
                        0.026770, 0.024045)   # Table 2: IIV CL = 41.2%, IIV V1 = 15.6%, correlation IIV V1 and CL = 43.6%
    etae_assay_osa_cc ~ 0.032236              # Table 2: IIV correction factor CSA = 18.1%; log(1 + 0.181^2)
    etale0   ~ 0.161033                       # Table 2: IIV baseline effect (thrombin potential) = 41.8% [Shr 15.7%]; log(1 + 0.418^2)
    etalec50 ~ 0.573462                       # Table 2: IIV EC50 (thrombin potential) = 88.0% [Shr 15.5%]; log(1 + 0.880^2)
    etalemax ~ 0.051112                       # Table 2: IIV Emax (thrombin potential) = 22.9% [Shr 17.5%]; log(1 + 0.229^2)

    # ---- Residual error ----
    propSd <- 0.112 ; label("Proportional residual SD for OSA-assayed FVIII activity (fraction)")   # Table 2: Proportional error OSA = 11.2%
    addSd  <- 4.15  ; label("Additive residual SD for OSA-assayed FVIII activity (IU/dL)")          # Table 2: Additive error OSA = 4.15 IU/dL

    addSd_thrombinPotential <- 8.62 ; label("Additive residual SD for normalized thrombin potential (% of normal pooled plasma)")   # Table 2: Additive error = 8.62% of NPP (RSE 12.2%)
  })

  model({
    # ---- PK ----
    ws <- WT / 75

    cl <- exp(lcl + etalcl) * ws^e_wt_cl * (VWF / 100)^e_vwf_cl *
      e_form_fviii_recomb_fulllength_cl^FORM_FVIII_RECOMB_FULLLENGTH * e_ada_pos_cl^ADA_POS
    vc <- exp(lvc + etalvc) * ws^e_wt_vc *
      e_form_fviii_recomb_fulllength_vc^FORM_FVIII_RECOMB_FULLLENGTH * e_ada_pos_vc^ADA_POS
    q  <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    fviii    <- central / vc
    fviii_bl <- 100 * FVIII_BL

    assay_csa <- 1 - ASSAY_OSA
    cf_csa    <- e_assay_osa_cc * exp(etae_assay_osa_cc)
    Cc <- fviii * cf_csa^assay_csa + fviii_bl

    propSd_eff <- propSd * e_assay_osa_propsd^assay_csa
    addSd_eff  <- addSd * e_assay_osa_addsd^assay_csa
    Cc ~ add(addSd_eff) + prop(propSd_eff)

    # ---- PD: normalized thrombin potential (% of normal pooled plasma) ----
    # Table 2 footer: E = Ebase + Emax * C^n / (EC50^n + C^n), with
    # Emax = theta_Emax * (Weight / 75)^-0.28 * 0.709^(mild hemophilia) and C
    # the one-stage FVIII activity including the endogenous baseline
    # (Equation 1).
    e0   <- exp(le0 + etale0)
    emax <- exp(lemax + etalemax) * ws^e_wt_emax *
      e_dis_hemophilia_a_mild_emax^DIS_HEMOPHILIA_A_MILD
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill)

    fviii_pd          <- fviii + fviii_bl
    edrug             <- emax * fviii_pd^hill / (ec50^hill + fviii_pd^hill)
    thrombinPotential <- e0 + edrug
    thrombinPotential ~ add(addSd_thrombinPotential)
  })
}
