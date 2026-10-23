Bukkems_2022_factorviii_thrombinPeak <- function() {
  description <- paste(
    "Sequential PK-PD model relating factor VIII (FVIII) activity to",
    "normalized thrombin peak height (% of normal pooled plasma) measured",
    "with the Nijmegen haemostasis assay, after an intravenous bolus of",
    "standard half-life FVIII concentrates in 24 children and adults with",
    "severe, moderate or mild hemophilia A (Bukkems 2022). The PK layer is",
    "the two-compartment allometric model of Bukkems 2022 (see",
    "modellib('Bukkems_2022_factorviii')); the PD layer is a",
    "proportional-baseline sigmoid Emax model, E = Ebase * (1 + Emax * C^n /",
    "(EC50^n + C^n)), driven by total one-stage FVIII activity (exogenous",
    "plus the patient's endogenous baseline), with baseline 15.6% of normal",
    "pooled plasma, Emax 7.05 times baseline, EC50 50.1 IU/dL and Hill",
    "coefficient 1.85. IIV on EC50 (55.1%) and Emax (37.3%); additive",
    "residual error 11.2% of normal pooled plasma.",
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
    concentration = "IU/dL (observation Cc is FVIII activity; the PD output thrombinPeak is normalized thrombin peak height as a percentage of normal pooled plasma)"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "PK layer only: allometric power scaling at a 75 kg reference, exponent 0.75 on CL and 1 on V1; Q and V2 not scaled (Bukkems 2022 Table 2 footer). PD cohort median 85 kg (range 11-134 kg; Table 1). No covariate was retained on the thrombin peak height PD parameters.",
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
      notes = "PK layer: multiplicative theta^flag effects 1.49 on CL and 1.14 on V1 (Table 2). 4 of 24 PD patients were NLTIA-positive. An NLTIA effect on the thrombin peak height EC50 (+30%) was found in forward inclusion but not retained after backward elimination (Discussion), so it is not encoded.",
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
    notes = "252 normalized thrombin peak height values from 24 of the 30 PK patients. Each NHA value is the mean of two measurements, normalized to the normal pooled plasma run alongside the patient sample. Sequential analysis: the individual empirical Bayes PK parameters from the Bukkems 2022 PK model were fixed as input to the PD estimation; the PK layer here is that population PK model, so simulations propagate PK variability as well."
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

    # ---- PD layer: normalized thrombin peak height (Table 2) ----
    le0   <- log(15.6) ; label("Baseline normalized thrombin peak height (% of normal pooled plasma)")              # Table 2: Baseline effect = 15.6% of NPP (RSE 18.8%)
    lemax <- log(7.05) ; label("Maximal FVIII effect on normalized thrombin peak height (factor of baseline)")     # Table 2: Maximal effect (factor of baseline) = 7.05 (RSE 33.6%)
    lec50 <- log(50.1) ; label("FVIII activity producing half the maximal thrombin peak height effect (IU/dL)")    # Table 2: EC50 = 50.1 IU/dL (RSE 24.4%)
    lhill <- log(1.85) ; label("Hill coefficient of the thrombin peak height Emax relationship (unitless)")        # Table 2: Hill coefficient = 1.85 (RSE 25.7%)

    # ---- Random effects: omega^2 = log(CV^2 + 1) (Table 2 abbreviations) ----
    etalcl + etalvc ~ c(0.156785,
                        0.026770, 0.024045)   # Table 2: IIV CL = 41.2%, IIV V1 = 15.6%, correlation IIV V1 and CL = 43.6%
    etae_assay_osa_cc ~ 0.032236              # Table 2: IIV correction factor CSA = 18.1%; log(1 + 0.181^2)
    etalec50 ~ 0.265130                       # Table 2: IIV EC50 (thrombin peak height) = 55.1% [Shr 12.5%]; log(1 + 0.551^2)
    etalemax ~ 0.130264                       # Table 2: IIV Emax (thrombin peak height) = 37.3% [Shr 16.8%]; log(1 + 0.373^2)

    # ---- Residual error ----
    propSd <- 0.112 ; label("Proportional residual SD for OSA-assayed FVIII activity (fraction)")   # Table 2: Proportional error OSA = 11.2%
    addSd  <- 4.15  ; label("Additive residual SD for OSA-assayed FVIII activity (IU/dL)")          # Table 2: Additive error OSA = 4.15 IU/dL

    addSd_thrombinPeak <- 11.2 ; label("Additive residual SD for normalized thrombin peak height (% of normal pooled plasma)")   # Table 2: Additive error = 11.2% of NPP (RSE 8.0%)
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

    # ---- PD: normalized thrombin peak height (% of normal pooled plasma) ----
    # Table 2 footer: E = Ebase * (1 + Emax * C^n / (EC50^n + C^n)), with C
    # the one-stage FVIII activity including the endogenous baseline
    # (Equation 1). Emax is a multiple of the baseline.
    e0   <- exp(le0)
    emax <- exp(lemax + etalemax)
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill)

    fviii_pd     <- fviii + fviii_bl
    edrug        <- emax * fviii_pd^hill / (ec50^hill + fviii_pd^hill)
    thrombinPeak <- e0 * (1 + edrug)
    thrombinPeak ~ add(addSd_thrombinPeak)
  })
}
