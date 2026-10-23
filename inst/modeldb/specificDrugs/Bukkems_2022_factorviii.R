Bukkems_2022_factorviii <- function() {
  description <- paste(
    "Two-compartment population PK model for factor VIII (FVIII) activity",
    "(IU/dL) after an intravenous bolus of standard half-life FVIII",
    "concentrates (plasma-derived, full-length recombinant or B-domain-deleted",
    "recombinant) in 30 children and adults with severe, moderate or mild",
    "hemophilia A (Bukkems 2022). CL and V1 are allometrically scaled to a",
    "75 kg reference (exponents 0.75 and 1); Q and V2 are not scaled. Typical",
    "values are CL 1.69 dL/h, V1 27.7 dL, Q 2.27 dL/h and V2 5.63 dL.",
    "Clearance falls with von Willebrand factor activity (power -0.52 at a",
    "100% reference) and is 1.27x higher for full-length recombinant products",
    "and 1.49x higher in patients with a positive Nijmegen low-titre inhibitor",
    "assay; V1 is 1.17x and 1.14x higher for the same two factors. The",
    "patient's lowest-ever measured FVIII activity is added to the prediction",
    "as an endogenous baseline. Chromogenic-assay (CSA) samples read 1.20x the",
    "one-stage-assay (OSA) value, with inter-individual variability (18.1%) on",
    "that correction factor. IIV on CL (41.2%) and V1 (15.6%) is correlated",
    "at 43.6%. Residual error is combined proportional plus additive, with",
    "assay-specific magnitudes.",
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
  units <- list(time = "h", dosing = "IU", concentration = "IU/dL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric power scaling at a 75 kg reference, exponent 0.75 on CL and 1 on V1; Q and V2 are not scaled (Bukkems 2022 Table 2 footer, final PK equations). Cohort median 75 kg (range 11-134 kg; Table 1). Supplementary Methods Eq. 3 describes a priori scaling of all parameters at 70 kg; the printed final equations and the 75 kg typical patient of Supplementary Figure S2 are followed here (see vignette Assumptions and deviations).",
      source_name = "Weight"
    ),
    VWF = list(
      description = "Von Willebrand factor ristocetin cofactor activity (VWF:Act)",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (VWF / 100)^-0.52 on CL (Table 2 footer; 'VWF exponent on CL' -0.52, bootstrap 95% CI -0.79 to -0.20). Measured at baseline of each PK profile; cohort median 100% (range 28-216%; Table 1). VWF was missing for 14 of 44 PK profiles; the paper estimated a separate parameter for the missing group (Supplementary Methods Eq. 6) but does not report its value, so VWF must be supplied for every simulated patient. Percent of normal pooled plasma, numerically equal to IU/dL.",
      source_name = "VWF:Act"
    ),
    FORM_FVIII_RECOMB_FULLLENGTH = list(
      description = "Full-length recombinant FVIII product indicator: 1 = full-length recombinant concentrate (Advate, Helixate, Kogenate), 0 = plasma-derived (Aafact, Haemate P / Humate P) or B-domain-deleted recombinant (NovoEight, Refacto) concentrate",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (plasma-derived or B-domain-deleted recombinant product)",
      notes = "Multiplicative theta^flag effects 1.27 on CL and 1.17 on V1 (Table 2 rows 'Full-length recombinant product on CL (%) 127' and '... on V1 (%) 117'; Results 3.2: 'a 27% higher clearance and a 17% higher central volume of distribution as compared to plasma-derived and B-domain deleted recombinant concentrates'). 30 of 44 PK profiles used a full-length recombinant product (Table 1). Per-profile (a patient can switch product between profiles).",
      source_name = "recombinant product"
    ),
    ADA_POS = list(
      description = "FVIII inhibitor status by the Nijmegen low-titre inhibitor assay (NLTIA): 1 = positive, 0 = negative. FVIII inhibitors are neutralising alloantibodies to the administered FVIII; mapped onto the canonical ADA_POS column per the neutralising-antibody alias in inst/references/covariate-columns.md.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (NLTIA negative)",
      notes = "Multiplicative theta^flag effects 1.49 on CL and 1.14 on V1 (Table 2 rows 'Positive NLTIA on CL (%) 149' and '... on V1 (%) 114'). 7 of 30 patients (23%) were NLTIA-positive, titres 0.01-0.15 NLTIU/mL (Table 1). The NLTIA described the inhibitor effect better than the Nijmegen Bethesda assay (Results 3.2). NLTIA was missing for 9 of 44 PK profiles; no missing-group value is reported. Per-profile.",
      source_name = "Positive NLTIA"
    ),
    ASSAY_OSA = list(
      description = "One-stage activated partial thromboplastin time clotting assay (OSA) indicator (1 = OSA, 0 = chromogenic substrate assay, CSA)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CSA) per the canonical column definition; in THIS model the structural PK parameters are on the OSA scale, so the assay correction and the CSA residual-error magnitudes are applied to the CSA (ASSAY_OSA = 0) rows",
      notes = "Table 2 footer: 'CSA activity = OSA activity * 1.20' (Correction factor CSA 1.20, RSE 3.5%, IIV 18.1%). Both assays were fitted simultaneously: 466 OSA and 386 CSA observations (Results 3.1). Residual error differs by assay (proportional 11.2% / additive 4.15 IU/dL for OSA; 10.5% / 4.28 IU/dL for CSA). Per-observation (per-row). Set ASSAY_OSA = 1 to simulate the one-stage readout.",
      source_name = "OSA / CSA"
    ),
    FVIII_BL = list(
      description = "Endogenous baseline FVIII activity: the lowest FVIII activity ever measured in the patient",
      units = "IU/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Results 3.1: 'the lowest FVIII activity level ever measured was regarded as the endogenous baseline and subtracted from the observed FVIII concentration in the modelling process'; Supplementary Methods Eq. 2: C_obs = (C_pred + baseline FVIII_i) * (1 + eps_prop) + eps_add. Supplied in IU/mL per the canonical column and converted x100 to the model's IU/dL inside model(). Severity bands (Methods 2.1): severe < 0.01, moderate 0.01-0.05, mild > 0.05-0.40 IU/mL. Time-fixed per patient. Set to 0 to simulate the drug-attributable FVIII only.",
      source_name = "baseline FVIII"
    )
  )

  compartmentData <- list(
    central = list(analyte = "factor VIII", units = "IU", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "factor VIII", units = "IU", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 30L,
    n_studies = 1L,
    age_range = "2-77 years",
    age_median = "40 years",
    weight_range = "11-134 kg",
    weight_median = "75 kg",
    sex_female_pct = NULL,
    race_ethnicity = NULL,
    disease_state = "Hemophilia A: 20 severe (FVIII < 1 IU/dL), 4 moderate (1-5 IU/dL) and 6 mild (> 5-40 IU/dL) patients (Table 1). 7 (23%) NLTIA-positive and 4 (13%) Nijmegen-Bethesda-assay-positive low-titre inhibitors. Median VWF activity 100% (range 28-216).",
    dose_range = "Intravenous bolus of a standard half-life FVIII concentrate, median 29.0 IU/kg (range 10.9-90.9). 44 PK profiles (1-4 per patient): Aafact 8, Advate 17, Haemate P / Humate P 3, Helixate 8, Kogenate 5, NovoEight 2, Refacto 1.",
    regions = "The Netherlands (Radboud University Medical Center, Nijmegen; single centre).",
    notes = "Data collected 1 February 2003 to 1 July 2013 (Dutch Trial Registry NL2808). Five children < 12 years, none aged 12-18. Sex not reported. Samples pre-dose and at 3, 5, 15, 30 min and 1, 3, 6, 9 and 24 h after the bolus. 466 OSA and 386 CSA FVIII activity values; 5.8% (OSA) and 3.4% (CSA) were below the 1 IU/dL quantification limit and excluded. Where no wash-out preceded the profile, the last three prophylactic doses were included in the dataset to account for residual FVIII. This PK model supplied the individual PK parameters for the three sequential PD models (Bukkems_2022_factorviii_thrombinPeak, _thrombinPotential, _plasminPeak)."
  )

  ini({
    # Structural PK parameters at the 75 kg reference (Table 2; volumes in dL
    # and clearances in dL/h, so central / vc is directly in IU/dL).
    lcl <- log(1.69) ; label("Clearance CL at the reference 75 kg (dL/h)")                 # Table 2: CL = 1.69 dL/h (RSE 11.1%)
    lvc <- log(27.7) ; label("Central volume V1 at the reference 75 kg (dL)")              # Table 2: V1 = 27.7 dL (RSE 5.8%)
    lq  <- log(2.27) ; label("Inter-compartmental clearance Q (dL/h)")                     # Table 2: Q = 2.27 dL/h (RSE 44.5%)
    lvp <- log(5.63) ; label("Peripheral volume V2 (dL)")                                  # Table 2: V2 = 5.63 dL (RSE 19.2%)

    # A priori allometric exponents (Supplementary Methods Eq. 3).
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of (WT / 75) on CL (unitless)")   # Table 2 footer: CL = theta_CL * (Weight / 75)^0.75 ...; Supplementary Eq. 3 fixes 0.75 for clearances
    e_wt_vc <- fixed(1)    ; label("Allometric exponent of (WT / 75) on V1 (unitless)")   # Table 2 footer: V1 = theta_V1 * (Weight / 75) ...; Supplementary Eq. 3 fixes 1 for volumes

    # Covariate effects (Table 2 footer, final PK equations).
    e_vwf_cl <- -0.52 ; label("Power exponent of (VWF / 100) on CL (unitless)")            # Table 2: VWF exponent on CL = -0.52 (RSE 26.6%); bootstrap 95% CI -0.79 to -0.20
    e_form_fviii_recomb_fulllength_cl <- 1.27 ; label("Multiplicative factor on CL for full-length recombinant products (unitless)")   # Table 2: Full-length recombinant product on CL (%) = 127 (RSE 10.3%)
    e_form_fviii_recomb_fulllength_vc <- 1.17 ; label("Multiplicative factor on V1 for full-length recombinant products (unitless)")   # Table 2: Full-length recombinant product on V1 (%) = 117 (RSE 6.8%)
    e_ada_pos_cl <- 1.49 ; label("Multiplicative factor on CL for NLTIA inhibitor-positive patients (unitless)")   # Table 2: Positive NLTIA on CL (%) = 149 (RSE 11.1%)
    e_ada_pos_vc <- 1.14 ; label("Multiplicative factor on V1 for NLTIA inhibitor-positive patients (unitless)")   # Table 2: Positive NLTIA on V1 (%) = 114 (RSE 3.9%)

    # Assay effects, applied to chromogenic-assay rows (ASSAY_OSA = 0).
    e_assay_osa_cc     <- 1.20   ; label("Typical multiplicative factor on predicted FVIII activity for CSA- vs OSA-assayed samples (unitless)")   # Table 2: Correction factor CSA = 1.20 (RSE 3.5%); footer: CSA activity = OSA activity * 1.20
    e_assay_osa_propsd <- 0.9375 ; label("Multiplicative factor on the proportional residual SD for CSA- vs OSA-assayed samples (unitless)") # Table 2: proportional error CSA 10.5% / OSA 11.2% = 0.9375
    e_assay_osa_addsd  <- 1.0313 ; label("Multiplicative factor on the additive residual SD for CSA- vs OSA-assayed samples (unitless)")     # Table 2: additive error CSA 4.28 / OSA 4.15 IU/dL = 1.0313

    # IIV reported as CV = sqrt(exp(omega^2) - 1) (Table 2 abbreviations), so
    # omega^2 = log(CV^2 + 1). Off-diagonal = 0.436 * sqrt(0.156785 * 0.024045).
    etalcl + etalvc ~ c(0.156785,
                        0.026770, 0.024045)   # Table 2: IIV CL = 41.2% [Shr 2.1%], IIV V1 = 15.6% [Shr 2.8%], correlation IIV V1 and CL = 43.6%
    etae_assay_osa_cc ~ 0.032236              # Table 2: IIV correction factor CSA = 18.1% [Shr 4.7%]; log(1 + 0.181^2)

    # Combined residual error; OSA (reference) magnitudes, with the CSA
    # magnitudes reached through the e_assay_osa_* multipliers.
    propSd <- 0.112 ; label("Proportional residual SD for OSA-assayed FVIII activity (fraction)")   # Table 2: Proportional error OSA = 11.2% (RSE 21.6%)
    addSd  <- 4.15  ; label("Additive residual SD for OSA-assayed FVIII activity (IU/dL)")          # Table 2: Additive error OSA = 4.15 IU/dL (RSE 14.9%)
  })

  model({
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

    # Exogenous FVIII activity on the one-stage-assay scale (IU/dL); the
    # endogenous baseline (FVIII_BL, IU/mL) is converted x100 to IU/dL.
    fviii    <- central / vc
    fviii_bl <- 100 * FVIII_BL

    # Chromogenic-assay rows: the 1.20 factor (with IIV) scales the exogenous
    # prediction, to which the measured endogenous baseline is then added.
    assay_csa <- 1 - ASSAY_OSA
    cf_csa    <- e_assay_osa_cc * exp(etae_assay_osa_cc)
    Cc <- fviii * cf_csa^assay_csa + fviii_bl

    propSd_eff <- propSd * e_assay_osa_propsd^assay_csa
    addSd_eff  <- addSd * e_assay_osa_addsd^assay_csa
    Cc ~ add(addSd_eff) + prop(propSd_eff)
  })
}
