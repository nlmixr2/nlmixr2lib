Valke_2024_factorviii <- function() {
  description <- paste(
    "Two-compartment population PK model for factor VIII (FVIII) activity",
    "(IU/dL) after a single intravenous bolus of plasma-derived von",
    "Willebrand factor / factor VIII concentrate (pdVWF/FVIII, Humate-P) in",
    "29 adults with severe hemophilia A (Valke 2024). The model was newly",
    "estimated after the earlier Bukkems 2022 popPK model over-predicted",
    "FVIII activity in this cohort (MPE 60.9%). CL, V1, Q and V2 are",
    "allometrically scaled to a 70 kg reference with exponents fixed at 0.75",
    "on clearances and 1 on volumes; typical values are CL 3.07 dL/h,",
    "V1 39.1 dL, Q 1.09 dL/h, V2 9.16 dL, giving a terminal half-life of",
    "12.7 h. Clearance is 1.53x higher in patients with a Nijmegen-modified",
    "Bethesda assay inhibitor. FVIII activity was measured by both the",
    "one-stage clotting assay (OSA, the reference scale of this model) and",
    "the chromogenic substrate assay (CSA); CSA samples read 0.939x the OSA",
    "value and carry their own residual-error magnitudes. IIV on CL (57.3%)",
    "and V1 (38.2%) is correlated at 77.0%. Residual error is combined",
    "proportional plus additive, with assay-specific magnitudes.",
    sep = " "
  )
  reference <- paste(
    "Valke LLFG, Cloesmeijer ME, Mansouritorghabeh H, Barteling W,",
    "Blijlevens NMA, Cnossen MH, Mathot RAA, Schols SEM, van Heerde WL.",
    "Pharmacokinetic-Pharmacodynamic Modelling in Hemophilia A: Relating",
    "Thrombin and Plasmin Generation to Factor VIII Activity After",
    "Administration of a VWF/FVIII Concentrate.",
    "Eur J Drug Metab Pharmacokinet. 2024 Mar;49(2):191-205.",
    "doi:10.1007/s13318-024-00876-6. PMID:38367174. PMCID:PMC10904421.",
    sep = " "
  )
  vignette <- "Valke_2024_factorviii_thrombin_plasmin"
  units <- list(time = "h", dosing = "IU", concentration = "IU/dL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "A priori allometric power scaling at a 70 kg reference with exponents fixed at 0.75 for CL and Q and 1 for V1 and V2 (Valke 2024 Supplementary Methods Eq. 1: theta_popPK = theta_pk * (Bodyweight / 70)^theta_exp). Cohort mean 62 kg (SD 7), median 62 kg (range 48-73 kg; Valke 2024 Table 1). Subject-level, time-fixed in this single-occasion study.",
      source_name        = "Bodyweight"
    ),
    ADA_POS = list(
      description        = "FVIII inhibitor status by the Nijmegen-modified Bethesda assay (NBA): 1 = inhibitor-positive (titer >= 0.60 NBU/mL), 0 = negative. FVIII inhibitors are neutralizing alloantibodies to administered FVIII; mapped onto the canonical ADA_POS column per the NAB-subset alias documented in inst/references/covariate-columns.md.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (NBA inhibitor negative)",
      notes              = "Valke 2024 Table 2 and Results 3.4: the only covariate retained in the final PK model. 1 of 29 patients (3%) was NBA-positive (titer 1.1 NBU/mL). A separate, more sensitive low-titer assay (Nijmegen low-titer inhibitor assay, NLTIA, cut-off 0.04 NLTIU/mL) was positive in 7 of 29 patients (24%) but had no significant effect on CL in this cohort and is therefore not carried as a covariate here (in the earlier Bukkems 2022 model, NLTIA positivity acted on both CL and V1). Subject-level, time-fixed. The effect is encoded in the multiplicative theta^flag form of Valke 2024 Supplementary Methods Eq. 6 with theta = 1.53 (Table 2 row 'Positive NBA on CL (%) 153'); see the vignette Assumptions and deviations for why the table's percentage is read as a factor of the reference rather than as a percentage increase.",
      source_name        = "NBA"
    ),
    ASSAY_OSA = list(
      description        = "One-stage activated partial thromboplastin time clotting assay (OSA) indicator (1 = OSA, 0 = chromogenic substrate assay, CSA)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CSA) per the canonical column definition; note that in THIS model the structural PK parameters are on the OSA scale, so the assay correction and the alternative residual-error magnitudes are applied to the CSA (ASSAY_OSA = 0) rows",
      notes              = "Valke 2024 Results 3.4: every plasma sample was assayed by both methods (258 OSA + 258 CSA observations from 29 patients) and both were fitted simultaneously. 'Samples measured with CSA were 0.939 times lower compared to samples measured with OSA' (Table 2: Correction factor CSA = 0.939, RSE 2.2%), so the model's typical FVIII activity is on the OSA scale and CSA-assayed rows are multiplied by 0.939. Residual error also differs by assay (Table 2: proportional 25.0% / additive 0.854 IU/dL for OSA; proportional 21.0% / additive 4.28 IU/dL for CSA). IIV on the correction factor was tested but did not significantly improve the fit. Per-observation (per-row) indicator. Set ASSAY_OSA = 1 to simulate the one-stage-assay readout.",
      source_name        = "Assay method"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "factor VIII", units = "IU", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "factor VIII", units = "IU", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 29L,
    n_studies      = 1L,
    age_range      = "19-53 years",
    age_median     = "27 years",
    weight_range   = "48-73 kg",
    weight_median  = "62 kg",
    weight_mean_sd = "62 (SD 7) kg",
    sex_female_pct = 0,
    race_ethnicity = NULL,
    disease_state  = "Severe hemophilia A (FVIII activity level < 1 IU/dL). Median pre-bolus FVIII activity < 1 IU/dL by both assays; 8 of 29 patients had a detectable pre-bolus CSA level (1-2 IU/dL) and 4 had a detectable OSA level. Mean baseline VWF activity 117% (SD 46).",
    inhibitor_status = "1 of 29 (3%) positive by the Nijmegen-modified Bethesda assay (titer 1.1 NBU/mL); 7 of 29 (24%) positive by the more sensitive Nijmegen low-titer inhibitor assay (titers 0.04-0.05 NLTIU/mL). All patients were retained in the analysis.",
    dose_range     = "Single intravenous bolus of pdVWF/FVIII concentrate (Humate-P, CSL Behring): median 1600 IU FVIII (IQR 1500-1700), i.e. 25.0 IU/kg (IQR 24.6-25.4). 72 h wash-out before dosing.",
    regions        = "Iran (Ghaem Hospital, Mashhad University of Medical Sciences); samples analysed at Radboud University Medical Center, Nijmegen, The Netherlands.",
    notes          = "Sub-study of the IMPALA study (Dutch Trial Register NL2808), enrolled 1 August 2011 to 20 December 2012. All patients male (hemophilia A is X-linked recessive); race / ethnicity not reported. Samples were drawn pre-bolus and at nine time points to 24 h. 258 OSA and 258 CSA FVIII activity observations entered the PK analysis; 5.2% (OSA) and 4.2% (CSA) of samples were below the assay detection limit -- mostly pre-dose -- and were excluded. The pre-bolus FVIII activity was treated as the endogenous baseline and SUBTRACTED from the observed activities during model development, so the model describes the exogenous FVIII activity increment only and carries no endogenous-baseline parameter. To anchor the terminal phase, the pre-bolus sample was additionally used as if it had been taken 72 h after dosing (the wash-out duration). Median observed FVIII half-life 10.6 h (IQR 8.3-12.9); the model's terminal half-life at 70 kg is 12.7 h. This is a replication / external-validation study of the Bukkems 2022 model (Br J Clin Pharmacol 88(6):2757-2768); the Bukkems parameter estimates reproduced in Valke 2024 Tables 2 and 3 for comparison are NOT encoded here."
  )

  ini({
    # Structural PK parameters at the reference body weight of 70 kg.
    # Valke 2024 Table 2, "External dataset" column. Volumes are reported in
    # dL and clearances in dL/h, so with dose in IU the observable
    # central / vc is directly in the paper's IU/dL unit.
    lcl <- log(3.07) ; label("Clearance CL at the reference 70 kg (dL/h)")                          # Valke 2024 Table 2: CL = 3.07 dL/h/70 kg (RSE 10%)
    lvc <- log(39.1) ; label("Central volume V1 at the reference 70 kg (dL)")                       # Valke 2024 Table 2: V1 = 39.1 dL/70 kg (RSE 7.7%)
    lq  <- log(1.09) ; label("Inter-compartmental clearance Q at the reference 70 kg (dL/h)")       # Valke 2024 Table 2: Q = 1.09 dL/h/70 kg (RSE 35.5%)
    lvp <- log(9.16) ; label("Peripheral volume V2 at the reference 70 kg (dL)")                    # Valke 2024 Table 2: V2 = 9.16 dL/70 kg (RSE 32.3%)

    # A priori allometric exponents, held fixed by the authors.
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of (WT / 70) on CL and Q (unitless)")       # Valke 2024 Supplementary Methods Eq. 1: "theta_exp is an exponent fixed at 0.75 for CL and Q"
    e_wt_vc <- fixed(1.00) ; label("Allometric exponent of (WT / 70) on V1 and V2 (unitless)")      # Valke 2024 Supplementary Methods Eq. 1: "... and 1 for V1 and V2"

    # Categorical covariate effect, multiplicative in the theta^flag form of
    # Valke 2024 Supplementary Methods Eq. 6.
    e_ada_pos_cl <- 1.53 ; label("Multiplicative factor on CL for NBA inhibitor-positive patients (unitless)")   # Valke 2024 Table 2: Positive NBA on CL (%) = 153

    # Bioanalytical assay effects, applied to CSA-assayed rows
    # (ASSAY_OSA = 0). Same theta^flag form (Supplementary Methods Eq. 6).
    e_assay_osa_cc     <- 0.939 ; label("Multiplicative factor on predicted FVIII activity for CSA- vs OSA-assayed samples (unitless)")   # Valke 2024 Table 2: Correction factor CSA = 0.939 (RSE 2.2%)
    e_assay_osa_propsd <- 0.840 ; label("Multiplicative factor on the proportional residual SD for CSA- vs OSA-assayed samples (unitless)") # Valke 2024 Table 2: proportional error CSA 21.0% / OSA 25.0% = 0.840
    e_assay_osa_addsd  <- 5.012 ; label("Multiplicative factor on the additive residual SD for CSA- vs OSA-assayed samples (unitless)")     # Valke 2024 Table 2: additive error CSA 4.28 / OSA 0.854 IU/dL = 5.012

    # Correlated exponential IIV on CL and V1 (Supplementary Methods Eq. 2:
    # theta_i = theta_popPK * exp(eta_i)). The reported IIV percentages are
    # taken as omega = CV, so omega^2 = CV^2, following the convention used
    # by the related FVIII popPK models already in nlmixr2lib
    # (Hazendonk_2016_factor_viii, Chelle_2019_factorviii_fanhdi,
    # Nestorov_2014_factorviii). Off-diagonal = 0.770 * sqrt(0.328329 * 0.145924).
    etalcl + etalvc ~ c(0.328329,
                        0.168542, 0.145924)   # Valke 2024 Table 2: IIV CL = 57.3% (RSE 15%) [Shr 1.4%], IIV V1 = 38.2% (RSE 15%) [Shr 0.1%], correlation IIV CL-V1 = 77.0%

    # Combined residual error. The values below are the OSA (reference)
    # magnitudes; the CSA magnitudes are reached through the
    # e_assay_osa_* multipliers in model().
    propSd <- 0.250 ; label("Proportional residual SD for OSA-assayed FVIII activity (fraction)")   # Valke 2024 Table 2: Proportional error OSA = 25.0% (RSE 6.6%)
    addSd  <- 0.854 ; label("Additive residual SD for OSA-assayed FVIII activity (IU/dL)")          # Valke 2024 Table 2: Additive error OSA = 0.854 IU/dL (RSE 26.5%)
  })

  model({
    # Allometric size scaling at the 70 kg reference (Supplementary Eq. 1).
    ws <- WT / 70

    # NBA inhibitor effect on CL, multiplicative (Supplementary Eq. 6).
    inh_eff_cl <- e_ada_pos_cl ^ ADA_POS

    # Individual PK parameters. IIV is carried on CL and V1 only; Q and V2
    # were estimated without IIV.
    cl <- exp(lcl + etalcl) * ws ^ e_wt_cl * inh_eff_cl
    vc <- exp(lvc + etalvc) * ws ^ e_wt_vc
    q  <- exp(lq)           * ws ^ e_wt_cl
    vp <- exp(lvp)          * ws ^ e_wt_vc

    # Micro-constants for the explicit two-compartment ODEs.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV disposition; the bolus enters the central
    # compartment directly.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # FVIII activity on the one-stage-assay (OSA) scale, in IU/dL: dose in IU
    # and vc in dL give central / vc directly in IU/dL. Note that the
    # observations were baseline-subtracted during model development, so this
    # is the exogenous FVIII increment above the patient's endogenous level.
    fviii <- central / vc

    # Chromogenic-assay (CSA) rows read 0.939x the one-stage value and carry
    # their own residual-error magnitudes. ASSAY_OSA = 1 leaves every
    # multiplier at 1; ASSAY_OSA = 0 applies each e_assay_osa_* factor once.
    assay_csa <- 1 - ASSAY_OSA
    Cc <- fviii * e_assay_osa_cc ^ assay_csa

    propSd_eff <- propSd * e_assay_osa_propsd ^ assay_csa
    addSd_eff  <- addSd  * e_assay_osa_addsd  ^ assay_csa
    Cc ~ add(addSd_eff) + prop(propSd_eff)
  })
}
