Valke_2024_factorviii_thrombinPotential <- function() {
  description <- paste(
    "Sequential PK-PD model relating factor VIII (FVIII) activity to",
    "normalized thrombin potential (% of normal pooled plasma) measured with",
    "the Nijmegen hemostasis assay, after a single intravenous bolus of",
    "plasma-derived von Willebrand factor / factor VIII concentrate",
    "(pdVWF/FVIII, Humate-P) in 29 adults with severe hemophilia A",
    "(Valke 2024). The PK layer is the two-compartment allometric model of",
    "Valke 2024 Table 2 (see modellib('Valke_2024_factorviii')); the PD layer",
    "is an additive-baseline Emax model, E = Ebase + Emax * C^n / (EC50^n +",
    "C^n), with baseline 21.9% of normal pooled plasma, Emax 65.3% of normal",
    "pooled plasma, EC50 1.93 IU/dL and the Hill coefficient fixed at 1. The",
    "very low EC50 means thrombin potential stays near normal for 24 h after",
    "a bolus even as FVIII activity falls. Inter-individual variability is",
    "carried on Emax (33.1%); residual error on the PD readout is additive",
    "(16.2% of normal pooled plasma).",
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
  units <- list(
    time          = "h",
    dosing        = "IU",
    concentration = "IU/dL (observation Cc is FVIII activity; the PD output thrombinPotential is normalized thrombin potential as a percentage of normal pooled plasma)"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "A priori allometric power scaling of the PK parameters at a 70 kg reference with exponents fixed at 0.75 for CL and Q and 1 for V1 and V2 (Valke 2024 Supplementary Methods Eq. 1). Cohort mean 62 kg (SD 7), median 62 kg (range 48-73 kg; Table 1). Body weight was tested as a covariate on the PD parameters and was not significant in this cohort (Valke 2024 Results 3.5), so it acts on the PK layer only. Note that the earlier Bukkems 2022 model did retain a body-weight coefficient of -0.28 on the thrombin potential Emax (Valke 2024 Table 3); that term is not part of the re-estimated model encoded here.",
      source_name        = "Bodyweight"
    ),
    ADA_POS = list(
      description        = "FVIII inhibitor status by the Nijmegen-modified Bethesda assay (NBA): 1 = inhibitor-positive (titer >= 0.60 NBU/mL), 0 = negative. Mapped onto the canonical ADA_POS column per the NAB-subset alias documented in inst/references/covariate-columns.md.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (NBA inhibitor negative)",
      notes              = "Valke 2024 Table 2: the only covariate retained in the final PK model; CL is 1.53x higher in NBA-positive patients (theta^flag form of Supplementary Methods Eq. 6). 1 of 29 patients (3%) was NBA-positive (titer 1.1 NBU/mL). Valke 2024 Supplementary Figure 12 shows the thrombin potential response of this patient. No covariate was retained on the PD parameters.",
      source_name        = "NBA"
    ),
    ASSAY_OSA = list(
      description        = "One-stage activated partial thromboplastin time clotting assay (OSA) indicator (1 = OSA, 0 = chromogenic substrate assay, CSA)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CSA) per the canonical column definition; in THIS model the structural PK parameters are on the OSA scale, so the assay correction and the alternative residual-error magnitudes are applied to the CSA (ASSAY_OSA = 0) rows",
      notes              = "Valke 2024 Results 3.4 and Table 2: FVIII activity was measured by both assays and CSA samples read 0.939x the OSA value; residual error also differs by assay. Affects the FVIII activity observation Cc only -- the PD layer is driven by the untransformed one-stage-scale FVIII activity, which is what the individual PK predictions supplied to the sequential PD estimation step represent. Per-observation (per-row). Set ASSAY_OSA = 1 to simulate the one-stage-assay readout.",
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
    sex_female_pct = 0,
    disease_state  = "Severe hemophilia A (FVIII activity level < 1 IU/dL). Baseline thrombin potential 280 nM-min (undetectable-574), i.e. 12% of normal pooled plasma (undetectable-29), against a healthy-control reference of 1580 (SD 199) nM-min.",
    dose_range     = "Single intravenous bolus of pdVWF/FVIII concentrate (Humate-P, CSL Behring): median 1600 IU FVIII (IQR 1500-1700), i.e. 25.0 IU/kg (IQR 24.6-25.4).",
    regions        = "Iran (Ghaem Hospital, Mashhad University of Medical Sciences); samples analysed at Radboud University Medical Center, Nijmegen, The Netherlands.",
    notes          = "Sub-study of the IMPALA study (Dutch Trial Register NL2808). All patients male. 285 normalized thrombin potential values were available for the PD analysis. The analysis was sequential: individual PK parameters from the Valke 2024 PK model were fixed as the input to the PD estimation step. This is a replication study of Bukkems 2022 (Br J Clin Pharmacol 88(6):2757-2768); the Bukkems thrombin potential model already predicted this external cohort adequately (MPE 7.46%, MAPE 18.6%), and the re-estimated parameters encoded here are the paper's own fit to the replication cohort. The Bukkems estimates reproduced in Valke 2024 Table 3 for comparison are NOT encoded here."
  )

  ini({
    # ---- PK layer: Valke 2024 Table 2, "External dataset" column ----
    lcl <- log(3.07) ; label("Clearance CL at the reference 70 kg (dL/h)")                          # Valke 2024 Table 2: CL = 3.07 dL/h/70 kg (RSE 10%)
    lvc <- log(39.1) ; label("Central volume V1 at the reference 70 kg (dL)")                       # Valke 2024 Table 2: V1 = 39.1 dL/70 kg (RSE 7.7%)
    lq  <- log(1.09) ; label("Inter-compartmental clearance Q at the reference 70 kg (dL/h)")       # Valke 2024 Table 2: Q = 1.09 dL/h/70 kg (RSE 35.5%)
    lvp <- log(9.16) ; label("Peripheral volume V2 at the reference 70 kg (dL)")                    # Valke 2024 Table 2: V2 = 9.16 dL/70 kg (RSE 32.3%)

    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of (WT / 70) on CL and Q (unitless)")       # Valke 2024 Supplementary Methods Eq. 1: exponent fixed at 0.75 for CL and Q
    e_wt_vc <- fixed(1.00) ; label("Allometric exponent of (WT / 70) on V1 and V2 (unitless)")      # Valke 2024 Supplementary Methods Eq. 1: exponent fixed at 1 for V1 and V2

    e_ada_pos_cl <- 1.53 ; label("Multiplicative factor on CL for NBA inhibitor-positive patients (unitless)")   # Valke 2024 Table 2: Positive NBA on CL (%) = 153

    e_assay_osa_cc     <- 0.939 ; label("Multiplicative factor on predicted FVIII activity for CSA- vs OSA-assayed samples (unitless)")     # Valke 2024 Table 2: Correction factor CSA = 0.939 (RSE 2.2%)
    e_assay_osa_propsd <- 0.840 ; label("Multiplicative factor on the proportional residual SD for CSA- vs OSA-assayed samples (unitless)") # Valke 2024 Table 2: proportional error CSA 21.0% / OSA 25.0% = 0.840
    e_assay_osa_addsd  <- 5.012 ; label("Multiplicative factor on the additive residual SD for CSA- vs OSA-assayed samples (unitless)")     # Valke 2024 Table 2: additive error CSA 4.28 / OSA 0.854 IU/dL = 5.012

    # ---- PD layer: Valke 2024 Table 3, "External dataset" column ----
    # Normalized thrombin potential, additive-baseline Emax model
    # (Valke 2024 Supplementary Methods Eqs. 8 + 10 and the formula block
    # printed beneath Table 3).
    le0    <- log(21.9) ; label("Baseline normalized thrombin potential (% of normal pooled plasma)")            # Valke 2024 Table 3: Baseline effect = 21.9% of NPP (RSE 13.6%)
    lemax  <- log(65.3) ; label("Maximal FVIII effect on normalized thrombin potential (% of normal pooled plasma)")  # Valke 2024 Table 3: Maximal effect (Emax) = 65.3% of NPP (RSE 6.7%)
    lec50  <- log(1.93) ; label("FVIII activity producing half the maximal thrombin potential effect (IU/dL)")   # Valke 2024 Table 3: EC50 = 1.93 IU/dL (RSE 46.4%)
    lhill  <- fixed(log(1)) ; label("Hill coefficient of the thrombin potential Emax relationship (unitless)")   # Valke 2024 Table 3: Hill coefficient = 1 FIX

    # ---- Random effects ----
    # PK: correlated exponential IIV on CL and V1 (Supplementary Eq. 2).
    # Reported IIV percentages are taken as omega = CV, so omega^2 = CV^2,
    # following the convention of the related FVIII popPK models already in
    # nlmixr2lib. Off-diagonal = 0.770 * sqrt(0.328329 * 0.145924).
    etalcl + etalvc ~ c(0.328329,
                        0.168542, 0.145924)   # Valke 2024 Table 2: IIV CL = 57.3%, IIV V1 = 38.2%, correlation IIV CL-V1 = 77.0%

    # PD: IIV was added to a single PD parameter only -- the one with the
    # largest drop in objective function value (Valke 2024 Results 3.5). For
    # the thrombin potential model that parameter is Emax.
    etalemax ~ 0.109561   # Valke 2024 Table 3: IIV on the maximal effect = 33.1% (RSE 37.0%) [Shr 4%]; omega^2 = 0.331^2

    # ---- Residual error ----
    propSd <- 0.250 ; label("Proportional residual SD for OSA-assayed FVIII activity (fraction)")   # Valke 2024 Table 2: Proportional error OSA = 25.0% (RSE 6.6%)
    addSd  <- 0.854 ; label("Additive residual SD for OSA-assayed FVIII activity (IU/dL)")          # Valke 2024 Table 2: Additive error OSA = 0.854 IU/dL (RSE 26.5%)

    addSd_thrombinPotential <- 16.2 ; label("Additive residual SD for normalized thrombin potential (% of normal pooled plasma)")   # Valke 2024 Table 3: Additive error = 16.2% of NPP (RSE 6.7%)
  })

  model({
    # ---- PK ----
    ws         <- WT / 70
    inh_eff_cl <- e_ada_pos_cl ^ ADA_POS

    cl <- exp(lcl + etalcl) * ws ^ e_wt_cl * inh_eff_cl
    vc <- exp(lvc + etalvc) * ws ^ e_wt_vc
    q  <- exp(lq)           * ws ^ e_wt_cl
    vp <- exp(lvp)          * ws ^ e_wt_vc

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # FVIII activity on the one-stage-assay (OSA) scale (IU/dL). This is the
    # exogenous increment above the endogenous baseline, which was subtracted
    # from the observations during model development.
    fviii <- central / vc

    assay_csa <- 1 - ASSAY_OSA
    Cc <- fviii * e_assay_osa_cc ^ assay_csa

    propSd_eff <- propSd * e_assay_osa_propsd ^ assay_csa
    addSd_eff  <- addSd  * e_assay_osa_addsd  ^ assay_csa
    Cc ~ add(addSd_eff) + prop(propSd_eff)

    # ---- PD: normalized thrombin potential (% of normal pooled plasma) ----
    # Valke 2024, formula block beneath Table 3:
    #   E = Ebase + Emax * C^n / (EC50^n + C^n)
    # i.e. the additive-baseline relation of Supplementary Eq. 10 wrapped
    # around the sigmoid Emax drug effect of Supplementary Eq. 8. The
    # asymptote is Ebase + Emax = 21.9 + 65.3 = 87.2% of normal pooled plasma.
    # Supplementary Eq. 10 is printed as "E = Ebase + (1 + Edrug)"; the bare
    # "1 +" is a typographical carry-over from the proportional form of
    # Supplementary Eq. 9 and is not encoded here (see vignette Assumptions
    # and deviations).
    e0   <- exp(le0)
    emax <- exp(lemax + etalemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    edrug              <- emax * fviii ^ hill / (ec50 ^ hill + fviii ^ hill)
    thrombinPotential  <- e0 + edrug
    thrombinPotential ~ add(addSd_thrombinPotential)
  })
}
