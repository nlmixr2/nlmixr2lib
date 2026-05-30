Hanberg_2018_meropenem <- function() {
  description <- paste(
    "Two-compartment IV population PK model for meropenem in critically ill",
    "adults receiving venovenous or venoarterial extracorporeal membrane",
    "oxygenation (ECMO) treatment, with simultaneous fitting of plasma",
    "concentrations (central compartment Ac/Vc) and free subcutaneous",
    "adipose-tissue (SCT) concentrations sampled by microdialysis (peripheral",
    "compartment Ap/Vp scaled by an estimated fraction unbound in tissue",
    "f_u,tissue = 0.79). Elimination clearance is a direct linear function of",
    "the patient's estimated creatinine clearance (eCLCr, Cockcroft-Gault,",
    "raw mL/min) via CL_i = CLfrac * eCLCr_i with CLfrac = 0.0460 L/h per",
    "(mL/min); 9 of 10 patients were also on continuous renal replacement",
    "therapy so eCLCr partly reflects the CRRT contribution (Hanberg 2018)."
  )
  reference <- paste(
    "Hanberg P, Obrink-Hansen K, Thorsted A, Bue M, Tottrup M, Friberg LE,",
    "Hardlei TF, Soballe K, Gjedsted J. (2018).",
    "Population pharmacokinetics of meropenem in plasma and subcutis from",
    "patients on extracorporeal membrane oxygenation treatment.",
    "Antimicrob Agents Chemother 62(5):e02390-17.",
    "doi:10.1128/AAC.02390-17"
  )
  vignette <- "Hanberg_2018_meropenem"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Estimated creatinine clearance by the Cockcroft-Gault formula",
        "(raw mL/min, NOT BSA-normalized)."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column eCLCr. Computed by the Cockcroft-Gault equation in",
        "raw mL/min (NOT BSA-normalized to mL/min/1.73 m^2). Stored under",
        "the canonical CRCL column per inst/references/covariate-columns.md",
        "(CRCL accepts raw mL/min when the source paper does not apply BSA",
        "normalization, with the per-model description recording the assay",
        "form; precedent: Delattre_2010_amikacin.R, Shekar_2014_meropenem.R).",
        "Patient population: median 59.5 mL/min, individual range 36.9-135.1",
        "(Table 1 and Table 3). All but 1 of 10 patients were simultaneously",
        "on continuous renal replacement therapy (CRRT), so eCLCr partly",
        "reflects the CRRT-modulated solute removal; Hanberg 2018 Discussion",
        "cautions that the proportionality between CL and eCLCr should not",
        "be directly extrapolated to settings without CRRT. Used as a direct",
        "linear covariate on CL: CL_i = e_crcl_cl * CRCL (units L/h, since",
        "the published slope e_crcl_cl = 0.0460 already has units L/h per",
        "(mL/min); no in-model unit conversion is needed)."
      ),
      source_name        = "eCLCr"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "30-69 years (median 56)",
    age_median     = "56 years",
    weight_range   = "55-134 kg (median 100.5)",
    weight_median  = "100.5 kg",
    sex_female_pct = 50,
    race_ethnicity = "Not reported (single-centre Danish national VV-ECMO referral hospital, presumed predominantly European)",
    disease_state  = paste(
      "Critically ill adults on venovenous or venoarterial ECMO treatment",
      "for severe heart and/or lung failure not responding to conventional",
      "treatment; underlying infections included influenza A virus pneumonia",
      "(n = 6), Streptococcus pneumoniae pneumonia (n = 3), Stenotrophomonas",
      "maltophilia, Chlamydia pneumoniae, Staphylococcus haemolyticus +",
      "enterococcus, group A streptococcus, and aspergillus galactomannan",
      "co-infection. SOFA scores ranged 4-15 (median 11.5) at inclusion."
    ),
    dose_range     = paste(
      "Meropenem 1 g (n = 7) or 2 g (n = 3) IV bolus infused over 5 min,",
      "every 8 hours. ECMO and meropenem treatment were initiated less than",
      "96 h prior to inclusion; mean of 7 (range 3-10) prior meropenem doses",
      "before the sampling on day 1, with additional samples on days 2, 4,",
      "and 6."
    ),
    regions        = paste(
      "Denmark, single-centre (Aarhus University Hospital, national",
      "VV-ECMO centre of Denmark, 25 annual VV-ECMO and 75 VA-ECMO runs).",
      "Study conducted February-May 2016."
    ),
    renal_function = paste(
      "Cockcroft-Gault eCLCr median 59.5 mL/min, individual range",
      "36.9-135.1 mL/min (Table 3). 9 of 10 patients were on continuous",
      "renal replacement therapy (CRRT) concurrently with meropenem;",
      "eCLCr partly reflects the CRRT contribution."
    ),
    ecmo_modalities = paste(
      "Mix of venovenous (VV) and venoarterial (VA) ECMO; centrifugal pump",
      "speed median 3,200 RPM (range 2,660-3,930), pump flow median 4.56 LPM",
      "(3.42-5.31), sweep gas median 4.75 LPM (3.0-7.5), heparin dose",
      "median 650 IE/h (200-2,000). All patients had ECMO treatment ongoing",
      "for a mean of 51.9 h (range 20-69 h) prior to inclusion."
    ),
    tissue_sampling = paste(
      "Subcutaneous adipose tissue (SCT) concentrations were obtained by",
      "microdialysis (MD) probe (CMA 63, 30 mm membrane, 20-kDa cutoff)",
      "inserted in the upper arm SCT, perfused at 2 uL/min with 0.9% NaCl",
      "containing 5 ug/mL cefuroxime as internal calibrator. Mean (SD)",
      "relative recovery was 18.5% (5.2)%, computed by retrodialysis by",
      "calibrator. Dialysates were sampled at 15-min, 30-min, and 60-min",
      "intervals over the day-1 dosing interval. MD data were modeled per",
      "Tunblad et al. 2008 as the dialysate concentration at the end of the",
      "collection interval (not averaged across the interval). The MD probe",
      "for one of the ten patients was malfunctioning and the dialysate",
      "data from that subject were excluded from the SCT-arm fit; plasma",
      "data from all 10 patients were included."
    ),
    notes          = paste(
      "Baseline demographics per Hanberg 2018 Table 1 (per-patient and",
      "medians). Quantification by UHPLC-UV (Agilent 1290 Infinity, C18",
      "column, 304 nm detection); LLOQ 0.5 ug/mL plasma, interrun",
      "imprecision 3.0% at 2.0 ug/mL. Plasma protein binding of meropenem",
      "is assumed negligible (< 2%) so the estimated fraction unbound in",
      "SCT (f_u,tissue = 0.79) can also be interpreted as the SCT-to-plasma",
      "AUC ratio (pAUC/cAUC = 0.79 across all patients)."
    )
  )

  ini({
    # Structural fixed-effect parameters from Hanberg 2018 Table 2
    # (Final-model column, with eCLCr as a covariate on CL).
    # Hanberg 2018 publishes CL as the linear relationship
    # CL_i = CLfrac * eCLCr_i with CLfrac = 0.0460 L/h per (mL/min); the
    # packaged model re-expresses the same equation as a typical CL at
    # the cohort-median eCLCr (lcl = log(CLfrac * 59.5) = log(2.737))
    # times an eCLCr-ratio centering term (CRCL / 59.5). This is the
    # standard nlmixr2lib "typical-value + covariate-ratio"
    # parameterization (matches the parameter-naming convention's
    # etalcl <-> lcl pairing rule) and is exactly equivalent to the
    # published slope form: 2.737 * (CRCL/59.5) = 0.0460 * CRCL. The
    # typical CL of 2.737 L/h matches the Discussion text
    # "clearance (2.74 liters/h for a patient with median eCLCr of
    # 59.5 ml/min)".
    lcl <- log(0.0460 * 59.5)
    label("Typical CL at cohort-median eCLCr 59.5 mL/min (L/h)")           # Hanberg 2018 Table 2 (final model) row "eCLCr (mL/min)" = 0.0460 (%RSE 6.7); typical CL = 0.0460 * 59.5 = 2.737 L/h (Discussion p. 7)
    lvc       <- log(8.31)
    label("Central volume Vc (L)")                                          # Hanberg 2018 Table 2 (final model) row "Vc" = 8.31 L, %RSE 9.0
    lq        <- log(8.52)
    label("Intercompartmental clearance Q (L/h)")                           # Hanberg 2018 Table 2 (final model) row "Q" = 8.52 L/h, %RSE 31
    lvp       <- log(6.99)
    label("Peripheral volume Vp (L)")                                       # Hanberg 2018 Table 2 (final model) row "Vp" = 6.99 L, %RSE 15

    # Fraction unbound in SCT (estimated). Encoded on logit scale per the
    # nlmixr2lib logit-fraction convention (logitfr / logitfu precedent in
    # Tsuji_2017_linezolid.R). The _sct suffix disambiguates from plasma
    # protein-binding fraction unbound (which the paper treats as negligible
    # at <2% for meropenem).
    logitfu_sct <- log(0.790 / (1 - 0.790))
    label("Logit-transformed fraction unbound in SCT (f_u,tissue = 0.790)") # Hanberg 2018 Table 2 (final model) row "fu" = 0.790, %RSE 4.6

    # Between-subject variability from Hanberg 2018 Table 2 (final-model
    # column, %CV with shrinkage in parentheses). omega^2 = log(CV^2 + 1)
    # for log-normal etas. No IIV on Vp or fu_sct per Table 2.
    etalcl ~ 0.03511 # log(0.189^2 + 1); Hanberg 2018 Table 2: BSV on CL = 18.9 %CV (SHR 0%)
    etalvc ~ 0.06112 # log(0.251^2 + 1); Hanberg 2018 Table 2: BSV on Vc = 25.1 %CV (SHR 7%)
    etalq  ~ 0.34853 # log(0.646^2 + 1); Hanberg 2018 Table 2: BSV on Q  = 64.6 %CV (SHR 6%)

    # Proportional residual error. Hanberg 2018 Table 2 reports a single
    # ERR row (19.9 %CV) applied to both plasma (Cc) and SCT (Csct)
    # observations. The nlmixr2 multi-endpoint syntax requires a distinct
    # parameter per output, so the single source value is replicated into
    # propSd (parent, applied to Cc) and propSd_Csct (suffix per the
    # multi-output convention, applied to Csct); both carry the same
    # initial value 0.199 to preserve the source's shared-residual
    # encoding. See the vignette Assumptions and deviations section for
    # the explicit note that these two parameters are operationally a
    # single estimated quantity in the source publication.
    propSd      <- 0.199
    label("Proportional residual error on Cc (fraction)")            # Hanberg 2018 Table 2 (final model) row "ERR" = 19.9 %CV (SHR 4%)
    propSd_Csct <- 0.199
    label("Proportional residual error on Csct (fraction)")          # Hanberg 2018 Table 2 (final model) row "ERR" = 19.9 %CV (SHR 4%); same shared estimate as propSd
  })
  model({
    # Individual PK parameters with log-normal IIV. CL is driven by the
    # individual eCLCr through the published linear relationship
    # CL_i = CLfrac * eCLCr_i (Hanberg 2018 Results, "covariate
    # analysis"); re-expressed equivalently as a typical CL at median
    # eCLCr times an eCLCr-ratio centering term.
    cl <- exp(lcl + etalcl) * (CRCL / 59.5)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp)

    # Back-transform fraction unbound in SCT from logit scale (bounds 0-1
    # respected without explicit constraint at simulation time).
    fu_sct <- 1 / (1 + exp(-logitfu_sct))

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment ODE system. Dose enters central as IV (5-min infusion
    # at simulation time via the rate column on dosing rows).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Observation variables.
    #   Cc   = plasma (central) total meropenem; protein binding is < 2% so
    #          this is operationally the free plasma concentration as well.
    #   Csct = subcutaneous-tissue (peripheral) free meropenem; the
    #          peripheral compartment carries total tissue drug, scaled by
    #          fu_sct to yield the free SCT concentration that MD samples.
    # Dose in mg, volumes in L -> Cc / Csct have units mg/L.
    Cc   <- central     / vc
    Csct <- peripheral1 / vp * fu_sct

    # Residual error -- per-endpoint proportional SD parameters seeded
    # from the single shared ERR row in Hanberg 2018 Table 2.
    Cc   ~ prop(propSd)
    Csct ~ prop(propSd_Csct)
  })
}
