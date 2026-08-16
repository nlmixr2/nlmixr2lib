vandenBerg_2025_genericMab_mbma <- function() {
  description <- paste0(
    "MBMA. Generic ('one model fits all mAbs') two-compartment population PK model for ",
    "immunoglobulin G monoclonal antibodies with a subcutaneous absorption depot and linear ",
    "clearance, built from the medians of a systematic meta-analysis of 160 published ",
    "two-compartment population PK models covering 69 marketed canonical mAbs. Every ",
    "structural parameter is held at the across-model median (CL 0.22 L/day, Vc 3.42 L, ",
    "Vp 2.68 L, Q 0.54 L/day, ka 0.25 /day, F 0.69) and the inter-individual variances are ",
    "held at the across-model median IIV reported in the final models (CL 31.9 CV%, ",
    "Vc 24 CV%, ka 54.8 CV%); no covariates and no meaningful residual error were ",
    "introduced. Unlike a regression MBMA, the random effects here are genuine ",
    "between-SUBJECT variability (medians of the published between-subject omegas), not a ",
    "study-level effect, so the model simulates individual concentration-time profiles; the ",
    "separate between-MODEL variability reported by the paper (CL 54.6 CV%, Vc 25 CV%, ",
    "Vp 73.9 CV%, Q 108 CV%, ka 34 CV%, F 61 CV%) is not encoded. Intended as an a priori ",
    "/ prior-informed starting point for a new mAb PK model or for fixing poorly identifiable ",
    "parameters; it deliberately excludes target-mediated and time-dependent elimination. ",
    "Values are van den Berg 2025 Figure 5 caption and Supplementary Section M5 (NONMEM 7.5.0 ",
    "control stream, all THETA/OMEGA FIX)."
  )

  reference <- paste(
    "van den Berg SPH, Adolfsen PEA, Dorlo TPC, Rispens T.",
    "Does one model fit all mAbs? An evaluation of population pharmacokinetic models.",
    "mAbs. 2025;17(1):2512217.",
    "doi:10.1080/19420862.2025.2512217",
    sep = " "
  )

  vignette <- "vandenBerg_2025_genericMab"

  units <- list(
    time          = "day",
    dosing        = "mg (intravenous into central, or subcutaneous into depot; the paper's illustrative simulation used a single 100 mg dose by either route)",
    concentration = "mg/L (equivalently ug/mL; central amount in mg divided by Vc in L, matching the NONMEM scaling S2 = VC in Supplementary Section M5)"
  )

  # No covariates are part of this model -- the generic model was built
  # deliberately without them (Methods, "Generic model simulation": "No RUV or
  # covariates were introduced"). The covariates below WERE characterised across
  # the 160 collected models and their median effects are reported by the paper,
  # but the authors did not fold them into the generic model, so they are
  # documentation only. The two weight-based covariate configurations the authors
  # DID simulate live in the sibling models
  # vandenBerg_2025_genericMab_allometricClVc_mbma and
  # vandenBerg_2025_genericMab_allometricAll_mbma.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight; the most frequently used size descriptor (on CL in 66.9% and on Vc in 71.2% of the 160 models)",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Median (IQR) estimated allometric exponent 0.62 (0.52; 0.81) on CL and",
        "0.56 (0.47; 0.76) on Vc; median (IQR) reference weight 70 (69; 75) kg",
        "(van den Berg 2025 Results 'Covariate model' and Figure 4a).",
        "Not in the generic model; see the two sibling allometric models."
      ),
      source_name = "WT"
    ),
    ALBUMIN = list(
      description = "Serum albumin; negatively associated with CL in 41.9% of the 160 models",
      units       = "g/L",
      type        = "continuous",
      notes       = paste(
        "Median (IQR) power coefficient on CL -0.90 (-1.1; -0.69) with a median (IQR)",
        "reference of 40 (39; 42.8) g/L (van den Berg 2025 Results 'Covariate model',",
        "Figure 4c). Screened and summarised but not included in the generic model."
      ),
      source_name = "ALB"
    ),
    SEXF = list(
      description = "Sex; a significant covariate in 37.5% of the 160 models",
      units       = NA_character_,
      type        = "categorical",
      reference_category = "female (SEXF = 1)",
      notes       = paste(
        "Median (IQR) CL increase in males vs females 19% (15%; 26%); the effect shrinks",
        "from 1.30 to 1.18 once a size descriptor is also in the model",
        "(van den Berg 2025 Figure 4b,e; conversion to a female reference in",
        "Supplementary Section M2). Not included in the generic model."
      ),
      source_name = "SEX"
    ),
    ADA_POS = list(
      description = "Anti-drug antibody status; a significant covariate on CL in 26 of the 160 models (13 mAbs)",
      units       = NA_character_,
      type        = "categorical",
      reference_category = "ADA-negative (ADA_POS = 0)",
      notes       = paste(
        "Median (IQR) multiplicative coefficient on CL for ADA-positive patients",
        "1.39 (1.19; 1.49), averaged per mAb (van den Berg 2025 Results 'Covariate model',",
        "Figure 4d). Not included in the generic model."
      ),
      source_name = "ADA"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "genericMab", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "genericMab", units = "mg", specimen = "serum", verified = FALSE),
    peripheral1 = list(analyte = "genericMab", units = "mg", specimen = "serum", verified = FALSE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 143094L,
    n_studies      = 160L,
    n_models       = 160L,
    n_mabs         = 69L,
    disease_state  = paste(
      "Mixed. The pooled evidence base spans every indication in which a marketed canonical",
      "IgG mAb has had a population PK model published -- oncology, rheumatology,",
      "gastroenterology, dermatology, neurology, ophthalmology and infectious disease --",
      "plus healthy donors, who were included in 37 of the 160 models (26 of 69 mAbs).",
      "Paediatric-only models were excluded (van den Berg 2025 Figure 1a)."
    ),
    dose_range     = paste(
      "Not summarised across the source models. The paper's illustrative simulation used a",
      "single 100 mg intravenous or subcutaneous dose in 500 virtual individuals",
      "(van den Berg 2025 Methods 'Generic model simulation')."
    ),
    regions        = "Global; the collected models were not stratified by region.",
    design         = paste(
      "Meta-analysis of published population PK model PARAMETER ESTIMATES, not of raw or",
      "digitised concentration data. IMGT/mAb-DB was queried for marketed canonical IgG mAbs",
      "and PubMed searched per mAb; 303 candidate models were reduced to 160 by excluding",
      "one-compartment models, paediatric models, subcutaneous models with bioavailability",
      "fixed to 100%, models without linear clearance, models fixing more than two of",
      "CL/Vc/Vp/Q, PBPK and non-population models, non-i.v./s.c. routes, FcRn-engineered mAbs,",
      "and models with insufficient reporting (van den Berg 2025 Methods 'PK models",
      "acquisition'; Figure 1a; Supplementary Table S1)."
    ),
    reference_subject = paste(
      "A typical individual with no covariates. Parameter estimates were pooled AS REPORTED,",
      "without normalisation to a common reference covariate value; the authors confirmed that",
      "normalising the 114 weight-covariate models to 70 kg did not shift the medians",
      "(van den Berg 2025 Methods 'Data analysis'; Figure S1A; Supplementary Section M3)."
    ),
    composition = paste(
      "By model (n = 160): IgG1 71.9%, IgG2 8.8%, IgG4 18.1%, IgG2/4 1.2%; chimeric 20.0%,",
      "humanised 35.0%, fully human 45.0%; i.v. only 67.5%, s.c. only 0.6%, i.v./s.c. 31.9%;",
      "nonlinear elimination present in 47.5% and time-dependent clearance in 19.4%",
      "(van den Berg 2025 Table 1)."
    ),
    notes = paste(
      "n_subjects is the SUM of the per-model patient counts across the 157 of 160 models that",
      "reported one (Supplementary Table S1 column N); patients are double-counted wherever",
      "several models were fitted to overlapping datasets, so it is an upper bound on distinct",
      "individuals rather than a cohort size. Median (IQR) samples per model 5250",
      "(1788; 10223), with a median (IQR) of 8.6 (5.7; 13.3) samples per patient",
      "(van den Berg 2025 Results 'Sample size matters'). The absorption parameters ka and F",
      "are medians over the 52 models built on both i.v. and s.c. data; the disposition",
      "parameters are medians over all 160."
    )
  )

  ini({
    # =====================================================================
    # Generic mAb model. All values are the across-model medians reported in
    # the van den Berg 2025 Figure 5 caption and re-tabulated at the head of
    # Supplementary Section M5, and every one of them is FIX in the paper's
    # NONMEM 7.5.0 control stream ($THETA ... FIX / $OMEGA ... FIX), so every
    # one is fixed() here. Nothing in this model was estimated from data by
    # these authors -- the "estimation" was the meta-analysis itself.
    #
    # The medians reproduce exactly from the deposited 160-model dataset
    # (Supplementary Table S1 / KMAB_A_2512217_SM8524.xlsx sheet "Models"):
    # CL 0.2200, Vc 3.4230, Vp 2.6800, Q 0.5390, ka 0.2485, F 0.6945.
    # =====================================================================

    # ----- Structural disposition parameters (medians over all 160 models) -----
    lcl <- fixed(log(0.22));  label("Clearance CL (L/day)")                          # van den Berg 2025 Figure 5c caption "CL = 0.22 L/d"; Suppl. M5 $THETA 0.22 FIX; Results: median (IQR) 0.22 (0.17; 0.29) L/d
    lvc <- fixed(log(3.42));  label("Central volume of distribution Vc (L)")          # van den Berg 2025 Figure 5c caption "VC = 3.42 L"; Suppl. M5 $THETA 3.42 FIX; Results: median (IQR) 3.42 (2.96; 3.99) L
    lvp <- fixed(log(2.68));  label("Peripheral volume of distribution Vp (L)")       # van den Berg 2025 Figure 5c caption "VP = 2.68 L"; Suppl. M5 $THETA 2.68 FIX; Results: median (IQR) 2.68 (2.09; 3.54) L
    lq  <- fixed(log(0.54));  label("Intercompartmental clearance Q (L/day)")         # van den Berg 2025 Figure 5c caption "Q = 0.54 L/d"; Suppl. M5 $THETA 0.54 FIX; Results: median (IQR) 0.54 (0.36; 0.84) L/d

    # ----- Subcutaneous absorption (medians over the 52 i.v./s.c. models) -----
    # The abstract prints ka's unit as "L/d", which is dimensionally impossible
    # for a first-order rate constant; the Figure 5c caption and the Suppl. M5
    # $THETA comment both give "/d", which is the value used here.
    lka     <- fixed(log(0.25));  label("First-order absorption rate constant ka (1/day)")  # van den Berg 2025 Figure 5c caption "ka = 0.25/d"; Suppl. M5 $THETA 0.25 FIX ; T_ka (/d); Results: median (IQR) 0.25 (0.18; 0.29) /d
    lfdepot <- fixed(log(0.69));  label("Subcutaneous bioavailability F (fraction)")        # van den Berg 2025 Figure 5c caption "F = 0.69"; Suppl. M5 $THETA 0.69 FIX and F1 = THETA(6); Results: median (IQR) 69 (61; 77) %

    # =====================================================================
    # Inter-individual variability. These are BETWEEN-SUBJECT variances --
    # the medians of the between-subject omegas reported by the source
    # models -- not the study-level random effect of a regression MBMA, so
    # they carry the standard etal<param> names. ini() takes variances; the
    # paper gives both the variance and the log-domain CV%, which are
    # mutually consistent via CV = sqrt(exp(omega^2) - 1):
    #   sqrt(exp(0.097) - 1) = 0.319, sqrt(exp(0.056) - 1) = 0.240,
    #   sqrt(exp(0.26)  - 1) = 0.545 (paper prints 54.8%, i.e. omega^2 = 0.26
    #   is the rounded value).
    # There is deliberately no IIV on Vp or Q: the paper's generic model
    # assigns random effects only to CL, Vc and ka (Suppl. M5 $PK / $OMEGA).
    # =====================================================================
    etalcl ~ fixed(0.097)   # van den Berg 2025 Figure 5c caption 'IIVCL = 31.9% (omega^2 = 0.097)'; Suppl. M5 $OMEGA 0.097 FIX ; IIV_CL; median (IQR) final-model IIV on CL 31.9 (28.2-38) CV% (Figure 3a)
    etalvc ~ fixed(0.056)   # van den Berg 2025 Figure 5c caption 'IIVVC = 24% (omega^2 = 0.056)'; Suppl. M5 $OMEGA 0.056 FIX ; IIV_V1; median (IQR) final-model IIV on Vc 24 (18-34.5) CV% (Figure 3a)
    etalka ~ fixed(0.26)    # van den Berg 2025 Figure 5c caption 'IIVka = 54.8% (omega^2 = 0.26)'; Suppl. M5 $OMEGA 0.26 FIX ; IIV_ka

    # =====================================================================
    # Residual error. The paper states "No RUV or covariates were
    # introduced"; the control stream implements that by fixing the
    # proportional error SD to 1e-5 with SIGMA = 1 FIX, which is what is
    # carried here so the encoding matches the authors' own simulation
    # exactly. It is a numerical placeholder, not a measurement-error
    # estimate -- refit it before using this model to fit data.
    # =====================================================================
    propSd <- fixed(0.00001); label("Proportional residual error SD (fraction; numerical placeholder, no RUV was modelled)")  # van den Berg 2025 Suppl. M5 $THETA 0.00001 FIX ; T_Prop_RE (sd), with $SIGMA 1 FIX; Methods 'Generic model simulation': "No RUV or covariates were introduced"
  })

  model({
    # Individual parameters: typical value times exp(eta), matching
    # Suppl. M5 $PK (CL = THETA(1)*EXP(ETA(1)) etc.). Vp and Q carry no
    # random effect.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    ka <- exp(lka + etalka)
    vp <- exp(lvp)
    q  <- exp(lq)

    # Micro-constants, as K20 / K23 / K32 in Suppl. M5 $PK.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system, as Suppl. M5 $DES (DEPOT / CENTRAL / PERI).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Subcutaneous bioavailability; intravenous doses go straight to
    # central and are unaffected (Suppl. M5 F1 = THETA(6)).
    f(depot) <- exp(lfdepot)

    # Central amount in mg over Vc in L gives mg/L, matching S2 = VC.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
