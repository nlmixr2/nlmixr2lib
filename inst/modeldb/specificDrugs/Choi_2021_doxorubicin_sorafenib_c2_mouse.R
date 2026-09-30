Choi_2021_doxorubicin_sorafenib_c2_mouse <- function() {
  description <- paste(
    "Preclinical (mouse; PK in male athymic nude mice, tumor growth in female",
    "SCID mice bearing orthotopic human 143B osteosarcoma xenografts).",
    "Doxorubicin (IV) plus sorafenib (oral) combination-therapy PK-PD model",
    "(Choi 2021 conventional model C2, the Koch 2009 interaction form with",
    "the interaction factor assigned to sorafenib): each drug's two-compartment",
    "PK model from the monotherapy fits drives a single Simeoni-type",
    "tumor-growth-inhibition model with Koch 2009 exponential-then-linear",
    "natural growth, in which the kill rate is the sum of each drug's",
    "monotherapy potency times its concentration with the sorafenib term",
    "multiplied by an estimated interaction factor, followed by a three-stage",
    "damaged-cell transit chain. Reported by Choi 2021 as the comparator",
    "for the new model D. Typical-value (naive-pooled SimBiology) fit with",
    "an exponential residual error and no between-animal variability."
  )
  reference <- paste(
    "Choi YH, Zhang C, Liu Z, Tu MJ, Yu AX, Yu AM.",
    "A Novel Integrated Pharmacokinetic-Pharmacodynamic Model to Evaluate",
    "Combination Therapy and Determine In Vivo Synergism.",
    "J Pharmacol Exp Ther. 2021;377(3):305-315. doi:10.1124/jpet.121.000584.",
    "Tumor-growth data and dosing regimen from Jian C, Tu MJ, Ho PY, et al.",
    "Oncotarget. 2017;8(19):30742-30755. doi:10.18632/oncotarget.16372.",
    sep = " "
  )
  vignette <- "Choi_2021_doxorubicin_sorafenib_mouse"

  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central_dox = list(analyte = "doxorubicin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_dox = list(analyte = "doxorubicin", units = "mg", specimen = "tissue", verified = TRUE),
    depot_sorafenib = list(analyte = "sorafenib", units = "mg", specimen = "administration site", verified = TRUE),
    central_sorafenib = list(analyte = "sorafenib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_sorafenib = list(analyte = "sorafenib", units = "mg", specimen = "tissue", verified = TRUE),
    cycling_cells = list(analyte = "proliferating tumor cells", units = "cm^3", specimen = "tumor", verified = TRUE),
    damaged_cells1 = list(analyte = "drug-damaged tumor cells", units = "cm^3", specimen = "tumor", verified = TRUE),
    damaged_cells2 = list(analyte = "drug-damaged tumor cells", units = "cm^3", specimen = "tumor", verified = TRUE),
    damaged_cells3 = list(analyte = "drug-damaged tumor cells", units = "cm^3", specimen = "tumor", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "mouse (PK: male athymic nude Foxn1nu; tumor growth: female CB17 SCID with orthotopic 143B osteosarcoma xenograft)",
    n_subjects = 23L,
    n_studies = 2L,
    age_range = "7 weeks at purchase (PK mice)",
    weight_range = "approximately 30 g (PK mice)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state = "PK: tumor-free mice. PD: orthotopic (intratibial) human 143B osteosarcoma cell-line-derived xenograft.",
    dose_range = "PD: doxorubicin 12 ug/mouse IV every other day x4 then 10 ug/mouse IV every other day x4, plus sorafenib 200 ug/mouse oral gavage on 7 days of cycle 1 then every other day x4 in cycle 2, per Jian 2017.",
    regions = "University of California Davis, USA (preclinical)",
    notes = paste(
      "PK: n = 6 male athymic nude mice per arm (doxorubicin IV, sorafenib",
      "IV, sorafenib oral; Choi 2021 Table 1). PD: 5 of the 7 combination-",
      "treated tumor-bearing mice from Jian 2017 were randomly chosen for model",
      "development and the remaining 2 used for verification (Choi 2021",
      "Methods; Table 3). Only w0, k1' and the interaction factor C_Sor",
      "were estimated; the PK, the natural-growth rates and each drug's",
      "potency k2 were held at their earlier estimates."
    )
  )

  ini({
    # --- Doxorubicin PK: Choi 2021 Table 1, 'Dox (iv)' column (held fixed) --
    # Rate constants are printed in 1/h; model time is days, so model()
    # multiplies each rate constant by 24.
    lvc_dox  <- fixed(log(0.0764)); label("Doxorubicin central volume Vc (L)")                             # Table 1 Dox (iv): Vc = 0.0764 L
    lvp_dox  <- fixed(log(3.52));   label("Doxorubicin peripheral volume Vp (L)")                          # Table 1 Dox (iv): Vp = 3.52 L
    lk12_dox <- fixed(log(5.86));   label("Doxorubicin central-to-peripheral rate constant k12 (1/h)")     # Table 1 Dox (iv): k12 = 5.86 1/h
    lk21_dox <- fixed(log(0.127));  label("Doxorubicin peripheral-to-central rate constant k21 (1/h)")     # Table 1 Dox (iv): k21 = 0.127 1/h
    lkel_dox <- fixed(log(1.24));   label("Doxorubicin elimination rate constant ke (1/h)")                # Table 1 Dox (iv): ke = 1.24 1/h

    # --- Sorafenib PK: Choi 2021 Table 1 (held fixed) -----------------------
    lka_sorafenib  <- fixed(log(0.640));  label("Sorafenib oral absorption rate constant ka (1/h)")               # Table 1 Sor (p.o.): ka = 0.640 1/h
    lvc_sorafenib  <- fixed(log(0.0112)); label("Sorafenib central volume Vc (L)")                                # Table 1 Sor (iv): Vc = 0.0112 L
    lvp_sorafenib  <- fixed(log(0.0120)); label("Sorafenib peripheral volume Vp (L)")                             # Table 1 Sor (iv): Vp = 0.0120 L
    lk12_sorafenib <- fixed(log(0.663));  label("Sorafenib central-to-peripheral rate constant k12 (1/h)")        # Table 1 Sor (iv): k12 = 0.663 1/h
    lk21_sorafenib <- fixed(log(0.620));  label("Sorafenib peripheral-to-central rate constant k21 (1/h)")        # Table 1 Sor (iv): k21 = 0.620 1/h
    lkel_sorafenib <- fixed(log(0.344));  label("Sorafenib elimination rate constant ke (1/h)")                   # Table 1 Sor (iv): ke = 0.344 1/h

    # --- Natural growth and drug potencies (held fixed) ---------------------
    # Table 3 footnote a: L0 and L1 fixed at the Table 2 control estimates.
    # Methods (model D): 'potency of each drug (k2A and k2B) remains
    # constant', i.e. the monotherapy k2 values from Table 2.
    ltumorExpGrowth     <- fixed(log(0.107));  label("Exponential-phase tumor growth rate L0 (1/day)")          # Table 3 footnote a / Table 2 Control: L0 = 0.107 1/day
    ltumorLinGrowth     <- fixed(log(0.148));  label("Linear-phase tumor growth rate L1 (cm^3/day)")            # Table 3 footnote a / Table 2 Control: L1 = 0.148 cm^3/day
    ldrugSlope_dox      <- fixed(log(13.8));   label("Doxorubicin potency k2A (L/mg/day)")                      # Table 2 Dox: k2 = 13.8 L/mg/day
    ldrugSlope_sorafenib <- fixed(log(0.0267)); label("Sorafenib potency k2B (L/mg/day)")                       # Table 2 Sor: k2 = 0.0267 L/mg/day

    # --- Combination-therapy parameters: Choi 2021 Table 3, 'Model C2' ------
    # Figure 1C(ii): the interaction factor psi_B multiplies the sorafenib kill term only.
    lrbase_tumor   <- log(0.0405); label("Initial tumor volume w0 (cm^3)")                              # Table 3 Model C2: w0 = 0.0405 cm^3 (12.3 %CV)
    ldamageTransit <- log(0.451); label("Combination damaged-cell transit rate constant k1' (1/day)")  # Table 3 Model C2: k'1 = 0.451 1/day (69.6 %CV)
    linteract_sorafenib <- log(1.56); label("Interaction factor C_Sor on the sorafenib kill term (unitless)")    # Table 3 Model C2: C_Sor = 1.56 (21.9 %CV)

    # --- Residual error ------------------------------------------------------
    # Exponential error model (Choi 2021 Methods); magnitudes derived as
    # sqrt(MSE). The two concentration SDs come from the monotherapy PK fits
    # and are held fixed; the tumor-volume SD is from the model C2 fit
    # (Table 3 MSE = SSE / (n - p): 2.64 / 0.120 = 22 = 25 - 3).
    expSd_dox          <- fixed(sqrt(0.290)); label("Exponential residual error on doxorubicin concentration (SD, log scale)")  # derived: sqrt(Table 1 Dox (iv) MSE = 0.290)
    expSd_sorafenib <- fixed(sqrt(0.337)); label("Exponential residual error on sorafenib concentration (SD, log scale)")    # derived: sqrt(Table 1 Sor (p.o.) MSE = 0.337)
    expSd_tumor_vol    <- sqrt(0.120);        label("Exponential residual error on tumor volume (SD, log scale)")              # derived: sqrt(Table 3 Model C2 MSE = 0.120)
  })

  model({
    # Choi 2021 prints the PK rate constants in 1/h; x24 converts to 1/day.
    vc_dox  <- exp(lvc_dox)
    vp_dox  <- exp(lvp_dox)
    k12_dox <- exp(lk12_dox) * 24
    k21_dox <- exp(lk21_dox) * 24
    kel_dox <- exp(lkel_dox) * 24

    ka_sorafenib  <- exp(lka_sorafenib) * 24
    vc_sorafenib  <- exp(lvc_sorafenib)
    vp_sorafenib  <- exp(lvp_sorafenib)
    k12_sorafenib <- exp(lk12_sorafenib) * 24
    k21_sorafenib <- exp(lk21_sorafenib) * 24
    kel_sorafenib <- exp(lkel_sorafenib) * 24

    tumorExpGrowth     <- exp(ltumorExpGrowth)
    tumorLinGrowth     <- exp(ltumorLinGrowth)
    drugSlope_dox      <- exp(ldrugSlope_dox)
    drugSlope_sorafenib <- exp(ldrugSlope_sorafenib)
    rbase_tumor        <- exp(lrbase_tumor)
    damageTransit      <- exp(ldamageTransit)
    interact_sorafenib <- exp(linteract_sorafenib)

    # --- Doxorubicin PK, IV (Choi 2021 Equations 1-2; amounts in mg) --------
    d/dt(central_dox)     <- -k12_dox * central_dox - kel_dox * central_dox + k21_dox * peripheral1_dox
    d/dt(peripheral1_dox) <-  k12_dox * central_dox - k21_dox * peripheral1_dox
    Cc_dox <- central_dox / vc_dox

    # --- Sorafenib PK, oral (Choi 2021 Equations 3-5; amounts in mg) --------
    d/dt(depot_sorafenib)       <- -ka_sorafenib * depot_sorafenib
    d/dt(central_sorafenib)     <-  ka_sorafenib * depot_sorafenib - k12_sorafenib * central_sorafenib -
      kel_sorafenib * central_sorafenib + k21_sorafenib * peripheral1_sorafenib
    d/dt(peripheral1_sorafenib) <-  k12_sorafenib * central_sorafenib - k21_sorafenib * peripheral1_sorafenib
    Cc_sorafenib <- central_sorafenib / vc_sorafenib

    # --- Combination tumor growth inhibition (Choi 2021 Figure 1C(ii)) -------
    # Kill rate = k2A * C_A + psi_B * k2B * C_B.
    tumor_vol <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3
    drugEffectCyclingCells <- drugSlope_dox * Cc_dox +
      interact_sorafenib * drugSlope_sorafenib * Cc_sorafenib

    d/dt(cycling_cells)  <- 2 * tumorExpGrowth * tumorLinGrowth * cycling_cells^2 /
      ((tumorLinGrowth + 2 * tumorExpGrowth * cycling_cells) * tumor_vol) -
      drugEffectCyclingCells * cycling_cells
    d/dt(damaged_cells1) <- drugEffectCyclingCells * cycling_cells - damageTransit * damaged_cells1
    d/dt(damaged_cells2) <- damageTransit * (damaged_cells1 - damaged_cells2)
    d/dt(damaged_cells3) <- damageTransit * (damaged_cells2 - damaged_cells3)

    cycling_cells(0)  <- rbase_tumor
    damaged_cells1(0) <- 0
    damaged_cells2(0) <- 0
    damaged_cells3(0) <- 0

    Cc_dox ~ lnorm(expSd_dox)
    Cc_sorafenib ~ lnorm(expSd_sorafenib)
    tumor_vol ~ lnorm(expSd_tumor_vol)
  })
}
