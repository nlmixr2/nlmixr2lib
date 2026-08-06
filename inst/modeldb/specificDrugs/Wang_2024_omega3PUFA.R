Wang_2024_omega3PUFA <- function() {
  description <- "One-compartment population PK model with first-order oral absorption for omega-3 polyunsaturated fatty acids (PUFA) in community-dwelling adults with type 2 diabetes, with an HDL-cholesterol power covariate on Ka, V and CL, coupled to a direct-effect inhibitory Imax model for glycosylated hemoglobin (HbA1c) (Wang 2024)"
  reference <- "Wang L, Huang X, Sun M, Zheng T, Zheng L, Lin X, Ruan J, Lin F. New light on omega-3 polyunsaturated fatty acids and diabetes debate: a population pharmacokinetic-pharmacodynamic modelling and intake threshold study. Nutr Diabetes. 2024;14(1):8. doi:10.1038/s41387-024-00262-w"
  vignette <- "Wang_2024_omega3PUFA"
  units <- list(time = "h", dosing = "g", concentration = "g/L")

  # Issue #482: what molecule each compartment holds, in what units, in what
  # biological matrix. Verified against Wang 2024 Methods (oral intake, plasma
  # omega-3 PUFA measured by ELISA) and the printed final-model equations.
  compartmentData <- list(
    depot   = list(analyte = "omega-3 PUFA", units = "g", specimen = "administration site", verified = TRUE),
    central = list(analyte = "omega-3 PUFA", units = "g", specimen = "plasma",              verified = TRUE)
  )

  covariateData <- list(
    HDLC = list(
      description        = "Serum high-density lipoprotein cholesterol",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. The only covariate retained in the final model.",
        "Enters Ka, V and CL as a power function normalised to the cohort mean of",
        "1.38 mmol/L (Wang 2024 Table 1: HDL 1.38 +/- 0.38 mmol/L, median 1.28,",
        "range 0.77-2.49). The Ka exponent (0.007) is numerically negligible but is",
        "retained because the paper prints it as part of the final model."
      ),
      source_name        = "HDL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 161L,
    n_studies      = 1L,
    study          = "Two-centre prospective community cohort, Fuzhou, China, 2020-2021 (ChiCTR2000036210)",
    age_range      = "43-85 years (mean 67.16, median 67)",
    height_range   = "136-184 cm (mean 159.91, median 159)",
    weight_range   = "39-102 kg (mean 63.52, median 63)",
    bmi_range      = "16.02-43.79 kg/m^2 (mean 24.82, median 24.61)",
    sex_female_pct = 60.2,
    race_ethnicity = "Not reported; recruited from two communities in Fuzhou, Fujian Province, China",
    disease_state  = "Type 2 diabetes, community-dwelling, aged 18 years or older on a regular diet",
    dose_range     = "Habitual dietary omega-3 PUFA intake 0.036-6.426 g/day (mean 1.14, median 0.774), quantified by food-frequency questionnaire",
    regions        = "China (Fuzhou, Fujian Province)",
    notes          = paste(
      "Baseline demographics from Wang 2024 Table 1. Baseline HbA1c 6.31 +/- 1.91 %",
      "(median 5.7, range 2-14.1); plasma omega-3 PUFA 0.0319 +/- 0.0091 g/L",
      "(median 0.0299, range 0.0107-0.0508). Sex 64 male / 97 female. A further 49",
      "patients with type 2 diabetes formed the external validation set (PK: MDPE",
      "11.23%, MAPE 29.00%, F20 38.78%, F30 81.63%; PK-PD: MDPE -6.87%, MAPE 15.32%,",
      "F20 69.39%, F30 85.71%). Estimation used FOCE-ELS in Phoenix NLME 1.2."
    )
  )

  ini({
    # ---- Structural PK: Wang 2024 Table 2, "Final model (CV%)" column ----------
    # The parenthesised value in that column is the relative standard error of the
    # estimate, not an IIV: the residual-error term sigma also carries one, and a
    # residual-error parameter has no inter-individual variability.
    lka <- log(1.175);   label("Absorption rate constant at the reference HDL-C of 1.38 mmol/L (Ka, 1/h)")               # Table 2: Ka = 1.175 (RSE 0.096%); printed final-model equation p. 3
    lvc <- log(26.151);  label("Apparent central volume of distribution at the reference HDL-C of 1.38 mmol/L (V/F, L)") # Table 2: V = 26.151 L (RSE 1.451%); printed final-model equation p. 3
    lcl <- log(0.411);   label("Apparent clearance at the reference HDL-C of 1.38 mmol/L (CL/F, L/h)")                   # Table 2: CL = 0.411 L/h (RSE 0.000%); printed final-model equation p. 3

    # ---- HDL-C power covariate exponents: printed final-model equations, p. 3 --
    # Ka = 1.175 * (HDL/1.38)^0.007  * exp(etaKa)
    # V  = 26.151 * (HDL/1.38)^-0.535 * exp(etaV)
    # CL = 0.411 * (HDL/1.38)^0.285  * exp(etaCL)
    e_hdlc_ka <-  0.007; label("HDL-C power exponent on Ka (unitless)")   # Wang 2024 p. 3, printed Ka equation
    e_hdlc_vc <- -0.535; label("HDL-C power exponent on V/F (unitless)")  # Wang 2024 p. 3, printed V equation (exponent is negative)
    e_hdlc_cl <-  0.285; label("HDL-C power exponent on CL/F (unitless)") # Wang 2024 p. 3, printed CL equation

    # ---- Direct-effect inhibitory Imax PD on HbA1c ----------------------------
    # Wang 2024 Table 2, "Final model (CV%)" column. Imax is an ABSOLUTE maximal
    # HbA1c reduction in percentage points, not a fraction of E0: only the
    # additive form HbA1c = E0 - Imax * Cc / (IC50 + Cc) reproduces the paper's
    # own simulation outputs (see the model() block and the vignette).
    le0    <- log(5.641); label("Baseline (drug-free) HbA1c, E0 (%)")                                     # Table 2: E0 = 5.641 % (RSE 2.42%)
    lec50  <- log(0.090); label("Plasma omega-3 PUFA giving half-maximal HbA1c reduction, IC50 (g/L)")    # Table 2: IC50 = 0.090 g/L (RSE 4.776%)
    limax  <- log(0.597); label("Maximal HbA1c reduction, Imax (HbA1c percentage points)")                # Table 2: Imax = 0.597 (RSE 8.262%)

    # ---- Inter-individual variability ----------------------------------------
    # Wang 2024 writes exp(etaKa), exp(etaV), exp(etaCL), exp(etaE0),
    # exp(etaIC50) and exp(etaImax) into the final-model equations, but never
    # reports a single omega. The eta structure is preserved so the published
    # model form is faithful; every variance is fixed at zero because no value
    # exists on disk to populate it. See vignette "Assumptions and deviations".
    etalka   ~ fixed(0)  # Wang 2024 p. 3: exp(etaKa) printed, magnitude never reported
    etalvc   ~ fixed(0)  # Wang 2024 p. 3: exp(etaV) printed, magnitude never reported
    etalcl   ~ fixed(0)  # Wang 2024 p. 3: exp(etaCL) printed, magnitude never reported
    etale0   ~ fixed(0)  # Wang 2024 p. 3: exp(etaE0) printed, magnitude never reported
    etalec50 ~ fixed(0)  # Wang 2024 p. 3: exp(etaIC50) printed, magnitude never reported
    etalimax ~ fixed(0)  # Wang 2024 p. 3: exp(etaImax) printed, magnitude never reported

    # ---- Residual error -------------------------------------------------------
    # Table 2 reports a single unlabelled "sigma" per sub-model with no error-model
    # statement. An additive interpretation is dimensionally impossible for the PK
    # observation (sigma = 0.354 against plasma levels of ~0.03 g/L), so both
    # sigmas are encoded as proportional.
    propSd       <- 0.354; label("Proportional residual error on plasma omega-3 PUFA (fraction)")  # Table 2 PK block: sigma = 0.354
    propSd_hba1c <- 0.354; label("Proportional residual error on HbA1c (fraction)")                # Table 2 PD block: sigma = 0.354
  })

  model({
    # HDL-C power covariate, normalised to the cohort mean HDL-C of 1.38 mmol/L
    # (Wang 2024 Table 1). The journal typesets the whole term
    # "<exponent> * exp(eta)" inside the superscript; that literal reading is
    # degenerate (it would abolish all IIV for any subject at HDL-C = 1.38, since
    # 1^x = 1), so the standard Phoenix NLME form theta * (cov/ref)^theta_cov *
    # exp(eta) is used. See the vignette "Assumptions and deviations".
    ka <- exp(lka + etalka) * (HDLC / 1.38)^e_hdlc_ka
    vc <- exp(lvc + etalvc) * (HDLC / 1.38)^e_hdlc_vc
    cl <- exp(lcl + etalcl) * (HDLC / 1.38)^e_hdlc_cl

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in g, volume in L -> g/L, matching the ELISA readout and the IC50 units
    Cc <- central / vc

    # Direct-effect inhibitory Imax on HbA1c. Wang 2024 prints the three
    # individual-parameter equations for E0, IC50 and Imax but never the effect
    # equation itself. The additive form below is the reading that reproduces the
    # paper's own reported simulation outputs: with a single lognormal spread on
    # E0 back-solved from the zero-dose intercepts of Figure 2, it recovers all
    # four Figure 2 target curves to within ~1-3 percentage points and both
    # quoted claims (0.4 g/day -> >95% of patients below 7.0%; 1.3 g/day ->
    # 90.08% below 6.5%). The alternative fractional reading
    # E0 * (1 - Imax * Cc / (IC50 + Cc)) drives attainment to ~100% at every dose
    # at or above 0.4 g/day and is inconsistent with Figure 2.
    e0    <- exp(le0    + etale0)
    ec50  <- exp(lec50  + etalec50)
    imax  <- exp(limax  + etalimax)

    hba1c <- e0 - imax * Cc / (ec50 + Cc)

    Cc    ~ prop(propSd)
    hba1c ~ prop(propSd_hba1c)
  })
}
