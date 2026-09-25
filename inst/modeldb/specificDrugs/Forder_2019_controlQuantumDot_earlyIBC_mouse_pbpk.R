Forder_2019_controlQuantumDot_earlyIBC_mouse_pbpk <- function() {
  description <- paste(
    "Preclinical (mouse, female, ~22 g). PBPK (flow-limited, perfusion-limited",
    "well-stirred organs; MATLAB ode23s) model for the biodistribution of a",
    "single intravenous dose of unconjugated control quantum dots (ConQD) in mice with early-stage inflammatory breast cancer (IBC).",
    "Eight flow-limited, well-stirred compartments: the healthy six plus a `breast` compartment in parallel with the other organs and a `tumor` housed inside the breast (Figure 1 middle). The tumor exchanges only with breast tissue, so breast flow to the tumor recirculates back to the breast (Eq. 17-20, as implemented in the supplementary MATLAB code).",
    "All volumes and flows are fixed physiology for a 22 g mouse (Table 1);",
    "the only fitted parameters are the kidney, liver, spleen, lung, tumor and breast",
    "partition coefficients (Table 3), estimated by relative residual sum of",
    "squares against day-4 fluorescence intensities (Table 2), with the",
    "`other` partition coefficient held at 1 and a negligible liver",
    "elimination term (kelim = 1e-6). States hold CONCENTRATION in the",
    "code-reported nM intensity scale; the dose (pmol) is divided by the",
    "plasma volume. No between-subject variability or residual error is",
    "reported, so the model is for typical-value simulation only.",
    sep = " "
  )
  reference <- paste(
    "Forder J, Smith M, Wagner M, Schaefer RJ, Gorky J, van Golen KL, Nohe A,",
    "Dhurjati P (2019). A physiologically-based pharmacokinetic model for",
    "targeting calcitriol-conjugated quantum dots to inflammatory breast",
    "cancer cells. Clinical and Translational Science 12(6):617-624.",
    "doi:10.1111/cts.12664.",
    "Structural equations and volumes/flows are taken from the supplementary",
    "MATLAB model code (CTS-12-617-s006) where it differs from the printed",
    "equations; see the vignette.",
    sep = " "
  )
  vignette <- "Forder_2019_calcitriolQuantumDots"
  units <- list(time = "h", dosing = "pmol", concentration = "nM")
  dosing <- "plasma"

  # Breast (mammary) tissue housing the early-stage tumor; no registered
  # canonical exists for a breast-tissue PBPK compartment.
  paper_specific_compartments <- c("breast")

  compartmentData <- list(
    plasma = list(
      analyte = "unconjugated (control) carboxyl quantum dot",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    kidney = list(
      analyte = "unconjugated (control) carboxyl quantum dot",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    liver = list(
      analyte = "unconjugated (control) carboxyl quantum dot",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    spleen = list(
      analyte = "unconjugated (control) carboxyl quantum dot",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    lung = list(
      analyte = "unconjugated (control) carboxyl quantum dot",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    other = list(
      analyte = "unconjugated (control) carboxyl quantum dot",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    tumor = list(
      analyte = "unconjugated (control) carboxyl quantum dot",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    breast = list(
      analyte = "unconjugated (control) carboxyl quantum dot",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (female, 13-16 weeks, ~22 g)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "13-16 weeks",
    weight_range = "~22.0 g (fixed model body weight)",
    sex_female_pct = 100,
    disease_state = "early-stage inflammatory breast cancer (localized SUM149 tumor within breast tissue)",
    dose_range = paste(
      "Single intravenous injection of ConQD; the model code starts plasma at",
      "40 nM (43.12 pmol in 1.078 mL)."
    ),
    notes = paste(
      "Partition coefficients were fitted (MATLAB fmincon, relative squared",
      "error) to mean day-4 (96 h) fluorescence pixel intensities from",
      "Schaefer 2012 (Forder 2019 Table 2, late-stage tumor column): kidney 12.54, liver 14.85, lung 8.22, spleen 26.43.",
      "No early-stage fluorescence data exist; the early-stage model was",
      "fitted to the late-stage data. Animal numbers are not reported."
    )
  )

  ini({
    # ================================================================
    # Physiology -- Forder 2019 Table 1 and supplementary MATLAB code
    # (CTS-12-617-s006). Literature constants, not estimated.
    # ================================================================
    v_plasma <- fixed(1.078); label("Plasma volume (mL)") # Table 1; code vp = 1.078
    v_kidney <- fixed(0.3674); label("Kidney volume (mL)") # Table 1; code vk = 0.3674
    v_liver <- fixed(1.2078); label("Liver volume (mL)") # Table 1 1.208; code vli = 1.2078
    v_spleen <- fixed(0.077); label("Spleen volume (mL)") # Table 1 0.077 (enlarged 0.154); code vs = 0.077
    v_lung <- fixed(0.1606); label("Lung volume (mL)") # Table 1; code vl = 0.1606
    v_other <- fixed(55.7572); label("Lumped other-tissue volume (mL)") # Table 1 55.76; code vo = 55.7572
    v_tumor <- fixed(0.08); label("Tumor volume (mL)") # Table 1; code vt = 0.08 (Schaefer 2012 tumor measurement)
    v_breast <- fixed(1.05); label("Breast (mammary) tissue volume (mL)") # Table 1; code vb = 1.05

    q_kidney <- fixed(85.7714); label("Kidney plasma flow (mL/h)") # Table 1 85.77; code qk = 85.7714
    q_liver <- fixed(132.8985); label("Hepatic arterial plasma flow (mL/h)") # Table 1 132.9; code qli = 132.8985
    q_spleen <- fixed(18.8509); label("Spleen plasma flow (mL/h)") # Table 1 18.85; code qs = 18.8509
    q_lung <- fixed(4.7127); label("Lung plasma flow (mL/h)") # Table 1 4.713; code ql = 4.7127
    q_other <- fixed(677.5191); label("Lumped other-tissue plasma flow (mL/h)") # Table 1 'Other' 677.5; code Tsm3QD() qo = 677.5191 (700.3091 minus the breast flow)
    q_tumor <- fixed(1.7664); label("Tumor plasma flow (mL/h)") # Table 1 1.766; code qt = 1.7664
    q_breast <- fixed(22.79); label("Breast plasma flow (mL/h)") # Table 1 22.79; code qb = 22.79

    # Eq. 3: kelim = 1.0e-6 per hour, from QD 705 (Lin 2008); code k1 = 0.000001.
    # In Eq. 12 it multiplies the liver concentration inside the 1/V_liver
    # bracket, so as written it acts as a clearance.
    kelim <- fixed(1e-6); label("Liver (biliary) elimination coefficient (1/h)") # Eq. 3

    # ================================================================
    # Partition coefficients -- Table 3, Early-stage tumor, ConQD column
    # (relative residual sum of squares fit; no uncertainty reported)
    # ================================================================
    lkp_kidney <- log(45.1); label("Log kidney:plasma partition coefficient (log unitless)") # Table 3, Early-stage tumor, ConQD
    lkp_liver <- log(53.4); label("Log liver:plasma partition coefficient (log unitless)") # Table 3, Early-stage tumor, ConQD
    lkp_spleen <- log(95.0); label("Log spleen:plasma partition coefficient (log unitless)") # Table 3, Early-stage tumor, ConQD
    lkp_lung <- log(29.6); label("Log lung:plasma partition coefficient (log unitless)") # Table 3, Early-stage tumor, ConQD
    lkp_tumor <- log(28.6); label("Log tumor:plasma partition coefficient (log unitless)") # Table 3, Early-stage tumor, ConQD
    lkp_breast <- log(2.67); label("Log breast:plasma partition coefficient (log unitless)") # Table 3, Early-stage tumor, ConQD
    lkp_other <- fixed(log(1)); label("Log other-tissue:plasma partition coefficient (log unitless)") # Results, Partition coefficients: assumed to be 1; code ro = 1
  })
  model({
    # States hold concentration (nM); each organ is flow limited:
    # dC/dt = Q / V * (C_in - C / P) (Eq. 6).
    kp_kidney <- exp(lkp_kidney)
    kp_liver <- exp(lkp_liver)
    kp_spleen <- exp(lkp_spleen)
    kp_lung <- exp(lkp_lung)
    kp_other <- exp(lkp_other)
    kp_tumor <- exp(lkp_tumor)
    kp_breast <- exp(lkp_breast)

    # Eq. 17 / supplementary code Tsm3QD(): the breast returns its own flow to
    # plasma; the tumor flow loops breast -> tumor -> breast.
    q_total <- q_kidney + q_liver + q_spleen + q_lung + q_other + q_breast # Eq. 18

    d/dt(plasma) <- (q_kidney * kidney / kp_kidney +
      (q_liver + q_spleen) * liver / kp_liver +
      q_lung * lung / kp_lung +
      q_other * other / kp_other +
      q_breast * breast / kp_breast -
      q_total * plasma) / v_plasma
    d/dt(kidney) <- q_kidney / v_kidney * (plasma - kidney / kp_kidney) # Eq. 9
    # Eq. 12: spleen outflow enters the liver; kelim sits inside the 1/V bracket.
    d/dt(liver) <- (q_liver * plasma + q_spleen * spleen / kp_spleen -
      (q_liver + q_spleen) * liver / kp_liver - kelim * liver) / v_liver
    d/dt(spleen) <- q_spleen / v_spleen * (plasma - spleen / kp_spleen) # Eq. 11
    d/dt(lung) <- q_lung / v_lung * (plasma - lung / kp_lung) # Eq. 10 (C_b = plasma)
    d/dt(other) <- q_other / v_other * (plasma - other / kp_other) # Eq. 13
    # Supplementary code Tsm3QD(): the tumor is fed by the breast VENOUS
    # concentration breast / kp_breast, and the breast loses q_tumor at that
    # same concentration. The printed Eq. 19-20 omit the "/ kp_breast" on both
    # terms; the code form is kept because it reproduces the Table 2 fit of the
    # Table 3 early-stage ConQD and CalQD partition coefficients exactly.
    d/dt(breast) <- (q_breast * plasma + q_tumor * tumor / kp_tumor -
      q_tumor * breast / kp_breast - q_breast * breast / kp_breast) / v_breast
    d/dt(tumor) <- q_tumor / v_tumor * (breast / kp_breast - tumor / kp_tumor)

    # Intravenous pulse into plasma: the dose (pmol) divided by the plasma
    # volume (mL) gives the initial plasma concentration (nM).
    f(plasma) <- 1 / v_plasma

    Cc <- plasma
  })
}
