Maximova_2022_acyclovir <- function() {
  description <- "One-compartment population PK model for intravenous acyclovir in 120 oncologic children (age 0-18 years; mean weight 32.4 kg) receiving acyclovir prophylaxis or treatment for HSV/VZV infection during allogeneic HSCT or high-intensity chemotherapy, developed in NONMEM 7.4 from 374 therapeutic-drug-monitoring plasma concentrations (paired peak and trough samples on up to three occasions). Clearance scales allometrically with body weight (fixed exponent 0.75, reference 27.8 kg) and as a power function of Schwartz eGFR (reference 209.4 mL/min/1.73 m^2); volume scales linearly with body weight. Inter-individual variability on CL and V, inter-occasion variability on CL (three occasions), and proportional residual error."
  reference <- paste(
    "Maximova N, Nistico D, Luci G, Simeone R, Piscianz E, Segat L, Barbi E,",
    "Di Paolo A. (2022).",
    "Population Pharmacokinetics of Intravenous Acyclovir in Oncologic",
    "Pediatric Patients.",
    "Front Pharmacol 13:865871.",
    "doi:10.3389/fphar.2022.865871"
  )
  vignette <- "Maximova_2022_acyclovir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "acyclovir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling on CL (fixed exponent 0.75) and V (fixed exponent 1), both referenced to 27.8 kg, the cohort median (Maximova 2022 Table 1; Table 2 footnote equations). The paper reports body weight per occasion (32.2, 30.9 and 32.8 kg mean at occasions 1-3), so WT may be supplied time-varying.",
      source_name = "WGT"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate by the original Schwartz (1976) formula, BSA-normalised.",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Schwartz formula eGFR (Maximova 2022 Methods, Study Design and Population: 'The Schwartz formula determined the estimated glomerular filtration rate (eGFR) for each patient (Schwartz et al., 1976)'). Enters CL as the power term (CRCL / 209.4)^1.627, centred on the cohort median 209.4 mL/min/1.73 m^2 (Table 1; mean 228.1 +/- 79.7). The paper's simulations stratify at 250 mL/min/1.73 m^2 for augmented renal clearance.",
      source_name = "eGFR"
    ),
    OCC = list(
      description = "Therapeutic-drug-monitoring occasion (1, 2 or 3).",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Each occasion was one peak/trough sample pair drawn at steady state; 120, 54 and 13 patients contributed a first, second and third occasion (Maximova 2022 POP/PK Modeling). Decomposed inside model() into indicators that select one of three inter-occasion etas on CL. Values outside 1-3 switch IOV off (typical-occasion simulation).",
      source_name = "OCC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 120L,
    n_studies = 1L,
    age_range = "0-18 years (inclusion criterion)",
    age_mean = "8.6 +/- 5.0 years (median 9.5)",
    weight_mean = "32.4 +/- 19.1 kg",
    weight_median = "27.8 kg",
    sex_female_pct = 39.2,
    disease_state = "Children with haematological malignancies undergoing allogeneic HSCT or high-intensity chemotherapy, receiving IV acyclovir for HSV/VZV prophylaxis (n = 94) or treatment (n = 26).",
    dose_range = "IV acyclovir every 6-8 h as a 60-min infusion; median starting daily dose 40.7 (range 15.6-136.7) mg/kg/day; first-occasion dose 334.0 +/- 158.2 mg (median 350 mg). Doses adjusted by TDM to keep Cmin > 0.5 mg/L and Cmax < 25 mg/L.",
    renal_function = "Schwartz eGFR 228.1 +/- 79.7 mL/min/1.73 m^2 (median 209.4); serum creatinine 0.361 +/- 0.206 mg/dL.",
    regions = "Single centre, IRCCS Burlo Garofolo, Trieste, Italy (2011-2020).",
    sampling_design = "374 plasma concentrations: 187 peak/trough pairs (30 min after the end of infusion and 10 min before the next dose) on 1-3 steady-state occasions (120, 54 and 13 patients). LC-MS/MS assay.",
    notes = "Demographics from Maximova 2022 Table 1 and Results (73 boys, 47 girls). Height 128.2 +/- 30.0 cm; BSA 1.058 +/- 0.426 m^2."
  )

  ini({
    # Structural parameters -- Maximova 2022 Table 2 'Final model' column.
    # Reference child: 27.8 kg, Schwartz eGFR 209.4 mL/min/1.73 m^2.
    lcl <- log(6.184); label("Clearance for a 27.8 kg child at eGFR 209.4 mL/min/1.73 m^2 (L/h)") # Table 2 CL = 6.184 L/h (SE 1.630)
    lvc <- log(18.942); label("Volume of distribution for a 27.8 kg child (L)") # Table 2 V = 18.942 L (SE 2.135)

    # Allometric exponents: printed in the Table 2 footnote equations
    # with no estimate, SE or bootstrap interval, so held fixed.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT on CL (unitless)") # Table 2 footnote (WGT/27.8)^0.75
    e_wt_vc <- fixed(1); label("Allometric exponent of WT on V (unitless)") # Table 2 footnote (WGT/27.8)^1

    # Power effect of Schwartz eGFR on CL.
    e_crcl_cl <- 1.627; label("Power exponent of eGFR on CL, (CRCL/209.4)^e_crcl_cl (unitless)") # Table 2 EGFR on CL = 1.627 (SE 0.269)

    # Inter-individual variability (log-normal), Table 2 reported as CV%;
    # omega^2 = log(1 + CV^2).
    etalcl ~ 0.18596 # Table 2 IIV CL = 45.2% CV -> log(1 + 0.452^2)
    etalvc ~ 0.27878 # Table 2 IIV V = 56.7% CV -> log(1 + 0.567^2)

    # Inter-occasion variability on CL, Table 2 IOV CL = 20.0% CV ->
    # log(1 + 0.20^2) = 0.03922; one eta per occasion with a shared variance.
    etaiov_cl_1 ~ 0.03922 # Table 2 IOV CL = 20.0% CV, occasion 1
    etaiov_cl_2 ~ fixed(0.03922) # Table 2 IOV CL = 20.0% CV, occasion 2 (shared variance)
    etaiov_cl_3 ~ fixed(0.03922) # Table 2 IOV CL = 20.0% CV, occasion 3 (shared variance)

    # Proportional residual error. Table 2 'Residual variability' = 0.181 is
    # read as the NONMEM SIGMA variance (SD = sqrt(0.181) = 0.4254); the
    # pcVPC band widths of Figure 4 match this reading (see vignette).
    propSd <- 0.4254; label("Proportional residual SD (fraction)") # Table 2 Residual variability = 0.181 (variance)
  })

  model({
    # Occasion indicators for the inter-occasion eta on CL.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3

    # Individual parameters (Table 2 footnote equations, with the printed
    # '+' signs read as multiplication and the eta superscripts as
    # exponentials -- see vignette Assumptions).
    cl <- exp(lcl + etalcl + iov_cl) * (CRCL / 209.4)^e_crcl_cl * (WT / 27.8)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 27.8)^e_wt_vc

    kel <- cl / vc
    d/dt(central) <- -kel * central

    # Plasma concentration; dose in mg and V in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
