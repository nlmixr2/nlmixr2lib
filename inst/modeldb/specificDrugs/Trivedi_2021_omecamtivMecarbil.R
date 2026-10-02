Trivedi_2021_omecamtivMecarbil <- function() {
  description <- paste(
    "Two-compartment oral population pharmacokinetic reduction of the",
    "Simcyp minimal-PBPK-with-single-adjusting-compartment (SAC) compound",
    "model for the cardiac myosin activator omecamtiv mecarbil in healthy",
    "adults (Trivedi 2021). Distribution was fitted by the authors as an",
    "ordinary two-compartment model to intravenous data (k12, k21 and the",
    "central volume Vc, with Vsac = Vc * k12 / k21) and then held fixed",
    "while first-order absorption ka and the apparent oral clearance",
    "CLpo were fitted to oral single- and multiple-dose data, so the",
    "Simcyp compound layer IS a two-compartment model with first-order",
    "absorption and needs no platform physiology to encode. Volumes are",
    "per kg body weight (Vss 3.8575 L/kg, Vsac 2.369 L/kg, systemic",
    "compartment Vss - Vsac); clearance is not weight-scaled. Inter-",
    "individual variability is the compound file's 30% CV on ka, CLpo and",
    "Vss (Supplementary Table S1 workbook), with Vsac held fixed as",
    "Simcyp does, so all Vss variability falls on the systemic volume.",
    "Elimination uses the apparent oral clearance with bioavailability 1;",
    "the typical-value reduction reproduces the medians of the 50",
    "deposited Simcyp virtual subjects (Cmax, Tmax, AUC0-48 after 50 mg)",
    "within 7.5%. The paper's main result, rosuvastatin (BCRP substrate)",
    "exposure with omecamtiv mecarbil co-administration, runs on a",
    "Simcyp library full-PBPK/ADAM rosuvastatin compound file with",
    "transporter kinetics and is NOT reproducible from this model.",
    sep = " "
  )
  reference <- paste(
    "Trivedi A, Sohn W, Kulkarni P, Jafarinasabian P, Zhang H, Spring M,",
    "Flach S, Abbasi S, Wahlstrom J, Lee E, Dutta S. (2021). Evaluation",
    "of drug-drug interaction potential between omecamtiv mecarbil and",
    "rosuvastatin, a BCRP substrate, with a clinical study in healthy",
    "subjects and using a physiologically-based pharmacokinetic model.",
    "Clin Transl Sci 14(6):2510-2520. doi:10.1111/cts.13118.",
    sep = " "
  )
  vignette <- "Trivedi_2021_omecamtivMecarbil"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear (exponent 1) scaling of both distribution volumes: Simcyp",
        "takes Vss and Vsac in L/kg (Trivedi 2021 Table 3; Supplementary",
        "Table S1 workbook 'Input Sheet'). Clearance is entered in L/h and",
        "is not weight-scaled. Table 3 converts the fitted litre values",
        "with a mean body weight of 83.79 kg; the deposited Simcyp",
        "Sim-Healthy Volunteers reference body weight is 80.7 kg."
      ),
      source_name = "body weight (Simcyp virtual individual)"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "omecamtiv mecarbil",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "omecamtiv mecarbil",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "omecamtiv mecarbil",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 81L,
    n_studies = 6L,
    age_range = "18-54 years across the six clinical studies (Table 1)",
    weight_mean = "83.79 kg (mean weight used to express the fitted volumes per kg, Table 3)",
    sex_female_pct = 34,
    race_ethnicity = paste(
      "Rosuvastatin DDI study (n = 14): 71.4% White, 21.4% Black, 7.1%",
      "multiple race, no Asian participants. Not reported for the other",
      "studies."
    ),
    disease_state = "Healthy adult volunteers.",
    dose_range = paste(
      "35 mg single intravenous dose (study 5, distribution fit); 50 mg",
      "single oral dose (study 1 and the rosuvastatin DDI study); 25 mg",
      "and 50 mg oral twice daily for 7-14 days (studies 2-4)."
    ),
    regions = "Not reported; Simcyp Sim-Healthy Volunteers virtual population.",
    notes = paste(
      "n_subjects sums the participants of the six studies in Trivedi",
      "2021 Table 1 (14 + 20 + 13 + 13 + 7 + 14); sex_female_pct is the",
      "participant-weighted proportion of women from the same table.",
      "Distribution parameters were fitted to study 5 (35 mg IV, n = 7)",
      "with SAAM II; ka and CLpo were fitted in Simcyp to studies 1 and 2;",
      "studies 3 and 4 were verification only. This is a PBPK analysis,",
      "not a pooled population-PK fit: the IIV below is the Simcyp",
      "compound-file default of 30% CV, not an estimated omega."
    )
  )

  ini({
    # Values are the as-run Simcyp V17 compound-file inputs for
    # omecamtiv mecarbil ('AMG 423_target formulation', Inhibitor 1) in
    # the Supplementary Table S1 workbook, sheet 'Input Sheet'. Each agrees
    # with the rounded Trivedi 2021 Table 3 value quoted alongside it.

    lka <- log(0.2234)
    label("First-order absorption rate constant ka (1/h)")
    # S1 workbook 'ka (1/h)' = 0.2234; Table 3 ka = 0.220 'Parameterized';
    # Results: 0.22 (95% CI 0.18 to 0.26) 1/h, fitted in Simcyp to studies 1 and 2.

    lcl <- log(10.81)
    label("Apparent oral clearance CLpo (CL/F, L/h)")
    # S1 workbook 'CL (po) (L/h)' = 10.81; Table 3 CLpo = 10.8 L/h 'Parameterized';
    # Results: 10.8 (95% CI 10.3 to 11.3) L/h. Not weight-scaled.

    lvss <- fixed(log(3.8575))
    label("Steady-state volume of distribution Vss per kg body weight (L/kg)")
    # S1 workbook 'Vss (L/kg)' = 3.8575; Table 3 Vss = 3.86 L/kg 'Calculated
    # (Vss = 323.3 L; mean weight 83.79 kg)'. Fixed from the IV compartmental fit.

    lvp <- fixed(log(2.369))
    label("Single adjusting compartment volume Vsac per kg body weight (L/kg)")
    # S1 workbook 'Volume [Vsac] (L/kg)' = 2.369; Table 3 Vsac = 2.37 L/kg
    # 'Calculated (Vsac = 198.5 L; mean weight 83.79 kg)'. Fixed from the IV fit.

    lk12 <- fixed(log(0.31955))
    label("Systemic-to-SAC transfer rate constant k12 (Simcyp SAC kin, 1/h)")
    # S1 workbook 'SAC kin (1/h)' = 0.31955; Table 3 k12 = 0.319 1/h; Results:
    # 0.319 (95% CI 0.186 to 0.452) 1/h from the SAAM II fit to 35 mg IV (study 5).

    lk21 <- fixed(log(0.2))
    label("SAC-to-systemic transfer rate constant k21 (Simcyp SAC kout, 1/h)")
    # S1 workbook 'SAC kout (1/h)' = 0.2; Table 3 k21 = 0.200 1/h; Results:
    # 0.200 (95% CI 0.136 to 0.265) 1/h from the SAAM II fit to 35 mg IV.

    # Simcyp compound-file variability, S1 workbook 'CV ka (%)' = 30,
    # 'CV CL (po) (%)' = 30 and 'CV Vss (%)' = 30. Encoded log-normally:
    # omega^2 = log(1 + 0.30^2) = 0.0861777. Vsac, k12 and k21 carry no CV
    # in the compound file. These are platform defaults, not estimates.
    etalka ~ fixed(0.0861777)
    etalcl ~ fixed(0.0861777)
    etalvss ~ fixed(0.0861777)

    # A Simcyp PBPK simulation has no residual-error model; zero keeps
    # this a simulation model without inventing a variance.
    propSd <- fixed(0)
    label("Proportional residual error SD (fraction; zero, no error model in the source)")
  })

  model({
    # 1. Individual parameters. Volumes scale linearly with body weight
    #    because Simcyp takes them in L/kg; CLpo is in L/h.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vss <- exp(lvss + etalvss)
    vsac <- exp(lvp)
    k12 <- exp(lk12)
    k21 <- exp(lk21)

    # 2. Systemic (central) compartment volume per kg. Simcyp samples Vss
    #    with a 30% CV but holds Vsac fixed, so Vsys = Vss - Vsac (typical
    #    1.4885 L/kg; Table 3 Vc = 1.49 L/kg). For the ~5% of draws with
    #    Vss below Vsac the platform keeps a small positive systemic volume
    #    instead; the 0.035 L/kg floor is the smallest systemic volume
    #    among the 50 virtual subjects of the S1 workbook ('Distribution -
    #    Vols', Inhibitor 1 Vsys minimum 0.0347 L/kg). The floor is a
    #    maintainer approximation of an unpublished platform rule and does
    #    not act at the typical value.
    vckg <- vss - vsac
    if (vckg < 0.035) {
      vckg <- 0.035
    }
    vc <- vckg * WT
    vp <- (vss - vckg) * WT

    # 3. Elimination micro-constant from the apparent oral clearance.
    #    With bioavailability 1 the oral AUC(0,inf) equals Dose / CLpo,
    #    which is how Simcyp's in vivo CLpo input is defined.
    kel <- cl / vc

    # 4. ODEs: first-order absorption (Simcyp fa = 1, no lag), systemic
    #    compartment and the SAC exchanging by the mass-based k12 / k21 of
    #    the SAAM II intravenous fit. Amounts in mg.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Observation: mg / L = ug/mL; x 1000 for ng/mL (Tables 4-5 units).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
