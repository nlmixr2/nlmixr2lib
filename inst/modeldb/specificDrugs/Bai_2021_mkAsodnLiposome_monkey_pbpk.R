Bai_2021_mkAsodnLiposome_monkey_pbpk <- function() {
  description <- paste(
    "Preclinical (macaque monkey). PBPK-derived reduced one-compartment",
    "intravenous model for MK-ASODN nanoliposomes, a nanoliposome-encapsulated",
    "20-mer antisense oligonucleotide (5'-CCCCGGGCCGCCCTTCTTCA, 6044.4 Da)",
    "against midkine mRNA, developed for hepatocellular carcinoma. The source",
    "paper built a 14-tissue perfusion-limited whole-body PBPK model in",
    "GastroPlus 8.0 using the software's default monkey physiology and",
    "Rodgers-Single (Lukacova) tissue-to-plasma partition coefficients; those",
    "organ volumes, blood flows and Kp values are platform internals that the",
    "publication does not print, so the whole-body structure is not reproduced",
    "here. What this file encodes is the plasma-level behaviour of that model:",
    "a single linear clearance of 0.1002 L/h/kg, which reproduces the paper's",
    "own predicted AUC0-inf on all three dose arms (11.5, 23, 46 mg/kg) to four",
    "significant figures, and a volume of 0.1326 L/kg from the terminal slope",
    "of the simulated curves in Figure 1, whose back-extrapolated intercept",
    "reproduces the paper's own tabulated predicted Cmax. The PBPK curve has a",
    "brief initial spike (under about 10 minutes) that a one-compartment",
    "reduction does not carry; the spike holds under 2% of the AUC. The PBPK",
    "model is linear, so it does not carry the less-than-dose-proportional",
    "exposure the observed monkey NCA showed. Clearance and volume scale",
    "linearly with body weight because the paper reports them per kg. The",
    "paper reports no interindividual or residual variability for the",
    "simulation, so no etas are declared and both residual error terms are",
    "fixed(0).",
    sep = " "
  )
  reference <- paste(
    "Bai H, Cheng Y, Che J. Pharmacokinetics and Disposition of Heparin-Binding",
    "Growth Factor Midkine Antisense Oligonucleotide Nanoliposomes in",
    "Experimental Animal Species and Prediction of Human Pharmacokinetics Using",
    "a Physiologically Based Pharmacokinetic Model.",
    "Front Pharmacol. 2021;12:769538. doi:10.3389/fphar.2021.769538.",
    "PMCID PMC8595129.",
    "Monkey study design is Methods 2.2 and 2.4; the PBPK method is Methods",
    "2.7; observed monkey NCA is Table 1; observed and PBPK-predicted Cmax and",
    "AUC0-inf are Table 4; observed and simulated monkey plasma profiles are",
    "Figure 1 (panels A-C).",
    sep = " "
  )
  vignette <- "Bai_2021_mkAsodnLiposome_pbpk"
  # The paper reports concentrations in ug/mL and AUC in ug.h/mL; mg/L is the
  # numerically identical form.
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales clearance and volume linearly (exponent 1), because the paper",
        "reports monkey clearance per kg (Table 1, mL/kg/h) and doses per kg.",
        "Reference 6 kg is the approximate monkey weight in Methods 2.2 ('each",
        "weighing approximately 6 kg'). With a mg/kg dose, concentrations are",
        "independent of the weight chosen.",
        sep = " "
      ),
      source_name = "BW"
    )
  )

  compartmentData <- list(
    central = list(analyte = "MK-ASODN", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "macaque monkey",
    n_subjects = 9,
    n_studies = 1,
    age_range = NA_character_,
    weight_range = "approximately 6 kg each (Methods 2.2)",
    sex_female_pct = 44.4,
    race_ethnicity = NA_character_,
    disease_state = "Healthy animals",
    dose_range = paste(
      "Single intravenous injection of 11.5, 23 or 46 mg/kg, 3 monkeys per",
      "dose group (Methods 2.4).",
      sep = " "
    ),
    regions = "China (Beijing)",
    notes = paste(
      "Five males and four females (Methods 2.2). Plasma sampled pre-dose and",
      "at 0.08-6 h post-dose (Methods 2.4). The PBPK simulation represents a",
      "single typical monkey using GastroPlus default monkey physiology.",
      sep = " "
    )
  )

  ini({
    # Clearance. Methods 2.7 states the monkey PBPK clearance input was the in
    # vivo clearance, but does not print the value used. It is recovered
    # exactly from Table 4: dose / predicted AUC0-inf = 11.5 / 114.77 =
    # 0.10020, 23 / 229.54 = 0.10020 and 46 / 459.09 = 0.10020 L/h/kg. The
    # three agree to five significant figures, which also proves the PBPK
    # model's clearance is dose-independent (linear). 0.1002 L/h/kg x 6 kg =
    # 0.6012 L/h.
    lcl <- fixed(log(0.6012)) ; label("Clearance for a 6 kg monkey (L/h)")  # Table 4, dose / predicted AUC0-inf = 0.1002 L/h/kg on all three arms, x 6 kg

    # Volume. Not printed (the paper gives only the observed NCA Vss in Table
    # 1). Digitised by the maintainers from the simulated (solid-line) curves
    # of Figure 1: a log-linear regression over 1.7-5.9 h gives k = 0.7557 /h
    # (panel A, 11.5 mg/kg) and 0.7555 /h (panel B, 23 mg/kg), so V = CL / k =
    # 0.1326 L/kg on both. The intercepts independently give dose / C0 =
    # 0.1327 and 0.1335 L/kg. Panel C (46 mg/kg; coarser axis) gives 0.1341.
    # 0.1326 L/kg x 6 kg = 0.7956 L.
    lvc <- fixed(log(0.7956)) ; label("Central volume of distribution for a 6 kg monkey (L)")  # Figure 1A-B simulated terminal slope, V = CL / k = 0.1326 L/kg, x 6 kg

    # Residual error. The paper simulated a single typical monkey and reports
    # no variability model, so both terms are fixed at zero rather than
    # invented.
    addSd <- fixed(0) ; label("Additive residual error (mg/L; not published by Bai 2021)")
    propSd <- fixed(0) ; label("Proportional residual error (fraction; not published by Bai 2021)")
  })

  model({
    # Linear per-kg scaling (Table 1 reports CL per kg); reference 6 kg
    # (Methods 2.2).
    cl <- exp(lcl) * (WT / 6)
    vc <- exp(lvc) * (WT / 6)

    kel <- cl / vc

    # Intravenous injection only (Methods 2.4); no absorption model.
    d/dt(central) <- -kel * central

    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
