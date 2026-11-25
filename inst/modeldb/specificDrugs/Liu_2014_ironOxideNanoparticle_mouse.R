Liu_2014_ironOxideNanoparticle_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, male athymic nude). Two-compartment (blood and liver)",
    "kinetic model of a lipid-coated superparamagnetic iron oxide (SPIO)",
    "nanoparticle after an IV bolus, fitted per animal to dynamic T2*-weighted MRI",
    "signal; separate rate-constant sets for 0.5, 3 and 10 mg Fe/kg and for",
    "clodronate-liposome or blank-liposome pretreatment before 10 mg Fe/kg.",
    "RECORD OF THE PUBLICATION ONLY, NOT USABLE FOR PREDICTION: the printed",
    "Table 3 rate constants cannot reproduce the paper's own fitted curve",
    "(Figure 1c) for any blood volume, and the blood volume and input fraction are",
    "unreported; see the vignette Errata."
  )
  reference <- paste(
    "Liu T, Choi H, Zhou R, Chen I-W. (2014). Quantitative evaluation of the",
    "reticuloendothelial system function with dynamic MRI. PLoS One 9(8):e103576.",
    "doi:10.1371/journal.pone.0103576."
  )
  vignette <- "Liu_2014_ironOxideNanoparticle"
  units <- list(time = "s", dosing = "mg/kg", concentration = "mg/L")

  covariateData <- list(
    DOSE_IRONOXIDENANOPARTICLE_MGKG = list(
      description = "SPIO nanoparticle IV bolus dose level of the animal's group, as iron",
      units = "mg Fe/kg",
      type = "continuous",
      reference_category = "10 mg Fe/kg (any value other than 0.5 or 3 takes the 10 mg Fe/kg rate constants)",
      notes = paste(
        "Selects the per-dose-group rate constants of Table 3 (0.5, 3 and 10 mg Fe/kg).",
        "The paper reports one parameter set per dose group, not a dose-exposure",
        "function, so the levels are categorical. Needed as a covariate because the",
        "parameter set cannot be recovered from the dose records."
      ),
      source_name = "dose group (Table 3 column)"
    ),
    CONMED_CLODRONATE = list(
      description = "Clodronate-liposome (clodrosome) macrophage-depletion pretreatment indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no clodronate liposome (saline)",
      notes = paste(
        "1 = clodrosome 50 mg/kg IV once daily for 2 days before the 10 mg Fe/kg SPIO",
        "dose (Methods 'Animals'). Selects the 'Clodrosome +10 mg Fe/kg' column of",
        "Table 3. Takes precedence over the dose level."
      ),
      source_name = "Clodrosome +10 mg Fe/kg (Table 3 column)"
    ),
    CONMED_PHOSPHATIDYLCHOLINE = list(
      description = "Blank phosphatidylcholine-liposome RES-blockade pretreatment indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no liposome pretreatment (saline)",
      notes = paste(
        "1 = 376 mg L-alpha-phosphatidylcholine/kg encapsulated liposome IV 1.5 h",
        "before the 10 mg Fe/kg SPIO dose (Methods 'Animals'). Selects the 'liposome",
        "+10 mg Fe/kg' column of Table 3. Takes precedence over the dose level."
      ),
      source_name = "liposome +10 mg Fe/kg (Table 3 column)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "SPIO nanoparticle (iron)", units = "mg/kg", specimen = "whole blood", verified = TRUE),
    liver = list(analyte = "SPIO nanoparticle (iron)", units = "mg/kg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "mouse (male athymic nude, 01b74)",
    n_subjects = 18L,
    n_studies = 1L,
    age_range = "young adult (not further specified)",
    sex_female_pct = 0,
    disease_state = "healthy; macrophage-depleted (clodronate liposome) or RES-blocked (blank liposome) subgroups",
    dose_range = "single IV bolus of 0.5, 3 or 10 mg Fe/kg SPIO nanoparticles over about 25 s, with a 100 uL saline flush",
    notes = paste(
      "24 mice in total. Kinetic modelling: 10 untreated mice across 0.5, 3 and 10 mg",
      "Fe/kg (one 10 mg Fe/kg mouse excluded for a high Ka; Figure 4 shows 3 animals",
      "per dose group), plus 4 clodrosome-pretreated and 4 liposome-pretreated mice at",
      "10 mg Fe/kg. The remaining 6 mice were used for ex vivo optical biodistribution.",
      "Body weight is not reported. Concentrations are delta-R2*/r2* from vena cava",
      "and liver regions of interest at 9.4 T, in mg Fe/L. Each animal was fitted",
      "separately by least squares in MATLAB SimBiology; Table 3 gives the group",
      "mean +/- SD of the per-animal estimates."
    )
  )

  ini({
    # Reference group: untreated mice at 10 mg Fe/kg (Table 3, which tests every
    # other column against this one). Rate constants are in 1/s.
    lk_central_liver <- log(378027.6); label("Blood-to-liver association rate constant Kin, 10 mg Fe/kg (1/s)") # Table 3: Kin 10 mg Fe/kg = 378027.6 +/- 9504.385 s^-1
    lk_liver_central <- log(0.365); label("Liver-to-blood dissociation rate constant Kout, 10 mg Fe/kg (1/s)") # Table 3: Kout 10 mg Fe/kg = 0.365 +/- 0.01 s^-1
    lkel <- log(0.025); label("Extrahepatic (blood) elimination rate constant Ke, 10 mg Fe/kg (1/s)") # Table 3: Ke 10 mg Fe/kg = 0.025 +/- 0.001 s^-1

    # Volumes. VL is fixed by the authors. Vp is named in the Methods and the
    # equations but no value is printed anywhere; 0.1 L/kg is the maintainers'
    # reading of Figure 1c (Cp about 100 mg/L just after 10 mg Fe/kg, so
    # Vp = Dose / Cp0). The printed rate constants are volume-independent in
    # their failure to reproduce Figure 1c (see the vignette Errata).
    lvc <- fixed(log(0.1)); label("Blood distribution volume Vp (L/kg)") # Not printed; maintainers' digitisation of Figure 1c Cp0 ~100 mg/L at 10 mg Fe/kg -> 10/100 L/kg
    lv_liver <- fixed(log(1e-3)); label("Liver compartment distribution volume VL (L/kg)") # Methods 'Pharmacokinetic Modeling': VL 'fixed at the reasonable value (1x10-3 L/kg)'

    # Group effects as log-ratios to the 10 mg Fe/kg column, written so that each
    # group's rate constant reproduces the printed Table 3 mean exactly.
    e_dose05_k_central_liver <- log(38.527 / 378027.6); label("Log-ratio of Kin, 0.5 vs 10 mg Fe/kg (unitless)") # Table 3: Kin 0.5 mg Fe/kg = 38.527 +/- 3.512 s^-1
    e_dose05_k_liver_central <- log(0.014 / 0.365); label("Log-ratio of Kout, 0.5 vs 10 mg Fe/kg (unitless)") # Table 3: Kout 0.5 mg Fe/kg = 0.014 +/- 0.001 s^-1
    e_dose05_kel <- log(0.008 / 0.025); label("Log-ratio of Ke, 0.5 vs 10 mg Fe/kg (unitless)") # Table 3: Ke 0.5 mg Fe/kg = 0.008 +/- 0.001 s^-1

    e_dose3_k_central_liver <- log(1750.480 / 378027.6); label("Log-ratio of Kin, 3 vs 10 mg Fe/kg (unitless)") # Table 3: Kin 3 mg Fe/kg = 1750.480 +/- 100 s^-1
    e_dose3_k_liver_central <- log(0.098 / 0.365); label("Log-ratio of Kout, 3 vs 10 mg Fe/kg (unitless)") # Table 3: Kout 3 mg Fe/kg = 0.098 +/- 0.01 s^-1
    e_dose3_kel <- log(0.015 / 0.025); label("Log-ratio of Ke, 3 vs 10 mg Fe/kg (unitless)") # Table 3: Ke 3 mg Fe/kg = 0.015 +/- 0.001 s^-1

    e_conmed_clodronate_k_central_liver <- log(4750.368 / 378027.6); label("Log-ratio of Kin, clodrosome + 10 vs 10 mg Fe/kg (unitless)") # Table 3: Kin Clodrosome +10 mg Fe/kg = 4750.368 +/- 83.616 s^-1
    e_conmed_clodronate_k_liver_central <- log(0.128 / 0.365); label("Log-ratio of Kout, clodrosome + 10 vs 10 mg Fe/kg (unitless)") # Table 3: Kout Clodrosome +10 mg Fe/kg = 0.128 +/- 0.056 s^-1
    e_conmed_clodronate_kel <- log(0.015 / 0.025); label("Log-ratio of Ke, clodrosome + 10 vs 10 mg Fe/kg (unitless)") # Table 3: Ke Clodrosome +10 mg Fe/kg = 0.015 +/- 0.001 s^-1

    e_conmed_phosphatidylcholine_k_central_liver <- log(1770.913 / 378027.6); label("Log-ratio of Kin, liposome + 10 vs 10 mg Fe/kg (unitless)") # Table 3: Kin liposome +10 mg Fe/kg = 1770.913 +/- 73.476 s^-1
    e_conmed_phosphatidylcholine_k_liver_central <- log(5.875 / 0.365); label("Log-ratio of Kout, liposome + 10 vs 10 mg Fe/kg (unitless)") # Table 3: Kout liposome +10 mg Fe/kg = 5.875 +/- 0.783 s^-1
    e_conmed_phosphatidylcholine_kel <- log(0.015 / 0.025); label("Log-ratio of Ke, liposome + 10 vs 10 mg Fe/kg (unitless)") # Table 3: Ke liposome +10 mg Fe/kg = 0.015 +/- 0.001 s^-1
  })
  model({
    # Group indicators. A pretreatment indicator overrides the dose level, since
    # both pretreatment groups received 10 mg Fe/kg (Methods 'Animals').
    pretreated <- 0
    if (CONMED_CLODRONATE == 1) pretreated <- 1
    if (CONMED_PHOSPHATIDYLCHOLINE == 1) pretreated <- 1
    dose05 <- 0
    if (DOSE_IRONOXIDENANOPARTICLE_MGKG == 0.5) dose05 <- 1
    dose3 <- 0
    if (DOSE_IRONOXIDENANOPARTICLE_MGKG == 3) dose3 <- 1
    dose05 <- dose05 * (1 - pretreated)
    dose3 <- dose3 * (1 - pretreated)
    clod <- CONMED_CLODRONATE
    lipo <- CONMED_PHOSPHATIDYLCHOLINE * (1 - CONMED_CLODRONATE)

    k_central_liver <- exp(lk_central_liver + e_dose05_k_central_liver * dose05 +
      e_dose3_k_central_liver * dose3 + e_conmed_clodronate_k_central_liver * clod +
      e_conmed_phosphatidylcholine_k_central_liver * lipo)
    k_liver_central <- exp(lk_liver_central + e_dose05_k_liver_central * dose05 +
      e_dose3_k_liver_central * dose3 + e_conmed_clodronate_k_liver_central * clod +
      e_conmed_phosphatidylcholine_k_liver_central * lipo)
    kel <- exp(lkel + e_dose05_kel * dose05 + e_dose3_kel * dose3 +
      e_conmed_clodronate_kel * clod + e_conmed_phosphatidylcholine_kel * lipo)
    vc <- exp(lvc)
    v_liver <- exp(lv_liver)

    # Equations e001-e002 are written on amounts: d(Cp*Vp)/dt and d(CL*VL)/dt,
    # with transfer fluxes Kin*Cp*Vp = Kin*central and Kout*CL*VL = Kout*liver.
    # The input term Ka*Dose is an IV bolus scaled by an input fraction Ka that
    # the authors say 'should be approximately 1' and never report; the bolus
    # dose enters central unscaled (Ka = 1).
    d/dt(central) <- -kel * central - k_central_liver * central + k_liver_central * liver
    d/dt(liver) <- k_central_liver * central - k_liver_central * liver

    # Blood (vena cava) and liver SPIO concentrations, mg Fe/L.
    Cc <- central / vc
    Cliver <- liver / v_liver
  })
}
