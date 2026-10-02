Agbo_2021_apomorphine <- function() {
  description <- "Joint parent-metabolite population PK model for apomorphine and apomorphine sulfate after apomorphine sublingual film (APL) or subcutaneous apomorphine in healthy subjects and patients with Parkinson's disease and OFF episodes (Agbo 2021). Apomorphine: two-compartment disposition, sublingual and subcutaneous absorption each through two transit compartments, all apomorphine elimination by first-order metabolism to apomorphine sulfate. Apomorphine sulfate: one-compartment disposition with first-order elimination, plus direct first-order input of the swallowed fraction of the sublingual dose. Sublingual bioavailability relative to subcutaneous (about 18%) decreases with sublingual dose (power, reference 20 mg); sublingual absorption rate decreases with contact time under the tongue (power, reference 2 min); apomorphine central volume increases with body weight (power, reference 69.3 kg) and is lower in study CTH-103; apomorphine sulfate volume is lower in women."
  reference <- "Agbo F, Crass RL, Chiu YY, Chapel S, Galluppi G, Blum D, Navia B. Population pharmacokinetic analysis of apomorphine sublingual film or subcutaneous apomorphine in healthy subjects and patients with Parkinson's disease. Clin Transl Sci. 2021;14(4):1464-1475. doi:10.1111/cts.13008"
  vignette <- "Agbo_2021_apomorphine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on apomorphine central volume V2/F, reference 69.3 kg (Results 'Final model'; Table 3 'Body weight on V2/F' = 1.53). Apparent clearance is NOT weight-dependent, so a heavier subject has a larger V2/F, a smaller k23 = CL/V2 and a lower Cmax and early AUC.",
      source_name = "WT"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "male (SEXF = 0)",
      notes = "Proportional shift on apomorphine sulfate volume V3/F: V3/F * (1 + e_sexf_vc_sulf * SEXF) (Table 3 'Female sex on V3/F' = -0.310; Results 'Final model': 'represent the proportional shift from the reference condition'; Methods Eq. 5). No effect on apomorphine itself.",
      source_name = "SEX"
    ),
    STUDY_CTH103 = list(
      description = "Study CTH-103 indicator (1 = subject enrolled in phase I study CTH-103, 0 = any other study)",
      units = "(binary)",
      type = "binary",
      reference_category = "other studies (STUDY_CTH103 = 0)",
      notes = "Proportional shift on apomorphine central volume V2/F: V2/F * (1 + e_study_cth103_vc * STUDY_CTH103) (Table 3 'Study CTH-103 on V2/F' = -0.555, Methods Eq. 5). Added at the base-model stage to account for very high peak concentrations observed in CTH-103, particularly after subcutaneous dosing (Results 'Base model'). Set to 0 for simulation of new subjects.",
      source_name = "STUDY"
    ),
    DOSE_APOMORPHINE_SL_MG = list(
      description = "Administered apomorphine sublingual film dose (mg) of the current administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the sublingual bioavailability F1, reference 20 mg (Table 3 'Dose of sublingual apomorphine on F1' = -0.206; Results: 'the power to which the ratio of the covariate value to the reference value (20 mg ...) is raised'). Supplied as a data column because model code cannot read the dose record's amt. Only scales f(depot), so its value on subcutaneous records has no effect; set it to the sublingual film dose on sublingual records.",
      source_name = "DOSE"
    ),
    DUR_SL_CONTACT = list(
      description = "Contact time of the apomorphine sublingual film under the tongue (min)",
      units = "min",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the sublingual absorption rate constant, reference 2 min (Table 3 'Contact time under the tongue for sublingual film on ka for sublingual administration' = -0.194; Results 'Final model'). Only scales the sublingual chain, so its value on subcutaneous records has no effect. Figure 2's typical patient uses 3 min.",
      source_name = "contact time under the tongue"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "apomorphine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "apomorphine", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "apomorphine", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "apomorphine", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "apomorphine", units = "mg", specimen = "administration site", verified = TRUE),
    transit4 = list(analyte = "apomorphine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "apomorphine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "apomorphine", units = "mg", specimen = "tissue", verified = TRUE),
    depot_sulf = list(
      analyte = "apomorphine (swallowed fraction of the sublingual dose, absorbed as apomorphine sulfate)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central_sulf = list(
      analyte = "apomorphine sulfate (apomorphine-equivalent amount)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 158L,
    n_studies = 9L,
    age_range = "21.0-88.1 years (healthy subjects 21.0-50.0; patients with PD 45.6-88.1)",
    age_mean = "27.2 years (healthy subjects); 65.0 years (patients with PD)",
    weight_range = "50.3-146.1 kg (healthy subjects 50.3-87.8; patients with PD 50.7-146.1)",
    weight_mean = "65.8 kg (healthy subjects); 80.0 kg (patients with PD)",
    weight_reference = "69.3 kg (approximate population median; reference for the V2/F weight effect)",
    sex_female_pct = 27.8,
    race_ethnicity = c(White = 42.4, Black = 1.9, Asian = 53.8, Other = 1.9),
    disease_state = "87 healthy subjects (predominantly Asian) and 71 patients with Parkinson's disease and OFF episodes (predominantly White)",
    dose_range = "Apomorphine sublingual film 10-35 mg and 50 mg (exploratory formulation K and commercial formulation O); subcutaneous apomorphine 2-5 mg (APOKYN or APO-go)",
    regions = "Not reported by region in the source",
    renal_function = "Creatinine clearance 86-155 mL/min (healthy subjects), 52-213 mL/min (patients with PD); 26 patients with PD had mild renal impairment",
    n_observations = "2485 apomorphine samples (158 subjects); 1182 apomorphine sulfate samples (101 subjects)",
    notes = "Pooled 5 phase I, 3 phase II and 1 phase III studies (Supplementary Table S1: CTH-103, -104, -105, -106, -107, -200, -201, -203, -301). Demographics from Table 2; female 44/158 and race counts White 67, Black 3, Asian 85, Other 3 summed across the two groups of Table 2."
  )

  ini({
    # Apomorphine disposition (Table 3; typical values at WT = 69.3 kg, study != CTH-103)
    lcl <- log(80.7); label("Apparent apomorphine clearance CL/F (L/h); k23 = CL/V2") # Table 3 'CL/F' = 80.7 L/h (RSE 25.0%)
    lvc <- log(438); label("Apparent apomorphine central volume V2/F at 69.3 kg (L)") # Table 3 'V2/F' = 438 L (RSE 8.39%)
    lk12 <- log(0.613); label("Distribution rate constant central -> peripheral k26 (1/h)") # Table 3 'k26' = 0.613 1/h (RSE 14.0%)
    lk21 <- log(0.00480); label("Distribution rate constant peripheral -> central k62 (1/h)") # Table 3 'k62' = 0.00480 1/h (RSE 15.9%)

    # Absorption (Table 3; Figure 1)
    lka <- log(6.58); label("Sublingual absorption / transit rate constant k17 = k78 = k82 at 2 min contact time (1/h)") # Table 3 'ka for sublingual administration' = 6.58 1/h (RSE 16.2%)
    lka2 <- log(17.6); label("Subcutaneous absorption / transit rate constant k49 = k9T10 = k10T2 (1/h)") # Table 3 'ka for subcutaneous administration' = 17.6 1/h (RSE 8.40%)
    lfdepot <- log(0.202); label("Sublingual bioavailability relative to subcutaneous, Biorsc, at 20 mg (fraction)") # Table 3 'Fraction absorbed relative to subcutaneous administration' = 0.202 (RSE 12.8%)
    lfdepot_sulf <- log(1 - 0.910); label("Swallowed fraction of the sublingual dose absorbed as apomorphine sulfate, F5 = 1 - Biosl (fraction)") # Table 3 'Fraction not swallowed and available for sublingual absorption' Biosl = 0.910 (RSE 5.86%); Figure 1: F5 = 1 - Biosl, F1 = Biorsc x Biosl

    # Apomorphine sulfate (Table 3)
    lvc_sulf <- log(1.42); label("Apparent apomorphine sulfate volume V3/F, male (L)") # Table 3 'V3/F' = 1.42 L (RSE 34.2%)
    lkel_sulf <- log(1.28); label("Apomorphine sulfate elimination rate constant k30 (1/h)") # Table 3 'k30' = 1.28 1/h (RSE 38.5%)
    lka_sulf <- log(0.205); label("Absorption rate constant of apomorphine sulfate from the gastrointestinal tract k53 (1/h)") # Table 3 'ka for apomorphine sulfate absorption from the gastrointestinal tract' = 0.205 1/h (RSE 87.9%)

    # Covariate effects (Table 3; Results 'Final model'; Methods Eqs. 4 and 5)
    e_dose_fdepot <- -0.206; label("Power exponent of sublingual dose on F1 (unitless; reference 20 mg)") # Table 3 'Dose of sublingual apomorphine on F1' = -0.206 (RSE 165%)
    e_slct_ka <- -0.194; label("Power exponent of contact time under the tongue on sublingual ka (unitless; reference 2 min)") # Table 3 'Contact time under the tongue for sublingual film on ka for sublingual administration' = -0.194 (RSE 61.5%)
    e_wt_vc <- 1.53; label("Power exponent of body weight on V2/F (unitless; reference 69.3 kg)") # Table 3 'Body weight on V2/F' = 1.53 (RSE 24.8%)
    e_study_cth103_vc <- -0.555; label("Proportional shift in V2/F for study CTH-103 (unitless)") # Table 3 'Study CTH-103 on V2/F' = -0.555 (RSE 8.59%)
    e_sexf_vc_sulf <- -0.310; label("Proportional shift in V3/F for female sex (unitless)") # Table 3 'Female sex on V3/F' = -0.310 (RSE 21.4%)

    # IIV: Methods Eq. 2 defines %CV = sqrt(omega^2) * 100, so omega^2 = (CV/100)^2
    etalcl ~ 0.128881 # Table 3 IIV 'k23' = 35.9 %CV; 0.359^2. k23 = CL/V2 * exp(eta) is CL_i = CL * exp(eta)
    etalvc ~ 0.173889 # Table 3 IIV 'V2/F' = 41.7 %CV; 0.417^2
    etalvc_sulf ~ 0.162409 # Table 3 IIV 'V3/F' = 40.3 %CV; 0.403^2
    etalka ~ 0.127449 # Table 3 IIV 'ka for sublingual administration' = 35.7 %CV; 0.357^2
    etalfdepot ~ 0.146689 # Table 3 IIV 'F1' = 38.3 %CV; 0.383^2

    # Residual error (Table 3 'Residual variability'; Results 'Base model': proportional for both analytes)
    propSd <- 0.443; label("Proportional residual error, apomorphine (fraction)") # Table 3 'Apomorphine' = 44.3% (RSE 3.26%)
    propSd_sulf <- 0.580; label("Proportional residual error, apomorphine sulfate (fraction)") # Table 3 'Apomorphine sulfate' = 58.0% (RSE 3.81%)
  })
  model({
    # Individual parameters
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) *
      (WT / 69.3)^e_wt_vc *
      (1 + e_study_cth103_vc * STUDY_CTH103)
    ka <- exp(lka + etalka) * (DUR_SL_CONTACT / 2)^e_slct_ka
    ka2 <- exp(lka2)
    vc_sulf <- exp(lvc_sulf + etalvc_sulf) * (1 + e_sexf_vc_sulf * SEXF)
    kel_sulf <- exp(lkel_sulf)
    ka_sulf <- exp(lka_sulf)
    fdepot_sulf <- exp(lfdepot_sulf)

    # Micro-constants: k23 (metabolism to the sulfate, the only apomorphine
    # elimination pathway), k26, k62
    kel <- cl / vc
    k12 <- exp(lk12)
    k21 <- exp(lk21)

    # Sublingual chain depot (CMT 1) -> transit1 (CMT 7) -> transit2 (CMT 8) -> central
    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ka * transit1
    d/dt(transit2) <- ka * transit1 - ka * transit2
    # Subcutaneous chain depot2 (CMT 4) -> transit3 (CMT 9) -> transit4 (CMT 10) -> central
    d/dt(depot2) <- -ka2 * depot2
    d/dt(transit3) <- ka2 * depot2 - ka2 * transit3
    d/dt(transit4) <- ka2 * transit3 - ka2 * transit4
    # Apomorphine central (CMT 2) and peripheral (CMT 6)
    d/dt(central) <- ka * transit2 +
      ka2 * transit4 -
      kel * central -
      k12 * central +
      k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    # Swallowed fraction of the sublingual dose (CMT 5) -> apomorphine sulfate (CMT 3)
    d/dt(depot_sulf) <- -ka_sulf * depot_sulf
    d/dt(central_sulf) <- kel * central +
      ka_sulf * depot_sulf -
      kel_sulf * central_sulf

    # Bioavailability (Figure 1): F1 = Biorsc x Biosl x dose effect; F4 = 1; F5 = 1 - Biosl
    f(depot) <- (1 - fdepot_sulf) *
      exp(lfdepot + etalfdepot) *
      (DOSE_APOMORPHINE_SL_MG / 20)^e_dose_fdepot
    f(depot_sulf) <- fdepot_sulf

    # mg / L -> ng/mL
    Cc <- 1000 * central / vc
    Cc_sulf <- 1000 * central_sulf / vc_sulf

    Cc ~ prop(propSd)
    Cc_sulf ~ prop(propSd_sulf)
  })
}
