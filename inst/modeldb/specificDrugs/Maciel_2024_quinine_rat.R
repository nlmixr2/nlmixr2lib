Maciel_2024_quinine_rat <- function() {
  description <- paste(
    "Preclinical (rat, male Wistar infected with Plasmodium berghei). Two-compartment",
    "population PK model for quinine in WHOLE BLOOD after a single 20 mg/kg IV bolus of",
    "free quinine or of quinine-loaded polymeric nanocapsules coated with polysorbate 80",
    "(NCP80, anionic) or prepared with Eudragit RS100 (NCEUD, cationic). Nanoencapsulation",
    "(either coating) lowers clearance and central volume; the peripheral volume is",
    "formulation-specific. Parameters are per kg body weight. See",
    "modellib('Maciel_2024_quinine_mouse') for the rat-to-mouse allometric translation",
    "and the parasitemia PK/PD model built on it."
  )
  reference <- paste(
    "Maciel TR, Funguetto-Ribeiro AC, Olivo LB, Teixeira FEG, Pacheco CO, de Araujo BV,",
    "Haas SE. (2024). Improved Malaria Therapy with Cationic Nanocapsules Demonstrated in",
    "Plasmodium berghei-Infected Rodents Using Whole Blood Surrogate Population PK/PD",
    "Modeling. Pharmaceutics 16(11):1369. doi:10.3390/pharmaceutics16111369."
  )
  vignette <- "Maciel_2024_quinine"

  # Table 1 labels CL and Q "mL/h" and the volumes "mL", but the values are per
  # kg: the dose entered as 20 (mg/kg) and the reported mean whole-blood
  # AUC0-inf of 127.8 ug*h/mL for free quinine (Results 3.1) gives
  # 20 / 127.8 = 0.156 L/h/kg, the Table 1 CL_pop exactly.
  units <- list(time = "h", dosing = "mg/kg", concentration = "mg/L")

  covariateData <- list(
    FORM_QUININE_NCP80 = list(
      description = "Quinine dose given as the polysorbate-80-coated poly(epsilon-caprolactone) nanocapsule (NCP80)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (free quinine solution, or NCEUD when FORM_QUININE_NCEUD = 1)",
      notes = paste(
        "Per-animal treatment-arm indicator. FORM_QUININE_NCP80 and FORM_QUININE_NCEUD",
        "are mutually exclusive; both 0 = free quinine. Their sum is the 'NC' indicator",
        "that Maciel 2024 uses on CL and V1 (either nanocapsule), while V2 takes a",
        "separate value for each coating (Table 1 and Table 2 footnote)."
      ),
      source_name = "NCP80"
    ),
    FORM_QUININE_NCEUD = list(
      description = "Quinine dose given as the Eudragit RS100 (cationic) nanocapsule (NCEUD)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (free quinine solution, or NCP80 when FORM_QUININE_NCP80 = 1)",
      notes = paste(
        "Per-animal treatment-arm indicator, mutually exclusive with FORM_QUININE_NCP80;",
        "both 0 = free quinine (Table 1)."
      ),
      source_name = "NCEUD"
    )
  )

  compartmentData <- list(
    central = list(analyte = "quinine", units = "mg/kg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "quinine", units = "mg/kg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "rat (male Wistar, Plasmodium berghei-infected)",
    n_subjects = 18L,
    n_studies = 1L,
    weight_range = "(not reported)",
    sex_female_pct = 0,
    disease_state = paste(
      "Plasmodium berghei malaria; PK studied on day 7 after infection. The erythrocyte",
      "partition coefficients were measured in blood from infected rats with 9.25 +/- 1.5%",
      "parasitemia (Methods 2.1)."
    ),
    dose_range = "Single 20 mg/kg quinine IV bolus as free quinine, NCP80 or NCEUD (n = 6 per group)",
    regions = "Brazil (Federal University of Pampa)",
    notes = paste(
      "197 whole-blood concentrations sampled at 0.083-24 h (Results 3.2). The whole-blood",
      "values were not measured directly: the plasma concentrations of Michels 2019",
      "(Int J Nanomedicine 14:10165, doi:10.2147/IJN.S227914) were multiplied by the",
      "in-vitro erythrocyte partition coefficients Kp = 3.74 (free), 8.92 (NCP80) and",
      "6.48 (NCEUD) (Methods 2.2). NONMEM 7.4, FOCE-I."
    )
  )

  ini({
    # Table 1 reports each formulation's typical value as its own estimate
    # (with its own RSE), so each stratum carries its own parameter. Table 2
    # ('Influence of Nanostructures in Rats') prints the same model as additive
    # shifts: CL_NC = CL_FREE - 0.0973, V1_NC = V1_FREE - 0.233,
    # V2_NCP80 = V2_FREE + 0.138, V2_NCEUD = V2_FREE + 0.0335, which agree
    # with these values to the printed precision.
    lcl_free <- log(0.156); label("Clearance of free quinine (L/h/kg)") # Table 1: CL_pop = 0.156 (RSE 9.0%)
    lcl_nc <- log(0.059); label("Clearance of nanoencapsulated quinine, NCP80 or NCEUD (L/h/kg)") # Table 1: CL_NC = 0.059 (RSE 15.5%)
    lvc_free <- log(0.261); label("Central volume of free quinine (L/kg)") # Table 1: V1_pop = 0.261 (RSE 10.2%)
    lvc_nc <- log(0.028); label("Central volume of nanoencapsulated quinine, NCP80 or NCEUD (L/kg)") # Table 1: V1_NC = 0.028 (RSE 11.5%)
    lq <- log(0.025); label("Intercompartmental clearance, all formulations (L/h/kg)") # Table 1: Q_pop = 0.025 (RSE 25.3%)
    lvp_free <- log(0.066); label("Peripheral volume of free quinine (L/kg)") # Table 1: V2_pop = 0.066 (RSE 49.9%)
    lvp_ncp80 <- log(0.204); label("Peripheral volume of NCP80 quinine (L/kg)") # Table 1: V2_NCP80 = 0.204 (RSE 17.2%)
    lvp_nceud <- log(0.099); label("Peripheral volume of NCEUD quinine (L/kg)") # Table 1: V2_NCEUD = 0.099 (RSE 39.7%)

    # Table 1 'Interindividual Variability (CV%)' column; the bracketed number
    # after each CV is read as its RSE. OMEGA = log(CV^2 + 1).
    etalcl ~ 0.00719902 # Table 1: IIV CL = 8.5% CV
    etalvc ~ 0.000168986 # Table 1: IIV V1 = 1.3% CV
    etalq ~ 0.591338 # Table 1: IIV Q = 89.8% CV
    etalvp ~ 1.48168 # Table 1: IIV V2 = 184.4% CV

    propSd <- 0.436; label("Proportional residual error (fraction)") # Table 1: Proportional model = 43.6 (RSE 7.3%)
  })

  model({
    # Either nanocapsule (the 'NC' covariate of Table 1)
    nc <- FORM_QUININE_NCP80 + FORM_QUININE_NCEUD

    cl <- exp(lcl_free * (1 - nc) + lcl_nc * nc + etalcl)
    vc <- exp(lvc_free * (1 - nc) + lvc_nc * nc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(
      lvp_free * (1 - nc) +
        lvp_ncp80 * FORM_QUININE_NCP80 +
        lvp_nceud * FORM_QUININE_NCEUD +
        etalvp
    )

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Equations 2-3; Equation 4 (the eliminated amount) is not carried.
    d/dt(central) <- k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Whole-blood quinine concentration (mg/L)
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
