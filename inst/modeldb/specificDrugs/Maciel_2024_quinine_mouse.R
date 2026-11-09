Maciel_2024_quinine_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, Plasmodium berghei-infected). Population PK/PD model for",
    "quinine given as free drug or as polysorbate-80-coated (NCP80) or Eudragit RS100",
    "(NCEUD) nanocapsules in the Peters four-day suppressive test. The whole-blood PK",
    "is a two-compartment model translated from rats by allometry (Table 2 power",
    "equations in body weight, with additive nanocapsule shifts on CL, V1 and V2) and",
    "carries no between-animal variability; the doses enter the central compartment",
    "with complete bioavailability. Percent parasitemia grows logistically towards a",
    "fixed ceiling and is killed by an Emax function of the whole-blood concentration,",
    "with a formulation-specific Emax. See modellib('Maciel_2024_quinine_rat') for the",
    "rat whole-blood popPK model the translation starts from."
  )
  reference <- paste(
    "Maciel TR, Funguetto-Ribeiro AC, Olivo LB, Teixeira FEG, Pacheco CO, de Araujo BV,",
    "Haas SE. (2024). Improved Malaria Therapy with Cationic Nanocapsules Demonstrated in",
    "Plasmodium berghei-Infected Rodents Using Whole Blood Surrogate Population PK/PD",
    "Modeling. Pharmaceutics 16(11):1369. doi:10.3390/pharmaceutics16111369."
  )
  vignette <- "Maciel_2024_quinine"

  # Volumes are L/kg and clearances L/h/kg, as in the rat model the allometry
  # was derived from; with a 20 mg/kg dose the free-quinine peak is about
  # 4.4 mg/L at WT = 0.025 kg, the Figure 5 level.
  units <- list(time = "h", dosing = "mg/kg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the Table 2 allometric equations P = a * WT^b directly, with no",
        "reference weight (a is the value at WT = 1 kg). The exponents are negative",
        "and steep, and the nanocapsule shifts are constants, so the model is only",
        "meaningful over the narrow mouse weight range it was built for: the",
        "nanocapsule central volume V1 = 0.00227 * WT^-2.06 - 3.695 falls to zero at",
        "WT = 0.0276 kg and is negative above it. The mouse weight is not reported;",
        "the vignette uses 0.025 kg."
      ),
      source_name = "WT"
    ),
    FORM_QUININE_NCP80 = list(
      description = "Quinine dose given as the polysorbate-80-coated poly(epsilon-caprolactone) nanocapsule (NCP80)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (free quinine solution, or NCEUD when FORM_QUININE_NCEUD = 1)",
      notes = paste(
        "Per-animal treatment-arm indicator, mutually exclusive with",
        "FORM_QUININE_NCEUD; both 0 = free quinine (also used for the saline controls,",
        "which receive no dose). The sum of the two indicators is the 'NC' covariate",
        "on CL and V1 (Table 2); V2 and Emax take a separate value per coating."
      ),
      source_name = "NCP80"
    ),
    FORM_QUININE_NCEUD = list(
      description = "Quinine dose given as the Eudragit RS100 (cationic) nanocapsule (NCEUD)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (free quinine solution, or NCP80 when FORM_QUININE_NCP80 = 1)",
      notes = "Per-animal treatment-arm indicator, mutually exclusive with FORM_QUININE_NCP80 (Tables 2 and 3).",
      source_name = "NCEUD"
    )
  )

  compartmentData <- list(
    central = list(analyte = "quinine", units = "mg/kg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "quinine", units = "mg/kg", specimen = "tissue", verified = TRUE),
    parasites = list(
      analyte = "Plasmodium berghei-infected erythrocytes",
      units = "% parasitemia",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (Plasmodium berghei-infected; strain not reported)",
    n_subjects = 20L,
    n_studies = 1L,
    weight_range = "(not reported)",
    disease_state = paste(
      "Plasmodium berghei malaria, Peters four-day suppressive test: infection on day 0",
      "(inoculum 1e8 infected erythrocytes), treatment from day 0.5 to day 3.5, parasitemia",
      "read from blood smears every other day until death or day 15 (Methods 2.1)."
    ),
    dose_range = "Quinine 20 mg/kg twice daily intraperitoneally from day 0.5 to day 3.5, as free quinine, NCP80 or NCEUD; saline controls",
    regions = "Brazil (Federal University of Pampa)",
    notes = paste(
      "Four groups of five mice (free quinine, NCP80, NCEUD, saline); 147 parasitemia",
      "observations (Results 3.4). The efficacy data are from Michels 2019 (Int J",
      "Nanomedicine 14:10165). No quinine was measured in the mice: the PD was fitted",
      "against the whole-blood profiles simulated with the allometric PK (Methods 2.4,",
      "Figure 5). NONMEM 7.4."
    )
  )

  ini({
    # Whole-blood PK, Table 2 'Allometric Scaling from Rats to Mice for Free QN':
    # P = a * WT^b with WT in kg. Each lX is log(a), the value at WT = 1 kg.
    lcl <- fixed(log(0.0112)); label("Allometric coefficient of free-quinine clearance (L/h/kg at WT = 1 kg)") # Table 2: CL_SPECIE = 0.0112 x WT^-1.43
    lvc <- fixed(log(0.00227)); label("Allometric coefficient of free-quinine central volume (L/kg at WT = 1 kg)") # Table 2: V1_SPECIE = 0.00227 x WT^-2.06
    lvp <- fixed(log(0.031)); label("Allometric coefficient of free-quinine peripheral volume (L/kg at WT = 1 kg)") # Table 2: V2_SPECIE = 0.031 x WT^-1.22
    lq <- fixed(log(0.0136)); label("Allometric coefficient of intercompartmental clearance (L/h/kg at WT = 1 kg)") # Table 2: Q_SPECIE = 0.0136 x WT^-1.93
    e_wt_cl <- fixed(-1.43); label("Allometric exponent of body weight on clearance (unitless)") # Table 2: CL_SPECIE exponent -1.43
    e_wt_vc <- fixed(-2.06); label("Allometric exponent of body weight on central volume (unitless)") # Table 2: V1_SPECIE exponent -2.06
    e_wt_vp <- fixed(-1.22); label("Allometric exponent of body weight on peripheral volume (unitless)") # Table 2: V2_SPECIE exponent -1.22
    e_wt_q <- fixed(-1.93); label("Allometric exponent of body weight on intercompartmental clearance (unitless)") # Table 2: Q_SPECIE exponent -1.93

    # Table 2 'Influence of Nanostructures in Mice': additive shifts on the
    # allometric free-quinine values.
    e_form_nc_cl <- fixed(-1.353); label("Additive shift in clearance for either nanocapsule (L/h/kg)") # Table 2 mice: CL_NC = CL_FREE - 1.353
    e_form_nc_vc <- fixed(-3.695); label("Additive shift in central volume for either nanocapsule (L/kg)") # Table 2 mice: V1_NC = V1_FREE - 3.695
    e_form_ncp80_vp <- fixed(4.494); label("Additive shift in peripheral volume for NCP80 (L/kg)") # Table 2 mice: V2_NCP80 = V2_FREE + 4.494
    e_form_nceud_vp <- fixed(0.274); label("Additive shift in peripheral volume for NCEUD (L/kg)") # Table 2 mice: V2_NCEUD = V2_FREE + 0.274

    # Parasite growth and drug effect, Equation 5 and Table 3.
    lkg <- log(0.012); label("Parasite growth rate constant (1/h)") # Table 3: kg_pop = 0.012 h-1 (RSE 5.8%)
    lnmax <- fixed(log(9.2)); label("Maximum parasitemia (percent)") # Table 3: Nmax (%) (Fix) = 9.2
    # Table 3 reports Emax once per formulation; Results 3.4 reads them as the
    # formulation values ('from 0.055 to 0.025 and 0.069'), not as increments.
    lemax_free <- log(0.055); label("Maximum parasite killing rate, free quinine (1/h)") # Table 3: Emax_pop = 0.055 (RSE 71.3%)
    lemax_ncp80 <- log(0.025); label("Maximum parasite killing rate, NCP80 (1/h)") # Table 3: Emax_NCP80 = 0.025 (RSE 77.9%)
    lemax_nceud <- log(0.069); label("Maximum parasite killing rate, NCEUD (1/h)") # Table 3: Emax_NCEUD = 0.069 (RSE 61.3%)
    lec50 <- log(3.59); label("Whole-blood quinine concentration at half-maximal killing (mg/L)") # Table 3: EC50_pop = 3.59 mg/L (RSE 104.2%)
    # Not reported. Back-projected to t = 0 through the Table 3 logistic growth
    # from the saline-control median of 2.5% at 24 h digitised from Figure 6a;
    # it reproduces the control medians of about 3.6, 5.0, 6.2 and 7.3% at
    # 72-216 h.
    lrbase <- fixed(log(2.0)); label("Parasitemia at infection, t = 0 (percent)") # figure-derived: Figure 6a Control panel, see comment above

    # Table 3 'Interindividual Variability (CV%)'; OMEGA = log(CV^2 + 1).
    etalemax ~ 6.4e-09 # Table 3: IIV Emax printed as 0.008 in the CV% column, taken as 0.008% CV
    etalec50 ~ 0.00347496 # Table 3: IIV EC50 = 5.9% CV

    addSd <- 0.933; label("Additive residual error (percent parasitemia)") # Table 3: Additive Error = -0.933 (RSE 10.7%), printed with a negative sign and 'mg/L'
    propSd <- 0.0934; label("Proportional residual error (fraction)") # Table 3: Proporcional Error (%) = 9.34 (RSE 55.2%)
  })

  model({
    # Either nanocapsule (the 'NC' covariate of Table 2)
    nc <- FORM_QUININE_NCP80 + FORM_QUININE_NCEUD

    cl <- exp(lcl) * WT^e_wt_cl + e_form_nc_cl * nc
    vc <- exp(lvc) * WT^e_wt_vc + e_form_nc_vc * nc
    vp <- exp(lvp) *
      WT^e_wt_vp +
      e_form_ncp80_vp * FORM_QUININE_NCP80 +
      e_form_nceud_vp * FORM_QUININE_NCEUD
    q <- exp(lq) * WT^e_wt_q

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Whole-blood quinine concentration (mg/L), the PD driver
    Cc <- central / vc

    kg <- exp(lkg)
    nmax <- exp(lnmax)
    emax <- exp(
      lemax_free * (1 - nc) +
        lemax_ncp80 * FORM_QUININE_NCP80 +
        lemax_nceud * FORM_QUININE_NCEUD +
        etalemax
    )
    ec50 <- exp(lec50 + etalec50)
    rbase <- exp(lrbase)

    # Equation 5: logistic growth minus Emax killing, N in percent parasitemia
    kkill <- emax * Cc / (ec50 + Cc)
    d/dt(parasites) <- kg * (1 - parasites / nmax) * parasites - kkill * parasites
    parasites(0) <- rbase

    parasites ~ add(addSd) + prop(propSd)
  })
}
