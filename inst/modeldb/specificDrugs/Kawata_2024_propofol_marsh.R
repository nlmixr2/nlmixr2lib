Kawata_2024_propofol_marsh <- function() {
  description <- "Extended Marsh model for propofol plasma and brain concentrations in Japanese adults and adolescents undergoing awake craniotomy (Kawata 2024; N = 29 model-building patients). Three-compartment intravenous PK with the Marsh micro-constants and weight-proportional central volume fixed, plus a brain compartment whose input carries a brain/plasma partition coefficient Kp: dCbrain/dt = k1e * (Kp * Cc - Cbrain). k1e and Kp were estimated; IIV on k10, k13 and Kp; proportional error on plasma and additive error on brain concentrations; no covariates retained beyond the Marsh weight scaling."
  reference <- paste(
    "Kawata M, Yonezawa A, Mineharu Y, Itohara K, Mizota T, Matsui Y, Kikuchi T,",
    "Yamao Y, Hattori EY, Hamada M, Hira D, Furukawa K, Miyamoto S, Terada T,",
    "Matsubara K, Arakawa Y. Development of extended pharmacokinetic models for propofol",
    "based on measured blood and brain concentrations. Sci Rep. 2024;14:6326.",
    "doi:10.1038/s41598-024-56863-z. PMCID: PMC10943190.",
    "Plasma PK fixed to the Marsh model (Marsh B, White M, Morton N, Kenny GN. Br J Anaesth. 1991;67:41-48. doi:10.1093/bja/67.1.41)."
  )
  vignette <- "Kawata_2024_propofol"
  units <- list(time = "min", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline. Enters only through the Marsh central volume, V1 = 0.228 L/kg * WT (linear, no reference weight). Supplementary code 1 'V1=THETA(6)*BW'.",
      source_name = "BW"
    )
  )

  covariatesDataExcluded <- list(
    HT = list(
      description = "Height",
      units = "m",
      type = "continuous",
      notes = "Paper reports height in m (Supplementary Table S1). Screened as a power effect on Q4 (k1e * V1); not retained (Supplementary Table S2, forward-inclusion delta OBJ -0.40)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened on Q4; not retained (Supplementary Table S2, delta OBJ -0.34)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on Q4 as (1 + theta * sex) with sex = 1 for females; not retained (Supplementary Table S2, delta OBJ 0.00)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on Q4; not retained (Supplementary Table S2, delta OBJ -1.75)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on Q4; not retained (Supplementary Table S2, delta OBJ -1.59)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Paper reports albumin in g/dL (Supplementary Table S1). Screened on Q4; not retained (Supplementary Table S2, delta OBJ -2.64)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened on Q4; significant in forward inclusion (delta OBJ -5.31) but removed in backward elimination (+5.31 < 7.88) (Supplementary Table S2)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened on Q4; not retained (Supplementary Table S2, delta OBJ -0.40). The Methods also name Cockcroft-Gault creatinine clearance."
    )
  )

  compartmentData <- list(
    central = list(analyte = "propofol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "propofol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "propofol", units = "mg", specimen = "plasma", verified = TRUE),
    brain = list(analyte = "propofol", units = "ug/g", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 29,
    n_studies = 1,
    age_range = "13-76 years",
    age_median = "49 years",
    weight_range = "44.0-92.8 kg",
    weight_median = "62.5 kg",
    sex_female_pct = 44.8,
    race_ethnicity = c(Asian = 100),
    disease_state = "Brain tumour (glioma, 26/29) or other brain lesion (epilepsy, cavernous hemangioma, tuberculoma) undergoing awake craniotomy under propofol anesthesia",
    dose_range = "Intravenous propofol infusion (target-controlled) until awakening; plasma (radial artery) and brain (lobectomy / uncapping margin) samples immediately before and 5, 15, 25 and 75 min after cessation",
    regions = "Japan (Kyoto University Hospital)",
    notes = "Model-building group = first 29 of 57 consecutively enrolled patients (May 2019 - June 2022); the remaining 28 formed an external validation group. Demographics from Supplementary Table S1 (median (range)); 13 of 29 female. 105 plasma and 80 brain concentrations. Japanese ethnicity inferred from the single-centre Kyoto enrolment."
  )

  ini({
    # Marsh three-compartment PK, fixed (Table 1 'Fixed'; Supplementary code 1 $THETA 1-6 FIXED)
    lkel <- fixed(log(0.119)); label("Elimination rate constant k10 (1/min)") # Table 1: k10 = 0.119 Fixed
    lk12 <- fixed(log(0.114)); label("Rate constant central -> peripheral1 k12 (1/min)") # Table 1: k12 = 0.114 Fixed
    lk13 <- fixed(log(0.0419)); label("Rate constant central -> peripheral2 k13 (1/min)") # Table 1: k13 = 0.0419 Fixed
    lk21 <- fixed(log(0.055)); label("Rate constant peripheral1 -> central k21 (1/min)") # Table 1: k21 = 0.055 Fixed
    lk31 <- fixed(log(0.0033)); label("Rate constant peripheral2 -> central k31 (1/min)") # Table 1: k31 = 0.0033 Fixed
    lvc <- fixed(log(0.228)); label("Central volume per kg body weight V1 (L/kg)") # Table 1: V1 = 0.228 L/kg Fixed

    # Brain compartment (estimated)
    lke0 <- log(0.014); label("Plasma-brain distribution rate constant k1e = ke0' (1/min)") # Table 1: k1e = 0.014 (34.8 %RSE)
    lkp_brain <- log(1.30); label("Brain/plasma partition coefficient Kp (unitless)") # Table 1: Kp = 1.30 (8.92 %RSE)

    # IIV: exponential; Table 1 reports CV% = 100 * sqrt(omega^2) (the same scale
    # Supplementary code 2 shows for the new model), so omega^2 = (CV/100)^2
    etalkel ~ 0.035344 # Table 1: IIV for k10 18.8 CV% -> 0.188^2
    etalk13 ~ 0.356409 # Table 1: IIV for k13 59.7 CV% -> 0.597^2
    etalkp_brain ~ 0.165649 # Table 1: IIV for Kp 40.7 CV% -> 0.407^2

    # Residual error (Supplementary code 1 $ERROR: CMT 1 Y=F*(1+EPS(1)), else Y=F+EPS(2))
    propSd <- 0.196; label("Proportional residual error on plasma Cc (fraction)") # Table 1: proportional error for plasma 19.6 CV%
    addSd_Cbrain <- 0.869; label("Additive residual error on brain Cbrain (ug/g)") # Table 1: additive error for brain 0.869 ug/g
  })

  model({
    kel <- exp(lkel + etalkel)
    k12 <- exp(lk12)
    k13 <- exp(lk13 + etalk13)
    k21 <- exp(lk21)
    k31 <- exp(lk31)
    vc <- exp(lvc) * WT
    ke0 <- exp(lke0)
    kp_brain <- exp(lkp_brain + etalkp_brain)

    d/dt(central) <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    Cc <- central / vc

    # Brain state is a concentration (NONMEM A(4) observed with no scale factor);
    # brain-tissue density taken as 1 g/mL so ug/mL plasma maps to ug/g brain.
    # Methods: dCbrain/dt = k1e * Kp * Ccentral - ke0' * Cbrain, with ke0' = k1e.
    d/dt(brain) <- ke0 * kp_brain * Cc - ke0 * brain
    Cbrain <- brain

    Cc ~ prop(propSd)
    Cbrain ~ add(addSd_Cbrain)
  })
}
