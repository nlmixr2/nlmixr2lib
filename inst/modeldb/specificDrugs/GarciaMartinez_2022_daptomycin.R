GarciaMartinez_2022_daptomycin <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous daptomycin in",
    "hospitalised adults with normal renal function or renal impairment",
    "(46 patients, 157 steady-state serum concentrations, 4-12 mg/kg q24h",
    "as a 30-min infusion; two Spanish hospitals). The ODE system describes",
    "UNBOUND daptomycin with linear elimination and linear peripheral",
    "distribution; the observed total serum concentration is reconstructed",
    "algebraically from a single-site saturable protein-binding equation",
    "Ctotal = Cu + Bmax * Cu / (KD + Cu). Unbound clearance scales with",
    "Cockcroft-Gault creatinine clearance by a power function centred on",
    "92.8 mL/min/1.73 m^2. Log-normal IIV on unbound clearance and unbound",
    "peripheral volume; residual error is additive on the log scale."
  )
  reference <- paste(
    "Garcia-Martinez T, Belles-Medall MD, Garcia-Cremades M,",
    "Ferrando-Piqueres R, Mangas-Sanjuan V, Merino-Sanjuan M. Population",
    "pharmacokinetic/pharmacodynamic modelling of daptomycin for schedule",
    "optimization in patients with renal impairment. Pharmaceutics.",
    "2022;14(10):2226. doi:10.3390/pharmaceutics14102226."
  )
  vignette <- "GarciaMartinez_2022_daptomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The ODE states hold unbound daptomycin (Figure 1: 'Unbound Central
  # concentration' and 'Unbound Peripheral concentration'); total drug is an
  # algebraic observable. The assay matrix is serum (Methods: 'the serum was
  # stored at -20 C until the drug concentrations were measured').
  compartmentData <- list(
    central = list(analyte = "daptomycin unbound", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "daptomycin unbound", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated by the Cockcroft-Gault formula, reported in mL/min/1.73 m^2",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Garcia-Martinez 2022 Methods: 'Creatinine clearance (CLCR) was",
        "estimated using the Cockroft-Gault formula' with the renal-function",
        "categories stated in mL/min/1.73 m^2; Table 1 median 93 (IQR",
        "50-136) mL/min/1.73 m^2. Measured on the day of drug sampling, so",
        "a single value per sampling occasion. Enters as (CRCL / 92.8)^0.19",
        "on unbound clearance (Equation 4); 92.8 is the median the paper",
        "centres on."
      ),
      source_name = "CLCR"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened in the univariate covariate analysis but not retained (Results 3.2.2: 'Additional continuous and categorical covariates were tested but did not result in a statistical (dOFV < 3.84) reduction of the OFV').",
      source_name = "age"
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained (Results 3.2.2). The paper reports height in m (Table 1 median 1.7 m).",
      source_name = "height"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained (Results 3.2.2). Weight enters the paper's simulations only through mg/kg dosing.",
      source_name = "body weight"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained (Results 3.2.2). Table 1 prints the unit as 'g/dL', a misprint for mg/dL (median 0.9).",
      source_name = "serum creatinine"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained (Results 3.2.2); Discussion: 'albumin levels could not be related to the binding mechanism'. Table 1 median 2.9 g/dL.",
      source_name = "albumin"
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      reference_category = "male",
      notes = "Screened but not retained (Results 3.2.2); only 3 of 46 patients were female.",
      source_name = "gender"
    ),
    CONMED_STATIN = list(
      description = "Statin co-administration",
      units = "(binary)",
      type = "binary",
      reference_category = "no statin",
      notes = "Screened as a categorical covariate but not retained (Methods; Results 3.2.2).",
      source_name = "statin co-administration"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 46L,
    n_studies = 1L,
    n_observations = 157L,
    age_range = "IQR 59-81 years",
    age_median = "68 years",
    weight_range = "IQR 65-85 kg",
    weight_median = "75 kg",
    sex_female_pct = 7,
    race_ethnicity = "not reported (Spanish hospital cohort)",
    disease_state = "Hospitalised adults receiving daptomycin for Gram-positive infection (S. aureus 47%, S. epidermidis 33% of the 36 patients with an isolate); renal replacement therapy excluded",
    renal_function = "Cockcroft-Gault CLCR median 93 (IQR 50-136) mL/min/1.73 m^2; >90: 35%, 60-89: 28%, 30-59: 30%, 15-29: 7%",
    dose_range = "4-12 mg/kg q24h as a 30-min IV infusion (median 675 mg, 9.1 mg/kg)",
    regions = "Spain (General University Hospital of Castellon; Arnau de Vilanova Hospital, Valencia)",
    notes = "Prospective multi-centre TDM study; five serial samples from day 4 (steady state): pre-dose, 0.5, 1-2 and 4-10 h after the end of infusion and before the next dose. Total daptomycin by HPLC-UV; unbound daptomycin was not measured. Demographics from Table 1."
  )

  ini({
    # Structural parameters refer to UNBOUND daptomycin (Results 3.2.2:
    # 'PK parameters are referred to unbound daptomycin concentrations').
    lcl <- log(6.98); label("Unbound clearance at CRCL = 92.8 mL/min/1.73 m^2 (L/h)") # Table 2 'CL (L/h)' = 6.98
    lvc <- log(0.95); label("Unbound central volume of distribution (L)") # Table 2 'V1 (L)' = 0.95
    lq <- log(1.96); label("Unbound intercompartmental clearance (L/h)") # Table 2 'Q (L/h)' = 1.96
    lvp <- log(21); label("Unbound peripheral volume of distribution (L)") # Table 2 'V2 (L/h)' = 21; unit misprinted as L/h in the table, L in Figure 1 caption

    # Single-site saturable binding, Equation 2: Ctotal = Cfree + Bmax * Cfree / (KD + Cfree).
    lbmax <- log(160); label("Maximal protein-binding capacity Bmax (mg/L)") # Table 2 'Bmax (mg/L)' = 160
    lkd <- log(3.56); label("Binding equilibrium dissociation constant KD (mg/L)") # Table 2 'KD (mg/L)' = 3.56 (Discussion's 'KD = 1.96' repeats the Q value)

    e_crcl_cl <- 0.19; label("Power exponent of CRCL / 92.8 on unbound clearance (unitless)") # Table 2 'CrCl on CL' = 0.19; Equation 4

    # IIV reported as % with exponential random effects; read as omega (SD
    # of eta) x 100, so variance = (pct / 100)^2.
    etalcl ~ 0.1024 # Table 2 IIV 'CL (%)' = 32; 0.32^2
    etalvp ~ 0.2209 # Table 2 IIV 'V2 (%)' = 47; 0.47^2

    expSd <- 0.22; label("Additive residual error on the log scale (SD)") # Table 2 'Additive on Log-scale (%)' = 22
  })

  model({
    cl <- exp(lcl + etalcl) * (CRCL / 92.8)^e_crcl_cl
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)
    bmax <- exp(lbmax)
    kd <- exp(lkd)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Unbound concentration from the ODE state; total (observed) concentration
    # from the saturable binding equation (Equation 2).
    Cunbound <- central / vc
    Cc <- Cunbound + bmax * Cunbound / (kd + Cunbound)

    Cc ~ lnorm(expSd)
  })
}
