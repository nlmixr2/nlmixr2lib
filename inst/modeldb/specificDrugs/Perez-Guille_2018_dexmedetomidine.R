"Perez-Guille_2018_dexmedetomidine" <- function() {
  description <- "Two-compartment IV population PK with sigmoidal Imax PD on heart rate (HR) and mean arterial pressure (MAP) fractional responses for dexmedetomidine in Mexican Mestizo children (2-18 y) undergoing ambulatory surgery, with a priori allometric scaling on CL and Q (exponent 0.75) and V1 and V2 (exponent 1) at a 70 kg reference weight (Perez-Guille et al. 2018, Tables 2 and 3, allometric model)"
  reference <- paste(
    "Perez-Guille MG, Toledo-Lopez A, Rivera-Espinosa L, Alemon-Medina R,",
    "Murata C, Lares-Asseff I, Chavez-Pacheco JL, Gomez-Garduno J,",
    "Zamora Gutierrez AL, Orozco-Galicia C, Ramirez-Morales K,",
    "Lugo-Goytia G.",
    "Population Pharmacokinetics and Pharmacodynamics of Dexmedetomidine",
    "in Children Undergoing Ambulatory Surgery.",
    "Anesth Analg. 2018;127(3):716-723.",
    "doi:10.1213/ANE.0000000000003413",
    sep = " "
  )
  vignette <- "Perez-Guille_2018_dexmedetomidine"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL") # Methods: dose 0.7 ug/kg single IV infusion over 10-15 min; plasma DEX measured by HPLC-ESI-MS/MS (5 pg/mL LLOQ); Table 2 CL in L/h, V1/V2 in L; Table 3 IC50 in ng/mL; with central in ug and vc in L, Cc = central/vc has units ug/L = ng/mL.

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "A priori allometric scaling per Perez-Guille 2018 Methods: CL and Q scale as (WT/70)^0.75 and V1 and V2 scale as (WT/70)^1; reference weight 70 kg. Cohort mean WT 43 kg (SD 19), age 2-18 years (Table 1). Age was tested as a separate covariate via a maturation model but was not retained (dBIC = 2.17; Results).",
      source_name        = "WT"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 30L,
    n_studies       = 1L,
    age_range       = "2-18 years",
    age_mean_sd     = "11 (SD 5) years",
    weight_range    = "Table 1 reports mean 43 (SD 19) kg; cohort weight range not separately tabulated",
    weight_mean_sd  = "43 (SD 19) kg",
    height_mean_sd  = "132 (SD 42) cm",
    sex_female_pct  = 30,
    race_ethnicity  = "Mexican Mestizo",
    disease_state   = "Healthy children with American Society of Anesthesiologists (ASA) physical status score of I (28/30) or II (2/30) scheduled for outpatient surgical procedures (urology, ORL, plastic, general). Exclusions: enzyme-inducing drugs, history of arrhythmias, delayed neurological development, malnutrition.",
    dose_range      = "Single IV infusion of dexmedetomidine 0.7 ug/kg over 10-15 min (exact infusion time recorded per patient and used as PK model input).",
    regions         = "Mexico (Instituto Nacional de Pediatria, Mexico City)",
    procedures      = "Urologic (6), ORL (12), plastic (8), general (4); mean duration of anesthesia 77 (SD 15) min.",
    co_medication   = "Standardised anaesthetic regimen: induction via inhaled sevoflurane + oxygen; fentanyl IV after IV access; IV propofol for endotracheal / laryngeal-mask insertion; maintenance with sevoflurane titrated to surgical stimulus.",
    sampling        = "Sparse: 2-5 venous samples per child (1 mL each) at randomly assigned times among 5, 10, 15, 20, 30, 45, 60, 90, 120, 180, 300, 420, and 600 min after end of infusion. HR and MAP recorded every 5 min during surgery and every 5 min for the first 30 min postoperatively, then hourly until discharge (Aldrete >8).",
    notes           = "Final fit by stochastic-EM + MCMC (Monolix v4.4). Allometric and proportional weight models both reported (Table 2); this entry encodes the allometric (footnote a) model. The PD model was fit sequentially with the PK parameters fixed at the individual EBE estimates (Methods)."
  )

  ini({
    # ============================================================
    # Population PK structural parameters (Table 2, allometric model)
    # ============================================================
    lcl <- log(20.8) ; label("Clearance at 70 kg (L/h)")                              # Table 2 (allometric model): theta_Cl = 20.8 L/h/70 kg (RSE 10%; bootstrap median 21.0, 95% CI 20.5-21.5)
    lvc <- log(21.9) ; label("Central volume V1 at 70 kg (L)")                        # Table 2: theta_V1 = 21.9 L/70 kg (RSE 17%; bootstrap median 21.8, 95% CI 19.9-24.1)
    lq  <- log(75.8) ; label("Intercompartmental clearance Q at 70 kg (L/h)")         # Table 2: theta_Q = 75.8 L/h/70 kg (RSE 12%; bootstrap median 74.1, 95% CI 65.7-80.4)
    lvp <- log(81.2) ; label("Peripheral volume V2 at 70 kg (L)")                     # Table 2: theta_V2 = 81.2 L/70 kg (RSE 9%; bootstrap median 81.7, 95% CI 79-83.5)

    # Allometric exponents -- fixed a priori per Methods (Anderson and Holford theoretical exponents)
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent on CL (unitless, FIXED)")     # Table 2: beta_Cl = 0.75 fixed
    e_wt_vc <- fixed(1)    ; label("Allometric exponent on V1 (unitless, FIXED)")     # Table 2: beta_V1 = 1.0 fixed
    e_wt_q  <- fixed(0.75) ; label("Allometric exponent on Q  (unitless, FIXED)")     # Table 2: beta_Q  = 0.75 fixed
    e_wt_vp <- fixed(1)    ; label("Allometric exponent on V2 (unitless, FIXED)")     # Table 2: beta_V2 = 1.0 fixed

    # ============================================================
    # Population PD: sigmoidal Imax model on HR fractional response
    # (Table 3, Heart rate panel; Methods Imax equation with So fixed to 1)
    # ============================================================
    limax_hr <- log(0.289) ; label("HR maximum fractional reduction Imax (unitless)") # Table 3 (HR): Imax = 28.9% (RSE 8%; bootstrap median 29.2, 95% CI 27-36.6); stored as fraction 0.289
    lic50_hr <- log(0.552) ; label("HR concentration giving 50% of Imax IC50 (ng/mL)")# Table 3 (HR): IC50 = 0.552 ng/mL (RSE 10%; bootstrap median 0.569, 95% CI 0.505-0.834)
    lhill_hr <- log(1.86)  ; label("HR Hill exponent gamma (unitless)")                # Table 3 (HR): gamma = 1.86 (RSE 17%; bootstrap median 1.65, 95% CI 1.19-2.0)

    # ============================================================
    # Population PD: sigmoidal Imax model on MAP fractional response
    # (Table 3, Mean arterial pressure panel)
    # ============================================================
    limax_map <- log(0.45)  ; label("MAP maximum fractional reduction Imax (unitless)") # Table 3 (MAP): Imax = 45% (RSE 17%; bootstrap median 48.3, 95% CI 16.1-52.6); stored as fraction 0.45
    lic50_map <- log(0.501) ; label("MAP concentration giving 50% of Imax IC50 (ng/mL)")# Table 3 (MAP): IC50 = 0.501 ng/mL (RSE 30%; bootstrap median 0.539, 95% CI 0.363-1.07)
    lhill_map <- log(1.26)  ; label("MAP Hill exponent gamma (unitless)")               # Table 3 (MAP): gamma = 1.26 (RSE 31%; bootstrap median 1.28, 95% CI 0.842-1.58)

    # ============================================================
    # Inter-individual variability
    # Table 2 / Table 3 column header reads "omega^2 (CV%)" but the
    # numerical values match the abstract's reported %BSV when read as
    # coefficients of variation (e.g., Table 2 lists 0.275 for CL, and
    # the abstract reports CL BSV = 27%). For a log-normal IIV the
    # internal variance is omega^2 = log(CV^2 + 1); the conversion is
    # documented next to each value below.
    # ============================================================
    etalcl ~ 0.072902     # Table 2 (allometric model): omega^2_Cl reported as CV 0.275; omega^2 = log(0.275^2 + 1) = 0.072902 (RSE 17%; bootstrap 0.275, 95% CI 0.259-0.289)
    etalvc ~ 0.039993     # Table 2: CV 0.202; omega^2 = log(0.202^2 + 1) = 0.039993 (RSE 52%; bootstrap 0.161, 95% CI 0.089-0.263)
    etalq  ~ 0.062044     # Table 2: CV 0.253; omega^2 = log(0.253^2 + 1) = 0.062044 (RSE 48%; bootstrap 0.248, 95% CI 0.186-0.289)
    etalvp ~ 0.046429     # Table 2: CV 0.218; omega^2 = log(0.218^2 + 1) = 0.046429 (RSE 35%; bootstrap 0.218, 95% CI 0.199-0.247)

    etalimax_hr ~ 0.004347 # Table 3 (HR): CV 0.066; omega^2 = log(0.066^2 + 1) = 0.004347 (RSE 73%; bootstrap 0.066, 95% CI 0.033-0.083)
    etalic50_hr ~ 0.009365 # Table 3 (HR): CV 0.097; omega^2 = log(0.097^2 + 1) = 0.009365 (RSE 87%; bootstrap 0.112, 95% CI 0.050-0.150)
    etalhill_hr ~ 0.134880 # Table 3 (HR): CV 0.38;  omega^2 = log(0.38^2  + 1) = 0.134880 (RSE 39%; bootstrap 0.350, 95% CI 0.290-0.395)

    etalimax_map ~ 0.001155 # Table 3 (MAP): CV 0.034; omega^2 = log(0.034^2 + 1) = 0.001155 (RSE 97%; bootstrap 0.040, 95% CI 0.029-0.05)
    etalic50_map ~ 0.014062 # Table 3 (MAP): CV 0.119; omega^2 = log(0.119^2 + 1) = 0.014062 (RSE 93%; bootstrap 0.115, 95% CI 0.072-0.151)
    etalhill_map ~ 0.056913 # Table 3 (MAP): CV 0.242; omega^2 = log(0.242^2 + 1) = 0.056913 (RSE 45%; bootstrap 0.210, 95% CI 0.116-0.266)

    # ============================================================
    # Residual error
    # PK proportional residual (Methods + Table 2 footnote: C_ij = P_ij * (1 + eps_p)).
    # PD residual reported as "Residual variability (%)" applied to the fractional
    # response; interpreted here as proportional residual (same column convention
    # as the PK panel). Documented in vignette Assumptions and deviations.
    # ============================================================
    propSd     <- 0.142 ; label("PK proportional residual SD on Cc (fraction)")        # Table 2 (allometric model): sigma_proportional = 14.2% (RSE 13%; bootstrap 14.5%, 95% CI 13.9-15.4)
    propSd_HR  <- 0.065 ; label("HR proportional residual SD on fractional response")  # Table 3 (HR): Residual variability = 6.5% (RSE 5%; bootstrap 6.5%, 95% CI 6.4-6.5)
    propSd_MAP <- 0.15  ; label("MAP proportional residual SD on fractional response") # Table 3 (MAP): Residual variability = 15% (RSE 4.6%; bootstrap 15.2%, 95% CI 14.9-15.3)
  })

  model({
    # 1. Individual PK parameters with a priori allometric weight scaling (reference 70 kg)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq  + etalq)  * (WT / 70)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp

    # 2. Individual PD parameters (sigmoidal Imax for HR and MAP)
    imax_hr <- exp(limax_hr + etalimax_hr)
    ic50_hr <- exp(lic50_hr + etalic50_hr)
    hill_hr <- exp(lhill_hr + etalhill_hr)

    imax_map <- exp(limax_map + etalimax_map)
    ic50_map <- exp(lic50_map + etalic50_map)
    hill_map <- exp(lhill_map + etalhill_map)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. Two-compartment IV disposition (dose into central as IV infusion)
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-   k12 * central        - k21 * peripheral1

    # 5. Plasma DEX concentration (central in ug, vc in L => Cc in ug/L = ng/mL)
    Cc <- central / vc

    # 6. PD response as fraction of baseline (So fixed to 1 per Methods):
    #    response_t / baseline = 1 - Imax * Cc^gamma / (IC50^gamma + Cc^gamma)
    HR  <- 1 - imax_hr  * Cc^hill_hr  / (ic50_hr^hill_hr  + Cc^hill_hr)
    MAP <- 1 - imax_map * Cc^hill_map / (ic50_map^hill_map + Cc^hill_map)

    # 7. Observation models
    Cc  ~ prop(propSd)
    HR  ~ prop(propSd_HR)
    MAP ~ prop(propSd_MAP)
  })
}
