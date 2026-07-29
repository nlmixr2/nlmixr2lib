Perrottet_2009_ganciclovir <- function() {
  description <- "Two-compartment population PK model for ganciclovir (administered as oral valganciclovir prodrug) in adult solid-organ transplant recipients (Perrottet 2009)"
  reference <- "Perrottet N, Csajka C, Pascual M, Manuel O, Lamoth F, Meylan P, Aubert JD, Venetz JP, Soccal P, Decosterd LA, Biollaz J, Buclin T. Population pharmacokinetics of ganciclovir in solid-organ transplant recipients receiving oral valganciclovir. Antimicrob Agents Chemother. 2009;53(7):3017-3023. doi:10.1128/AAC.00836-08"
  vignette <- "Perrottet_2009_ganciclovir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; recorded on each sampling occasion. Linear ratio (WT / 70) on V1. Cohort median 72 kg, range 46-115 kg (Perrottet 2009 Table 1).",
      source_name        = "BW"
    ),
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used to de-index MDRD-eGFR (mL/min/1.73 m^2) into GFR (L/h) inside model(). Reference BSA 1.73 m^2 (Perrottet 2009 simulation cohort).",
      source_name        = "BSA"
    ),
    CRCL = list(
      description        = "MDRD-estimated glomerular filtration rate (four-variable MDRD formula)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. The four-variable MDRD formula GFR = 175 * (Crs/88.4)^-1.154 * age^-0.203 * 0.742^SEXF returns mL/min/1.73 m^2; model() converts to L/h using individual BSA (Perrottet 2009 Methods, structural model and Discussion).",
      source_name        = "GFR_MDRD"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = male",
      notes              = "Multiplicative effects on both CL (theta_female = 1.21, females 21% higher) and V1 (theta_female = 0.78, females 22% lower). The CL sex effect is additional to the implicit sex correction inside the MDRD formula.",
      source_name        = "sex"
    ),
    TX_HEART = list(
      description        = "Heart-transplant recipient indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = non-heart solid-organ graft (kidney, lung, or liver)",
      notes              = "Time-fixed per subject. Perrottet 2009 cohort: 10/65 patients (15%). Applied as a multiplicative effect on the CL-vs-GFR slope (theta_heart = 0.86 vs theta_kidney = 1.68).",
      source_name        = "graft_type"
    ),
    TX_LUNG = list(
      description        = "Lung-transplant recipient indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = non-lung solid-organ graft (kidney, heart, or liver)",
      notes              = "Time-fixed per subject. Perrottet 2009 cohort: 12/65 patients (18%). Pooled with liver recipients under a shared CL slope (theta_lung/liver = 1.17 vs theta_kidney = 1.68).",
      source_name        = "graft_type"
    ),
    TX_LIVER = list(
      description        = "Liver-transplant recipient indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = non-liver solid-organ graft (kidney, heart, or lung)",
      notes              = "Time-fixed per subject. Perrottet 2009 cohort: 2/65 patients (3%). Pooled with lung recipients under a shared CL slope (theta_lung/liver = 1.17 vs theta_kidney = 1.68).",
      source_name        = "graft_type"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 65,
    n_studies      = 1,
    age_range      = "18-70 years (median 55)",
    weight_range   = "46-115 kg (median 72)",
    sex_female_pct = 31,
    disease_state  = "Adult solid-organ transplant recipients receiving oral valganciclovir for CMV prophylaxis or treatment, or intravenous ganciclovir for CMV disease.",
    graft_breakdown = "Kidney 41 (63%), heart 10 (15%), lung 12 (18%), liver 2 (3%).",
    serostatus     = "CMV D+/R- 22 (34%), D+/R+ 28 (43%), D-/R+ 15 (23%).",
    renal_function = "Wide range; cohort serum creatinine 29-691 umol/L (median 108); Cockroft-Gault GFR 10-170 mL/min reported in Discussion.",
    dose_range     = "Oral valganciclovir 450 mg or 900 mg QD (prophylaxis) up to 900 mg BID (treatment), each adjusted to renal function. Intravenous ganciclovir 5 mg/kg Q12H in two patients.",
    regions        = "Switzerland (single-centre study: University Hospitals of Lausanne and Geneva).",
    n_samples      = 437,
    notes          = "Perrottet 2009 Table 1. Comorbidities included cardiopathy (19/65), overweight (9/65), and cystic fibrosis (1/65). Concomitant medications evaluated as covariates (none retained): calcineurin inhibitors, mycophenolate, cotrimoxazole, OAT inhibitors."
  )

  ini({
    # Structural parameters (Perrottet 2009 Table 3, final model with IOV)
    # CL parameterisation: CL = theta_graft x GFR(L/h) x theta_female^SEXF
    #   theta_kidney = 1.68, theta_heart = 0.86, theta_lung/liver = 1.17 (L/h per L/h GFR)
    #   theta_female = 1.21 (males: SEXF = 0; females: SEXF = 1)
    # V1 parameterisation: V1 = theta_BW x (BW/70) x theta_female^SEXF
    #   theta_BW = 24 L (70-kg male reference); theta_female = 0.78
    lka              <- log(0.56);          label("Absorption rate constant ka (1/h)")                              # Table 3
    lcl              <- log(1.68);          label("CL slope on GFR for a male kidney recipient (L/h per L/h GFR)")  # Table 3 theta_kidney
    e_tx_heart_cl    <- log(0.86 / 1.68);   label("CL multiplier (log): heart vs kidney recipients")                # Table 3 theta_heart / theta_kidney
    e_tx_lungliver_cl <- log(1.17 / 1.68); label("CL multiplier (log): lung or liver vs kidney recipients (Perrottet 2009 pools both groups under one CL slope)")  # Table 3 theta_lung/liver / theta_kidney
    e_sexf_cl        <- log(1.21);          label("CL multiplier (log): female vs male recipients")                 # Table 3 theta_female on CL
    lvc              <- log(24);            label("Central volume V1 for a 70-kg male reference (L)")               # Table 3 theta_BW
    e_sexf_vc        <- log(0.78);          label("V1 multiplier (log): female vs male recipients")                 # Table 3 theta_female on V1
    lq               <- log(4.1);           label("Intercompartmental clearance Q (L/h)")                           # Table 3
    lvp              <- log(22);            label("Peripheral volume V2 (L)")                                       # Table 3
    lfdepot          <- fixed(log(0.6));    label("Oral ganciclovir-equivalent bioavailability from valganciclovir (fraction)")  # Methods (Structural model): fixed at 0.6 per refs 6, 10, 19, 24

    # Inter-patient variability (omega^2 = log(CV^2 + 1))
    # CL: 26% CV  -> omega^2 = log(0.26^2 + 1) = 0.06544
    # V1: 20% CV  -> omega^2 = log(0.20^2 + 1) = 0.03922
    etalcl ~ 0.06544          # Table 3, 26% CV interpatient on CL
    etalvc ~ 0.03922          # Table 3, 20% CV interpatient on V1

    # Residual error: proportional only after IOV inclusion
    # 21% CV proportional residual error (Table 3, sigma_prop)
    propSd <- 0.21;            label("Proportional residual error (fraction)")                                       # Table 3 sigma_prop
  })

  model({
    # 1. Derived covariate terms
    # De-index BSA-normalised MDRD-eGFR (mL/min/1.73 m^2) into individual GFR (L/h)
    # using each subject's BSA. For BSA = 1.73 m^2 the conversion reduces to CRCL * 0.06.
    gfr_lh <- CRCL * BSA / 1.73 * 60 / 1000

    # 2. Individual parameters
    ka <- exp(lka)
    cl <- exp(lcl
              + e_tx_heart_cl * TX_HEART
              + e_tx_lungliver_cl * (TX_LUNG + TX_LIVER)
              + e_sexf_cl * SEXF
              + etalcl) * gfr_lh
    vc <- exp(lvc + e_sexf_vc * SEXF + etalvc) * (WT / 70)
    q  <- exp(lq)
    vp <- exp(lvp)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. ODE system (two-compartment with first-order absorption from depot)
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Bioavailability anchor on the depot (oral valganciclovir prodrug)
    f(depot) <- exp(lfdepot)

    # 6. Observation and error
    # Dose in mg, volume in L -> mg/L
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
