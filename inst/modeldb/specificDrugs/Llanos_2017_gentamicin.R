Llanos_2017_gentamicin <- function() {
  description <- "Two-compartment population PK model of gentamicin in pediatric oncology patients with febrile neutropenia (Llanos-Paez 2017)"
  reference <- "Llanos-Paez CC, Staatz CE, Lawson R, Hennig S. A Population Pharmacokinetic Model of Gentamicin in Pediatric Oncology Patients To Facilitate Personalized Dosing. Antimicrob Agents Chemother. 2017;61(8):e00205-17. doi:10.1128/AAC.00205-17"
  vignette <- "Llanos_2017_gentamicin"
  units <- list(time = "hr", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Reference 70 kg. Drives the FFM allometric exponent of 0.75 inside GFRmat (CL pathway) and the FFM/70 multipliers on V1, V2 (linear) and Q (^0.75). Fat-free mass should be precomputed via the Janmahasatian et al. 2005 formula from total body weight, height, and sex.",
      source_name        = "FFM"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age in weeks / 4.35 + postnatal age in months)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. The Llanos 2017 paper uses postmenstrual age in WEEKS in its Rhodin et al. 2009 maturation function (TM50 = 55.4 weeks, Hill = 3.33). Canonical PAGE in nlmixr2lib is in months, so the model converts internally as PMA_weeks = PAGE * 4.35 before evaluating the Hill equation.",
      source_name        = "PMA"
    ),
    CREAT = list(
      description        = "Observed serum creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-patient measurement. Per Llanos 2017 Methods, values below the lowest reportable assay limit (< 30 umol/L) were imputed to the age- and sex-expected mean creatinine (CREAT_REF) in the model-building dataset; downstream users may apply the same imputation. Effect on CL is the inverse power model (CREAT_REF / CREAT)^0.55, so a higher individual creatinine reduces CL.",
      source_name        = "Scri"
    ),
    CREAT_REF = list(
      description        = "Age- and sex-expected mean serum creatinine for an individual, computed externally per Ceriotti et al. 2008 (Clin Chem 54:559-566)",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-patient reference value used to standardize CREAT. Time-varying through age dependence. The Ceriotti 2008 reference table is age- and (in adolescents) sex-stratified; the user precomputes CREAT_REF externally and supplies it as a covariate column. Same units as CREAT so the ratio is dimensionless.",
      source_name        = "Scrmean"
    )
  )

  population <- list(
    n_subjects     = 423L,
    n_studies      = 1L,
    age_range      = "Postnatal age 0.2-18.2 years (median 5.18 years); postmenstrual age 50.9-985 weeks (median 309 weeks)",
    weight_range   = "Total body weight 4.8-102.8 kg (median 19.4 kg); fat-free mass 3.4-72.6 kg (median 15.7 kg)",
    sex_female_pct = 48.2,
    race_ethnicity = "Not reported",
    disease_state  = "Pediatric oncology (febrile neutropenia in 88%, fever-only neutropenia in 12%); 24% had received prior nephrotoxic chemotherapy (cisplatin or carboplatin) within 6 months prior to gentamicin",
    dose_range     = "7.5 mg/kg once daily for patients < 10 years; 6 mg/kg once daily for patients >= 10 years; 30-min IV infusion (per local hospital guidelines)",
    regions        = "Australia (single center: The Lady Cilento Children's Hospital, Brisbane)",
    notes          = "Retrospective therapeutic drug monitoring data 2008-2013 (model-building) and 2014-2015 (n = 52, external evaluation). 2,422 gentamicin concentrations from 423 patients, sampled 0.5-36.0 h after end of infusion. 15% of measurements were below the lower limit of quantitation and replaced by LLOQ/2."
  )

  ini({
    # Structural parameters from Llanos 2017 Table 2 (Final model column).
    # Typical values are normalized to a 70 kg adult-equivalent reference subject
    # via the FFM/70 allometric scaling and the GFRmat asymptote 112 mL/min.
    lcl <- log(5.77); label("Clearance (CL, L/h per 70 kg reference)")                                  # Llanos 2017 Table 2
    lvc <- log(21.6); label("Central volume of distribution (V1, L per 70 kg reference)")              # Llanos 2017 Table 2
    lvp <- log(13.8); label("Peripheral volume of distribution (V2, L per 70 kg reference)")           # Llanos 2017 Table 2
    lq  <- log(0.62); label("Intercompartmental clearance (Q, L/h per 70 kg reference)")               # Llanos 2017 Table 2

    # Covariate effect coefficients (Llanos 2017 Table 2 footer; Methods Eq. 3-4).
    e_creat_cl  <- 0.55; label("Power exponent of CREAT_REF/CREAT on CL (theta_Scr, unitless)")        # Llanos 2017 Table 2
    e_ffm_cl_q  <- 0.75; label("Shared FFM/70 allometric exponent on CL (via GFRmat) and on Q (unitless)") # Llanos 2017 Methods Eq. 4
    e_ffm_vc_vp <- 1;    label("Shared FFM/70 allometric exponent on V1 and V2 (unitless)")            # Llanos 2017 Table 2 footer
    pma_hill    <- 3.33; label("Hill coefficient for postmenstrual-age GFR maturation (unitless)")     # Llanos 2017 Methods Eq. 3 (Rhodin 2009)
    pma_tm50    <- 55.4; label("Postmenstrual age at half-maximal GFR maturation (weeks)")             # Llanos 2017 Methods Eq. 3 (Rhodin 2009)
    gfr_max     <- 112;  label("Asymptotic mature glomerular filtration rate (mL/min, 70 kg adult)")   # Llanos 2017 Methods Eq. 3 (Rhodin 2009)
    gfr_ref     <- 100;  label("GFRmat reference value used to normalize CL covariate effect (mL/min)")# Llanos 2017 Table 2 footer

    # Inter-individual variability (omega^2 = log(CV^2 + 1) for log-normal).
    # Llanos 2017 Table 2 final-model BSV: CL 16.0% CV, V1 21.5% CV, V2 62.4% CV.
    # CL and V1 are correlated (correlation 69.2%); V2 is independent.
    # cov(CL, V1) = 0.692 * sqrt(0.02527 * 0.04518) = 0.02338.
    etalcl + etalvc ~ c(0.02527,
                        0.02338, 0.04518) # Llanos 2017 Table 2 (BSV CL, cov CL:V1, BSV V1)
    etalvp ~ 0.32870                      # Llanos 2017 Table 2 (BSV V2, CV 62.4%)

    # Residual error (Llanos 2017 Table 2 final model).
    propSd <- 0.275; label("Proportional residual error (SD, fraction)")                                # Llanos 2017 Table 2
    addSd  <- 0.04;  label("Additive residual error (SD, mg/L)")                                        # Llanos 2017 Table 2
  })
  model({
    # Convert canonical PAGE (months) to postmenstrual age in weeks for the
    # Llanos / Rhodin maturation equation (TM50 in weeks).
    pma_wk <- PAGE * 4.35

    # GFR maturation function (Llanos 2017 Methods Eq. 3, after Rhodin et al. 2009).
    # FFM allometric scaling is baked into GFRmat itself; do NOT multiply CL by
    # (FFM/70)^0.75 again outside this term.
    gfr_mat <- ((FFM / 70)^e_ffm_cl_q) *
               (pma_wk^pma_hill / (pma_tm50^pma_hill + pma_wk^pma_hill)) *
               gfr_max

    # Inverse-power Scr standardization on CL (Llanos 2017 Methods).
    scr_eff <- (CREAT_REF / CREAT)^e_creat_cl

    # Individual PK parameters (Llanos 2017 Table 2 footer).
    cl <- exp(lcl + etalcl) * (gfr_mat / gfr_ref) * scr_eff
    vc <- exp(lvc + etalvc) * (FFM / 70)^e_ffm_vc_vp
    vp <- exp(lvp + etalvp) * (FFM / 70)^e_ffm_vc_vp
    q  <- exp(lq)           * (FFM / 70)^e_ffm_cl_q

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Concentration: dose in mg, volume in L -> mg/L.
    Cc <- central / vc

    Cc ~ prop(propSd) + add(addSd)
  })
}
