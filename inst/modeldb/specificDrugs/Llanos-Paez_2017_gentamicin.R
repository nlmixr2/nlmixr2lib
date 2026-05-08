`Llanos-Paez_2017_gentamicin` <- function() {
  description <- "Two-compartment population PK model for gentamicin in pediatric oncology patients (Llanos-Paez 2017 AAC) extended with a renal-cortex accumulation compartment and an Emax model of relative renal-function reduction (Llanos-Paez 2017 AAPS J)."
  reference <- paste(
    "Llanos-Paez CC, Staatz CE, Hennig S.",
    "Balancing Antibacterial Efficacy and Reduction in Renal Function to Optimise",
    "Initial Gentamicin Dosing in Paediatric Oncology Patients.",
    "AAPS J. 2017;20(1):14. doi:10.1208/s12248-017-0173-6.",
    "PK structure adapted from Llanos-Paez CC, Staatz CE, Lawson R, Hennig S.",
    "A Population Pharmacokinetic Model of Gentamicin in Pediatric Oncology",
    "Patients To Facilitate Personalized Dosing.",
    "Antimicrob Agents Chemother. 2017;61(8):e00205-17. doi:10.1128/AAC.00205-17.",
    sep = " "
  )
  vignette <- "Llanos-Paez_2017_gentamicin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Linear allometric scaling on V1 and V2 with reference 70 kg;",
        "power-0.75 scaling on Q and on the GFR-maturation factor used in CL,",
        "per Llanos-Paez 2017 AAC Table 2 footnote."
      ),
      source_name        = "FFM"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age in weeks / 4.35 + postnatal age in months)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Drives the Hill-type GFR-maturation function. The source",
        "paper reports PMA in weeks; this model converts the canonical PAGE",
        "(months) back to weeks via PMA_wk = PAGE * 4.35 so the published",
        "55.4-week half-maturation and exponent of 3.33 apply unchanged."
      ),
      source_name        = "PMA (weeks)"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives a power-form effect on CL: (CREAT_REF / CREAT)^0.55 with",
        "CREAT_REF fixed at 34 umol/L (the typical-patient value reported in",
        "Llanos-Paez 2018 Methods, 'Evaluation of Gentamicin Accumulation on",
        "the Long-Term Effect on Patients Renal Function'). The original",
        "Llanos-Paez 2017 AAC parameterisation uses (Scr_mean_i / Scr_i)^0.55",
        "where Scr_mean_i is per-subject expected physiological creatinine",
        "from the Ceriotti et al. (2008) age/sex formula; the single",
        "fixed-reference simplification is documented as a deviation in the",
        "validation vignette."
      ),
      source_name        = "Scr"
    )
  )

  population <- list(
    n_subjects     = 475,
    n_studies      = 1,
    age_range      = "0.2-18.2 years (postnatal age)",
    age_median     = "5.2 years (postnatal age)",
    weight_range   = "4.5-121.0 kg",
    weight_median  = "19.5 kg",
    ffm_range      = "3.7-64.8 kg",
    ffm_median     = "15.3 kg",
    pma_range      = "Llanos-Paez 2017 AAC building cohort spans PMA 50.9-985 weeks (median 309 weeks); Llanos-Paez 2018 reuses the same population.",
    sex_female_pct = "48% female (Llanos-Paez 2017 AAC building cohort, 204/423); sex distribution of the 52-subject external evaluation cohort 58% female; Llanos-Paez 2018 does not re-report sex separately for the combined 475-patient cohort.",
    race_ethnicity = "Not reported.",
    disease_state  = "Pediatric oncology patients (febrile or fever-only neutropenia). Eighty-eight percent of the building cohort had febrile neutropenia and 12% had fever-only neutropenia.",
    dose_range     = "Local clinical-practice initial dose 7.5 mg/kg/q24h (<10 yrs) or 6 mg/kg/q24h (>=10 yrs), administered as a 30-min IV infusion. Llanos-Paez 2018 evaluates initial doses 7.1, 9.5, 10.8, 12.8, and 14.6 mg/kg/q24h.",
    regions        = "Australia (Lady Cilento Children's Hospital, Brisbane).",
    notes          = paste(
      "475-patient pediatric oncology cohort comprising the 423 model-development",
      "subjects and 52-subject external evaluation cohort from Llanos-Paez 2017",
      "AAC. Llanos-Paez 2018 reuses the AAC PK model and parameter estimates",
      "and extends the model with a renal-cortex accumulation compartment and",
      "an Emax renal-function reduction model (Eqs 10-14 and Table I).",
      "Demographic ranges are reported across the two papers: Llanos-Paez 2017",
      "AAC Table 1 and Llanos-Paez 2018 'Patients and PK Model' Results section."
    )
  )

  ini({
    # ===== Structural PK (Llanos-Paez 2017 AAC Table 2, final model) =====
    # Reference subject: 70 kg fat-free mass at full GFR maturation with Scr = Scr_mean_ref.
    lcl <- log(5.77);  label("Typical CL at FFM = 70 kg, GFR_mat = 100, CREAT = 34 umol/L (L/h)")  # AAC Table 2 final model: CL = 5.77 L/h/70 kg
    lvc <- log(21.6);  label("Typical V1 (central volume of distribution) at FFM = 70 kg (L)")    # AAC Table 2 final model: V1 = 21.6 L/70 kg
    lq  <- log(0.62);  label("Typical Q (intercompartmental clearance) at FFM = 70 kg (L/h)")     # AAC Table 2 final model: Q = 0.62 L/h/70 kg
    lvp <- log(13.8);  label("Typical V2 (peripheral volume of distribution) at FFM = 70 kg (L)") # AAC Table 2 final model: V2 = 13.8 L/70 kg

    # Allometric / maturation / covariate exponents (Llanos-Paez 2017 AAC Table 2 footnote)
    e_ffm_cl_q  <- 0.75;  label("Shared FFM allometric exponent on Q and on the GFR_mat factor")
    e_creat_cl  <- 0.55;  label("Power exponent on (CREAT_REF / CREAT) for CL")  # AAC Table 2 final model: theta_Scr = 0.55

    # GFR maturation Hill function constants (Llanos-Paez 2017 AAC Table 2 footnote)
    pma50_gfr     <- 55.4;  label("Postmenstrual age at half-maximal GFR maturation (weeks)")
    h_pma_gfr     <- 3.33;  label("Hill coefficient on PMA for GFR maturation (unitless)")
    gfr_max_adult <- 112;   label("Adult-equivalent maximum GFR at full maturation (mL/min)")
    creat_ref     <- 34;    label("Reference serum creatinine for the (CREAT_REF / CREAT) ratio (umol/L)")  # Llanos-Paez 2018 'Evaluation of Gentamicin Accumulation' Methods: typical-patient Scr = 34 umol/L
    gfr_ref_cl    <- 100;   label("Reference GFR_mat used in the CL covariate model (mL/min)")

    # ===== Renal-cortex accumulation (Llanos-Paez 2018 Eq 10, Table I) =====
    lvmax_rc   <- log(1.0);    label("Maximum gentamicin uptake rate into renal cortex (mg/h)")  # Llanos-Paez 2018 Table I: V_max = 1.0 mg/h (citing Giuliano et al. 1986, ref 27)
    lkm_rc     <- log(15.0);   label("Michaelis constant for renal-cortex uptake (mg)")           # Llanos-Paez 2018 Table I: k_m = 15.0 mg (citing Rougier et al. 2003, ref 10)
    lkreabs_rc <- log(0.0693); label("Tubular reabsorption rate from renal cortex (1/h)")         # Llanos-Paez 2018 Table I: k_reabs = 0.0693 1/h (citing Croes et al. 2011, ref 11)

    # ===== Renal-function Emax model (Llanos-Paez 2018 Eq 11, 13; Table I) =====
    emax_erf <- 100;   label("Maximum E_RF (renal-function-reduction utility-scale variable, %)")  # Llanos-Paez 2018 Table I: E_max = 100 (ref 10)
    arc50    <- 50;    label("Renal-cortex amount producing E_RF = E_max / 2 (mg)")                # Llanos-Paez 2018 Table I: A_RC50 = 50 mg (ref 10)
    h_arc    <- 5;     label("Hill coefficient on A_Rc for E_RF (unitless)")                       # Llanos-Paez 2018 Table I: gamma = 5 (citing Croes et al., ref 11)
    erf50    <- 33.5;  label("E_RF value at which renal function further declines by half-max")    # Llanos-Paez 2018 Table I: E_GFR50 = 33.5 (ref 10)
    h_erf    <- 5.5;   label("Hill coefficient on E_RF for the renal-function reduction (unitless)")  # Llanos-Paez 2018 Table I: delta = 5.5 (ref 10)
    f_gfrmax <- 0.41;  label("Maximum fractional reduction in renal function (GFR_max / GFR_obs)")  # Llanos-Paez 2018 Table I: GFR_max = 41% of GFR_obs (ref 10)

    # ===== IIV (Llanos-Paez 2017 AAC Table 2 final model) =====
    # CV%-to-variance: omega^2 = log(1 + CV^2)
    # CL CV 16.0% -> 0.02528, V1 CV 21.5% -> 0.04525, corr(CL,V1) = 0.692,
    #   cov = 0.692 * sqrt(0.02528 * 0.04525) = 0.02340
    etalcl + etalvc ~ c(0.02528,
                        0.02339, 0.04519)  # AAC Table 2: BSV CL 16.0%, BSV V1 21.5%, corr 69.2%
    etalq  ~ 0.03119   # AAC Table 2: BSV Q  17.8% CV  -> log(1 + 0.178^2)
    etalvp ~ 0.32885   # AAC Table 2: BSV V2 62.4% CV  -> log(1 + 0.624^2)

    # ===== Residual error =====
    propSd <- 0.275;  label("Proportional residual error (fraction)")  # AAC Table 2: 27.5%
    addSd  <- 0.04;   label("Additive residual error (mg/L)")          # AAC Table 2: 0.04 mg/L
  })

  model({
    # ----- Derived covariate terms -----
    # Convert canonical PAGE (months) back to source-paper PMA (weeks) so the
    # published 55.4-week half-maturation and 3.33 Hill coefficient apply unchanged.
    pma_wk <- PAGE * 4.35

    # GFR maturation (Llanos-Paez 2017 AAC Table 2 footnote):
    #   GFR_mat (mL/min) = 112 * (FFM/70)^0.75 * PMA_wk^3.33 / (55.4^3.33 + PMA_wk^3.33)
    gfr_mat <- gfr_max_adult * (FFM / 70)^e_ffm_cl_q *
      pma_wk^h_pma_gfr / (pma50_gfr^h_pma_gfr + pma_wk^h_pma_gfr)

    # Creatinine effect on CL (Llanos-Paez 2017 AAC Table 2 footnote: (Scr_mean / Scr)^0.55).
    # See covariateData[[CREAT]]$notes and the vignette deviations: per-subject Scr_mean
    # from Ceriotti 2008 is replaced by a fixed cohort-typical reference (creat_ref).
    creat_factor <- (creat_ref / CREAT)^e_creat_cl

    # ----- Individual PK parameters -----
    cl <- exp(lcl + etalcl) * (gfr_mat / gfr_ref_cl) * creat_factor
    vc <- exp(lvc + etalvc) * (FFM / 70)
    q  <- exp(lq  + etalq)  * (FFM / 70)^e_ffm_cl_q
    vp <- exp(lvp + etalvp) * (FFM / 70)

    # ----- Renal-cortex parameters -----
    vmax_rc   <- exp(lvmax_rc)
    km_rc     <- exp(lkm_rc)
    kreabs_rc <- exp(lkreabs_rc)

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- ODE system (Llanos-Paez 2018 Fig 2 schematic; Eq 10) -----
    # Two-compartment IV PK + saturable uptake into the renal cortex with first-order
    # tubular reabsorption back out of the cortex. IV doses go directly to `central`.
    d/dt(central)      <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1)  <-  k12 * central - k21 * peripheral1
    d/dt(renal_cortex) <- -kreabs_rc * renal_cortex + vmax_rc * central / (km_rc + central)

    # ----- Outputs -----
    # Plasma concentration: dose in mg, vc in L -> mg/L
    Cc <- central / vc

    # Renal-function utility-scale variable (Llanos-Paez 2018 Eq 11):
    #   E_RF(t) = E_max * A_Rc^gamma / (A_RC50^gamma + A_Rc^gamma)
    erf <- emax_erf * renal_cortex^h_arc / (arc50^h_arc + renal_cortex^h_arc)

    # Relative reduction in renal function (Llanos-Paez 2018 Eqs 12-14, simplified):
    #   rDeltaRF = (RF_0 - RF_new) / RF_0 = (GFR_obs - GFR_new) / GFR_obs
    #            = f_gfrmax * E_RF^delta / (E_RF50^delta + E_RF^delta)
    # Independent of GFR_obs because GFR_max is parameterised as a fixed fraction
    # of GFR_obs (Llanos-Paez 2018 Table I, GFR_max = 41% of GFR_obs).
    rDeltaRF <- f_gfrmax * erf^h_erf / (erf50^h_erf + erf^h_erf)

    # ----- Error model -----
    Cc ~ add(addSd) + prop(propSd)
  })
}
