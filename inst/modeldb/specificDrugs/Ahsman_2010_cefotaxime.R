Ahsman_2010_cefotaxime <- function() {
  description <- paste0(
    "One-compartment population PK model for cefotaxime (CTX) and its active ",
    "metabolite desacetylcefotaxime (DACT) in critically ill neonates and ",
    "infants on extracorporeal membrane oxygenation (ECMO). IV bolus parent ",
    "with first-order elimination; the metabolite is generated 1:1 from ",
    "parent elimination on a CTX-equivalent mass basis (the source paper ",
    "converted observed DACT concentrations to CTX equivalents by the ",
    "molecular weight ratio Mr_CTX / Mr_DACT = 455.5 / 413.4 before fitting, ",
    "and assumed a conversion fraction FDACT/CTX = 1). Parent CL is scaled by ",
    "body weight (WT) via a power model centred at 3.5 kg; both parent and ",
    "metabolite CL include a power covariate effect of time after ECMO ",
    "decannulation (T_POST_ECMO) centred at 100 h that is removed during ",
    "ECMO (T_POST_ECMO = 0). Metabolite CL additionally includes a power ",
    "covariate effect of continuous venovenous hemofiltration flow Q_CVVH ",
    "centred at 193 mL/min that is removed when CVVH is not running ",
    "(Q_CVVH = 0). Volumes (Vc, Vc_dact) carry no covariate effects."
  )
  reference <- paste(
    "Ahsman MJ, Wildschut ED, Tibboel D, Mathot RA.",
    "Pharmacokinetics of cefotaxime and desacetylcefotaxime in infants",
    "during extracorporeal membrane oxygenation.",
    "Antimicrob Agents Chemother. 2010;54(5):1734-1741.",
    "doi:10.1128/AAC.01696-09.",
    sep = " "
  )
  vignette <- "Ahsman_2010_cefotaxime"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject in the source analysis (the most recent weight before ECMO was used for both dose calculation and PK analysis). Power covariate on CTX CL centred at 3.5 kg (cohort median, Table 1).",
      source_name        = "WT"
    ),
    Q_CVVH = list(
      description        = "Continuous venovenous hemofiltration (CVVH) circuit flow rate when the CVVH filter is active in the ECMO circuit",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; 0 when CVVH is not running. Power covariate on DACT CL centred at 193 mL/min (cohort median during CVVH operation, Table 1). When Q_CVVH = 0 (no CVVH active) the covariate effect is removed from the equation (multiplier = 1), per the source paper.",
      source_name        = "QCVVH"
    ),
    T_POST_ECMO = list(
      description        = "Time after ECMO decannulation (end of extracorporeal circulation)",
      units              = "hour",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; 0 before and during ECMO, becomes the wall-clock hours elapsed since decannulation thereafter. Power covariate on both CTX CL and DACT CL centred at 100 h (Appendix). When T_POST_ECMO = 0 the covariate effect is removed from the equation (multiplier = 1), per the source paper.",
      source_name        = "tEND"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 37L,
    n_studies      = 1L,
    age_range      = "postnatal age 0.67-199 days at ECMO start (median 3.3 days); gestational age 34-42 weeks at birth (median 37 weeks)",
    age_median     = "3.3 days postnatal",
    weight_range   = "2.0-6.2 kg at ECMO start",
    weight_median  = "3.5 kg",
    sex_female_pct = 51.4,
    race_ethnicity = NULL,
    disease_state  = "Critically ill neonates and infants on ECMO. Primary diagnoses: meconium aspiration syndrome (46%), congenital diaphragmatic hernia (22%), pulmonary hypertension other causes (14%), congenital heart defects (11%), other including sepsis and viral infections (7%).",
    dose_range     = "Standard IV bolus regimen: 50 mg/kg BID (PNA < 1 week, 65%), 50 mg/kg TID (PNA 1-4 weeks, 19%), 37.5 mg/kg QID (PNA > 4 weeks, 8%); small subgroup on 25 mg/kg BID (5%) or 37.5 mg/kg TID (3%) per clinical discretion.",
    regions        = "Single centre, Sophia Children's Hospital, Erasmus University Medical Center, Rotterdam, Netherlands.",
    notes          = "December 2006-June 2009. ECMO modality 54% VV / 46% VA (n = 41 ECMO runs across 37 patients; 4 patients had 2 runs each). Median ECMO duration 108 h (range 16-374 h). 30/37 patients underwent CVVH at median flow 193 mL/min (range 100-350 mL/min). Hypothermic (24 deg C) n = 2; normothermic (36 deg C) n = 35. Serum chemistry medians: albumin 31 g/L (21-40), serum creatinine 32 umol/L (19-69). Survival 25/37. 392 plasma samples (median 10 per patient, range 1-17). Cefotaxime administered IV into the extracorporeal line distal to the oxygenator. DACT concentrations were converted to CTX equivalents at dataset-build time via the molecular weight ratio Mr_CTX / Mr_DACT = 455.5 / 413.4."
  )

  ini({
    # ----------------------------------------------------------------------
    # Structural parameters -- Table 2 final-model column.
    # CTX (parent) and DACT (metabolite) each described by a one-compartment
    # model with first-order elimination; the metabolite mass-balance assumes
    # 1:1 conversion of CTX to DACT (FDACT/CTX = 1, footnote a of Table 2),
    # with concentrations carried in CTX-equivalent units (the source paper
    # converted observed DACT by Mr_CTX / Mr_DACT = 455.5 / 413.4 before
    # fitting). Reference weight 3.5 kg = cohort median (Table 1).
    # ----------------------------------------------------------------------
    lvc       <- log(1.82);  label("CTX central volume of distribution Vc (L)")                                  # Table 2
    lcl       <- log(0.36);  label("CTX clearance CL at reference WT = 3.5 kg, T_POST_ECMO = 0 (L/h)")           # Table 2
    lvc_dact  <- log(11.0);  label("DACT apparent central volume of distribution Vc_dact in CTX equivalents (L); assumes FDACT/CTX = 1") # Table 2
    lcl_dact  <- log(1.46);  label("DACT apparent clearance CL_dact at reference Q_CVVH = 193 mL/min, T_POST_ECMO = 0 in CTX equivalents (L/h); assumes FDACT/CTX = 1") # Table 2

    # ----------------------------------------------------------------------
    # Covariate effects -- Table 2 "Covariate effects" rows and Appendix
    # equations. All are power exponents centred at the reference values
    # carried by the corresponding covariate column. The CVVH and T_POST_ECMO
    # effects are removed from the equation when the covariate is zero
    # (Q_CVVH = 0 means no CVVH active; T_POST_ECMO = 0 means before or
    # during ECMO); see model() for the conditional implementation.
    # ----------------------------------------------------------------------
    e_wt_cl                <- 0.56; label("Power exponent of WT on CTX CL (reference 3.5 kg, unitless)")                              # Table 2
    e_q_cvvh_cl_dact       <- 0.72; label("Power exponent of Q_CVVH on DACT CL (reference 193 mL/min, unitless)")                    # Table 2
    e_t_post_ecmo_cl       <- 0.16; label("Power exponent of T_POST_ECMO on CTX CL (reference 100 h, unitless)")                     # Table 2
    e_t_post_ecmo_cl_dact  <- 0.53; label("Power exponent of T_POST_ECMO on DACT CL (reference 100 h, unitless)")                    # Table 2

    # ----------------------------------------------------------------------
    # Inter-individual variability -- Table 2 IIV rows.
    # Variances on the log scale computed from the reported %CV via
    # omega^2 = log(CV^2 + 1). The source paper additionally estimated
    # an omega-block correlation structure between CL and V (correlation
    # range 70.6% for CL_DACT-Vc_dact to 90.8% for Vc-Vc_dact); only those
    # two boundary correlations are reported, the full 4x4 correlation
    # matrix is not given. IIV is therefore encoded here as diagonal
    # (independent etas) and the omitted correlations are documented in the
    # vignette Assumptions and deviations section.
    # ----------------------------------------------------------------------
    etalvc      ~ 0.118064; label("IIV variance on log Vc (35.4% CV)")              # Table 2
    etalcl      ~ 0.122502; label("IIV variance on log CL (36.1% CV)")              # Table 2
    etalvc_dact ~ 0.305721; label("IIV variance on log Vc_dact (59.8% CV)")         # Table 2
    etalcl_dact ~ 0.234436; label("IIV variance on log CL_dact (51.4% CV)")         # Table 2

    # ----------------------------------------------------------------------
    # Residual error -- Table 2 "Residual variability" rows.
    # The source paper estimated TWO proportional residual error terms, one
    # for samples with t_dose < 1 h (69.4% CV) and one for samples with
    # t_dose > 1 h (32.7% CV); the early-time term was introduced to absorb
    # extra variability from nurse-recorded dose times vs the actual
    # injection time, not as a structural feature of CTX or DACT PK. The
    # tabulated values are identical between the CTX and DACT columns
    # (shared proportional variance across the two outputs). For the
    # library model we keep the t_dose > 1 h proportional residual (the
    # representative steady-state value, 32.7% CV) and document the
    # collapsed time-window structure in the vignette.
    # ----------------------------------------------------------------------
    propSd      <- 0.327; label("Proportional residual error on CTX (fraction, t_dose > 1 h)")   # Table 2
    propSd_dact <- 0.327; label("Proportional residual error on DACT (fraction, t_dose > 1 h)") # Table 2
  })

  model({
    # ------------------------------------------------------------------
    # Reference (centring) values for the covariate effects -- Table 2
    # remarks and Appendix.
    # ------------------------------------------------------------------
    ref_wt        <- 3.5    # kg, cohort median weight at ECMO start (Table 1)
    ref_q_cvvh    <- 193    # mL/min, cohort median CVVH flow (Table 1)
    ref_t_post_ecmo <- 100  # h, centring constant for T_POST_ECMO power covariate (Appendix)

    # ------------------------------------------------------------------
    # Individual parameters -- exp(lX + etalX) * covariate multipliers.
    # The Q_CVVH and T_POST_ECMO power terms are gated to 1 when the
    # covariate is zero so the (0/ref)^exponent expression does not
    # collapse to zero (the source paper states "when Q_CVVH = 0" and
    # "when t_END = 0" the corresponding covariate effect is removed
    # from the equation, i.e. multiplier = 1).
    # ------------------------------------------------------------------
    cov_wt_cl              <- (WT / ref_wt)^e_wt_cl
    cov_t_post_ecmo_cl     <- ifelse(T_POST_ECMO > 0, (T_POST_ECMO / ref_t_post_ecmo)^e_t_post_ecmo_cl,      1)
    cov_t_post_ecmo_cl_d   <- ifelse(T_POST_ECMO > 0, (T_POST_ECMO / ref_t_post_ecmo)^e_t_post_ecmo_cl_dact, 1)
    cov_q_cvvh_cl_dact     <- ifelse(Q_CVVH       > 0, (Q_CVVH       / ref_q_cvvh)^e_q_cvvh_cl_dact,            1)

    vc      <- exp(lvc      + etalvc)
    cl      <- exp(lcl      + etalcl)      * cov_wt_cl * cov_t_post_ecmo_cl
    vc_dact <- exp(lvc_dact + etalvc_dact)
    cl_dact <- exp(lcl_dact + etalcl_dact) * cov_q_cvvh_cl_dact * cov_t_post_ecmo_cl_d

    # ------------------------------------------------------------------
    # Micro-constants. AMTCMT1 / AMTCMT2 in the Appendix correspond to
    # central (CTX) and central_dact (DACT) respectively.
    # ------------------------------------------------------------------
    kel      <- cl      / vc
    kel_dact <- cl_dact / vc_dact

    # ------------------------------------------------------------------
    # ODEs (Appendix).
    # IV bolus dosing on `central` (no depot). The CTX -> DACT transfer is
    # carried as a flux equal to CL_CTX * Cp_CTX = kel * central, consistent
    # with the assumed 1:1 (CTX-equivalent) conversion FDACT/CTX = 1.
    # ------------------------------------------------------------------
    d/dt(central)      <- -kel * central
    d/dt(central_dact) <-  kel * central - kel_dact * central_dact

    # ------------------------------------------------------------------
    # Observations and residual error.
    # Cc is the cefotaxime plasma concentration in mg/L; Cc_dact is the
    # desacetylcefotaxime plasma concentration in mg/L of CTX equivalents
    # (the source paper converted observed DACT by Mr_CTX / Mr_DACT =
    # 455.5 / 413.4 before fitting, so the modeled mass and the modeled
    # Vc_dact / CL_dact are both in CTX-equivalent units). To recover the
    # actual DACT concentration in DACT mass units, multiply Cc_dact by
    # Mr_DACT / Mr_CTX = 413.4 / 455.5.
    # ------------------------------------------------------------------
    Cc      <- central      / vc
    Cc_dact <- central_dact / vc_dact

    Cc      ~ prop(propSd)
    Cc_dact ~ prop(propSd_dact)
  })
}
