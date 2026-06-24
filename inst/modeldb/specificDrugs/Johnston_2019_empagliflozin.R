Johnston_2019_empagliflozin <- function() {
  description <- paste0(
    "Exposure-response (PD-only) model for the effect of empagliflozin on ",
    "HbA1c in patients with type 1 diabetes mellitus (T1DM) on background ",
    "insulin therapy (M-EASE-2; Johnston 2019). A direct-response Emax ",
    "function of individual steady-state empagliflozin AUC (AUC_EMPA, ",
    "supplied as a per-subject covariate column from an upstream popPK ",
    "analysis -- Mondick 2018 plus EASE-2 / EASE-3 data-on-file) reduces ",
    "the model-predicted baseline HbA1c, with an additional linear placebo ",
    "drift over time. Full covariate model on baseline HbA1c, Emax, and ",
    "placebo (sex, insulin delivery type, body weight, eGFR, baseline ",
    "insulin daily dose, and -- on Emax only -- baseline HbA1c)."
  )
  reference <- paste(
    "Johnston CK, Riggs MM, Marquard J, Soleymanlou N, Nock V,",
    "Liesenfeld K-H. M-EASE-2: A Modelling and Simulation Study",
    "Conducted to Further Characterise the Efficacy of Low-dose",
    "Empagliflozin as Adjunctive to InSulin ThErapy (M-EASE) in Type 1",
    "Diabetes Mellitus. American Diabetes Association 79th Scientific",
    "Sessions, 2019; poster 1198-P.",
    "doi:10.2337/db19-1198-p.",
    "PDF: https://metrumrg.com/wp-content/uploads/Pubs/M-EASE-2-page2019-johnston.pdf.",
    "Upstream popPK structure cited by the source authors:",
    "Mondick J et al. J Clin Pharmacol. 2018;58:640-649",
    "(updated internally with EASE-2 / EASE-3 data on file to compute",
    "the individual AUCss inputs to this PD model)."
  )
  vignette <- "Johnston_2019_empagliflozin"
  units <- list(
    time          = "hour",
    dosing        = "n/a (no drug dosing events; empagliflozin exposure enters as the per-subject AUC_EMPA covariate from an upstream popPK)",
    concentration = "% HbA1c (NGSP; observation output -- not a drug concentration)",
    AUC_EMPA      = "nmol*h/L (empagliflozin steady-state AUC over the q24h dosing interval)"
  )

  covariateData <- list(
    AUC_EMPA = list(
      description        = "Steady-state empagliflozin AUC over the q24h dosing interval supplied as a per-subject (time-fixed) drug-exposure covariate from an upstream popPK analysis (Mondick 2018, updated with EASE-2 / EASE-3 data on file).",
      units              = "nmol*h/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-subject (time-fixed) steady-state AUC of empagliflozin. The",
        "source authors generate AUCss from individual empirical Bayes",
        "estimates of an upstream popPK model (Mondick 2018 plus EASE-2 /",
        "EASE-3 data on file) and pass them to the PD model as a static",
        "covariate. AUC_EMPA = 0 reproduces the placebo arm (drug-effect",
        "term vanishes). Reference AUCss for the 2.5 mg / 10 mg / 25 mg",
        "QD dose groups is approximately AUC50 = 498 nmol*h/L for 2.5 mg",
        "(back-calculated from the simulated half-maximal -0.29% HbA1c",
        "change reported in the poster Results)."
      ),
      source_name        = "AUCss"
    ),
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline (pre-treatment) body weight; reference 82 kg (Johnston 2019 Results / reference-patient description). Enters baseline HbA1c and Emax via power form `(WT / 82)^e_wt_<param>`.",
      source_name        = "WTB"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate (BSA-normalised eGFR)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline (pre-treatment) eGFR; reference 98 mL/min/1.73 m^2 (Johnston 2019 reference-patient description). The source paper does not specify which creatinine-based eGFR equation was used (MDRD or CKD-EPI). Enters baseline HbA1c and Emax via power form `(CRCL / 98)^e_crcl_<param>`.",
      source_name        = "eGFR"
    ),
    HBA1C = list(
      description        = "Baseline (pre-treatment) HbA1c (per-subject, time-fixed)",
      units              = "% (NGSP)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline HbA1c as a per-subject covariate on Emax in the full-random-effects covariate model (Johnston 2019 Equation 1; reference 8.14 %, the model-estimated typical baseline at the reference patient). NOT the same column as the time-course HbA1c observations: this is the per-subject baseline anchor used inside the covariate equation for Emax; observations are the model-predicted longitudinal HbA1c trajectory. Enters via power form `(HBA1C / 8.14)^e_hba1c_emax`.",
      source_name        = "HbA1c (baseline)"
    ),
    INSDOSE_BL = list(
      description        = "Baseline total daily insulin dose normalised to body weight",
      units              = "U/kg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject (time-fixed) total daily exogenous insulin dose at baseline, normalised to body weight (U/kg). Reference 0.660 U/kg (Johnston 2019 reference-patient description). Enters baseline HbA1c and Emax via power form `(INSDOSE_BL / 0.660)^e_insdose_bl_<param>`. Distinct from `INS_BL` (plasma insulin concentration) -- this is an administered dose level, not a measured biomarker.",
      source_name        = "IDB"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male reference)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; the reference patient is male per Johnston 2019 Results)",
      notes              = "Encoded as multiplier `e_sexf_<param>^SEXF` on baseline HbA1c, Emax, and placebo rate (i.e., the Table 2 multiplier value applies when SEXF = 1, and 1.0 when SEXF = 0). 391 males and 405 females in the modelled cohort.",
      source_name        = "Sex"
    ),
    INSDT_CSII = list(
      description        = "Insulin delivery type indicator (1 = continuous subcutaneous insulin infusion (CSII), 0 = multiple daily injections (MDI) reference)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (MDI; the reference patient is on MDI per Johnston 2019 Results)",
      notes              = "Encoded as multiplier `e_insdt_csii_<param>^INSDT_CSII` on baseline HbA1c, Emax, and placebo rate. Follows the `<COLUMN>_<LEVEL>` decomposition pattern used for `REGI_BID` and `TUMTP_*` categoricals: a future paper that retains the same INSDT covariate but with more than two levels would register sibling canonicals (e.g., `INSDT_CGM`) rather than overload this name. M-EASE-2 distinguishes only CSII vs MDI.",
      source_name        = "INSDT"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 796L,
    n_subjects_male = 391L,
    n_subjects_female = 405L,
    n_studies       = 2L,
    studies         = "EASE-2 (52-week phase 3 trial, empagliflozin 10 and 25 mg QD treatment arms; primary development data) and EASE-1 (4-week phase 2 dose-finding study, empagliflozin 2.5, 10, 25 mg QD). External qualification used EASE-3 (out-of-sample 26-week phase 3 trial that included a 2.5 mg QD arm).",
    age_range       = "21-69 years (95th-percentile interval at baseline)",
    weight_range    = "55-125 kg (95th-percentile interval at baseline)",
    sex_female_pct  = 100 * 405 / 796,  # = 50.88 %
    disease_state   = paste0(
      "Adults with type 1 diabetes mellitus (T1DM) on background insulin ",
      "therapy. Baseline HbA1c 7.2-9.5 % (95th-percentile interval). ",
      "Reference patient: male, multiple-daily-injection (MDI) insulin, ",
      "total daily insulin dose 0.660 U/kg, HbA1c 8.1 %, eGFR 98 mL/min/",
      "1.73 m^2, body weight 82 kg."
    ),
    dose_range      = "Empagliflozin 0 (placebo), 2.5 (simulated only -- not directly studied in EASE-2), 10, 25 mg QD adjunctive to insulin; per-subject AUCss supplied as the AUC_EMPA covariate.",
    egfr_range      = "57-127 mL/min/1.73 m^2 (95th-percentile interval at baseline)",
    regions         = "Multi-national (EASE-2 / EASE-1 trials).",
    n_observations  = "Not reported in the poster.",
    notes           = paste0(
      "Pop-PD analysis using Markov Chain Monte Carlo Bayesian estimation ",
      "(NONMEM v7.4). AUC50 was estimated under an informative prior ",
      "derived from a Type 2 diabetes exposure-response analysis (Baron ",
      "2016, ref 6 in the poster) because the M-EASE-2 cohort had limited ",
      "data at empagliflozin 2.5 mg. Sensitivity analyses varied the ",
      "prior's variance and location (Johnston 2019 Figure 3). The ",
      "primary objective was to simulate the placebo-corrected HbA1c ",
      "change from baseline up to 52 weeks for a hypothetical 2.5 mg QD ",
      "arm in the EASE-2 population. Source: Johnston 2019 Results / ",
      "Data-study-population section."
    )
  )

  ini({
    # ---- Structural parameters (Johnston 2019 Table 2 'Estimate' column;
    #      reference patient: male, MDI, 0.660 U/kg, 8.14 % HbA1c,
    #      98 mL/min/1.73 m^2, 82 kg) ----
    lhba1c0      <- log(8.14)   ; label("Typical baseline HbA1c at the reference patient (% NGSP)")             # Johnston 2019 Table 2 (95% CI 8.07, 8.22)
    lauc50       <- log(498)    ; label("Typical AUC50 -- empagliflozin AUCss at half-maximal effect (nmol*h/L)") # Johnston 2019 Table 2 (95% CI 296, 819); informative prior from Baron 2016 T2DM analysis
    lemax        <- log(0.579)  ; label("Typical maximum HbA1c reduction Emax at the reference patient (%)")     # Johnston 2019 Table 2 (95% CI 0.491, 0.678)
    lkplacebo    <- log(2.61e-5); label("Typical placebo HbA1c drift rate per hour at the reference patient (%/h)") # Johnston 2019 Table 2 (95% CI 1.96e-5, 3.29e-5)

    # ---- Covariate effects on baseline HbA1c (Johnston 2019 Table 2;
    #      power exponents for continuous, multiplier^indicator for binary) ----
    e_wt_hba1c0          <- -0.0311 ; label("Power exponent for body weight on baseline HbA1c (unitless)")                                       # Johnston 2019 Table 2 (95% CI -0.0612, -0.00102)
    e_crcl_hba1c0        <-  0.0123 ; label("Power exponent for eGFR on baseline HbA1c (unitless)")                                              # Johnston 2019 Table 2 (95% CI -0.0157, 0.0403)
    e_insdose_bl_hba1c0  <-  0.0141 ; label("Power exponent for baseline insulin daily dose on baseline HbA1c (unitless)")                       # Johnston 2019 Table 2 (95% CI -0.00425, 0.0326)
    e_sexf_hba1c0        <-  0.988  ; label("Female / male ratio multiplier on baseline HbA1c (applied as e_sexf_hba1c0^SEXF)")                  # Johnston 2019 Table 2 (95% CI 0.977, 1.00)
    e_insdt_csii_hba1c0  <-  1.00   ; label("CSII / MDI ratio multiplier on baseline HbA1c (applied as e_insdt_csii_hba1c0^INSDT_CSII)")          # Johnston 2019 Table 2 (95% CI 0.988, 1.01)

    # ---- Covariate effects on Emax (Johnston 2019 Table 2) ----
    e_wt_emax          <-  0.0555 ; label("Power exponent for body weight on Emax (unitless)")                                  # Johnston 2019 Table 2 (95% CI -0.351, 0.458)
    e_crcl_emax        <-  0.504  ; label("Power exponent for eGFR on Emax (unitless)")                                         # Johnston 2019 Table 2 (95% CI 0.116, 0.917)
    e_insdose_bl_emax  <-  0.0552 ; label("Power exponent for baseline insulin daily dose on Emax (unitless)")                  # Johnston 2019 Table 2 (95% CI -0.190, 0.300)
    e_hba1c_emax       <-  0.999  ; label("Power exponent for baseline HbA1c on Emax (unitless)")                               # Johnston 2019 Table 2 (95% CI -0.358, 2.33)
    e_sexf_emax        <-  0.984  ; label("Female / male ratio multiplier on Emax (applied as e_sexf_emax^SEXF)")               # Johnston 2019 Table 2 (95% CI 0.827, 1.17)
    e_insdt_csii_emax  <-  0.880  ; label("CSII / MDI ratio multiplier on Emax (applied as e_insdt_csii_emax^INSDT_CSII)")       # Johnston 2019 Table 2 (95% CI 0.737, 1.04)

    # ---- Covariate effects on placebo rate (Johnston 2019 Table 2;
    #      only sex and insulin delivery type are reported) ----
    e_sexf_kplacebo       <- 0.727 ; label("Female / male ratio multiplier on placebo rate (applied as e_sexf_kplacebo^SEXF)")              # Johnston 2019 Table 2 (95% CI 0.534, 0.971)
    e_insdt_csii_kplacebo <- 1.47  ; label("CSII / MDI ratio multiplier on placebo rate (applied as e_insdt_csii_kplacebo^INSDT_CSII)")     # Johnston 2019 Table 2 (95% CI 1.10, 1.99)

    # ---- IIV (Johnston 2019 narrative under 'Model development':
    #      'Inter-individual variance (CV%) for baseline HbA1c and Emax
    #      were 7.2 % and 38 %, respectively'.
    #      log-normal convention: omega^2 = log(1 + CV^2). ----
    etalhba1c0 ~ 0.005177  # omega^2 = log(1 + 0.072^2) = 0.005177 (CV 7.2 % on baseline HbA1c)
    etalemax   ~ 0.135108  # omega^2 = log(1 + 0.38^2)  = 0.135108 (CV 38 % on Emax)

    # ---- Residual error (Johnston 2019 narrative under 'Model development':
    #      'the proportional and additive residual variability estimates
    #      (CV% and standard deviation) were 4.6 % and 0.11, respectively'). ----
    propSd <- 0.046 ; label("Proportional residual error on HbA1c (fraction)")  # Johnston 2019 narrative (4.6 % CV)
    addSd  <- 0.11  ; label("Additive residual error on HbA1c (% NGSP)")        # Johnston 2019 narrative (SD 0.11 %)
  })

  model({
    # ---- Reference (typical-patient) covariate centring values ----
    # All from Johnston 2019 Results / 'Data study population' section: the
    # 'reference patient' description (male, MDI, 0.660 U/kg, 8.14 % HbA1c,
    # 98 mL/min/1.73 m^2, 82 kg).
    ref_wt        <- 82
    ref_crcl      <- 98
    ref_hba1c     <- 8.14
    ref_insdose   <- 0.660

    # ---- Individual baseline HbA1c (Johnston 2019 Equation 1, first line):
    # Baseline_HbA1c,i = theta_a
    #                  * prod_m (cov_m / ref_m)^theta_(m+a)        [continuous]
    #                  * prod_p (theta_(p+m+a))^cov_p              [binary]
    #                  * exp(eta_z).
    baseline_hba1c_i <- exp(lhba1c0 + etalhba1c0) *
      (WT         / ref_wt     )^e_wt_hba1c0 *
      (CRCL       / ref_crcl   )^e_crcl_hba1c0 *
      (INSDOSE_BL / ref_insdose)^e_insdose_bl_hba1c0 *
      e_sexf_hba1c0^SEXF *
      e_insdt_csii_hba1c0^INSDT_CSII

    # ---- Individual Emax (Johnston 2019 Equation 1, third line):
    # Emax_i = theta_c
    #        * prod_m (cov_m / ref_m)^theta_(m+n)
    #        * prod_p (theta_(p+m+n))^cov_p
    #        * exp(eta_s).
    # The 'baseline HbA1c on Emax' covariate uses the OBSERVED per-subject
    # baseline (the HBA1C covariate column), not the latent baseline_hba1c_i
    # state above; in a typical FREM parameterisation the source paper
    # estimates the effect of the observed baseline on Emax.
    emax_i <- exp(lemax + etalemax) *
      (WT         / ref_wt     )^e_wt_emax *
      (CRCL       / ref_crcl   )^e_crcl_emax *
      (INSDOSE_BL / ref_insdose)^e_insdose_bl_emax *
      (HBA1C      / ref_hba1c  )^e_hba1c_emax *
      e_sexf_emax^SEXF *
      e_insdt_csii_emax^INSDT_CSII

    # ---- Individual placebo rate (Johnston 2019 Equation 1, fourth line
    # combined with Table 2 placebo-row covariate multipliers; the poster
    # equation lists Placebo = theta_d but Table 2 shows non-null SEXF and
    # INSDT_CSII multipliers on the placebo, encoded here as multiplicative
    # ratios). No IIV on placebo (not reported by the source). ----
    kplacebo_i <- exp(lkplacebo) *
      e_sexf_kplacebo^SEXF *
      e_insdt_csii_kplacebo^INSDT_CSII

    # ---- Placebo-arm HbA1c trajectory (canonical compartment hba1c_placebo).
    # State is initialised at the individual baseline HbA1c and drifts
    # linearly at kplacebo_i (%/h), i.e. encodes 'Baseline + Placebo * TIME'
    # from Johnston 2019 Equation 1. With AUC_EMPA = 0 this trajectory IS
    # the placebo-arm observation. ----
    d/dt(hba1c_placebo) <- kplacebo_i
    hba1c_placebo(0)    <- baseline_hba1c_i

    # ---- Direct-response Emax drug effect (Johnston 2019 Equation 1,
    # second line). AUC_EMPA is supplied as a static covariate from the
    # upstream popPK; the effect is instantaneous (no time delay). ----
    drug_effect <- emax_i * AUC_EMPA / (exp(lauc50) + AUC_EMPA)

    # ---- Drug-arm HbA1c trajectory (canonical hba1c_drug) = placebo-arm
    # trajectory minus the AUCss-driven Emax reduction. Combined
    # additive + proportional residual error per Johnston 2019 narrative. ----
    hba1c_drug <- hba1c_placebo - drug_effect
    hba1c_drug ~ add(addSd) + prop(propSd)
  })
}
