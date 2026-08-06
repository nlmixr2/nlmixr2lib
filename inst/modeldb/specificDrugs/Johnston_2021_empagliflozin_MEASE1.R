Johnston_2021_empagliflozin_MEASE1 <- function() {
  description <- paste(
    "Semi-mechanistic direct-effect exposure-response (PD-only) model for",
    "empagliflozin on total daily insulin dose (TDID), mean daily glucose",
    "(MDG), and glycated hemoglobin (HbA1c) in adults with type 1 diabetes",
    "(T1D) (M-EASE-1; Johnston 2021). Individual steady-state empagliflozin",
    "AUC (AUC_EMPA, supplied as a per-subject covariate column from the",
    "upstream Johnston 2021 popPK; see Johnston_2021_empagliflozin_popPK)",
    "drives (i) a proportional-Emax reduction in TDID, (ii) a direct-Emax",
    "reduction in MDG offset by a linear time-dependent placebo drift and",
    "coupled to the TDID reduction via a power-form insulin effect, and",
    "(iii) HbA1c as a power-form function of MDG relative to baseline.",
    "Covariates: sex, body weight, eGFR, insulin delivery type (CSII vs",
    "MDI), and baseline HbA1c."
  )
  reference <- paste(
    "Johnston CK, Eudy-Byrne RJ, Elmokadem A, Nock V, Marquard J,",
    "Soleymanlou N, Riggs MM, Liesenfeld K-H.",
    "A Model-Informed Drug Development (MIDD) Approach for a Low Dose",
    "of Empagliflozin in Patients with Type 1 Diabetes.",
    "Pharmaceutics. 2021;13(4):485. doi:10.3390/pharmaceutics13040485."
  )
  vignette <- "Johnston_2021_empagliflozin"
  units <- list(
    time          = "h",
    dosing        = "n/a (no drug dosing events; empagliflozin exposure enters as the per-subject AUC_EMPA covariate from the upstream popPK)",
    concentration = "n/a (multi-output PD-only model; the observations are HbA1c in % NGSP, MDG in mg*day/dL, and TDID in IU/kg)",
    AUC_EMPA      = "nmol*h/L"
  )
  covariateData <- list(
    AUC_EMPA = list(
      description        = "Steady-state empagliflozin AUC over the q24h dosing interval supplied as a per-subject (time-fixed) drug-exposure covariate from the upstream Johnston 2021 popPK.",
      units              = "nmol*h/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-subject (time-fixed) steady-state AUC of empagliflozin. The",
        "source authors generate AUCss from individual empirical Bayes",
        "estimates of the upstream popPK (Johnston 2021 Table S2; see",
        "Johnston_2021_empagliflozin_popPK) and pass them to the M-EASE-1",
        "PD model as a static covariate. AUC_EMPA = 0 recovers the placebo",
        "arm (the empagliflozin-Emax terms vanish). Reported AUC50 values",
        "in the M-EASE-1 fit: 110 nmol*h/L (TDID) and 370 nmol*h/L (MDG)",
        "(Table S3 / Table 3)."
      ),
      source_name        = "AUCss"
    ),
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effects centred at 82 kg on baseline TDID, baseline HbA1c, and Emax_MDG (Table S3 / equations page 13).",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate (BSA-normalised)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effects centred at 99 mL/min/1.73 m^2 on baseline TDID and Emax_MDG (Table S3 / equations page 13). Renamed from the source column eGFR to the canonical CRCL per covariate-columns.md.",
      source_name        = "eGFR"
    ),
    HBA1C = list(
      description        = "Baseline (pre-treatment) HbA1c (per-subject, time-fixed)",
      units              = "% (NGSP)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline HbA1c as a per-subject covariate on baseline TDID (Table S3 TDID equation; reference 8.1 %). NOT the same as the time-course HbA1c observations: this is the per-subject baseline anchor used inside the TDID equation; observations are the model-predicted longitudinal HbA1c trajectory.",
      source_name        = "Base.HbA1c"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male reference)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; the reference patient is male per Table S3 M-EASE-1 header)",
      notes              = "Multiplicative effects on baseline TDID, baseline HbA1c, and Emax_MDG (Table S3 / equations page 13). Encoded as `<multiplier>^SEXF`.",
      source_name        = "Sex (1 = female; 0 = male)"
    ),
    INSDT_CSII = list(
      description        = "Insulin delivery type indicator (1 = continuous subcutaneous insulin infusion (CSII), 0 = multiple daily injections (MDI) reference)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (MDI; the reference patient is on MDI per Johnston 2021 population description)",
      notes              = "Multiplicative effect on Emax_MDG (Table S3 / equations page 13). Follows the `<COLUMN>_<LEVEL>` decomposition pattern used elsewhere in the registry.",
      source_name        = "INSDT"
    )
  )

  population <- list(
    species             = "human",
    n_subjects          = 796L,
    n_subjects_empa     = 534L,
    n_subjects_placebo  = 262L,
    n_hba1c_observations = 4824L,
    n_tdid_observations  = 189182L,
    n_mdg_observations   = 4243L,
    n_studies           = 2L,
    studies             = "EASE-1 (phase II, 4-week) and EASE-2 (phase III, 52-week). External evaluation on EASE-3 (phase III, 26-week, out-of-sample).",
    age_range           = "21-69 years (95th-percentile interval at baseline)",
    weight_range        = "55-125 kg (95th-percentile interval at baseline)",
    egfr_range          = "57-127 mL/min/1.73 m^2 (95th-percentile interval at baseline)",
    hba1c_range         = "7.2-9.5 % (95th-percentile interval at baseline)",
    disease_state       = paste0(
      "Adults with type 1 diabetes mellitus (T1D) on background insulin ",
      "therapy. Reference patient (Table S3 M-EASE-1 header): male, eGFR ",
      "99 mL/min/1.73 m^2, body weight 82 kg, baseline HbA1c 8.1 %, ",
      "cumulative MDG over 24 h = 4266 mg*day/dL."
    ),
    dose_range          = "Empagliflozin 0 (placebo), 2.5, 10, 25 mg QD (as the empagliflozin arms in EASE-1 and EASE-2); AUC_EMPA covariate supplies the exposure.",
    regions             = "Multi-national (EASE-1 / EASE-2 trials).",
    notes               = paste0(
      "Semi-mechanistic stepwise-fit model developed in NONMEM 7.4 with ",
      "FOCE-I (Section 2.2). Step 1 fits TDID as a proportional-Emax ",
      "function of AUC_EMPA; Step 2 fits MDG using the derived TDID and a ",
      "placebo-only-fit MDG_t0, AUC50_MDG, and PBO_MDG (fixed); Step 3 ",
      "fits HbA1c as a power-form function of MDG relative to baseline. ",
      "This file encodes the final integrated model; the fixed-from-placebo ",
      "MDG parameters (MDG_t0, AUC50_MDG, PBO_MDG) and the fixed-from-",
      "placebo TDID AUC50 are all wrapped in `fixed()` per the source's ",
      "sequential-fit strategy (Table S3 footnotes a, c)."
    )
  )

  ini({
    # ---- TDID model (Table S3 M-EASE-1 TDID block; reference patient at ----
    # ---- baseline TDID = 0.657 IU/kg for male, WT 82 kg, eGFR 99, HbA1c 8.1) ----
    ltdid_t0        <- log(0.657); label("Typical baseline TDID at the reference patient (IU/kg; log-scale)")  # Table S3 theta_1 = 0.657 (95% CI 0.617, 0.7)
    e_wt_tdid_t0    <-  0.317;     label("Power exponent of (WT/82) on baseline TDID")                          # Table S3 theta_2 = 0.317
    sexf_tdid_t0    <-  0.96;      label("Multiplicative effect of female sex on baseline TDID (ref = male)")  # Table S3 theta_3 = 0.96
    e_egfr_tdid_t0  <-  0.145;     label("Power exponent of (CRCL/99) on baseline TDID")                        # Table S3 theta_4 = 0.145
    e_hba1c_tdid_t0 <-  0.368;     label("Power exponent of (HBA1C/8.1) on baseline TDID")                      # Table S3 theta_5 = 0.368

    emax_tdid       <- 0.186;                label("Typical Emax_TDID -- maximal proportional reduction in TDID (fraction)")  # Table S3 theta_11 / Table 3 = 0.186 (95% CI 0.145, 0.238)
    lauc50_tdid     <- fixed(log(110));      label("AUC50 for TDID (nmol*h/L; log-scale, from placebo-only fit)")       # Table S3 theta_12 = 110 nmol*h/L (fixed; footnote a "Estimated from placebo only data and fixed")

    # ---- MDG model (Table S3 M-EASE-1 MDG block; reference patient MDG_t0 = 4160 mg*day/dL) ----
    lmdg_t0         <- log(4160);            label("Typical baseline MDG (mg*day/dL cumulative over 24 h; log-scale)")          # Table S3 theta_1 = 4160
    ins_mdg_eff     <- fixed(-0.261);        label("Power exponent for (TDID/TDID_t0) on MDG time course (unitless)")   # Table S3 theta_3 = -0.261 (INS_MDG effect; fixed per footnote c)
    pbo_mdg_val     <- fixed(0.0136);        label("Placebo drift rate on MDG (mg*day/dL per hour)")                     # Table S3 theta_3 = 0.0136 (PBO_MDG; fixed per footnote c)
    lauc50_mdg      <- fixed(log(370));      label("AUC50 for MDG (nmol*h/L; log-scale, from placebo-only fit)")          # Table S3 theta_4 = 370 nmol*h/L (fixed per footnote c)
    lemax_mdg       <- log(634);             label("Typical Emax_MDG -- maximal MDG reduction (mg*day/dL; log-scale)")          # Table S3 theta_5 = 634 mg*day/dL (95% CI 534, 753)
    e_wt_emax_mdg   <- -0.113;               label("Power exponent of (WT/82) on Emax_MDG")                                     # Table S3 theta_6 = -0.113
    sexf_emax_mdg   <-  1.09;                label("Multiplicative effect of female sex on Emax_MDG (ref = male)")             # Table S3 theta_7 = 1.09
    e_egfr_emax_mdg <-  0.0707;              label("Power exponent of (CRCL/99) on Emax_MDG")                                   # Table S3 theta_8 = 0.0707
    insdt_emax_mdg  <-  0.995;               label("Multiplicative effect of CSII (vs MDI) on Emax_MDG")                        # Table S3 theta_9 = 0.995

    # ---- HbA1c model (Table S3 M-EASE-1 HbA1c block; reference patient at ----
    # ---- HbA1c_t0 = 8.15 % for male, WT 82 kg) ----
    lhba1c_t0       <- log(8.15); label("Typical baseline HbA1c at the reference patient (%; log-scale)")  # Table S3 theta_1 = 8.15 (95% CI 8.09, 8.21)
    sexf_hba1c_t0   <- 0.99;      label("Multiplicative effect of female sex on baseline HbA1c (ref = male)")  # Table S3 theta_2 = 0.99
    e_wt_hba1c_t0   <- -0.0258;   label("Power exponent of (WT/82) on baseline HbA1c")                          # Table S3 theta_3 = -0.0258
    gamma_mdgeff    <- 0.487;     label("Typical MDGEFF -- power exponent on (MDG/MDG_t0) driving HbA1c")      # Table S3 theta_4 = 0.487 (95% CI 0.445, 0.532)

    # ---- Inter-individual variability (Table S3 M-EASE-1) ----
    # TDID block: 2x2 correlated (baseline TDID + Emax_TDID)
    #   omega^2_TDIDBASE = 0.0974 (32.0% CV)
    #   Cov(TDIDBASE, TDIDEmax) = 0.00579 (rho = 0.0215)
    #   omega^2_TDIDEmax = 0.554 (86.0% CV)
    # Emax_TDID IIV is encoded here as log-normal on the fraction directly
    # (etaemax_tdid multiplies emax_tdid). The source uses a logit-scale
    # parameterisation with an additional EASE-2 study-specific offset; see
    # the Errata in the vignette Assumptions and deviations section.
    etaltdid_t0 + etaemax_tdid ~ c(0.0974,
                                   0.00579, 0.554)

    # MDG block: two independent etas
    etalmdg_t0   ~ 0.009   # omega^2_MDGt0 = 0.009 (9.51% CV)
    etalemax_mdg ~ 0.0744  # omega^2_MDGEmax = 0.0744 (27.8% CV)

    # HbA1c block: 2x2 correlated (baseline HbA1c + MDGEFF exponent).
    # Source parameterisation: gamma_i = gamma * exp(etagamma_mdgeff) so
    # etagamma_mdgeff multiplies the reported gamma value log-normally.
    etalhba1c_t0 + etagamma_mdgeff ~ c(0.00437,
                                       0.0106,  0.461)

    # ---- Residual error (Table S3; per-endpoint proportional + additive) ----
    # TDID:   sigma^2_prop = 0.0239 -> CV = 15.6% ; sigma^2_add = 0.001 -> SD = 0.0316 IU/kg
    # MDG:    sigma^2_prop = 0.0254 -> CV = 16.0% ; sigma^2_add = 0.001 -> SD = 0.0316 mg*day/dL
    # HbA1c:  sigma^2_prop = 0.00218 -> CV = 4.67% (proportional only reported)
    propSd_tdid  <- 0.156;  label("Proportional residual error on TDID (fraction; 15.6% CV, sqrt(exp(0.0239)-1))")  # Table S3
    addSd_tdid   <- 0.0316; label("Additive residual error on TDID (IU/kg; SD = sqrt(0.001))")                       # Table S3
    propSd_mdg   <- 0.160;  label("Proportional residual error on MDG (fraction; 16.0% CV, sqrt(exp(0.0254)-1))")   # Table S3
    addSd_mdg    <- 0.0316; label("Additive residual error on MDG (mg*day/dL; SD = sqrt(0.001) -- essentially negligible relative to typical MDG ~4000)")  # Table S3
    propSd_hba1c <- 0.0467; label("Proportional residual error on HbA1c (fraction; 4.67% CV, sqrt(exp(0.00218)-1))")  # Table S3
  })

  model({
    # ---- Reference values (Table S3 M-EASE-1 header + equations page 13) ----
    ref_wt    <- 82   # kg
    ref_egfr  <- 99   # mL/min/1.73 m^2
    ref_hba1c <- 8.1  # % (NGSP; the /8.1 normalisation in the TDID equation)

    # ---- Step 1: Individual baseline TDID (Table S3 M-EASE-1 TDID equation) ----
    # TDID_t0,i = exp(theta_1) * theta_3^SEXF * (WT/82)^theta_2
    #             * (eGFR/99)^theta_4 * (HbA1c/8.1)^theta_5 * exp(eta_1)
    tdid_t0_i <- exp(ltdid_t0 + etaltdid_t0) *
                 sexf_tdid_t0^SEXF *
                 (WT / ref_wt)^e_wt_tdid_t0 *
                 (CRCL / ref_egfr)^e_egfr_tdid_t0 *
                 (HBA1C / ref_hba1c)^e_hba1c_tdid_t0

    # Emax_TDID individual (log-normal IIV on the reported fraction).
    emax_tdid_i <- emax_tdid * exp(etaemax_tdid)

    # TDID(t) -- proportional Emax reduction driven by AUC_EMPA. Steady-state
    # w.r.t. time because AUC_EMPA is time-invariant per subject; the EASE-1
    # first-week INC multiplier is omitted from this extraction (see Errata).
    auc50_tdid <- exp(lauc50_tdid)
    tdid <- tdid_t0_i * (1 - emax_tdid_i * AUC_EMPA / (auc50_tdid + AUC_EMPA))

    # ---- Step 2: Individual baseline MDG and Emax_MDG (Table S3 MDG block) ----
    mdg_t0_i    <- exp(lmdg_t0 + etalmdg_t0)

    emax_mdg_i  <- exp(lemax_mdg + etalemax_mdg) *
                   sexf_emax_mdg^SEXF *
                   (WT / ref_wt)^e_wt_emax_mdg *
                   (CRCL / ref_egfr)^e_egfr_emax_mdg *
                   insdt_emax_mdg^INSDT_CSII

    auc50_mdg <- exp(lauc50_mdg)

    # MDG(t) -- Equation 2:
    #   MDG = MDG_t0 * (TDID/TDID_t0)^INS_MDG_EFF + PBO_MDG*TIME
    #         - Emax_MDG * AUC / (AUC50_MDG + AUC)
    # `time` is model time in hours; the placebo term accumulates linearly.
    mdg <- mdg_t0_i * (tdid / tdid_t0_i)^ins_mdg_eff +
           pbo_mdg_val * time -
           emax_mdg_i * AUC_EMPA / (auc50_mdg + AUC_EMPA)

    # ---- Step 3: Individual baseline HbA1c and gamma (Table S3 HbA1c block) ----
    hba1c_t0_i <- exp(lhba1c_t0 + etalhba1c_t0) *
                  sexf_hba1c_t0^SEXF *
                  (WT / ref_wt)^e_wt_hba1c_t0

    # gamma_i = gamma * exp(eta_gamma)  (Equation 3: (MDG/MDG_t0)^(theta_4*exp(eta_2)))
    gamma_i <- gamma_mdgeff * exp(etagamma_mdgeff)

    # HbA1c(t) -- Equation 3
    hba1c <- hba1c_t0_i * (mdg / mdg_t0_i)^gamma_i

    # ---- Multi-output residual error ----
    tdid  ~ prop(propSd_tdid)  + add(addSd_tdid)
    mdg   ~ prop(propSd_mdg)   + add(addSd_mdg)
    hba1c ~ prop(propSd_hba1c)
  })
}
