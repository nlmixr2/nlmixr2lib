DeJongh_2014_antithrombin <- function() {
  description <- paste(
    "Two-compartment population PK model with linear elimination from",
    "the central compartment for recombinant human antithrombin (rhAT)",
    "in adults with hereditary antithrombin deficiency (HD) undergoing",
    "perioperative or peripartum treatment (DeJongh 2014). Delivery",
    "(pregnancy) status approximately doubles CL and Vss versus",
    "non-pregnant surgery patients. Baseline AT activity is modelled",
    "as a subject-level endogenous parameter (Dansirikul-type baseline",
    "handling) that adds to the drug-derived AT activity in the",
    "observation."
  )
  reference <- paste(
    "DeJongh J, Frieling J, Lowry S, Drenth H-J.",
    "Pharmacokinetics of recombinant human antithrombin in delivery",
    "and surgery patients with hereditary antithrombin deficiency.",
    "Clin Appl Thromb Hemost. 2014;20(4):355-364.",
    "doi:10.1177/1076029613516188.",
    sep = " "
  )
  vignette <- "DeJongh_2014_antithrombin"
  units <- list(
    time          = "h",
    dosing        = "IU",
    concentration = "% of normal"
  )

  covariateData <- list(
    PREG = list(
      description        = "Pregnancy status (delivery-vs-surgery cohort indicator)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-pregnant surgery-cohort patient; canonical reference)",
      notes              = paste(
        "1 = delivery patient (peripartum pregnant woman treated for",
        "peripartum VTE prevention); 0 = surgery patient (non-pregnant",
        "adult receiving perioperative treatment). The DeJongh 2014",
        "cohort has no pregnant non-delivery women, so PREG functions",
        "as the delivery-vs-surgery cohort indicator; the source paper",
        "calls the covariate 'delivery'. Time-fixed per subject.",
        "Discussion attributes the CL and Vss increase in delivery",
        "patients to pregnancy physiology (intra- and extra-vascular",
        "volume expansion, altered hepatic clearance) rather than to",
        "the delivery event itself."
      ),
      source_name        = "DELIVERY"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 47L,
    n_studies      = 3L,
    age_range      = "adults (specific range not reported by the paper)",
    weight_range   = "not reported (simulations used typical WT = 76 kg with variance 204 kg^2)",
    sex_female_pct = 87,
    disease_state  = paste(
      "Hereditary antithrombin (AT) deficiency (HD).",
      "Development: 15 non-high-risk HD patients receiving a single",
      "IV bolus dose of rhAT (AT III 009-00; 50 or 100 IU/kg;",
      "female:male 13:2). External validation: 32 HD patients in",
      "high-risk situations (21 delivery peripartum, 11 surgery",
      "perioperative) across AT III 01002 (NCT00056550) and",
      "AT HD 012-04 (NCT00110513). Mean baseline AT activity",
      "46-54% of normal across cohorts."
    ),
    dose_range     = paste(
      "Single-dose PK: 50 or 100 IU/kg IV bolus (15-min infusion).",
      "Clinical trials: individualized loading (15-min infusion)",
      "followed immediately by continuous infusion maintenance;",
      "median total dose 523-1194 IU/kg (surgery cohorts) and",
      "726-973 IU/kg (delivery cohorts) over 3-19 days."
    ),
    regions        = "European Union (development PK study AT III 009-00)",
    notes          = paste(
      "Baseline demographics from paper Methods 'Studies and",
      "Patients' and Table 1. NONMEM versions V.2 and 7.2 (Icon)",
      "with FOCE-I. Proportional residual error model. Population",
      "PK model developed on AT III 009-00 (n=15), externally",
      "validated on AT III 01002 (n=14; delivery covariate",
      "identified) and AT HD 012-04 (n=18)."
    )
  )

  ini({
    # Structural parameters -- surgery-cohort reference (PREG = 0). The delivery
    # cohort (PREG = 1) is a log-additive shift on CL and Vss; see below.
    # DeJongh 2014 Table 3. The paper parameterises volume as Vss (steady-state
    # total volume) plus Vr = Vss/Vcentral - 1 (peripheral-to-central ratio) to
    # improve NONMEM numerical stability under IIV on Vss. Vr is estimated as a
    # single typical value common to both cohorts (no IIV, no PREG effect), so
    # the reparameterisation collapses back to canonical (Vc, Vp) with a fixed
    # Vp/Vc ratio: Vc = Vss/(1 + Vr) and Vp = Vss - Vc = Vss * Vr/(1 + Vr).
    lcl        <- log(0.665)                    ; label("Clearance for surgery-cohort reference (L/h)")                    # Table 3: CL_surgery = 0.665 L/h (SE 0.0493)
    lvc        <- log(7.72 / (1 + 1.51))        ; label("Central volume of distribution for surgery reference (L)")        # Table 3: Vss_surgery = 7.72 L (SE 1.26) and Vr = 1.51 (SE 0.331); Vc = Vss/(1+Vr) = 3.0757 L
    lvp        <- log(7.72 * 1.51 / (1 + 1.51)) ; label("Peripheral volume of distribution for surgery reference (L)")     # Table 3: Vp = Vss - Vc = Vss * Vr/(1+Vr) = 4.6443 L (surgery)
    lq         <- log(0.613)                    ; label("Intercompartmental clearance (L/h; shared between cohorts)")      # Table 3: Q = 0.613 L/h (SE 0.646; large SE consistent with the paper's note that Q was reparameterised as a multiple of CL for numerical stability)
    lrbase     <- log(44.7)                     ; label("Baseline AT activity typical value (% of normal)")                # Table 3: ATBL = 44.7% (SE 3.03)

    # Delivery-cohort covariate effects. Table 3 reports CL_delivery = 1.38 L/h
    # and Vss_delivery = 14.3 L, giving log-additive shifts of log(1.38/0.665)
    # and log(14.3/7.72) respectively. Because Vr is common to both cohorts,
    # the Vss shift propagates identically to Vc and Vp; the canonical shared-
    # exponent naming convention e_<cov>_vc_vp captures this.
    e_preg_cl    <- log(1.38 / 0.665)           ; label("PREG log-additive effect on CL (unitless)")                       # Table 3: CL_delivery / CL_surgery = 1.38 / 0.665 = 2.075
    e_preg_vc_vp <- log(14.3 / 7.72)            ; label("Shared PREG log-additive effect on Vc and Vp (unitless)")         # Table 3: Vss_delivery / Vss_surgery = 14.3 / 7.72 = 1.852

    # IIV: variances (omega^2) reported directly in Table 3. CL and Vss share
    # a 2 x 2 covariance block; ATBL has an independent variance. The single
    # Vss random effect propagates uniformly to Vc and Vp because Vr has no
    # IIV, so etalvc is reused on both vc and vp inside model().
    etalcl + etalvc ~ c(0.0676,
                        0.0395, 0.0521)                                                                                    # Table 3: omega^2_CL = 0.0676 (SE 0.0205), cov(CL,Vss) = 0.0395 (SE 0.0175), omega^2_Vss = 0.0521 (SE 0.026)
    etalrbase       ~ 0.0519                                                                                               # Table 3: omega^2_ATBL = 0.0519 (SE 0.0221)

    # Residual error: proportional model. Table 3 reports sigma^2 = 0.0289
    # (SE 0.00543), so propSd = sqrt(0.0289) = 0.17.
    propSd     <- 0.17                          ; label("Proportional residual error (fraction of predicted AT activity)")  # Table 3: sigma^2 = 0.0289 (SE 0.00543) -> propSd = sqrt(0.0289) = 0.17
  })

  model({
    # 1. Individual PK parameters. etalvc represents the single Vss random effect
    #    that propagates equally to Vc and Vp (Vr is a fixed typical value).
    cl <- exp(lcl + e_preg_cl    * PREG + etalcl)
    vc <- exp(lvc + e_preg_vc_vp * PREG + etalvc)
    vp <- exp(lvp + e_preg_vc_vp * PREG + etalvc)
    q  <- exp(lq)

    # 2. Subject-level baseline AT activity (Dansirikul-type endogenous baseline
    #    parameter). This adds to the drug-derived AT activity in the observation.
    rbase <- exp(lrbase + etalrbase)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. Two-compartment IV model. `central` and `peripheral1` hold antithrombin
    #    amount in IU. Loading is administered as a 15-min bolus infusion into
    #    `central`; maintenance is continuous infusion into `central`.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Observation: AT activity (% of normal) = subject baseline + drug-derived
    #    activity. Unit conversion factor 0.1 %/(IU/L) applies the WHO potency
    #    definition of the antithrombin International Unit: 1 IU = amount of AT
    #    in 1 mL of pooled normal human plasma, so 1 IU/mL = 100% activity and
    #    1 IU/L = 0.1% activity. Vc is in L and `central` in IU; central/vc is
    #    IU/L; multiplied by 0.1 gives % of normal. The paper does not state
    #    this conversion factor explicitly but it is required for the paper's
    #    L-valued Vss to reproduce the paper's %-valued observations and the
    #    published incremental-recovery point estimates (Table 2:
    #    IR ~ 2 %/(IU/kg) matches Vc / body-weight and this conversion factor).
    Cc <- rbase + (central / vc) * 0.1
    Cc ~ prop(propSd)
  })
}
