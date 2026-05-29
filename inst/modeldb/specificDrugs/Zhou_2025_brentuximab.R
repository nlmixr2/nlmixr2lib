Zhou_2025_brentuximab <- function() {
  description <- "Coupled population PK model for brentuximab vedotin antibody-drug conjugate (ADC) and its released payload monomethyl auristatin E (MMAE) in pediatric patients (5-18 years) with relapsed/refractory or newly diagnosed Hodgkin lymphoma or systemic anaplastic large-cell lymphoma (Zhou 2025). ADC is described by a linear 3-compartment model with first-order elimination; MMAE by a 2-compartment model with first-order elimination. ADC -> MMAE flux is the sum of (a) a one-time saturable target-binding flux Kd*Target*ADC (initial Target = 1 unitless, irreversibly depleted) and (b) a proteolytic flux FM*exp(-ALFM*tad)*K10*ADC where the conversion fraction declines as a function of time after the most recent dose. Both fluxes accumulate in an intermediate Lag compartment that empties to MMAE central with rate Klag. Final-model parameter values come from Zhou 2025 supplementary Tables S1 (ADC) and S2 (MMAE); equations come from the NONMEM control streams in Zhou 2025 Supplementary Methods."
  reference <- "Zhou X, Mould DR, Gore L, Bai X, Gupta N. Optimizing Brentuximab Vedotin Dosing in Pediatric Patients with Advanced Hodgkin Lymphoma: A Population Pharmacokinetic and Exposure-Response Analysis. Clin Pharmacol Ther. 2025;117(6):1803-1810. doi:10.1002/cpt.3629. PMID 40095373."
  vignette <- "Zhou_2025_brentuximab"
  paper_specific_compartments <- c("lag")

  units <- list(
    time          = "hour",
    dosing        = "umol",
    concentration = "umol/L"
  )

  covariateData <- list(
    BSA = list(
      description        = "Baseline body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power effects on ADC CL (exponent 1.38), ADC V3 = vp2 (exponent 1.96), MMAE CL (exponent 0.772), and MMAE central volume VM (exponent 0.546). Reference 1.8 m^2 (NONMEM normalization NBSA = BSA / 1.8 from Zhou 2025 supplement control stream). Pediatric study median BSA 1.52 m^2; range 0.79-2.03 m^2.",
      source_name        = "BSA"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power effects on ADC CL (exponent -0.776), MMAE CL (exponent -0.0805), and MMAE Kd binding rate (exponent -4.11). Reference 40 g/L (NONMEM normalization NALB = ALB / 40 from Zhou 2025 supplement control stream).",
      source_name        = "ALB"
    ),
    CREAT = list(
      description        = "Baseline serum creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power effect on MMAE CL only (exponent -0.0952). Reference 45.689 umol/L (NONMEM normalization NCREAT = CREAT / 45.689 from Zhou 2025 supplement control stream). Bilirubin column in Zhou 2025 Table 1 reports umol/L; the supplement control stream uses the same 'umol/L' convention for CREAT.",
      source_name        = "CREAT"
    ),
    TUMSZ = list(
      description        = "Baseline tumor linear diameter (sum of linear diameters of target lesions)",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline. Power effect on ADC CL only (exponent 0.12). Reference 41 mm (NONMEM normalization NLDIAM = LDIAM / 41 from Zhou 2025 supplement control stream). Distinct from the 'sum of tumor area, mm^2' column in Zhou 2025 Table 1 (median 1581-1639 mm^2 for the sum-of-areas convention).",
      source_name        = "LDIAM"
    ),
    TUMTP_CHL = list(
      description        = "Classical Hodgkin lymphoma tumor-type indicator (1 = HL, 0 = non-HL)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 1,
      notes              = "Time-fixed. Used via the derived NONHL = 1 - TUMTP_CHL indicator in the model: multiplicative effect 0.509 on ADC Q2 for non-HL patients (~50% reduction relative to HL); multiplicative effect 0.296 on MMAE central volume VM for non-HL; multiplicative effect 0.884 on the ADC->MMAE conversion-decay rate ALFM for non-HL. Reference category here is HL (TUMTP_CHL = 1) because that is the more common group in this pediatric cohort and the paper anchors typical-value parameters to HL patients. Non-HL in Zhou 2025 = systemic anaplastic large-cell lymphoma (sALCL).",
      source_name        = "DIS"
    ),
    ADA_POS = list(
      description        = "Anti-drug antibody positivity indicator (ever positive)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed once positive (Zhou 2025 supplement: 'ATAPOS and ADA are filled as once positive then positive for the rest of the study'). Multiplicative effect 2.6 on ADC CL (~2.6-fold higher CL for ADA-positive subjects) and 0.696 on MMAE CL (~30% lower MMAE CL for ADA-positive subjects). Pediatric ADA incidence was low in these two studies, so the effect estimate is informed by a small subgroup.",
      source_name        = "ATAPOS"
    ),
    CONMED_AVD = list(
      description        = "Brentuximab vedotin + AVD chemotherapy combination indicator (1 = A+AVD, 0 = single-agent BV)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Time-fixed per study (Study 1 = single-agent BV, CONMED_AVD = 0; Study 2 = BV + AVD, CONMED_AVD = 1). Multiplicative effect 2.12 on ADC CL — patients on the A+AVD regimen have ~2.1-fold higher ADC CL relative to single-agent BV. The Zhou 2025 NONMEM dataset encodes this as DOX (doxorubicin / adriamycin administration flag); doxorubicin is given on the same days as vinblastine and dacarbazine in the AVD backbone, so the DOX column is equivalent to the A+AVD regimen indicator.",
      source_name        = "DOX"
    )
  )

  population <- list(
    n_subjects     = 95L,
    n_studies      = 2L,
    age_range      = "5-18 years (study 1: 7-18; study 2: 5 to <18)",
    age_median     = "14 years (both studies)",
    weight_range   = "18.8-87.0 kg",
    weight_median  = "49 kg (study 2) / 49.9 kg (study 1)",
    bsa_range      = "0.79-2.03 m^2",
    bsa_median     = "1.52 m^2",
    sex_female_pct = 41.1,
    race_ethnicity = c(White = 68.4, Black = 12.6, Asian = 5.3, Other = 13.7),
    disease_state  = "Pediatric (study 1) relapsed/refractory systemic anaplastic large-cell lymphoma (sALCL) or Hodgkin lymphoma (HL); (study 2) advanced-stage CD30+ newly diagnosed classical HL.",
    dose_range     = "Study 1 (C25002): single-agent BV 1.4-1.8 mg/kg IV every 3 weeks. Study 2 (C25004): BV 48 mg/m^2 IV every 2 weeks combined with doxorubicin 25 mg/m^2 + vinblastine 6 mg/m^2 + dacarbazine 375 mg/m^2.",
    regions        = "United States (children's hospital network).",
    study_phase    = "Two open-label phase I/II studies (NCT01492088 single-agent dose escalation; NCT02979522 BSA-based BV + AVD).",
    n_observations = "9479 ADC + MMAE concentration records (2608 from study 1, 6871 from study 2).",
    reference_subject = "BSA 1.8 m^2, ALB 40 g/L, CREAT 45.689 umol/L, LDIAM 41 mm, TUMTP_CHL 1 (HL), ADA_POS 0, CONMED_AVD 0.",
    notes          = "Baseline characteristics from Zhou 2025 Table 1 (per-study); pooled n = 95 patients. Albumin range 23-51 g/L (median 39-40); creatinine clearance Cockcroft 85.5-301.4 mL/min (median 132.3-165.3); sum-of-tumor-area median 1581-1639 mm^2. The NONMEM reference values 1.8 m^2 (BSA) and 40 g/L (ALB) are adult/normal-range anchors rather than pediatric medians."
  )

  ini({
    # ADC structural parameters (Zhou 2025 Table S1; 3-compartment linear with
    # first-order elimination from central; ADVAN11 TRANS4 in NONMEM).
    lcl  <- log(0.0208); label("ADC clearance (CL, L/hr)")                                     # Zhou 2025 Table S1: 0.0208 (16.3% RSE)
    lvc  <- log(2.54);   label("ADC central volume (V1, L)")                                   # Zhou 2025 Table S1: 2.54 (0.7% RSE)
    lq   <- log(0.0192); label("ADC inter-compartmental clearance to peripheral 1 (Q2, L/hr)") # Zhou 2025 Table S1: 0.0192 (14.2% RSE)
    lvp  <- log(97.1);   label("ADC peripheral volume 1 (V2, L)")                              # Zhou 2025 Table S1: 97.1 (18.9% RSE)
    lq2  <- log(0.0865); label("ADC inter-compartmental clearance to peripheral 2 (Q3, L/hr)") # Zhou 2025 Table S1: 0.0865 (7.3% RSE)
    lvp2 <- log(3.39);   label("ADC peripheral volume 2 (V3, L)")                              # Zhou 2025 Table S1: 3.39 (15.8% RSE)

    # ADC covariate effects (Zhou 2025 Table S1; reference values from the
    # NONMEM control stream NBSA = BSA/1.8, NALB = ALB/40, NLDIAM = LDIAM/41).
    e_bsa_cl    <- 1.38;   label("Power exponent of (BSA/1.8) on ADC CL (unitless)")               # Zhou 2025 Table S1: 1.38 (23.2% RSE)
    e_alb_cl    <- -0.776; label("Power exponent of (ALB/40) on ADC CL (unitless)")                # Zhou 2025 Table S1: -0.776 (12.3% RSE)
    e_tumsz_cl  <- 0.12;   label("Power exponent of (TUMSZ/41) on ADC CL (unitless)")              # Zhou 2025 Table S1: 0.12 (17.3% RSE)
    e_nonhl_q2  <- 0.509;  label("Power-form multiplier of non-HL on ADC Q2: q2 *= e_nonhl_q2^(1 - TUMTP_CHL)") # Zhou 2025 Table S1: 0.509 (32.2% RSE)
    e_bsa_vp2   <- 1.96;   label("Power exponent of (BSA/1.8) on ADC V3 = vp2 (unitless)")         # Zhou 2025 Table S1: 1.96 (20.8% RSE)
    e_ada_cl    <- 2.6;    label("Power-form multiplier of ADA positivity on ADC CL: cl *= e_ada_cl^ADA_POS") # Zhou 2025 Table S1: 2.6 (6.4% RSE)
    e_avd_cl    <- 2.12;   label("Power-form multiplier of A+AVD coadministration on ADC CL: cl *= e_avd_cl^CONMED_AVD") # Zhou 2025 Table S1: 2.12 (18.7% RSE) — supplement table label of 'theta13' for this row is a typo (control stream confirms theta14)

    # MMAE structural parameters (Zhou 2025 Table S2; 2-compartment linear with
    # an upstream Target binding pool and Lag compartment fed by ADC, ADVAN13
    # custom ODE in NONMEM).
    lcl_mmae   <- log(0.794);    label("MMAE clearance (CLM, L/hr)")                                    # Zhou 2025 Table S2: 0.794 (2.0% RSE)
    lvc_mmae   <- log(20.1);     label("MMAE central volume (VM, L)")                                   # Zhou 2025 Table S2: 20.1 (0.7% RSE)
    lq_mmae    <- log(0.628);    label("MMAE inter-compartmental clearance (QM, L/hr)")                 # Zhou 2025 Table S2: 0.628 (1.2% RSE)
    lvp_mmae   <- log(2.74);     label("MMAE peripheral volume (VMP, L)")                               # Zhou 2025 Table S2: 2.74 (0.7% RSE)
    lkd_mmae   <- log(0.0186);   label("MMAE binding rate constant (Kd, 1/hr)")                         # Zhou 2025 Table S2: 0.0186 (0.4% RSE)
    lalfm_mmae <- log(0.00462);  label("Decay rate of ADC->MMAE proteolytic-conversion fraction (ALFM, 1/hr)") # Zhou 2025 Table S2: 0.00462 (0.7% RSE)
    lklag_mmae <- log(60.8);     label("Lag-compartment empty rate constant (Klag, 1/hr)")              # Zhou 2025 Table S2: 60.8 (1.3% RSE)
    # FM (fraction metabolized) is fixed to 1 in Zhou 2025 Table S2 — encoded
    # as a literal constant in model() rather than an estimated parameter.

    # MMAE covariate effects (Zhou 2025 Table S2; reference values from the
    # NONMEM control stream NCREAT = CREAT/45.689, NALB = ALB/40, NBSA = BSA/1.8).
    e_creat_cl_mmae    <- -0.0952; label("Power exponent of (CREAT/45.689) on MMAE CL (unitless)")          # Zhou 2025 Table S2: -0.0952 (9.3% RSE)
    e_alb_cl_mmae      <- -0.0805; label("Power exponent of (ALB/40) on MMAE CL (unitless)")                # Zhou 2025 Table S2: -0.0805 (27.5% RSE)
    e_bsa_cl_mmae      <- 0.772;   label("Power exponent of (BSA/1.8) on MMAE CL (unitless)")               # Zhou 2025 Table S2: 0.772 (8.4% RSE)
    e_ada_cl_mmae      <- 0.696;   label("Power-form multiplier of ADA positivity on MMAE CL: clm *= e_ada_cl_mmae^ADA_POS") # Zhou 2025 Table S2: 0.696 (5.3% RSE)
    e_nonhl_vc_mmae    <- 0.296;   label("Power-form multiplier of non-HL on MMAE VM: vcm *= e_nonhl_vc_mmae^(1 - TUMTP_CHL)") # Zhou 2025 Table S2: 0.296 (15.1% RSE)
    e_bsa_vc_mmae      <- 0.546;   label("Power exponent of (BSA/1.8) on MMAE VM (unitless)")               # Zhou 2025 Table S2: 0.546 (9.7% RSE)
    e_nonhl_alfm_mmae  <- 0.884;   label("Power-form multiplier of non-HL on ALFM: alfm *= e_nonhl_alfm_mmae^(1 - TUMTP_CHL)") # Zhou 2025 Table S2: 0.884 (17.7% RSE)
    e_alb_kd_mmae      <- -4.11;   label("Power exponent of (ALB/40) on MMAE Kd (unitless)")                # Zhou 2025 Table S2: -4.11 (3.3% RSE)

    # IIV (log-normal). %CV from Zhou 2025 Tables S1 and S2; converted via
    # omega^2 = log(CV^2 + 1). Off-diagonal correlations are not reported in
    # the supplement, so the IIV block is diagonal here (the NONMEM control
    # stream estimates a 3x3 OMEGA BLOCK on ADC CL/Q2/V3 but only the
    # diagonal CV% is reported in Table S1 — see the vignette's Assumptions
    # and deviations section).
    etalcl  ~ 0.21        # Zhou 2025 Table S1: CL 48.4% CV; omega^2 = log(0.484^2 + 1) = 0.211
    etalq   ~ 0.353       # Zhou 2025 Table S1: Q2 65.0% CV; omega^2 = log(0.65^2 + 1)  = 0.353
    etalvp2 ~ 0.31        # Zhou 2025 Table S1: V3 60.3% CV; omega^2 = log(0.603^2 + 1) = 0.310

    etalcl_mmae   ~ 0.228      # Zhou 2025 Table S2: CLM 50.6% CV; omega^2 = log(0.506^2 + 1) = 0.228
    etalvc_mmae   ~ 0.494      # Zhou 2025 Table S2: VM 79.9% CV;  omega^2 = log(0.799^2 + 1) = 0.494
    etalkd_mmae   ~ 1.151      # Zhou 2025 Table S2: Kd 147% CV;   omega^2 = log(1.47^2 + 1)  = 1.151
    etalalfm_mmae ~ 0.582      # Zhou 2025 Table S2: ALFM 88.8% CV; omega^2 = log(0.888^2 + 1) = 0.582

    # Residual variability (LTBS in NONMEM = proportional in linear space for
    # nlmixr2; %CV from Tables S1 and S2 is the proportional SD on the linear
    # concentration scale).
    propSd       <- 0.321; label("Proportional residual error on ADC Cc (fraction)")    # Zhou 2025 Table S1: 32.1% CV (0.7% RSE)
    propSd_mmae    <- 0.375; label("Proportional residual error on MMAE Cc_mmae (fraction)") # Zhou 2025 Table S2: 37.5% CV (1.2% RSE)
  })

  model({
    # 1. Derived covariate terms. Reference values: BSA 1.8 m^2, ALB 40 g/L,
    # CREAT 45.689 umol/L, TUMSZ 41 mm (from Zhou 2025 supplement NONMEM
    # control stream). NONHL = 1 - TUMTP_CHL puts HL patients in the
    # reference (effect = 1) and non-HL (sALCL) patients on the multiplier.
    nbsa   <- BSA / 1.8
    nalb   <- ALB / 40
    ncreat <- CREAT / 45.689
    ntumsz <- TUMSZ / 41
    nonhl  <- 1 - TUMTP_CHL

    # 2. Individual ADC parameters
    cl_adc <- exp(lcl + etalcl) *
              nbsa^e_bsa_cl *
              nalb^e_alb_cl *
              ntumsz^e_tumsz_cl *
              e_ada_cl^ADA_POS *
              e_avd_cl^CONMED_AVD
    v1_adc <- exp(lvc)
    q2_adc <- exp(lq + etalq) * e_nonhl_q2^nonhl
    v2_adc <- exp(lvp)
    q3_adc <- exp(lq2)
    v3_adc <- exp(lvp2 + etalvp2) * nbsa^e_bsa_vp2

    # 2b. Individual MMAE parameters
    cl_mmae <- exp(lcl_mmae + etalcl_mmae) *
               ncreat^e_creat_cl_mmae *
               nalb^e_alb_cl_mmae *
               nbsa^e_bsa_cl_mmae *
               e_ada_cl_mmae^ADA_POS
    vc_mmae <- exp(lvc_mmae + etalvc_mmae) * e_nonhl_vc_mmae^nonhl * nbsa^e_bsa_vc_mmae
    qm      <- exp(lq_mmae)
    vp_mmae <- exp(lvp_mmae)
    kd      <- exp(lkd_mmae + etalkd_mmae) * nalb^e_alb_kd_mmae
    alfm    <- exp(lalfm_mmae + etalalfm_mmae) * e_nonhl_alfm_mmae^nonhl
    klag    <- exp(lklag_mmae)

    # 3. Micro-constants
    k10 <- cl_adc  / v1_adc
    k12 <- q2_adc  / v1_adc
    k21 <- q2_adc  / v2_adc
    k13 <- q3_adc  / v1_adc
    k31 <- q3_adc  / v3_adc
    k40 <- cl_mmae / vc_mmae
    k47 <- qm      / vc_mmae
    k74 <- qm      / vp_mmae

    # 4. Time-varying ADC->MMAE proteolytic-conversion fraction (Zhou 2025
    # supplement $DES: FMT = FM * exp(-ALFM * (T - DTIME))). FM is fixed at
    # 1 in Zhou 2025 Table S2. tad() returns time after the most recent dose,
    # which mirrors NONMEM's DTIME bookkeeping (resets at each dose).
    fmt <- exp(-alfm * tad())

    # 5. ODE system (Zhou 2025 supplement $DES, mapped from NONMEM
    # compartments to nlmixr2 names: A(1)=central, A(2)=peripheral1,
    # A(3)=peripheral2, A(4)=central_mmae, A(5)=target, A(6)=lag,
    # A(7)=peripheral1_mmae).
    target(0) <- 1
    d/dt(central)          <- -(k10 + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1)      <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2)      <-  k13 * central - k31 * peripheral2
    d/dt(target)           <- -kd * target * central
    d/dt(lag)              <-  kd * target * central + fmt * k10 * central - klag * lag
    d/dt(central_mmae)     <-  klag * lag - (k40 + k47) * central_mmae + k74 * peripheral1_mmae
    d/dt(peripheral1_mmae) <-  k47 * central_mmae - k74 * peripheral1_mmae

    # 6. Observations and residual error. The Zhou 2025 NONMEM dataset uses
    # molar units (AMT in umol, DV in umol/L = uM) — see the supplement
    # control stream comment "AMT IN UM ; DVUM IS DV IN UM". Doses provided
    # in mg must be converted before dosing this model: umol = mg / MW_kDa
    # (MW_BV approx 153.4 kDa for the ADC; MW_MMAE approx 0.718 kDa for the
    # payload). The binding equation kd * target * central is the only
    # term whose absolute parameter value depends on the amount units, so
    # using mass units here without rescaling kd would corrupt the
    # ADC -> MMAE conversion magnitude.
    Cc      <- central      / v1_adc
    Cc_mmae <- central_mmae / vc_mmae

    Cc      ~ prop(propSd)
    Cc_mmae ~ prop(propSd_mmae)
  })
}
