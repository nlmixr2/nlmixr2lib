Yin_2020_pexidartinib <- function() {
  description <- "Two-compartment population PK model for oral pexidartinib (CSF1R/KIT/FLT3 inhibitor) in healthy subjects and adult patients with tenosynovial giant cell tumour (TGCT) or other advanced solid tumours (Yin 2020). Absorption is sequential zero-order deposition into a depot (duration D1, lag time ALAG1) followed by first-order absorption (KA) into the central compartment, with linear elimination from central. Apparent clearance CL/F scales allometrically on (WT/80)^0.75 and is additionally modified by piecewise power effects of CRCL (active only when CRCL < 90 mL/min), AST (active only when AST > 80 U/L), and total bilirubin (active only when TBILI > 20.5 umol/L), plus multiplicative effects for Asian race (1.27x), healthy-participant cohort (1.26x; the Phase 1 healthy-subject studies), and female sex (0.869x). Apparent central and peripheral volumes Vc/F and Vp/F scale on (WT/80)^1; apparent inter-compartmental clearance Q/F scales on (WT/80)^0.75. Relative bioavailability of the Phase 1 formulation is fixed at 0.855 vs the Phase 3 / commercial reference formulation. Inter-individual variability is a 3x3 block on log(CL,Vc,Vp), independent diagonals on log(KA) and log(Q), and a Phase-1-formulation-specific IIV on the F1 bioavailability anchor. The published inter-occasion variability (5 occasions on KA, 10 occasions on F1) is not encoded structurally here (following the Andrews 2017 / Brooks 2021 tacrolimus precedent for the model-library use case where no operational occasion column is defined). Residual error is proportional with separate magnitudes for patient samples (29.7% CV) and healthy-subject samples (19.6% CV), switched per-subject by the DIS_HEALTHY indicator."
  reference <- "Yin O, Kang J, Knebel W, Zahir H, van de Sande M, Tap WD, Gelderblom H, Stacchiotti S, Greenberg J, Shuster D, Wagner AJ. Population Pharmacokinetic Analysis of Pexidartinib in Healthy Subjects and Patients With Tenosynovial Giant Cell Tumor or Other Solid Tumors. J Clin Pharmacol. 2021 Apr;61(4):480-492. doi:10.1002/jcph.1753. PDF on disk: ACoP 2019 poster of the same analysis (metrum_nd_pexidartinib_healthy.pdf)."
  vignette <- "Yin_2020_pexidartinib"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 80 kg (Yin 2020 Table 2 / Figure 2 caption: 'median value of 80 kg'). Allometric power scaling: (WT/80)^0.75 on CL/F and Q/F (fixed at theoretical); (WT/80)^1.0 on Vc/F and Vp/F (fixed at theoretical). Time-fixed per subject in the source analysis (baseline weight).",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; the reference subject in Yin 2020 Figure 2 caption is male)",
      notes              = "Multiplicative effect on CL/F: exp(theta14) = 0.869 when SEXF = 1 (female), 1.0 when SEXF = 0 (male). Yin 2020 Table 2 reports the female effect as 0.869 (95% CI 0.808-0.934).",
      source_name        = "Female"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator, 1 = Asian, 0 = non-Asian (the Yin 2020 reference race category)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian; the reference subject in Yin 2020 Figure 2 caption is non-Asian)",
      notes              = "Multiplicative effect on CL/F: exp(theta10) = 1.27 when RACE_ASIAN = 1 (Asian), 1.0 when RACE_ASIAN = 0. Yin 2020 Table 2 reports the Asian effect as 1.27 (95% CI 1.05-1.54). Asian subjects (N = 8) were a small subgroup of the pooled 375-subject analysis dataset; the effect was retained but Yin 2020 Results explicitly notes the wide 95% confidence interval.",
      source_name        = "Asian"
    ),
    CRCL = list(
      description        = "Baseline creatinine clearance (calculated)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 90 mL/min. Piecewise power effect on CL/F: factor = (CRCL/90)^(-0.0941) when CRCL < 90 mL/min, factor = 1 when CRCL >= 90 mL/min. Yin 2020 parameterises the renal-function effect as ( CRCL <90 / 90 )^theta9 (Table 2 row 3); the floor at CRCL = 90 means the effect is inactive at and above the reference. Yin 2020 Results explicitly notes the limited evaluation of CRCL effect because the analysis dataset contained a relatively narrow CRCL range and a small number of renally impaired subjects. Time-fixed per subject in the source analysis (baseline).",
      source_name        = "CRCL"
    ),
    AST = list(
      description        = "Baseline aspartate aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 80 U/L. Piecewise power effect on CL/F: factor = (AST/80)^0.0709 when AST > 80 U/L, factor = 1 when AST <= 80 U/L. Yin 2020 parameterises the hepatic-enzyme effect as ( AST >80 / 80 )^theta11 (Table 2 row 5); the floor at AST = 80 means the effect is inactive at and below the reference. Time-fixed per subject in the source analysis (baseline).",
      source_name        = "AST"
    ),
    TBILI = list(
      description        = "Baseline total bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 20.5 umol/L (Yin 2020 Table 2 row 6). Piecewise power effect on CL/F: factor = (TBILI/20.5)^0.244 when TBILI > 20.5 umol/L, factor = 1 when TBILI <= 20.5 umol/L. Yin 2020 parameterises the hyperbilirubinaemia effect as ( TBIL >20.5 / 20.5 )^theta12 (Table 2 row 6); the floor at TBILI = 20.5 means the effect is inactive at and below the reference. Source paper reports total bilirubin in SI umol/L (the canonical TBILI unit; no inline conversion needed). Time-fixed per subject (baseline).",
      source_name        = "TBIL"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator, 1 = healthy subject (Phase 1 clinical pharmacology study), 0 = patient with TGCT or other solid tumour",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient cohort; the reference subject in Yin 2020 Figure 2 caption is a patient)",
      notes              = "Time-fixed per subject. Yin 2020 codes this indicator as 'StHT' (study with healthy subjects). Multiplicative effect on CL/F: exp(theta13) = 1.26 when DIS_HEALTHY = 1 (healthy subject), 1.0 when DIS_HEALTHY = 0 (patient). Yin 2020 Table 2 reports the StHT effect as 1.26 (95% CI 1.16-1.36); Results states this corresponds to a 21% lower steady-state AUC0-24 in healthy subjects vs patients. The residual-error magnitude also switches with DIS_HEALTHY: 29.7% CV for patients vs 19.6% CV for healthy subjects (Table 2 Sigma rows). The N = 375 analysis dataset pooled 159 healthy subjects (seven Phase 1 clinical pharmacology studies U114-U121) with 216 patients (PLX108-01 and ENLIVEN).",
      source_name        = "StHT"
    ),
    FORM_PEX_PHASE1 = list(
      description        = "Pexidartinib formulation indicator, 1 = subject received the Phase 1 clinical formulation (used in PLX108-01 and U114), 0 = subject received the Phase 3 / commercial formulation (used in U116-U121 and ENLIVEN)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Phase 3 / commercial formulation; the typical-value bioavailability reference at F1 = 1.0)",
      notes              = "Per-dose-occasion indicator (a subject can in principle carry both values across studies, though within a single study the formulation is fixed). Fixed-effect multiplier on the depot bioavailability: F1_Phase1 = 0.855 (held fixed in Yin 2020 final model) when FORM_PEX_PHASE1 = 1, F1 = 1.0 when FORM_PEX_PHASE1 = 0. The Phase-1-formulation IIV (Omega 6.6 = 0.101, 32.6% CV) applies only when FORM_PEX_PHASE1 = 1 and has no effect when FORM_PEX_PHASE1 = 0 (Phase 3 / commercial). For typical simulations of the approved 800 mg/day TGCT regimen used in ENLIVEN, set FORM_PEX_PHASE1 = 0 for every subject. Yin 2020 Results: 'A formulation effect on F1 was fixed in the model to account for a 17% higher observed pexidartinib exposure with the phase 3 formulation compared with the phase 1 formulation' (1 / 0.855 = 1.17).",
      source_name        = "FORM"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 375L,
    n_observations  = 8430L,
    n_studies       = 9L,
    age_range       = "adults (>= 18 years); paper does not report a single pooled-cohort numeric range",
    weight_range    = "5th, 25th, 75th, 95th percentile and median 80 kg per Yin 2020 Figure 2 caption (53, ~, ~, ~, 80 kg); the explicit low percentile is 53 kg (~36% AUC increase vs the 80 kg reference)",
    weight_median   = "80 kg (pooled-cohort median; Yin 2020 Figure 2 caption)",
    sex_female_pct  = NA_real_,
    race_ethnicity  = c(Asian = round(8 / 375 * 100, 1), `non-Asian` = round(367 / 375 * 100, 1)),
    disease_state   = "Pooled cohort of healthy subjects and adult patients with tenosynovial giant cell tumour (TGCT) or other advanced solid tumours. 159 healthy subjects enrolled in seven Phase 1 clinical pharmacology studies (relative bioavailability, dose proportionality, drug-drug interaction with itraconazole / rifampin / esomeprazole, and food effect; doses 200-2400 mg single doses). 132 patients in Study PLX108-01 (Phase 1 dose-ranging in TGCT or other solid tumour; 200-1200 mg/day). 84 patients in Study PLX108-10 ENLIVEN (Phase 3 in TGCT; Part 1: 1000 mg/day for 2 weeks then 800 mg/day; Part 2: 800 mg/day).",
    dose_range      = "Single doses 200-2400 mg in healthy subjects; multiple doses 200-1200 mg/day in PLX108-01; 800 mg/day (400 mg BID) in ENLIVEN after the 2-week 1000 mg/day Part 1 lead-in.",
    regions         = "International (PLX108-01 and ENLIVEN multi-regional)",
    sampling_window = "Healthy-subject studies: serial PK samples up to 144 or 192 hours post-dose. PLX108-01: 5-6 PK samples per patient on Cycle 1 Day 1 and Cycle 1 Day 15, predose samples on Days 1, 8, and 16 of Cycle 1. ENLIVEN: 7 PK samples on Cycle 1 Day 15, random samples on Cycle 3 Day 1 and Cycle 5 Day 1.",
    assay           = "Validated liquid chromatography-tandem mass spectrometry method; LLOQ 2.5 ng/mL.",
    iov_structure   = "Yin 2020 Table 2 reports inter-occasion variability (IOV) on KA across 5 occasions (Omega 7.7 = 1.83, 229% CV) and on F1 across 10 occasions (Omega 12.12 = 0.0652, 25.9% CV). This model file does NOT encode the IOV structurally -- the source paper does not define an operational occasion column for the model-library use case, and the nlmixr2lib convention (Andrews 2017 / Brooks 2021 tacrolimus precedent) is to omit IOV when no occasion mapping is defined; see vignette Assumptions and deviations.",
    notes           = "Pooled population PK analysis of nine clinical studies (375 subjects, 8430 PK samples). Reference subject in Yin 2020 Figure 2 (covariate forest plot): male, non-Asian, patient cohort, WT 80 kg, CRCL >= 90 mL/min, AST <= 80 U/L, TBIL <= 20.5 umol/L."
  )

  ini({
    # Structural fixed effects -- Yin 2020 Table 2 final-model "Estimate" column.
    # Per the table footnote: "Estimates of theta modeled in the log domain were
    # exponentiated and are reported in the table" -- so the linear-scale
    # reference values are taken directly from Table 2 and re-logged here for the
    # `l` (log) prefix convention.
    lcl   <- log(5.83) ; label("Apparent oral clearance CL/F (L/hr) at the reference subject")              # Yin 2020 Table 2: CL/F exp(theta1) = 5.83 L/hr (95% CI 5.43-6.27)
    lvc   <- log(98.0) ; label("Apparent central volume Vc/F (L) at WT = 80 kg")                            # Yin 2020 Table 2: Vc/F exp(theta2) = 98.0 L (95% CI 90.0-107)
    lvp   <- log(116)  ; label("Apparent peripheral volume Vp/F (L) at WT = 80 kg")                         # Yin 2020 Table 2: Vp/F exp(theta3) = 116 L (95% CI 106-128)
    lq    <- log(20.7) ; label("Apparent inter-compartmental clearance Q/F (L/hr) at WT = 80 kg")           # Yin 2020 Table 2: Q/F exp(theta4) = 20.7 L/hr (95% CI 17.9-23.8)
    lka   <- log(6.82) ; label("First-order absorption rate constant KA (1/hr)")                            # Yin 2020 Table 2: KA exp(theta5) = 6.82 1/hr (95% CI 5.09-9.14)
    ltlag <- log(0.387); label("Absorption lag time ALAG1 (hr)")                                            # Yin 2020 Table 2: ALAG1 exp(theta6) = 0.387 hr (95% CI 0.385-0.390)
    ld1   <- log(1.22) ; label("Duration of zero-order deposition D1 (hr)")                                 # Yin 2020 Table 2: D1 exp(theta7) = 1.22 hr (95% CI 1.20-1.25)

    # Phase 1 formulation bioavailability anchor (Phase 3 = 1.0 reference).
    # Yin 2020 Results: "A formulation effect on F1 was fixed in the model to
    # account for a 17% higher observed pexidartinib exposure with the phase 3
    # formulation compared with the phase 1 formulation" (1 / 0.855 = 1.17).
    lfdepot <- fixed(log(0.855)); label("Log bioavailability of depot for Phase 1 formulation, relative to Phase 3 reference (unitless)") # Yin 2020 Table 2: F1 Phase1 exp(theta8) = 0.855 Fixed

    # Allometric exponents on body weight (Yin 2020 Table 2; theoretical
    # fixed exponents 0.75 on clearances and 1 on volumes, identical to the
    # standard size-scaling convention; reported in the table as the trailing
    # superscripts in the parameter equations).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT/80) on CL/F (unitless; fixed at theory value)")  # Yin 2020 Table 2 row 1 trailing exponent
    e_wt_vc <- fixed(1)   ; label("Allometric exponent of (WT/80) on Vc/F (unitless; fixed at theory value)")  # Yin 2020 Table 2 row 8 trailing exponent
    e_wt_vp <- fixed(1)   ; label("Allometric exponent of (WT/80) on Vp/F (unitless; fixed at theory value)")  # Yin 2020 Table 2 row 9 trailing exponent
    e_wt_q  <- fixed(0.75); label("Allometric exponent of (WT/80) on Q/F (unitless; fixed at theory value)")   # Yin 2020 Table 2 row 10 trailing exponent

    # Covariate effects on CL/F. The table reports raw power-exponent estimates
    # for the continuous-covariate piecewise effects (CRCL, AST, TBILI) and
    # already-exponentiated multiplicative estimates for the binary-covariate
    # effects (Asian, StHT/healthy, Female).
    e_crcl_cl   <- -0.0941; label("Piecewise power exponent of (CRCL/90) on CL/F, active only when CRCL < 90 mL/min")             # Yin 2020 Table 2: theta9 = -0.0941 (95% CI -0.402, 0.214)
    e_asian_cl  <-  1.27  ; label("Multiplicative effect on CL/F for Asian race (RACE_ASIAN = 1)")                                  # Yin 2020 Table 2: exp(theta10) = 1.27 (95% CI 1.05-1.54)
    e_ast_cl    <-  0.0709; label("Piecewise power exponent of (AST/80) on CL/F, active only when AST > 80 U/L")                    # Yin 2020 Table 2: theta11 = 0.0709 (95% CI -0.180, 0.322)
    e_tbili_cl  <-  0.244 ; label("Piecewise power exponent of (TBILI/20.5) on CL/F, active only when TBILI > 20.5 umol/L")         # Yin 2020 Table 2: theta12 = 0.244 (95% CI 0.183-0.306)
    e_healthy_cl<-  1.26  ; label("Multiplicative effect on CL/F for healthy-subject cohort (DIS_HEALTHY = 1)")                    # Yin 2020 Table 2: exp(theta13) = 1.26 (95% CI 1.16-1.36)
    e_female_cl <-  0.869 ; label("Multiplicative effect on CL/F for female sex (SEXF = 1)")                                       # Yin 2020 Table 2: exp(theta14) = 0.869 (95% CI 0.808-0.934)

    # Inter-individual variability. Yin 2020 Table 2 reports the Omega matrix
    # entries directly (variances on the log-eta scale; off-diagonal covariances
    # with derived correlations also tabulated). Variances verified via
    # %CV = sqrt(exp(omega^2) - 1) against the table's %CV column.
    #   Omega 1.1 CL/F = 0.0860 (%CV = 30.0)
    #   Omega 2.1 cov(Vc/F, CL/F) = 0.0774 (corr = 0.504)
    #   Omega 2.2 Vc/F = 0.274 (%CV = 56.1)
    #   Omega 3.1 cov(Vp/F, CL/F) = 0.0149 (corr = 0.110)
    #   Omega 3.2 cov(Vp/F, Vc/F) = -0.0467 (corr = -0.193)
    #   Omega 3.3 Vp/F = 0.213 (%CV = 48.8)
    #   Omega 4.4 Q/F = 0.406 (%CV = 70.8) [independent]
    #   Omega 5.5 KA = 1.31 (%CV = 165) [independent]
    #   Omega 6.6 Phase 1 F1 = 0.101 (%CV = 32.6) [independent; gated by FORM_PEX_PHASE1]
    etalcl + etalvc + etalvp ~ c(0.0860,
                                  0.0774, 0.274,
                                  0.0149, -0.0467, 0.213)
    etalq        ~ 0.406  # Yin 2020 Table 2 Omega 4.4 Q/F (independent diagonal)
    etalka       ~ 1.31   # Yin 2020 Table 2 Omega 5.5 KA (independent diagonal)
    etalfdepot   ~ 0.101  # Yin 2020 Table 2 Omega 6.6 Phase 1 F1 (independent diagonal; gated by FORM_PEX_PHASE1 in model())

    # Residual error -- proportional with population-specific magnitudes
    # switched per-subject by DIS_HEALTHY.
    #   Sigma 1.1 prop pat = 0.0883 (%CV = 29.7) -> propSd = sqrt(0.0883) = 0.297
    #   Sigma 2.2 prop ht  = 0.0384 (%CV = 19.6) -> propSd = sqrt(0.0384) = 0.196
    propSd_patient <- 0.297; label("Proportional residual SD for patient samples (fraction)")    # Yin 2020 Table 2: Sigma 1.1 prop pat = 0.0883 (%CV 29.7)
    propSd_healthy <- 0.196; label("Proportional residual SD for healthy-subject samples (fraction)")  # Yin 2020 Table 2: Sigma 2.2 prop ht = 0.0384 (%CV 19.6)
  })

  model({
    # Body-weight scaling reference: 80 kg pooled-cohort median (Yin 2020
    # Figure 2 caption).
    wt80 <- WT / 80

    # Piecewise covariate effects on CL/F. Yin 2020 uses the conditional
    # power-form notation (CRCL<90/90)^theta9, (AST>80/80)^theta11, and
    # (TBIL>20.5/20.5)^theta12 -- shorthand for: when the inequality holds,
    # use (value/reference)^theta; otherwise use 1. Encoded via the
    # Andrews_2017_tacrolimus precedent: raise the power-form factor to the
    # 0/1 indicator so that anything^0 = 1 selects the no-effect branch.
    crcl_below <- (CRCL < 90)
    f_crcl     <- ((CRCL / 90) ^ e_crcl_cl) ^ crcl_below

    ast_above  <- (AST > 80)
    f_ast      <- ((AST / 80) ^ e_ast_cl) ^ ast_above

    tbili_above <- (TBILI > 20.5)
    f_tbili     <- ((TBILI / 20.5) ^ e_tbili_cl) ^ tbili_above

    # Binary-covariate multiplicative effects on CL/F. Multiplier = e_<cov>_cl
    # when the indicator = 1, else 1.0. Encoded as 1 + (multiplier - 1) *
    # indicator (Andrews_2017_tacrolimus precedent).
    f_asian   <- 1 + (e_asian_cl   - 1) * RACE_ASIAN
    f_healthy <- 1 + (e_healthy_cl - 1) * DIS_HEALTHY
    f_female  <- 1 + (e_female_cl  - 1) * SEXF

    # Individual PK parameters with the Yin 2020 final-model covariate
    # equations (Table 2).
    cl   <- exp(lcl + etalcl) * wt80 ^ e_wt_cl * f_crcl * f_asian * f_ast * f_tbili * f_healthy * f_female
    vc   <- exp(lvc + etalvc) * wt80 ^ e_wt_vc
    vp   <- exp(lvp + etalvp) * wt80 ^ e_wt_vp
    q    <- exp(lq  + etalq)  * wt80 ^ e_wt_q
    ka   <- exp(lka + etalka)
    tlag <- exp(ltlag)
    d1   <- exp(ld1)

    # Bioavailability of the depot compartment.
    # - Phase 3 / commercial formulation (FORM_PEX_PHASE1 = 0): F1 = 1.0,
    #   with no IIV (etalfdepot has no effect because the eta term is
    #   multiplied by FORM_PEX_PHASE1).
    # - Phase 1 formulation (FORM_PEX_PHASE1 = 1): F1 = 0.855 * exp(etalfdepot)
    #   where etalfdepot ~ N(0, 0.101). The typical-value 0.855 anchor is
    #   held fixed (lfdepot <- fixed(log(0.855)) in ini()).
    fd_typ <- (1 - FORM_PEX_PHASE1) + FORM_PEX_PHASE1 * exp(lfdepot)
    fd     <- fd_typ * exp(etalfdepot * FORM_PEX_PHASE1)

    # Two-compartment linear elimination from central.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Sequential zero-order deposition (duration d1, lag time tlag) into the
    # depot followed by first-order absorption (ka) into central. Bioavailability
    # fd multiplies the dose entering depot.
    alag(depot) <- tlag
    dur(depot)  <- d1
    f(depot)    <- fd

    # Plasma pexidartinib concentration: amount in mg, volume in L -> mg/L;
    # multiply by 1000 to report ng/mL to align with the bioanalytical LLOQ
    # (2.5 ng/mL) and clinical convention.
    Cc <- central / vc * 1000

    # Per-subject population-conditional proportional residual error switched
    # by DIS_HEALTHY.
    propSd <- propSd_patient * (1 - DIS_HEALTHY) + propSd_healthy * DIS_HEALTHY
    Cc ~ prop(propSd)
  })
}
