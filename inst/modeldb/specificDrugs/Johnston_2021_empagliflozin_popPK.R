Johnston_2021_empagliflozin_popPK <- function() {
  description <- paste(
    "Population pharmacokinetic (PopPK) model for empagliflozin in adults",
    "with type 1 diabetes (T1D) enrolled in the EASE clinical program",
    "(pooled EASE-1 phase II, EASE-2 phase III, and EASE-3 phase III).",
    "Two-compartment model with first-order elimination and a sequential",
    "zero-order (duration D1) plus first-order (ka) oral absorption with",
    "a depot lag time (ALAG). Full covariate power / multiplicative model",
    "on CL/F (age via TPRO/AP/eGFR/TDID/WT + sex + smoking), V2/F, V3/F,",
    "Q/F, and ka (Johnston 2021 Table S2). The PopPK was adapted from a",
    "prior T2D empagliflozin analysis (Baron 2016) and re-fit on the",
    "EASE-1 + EASE-3 dataset with EASE-2 as external evaluation, then",
    "re-estimated on the combined EASE-1/-2/-3 data to generate the",
    "individual AUCss inputs to the M-EASE-1 and M-EASE-2 exposure-response",
    "models. This file encodes the final full-covariate model."
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
    dosing        = "mg empagliflozin (oral, once daily)",
    concentration = "Cc in nmol/L (converted from mg/L via MW 450.91 g/mol)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "empagliflozin popPK", units = NA_character_, specimen = "administration site", verified = FALSE),
    central     = list(analyte = "empagliflozin popPK", units = NA_character_, specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "empagliflozin popPK", units = NA_character_, specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effects centred at 44 years on V2/F, V3/F, and ka (Table S2 / equations page 10).",
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effects centred at 70 kg on CL/F, V2/F, Q/F, and V3/F (Table S2).",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate (BSA-normalised)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-form effect centred at 99 mL/min/1.73 m^2 on CL/F (Table S2 theta_25).",
        "Renamed from the source column eGFR to the canonical CRCL per",
        "covariate-columns.md. The paper does not specify which eGFR equation",
        "(MDRD vs CKD-EPI) was used."
      ),
      source_name        = "eGFR"
    ),
    TPRO = list(
      description        = "Total serum protein at baseline",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effects centred at 68 g/L on CL/F, V2/F, and V3/F (Table S2). Reference value read from the equations on page 10 (`TPRO_i(g/L)/68(g/L)`).",
      source_name        = "TPRO"
    ),
    ALP = list(
      description        = "Alkaline phosphatase at baseline",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect centred at 73 U/L on CL/F (Table S2 theta_24). Reference value read from the equation `AP_i(IU/L)/73(IU/L)` on page 10 (the Figure S1 caption prints `AP = 73 IU/kg` but the model equation uses IU/L; treated as a caption typo). Renamed from the source column AP to the canonical ALP per covariate-columns.md.",
      source_name        = "AP"
    ),
    INSDOSE_BL = list(
      description        = "Baseline total daily insulin dose normalised to body weight",
      units              = "U/kg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect centred at 0.6 IU/kg on CL/F (Table S2 theta_26). Reference value read from the equation `TDID_i(IU/kg)/0.6(IU/kg)` on page 10.",
      source_name        = "TDID"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male reference)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Multiplicative effects on CL/F, V2/F, V3/F, and ka (Table S2 theta_8 / theta_9 / theta_10 / theta_11). Encoded as `<multiplier>^SEXF` so SEXF=0 recovers the male baseline.",
      source_name        = "Sex (1 = female; 0 = male)"
    ),
    SMOKE_CURRENT = list(
      description        = "Current-smoker indicator (paired with SMOKE_NEVER)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = former smoker (when paired with SMOKE_NEVER = 0)",
      notes              = paste(
        "Johnston 2021 Table S2 reports smoking as a 3-level categorical with",
        "NON-smoker (=never smoker) as the reference (CL multipliers",
        "theta_12 = 1.02 for ex-smoker, theta_13 = 1.08 for current smoker).",
        "The canonical SMOKE_NEVER + SMOKE_CURRENT pairing uses FORMER smoker",
        "as the reference (both indicators = 0); the model recodes the",
        "published coefficients so the reference becomes former smoker",
        "(same recoding pattern used in Baron_2016_empagliflozin)."
      ),
      source_name        = "Cur-Smoker (1 = current; 0 otherwise)"
    ),
    SMOKE_NEVER = list(
      description        = "Never-smoker indicator (paired with SMOKE_CURRENT)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = former smoker (when paired with SMOKE_CURRENT = 0)",
      notes              = "See SMOKE_CURRENT notes for the recoding from Johnston 2021's never-as-reference encoding to the canonical former-as-reference pairing.",
      source_name        = "Ex-Smoker (1 = ex/former; 0 otherwise)"
    )
  )

  population <- list(
    species             = "human",
    n_subjects          = 1241L,
    n_subjects_male     = 614L,
    n_subjects_female   = 627L,
    n_studies           = 3L,
    studies             = "EASE-1 (phase II, 4-week, 75 patients with T1D, once-daily empagliflozin 2.5 / 10 / 25 mg or placebo), EASE-2 (phase III, 52-week, 721 patients with T1D, empagliflozin 10 / 25 mg or placebo), EASE-3 (phase III, 26-week, 948 patients with T1D, empagliflozin 2.5 / 10 / 25 mg or placebo). Table 1.",
    age_range           = "21-69 years (95th-percentile interval across studies at baseline)",
    weight_range        = "52-126 kg (95th-percentile interval across studies at baseline)",
    egfr_range          = "55-129 mL/min/1.73 m^2 (95th-percentile interval across studies at baseline)",
    hba1c_range         = "7.2-9.6 % (95th-percentile interval across studies at baseline)",
    insulin_range       = "0.36-1.34 IU/kg total daily insulin dose (95th-percentile interval across studies at baseline)",
    sex_female_pct      = 100 * 627 / 1241,
    disease_state       = paste0(
      "Adults with type 1 diabetes mellitus (T1D) on background insulin ",
      "therapy. Reference subject for the full covariate model: male, ",
      "nonsmoker (never smoker), total daily insulin dose 0.6 IU/kg, ",
      "alkaline phosphatase 73 IU/L, total protein 68 g/L, eGFR 99 ",
      "mL/min/1.73 m^2, weight 70 kg, age 44 years (Figure S1 caption + ",
      "equations page 10)."
    ),
    dose_range          = "Empagliflozin 0 (placebo), 2.5, 10, 25 mg once daily (oral).",
    regions             = "Multi-national (EASE-1 / -2 / -3 trials).",
    notes               = paste0(
      "PopPK was estimated with NONMEM 7.3 using FOCE-I (Section 2.2). The ",
      "model was adapted from the prior T2D empagliflozin PopPK of Baron ",
      "2016; initial development used EASE-1 and EASE-3 with EASE-2 as an ",
      "external evaluation dataset, then the best-performing model was ",
      "re-estimated on the combined dataset to produce individual AUCss ",
      "estimates for the exposure-response analyses. No covariates were ",
      "found to have clinically relevant impacts on exposure (Section 3.1; ",
      "Figure S1 forest plot)."
    )
  )

  ini({
    # ---- Structural PK parameters (Johnston 2021 Table S2; reference subject ----
    # ---- as summarised in Figure S1 caption: male, nonsmoker, TDID 0.6 IU/kg, ----
    # ---- AP 73 IU/L, TPRO 68 g/L, eGFR 99, WT 70 kg, age 44) ----
    # Table S2 theta_1 = 11.2 L/h is CL/F for the paper's reference subject, who is a
    # NEVER smoker. The smoking multipliers below are recoded onto the canonical
    # FORMER-smoker reference, so the base value must be rebased onto a former smoker
    # too (x theta_12 = 1.02). Without that, cl_smoke_never = 1/1.02 is applied to a
    # base that already represents a never smoker and every never-smoker prediction
    # comes out 2% low (11.2 / 1.02 = 10.98 L/h). With the rebasing, all three strata
    # match the paper: never = 11.424 / 1.02 = 11.2 (Table S2 theta_1);
    # former = 11.424 = 11.2 x 1.02; current = 11.424 x 1.08 / 1.02 = 12.096 = 11.2 x 1.08.
    lcl   <- log(11.2 * 1.02);  label("Apparent oral clearance CL/F (L/h), rebased to the canonical former-smoker reference")  # Table S2 theta_1 = 11.2 L/h (never smoker) x theta_12 = 1.02 (ex-smoker vs never)
    lvc   <- log(1.69);  label("Apparent central volume V2/F (L) for the reference subject (renamed V2 -> vc per canonical convention)")  # Table S2 theta_2 = 1.69 L
    lq    <- log(6.14);  label("Apparent inter-compartmental clearance Q/F (L/h) for the reference subject")  # Table S2 theta_3 = 6.14 L/h
    lvp   <- log(82.2);  label("Apparent peripheral volume V3/F (L) for the reference subject (renamed V3 -> vp per canonical convention)") # Table S2 theta_4 = 82.2 L
    lka   <- log(0.233); label("First-order absorption rate constant ka (1/h) for the reference subject")     # Table S2 theta_5 = 0.233 1/h
    ld1   <- log(0.623); label("Duration of zero-order input into the depot D1 (h)")                          # Table S2 theta_6 = 0.623 h
    lalag <- log(0.135); label("Absorption lag time on the depot compartment ALAG (h)")                       # Table S2 theta_7 = 0.135 h

    # ---- Categorical (multiplicative) covariate effects (Table S2) ----
    cl_female         <- 0.892;  label("Multiplicative effect of female sex on CL/F (theta_8; ref = male)")           # Table S2 theta_8 = 0.892
    vc_female         <- 0.986;  label("Multiplicative effect of female sex on V2/F (theta_9; ref = male)")           # Table S2 theta_9 = 0.986
    vp_female         <- 0.762;  label("Multiplicative effect of female sex on V3/F (theta_10; ref = male)")          # Table S2 theta_10 = 0.762
    ka_female         <- 1.05;   label("Multiplicative effect of female sex on ka (theta_11; ref = male)")            # Table S2 theta_11 = 1.05
    # Smoking (recoded from Johnston 2021's never-as-reference to canonical
    # SMOKE_CURRENT + SMOKE_NEVER pair with former-as-reference; same recoding
    # pattern as Baron_2016_empagliflozin).
    cl_smoke_never    <- 1 / 1.02;  label("Multiplicative effect of never-smoker (canonical vs former ref) on CL/F = 1 / theta_12_paper")      # Recoded from Table S2 theta_12 = 1.02 (ex-smoker vs never)
    cl_smoke_current  <- 1.08 / 1.02; label("Multiplicative effect of current-smoker (canonical vs former ref) on CL/F = theta_13_paper / theta_12_paper")  # Recoded from Table S2 theta_13 = 1.08 (current smoker vs never) and theta_12 = 1.02

    # ---- Continuous (power-form) covariate effects on CL/F (Table S2) ----
    e_wt_cl        <-  0.394;  label("Power exponent of (WT/70) on CL/F")            # Table S2 theta_17 = 0.394
    e_tpro_cl      <- -0.245;  label("Power exponent of (TPRO/68) on CL/F")          # Table S2 theta_21 = -0.245
    e_ap_cl        <- -0.0541; label("Power exponent of (ALP/73) on CL/F")           # Table S2 theta_24 = -0.0541
    e_egfr_cl      <-  0.271;  label("Power exponent of (CRCL/99) on CL/F")          # Table S2 theta_25 = 0.271
    e_tdid_cl      <-  0.0469; label("Power exponent of (INSDOSE_BL/0.6) on CL/F")   # Table S2 theta_26 = 0.0469

    # ---- Continuous (power-form) covariate effects on V2/F (Table S2) ----
    e_age_vc       <- -1.54;   label("Power exponent of (AGE/44) on V2/F (now vc)")   # Table S2 theta_14 = -1.54
    e_tpro_vc      <- -4.27;   label("Power exponent of (TPRO/68) on V2/F (now vc)")  # Table S2 theta_22 = -4.27
    e_wt_vc        <-  2.57;   label("Power exponent of (WT/70) on V2/F (now vc)")    # Table S2 theta_18 = 2.57

    # ---- Continuous (power-form) covariate effects on V3/F (Table S2) ----
    e_age_vp       <-  0.190;  label("Power exponent of (AGE/44) on V3/F (now vp)")   # Table S2 theta_15 = 0.190
    e_tpro_vp      <- -0.381;  label("Power exponent of (TPRO/68) on V3/F (now vp)")  # Table S2 theta_23 = -0.381
    e_wt_vp        <-  0.414;  label("Power exponent of (WT/70) on V3/F (now vp)")    # Table S2 theta_20 = 0.414

    # ---- Continuous (power-form) covariate effect on Q/F (Table S2) ----
    e_wt_q         <-  1.11;   label("Power exponent of (WT/70) on Q/F")              # Table S2 theta_19 = 1.11

    # ---- Continuous (power-form) covariate effect on ka (Table S2) ----
    e_age_ka       <-  0.0419; label("Power exponent of (AGE/44) on ka")              # Table S2 theta_16 = 0.0419

    # ---- PK inter-individual variability (Table S2; log-normal on-diagonal blocks) ----
    # CL/F, Q/F, and V3/F share a full 3x3 covariance block (Cov: CL-Q, CL-V3, Q-V3).
    # ka is independent. No IIV was estimated on V2/F, D1, or ALAG.
    etalcl + etalq + etalvp ~ c( 0.0644,
                                -0.0614,  0.0952,
                                 0.0467, -0.0806,  0.190)
    # Cov(CL,Q)  = -0.0614 (rho = -0.784); Cov(CL,V3) = 0.0467 (rho = 0.422);
    # Cov(Q,V3)  = -0.0806 (rho = -0.599).
    # Diagonals: omega^2_CL = 0.0644 (25.8% CV), omega^2_Q = 0.0952 (31.6% CV),
    # omega^2_V3 = 0.190 (45.7% CV).
    etalka ~ 0.0258            # omega^2_ka = 0.0258 (16.2% CV)

    # ---- Residual error (Table S2; proportional, log-normal convention) ----
    # Table S2 reports two study-specific proportional residuals:
    #   Studies-3 (EASE-3): omega^2 = 0.128 (37.0% CV) -- larger dataset, primary
    #   Studies-1 (EASE-1): omega^2 = 0.0796 (28.8% CV)
    # This file encodes the EASE-3 residual (37.0% CV, sqrt(0.128) = 0.358 on
    # linear scale) as the primary; the EASE-1 alternative is documented in the
    # vignette Assumptions and deviations section. No EASE-2 residual was
    # reported (EASE-2 was originally the external-eval dataset; the final
    # combined re-estimation re-used the EASE-3 residual).
    propSd <- 0.370; label("Proportional residual error on Cc (fraction; 37.0% CV, sqrt(omega^2 = 0.128) from log-normal convention)")  # Table S2 sigma^2_EASE3 = 0.128 -> CV = sqrt(exp(0.128)-1) = 37.0%
  })

  model({
    # ---- Empagliflozin molecular weight (BI 10773; CAS 864070-44-0; C23H27ClO7 -> 450.91 g/mol) ----
    mw_empa     <- 450.91
    nmol_per_mg <- 1e6 / mw_empa

    # ---- Reference values for centred covariate effects (Table S2 equations page 10) ----
    ref_age  <- 44    # years
    ref_wt   <- 70    # kg
    ref_egfr <- 99    # mL/min/1.73 m^2
    ref_tpro <- 68    # g/L
    ref_ap   <- 73    # U/L
    ref_tdid <- 0.6   # IU/kg

    # ---- Individual PK parameters (Table S2 equations page 10) ----
    cl <- exp(lcl + etalcl) *
          cl_female^SEXF *
          cl_smoke_never^SMOKE_NEVER *
          cl_smoke_current^SMOKE_CURRENT *
          (WT / ref_wt)^e_wt_cl *
          (TPRO / ref_tpro)^e_tpro_cl *
          (ALP / ref_ap)^e_ap_cl *
          (CRCL / ref_egfr)^e_egfr_cl *
          (INSDOSE_BL / ref_tdid)^e_tdid_cl

    vc <- exp(lvc) *
          vc_female^SEXF *
          (AGE / ref_age)^e_age_vc *
          (TPRO / ref_tpro)^e_tpro_vc *
          (WT / ref_wt)^e_wt_vc

    vp <- exp(lvp + etalvp) *
          vp_female^SEXF *
          (AGE / ref_age)^e_age_vp *
          (TPRO / ref_tpro)^e_tpro_vp *
          (WT / ref_wt)^e_wt_vp

    q  <- exp(lq + etalq) *
          (WT / ref_wt)^e_wt_q

    ka <- exp(lka + etalka) *
          ka_female^SEXF *
          (AGE / ref_age)^e_age_ka

    d1        <- exp(ld1)
    alag_dpt  <- exp(lalag)

    # ---- Micro-rate constants (two-compartment linear) ----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- Sequential zero-order + first-order oral absorption with lag ----
    # Dose enters the depot as a zero-order infusion of duration D1 (delivered
    # via `rate = -2` on the dose event so rxode2 uses the model-specified
    # `dur(depot)`) starting after the lag time ALAG. The depot then drains
    # via first-order ka into central.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    dur(depot)  <- d1
    alag(depot) <- alag_dpt

    # ---- Observation: empagliflozin plasma concentration in nmol/L ----
    # central is in mg (from oral dosing), vc in L -> mg/L; multiply by 1e6/MW
    # to convert to nmol/L (matching the paper's AUC unit nmol*h/L).
    Cc <- central / vc * nmol_per_mg
    Cc ~ prop(propSd)
  })
}
