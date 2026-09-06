Ahn_2023_vupanorsen <- function() {
  description <- "Population PK and PK/PD model for vupanorsen (PF-07285557), a GalNAc3-conjugated 2'-O-methoxyethyl antisense oligonucleotide targeting ANGPTL3 mRNA, pooled across two phase I and two phase II studies (Ahn 2023). Two-compartment disposition with first-order subcutaneous absorption; allometric body weight on all disposition parameters (88 kg reference), Asian race on CL/F and Vc/F, female sex and anti-drug-antibody positivity on CL/F, and a 160 mg dose-level effect on Q/F. Three simultaneously fitted indirect-response endpoints (ANGPTL3, triglycerides, non-HDL-cholesterol) in which the predicted peripheral-compartment concentration inhibits the zero-order production rate, with study-population factors on baseline and on potency."
  reference <- "Ahn JE, Terra SG, Liu J. A population pharmacokinetic and pharmacokinetic-pharmacodynamic analysis of vupanorsen from phase I and phase II studies. CPT Pharmacometrics Syst Pharmacol. 2023;12(7):988-1000. doi:10.1002/psp4.12969"
  vignette <- "Ahn_2023_vupanorsen"

  # ANGPTL3, triglyceride and non-HDL-cholesterol turnover pools. The register
  # already carries `ldl` and `hc24` as canonical lipid / sterol PD outputs, and
  # `tg` appears as a paper-specific state in `Ousey_2026_plozasiran.R`; all
  # three are declared paper-specific here pending canonicalisation.
  paper_specific_compartments <- c("angptl3", "tg", "nonhdlc")

  # The magnitude of the vupanorsen residual error differs between the phase I
  # and the phase II studies, so Table 2 reports two combined additive +
  # proportional error models rather than one (Methods "Base PK model"). The
  # applicable pair is selected inside `model()` from the study indicators.
  paper_specific_residual_sds <- c(
    "propSdPhase1", "addSdPhase1", "propSdPhase2", "addSdPhase2"
  )

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling with FIXED exponents 0.75 on CL/F and Q/F and 1 on Vc/F and Vp/F, normalised to an 88 kg reference (Ahn 2023 Equations 11-14; Methods 'Base PK model' states the allometry constants were fixed). The 88 kg constant is the value printed in the equations; the pooled cohort mean was 89.15 kg (Table 1).",
      source_name        = "BWT"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = non-Asian (White, Black, or missing)",
      notes              = "Retained on both CL/F and Vc/F in the final model (Ahn 2023 Equations 11-12). Asian participants were 9.3% of the pooled cohort (42/451, Table 1), all 12 phase I Japanese volunteers plus 30 in the Western / phase II studies. The Discussion notes the bootstrapped 95% CI for the CL/F effect included the null value of 1 while the parametric CI did not; the covariate was retained because it was still of interest.",
      source_name        = "Asian"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = male",
      notes              = "Retained on CL/F only (Ahn 2023 Equation 11). Female participants were 40.6% of the pooled cohort (183/451, Table 1). The effect is after accounting for body weight, which is already in the model allometrically (Discussion).",
      source_name        = "Female"
    ),
    ADA_POS = list(
      description        = "Anti-drug antibody positive indicator, stationary (subject-level) definition",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = ADA-negative or not assessed",
      notes              = "Stationary definition: positive if any post-treatment ADA sample was positive, negative if none were (Methods 'Full PK model'). ADA was not assessed in the phase I studies and those subjects were assumed ADA-negative, because the median onset of treatment-emergent ADA was at least 164 days in the phase IIa study (Methods 'Missing data and imputations'). A time-varying ADA definition was also constructed but the stationary one was used for the model. Retained on CL/F (Ahn 2023 Equation 11).",
      source_name        = "ADAP"
    ),
    DOSE = list(
      description        = "Vupanorsen dose level of the treatment arm",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used only to build the binary 160 mg indicator that multiplies Q/F (Ahn 2023 Equation 13); all other dose levels share the reference Q/F. The effect captures a greater-than-dose-proportional exposure increase seen only in the 160 mg cohort of the Japanese phase I study, which the authors attribute to possible partial saturation of asialoglycoprotein-receptor-mediated hepatic uptake (Methods 'Full PK model'; Discussion). Supply the arm's nominal dose level in mg; the model tests DOSE == 160.",
      source_name        = "160mg"
    ),
    STUDY_PHASE2A = list(
      description        = "Phase IIa study indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = phase I studies (the PD reference population) or the phase IIb study",
      notes              = "NCT03371355, dose-finding in patients with hypertriglyceridemia, type 2 diabetes and nonalcoholic fatty liver disease (N = 105, Table 1). Carries its own baseline factor for each of the three PD endpoints (Ahn 2023 Equation 17, Table 3), and shares the single phase II potency factor with STUDY_PHASE2B (Equation 16). PK parameters do not depend on it.",
      source_name        = "Phase IIa"
    ),
    STUDY_PHASE2B = list(
      description        = "Phase IIb (TRANSLATE-TIMI 70) study indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = phase I studies (the PD reference population) or the phase IIa study",
      notes              = "NCT04516291, TRANSLATE-TIMI 70, dose-ranging in statin-treated patients with dyslipidemia (N = 286, Table 1). Carries its own baseline factor for each of the three PD endpoints (Ahn 2023 Equation 17, Table 3), and shares the single phase II potency factor with STUDY_PHASE2A (Equation 16). PK parameters do not depend on it.",
      source_name        = "Phase IIb"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Baseline age. Screened on CL/F in the full PK model as a power term (Ahn 2023 Equation 6, normalised to 60 years) and removed from the final model.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Judged not clinically important: the bootstrapped estimate lay close to the null value of 0 for a power exponent and the 95% CI was wide and included it (Results 'Full PK model'; Figure 2). Removing age, eGFR and Black race together raised the objective function by only about 2.5. The point estimate is reported graphically in Figure 2 only, so no usable value is available even for documentation. Cohort mean 59.5 years (SD 9.9), range 21-87 (Table 1)."
    ),
    CRCL = list(
      description        = "Baseline estimated glomerular filtration rate. Screened on CL/F in the full PK model as a power term (Ahn 2023 Equation 6, normalised to 94 mL/min/1.73 m^2) and removed from the final model.",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Judged not clinically important on the same grounds as AGE (Results 'Full PK model'; Figure 2). The paper does not state which estimating equation produced the eGFR values. Point estimate reported graphically in Figure 2 only. Cohort mean 91.30 mL/min/1.73 m^2 (SD 17.18), range 30.0-136.8 (Table 1)."
    ),
    RACE_BLACK = list(
      description        = "Black race indicator. Screened on CL/F in the full PK model (Ahn 2023 Equation 6) and removed from the final model.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = non-Black",
      notes              = "Judged not clinically important: the estimate lay close to the null value of 1 and the 95% CI was wide and included it (Results 'Full PK model'; Figure 2). Point estimate reported graphically in Figure 2 only. Black participants were 5.1% of the pooled cohort (23/451, Table 1). Note that RACE_ASIAN, screened in the same equation, WAS retained."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte  = "vupanorsen",
      units    = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte  = "vupanorsen",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte  = "vupanorsen",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    angptl3 = list(
      analyte  = "angiopoietin-like 3 protein",
      units    = "ng/mL",
      specimen = "serum",
      verified = TRUE
    ),
    tg = list(
      analyte  = "triglycerides",
      units    = "mg/dL",
      specimen = "serum",
      verified = TRUE
    ),
    nonhdlc = list(
      analyte  = "non-high-density-lipoprotein cholesterol",
      units    = "mg/dL",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 451,
    n_studies      = 4,
    age_range      = "21-87 years",
    age_median     = "mean 59.5 years (SD 9.9)",
    weight_range   = "52.0-138.0 kg",
    weight_median  = "mean 89.15 kg (SD 17.00)",
    sex_female_pct = 40.6,
    race_ethnicity = c(White = 84.7, Black = 5.1, Asian = 9.3, Missing = 0.9),
    disease_state  = "elevated triglycerides in otherwise healthy volunteers (phase I); hypertriglyceridemia with type 2 diabetes and nonalcoholic fatty liver disease (phase IIa); statin-treated dyslipidemia (phase IIb)",
    dose_range     = "20-120 mg subcutaneous single dose and 10-60 mg weekly in phase I; every-2-week and every-4-week regimens up to 160 mg in phase II",
    regions        = "Western (phase I NCT02709850, phase IIa NCT03371355, phase IIb NCT04516291) and Japan (phase I NCT04459767)",
    renal_function = "baseline eGFR mean 91.30 mL/min/1.73 m^2 (SD 17.18), range 30.0-136.8",
    co_medication  = "baseline statin use in 74.1% of the pooled cohort (334/451); all 286 phase IIb participants were statin-treated",
    notes          = "Baseline demographics from Ahn 2023 Table 1 (N = 451 including placebo recipients). The PK analysis used 2531 concentrations from 364 vupanorsen-treated participants plus one placebo participant with two quantifiable concentrations; the PD analysis used 3312 ANGPTL3, 3551 triglyceride and 3551 non-HDL-cholesterol observations from 451 participants (Results 'Observed data'). Observed baseline means were ANGPTL3 103.65 ng/mL (SD 35.11), triglycerides 248.6 mg/dL (SD 148.5) and non-HDL-cholesterol 149.72 mg/dL (SD 39.35) (Table S3); the model-estimated typical baselines in Table 3 differ from these because they are conditioned on the study-population factors."
  )

  ini({
    # ---- Structural PK, reference 88 kg body weight (Table 2) --------------
    lka <- log(0.737)     ; label("Absorption rate constant (Ka, 1/h)")                       # Table 2: Ka 0.737 /h
    lcl <- log(34.5)      ; label("Apparent clearance (CL/F, L/h)")                           # Table 2: CL/F 34.5 L/h
    lvc <- log(318)       ; label("Apparent central volume of distribution (Vc/F, L)")        # Table 2: Vc/F 318 L
    lq  <- log(8.49)      ; label("Apparent intercompartmental clearance (Q/F, L/h)")         # Table 2: Q/F 8.49 L/h
    lvp <- log(12100)     ; label("Apparent peripheral volume of distribution (Vp/F, L)")     # Table 2: Vp/F 12,100 L

    # ---- Allometric exponents ---------------------------------------------
    # Methods "Base PK model": "BWT was included with fixed allometry constants
    # (i.e., 0.75 and 1 for clearance and volume parameters, respectively)".
    # A single exponent is shared by CL/F and Q/F, and a single exponent by
    # Vc/F and Vp/F, per Equations 11-14.
    e_wt_cl_q  <- fixed(0.75) ; label("Allometric exponent of (WT/88) on CL/F and Q/F (unitless)")   # Methods "Base PK model"; Equations 11, 13
    e_wt_vc_vp <- fixed(1)    ; label("Allometric exponent of (WT/88) on Vc/F and Vp/F (unitless)")  # Methods "Base PK model"; Equations 12, 14

    # ---- Covariate effects on PK (multiplicative, applied as theta^indicator)
    e_race_asian_cl  <- 0.700 ; label("Multiplicative effect of Asian race on CL/F (unitless)")            # Table 2: Asian on CL/F 0.700
    e_sexf_cl        <- 0.820 ; label("Multiplicative effect of female sex on CL/F (unitless)")            # Table 2: Female on CL/F 0.820
    e_ada_pos_cl     <- 0.379 ; label("Multiplicative effect of ADA positivity on CL/F (unitless)")        # Table 2: ADAP on CL/F 0.379
    e_race_asian_vc  <- 0.519 ; label("Multiplicative effect of Asian race on Vc/F (unitless)")            # Table 2: Asian on Vc/F 0.519
    e_dose_160mg_q   <- 0.570 ; label("Multiplicative effect of the 160 mg dose level on Q/F (unitless)")  # Table 2: Dose 160 mg on Q/F 0.570

    # ---- PK interindividual variability -----------------------------------
    # Table 2 reports these under the heading "IIV", i.e. interindividual
    # VARIANCES (omega^2); they are used here unchanged.
    etalka ~ 0.318   # Table 2: IIV Ka 0.318
    etalcl ~ 0.445   # Table 2: IIV CL/F 0.445
    etalvc ~ 0.410   # Table 2: IIV Vc/F 0.410
    etalq  ~ 0.0883  # Table 2: IIV Q/F 0.0883
    etalvp ~ 0.709   # Table 2: IIV Vp/F 0.709

    # ---- PK residual error -------------------------------------------------
    # Equation 2 is a combined additive + proportional model on the linear
    # scale. Table 2 reports the residual unexplained VARIANCES; the entries
    # below are their square roots. Separate magnitudes were estimated for the
    # phase I and phase II studies (Methods "Base PK model").
    propSdPhase1 <- 0.3577709  ; label("Proportional residual SD, phase I studies (fraction)")     # Table 2: RUV prop phase I 0.128 -> sqrt(0.128)
    addSdPhase1  <- 0.06300794 ; label("Additive residual SD, phase I studies (ng/mL)")            # Table 2: RUV add phase I 0.00397 -> sqrt(0.00397)
    propSdPhase2 <- 0.4774935  ; label("Proportional residual SD, phase II studies (fraction)")    # Table 2: RUV prop phase II 0.228 -> sqrt(0.228)
    addSdPhase2  <- 0.03449638 ; label("Additive residual SD, phase II studies (ng/mL)")           # Table 2: RUV add phase II 0.00119 -> sqrt(0.00119)

    # ---- PD: typical baselines in the phase I reference population (Table 3)
    lrbase_angptl3 <- log(105) ; label("Typical baseline ANGPTL3 (ng/mL)")                   # Table 3: Baseline ANGPTL3 105 ng/mL
    lrbase_tg      <- log(186) ; label("Typical baseline triglycerides (mg/dL)")             # Table 3: Baseline TG 186 mg/dL
    lrbase_nonhdlc <- log(175) ; label("Typical baseline non-HDL-cholesterol (mg/dL)")       # Table 3: Baseline non-HDL-C 175 mg/dL

    # ---- PD: study-population baseline factors (Equation 17, FB_Patients) ---
    e_study_phase2a_rbase_angptl3 <- 0.978 ; label("Phase IIa factor on baseline ANGPTL3 (unitless)")                  # Table 3: Phase IIa factor in ANGPTL3 baseline 0.978
    e_study_phase2b_rbase_angptl3 <- 0.897 ; label("Phase IIb factor on baseline ANGPTL3 (unitless)")                  # Table 3: Phase IIb factor in ANGPTL3 baseline 0.897
    e_study_phase2a_rbase_tg      <- 1.46  ; label("Phase IIa factor on baseline triglycerides (unitless)")            # Table 3: Phase IIa factor in TG baseline 1.46
    e_study_phase2b_rbase_tg      <- 1.17  ; label("Phase IIb factor on baseline triglycerides (unitless)")            # Table 3: Phase IIb factor in TG baseline 1.17
    e_study_phase2a_rbase_nonhdlc <- 0.840 ; label("Phase IIa factor on baseline non-HDL-cholesterol (unitless)")      # Table 3: Phase IIa factor for non-HDL-C baseline 0.840
    e_study_phase2b_rbase_nonhdlc <- 0.768 ; label("Phase IIb factor on baseline non-HDL-cholesterol (unitless)")      # Table 3: Phase IIb factor for non-HDL-C baseline 0.768

    # ---- PD: potency in the phase I reference population (Table 3) ---------
    # IC50 values are peripheral-compartment vupanorsen concentrations.
    lec50_angptl3 <- log(0.929) ; label("IC50 for ANGPTL3 inhibition, phase I population (ng/mL)")                  # Table 3: IC50 ANGPTL3 0.929 ng/mL
    lec50_tg      <- log(0.741) ; label("IC50 for triglyceride inhibition, phase I population (ng/mL)")             # Table 3: IC50 TG 0.741 ng/mL
    lec50_nonhdlc <- log(6.10)  ; label("IC50 for non-HDL-cholesterol inhibition, phase I population (ng/mL)")      # Table 3: IC50 non-HDL-C 6.10 ng/mL

    # ---- PD: phase II potency factors (Equation 16, FP_Patients) -----------
    # One factor per end point, shared by the phase IIa and phase IIb studies
    # ("the drug potencies be different between phase I and phase II
    # populations", Results "Population PK/PD model results").
    e_study_phase2_ec50_angptl3 <- 3.60 ; label("Phase II factor on the ANGPTL3 IC50 (unitless)")                 # Table 3: Phase II factor for ANGPTL3 IC50 3.60
    e_study_phase2_ec50_tg      <- 5.13 ; label("Phase II factor on the triglyceride IC50 (unitless)")            # Table 3: Phase II factor for TG IC50 5.13
    e_study_phase2_ec50_nonhdlc <- 3.87 ; label("Phase II factor on the non-HDL-cholesterol IC50 (unitless)")     # Table 3: Phase II factor for non-HDL-C IC50 3.87

    # ---- PD: maximum inhibition (Table 3 footnote) -------------------------
    # "Imax fixed to 1, 0.815, and 0.690, respectively, based on the previous
    # developed pharmacokinetic/pharmacodynamic model parameter estimates using
    # healthy participants."
    limax_angptl3 <- fixed(log(1))     ; label("Maximum fractional inhibition of ANGPTL3 production (unitless)")                 # Table 3 footnote: Imax ANGPTL3 1
    limax_tg      <- fixed(log(0.815)) ; label("Maximum fractional inhibition of triglyceride production (unitless)")            # Table 3 footnote: Imax TG 0.815
    limax_nonhdlc <- fixed(log(0.690)) ; label("Maximum fractional inhibition of non-HDL-cholesterol production (unitless)")     # Table 3 footnote: Imax non-HDL-C 0.690

    # ---- PD: turnover rate constants (Table 3 footnote) --------------------
    # "Kout values for ANGPTL3, TG, and non-HDC-C were fixed to 0.0134, 0.0107,
    # and 0.0072 h-1, respectively, ... based on the previous developed
    # pharmacokinetic/pharmacodynamic model parameter estimates using healthy
    # participants."
    lkout_angptl3 <- fixed(log(0.0134)) ; label("First-order turnover rate constant for ANGPTL3 (1/h)")                 # Table 3 footnote: Kout ANGPTL3 0.0134 /h
    lkout_tg      <- fixed(log(0.0107)) ; label("First-order turnover rate constant for triglycerides (1/h)")           # Table 3 footnote: Kout TG 0.0107 /h
    lkout_nonhdlc <- fixed(log(0.0072)) ; label("First-order turnover rate constant for non-HDL-cholesterol (1/h)")     # Table 3 footnote: Kout non-HDL-C 0.0072 /h

    # ---- PD: sigmoidicity (Equation 3; Table 3) ----------------------------
    # "gamma, sigmoidicity (slope) factor (estimated for TG and non-HDL-C,
    # fixed to 1 for ANGPTL3)".
    lhill_angptl3 <- fixed(log(1)) ; label("Sigmoidicity factor for the ANGPTL3 inhibition function (unitless)")                 # Equation 3 text: gamma fixed to 1 for ANGPTL3
    lhill_tg      <- log(0.627)    ; label("Sigmoidicity factor for the triglyceride inhibition function (unitless)")            # Table 3: Gamma TG 0.627
    lhill_nonhdlc <- log(0.631)    ; label("Sigmoidicity factor for the non-HDL-C inhibition function (unitless)")               # Table 3: Gamma non-HDL-C 0.631

    # ---- PD interindividual variability (variances, as reported) -----------
    etalrbase_angptl3 ~ 0.0697  # Table 3: IIV ANGPTL3 baseline 0.0697
    etalec50_angptl3  ~ 0.814   # Table 3: IIV ANGPTL3 IC50 0.814
    etalrbase_tg      ~ 0.132   # Table 3: IIV TG baseline 0.132
    etalec50_tg       ~ 1.78    # Table 3: IIV TG IC50 1.78
    etalrbase_nonhdlc ~ 0.0577  # Table 3: IIV non-HDL-C baseline 0.0577
    etalec50_nonhdlc  ~ 1.79    # Table 3: IIV non-HDL-C IC50 1.79

    # ---- PD residual error (Equation 4, proportional) ----------------------
    # Table 3 reports residual unexplained VARIANCES; entries are square roots.
    propSd_angptl3 <- 0.2039608 ; label("Proportional residual SD for ANGPTL3 (fraction)")                  # Table 3: RUV ANGPTL3 0.0416 -> sqrt(0.0416)
    propSd_tg      <- 0.2437212 ; label("Proportional residual SD for triglycerides (fraction)")            # Table 3: RUV TG 0.0594 -> sqrt(0.0594)
    propSd_nonhdlc <- 0.1276715 ; label("Proportional residual SD for non-HDL-cholesterol (fraction)")      # Table 3: RUV non-HDL-C 0.0163 -> sqrt(0.0163)
  })

  model({
    # ---- 0. Unit conversion ------------------------------------------------
    # Doses are in mg and the apparent volumes in L, so amount / volume is
    # mg/L. The paper reports vupanorsen concentrations, the PD IC50 values and
    # the PK residual error in ng/mL, and 1 mg/L = 1000 ng/mL.
    ngmlPerMgl <- 1000

    # ---- 1. Derived covariate indicators -----------------------------------
    # Equation 13: the theta_160mg factor applies to the 160 mg dose level only.
    dose160 <- (DOSE > 159.5) * (DOSE < 160.5)
    # Equation 16: a single potency factor applies to both phase II studies
    # (the "patient" populations). Phase I is the reference. The two study
    # indicators are mutually exclusive, so their sum is the patient indicator.
    isPhase2 <- STUDY_PHASE2A + STUDY_PHASE2B

    # ---- 2. Individual PK parameters (Equations 11-15) ---------------------
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 88)^e_wt_cl_q *
      e_race_asian_cl^RACE_ASIAN * e_ada_pos_cl^ADA_POS * e_sexf_cl^SEXF
    vc <- exp(lvc + etalvc) * (WT / 88)^e_wt_vc_vp *
      e_race_asian_vc^RACE_ASIAN
    q  <- exp(lq + etalq) * (WT / 88)^e_wt_cl_q * e_dose_160mg_q^dose160
    vp <- exp(lvp + etalvp) * (WT / 88)^e_wt_vc_vp

    # ---- 3. Micro-constants ------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 4. PK disposition -------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ---- 5. Concentrations -------------------------------------------------
    Cc <- ngmlPerMgl * central / vc
    # C3 is the predicted PERIPHERAL-compartment concentration, which is what
    # drives the PD model: "Peripheral concentration was preferred to central
    # concentration as the site of vupanorsen action is known to be the liver"
    # (Methods "Population PK/PD model development"). Clamped at zero so that
    # solver round-off cannot raise a negative base to the fractional
    # sigmoidicity exponents estimated for TG and non-HDL-C.
    C3 <- max(ngmlPerMgl * peripheral1 / vp, 0.0)

    # ---- 6. PD parameters (Equations 16-17) --------------------------------
    kout_angptl3 <- exp(lkout_angptl3)
    kout_tg      <- exp(lkout_tg)
    kout_nonhdlc <- exp(lkout_nonhdlc)

    imax_angptl3 <- exp(limax_angptl3)
    imax_tg      <- exp(limax_tg)
    imax_nonhdlc <- exp(limax_nonhdlc)

    hill_angptl3 <- exp(lhill_angptl3)
    hill_tg      <- exp(lhill_tg)
    hill_nonhdlc <- exp(lhill_nonhdlc)

    # Equation 17: Kin = FB_Patients * Baseline_i * Kout, so the drug-free
    # steady state of each pool is FB_Patients * Baseline_i.
    rbase_angptl3 <- exp(lrbase_angptl3 + etalrbase_angptl3) *
      e_study_phase2a_rbase_angptl3^STUDY_PHASE2A *
      e_study_phase2b_rbase_angptl3^STUDY_PHASE2B
    rbase_tg <- exp(lrbase_tg + etalrbase_tg) *
      e_study_phase2a_rbase_tg^STUDY_PHASE2A *
      e_study_phase2b_rbase_tg^STUDY_PHASE2B
    rbase_nonhdlc <- exp(lrbase_nonhdlc + etalrbase_nonhdlc) *
      e_study_phase2a_rbase_nonhdlc^STUDY_PHASE2A *
      e_study_phase2b_rbase_nonhdlc^STUDY_PHASE2B

    kin_angptl3 <- rbase_angptl3 * kout_angptl3
    kin_tg      <- rbase_tg      * kout_tg
    kin_nonhdlc <- rbase_nonhdlc * kout_nonhdlc

    # Equation 16: the phase II potency factor multiplies IC50 inside the
    # sigmoid, i.e. the denominator term is (IC50_i * FP_Patients)^gamma.
    ec50_angptl3 <- exp(lec50_angptl3 + etalec50_angptl3) *
      e_study_phase2_ec50_angptl3^isPhase2
    ec50_tg <- exp(lec50_tg + etalec50_tg) *
      e_study_phase2_ec50_tg^isPhase2
    ec50_nonhdlc <- exp(lec50_nonhdlc + etalec50_nonhdlc) *
      e_study_phase2_ec50_nonhdlc^isPhase2

    # ---- 7. Indirect-response PD system (Equation 16) ----------------------
    inh_angptl3 <- imax_angptl3 * C3^hill_angptl3 /
      (ec50_angptl3^hill_angptl3 + C3^hill_angptl3)
    inh_tg <- imax_tg * C3^hill_tg /
      (ec50_tg^hill_tg + C3^hill_tg)
    inh_nonhdlc <- imax_nonhdlc * C3^hill_nonhdlc /
      (ec50_nonhdlc^hill_nonhdlc + C3^hill_nonhdlc)

    d/dt(angptl3) <- kin_angptl3 * (1 - inh_angptl3) - kout_angptl3 * angptl3
    d/dt(tg)      <- kin_tg      * (1 - inh_tg)      - kout_tg      * tg
    d/dt(nonhdlc) <- kin_nonhdlc * (1 - inh_nonhdlc) - kout_nonhdlc * nonhdlc

    angptl3(0) <- rbase_angptl3
    tg(0)      <- rbase_tg
    nonhdlc(0) <- rbase_nonhdlc

    # ---- 8. Observations and residual error --------------------------------
    # Equation 2, with the phase I / phase II residual magnitudes selected by
    # the study indicators (Methods "Base PK model").
    propSdCc <- propSdPhase1 * (1 - isPhase2) + propSdPhase2 * isPhase2
    addSdCc  <- addSdPhase1  * (1 - isPhase2) + addSdPhase2  * isPhase2

    Cc      ~ add(addSdCc) + prop(propSdCc)
    angptl3 ~ prop(propSd_angptl3)
    tg      ~ prop(propSd_tg)
    nonhdlc ~ prop(propSd_nonhdlc)
  })
}
