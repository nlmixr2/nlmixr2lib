Bhagunde_2026_lecanemab_gfap <- function() {
  description <- paste(
    "Indirect-response PK/PD model with a disease-progression term for",
    "plasma glial fibrillary acidic protein (GFAP) in subjects with early",
    "Alzheimer's disease receiving the anti-protofibril monoclonal",
    "antibody lecanemab. Unlike the companion Abeta42/40 and p-tau181",
    "models, the drug effect on GFAP is driven not by lecanemab",
    "concentration but by the relative reduction in brain amyloid plaque",
    "from baseline, which had the lowest objective function value of the",
    "structures tested:",
    "dGFAP/dt = KinG * (1 - SLP * (Plaque_baseline - Plaque(t)) /",
    "Plaque_baseline + DP * TIME) - KoutG * GFAP, with",
    "KoutG = KinG / BGFAP at baseline. The amyloid-plaque state is the",
    "semi-mechanistic amyloid-PET turnover model of Bhagunde 2026",
    "(CPT:PSP), in which lecanemab stimulates plaque elimination linearly",
    "in serum concentration, and lecanemab serum concentration comes from",
    "the two-compartment population PK model of Majid 2024. Both upstream",
    "layers are held fixed because the present paper did not re-estimate",
    "them. Covariates on baseline GFAP are power terms for age, body",
    "weight and baseline Mini-Mental State Examination score plus a",
    "study-effect term for Study 201 (an assay batch effect). Variability",
    "is exponential on BGFAP and KinG, additive on the logit of SLP, and",
    "residual error is proportional. Fit by Monolix SAEM to 3098 plasma",
    "GFAP observations from 736 subjects pooled across the lecanemab",
    "phase 2 Study 201 and phase 3 Clarity AD (Study 301) trials.",
    sep = " "
  )
  reference <- paste(
    "Bhagunde P, Penner N, Willis BA, Bell R, Sachdev P, Charil A,",
    "Irizarry MC, Hersch S, Reyderman L (2026).",
    "Pharmacokinetic/pharmacodynamic analyses of plasma pathophysiology",
    "biomarkers in subjects with early Alzheimer's disease following",
    "lecanemab treatment.",
    "Alzheimer's & Dementia: Translational Research & Clinical",
    "Interventions 12(2):e70246. doi:10.1002/trc2.70246.",
    "Amyloid-plaque driver (baseline on the logit scale, Kout, DESLP and",
    "their covariate effects) taken from the cited upstream model:",
    "Bhagunde P, Penner N, Willis BA, et al. (2026) Brain amyloid plaque",
    "levels affect clinical progression in Alzheimer disease: assessment",
    "of amyloid PET and change in CDR-SB utilizing semi-mechanistic",
    "model. CPT Pharmacometrics Syst Pharmacol 15(2):e70173.",
    "doi:10.1002/psp4.70173.",
    "Lecanemab population PK driver (CL, V1, V2, Q and their covariate",
    "effects) taken from the cited upstream model:",
    "Majid O, Cao Y, Willis BA, et al. (2024) Population pharmacokinetics",
    "and exposure-response analyses of safety (ARIA-E and isolated ARIA-H)",
    "of lecanemab in subjects with early Alzheimer's disease.",
    "CPT Pharmacometrics Syst Pharmacol 13(12):2111-2123.",
    "doi:10.1002/psp4.13224.",
    sep = " "
  )
  vignette <- "Bhagunde_2026_lecanemab_biomarkers"

  paper_specific_compartments <- c("plaque", "gfap")

  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(
      analyte = "lecanemab", units = "mg", specimen = "serum", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "lecanemab", units = "mg", specimen = "serum", verified = TRUE
    ),
    plaque = list(
      analyte = "brain amyloid-beta plaque burden (amyloid PET)",
      units = "CL (Centiloid)", specimen = "not applicable", verified = TRUE
    ),
    gfap = list(
      analyte = "glial fibrillary acidic protein (GFAP)",
      units = "pg/mL", specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters twice. (1) On lecanemab CL and V1 as a power term",
        "normalised to 72 kg (Majid 2024 Table 1 equations). (2) On",
        "baseline plasma GFAP as a power term (WT/72)^-0.66 (Bhagunde",
        "2026 Table 1). The 72 kg normalisation is confirmed by the",
        "Bhagunde 2026 Figure 4C weight forest panel, whose plotted fold",
        "changes (50 kg -> 1.27, 100 kg -> 0.805) imply reference weights",
        "of 71.8 and 72.0 kg under the tabulated exponent -0.66.",
        sep = " "
      ),
      source_name        = "Bodyweight"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power term on lecanemab CL normalised to 43 g/L",
        "(Majid 2024 Table 1 equation).",
        sep = " "
      ),
      source_name        = "ALB"
    ),
    SEXF = list(
      description        = "Sex",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "1 = female, 0 = male. Acts on lecanemab CL (ratio 0.791) and V1",
        "(ratio 0.868) in the Majid 2024 popPK model. Sex entered the",
        "Bhagunde 2026 GFAP forward-addition step but was eliminated",
        "during backward elimination (Table S5).",
        sep = " "
      ),
      source_name        = "SEX"
    ),
    ADA_POS = list(
      description        = "Anti-drug antibody status",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = paste(
        "1 = ADA-positive. Ratio 1.13 on lecanemab CL (Majid 2024",
        "Table 1). Time-varying in the source popPK analysis.",
        sep = " "
      ),
      source_name        = "ADA"
    ),
    RACE_JAPANESE = list(
      description        = "Japanese-heritage race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Japanese)",
      notes              = paste(
        "1 = Japanese. Acts on lecanemab V1 (ratio 0.920) and V2 (ratio",
        "0.671) in the Majid 2024 popPK model.",
        sep = " "
      ),
      source_name        = "JPN"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters twice. (1) On baseline plasma GFAP as a power term",
        "(AGE/72)^1.03 (Bhagunde 2026 Table 1). (2) On the amyloid-plaque",
        "drug-effect slope DESLP as a power term (AGE/72)^3.10 (Bhagunde",
        "2026 CPT:PSP Table 1). Neither source paper prints the",
        "normalisation age; 72 years is the reference-subject age stated",
        "in Bhagunde 2026 Section 3.2.1 and is used for both. See the",
        "vignette Assumptions and deviations section.",
        sep = " "
      ),
      source_name        = "AGE"
    ),
    SCORE_MMSE = list(
      description        = "Baseline Mini-Mental State Examination total score",
      units              = "(SCORE_MMSE units, 0-30 score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power term on baseline plasma GFAP, (SCORE_MMSE/26)^-0.63",
        "(Bhagunde 2026 Table 1). The paper does not print the",
        "normalisation value; 26 is recovered from the Bhagunde 2026",
        "Figure 4C baseline-MMSE forest panel (22 units -> 1.11 implies",
        "25.96; 29 units -> 0.929 implies 25.80). See the vignette",
        "Assumptions and deviations section.",
        sep = " "
      ),
      source_name        = "BMMSE"
    ),
    APOE4_CARRIER = list(
      description        = "APOE-epsilon4 carrier status",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-carrier)",
      notes              = paste(
        "1 = carries at least one APOE-epsilon4 allele, 0 = non-carrier.",
        "Acts on the logit-scale baseline amyloid plaque burden in the",
        "upstream amyloid-PET model as a multiplicative ratio 0.629 on",
        "the logit-scale typical value (Bhagunde 2026 CPT:PSP Table 1),",
        "which raises the typical baseline plaque from 61.5 to 82.7",
        "Centiloids and reproduces that paper's stated 62 and 83",
        "Centiloids. APOE4 was not retained as a covariate on baseline",
        "GFAP itself (Bhagunde 2026 Table S5).",
        sep = " "
      ),
      source_name        = "APOE4 carrier"
    ),
    STUDY_LEC201 = list(
      description        = "Lecanemab Study 201 (phase 2) cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Clarity AD / Study 301)",
      notes              = paste(
        "1 = subject enrolled in the lecanemab phase 2 Study 201,",
        "0 = subject enrolled in the phase 3 Clarity AD (Study 301)",
        "trial. Log-linear effect 0.388 on baseline plasma GFAP",
        "(Bhagunde 2026 Table 1). The authors attribute the effect to an",
        "assay batch effect: 'GFAP samples from Study 201 and Clarity AD",
        "were analyzed at different times' (Section 3.2.3). Only the GFAP",
        "model carries this covariate.",
        sep = " "
      ),
      source_name        = "Study"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 736,
    n_studies      = 2,
    n_observations = 3098,
    disease_state  = paste(
      "Early Alzheimer's disease (mild cognitive impairment due to",
      "Alzheimer's disease, or mild Alzheimer's disease dementia) with",
      "confirmed brain amyloid pathology.",
      sep = " "
    ),
    dose_range     = paste(
      "Lecanemab IV 2.5 mg/kg biweekly, 5 mg/kg monthly, 5 mg/kg",
      "biweekly, 10 mg/kg monthly, 10 mg/kg biweekly, or placebo",
      "(Bhagunde 2026 Table S1).",
      sep = " "
    ),
    regions        = "Global (North America, Europe, Asia including Japan)",
    notes          = paste(
      "The GFAP analysis dataset was restricted to subjects who also had",
      "amyloid PET data, because the drug effect is driven by relative",
      "change in amyloid plaque. Plasma GFAP was evaluated for only a",
      "subset of Study 201 subjects (about 35 per dose group in the Core",
      "phase) and was measured with the Quanterix Simoa GFAP Discovery",
      "Kit (Bhagunde 2026 Appendix A). Subject and observation counts are",
      "Bhagunde 2026 Table S1.",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------- #
    # Layer 1 -- lecanemab population PK (two-compartment, zero-order   #
    # IV input, linear elimination). All values FIXED: not re-estimated #
    # by Bhagunde 2026. Source: Majid 2024 Table 1 and its equations.   #
    # ---------------------------------------------------------------- #
    lcl <- fixed(log(0.0154)); label("Lecanemab clearance at the reference subject (CL, L/h)")        # Majid 2024 Table 1: CL = 0.0154 L/h (RSE 1.60%)
    lvc <- fixed(log(3.24));   label("Lecanemab central volume at the reference subject (V1, L)")     # Majid 2024 Table 1: V1 = 3.24 L (RSE 0.799%)
    lvp <- fixed(log(2.00));   label("Lecanemab peripheral volume at the reference subject (V2, L)")  # Majid 2024 Table 1: V2 = 2.00 L (RSE 4.09%)
    lq  <- fixed(log(0.00718)); label("Lecanemab intercompartmental clearance (Q, L/h)")              # Majid 2024 Table 1: Q = 0.00718 L/h (RSE 4.23%)

    e_wt_cl            <- fixed(0.353);  label("Power exponent for body weight on CL, (WT/72)^e (unitless)")            # Majid 2024 Table 1
    e_alb_cl           <- fixed(-0.374); label("Power exponent for serum albumin on CL, (ALB/43)^e (unitless)")         # Majid 2024 Table 1
    e_sexf_cl          <- fixed(0.791);  label("Ratio of CL in females to males, applied as e^SEXF (unitless)")         # Majid 2024 Table 1
    e_ada_pos_cl       <- fixed(1.13);   label("Ratio of CL in ADA-positive to ADA-negative subjects, applied as e^ADA_POS (unitless)")  # Majid 2024 Table 1
    e_wt_vc            <- fixed(0.513);  label("Power exponent for body weight on V1, (WT/72)^e (unitless)")            # Majid 2024 Table 1
    e_sexf_vc          <- fixed(0.868);  label("Ratio of V1 in females to males, applied as e^SEXF (unitless)")         # Majid 2024 Table 1
    e_race_japanese_vc <- fixed(0.920);  label("Ratio of V1 in Japanese to non-Japanese subjects, applied as e^RACE_JAPANESE (unitless)")  # Majid 2024 Table 1
    e_race_japanese_vp <- fixed(0.671);  label("Ratio of V2 in Japanese to non-Japanese subjects, applied as e^RACE_JAPANESE (unitless)")  # Majid 2024 Table 1

    # ---------------------------------------------------------------- #
    # Layer 2 -- brain amyloid plaque (amyloid PET) turnover model.     #
    # dPlaque/dt = Kin - Plaque * Kout * (1 + DESLP * Conc), printed in #
    # Bhagunde 2026 Section 2.2.2. Parameter values from the upstream    #
    # Bhagunde 2026 CPT:PSP paper, which the present paper cites as the #
    # source of the plaque predictions. All FIXED (not re-estimated).   #
    # ---------------------------------------------------------------- #
    logitrbase_plaque <- fixed(-1.12);        label("Typical baseline amyloid plaque burden on the logit scale, Plaque_0 = 250 * plogis(logit) Centiloid")  # Bhagunde 2026 CPT:PSP Table 1: baseline amyloid plaque, logit scale = -1.12 (RSE 5.35%). 250 CL is the stated maximum possible plaque
    lkout_plaque      <- fixed(log(0.0572));  label("First-order amyloid plaque elimination rate constant (Kout, 1/year)")  # Bhagunde 2026 CPT:PSP Table 1: Kout = 0.0572 /year (RSE 24.4%), t1/2 = 12.1 years
    lslope_plaque     <- fixed(log(0.154));   label("Linear slope of the lecanemab-concentration effect on plaque elimination (DESLP, per ug/mL)")  # Bhagunde 2026 CPT:PSP Table 1: Drug effect (DESLP) = 0.154 (RSE 23.7%)

    e_apoe4_carrier_logitrbase_plaque <- fixed(0.629); label("Ratio applied to the logit-scale baseline plaque for APOE4 carriers, e^APOE4_CARRIER (unitless)")  # Bhagunde 2026 CPT:PSP Table 1: APOE4 carrier on baseline amyloid PET (ratio) = 0.629 (RSE 6.36%)
    e_age_slope_plaque                <- fixed(3.10);  label("Power exponent for age on DESLP, (AGE/72)^e (unitless)")  # Bhagunde 2026 CPT:PSP Table 1: Age on drug effect = 3.10 (RSE 7.82%)

    # ---------------------------------------------------------------- #
    # Layer 3 -- plasma GFAP indirect response with disease progression #
    # (Bhagunde 2026 Section 2.2.2 equation and Table 1)                #
    # ---------------------------------------------------------------- #
    lrbase <- log(315); label("Baseline plasma GFAP at the reference subject (BGFAP, pg/mL)")                   # Bhagunde 2026 Table 1: Baseline GFAP = 315 pg/mL (RSE 1.66%)
    lkin   <- log(5.6); label("Zero-order plasma GFAP production rate at the reference subject (KinG, pg/mL/h)")  # Bhagunde 2026 Table 1: KinG = 5.6 pg/mL/h (RSE 23.5%), KoutG is derived as KinG / BGFAP
    logitslope_gfap <- log(0.237 / (1 - 0.237)); label("Logit of the slope of the relative-plaque-reduction effect on GFAP production (SLP, unitless in (0,1))")  # Bhagunde 2026 Table 1: SLP (logit) = 0.237 (RSE 8.33%), logit(0.237) = -1.16920, matching the additive-on-logit IIV
    lprog  <- log(0.0307); label("Linear disease-progression drift on GFAP production (DP, 1/year)")            # Bhagunde 2026 Table 1: DP = 0.0307 /year (RSE 0.279%)

    e_age_rbase          <- 1.03;  label("Power exponent for age on baseline GFAP, (AGE/72)^e (unitless)")                    # Bhagunde 2026 Table 1: Age on baseline GFAP (power) = 1.03 (RSE 12.1%)
    e_wt_rbase           <- -0.66; label("Power exponent for body weight on baseline GFAP, (WT/72)^e (unitless)")             # Bhagunde 2026 Table 1: Weight on baseline GFAP (power) = -0.66 (RSE 10.2%)
    e_mmse_rbase         <- -0.63; label("Power exponent for baseline MMSE on baseline GFAP, (SCORE_MMSE/26)^e (unitless)")   # Bhagunde 2026 Table 1: Baseline MMSE on baseline GFAP (power) = -0.63 (RSE 25.8%)
    e_study_lec201_rbase <- 0.388; label("Log-linear effect of Study 201 enrolment on baseline GFAP, exp(e * STUDY_LEC201) (unitless)")  # Bhagunde 2026 Table 1: Study 201 on baseline GFAP = 0.388 (RSE 9.01%), assay batch effect per Section 3.2.3

    # Lecanemab PK interindividual variability, fixed from Majid 2024
    # (omega^2 = log(CV^2 + 1)).
    etalcl + etalvc ~ fixed(c(0.114936,
                              0.005934, 0.014774))   # Majid 2024 Table 1: IIV CL 34.9 %CV, IIV V1 12.2 %CV, correlation CL ~ V1 R = 0.144
    etalvp          ~ fixed(0.639175)                # Majid 2024 Table 1: IIV V2 94.6 %CV

    # Plasma GFAP interindividual variability (estimated here)
    etalrbase          ~ 0.136213   # Bhagunde 2026 Table 1: IIV baseline GFAP 38.2 %CV -> omega^2 = log(1 + 0.382^2)
    etalkin            ~ 1.018451   # Bhagunde 2026 Table 1: IIV KinG 133 %CV -> omega^2 = log(1 + 1.33^2)
    etalogitslope_gfap ~ 0.866761   # Bhagunde 2026 Table 1: SLP additive on logit, SD 0.931 -> variance 0.931^2

    propSd <- 0.193; label("Proportional residual error on plasma GFAP (fraction)")  # Bhagunde 2026 Table 1: Proportional error 0.193 (RSE 1.56%), reported under a percent-CV heading but the value is the SD, i.e. 19.3 percent -- see the vignette Errata
  })

  model({
    # Constants and reference (centring) values
    hrs_per_year  <- 8766   # 365.25 days/year x 24 h/day; converts rate constants reported per year to the model time unit (h)
    wt_ref        <- 72     # kg;  Majid 2024 popPK equations normalise WT to 72 kg
    alb_ref       <- 43     # g/L; Majid 2024 popPK equations normalise ALB to 43 g/L
    age_ref       <- 72     # years; Bhagunde 2026 Section 3.2.1 reference subject
    mmse_ref      <- 26     # MMSE units; normalisation recovered from Bhagunde 2026 Figure 4C (see covariateData notes)
    plaque_max    <- 250    # Centiloid; Bhagunde 2026 CPT:PSP: "250 represents the maximum possible amyloid plaque"

    # --- Layer 1: individual lecanemab PK (Majid 2024 Table 1 equations)
    cl <- exp(lcl + etalcl) * (WT / wt_ref)^e_wt_cl * (ALB / alb_ref)^e_alb_cl *
      e_sexf_cl^SEXF * e_ada_pos_cl^ADA_POS
    vc <- exp(lvc + etalvc) * (WT / wt_ref)^e_wt_vc * e_sexf_vc^SEXF *
      e_race_japanese_vc^RACE_JAPANESE
    vp <- exp(lvp + etalvp) * e_race_japanese_vp^RACE_JAPANESE
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Serum lecanemab concentration -- the driver of plaque removal. This is an
    # individual prediction (IPRED), matching the paper's use of post hoc
    # individual PK parameters; it carries no residual error of its own.
    Cc <- central / vc

    # --- Layer 2: brain amyloid plaque turnover
    # BSL = 250 * exp(phi) / (1 + exp(phi)); phi = TVphi * COV
    phi_plaque   <- logitrbase_plaque * e_apoe4_carrier_logitrbase_plaque^APOE4_CARRIER
    rbase_plaque <- plaque_max * exp(phi_plaque) / (1 + exp(phi_plaque))
    kout_plaque  <- exp(lkout_plaque) / hrs_per_year
    slope_plaque <- exp(lslope_plaque) * (AGE / age_ref)^e_age_slope_plaque
    kin_plaque   <- rbase_plaque * kout_plaque   # kin = BSL * kout

    d/dt(plaque) <- kin_plaque - kout_plaque * (1 + slope_plaque * Cc) * plaque
    plaque(0)    <- rbase_plaque

    # Relative reduction in amyloid plaque from baseline -- the GFAP driver
    plaque_red <- (rbase_plaque - plaque) / rbase_plaque

    # --- Layer 3: plasma GFAP indirect response with disease progression
    rbase <- exp(lrbase + etalrbase) *
      (AGE / age_ref)^e_age_rbase *
      (WT / wt_ref)^e_wt_rbase *
      (SCORE_MMSE / mmse_ref)^e_mmse_rbase *
      exp(e_study_lec201_rbase * STUDY_LEC201)
    kin  <- exp(lkin + etalkin)
    kout <- kin / rbase   # Bhagunde 2026: "KinG = BGFAP * KoutG (at baseline)"
    # Kept on its own line so rxode2 recognises the mu-referenced expression
    logit_slp <- logitslope_gfap + etalogitslope_gfap
    slp  <- 1 / (1 + exp(-logit_slp))
    prog <- exp(lprog) / hrs_per_year

    # dGFAP/dt = KinG * (1 - SLP * relative plaque reduction + DP * TIME) - KoutG * GFAP
    d/dt(gfap) <- kin * (1 - slp * plaque_red + prog * t) - kout * gfap
    gfap(0)    <- rbase

    gfap ~ prop(propSd)
  })
}
