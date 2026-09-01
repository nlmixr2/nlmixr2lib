Bhagunde_2026_lecanemab_abeta4240 <- function() {
  description <- paste(
    "Indirect-response PK/PD model for the absolute plasma amyloid-beta",
    "42/40 ratio in subjects with early Alzheimer's disease receiving the",
    "anti-protofibril monoclonal antibody lecanemab. Serum lecanemab",
    "concentration Cc drives a linear stimulation of the Abeta42/40",
    "production rate: dR/dt = Kin * (1 + SLOPE * Cc) - Kout * R, with",
    "Kin = baseline * Kout so the untreated pool sits at its baseline.",
    "Covariate effects (all log-linear) are age, APOE4-carrier status,",
    "Japanese race and female sex on the baseline ratio, and age and body",
    "weight on SLOPE. Interindividual variability is exponential on",
    "baseline and SLOPE with an estimated correlation, and residual error",
    "is proportional. Fit by NONMEM FOCE-I to 12,468 plasma Abeta42/40",
    "observations from 1994 subjects pooled across the lecanemab phase 2",
    "Study 201 (Core and open-label extension) and phase 3 Clarity AD",
    "(Study 301, Core and open-label extension) trials. The lecanemab",
    "serum-concentration driver is the two-compartment population PK model",
    "of Majid 2024 (zero-order IV input, linear elimination from the",
    "central compartment); its parameters are held fixed here because the",
    "present paper did not re-estimate them.",
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

  paper_specific_compartments <- c("abeta4240")

  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(
      analyte = "lecanemab", units = "mg", specimen = "serum", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "lecanemab", units = "mg", specimen = "serum", verified = TRUE
    ),
    abeta4240 = list(
      analyte = "amyloid-beta 42/40 ratio", units = "ratio (unitless)",
      specimen = "plasma", verified = TRUE
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
        "normalised to 72 kg (Majid 2024 Table 1 equations). (2) On the",
        "Abeta42/40 drug-effect SLOPE as a log-linear (exponential) term",
        "centred at 72 kg. The 72 kg centring value is the reference",
        "weight printed in the Majid 2024 popPK equations and is",
        "independently recovered to 71.8-72.0 kg from the Bhagunde 2026",
        "Figure 4C body-weight forest panel; see the vignette Assumptions",
        "and deviations section.",
        sep = " "
      ),
      source_name        = "BW / Baseline Weight"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power term on lecanemab CL normalised to 43 g/L",
        "(Majid 2024 Table 1 equation: CL = 0.0154 * (BW/72)^0.353 *",
        "(ALB/43)^-0.374 * 0.791^SEX * 1.13^ADA).",
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
        "(ratio 0.868) in the Majid 2024 popPK model, and on the",
        "Abeta42/40 baseline ratio (log-linear coefficient -0.00888) in",
        "the Bhagunde 2026 PD model. The source papers write SEX with the",
        "same 0 = male / 1 = female orientation, so no value",
        "transformation is required.",
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
        "Table 1). Time-varying in the source popPK analysis; may be",
        "supplied as a time-varying column.",
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
        "0.671) in the Majid 2024 popPK model, and on the Abeta42/40",
        "baseline ratio (log-linear coefficient 0.0346) in the",
        "Bhagunde 2026 PD model.",
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
        "Log-linear (exponential) covariate on both the Abeta42/40",
        "baseline ratio and SLOPE, centred at 72 years -- the reference",
        "subject age stated in Bhagunde 2026 Section 3.2.1 ('Compared to",
        "a reference subject aged 72 years'). Centring at 72 years",
        "reproduces the Figure 4A age forest panel (83 years -> 1.24 vs",
        "the plotted 1.21; 57 years -> 0.75 vs the plotted 0.77).",
        sep = " "
      ),
      source_name        = "AGE"
    ),
    APOE4_CARRIER = list(
      description        = "APOE-epsilon4 carrier status",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-carrier)",
      notes              = paste(
        "1 = carries at least one APOE-epsilon4 allele (heterozygous or",
        "homozygous), 0 = non-carrier. Log-linear coefficient 0.0138 on",
        "the Abeta42/40 baseline ratio (Bhagunde 2026 Table 1). The",
        "source paper tests carrier-vs-noncarrier only and does not",
        "separate heterozygotes from homozygotes, so the binary",
        "APOE4_CARRIER encoding is used rather than APOE4_COUNT.",
        sep = " "
      ),
      source_name        = "APOE4 carrier"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1994,
    n_studies      = 2,
    n_observations = 12468,
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
      "Pooled Core and open-label-extension data from the lecanemab phase",
      "2 Study 201 (856 randomised across six treatment groups, 18-month",
      "Core, 9-59 month off-treatment gap, then open-label 10 mg/kg",
      "biweekly) and the phase 3 Clarity AD / Study 301 trial (18-month",
      "Core, 1:1 placebo vs 10 mg/kg biweekly, plus open-label",
      "extension). Subject and observation counts per biomarker are",
      "Bhagunde 2026 Table S1; demographic detail is in the study-design",
      "publications cited by Bhagunde 2026 Appendix A and is not",
      "tabulated in the present paper.",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------- #
    # Lecanemab population PK -- two-compartment, zero-order IV input,  #
    # linear elimination. All values are FIXED: Bhagunde 2026 states    #
    # the population PK model has been published elsewhere, and uses     #
    # post hoc individual PK parameters as the PD driver, so none of    #
    # these was re-estimated in the present paper. Source of every      #
    # value: Majid 2024 Table 1 and the equations printed beneath it.   #
    # ---------------------------------------------------------------- #
    lcl <- fixed(log(0.0154)); label("Lecanemab clearance at the reference subject (CL, L/h)")           # Majid 2024 Table 1: CL = 0.0154 L/h (RSE 1.60%)
    lvc <- fixed(log(3.24));   label("Lecanemab central volume at the reference subject (V1, L)")        # Majid 2024 Table 1: V1 = 3.24 L (RSE 0.799%)
    lvp <- fixed(log(2.00));   label("Lecanemab peripheral volume at the reference subject (V2, L)")     # Majid 2024 Table 1: V2 = 2.00 L (RSE 4.09%)
    lq  <- fixed(log(0.00718)); label("Lecanemab intercompartmental clearance (Q, L/h)")                 # Majid 2024 Table 1: Q = 0.00718 L/h (RSE 4.23%)

    e_wt_cl            <- fixed(0.353);  label("Power exponent for body weight on CL, (WT/72)^e (unitless)")            # Majid 2024 Table 1: Weight ~ CL exponent 0.353
    e_alb_cl           <- fixed(-0.374); label("Power exponent for serum albumin on CL, (ALB/43)^e (unitless)")         # Majid 2024 Table 1: Albumin ~ CL exponent -0.374
    e_sexf_cl          <- fixed(0.791);  label("Ratio of CL in females to males, applied as e^SEXF (unitless)")         # Majid 2024 Table 1: Females ~ CL ratio to males 0.791
    e_ada_pos_cl       <- fixed(1.13);   label("Ratio of CL in ADA-positive to ADA-negative subjects, applied as e^ADA_POS (unitless)")  # Majid 2024 Table 1: ADA positive ~ CL ratio 1.13
    e_wt_vc            <- fixed(0.513);  label("Power exponent for body weight on V1, (WT/72)^e (unitless)")            # Majid 2024 Table 1: Weight ~ V1 exponent 0.513
    e_sexf_vc          <- fixed(0.868);  label("Ratio of V1 in females to males, applied as e^SEXF (unitless)")         # Majid 2024 Table 1: Females ~ V1 ratio to males 0.868
    e_race_japanese_vc <- fixed(0.920);  label("Ratio of V1 in Japanese to non-Japanese subjects, applied as e^RACE_JAPANESE (unitless)")  # Majid 2024 Table 1: Japanese ethnicity ~ V1 ratio 0.920
    e_race_japanese_vp <- fixed(0.671);  label("Ratio of V2 in Japanese to non-Japanese subjects, applied as e^RACE_JAPANESE (unitless)")  # Majid 2024 Table 1: Japanese ethnicity ~ V2 ratio 0.671

    # ---------------------------------------------------------------- #
    # Plasma Abeta42/40 indirect-response PD model                      #
    # (Bhagunde 2026 Section 2.2.2 equation and Table 1)                #
    # ---------------------------------------------------------------- #
    lrbase <- log(0.0864);    label("Baseline plasma Abeta42/40 ratio at the reference subject (unitless ratio)")   # Bhagunde 2026 Table 1: Baseline plasma Abeta42/40 ratio = 0.0864 (RSE 0.272%)
    lkout  <- log(1.51);      label("First-order degradation rate constant of the plasma Abeta42/40 pool (Kout, 1/year)")  # Bhagunde 2026 Table 1: Kout = 1.51 /year (RSE 9.63%), t1/2 = ln(2)/1.51 = 0.46 year, matching the paper approximately 0.5 years
    lslope <- log(0.000704);  label("Linear slope of the lecanemab-concentration effect on Abeta42/40 production (SLOPE, per ug/mL)")  # Bhagunde 2026 Table 1: Slope for lecanemab exposure effect = 0.000704 per ug/mL (RSE 3.89%)

    e_age_rbase            <- -0.00181; label("Log-linear coefficient for age on baseline Abeta42/40, exp(e * (AGE - 72)) (per year)")            # Bhagunde 2026 Table 1: Age on baseline (exponential) = -0.00181 (RSE 15.4%)
    e_apoe4_carrier_rbase  <- 0.0138;   label("Log-linear coefficient for APOE4-carrier status on baseline Abeta42/40, exp(e * APOE4_CARRIER) (unitless)")  # Bhagunde 2026 Table 1: APOE4 carrier on baseline = 0.0138 (RSE 26.9%)
    e_race_japanese_rbase  <- 0.0346;   label("Log-linear coefficient for Japanese race on baseline Abeta42/40, exp(e * RACE_JAPANESE) (unitless)")  # Bhagunde 2026 Table 1: Japanese (vs others) on baseline = 0.0346 (RSE 14.4%)
    e_sexf_rbase           <- -0.00888; label("Log-linear coefficient for female sex on baseline Abeta42/40, exp(e * SEXF) (unitless)")            # Bhagunde 2026 Table 1: Females (vs males) on baseline = -0.00888 (RSE 38.3%)
    e_age_slope            <- 0.0193;   label("Log-linear coefficient for age on SLOPE, exp(e * (AGE - 72)) (per year)")                          # Bhagunde 2026 Table 1: Age on slope (exponential) = 0.0193 (RSE 19.4%)
    e_wt_slope             <- -0.00602; label("Log-linear coefficient for body weight on SLOPE, exp(e * (WT - 72)) (per kg)")                     # Bhagunde 2026 Table 1: Weight on slope (exponential) = -0.00602 (RSE 25.4%)

    # Lecanemab PK interindividual variability, fixed from Majid 2024
    # (omega^2 = log(CV^2 + 1). Majid 2024 Table 1 footnote gives
    #  CV percent = square root of variance x 100 for the residual error and
    #  the standard exponential-IIV %CV definition for the etas).
    etalcl + etalvc ~ fixed(c(0.114936,
                              0.005934, 0.014774))   # Majid 2024 Table 1: IIV CL 34.9 %CV, IIV V1 12.2 %CV, correlation CL ~ V1 R = 0.144
    etalvp          ~ fixed(0.639175)                # Majid 2024 Table 1: IIV V2 94.6 %CV

    # Plasma Abeta42/40 interindividual variability (estimated here)
    etalrbase + etalslope ~ c(0.005185,
                              -0.026570, 0.414814)   # Bhagunde 2026 Table 1: IIV baseline 7.21 %CV, IIV slope 71.7 %CV, correlation baseline ~ slope R = -0.573

    propSd <- 0.0772; label("Proportional residual error on the plasma Abeta42/40 ratio (fraction)")  # Bhagunde 2026 Table 1: Proportional residual variability 7.72 %CV (RSE 2.12%)
  })

  model({
    # Constants and reference (centring) values
    hrs_per_year <- 8766   # 365.25 days/year x 24 h/day; converts the PD rate constants, reported per year, to the model time unit (h)
    wt_ref       <- 72     # kg;   Majid 2024 popPK equations normalise WT to 72 kg
    alb_ref      <- 43     # g/L;  Majid 2024 popPK equations normalise ALB to 43 g/L
    age_ref      <- 72     # years; Bhagunde 2026 Section 3.2.1 reference subject

    # Individual lecanemab PK parameters (Majid 2024 Table 1 equations)
    cl <- exp(lcl + etalcl) * (WT / wt_ref)^e_wt_cl * (ALB / alb_ref)^e_alb_cl *
      e_sexf_cl^SEXF * e_ada_pos_cl^ADA_POS
    vc <- exp(lvc + etalvc) * (WT / wt_ref)^e_wt_vc * e_sexf_vc^SEXF *
      e_race_japanese_vc^RACE_JAPANESE
    vp <- exp(lvp + etalvp) * e_race_japanese_vp^RACE_JAPANESE
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Individual plasma Abeta42/40 PD parameters
    rbase <- exp(lrbase + etalrbase) *
      exp(e_age_rbase * (AGE - age_ref) +
            e_apoe4_carrier_rbase * APOE4_CARRIER +
            e_race_japanese_rbase * RACE_JAPANESE +
            e_sexf_rbase * SEXF)
    slope <- exp(lslope + etalslope) *
      exp(e_age_slope * (AGE - age_ref) + e_wt_slope * (WT - wt_ref))
    kout  <- exp(lkout) / hrs_per_year
    kin   <- rbase * kout    # Bhagunde 2026: Kin is set so the untreated pool sits at baseline

    # Lecanemab disposition
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Serum lecanemab concentration -- the PD driver. This is an individual
    # prediction (IPRED), matching the paper's use of post hoc individual PK
    # parameters; it therefore carries no residual error of its own.
    Cc <- central / vc

    # Plasma Abeta42/40 indirect response with linear stimulation of production
    # dR/dt = Kin * (1 + SLOPE * Conc) - R(t) * Kout   (Bhagunde 2026 Section 2.2.2)
    d/dt(abeta4240) <- kin * (1 + slope * Cc) - kout * abeta4240
    abeta4240(0)    <- rbase

    abeta4240 ~ prop(propSd)
  })
}
