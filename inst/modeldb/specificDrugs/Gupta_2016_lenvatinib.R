Gupta_2016_lenvatinib <- function() {
  description <- "Three-compartment population PK model for lenvatinib in healthy subjects and patients with cancer (Gupta 2016). Simultaneous first-order plus zero-order oral absorption into the central compartment, linear elimination, and covariate effects of body weight (allometric on CL/F and Q/F with exponent 0.75 and linear on V/F), CYP3A4 inducers (+30 percent on CL/F), CYP3A4 inhibitors (-7.8 percent on CL/F), serum albumin < 30 g/L (-16.3 percent on CL/F), alkaline phosphatase > ULN (-11.7 percent on CL/F), healthy-subject cohort (+15 percent on CL/F vs cancer patients), and capsule vs tablet formulation (relative bioavailability 0.896)."
  reference <- "Gupta A, Jarzab B, Capdevila J, Shumaker R, Hussein Z. Population pharmacokinetic analysis of lenvatinib in healthy subjects and patients with cancer. Br J Clin Pharmacol. 2016 Jun;81(6):1124-33. doi:10.1111/bcp.12907"
  vignette <- "Gupta_2016_lenvatinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline; reported in kg).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric exponent 0.75 on CL/F, Q1/F, and Q2/F; linear allometry (exponent 1) on V1/F, V2/F, and V3/F per Gupta 2016 Table 2. Reference body weight 75 kg (Table 1 baseline median).",
      source_name        = "WGT"
    ),
    ALB = list(
      description        = "Serum albumin concentration (baseline or time-varying; the source paper uses the last available value if multiple assessments exist).",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Binarized inline as alb_low <- (ALB < 30) per Gupta 2016 Table 2 footnote: ALB indicator = 0 if ALB >= 30 g/L, 1 if ALB < 30 g/L. Multiplicative power-form effect on CL/F: 0.837^alb_low (-16.3 percent when alb_low = 1).",
      source_name        = "ALB"
    ),
    ALP = list(
      description        = "Serum alkaline phosphatase activity.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Binarized inline as alp_high <- (ALP > 120) per Gupta 2016 Table 2 footnote: ALP indicator = 0 if ALP/ULN <= 1, 1 if ALP/ULN > 1. Multiplicative power-form effect on CL/F: 0.883^alp_high (-11.7 percent when alp_high = 1). The 120 U/L ULN used here is a representative adult cutoff; downstream users should supply ALP_HIGH directly or adjust the inline threshold to match their site's ULN.",
      source_name        = "ALP"
    ),
    CYP3A4_IND = list(
      description        = "Concomitant CYP3A4 inducer coadministration indicator (1 = any CYP3A4 inducer during the study; 0 = no concomitant CYP3A4 inducer).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CYP3A4 inducer).",
      notes              = "Multiplicative power-form effect on CL/F: 1.30^CYP3A4_IND (+30 percent when 1). The Gupta dataset pools any concomitant CYP3A4 inducer reported in the per-subject medication log; 19 of 779 subjects (2.4 percent) were positive.",
      source_name        = "INDU"
    ),
    CYP3A4_INH = list(
      description        = "Concomitant CYP3A4 inhibitor coadministration indicator (1 = any CYP3A4 inhibitor during the study; 0 = no concomitant CYP3A4 inhibitor).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CYP3A4 inhibitor).",
      notes              = "Multiplicative power-form effect on CL/F: 0.922^CYP3A4_INH (-7.8 percent when 1). The Gupta dataset pools any concomitant CYP3A4 inhibitor reported in the per-subject medication log; 49 of 779 subjects (6.3 percent) were positive.",
      source_name        = "INHIB"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-subject cohort indicator (1 = healthy subject from a phase 1 clinical pharmacology study; 0 = cancer patient from a phase 1, 2, or 3 study).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (cancer patient; reference cohort is the pooled solid-tumor / thyroid-cancer cohort across phase 1-3 studies).",
      notes              = "Multiplicative power-form effect on CL/F: 1.15^DIS_HEALTHY (+15 percent when 1). Reflects the systematic CL/F difference between phase 1 clinical pharmacology subjects and the cancer-patient pool.",
      source_name        = "TM"
    ),
    FORM_CAPSULE = list(
      description        = "Capsule vs tablet formulation indicator (1 = capsule; 0 = tablet).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tablet; F fixed at 1 for the tablet reference arm in Gupta 2016).",
      notes              = "Multiplicative effect on bioavailability of both the first-order and the zero-order absorption routes: 0.896 (capsule) vs 1 (tablet) per Gupta 2016 Table 2. The 30.2 percent CV IIV on F1 (etalfcap) applies only to the capsule arm; tablet subjects have F = 1 with no eta contribution. Reference category in this model is tablet rather than the solution comparator used for itraconazole capsules in Hennig 2006/2007.",
      source_name        = "FORM"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 779,
    n_studies      = 15,
    age_range      = "18-89 years",
    age_median     = "55 years",
    weight_range   = "32.6-177.5 kg",
    weight_median  = "75 kg",
    sex_female_pct = 44.0,
    race_ethnicity = c(
      White            = 70.2,
      Japanese         = 11.7,
      Black            = 9.4,
      Other            = 6.3,
      Hispanic         = 0.8,
      Asian            = 0.6,
      Native_Hawaiian  = 0.6,
      American_Indian  = 0.4
    ),
    disease_state  = "Pooled cohort: healthy adults (n = 196, 25.2 percent), differentiated thyroid cancer (DTC; n = 327, 42.0 percent), medullary thyroid cancer (MTC; n = 56, 7.2 percent), anaplastic thyroid cancer (ATC; n = 9, 1.2 percent), and other solid tumors (n = 191, 24.5 percent).",
    dose_range     = "3.2-32 mg oral lenvatinib, mainly once daily, predominantly as tablets or capsules.",
    regions        = "Multiregional (15 pooled studies including phase 3 SELECT trial in RR-DTC).",
    n_observations = "10 265 plasma concentrations across 15 studies: 5 077 from phase 1 healthy-subject full profiles, 3 192 from phase 1 patient profiles plus sparse samples, 354 from phase 2 thyroid-cancer studies, and 1 642 from phase 3 SELECT (DTC).",
    ecog_distribution = "ECOG 0: 32.5 percent; ECOG 1: 26.7 percent; ECOG 2: 2.4 percent; ECOG 3: 0.1 percent; missing 38.3 percent.",
    notes          = "Demographic counts and ranges reproduced from Gupta 2016 Table 1. Missing race subjects are encoded as 0 for every minority RACE_* indicator (treated as the White reference). 56.0 percent male (436/779) and 44.0 percent female (343/779)."
  )

  ini({
    # Structural parameters - reference weight 75 kg, tablet reference (F = 1).
    lcl        <- log(6.56);   label("Apparent clearance, CL/F, at 75 kg in cancer patients on tablet without CYP3A4 inducers/inhibitors and with normal ALB / ALP (L/h)")  # Gupta 2016 Table 2: CL/F = 6.56 L/h
    lvc        <- log(49.3);   label("Apparent central volume, V1/F, at 75 kg (L)")                                                                                          # Gupta 2016 Table 2: V1/F = 49.3 L
    lvp        <- log(30.7);   label("Apparent peripheral 1 volume, V2/F, at 75 kg (L)")                                                                                     # Gupta 2016 Table 2: V2/F = 30.7 L
    lvp2       <- log(37.1);   label("Apparent peripheral 2 volume, V3/F, at 75 kg (L)")                                                                                     # Gupta 2016 Table 2: V3/F = 37.1 L
    lq         <- log(3.52);   label("Apparent inter-compartmental clearance, Q1/F, at 75 kg (L/h)")                                                                         # Gupta 2016 Table 2: Q1/F = 3.52 L/h
    lq2        <- log(0.769);  label("Apparent inter-compartmental clearance, Q2/F, at 75 kg (L/h)")                                                                         # Gupta 2016 Table 2: Q2/F = 0.769 L/h
    lka        <- log(1.02);   label("First-order absorption rate constant (1/h)")                                                                                           # Gupta 2016 Table 2: Ka = 1.02 1/h
    lduration  <- log(1.22);   label("Duration of zero-order absorption, D1 (h)")                                                                                            # Gupta 2016 Table 2: D1 = 1.22 h
    lfcap      <- log(0.896);  label("Relative bioavailability of capsule vs tablet (unitless)")                                                                             # Gupta 2016 Table 2: F1 = 0.896

    # Allometric exponents held fixed at the paper's reported integer / canonical values.
    allo_cl    <- fixed(0.75); label("Allometric exponent on CL/F and Q/F (unitless)")                                                                                       # Gupta 2016 Table 2: WGT/75 raised to 0.75 in CL/F and Q1/F, Q2/F covariate equations
    allo_v     <- fixed(1.0);  label("Allometric exponent on V1/F, V2/F, V3/F (unitless)")                                                                                   # Gupta 2016 Table 2: WGT/75 raised to 1 (linear) in V1/F, V2/F, V3/F covariate equations

    # Covariate effects on CL/F expressed on the log scale: source THETAs
    # are entered as exp(e_<cov>_cl * <cov>) so that the linear-space
    # multiplier reproduces the source-paper power form theta^cov.
    e_cyp3a4_ind_cl <- log(1.30);  label("Log-effect of CYP3A4 inducer on CL/F (unitless)")  # Gupta 2016 Table 2: theta_INDU = 1.30 (+30 percent on CL/F)
    e_cyp3a4_inh_cl <- log(0.922); label("Log-effect of CYP3A4 inhibitor on CL/F (unitless)")  # Gupta 2016 Table 2: theta_INHIB = 0.922 (-7.8 percent on CL/F)
    e_alb_cl        <- log(0.837); label("Log-effect of low albumin (< 30 g/L) on CL/F (unitless)")  # Gupta 2016 Table 2: theta_ALB = 0.837 (-16.3 percent on CL/F)
    e_alp_cl        <- log(0.883); label("Log-effect of high alkaline phosphatase (> ULN) on CL/F (unitless)")  # Gupta 2016 Table 2: theta_ALP = 0.883 (-11.7 percent on CL/F)
    e_dis_healthy_cl <- log(1.15); label("Log-effect of healthy-subject cohort on CL/F (unitless)")  # Gupta 2016 Table 2: theta_TM = 1.15 (+15 percent on CL/F)

    # IIV - reported as percent CV in Gupta 2016 Table 2 (footnote a:
    # %CV ~ sqrt(variance) * 100). Variances stored here use the
    # log-normal conversion omega^2 = log(CV^2 + 1) so that the simulated
    # log-scale eta reproduces the linear-space CV when exponentiated.
    etalcl       ~ 0.06296    # Gupta 2016 Table 2: 25.5 percent CV -> log(0.255^2 + 1) = 0.06296
    etalvc       ~ 0.05074    # Gupta 2016 Table 2: 22.8 percent CV -> log(0.228^2 + 1) = 0.05074
    etalvp       ~ 0.14164    # Gupta 2016 Table 2: 39.0 percent CV -> log(0.390^2 + 1) = 0.14164
    etalvp2      ~ 0.08788    # Gupta 2016 Table 2: 30.3 percent CV -> log(0.303^2 + 1) = 0.08788
    etalka       ~ 0.26236    # Gupta 2016 Table 2: 54.8 percent CV -> log(0.548^2 + 1) = 0.26236
    etalduration ~ 0.46258    # Gupta 2016 Table 2: 76.7 percent CV -> log(0.767^2 + 1) = 0.46258
    etalfcap     ~ 0.08732    # Gupta 2016 Table 2: 30.2 percent CV -> log(0.302^2 + 1) = 0.08732 (applies to the capsule arm only)

    # Residual error - Gupta 2016 Table 2 reports a stratified
    # residual model (proportional 16.9 percent for clinical pharmacology
    # studies, 33.3 percent for patient studies, and a separate
    # 48.1 percent proportional + 7.19 ng/mL additive arm for TAD < 2 h).
    # The simplified combined form below uses the patient-study
    # proportional plus the early-absorption additive term as a single
    # combined error model suitable for simulation. The vignette
    # "Assumptions and deviations" section quotes the alternative arms.
    propSd <- 0.333;  label("Proportional residual error (fraction)")  # Gupta 2016 Table 2: 33.3 percent CV (patient studies)
    addSd  <- 7.19;   label("Additive residual error (ng/mL)")          # Gupta 2016 Table 2: 7.19 ng/mL (TAD < 2 h additive arm)
  })

  model({
    # Reference covariate values - Gupta 2016 Table 1 baseline medians.
    ref_wt   <- 75
    alp_uln  <- 120  # Representative adult ALP upper limit of normal (U/L); see covariateData[[ALP]]$notes.

    # Binary indicators derived from continuous lab values - match the
    # paper's NONMEM coding (Table 2 footnote).
    alb_low  <- (ALB < 30)
    alp_high <- (ALP > alp_uln)

    # Individual structural parameters.
    ka       <- exp(lka + etalka)
    duration <- exp(lduration + etalduration)

    cl <- exp(lcl + etalcl) *
          (WT / ref_wt)^allo_cl *
          exp(e_cyp3a4_ind_cl  * CYP3A4_IND) *
          exp(e_cyp3a4_inh_cl  * CYP3A4_INH) *
          exp(e_alb_cl         * alb_low) *
          exp(e_alp_cl         * alp_high) *
          exp(e_dis_healthy_cl * DIS_HEALTHY)

    vc  <- exp(lvc  + etalvc)  * (WT / ref_wt)^allo_v
    vp  <- exp(lvp  + etalvp)  * (WT / ref_wt)^allo_v
    vp2 <- exp(lvp2 + etalvp2) * (WT / ref_wt)^allo_v

    q   <- exp(lq)  * (WT / ref_wt)^allo_cl
    q2  <- exp(lq2) * (WT / ref_wt)^allo_cl

    # Formulation-effect bioavailability: 1 for the tablet reference arm,
    # 0.896 * exp(eta) for the capsule arm. The 30.2 percent CV IIV on F1
    # therefore applies only to capsule subjects.
    formul_bio <- (1 - FORM_CAPSULE) + FORM_CAPSULE * exp(lfcap + etalfcap)

    # ODE system. The simultaneous first- and zero-order absorption model
    # routes the dose into two compartments in parallel: a first-order
    # depot with rate ka, and a direct-to-central zero-order infusion
    # with duration D1. The 50:50 dose split between the two routes is a
    # documented assumption (the source paper does not disclose the
    # explicit fraction); see vignette Assumptions and deviations.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (cl + q + q2) / vc * central +
                                       q / vp  * peripheral1 +
                                       q2 / vp2 * peripheral2
    d/dt(peripheral1) <-  q  / vc * central - q  / vp  * peripheral1
    d/dt(peripheral2) <-  q2 / vc * central - q2 / vp2 * peripheral2

    # Zero-order direct-to-central infusion duration.
    dur(central) <- duration

    # Bioavailability split: 0.5 of the administered dose enters via the
    # first-order route (depot) and 0.5 enters via the zero-order route
    # (central, with duration D1). Both routes share the formulation
    # effect (formul_bio).
    f(depot)   <- 0.5 * formul_bio
    f(central) <- 0.5 * formul_bio

    Cc <- central / vc * 1000  # 1 mg/L = 1000 ng/mL (dose in mg, volume in L).
    Cc ~ prop(propSd) + add(addSd)
  })
}
