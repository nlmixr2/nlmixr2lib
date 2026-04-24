Li_2019_abatacept <- function() {
  description <- "Two-compartment population PK model for abatacept (CTLA4-Ig Fc-fusion) in adults with rheumatoid arthritis (Li 2019), with first-order SC absorption, zero-order IV infusion support, first-order linear elimination, logit-scale SC bioavailability, full-block IIV on CL/VC/Q/VP, and a KA parameterisation that enforces KA > k_el."
  reference <- "Li X, Roy A, Murthy B. Population Pharmacokinetics and Exposure-Response Relationship of Intravenous and Subcutaneous Abatacept in Patients With Rheumatoid Arthritis. J Clin Pharmacol. 2019 Feb;59(2):245-257. doi:10.1002/jcph.1308"
  vignette <- "Li_2019_abatacept"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 70 kg. Power effect on CL (exp 0.651), VC (exp 0.452), and VP (exp 0.457) per Li 2019 Table 1B and the final-model covariate equation in Methods/Results.",
      source_name        = "BWT"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 50 years. Power effect on CL (exp -0.186) per Li 2019 Table 1B.",
      source_name        = "AGE"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 4.0 g/dL. Li 2019 Methods states 'baseline albumin of 4.0 mg/dL' which is a publication unit typo (normal human albumin is ~4 g/dL; 4 mg/dL is physiologically impossible). Coded as 4.0 g/dL here; see vignette Errata. Power effect on CL (exp -0.687) per Li 2019 Table 1B.",
      source_name        = "ALB"
    ),
    CRCL = list(
      description        = "Calculated glomerular filtration rate (BSA-normalized)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 90 mL/min/1.73 m^2. Li 2019 uses 'cGFR' (calculated GFR, BSA-normalized); mapped to the canonical CRCL (which is defined to accept either MDRD-estimated eGFR or BSA-normalized measured CrCl). Power effect on CL (exp 0.162) per Li 2019 Table 1B.",
      source_name        = "cGFR"
    ),
    SWOL_28JOINT = list(
      description        = "Baseline swollen joint count (28-joint scale)",
      units              = "count (0-28)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 16. Applied as ((SWOL_28JOINT + 1)/(16 + 1))^0.0965 on CL to avoid the zero-count edge case per the Li 2019 final-model covariate equation. Not clinically relevant per Li 2019.",
      source_name        = "SWOL"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; reference subject is male per Li 2019 Methods)",
      notes              = "Li 2019 codes SEX = 0 for males (reference) and SEX = 1 for females; the CL ~ SEX coefficient is -0.0722 (CL lower in females). SEXF is a direct one-to-one rename of the Li 2019 SEX column.",
      source_name        = "SEX"
    ),
    CONMED_NSAID = list(
      description        = "Concomitant NSAID use at baseline",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on concomitant NSAIDs; typical patient)",
      notes              = "Exponential effect on CL: CL multiplied by exp(0.0640 * CONMED_NSAID) per Li 2019 final-model covariate equation. Patients on NSAIDs have ~6.6% higher CL (not clinically relevant per Li 2019). Baseline-only in Li 2019.",
      source_name        = "NSAID"
    ),
    FORM_ABA_PHASE2 = list(
      description        = "Abatacept SC formulation indicator, 1 = phase-2 SC formulation, 0 = phase-3 (commercial) SC formulation",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (phase-3 commercial SC formulation)",
      notes              = "Model-specific covariate: formulation indicators are kept non-canonical per the nlmixr2lib global policy unless they clearly generalize across multiple drugs. Added on the logit scale to the typical F: logit_F = logit_F_TV + CONMED_FORM_ABA_PHASE2 * (-1.16). The phase-2 SC formulation had a different pH and a lower absolute bioavailability (~0.56) than the phase-3/commercial formulation (~0.81, the reference). Set to 0 for routine simulation of the commercial 125 mg SC regimen.",
      source_name        = "FORM"
    )
  )

  population <- list(
    n_subjects     = 2244L,
    n_observations = 10382L,
    n_studies      = 11L,
    disease_state  = "Rheumatoid arthritis (adult); 4 phase 2 and 7 phase 3 studies pooled.",
    dose_range     = "IV abatacept 0.5-10 mg/kg Q4W (6 studies) and SC abatacept 75-200 mg QW (4 studies); one study (ACQUIRE) included both IV and SC. Approved regimens are weight-tiered ~10 mg/kg IV Q4W and fixed 125 mg SC QW.",
    regions        = "Multi-regional (11 pooled global phase 2 and phase 3 studies).",
    reference_subject = "50-year-old male, BWT 70 kg, baseline albumin 4.0 g/dL (publication typo: 'mg/dL'), cGFR 90 mL/min/1.73 m^2, swollen joint count 16, not on concomitant NSAIDs, phase-3 SC formulation; reference values approximate the median (continuous) or mode (categorical) of the popPK dataset.",
    notes          = "Baseline demographics per Li 2019 Methods Data and Study Populations (Table S3 referenced for full demographic summary). After exclusion of samples missing dose/sample information and below-LLOQ concentrations, 10 382/13 610 (76.3%) samples from 2244 patients entered the analysis dataset. LLOQ of the validated ELISA was 1.0 ng/mL. SC phase-2 formulation (different pH) was replaced by the phase-3 commercial formulation to improve product stability."
  )

  ini({
    # Structural PK parameters - Li 2019 Table 1A final-model estimates. Paper
    # reports CL, Q in L/h and KA in 1/h; values below convert to per-day
    # (x 24) to match the nlmixr2lib convention of time in days.
    lcl         <- log(0.0204 * 24);    label("Clearance CL (L/day) at reference covariates")                  # Li 2019 Table 1A: CL_TV,ref = 0.0204 L/h
    lvc         <- log(3.27);           label("Central volume VC (L) at reference covariates")                 # Li 2019 Table 1A: VC_TV,ref = 3.27 L
    lq          <- log(0.0265 * 24);    label("Inter-compartmental clearance Q (L/day)")                       # Li 2019 Table 1A: Q_TV,ref = 0.0265 L/h
    lvp         <- log(4.26);           label("Peripheral volume VP (L) at reference covariates")              # Li 2019 Table 1A: VP_TV,ref = 4.26 L
    lka         <- log(0.00305 * 24);   label("Typical absorption rate KA_TV (1/day)")                         # Li 2019 Table 1A: KA_TV = 0.00305 1/h
    logitfdepot <- 1.42;                label("Logit-scale typical SC bioavailability F_TV,ref (unitless)")     # Li 2019 Table 1A: F_TV,ref = 1.42 (logit scale); F_abs = 1/(1+exp(-1.42)) ~= 0.805

    # Covariate effects - Li 2019 Table 1B; see vignette Errata regarding the
    # 'Covariates of the full PPK model' caption. Values below are the
    # published estimates for the 10 covariate effects retained in the final
    # model (MTX was dropped by backward elimination and is not included in
    # the Methods/Results final-model equation).
    e_wt_cl    <-  0.651;  label("Power exponent of (WT/70 kg) on CL (unitless)")                              # Li 2019 Table 1B: CL~BWT = 0.651
    e_wt_vc    <-  0.452;  label("Power exponent of (WT/70 kg) on VC (unitless)")                              # Li 2019 Table 1B: VC~BWT = 0.452
    e_wt_vp    <-  0.457;  label("Power exponent of (WT/70 kg) on VP (unitless)")                              # Li 2019 Table 1B: VP~BWT = 0.457
    e_age_cl   <- -0.186;  label("Power exponent of (AGE/50 yr) on CL (unitless)")                             # Li 2019 Table 1B: CL~AGE = -0.186
    e_alb_cl   <- -0.687;  label("Power exponent of (ALB/4.0 g/dL) on CL (unitless)")                          # Li 2019 Table 1B: CL~ALB = -0.687 (reference 4.0 g/dL; paper typo 'mg/dL')
    e_crcl_cl  <-  0.162;  label("Power exponent of (CRCL/90 mL/min/1.73m^2) on CL (unitless)")                # Li 2019 Table 1B: CL~cGFR = 0.162
    e_swol_cl  <-  0.0965; label("Power exponent of ((SWOL_28JOINT+1)/(16+1)) on CL (unitless)")               # Li 2019 Table 1B: CL~SWOL = 0.0965
    e_sexf_cl  <- -0.0722; label("Exponential coefficient on CL for female sex (SEXF=1; unitless)")             # Li 2019 Table 1B: CL~SEX = -0.0722
    e_nsaid_cl <-  0.0640; label("Exponential coefficient on CL for concomitant NSAID use (unitless)")          # Li 2019 Table 1B: CL~NSAID = 0.0640
    e_form_f   <- -1.16;   label("Additive coefficient on logit-F for SC phase-2 formulation (unitless)")       # Li 2019 Table 1B: F~FORM = -1.16

    # Inter-individual variability - Li 2019 Table 1A reports variances
    # (omega^2) with SDs in parentheses. Full-block matrix on CL, VC, Q, VP;
    # independent IIV on KA and on logit-F.
    # Block lower-triangle (row-major), Li 2019 Table 1A values:
    #   (1,1) var(ZCL)=0.0991
    #   (2,1) cov(ZCL,ZVC)=0.0412, (2,2) var(ZVC)=0.0632
    #   (3,1) cov(ZCL,ZQ)=0.0952, (3,2) cov(ZVC,ZQ)=0.0407, (3,3) var(ZQ)=0.429
    #   (4,1) cov(ZCL,ZVP)=0.0910, (4,2) cov(ZVC,ZVP)=0.0675, (4,3) cov(ZQ,ZVP)=0.280, (4,4) var(ZVP)=0.377
    etalcl + etalvc + etalq + etalvp ~
      c(0.0991,
        0.0412, 0.0632,
        0.0952, 0.0407, 0.429,
        0.0910, 0.0675, 0.280, 0.377)
    etalka         ~ 1.63    ; label("IIV variance on the relative log-KA (Li 2019 Table 1A var(ZKA); independent)")
    etalogitfdepot ~ 0.710   ; label("IIV variance on logit-F (Li 2019 Table 1A var(ZF); independent)")

    # Residual error - Li 2019 Table 1A combined proportional + additive.
    propSd <- 0.215; label("Proportional residual error (fraction)")                     # Li 2019 Table 1A: theta_PROP = 0.215
    addSd  <- 0.341; label("Additive residual error (mg/L = ug/mL)")                     # Li 2019 Table 1A: theta_ADD = 0.341 ug/mL
  })

  model({
    # Individual PK parameters with Li 2019 final-model covariate equations
    # (Methods text and Table 1B). Reference subject: 50-year-old male, BWT
    # 70 kg, ALB 4.0 g/dL (paper typo 'mg/dL'), cGFR 90 mL/min/1.73 m^2,
    # swollen joint count 16, not on NSAIDs, phase-3 SC formulation.
    cl <- exp(lcl + etalcl) *
          (WT / 70)^e_wt_cl *
          (AGE / 50)^e_age_cl *
          (ALB / 4.0)^e_alb_cl *
          (CRCL / 90)^e_crcl_cl *
          ((SWOL_28JOINT + 1) / (16 + 1))^e_swol_cl *
          exp(SEXF * e_sexf_cl + CONMED_NSAID * e_nsaid_cl)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp

    # KA is parameterised per Li 2019 to ensure KA > k_el (flip-flop
    # prevention): KA_i = KA_TV * exp(etaKA) + CL_i / VC_i.
    kel <- cl / vc
    ka  <- exp(lka + etalka) + kel

    # Logit-scale SC bioavailability with independent IIV; the phase-2
    # formulation indicator is additive on the logit scale. F_abs is bounded
    # in (0, 1) via the inverse logit.
    logit_f <- logitfdepot + etalogitfdepot + FORM_ABA_PHASE2 * e_form_f
    fdepot  <- 1 / (1 + exp(-logit_f))

    # Two-compartment PK; IV doses go directly to central, SC doses to depot
    # with bioavailability fdepot.
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # Concentration in mg/L (= ug/mL), matching the addSd unit.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
