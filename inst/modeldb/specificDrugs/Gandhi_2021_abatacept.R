Gandhi_2021_abatacept <- function() {
  description <- "Two-compartment population PK model for abatacept (CTLA4-Ig Fc-fusion) pooled across adults with rheumatoid arthritis and patients aged 2-17 years with polyarticular juvenile idiopathic arthritis (Gandhi 2021), with first-order SC absorption, zero-order IV infusion support, first-order linear elimination, logit-scale SC bioavailability with disease/age/weight covariates, and a KA parameterisation that enforces KA > k_el."
  reference <- "Gandhi V, Sun H, Subramanian K, Lon HK, Roy A. Model-Based Selection and Recommendation for Subcutaneous Abatacept Dose in Patients With Polyarticular Juvenile Idiopathic Arthritis. J Clin Pharmacol. 2021 May;61(5):651-661. doi:10.1002/jcph.1781"
  vignette <- "Gandhi_2021_abatacept"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 68 kg (Gandhi 2021 Figure 1 caption: reference patient weighs 68 kg). Power effect on CL (exp 0.706), VC (exp 0.603), VP (exp 0.575), and on logit-F (slope -0.506).",
      source_name        = "BWT"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 49 years (Gandhi 2021 Figure 1 caption: 49-year-old reference patient). Power effect on VC (exp 0.114) and on logit-F (slope 0.487). Not retained on CL in the final model.",
      source_name        = "AGE"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 4.1 g/dL (Gandhi 2021 Figure 1 caption). Power effect on CL (exp -0.722).",
      source_name        = "ALB"
    ),
    CRCL = list(
      description        = "Calculated glomerular filtration rate (BSA-normalized)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 99.18 mL/min/1.73 m^2 (Gandhi 2021 Figure 1 caption). Gandhi 2021 uses 'cGFR' (calculated GFR, BSA-normalized); mapped to the canonical CRCL (which accepts either MDRD-estimated eGFR or BSA-normalized measured CrCl). Power effect on CL (exp 0.259).",
      source_name        = "cGFR"
    ),
    SWOL_28JOINT = list(
      description        = "Baseline swollen joint count (28-joint scale)",
      units              = "count (0-28)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 15 (Gandhi 2021 Figure 1 caption). Applied as ((SWOL_28JOINT + 1)/(15 + 1))^0.0742 on CL to avoid the zero-count edge case (matches the shifted-power form used in Li 2019 from the same author group). The paper does not explicitly state the joint-count scale; the 28-joint scale is assumed because (1) the same author group used 28-joint counts in the prior Li 2019 RA-only analysis and (2) the reference value 15 is consistent with the 28-joint scale (cap 28).",
      source_name        = "SJC"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; reference subject is male per Gandhi 2021 Methods 'except for sex, for which male was used as the reference')",
      notes              = "Gandhi 2021 Methods: 'for sex, male was used as the reference'. SEXF = 1 for female, 0 for male; the published 'Exponent of male sex on CL = 0.0674' is applied as `exp(SEXF * 0.0674)` so females have ~7% higher CL than males in the pooled RA + pJIA dataset (sign opposite to Li 2019's RA-only -0.0722; the dataset expansion to pJIA produced the sign flip). Convention follows the existing Li 2019 implementation in nlmixr2lib.",
      source_name        = "SEX"
    ),
    CONMED_NSAID = list(
      description        = "Concomitant NSAID use at baseline",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on concomitant NSAIDs; the reference subject is not on NSAIDs per Gandhi 2021 Figure 1 caption)",
      notes              = "Exponential effect on CL: CL multiplied by exp(0.102 * CONMED_NSAID) per Gandhi 2021 Table 2 'Exponent of NSAID on CL = 0.102'. Patients on NSAIDs have ~10.7% higher CL (not clinically relevant per Gandhi 2021).",
      source_name        = "NSAID"
    ),
    DIS_PJIA = list(
      description        = "Polyarticular juvenile idiopathic arthritis disease-state indicator, 1 = pJIA, 0 = adult RA (or other non-pJIA)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult RA; the typical reference subject in Gandhi 2021 is an RA patient)",
      notes              = "Gandhi 2021 codes a pJIA-vs-RA disease indicator (JIA = 1 if pJIA, 0 if RA). The covariate enters the bioavailability model only (additive on the logit scale: `logit_F = logit_F_TV + 3.08 * DIS_PJIA + ...`); pJIA patients have substantially higher SC bioavailability than adult RA patients (logit_F shifts from 1.208 at the reference to ~4.288 in pJIA, i.e. F_abs ~ 0.770 -> ~0.987). Disease (pJIA vs RA) was tested but did NOT have a clinically relevant effect on abatacept CL (Gandhi 2021 Results paragraph 1 of Population PK Analysis section).",
      source_name        = "JIA"
    )
  )

  population <- list(
    n_subjects     = 2616L,
    n_observations = 12759L,
    n_studies      = 13L,
    age_range      = "2-17 years (pJIA cohort) plus adult RA cohort",
    weight_range   = "10-100+ kg (simulated range; actual study weight range not separately tabulated in PMC text)",
    disease_state  = "Pooled adult rheumatoid arthritis (n = 2213) and polyarticular juvenile idiopathic arthritis (pJIA; n = 403, ages 2-17 years).",
    dose_range     = "IV abatacept 0.5-10 mg/kg Q4W and SC abatacept 75-200 mg QW (RA studies); approved IV regimen ~10 mg/kg Q4W; weight-tiered SC pJIA regimen 50 mg (<25 kg), 87.5 mg (25-<50 kg), 125 mg (>=50 kg) QW; approved adult SC regimen 125 mg QW.",
    regions        = "Multi-regional (13 pooled phase 2/3 studies).",
    reference_subject = "49-year-old female RA patient, BWT 68 kg, baseline albumin 4.1 g/dL, calculated GFR 99.18 mL/min/1.73 m^2, swollen joint count 15, not on NSAIDs (Gandhi 2021 Figure 1 caption). For the typical-value parameter values (CL_TYP, etc.), the reference is male per the Gandhi 2021 Methods text 'for sex, male was used as the reference'; the reference patient's female designation in Figure 1 is the visualization baseline used for the covariate-effect forest plot.",
    notes          = "Pooled population PK dataset of 12 759 abatacept serum concentrations (RA, 9420; pJIA, 3339) from 2616 patients (RA, 2213; pJIA, 403) drawn from 13 phase 2/3 studies (Table S1, not embedded in PMC). Concentrations were quantified by validated ELISA (LLOQ 1.0 ng/mL); 8.7% of samples below LLOQ were excluded. The study population spans intravenous abatacept (RA + pJIA) and subcutaneous abatacept (RA + pJIA), so the analysis dataset includes mixed routes per subject in some studies. Gandhi 2021 Methods Data and Study Populations."
  )

  ini({
    # Structural PK parameters - Gandhi 2021 Table 2 final-model estimates.
    # Paper reports CL, Q in L/h and KA in 1/h; values below convert to per-day
    # (x 24) to match the nlmixr2lib convention of time in days. The KA units
    # column header in Table 2 is mis-typed as 'L/h' in the paper; KA is a
    # rate constant, treated here as 1/h consistent with the rest of Methods.
    lcl         <- log(0.0179 * 24);    label("Clearance CL (L/day) at reference covariates")                  # Gandhi 2021 Table 2: CL_TV = 0.0179 L/h
    lvc         <- log(3.29);           label("Central volume VC (L) at reference covariates")                 # Gandhi 2021 Table 2: VC_TV = 3.29 L
    lq          <- log(0.0231 * 24);    label("Inter-compartmental clearance Q (L/day)")                       # Gandhi 2021 Table 2: Q_TV = 0.0231 L/h
    lvp         <- log(3.67);           label("Peripheral volume VP (L) at reference covariates")              # Gandhi 2021 Table 2: VP_TV = 3.67 L
    lka         <- log(0.00521 * 24);   label("Typical relative absorption rate KA_TV (1/day)")                 # Gandhi 2021 Table 2: KA_TV = 0.00521 1/h (paper 'L/h' is a unit typo for the rate constant)
    logitfdepot <- 1.20831;             label("Logit-scale typical SC bioavailability F_TV,ref (unitless)")     # Gandhi 2021 Table 2: SC F_TV = 0.770; logit(0.770) = 1.20831

    # Covariate effects - Gandhi 2021 Table 2. Reference covariates per Figure 1
    # caption: 49-year-old, BWT 68 kg, ALB 4.1 g/dL, cGFR 99.18 mL/min/1.73 m^2,
    # SJC 15, no NSAIDs, RA (DIS_PJIA = 0), male sex (per Methods).
    e_wt_cl    <-  0.706;  label("Power exponent of (WT/68 kg) on CL (unitless)")                              # Gandhi 2021 Table 2: 'Power of body weight on CL' = 0.706
    e_wt_vc    <-  0.603;  label("Power exponent of (WT/68 kg) on VC (unitless)")                              # Gandhi 2021 Table 2: 'Power of body weight on VC' = 0.603
    e_wt_vp    <-  0.575;  label("Power exponent of (WT/68 kg) on VP (unitless)")                              # Gandhi 2021 Table 2: 'Power of body weight on VP' = 0.575
    e_wt_f     <- -0.506;  label("Slope of log(WT/68 kg) on logit-F (unitless)")                                # Gandhi 2021 Table 2: 'Power of body weight on bioavailability' = -0.506
    e_age_vc   <-  0.114;  label("Power exponent of (AGE/49 yr) on VC (unitless)")                              # Gandhi 2021 Table 2: 'Power of age on VC' = 0.114
    e_age_f    <-  0.487;  label("Slope of log(AGE/49 yr) on logit-F (unitless)")                               # Gandhi 2021 Table 2: 'Power of age on bioavailability' = 0.487
    e_alb_cl   <- -0.722;  label("Power exponent of (ALB/4.1 g/dL) on CL (unitless)")                          # Gandhi 2021 Table 2: 'Power of albumin on CL' = -0.722
    e_crcl_cl  <-  0.259;  label("Power exponent of (CRCL/99.18 mL/min/1.73m^2) on CL (unitless)")             # Gandhi 2021 Table 2: 'Power of GFR on CL' = 0.259
    e_swol_cl  <-  0.0742; label("Power exponent of ((SWOL_28JOINT+1)/(15+1)) on CL (unitless)")               # Gandhi 2021 Table 2: 'Power of SJC on CL' = 0.0742
    e_sexf_cl  <-  0.0674; label("Exponential coefficient on CL for female sex (SEXF=1; unitless)")             # Gandhi 2021 Table 2: 'Exponent of male sex on CL' = 0.0674 (male = reference per Methods)
    e_nsaid_cl <-  0.102;  label("Exponential coefficient on CL for concomitant NSAID use (unitless)")          # Gandhi 2021 Table 2: 'Exponent of NSAID on CL' = 0.102
    e_jia_f    <-  3.08;   label("Additive coefficient on logit-F for pJIA (DIS_PJIA=1; unitless)")             # Gandhi 2021 Table 2: 'Exponent of JIA on bioavailability' = 3.08

    # Inter-individual variability - Gandhi 2021 Table 2 reports IIV variances
    # in the 'Estimate' column (the same column carries the SIGMA variances
    # for the residual error). The paper does not report a full block; IIVs
    # are treated as independent here, matching the Table 2 layout where each
    # ETA is reported with its own variance and shrinkage.
    etalka         ~ 1.11    ; label("IIV variance on the relative log-KA (Gandhi 2021 Table 2; shrinkage 74.9%)")
    etalvc         ~ 0.0464  ; label("IIV variance on log-VC (Gandhi 2021 Table 2; shrinkage 61.5%)")
    etalcl         ~ 0.0637  ; label("IIV variance on log-CL (Gandhi 2021 Table 2; shrinkage 14.3%)")
    etalvp         ~ 0.154   ; label("IIV variance on log-VP (Gandhi 2021 Table 2; shrinkage 54.7%)")
    etalogitfdepot ~ 0.516   ; label("IIV variance on logit-F (Gandhi 2021 Table 2; shrinkage 49.2%)")

    # Residual error - Gandhi 2021 Table 2 reports SIGMA matrix variances in
    # the same 'Estimate' column as the IIV variances (operator decision; see
    # vignette Errata for the rationale). The proportional and additive
    # variances are converted to standard deviations by sqrt() for nlmixr2's
    # add() / prop() conventions.
    propSd <- sqrt(0.0615); label("Proportional residual error (fraction)")                # Gandhi 2021 Table 2: SIGMA_PROP = 0.0615 (variance); SD = 0.248
    addSd  <- sqrt(0.00134); label("Additive residual error (mg/L = ug/mL)")                # Gandhi 2021 Table 2: SIGMA_ADD = 0.00134 (variance, mg^2/L^2); SD = 0.0366 mg/L
  })

  model({
    # Individual PK parameters with Gandhi 2021 final-model covariate equations
    # (Table 2). Reference subject: 49-year-old male, BWT 68 kg, ALB 4.1 g/dL,
    # cGFR 99.18 mL/min/1.73 m^2, swollen joint count 15, not on NSAIDs, RA
    # (DIS_PJIA = 0).
    cl <- exp(lcl + etalcl) *
          (WT / 68)^e_wt_cl *
          (ALB / 4.1)^e_alb_cl *
          (CRCL / 99.18)^e_crcl_cl *
          ((SWOL_28JOINT + 1) / (15 + 1))^e_swol_cl *
          exp(SEXF * e_sexf_cl + CONMED_NSAID * e_nsaid_cl)
    vc <- exp(lvc + etalvc) * (WT / 68)^e_wt_vc * (AGE / 49)^e_age_vc
    q  <- exp(lq)
    vp <- exp(lvp + etalvp) * (WT / 68)^e_wt_vp

    # KA is parameterised per Gandhi 2021 Methods (Equation S2) to ensure
    # KA > k_el (flip-flop prevention): KA_i = KA_TV * exp(etaKA) + CL_i / VC_i.
    kel <- cl / vc
    ka  <- exp(lka + etalka) + kel

    # Logit-scale SC bioavailability with independent IIV. Disease (pJIA),
    # body weight, and age covariates are added on the logit scale per Gandhi
    # 2021 Methods (Supplementary Equations S1 and S3). F_abs is bounded in
    # (0, 1) via the inverse logit.
    logit_f <- logitfdepot + etalogitfdepot +
               DIS_PJIA * e_jia_f +
               e_wt_f  * log(WT  / 68) +
               e_age_f * log(AGE / 49)
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
