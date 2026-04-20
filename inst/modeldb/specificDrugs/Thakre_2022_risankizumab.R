Thakre_2022_risankizumab <- function() {
  description <- "Two-compartment population PK model of risankizumab (anti-IL-23 mAb) with first-order SC absorption in patients with active psoriatic arthritis (Thakre 2022)"
  reference <- "Thakre N, D'Cunha R, Goebel A, Liu W, Pang Y, Suleiman AA. Population Pharmacokinetics and Exposure-Response Analyses for Risankizumab in Patients with Active Psoriatic Arthritis. Rheumatol Ther. 2022;9(6):1587-1603. doi:10.1007/s40744-022-00495-0"
  vignette <- "Thakre_2022_risankizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL and Vc; normalized as WT/70 per Thakre 2022 Eq. 2 and Eq. 3 (reference: 70 kg). Source column 'WTKG' renamed to canonical WT.",
      source_name        = "WTKG"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL; normalized as AGE/52 per Thakre 2022 Eq. 2 (reference: 52 years).",
      source_name        = "AGE"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL; normalized as ALB/45 per Thakre 2022 Eq. 2 (reference: 45 g/L). Source uses SI units (g/L); convert g/dL to g/L by x10 if needed.",
      source_name        = "ALB"
    ),
    CREAT = list(
      description        = "Baseline serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL; normalized as CREAT/70.73 per Thakre 2022 Eq. 2 (reference: 70.73 umol/L, equivalent to 0.8 mg/dL). Source column 'CRE' renamed to canonical CREAT.",
      source_name        = "CRE"
    ),
    CRP = list(
      description        = "Baseline high-sensitivity C-reactive protein (hs-CRP assay; baseline, time-fixed per subject)",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL; normalized as CRP/5.21 per Thakre 2022 Eq. 2 (reference: 5.21 mg/L). Source column 'CRPHS' (high-sensitivity CRP) maps to the canonical general-scope CRP covariate; the assay type (hs-CRP) is documented here in the covariateData entry rather than via a separate hsCRP canonical.",
      source_name        = "CRPHS"
    )
  )

  population <- list(
    n_subjects     = 1527L,
    n_studies      = 5L,
    age_range      = "Adults, median reference ~52 years (range not published in main paper)",
    weight_range   = "Reference 70 kg; PsA phase-3 populations typically 40-160 kg",
    sex_female_pct = NA_real_,
    race_ethnicity = "Predominantly White, as typical in global PsA phase-2/3 trials (details in source supplement).",
    disease_state  = "Active psoriatic arthritis (PsA), pooled across one phase 1 healthy-participant study, one phase 2 dose-ranging study with open-label extension, and two pivotal phase 3 studies (KEEPsAKE 1 [NCT03675308] and KEEPsAKE 2 [NCT03671148]).",
    dose_range     = "Risankizumab 18-300 mg SC and 0.01-5 mg/kg IV single-dose (phase 1), plus 150 mg SC at weeks 0 and 4 and every 12 weeks thereafter (phase 2/3 clinical regimen).",
    regions        = "Global (multi-regional).",
    notes          = "Data set: 3631 concentration measurements from 1527 individuals (study 1: n=67, studies 2/3: n=177, study 4: n=391, study 5: n=892). ADA titers on CL and WT on V2 were evaluated but did not meet retention criteria and are not in the final model. Reference covariate values (used as the normalizers in the covariate power terms): WT = 70 kg, AGE = 52 years, ALB = 45 g/L, CREAT = 70.73 umol/L, CRP = 5.21 mg/L (hs-CRP assay)."
  )

  ini({
    # Structural parameters from Thakre 2022 Table 1 (final PsA population PK model).
    # Typical values are for the reference subject (70 kg, age 52, ALB 45 g/L,
    # CREAT 70.73 umol/L, CRP 5.21 mg/L [hs-CRP]).
    lcl     <- log(0.248); label("Clearance (CL, L/day)")                         # Thakre 2022 Table 1
    lvc     <- log(4.71);  label("Central volume of distribution (Vc / V1, L)")   # Thakre 2022 Table 1
    lq      <- log(0.839); label("Intercompartmental clearance (Q, L/day)")       # Thakre 2022 Table 1
    lvp     <- log(4.26);  label("Peripheral volume of distribution (Vp / V2, L)")# Thakre 2022 Table 1
    lka     <- log(0.218); label("First-order SC absorption rate (ka, 1/day)")    # Thakre 2022 Table 1
    lfdepot <- log(0.835); label("Absolute SC bioavailability (F)")                # Thakre 2022 Table 1

    # Covariate effects from Thakre 2022 Table 1 and Eq. 2 / Eq. 3. All are
    # power exponents on covariates normalized to their reference values.
    e_wt_cl    <-  0.869;  label("Power exponent of body weight on CL (unitless)")        # Thakre 2022 Table 1
    e_alb_cl   <- -0.703;  label("Power exponent of serum albumin on CL (unitless)")      # Thakre 2022 Table 1
    e_creat_cl <- -0.201;  label("Power exponent of serum creatinine on CL (unitless)")   # Thakre 2022 Table 1
    e_crp_cl   <-  0.0471; label("Power exponent of CRP on CL (unitless)")                # Thakre 2022 Table 1
    e_age_cl   <- -0.138;  label("Power exponent of age on CL (unitless)")                # Thakre 2022 Table 1
    e_wt_vc    <-  1.46;   label("Power exponent of body weight on Vc (unitless)")        # Thakre 2022 Table 1

    # Inter-individual variability. Table 1 reports variances directly
    # (omega^2 on the log scale) and the implied %CV = sqrt(exp(omega^2) - 1) * 100.
    # CL and V1 have a reported covariance of 0.0836 (correlation ~65.8%).
    etalcl + etalvc ~ c(0.0943,
                        0.0836, 0.171) # Thakre 2022 Table 1 (variance CL, covariance CL:V1, variance V1)
    etalka ~ 0.164                     # Thakre 2022 Table 1 (variance on ka)

    # Residual error. Table 1 reports the variance of the proportional residual
    # error (0.0382); nlmixr2 prop(propSd) expects an SD, so store sqrt(variance).
    propSd <- sqrt(0.0382); label("Proportional residual error (SD, fraction)")   # Thakre 2022 Table 1
  })
  model({
    # Individual PK parameters. Reference subject: 70 kg, age 52 years,
    # ALB 45 g/L, CREAT 70.73 umol/L, CRP 5.21 mg/L [hs-CRP]. Covariate forms per
    # Thakre 2022 Eq. 2 (CL) and Eq. 3 (V1) are power models normalized to
    # each covariate's reference value.
    cl <- exp(lcl + etalcl) *
      (WT    / 70)^e_wt_cl *
      (ALB   / 45)^e_alb_cl *
      (CREAT / 70.73)^e_creat_cl *
      (CRP   / 5.21)^e_crp_cl *
      (AGE   / 52)^e_age_cl

    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    vp <- exp(lvp)
    q  <- exp(lq)
    ka <- exp(lka + etalka)

    # Two-compartment model with first-order SC absorption.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                 k12 * central - k21 * peripheral1

    # Absolute SC bioavailability applied to the depot dose.
    f(depot) <- exp(lfdepot)

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL.
    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
