Wahlby_2004_gentamicin <- function() {
  description <- "Two-compartment population PK model for intravenous gentamicin in 210 cancer patients, demonstrating Wahlby 2004's extended covariate-model formulation in which the within- and between-subject components of a time-varying covariate are entered as separate model terms. Final-model clearance depends on baseline creatinine clearance (CRCL_BASE) and the time-varying delta-from-baseline (CRCL - CRCL_BASE); central volume depends on baseline body surface area (BSA_BASE) and the time-varying albumin (ALB) ratio (ALB/34)^-0.41. Underlying structural PK (compartment count, IIV/residual-error structure) follows Rosario et al. 1998 (Br J Clin Pharmacol 46:229-236)."
  reference <- paste(
    "Wahlby U, Thomson AH, Milligan PA, Karlsson MO.",
    "Models for time-varying covariates in population pharmacokinetic-pharmacodynamic analysis.",
    "Br J Clin Pharmacol 2004;58(4):367-377.",
    "doi:10.1111/j.1365-2125.2004.02170.x.",
    "Structural PK model carried from Rosario MC, Thomson AH, Jodrell DI, Sharp CA, Elliott HL.",
    "Population pharmacokinetics of gentamicin in patients with cancer.",
    "Br J Clin Pharmacol 1998;46(3):229-236; the present entry encodes the Wahlby 2004 final-model parameter estimates (Table 5, Final-Model column).",
    sep = " "
  )
  vignette <- "Wahlby_2004_time_varying_covariates"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance, time-varying within an individual (raw mL/min, NOT BSA-normalized; reflects Wahlby 2004's CLC measurement convention).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Wahlby 2004 reports CLC as raw mL/min (Table 1). The canonical CRCL register entry is for BSA-normalized mL/min/1.73 m^2; this model preserves the source unit (mL/min) and documents it here. Population median (overall) = 76.9 mL/min, mean = 81.9, range 14.8-183 (Table 1). Centering value in Eq 4: 81 mL/min.",
      source_name        = "CLC"
    ),
    CRCL_BASE = list(
      description        = "Per-subject baseline creatinine clearance, time-fixed (= CRCL value at the first observation per subject).",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject constant across time. Median across 210 subjects = 71.7 mL/min (Table 1, BCLC row). Centering value for the BCLC effect on CL is the BCLC median (71.7 mL/min). The Wahlby 2004 final-model paper estimate theta_BCLC = 0.0098 per mL/min was reported relative to centering at the BCLC median per their Methods description ('BCOV - BCOV_median').",
      source_name        = "BCLC"
    ),
    BSA_BASE = list(
      description        = "Per-subject baseline body surface area, time-fixed (BBSA in the source paper).",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject constant. Median across 210 subjects = 1.7 m^2 (Table 1, BBSA row). Replaces time-varying BSA in the final-model V1 equation because BBSA was the better predictor (delta-BSA was imprecisely estimated, RSE 110%, per Results section).",
      source_name        = "BBSA"
    ),
    ALB = list(
      description        = "Serum albumin, time-varying within an individual.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Wahlby 2004 reports albumin in g/L (Table 1: mean 33.5, median 34.0, range 14-53). Centering value in the V1 model: ALB_ref = 34 g/L (Table 5 footnote a). The final model uses the time-varying ALB directly (not BALB / DALB) because the BCOV/DCOV split for ALB did not improve fit (Results section).",
      source_name        = "ALB"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 210L,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    disease_state  = "Adult cancer patients receiving intravenous gentamicin. Wahlby 2004 re-analyses the gentamicin cohort originally reported by Rosario MC, Thomson AH, Jodrell DI, Sharp CA, Elliott HL (Br J Clin Pharmacol 1998;46(3):229-236).",
    dose_range     = "Intravenous gentamicin; specific dose ranges and administration details are reported by Rosario 1998 and are not reproduced in Wahlby 2004.",
    n_observations = 574L,
    n_covariate_measurements_clc = 576L,
    n_covariate_measurements_alb = 439L,
    n_covariate_measurements_bsa = 248L,
    follow_up      = "Patients followed for one to five courses (median one), over a maximum of 25 months; on average data were available for the first 2 days per course, with 1-9 samples per patient per occasion.",
    regions        = NA_character_,
    notes          = "Demographic detail (age, weight, sex, race) is not reported in Wahlby 2004; consult Rosario 1998 for the underlying cohort description."
  )

  ini({
    # Structural typical-value parameters (Wahlby 2004 Table 5, Final-Model column)
    lcl <- log(4.1);  label("Population typical CL (L/h)")                         # Table 5: theta_CL = 4.1 L/h (RSE 2.3%)
    lvc <- log(8.63); label("Population typical V1 per unit BSA_BASE (L/m^2)")     # Table 5: theta_V1/BBSA = 8.63 L/m^2 (RSE 2.6%); enters as vc = 8.63 * BSA_BASE * (ALB/34)^e_alb_vc
    lq  <- log(1.21); label("Intercompartmental clearance Q (L/h)")                # Table 5: theta_Q = 1.21 L/h (RSE 19%)
    lvp <- log(8.12); label("Peripheral volume V2 (L)")                            # Table 5: theta_V2 = 8.12 L (RSE 16%)

    # Covariate effects on CL (Wahlby 2004 Eq 2 form with BCOV/DCOV split; Table 5 Final-Model column)
    e_bcrcl_cl <- 0.0098; label("Fractional change in CL per mL/min deviation of CRCL_BASE from 71.7") # Table 5: theta_BCLC = 0.0098 (RSE 5.5%); centered at BCLC median = 71.7
    e_dcrcl_cl <- 0.0068; label("Fractional change in CL per mL/min of (CRCL - CRCL_BASE)")            # Table 5: theta_DCLC = 0.0068 (RSE 23%)

    # Covariate effects on V1 (Eq 5 with final-model encoding)
    e_alb_vc <- -0.41;    label("Power exponent of (ALB / 34) on V1 (unitless)")     # Table 5 footnote a: V1 ~ [ALB/34]^theta_ALB, point estimate -0.41 (RSE 26%)

    # Inter-individual variability (Wahlby 2004 Table 5 Final-Model column, omega^2 form)
    # Methods note exponential parameter-variability models; reported omegas are SDs of the log-eta.
    # omega_CL = 0.22 (RSE 15%, relative to omega^2), omega_Q = 0.61 (RSE 56%, relative to omega^2).
    etalcl ~ 0.22^2  # Table 5: omega_CL = 0.22; variance = 0.0484
    etalq  ~ 0.61^2  # Table 5: omega_Q  = 0.61; variance = 0.3721

    # Residual error (Wahlby 2004 Table 5 Final-Model column)
    # Methods note "proportional and/or additive model for the residual variability".
    # sigma_prop = 0.17 (fraction; RSE 16% relative to sigma^2)
    # sigma_add  = 0.22 (mg/L; RSE 40% relative to sigma^2)
    propSd <- 0.17; label("Proportional residual error (fraction)")  # Table 5: sigma_prop = 0.17
    addSd  <- 0.22; label("Additive residual error (mg/L)")          # Table 5: sigma_add  = 0.22
  })

  model({
    # Delta-from-baseline derivation (the Wahlby 2004 DCOV term computed inline from
    # the user-supplied time-varying CRCL and per-subject CRCL_BASE columns).
    dcrcl <- CRCL - CRCL_BASE

    # Individual parameters with the Wahlby 2004 extended-covariate form on CL
    # (Eq 2): P = theta_p * [1 + theta_BCOV * (BCOV - BCOV_median) + theta_DCOV * DCOV].
    cl <- exp(lcl + etalcl) * (1 + e_bcrcl_cl * (CRCL_BASE - 71.7) + e_dcrcl_cl * dcrcl)

    # Central volume of distribution: V1 = (V1/BBSA) * BBSA * (ALB/34)^e_alb_vc per
    # Wahlby 2004 Eq 5 with BSA replaced by BSA_BASE in the final-model encoding
    # (Table 5 Final-Model column). Table 5 footnote a documents the power form.
    vc <- exp(lvc) * BSA_BASE * (ALB / 34)^e_alb_vc

    q  <- exp(lq + etalq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV gentamicin PK; doses go directly to central. The library
    # model does not hard-code an infusion duration; users specify rate / dur per
    # dose in their event table.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
