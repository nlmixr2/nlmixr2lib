Leger_2004_topotecan <- function() {
  description <- "Two-compartment population PK model for oral and intravenous topotecan in adult cancer patients, with first-order absorption + lag time for the oral route, additive linear creatinine-clearance plus linear-ordinal WHO performance-status effects on CL, and linear body-weight effect on the central volume of distribution (Leger 2004)"
  reference <- "Leger F, Loos WJ, Fourcade J, Bugat R, Goffinet M, Mathijssen RHJ, Verweij J, Sparreboom A, Chatelut E. Factors affecting pharmacokinetic variability of oral topotecan: a population analysis. Br J Cancer. 2004;90(2):343-347. doi:10.1038/sj.bjc.6601469"
  vignette <- "Leger_2004_topotecan"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear (not allometric) scaling on central volume: V1 = e_wt_vc * WT (Leger 2004 Table 3 final covariate model: V1 = 0.58 * body weight). Cohort mean 70 kg (range 42-117 kg) per Leger 2004 Table 1.",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Computed via the Cockcroft-Gault equation in raw mL/min (NOT BSA-normalized to mL/min/1.73 m^2). Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form). Cohort mean 80 mL/min (range 33-167) per Leger 2004 Table 1. In model() the value is converted from mL/min to L/h via the constant 0.06 (= 60 min/h divided by 1000 mL/L) because Leger 2004 Table 3 footnote states the structural formula uses CrCl in L/h.",
      source_name        = "CRCL"
    ),
    WHO_PS = list(
      description        = "World Health Organization performance status (integer score)",
      units              = "(integer score 0-3 in this cohort)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-ordinal effect on CL: typical CL is multiplied by (1 - 0.12 * WHO_PS), so a WHO_PS = 2 patient has 24% lower CL than a WHO_PS = 0 patient (Leger 2004 Discussion). Distribution in cohort (Leger 2004 Table 1): WHO PS 0 / 1 / 2 / 3 = 79 / 98 / 11 / 2 patients (N = 190).",
      source_name        = "PS"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 190L,
    n_studies      = 6L,
    age_range      = "18-76 years",
    age_median     = "55 years",
    weight_range   = "42-117 kg",
    weight_median  = "70 kg",
    bsa_range      = "1.36-2.44 m^2 (DuBois formula)",
    bsa_median     = "1.80 m^2",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported (French and Dutch oncology centres)",
    disease_state  = "Adult cancer patients (predominantly ovarian cancer and other solid tumours; some studies combined with cisplatin)",
    dose_range     = "Oral 0.15-2.7 mg/m^2/day or IV infusion over 30 min 0.2-2.4 mg/m^2/day for 5-21 consecutive days",
    regions        = "France (Toulouse, Institut Claudius-Regaud) and the Netherlands (Rotterdam, Erasmus MC - Daniel den Hoed Cancer Center)",
    renal_function = "Cockcroft-Gault creatinine clearance mean 80 mL/min (range 33-167); serum creatinine mean 87 umol/L (range 41-162)",
    performance_status = "WHO PS distribution 0 / 1 / 2 / 3 = 79 / 98 / 11 / 2",
    notes          = "Baseline demographics per Leger 2004 Table 1. Pooled data from five separate clinical trials (173 patients) plus 17 patients with drug monitoring for other reasons, total 190 patients. IV (n=72) administered over 30-min infusions; oral (n=118) administered as gelatin capsules. Total 2064 plasma samples (median 15 per patient in cycle 1)."
  )

  ini({
    # Structural parameters at a typical patient (WHO PS 0, CrCl 80 mL/min = 4.8 L/h,
    # WT 70 kg). Leger 2004 Table 3 final covariate model (cycle 1 fit):
    #   CL (L/h) = (theta1 + theta2 * CrCl_Lh) * (1 - theta3 * WHO_PS)
    #   V1 (L)   = theta4 * body_weight
    #   V2 (L)   = theta5
    #   Q  (L/h) = theta6
    #   F  (--)  = theta7
    #   Ka (1/h) = theta9
    #   Tlag (h) = theta10
    lcl     <- log(12.8); label("Non-renal / intercept component of CL (L/h)")  # Leger 2004 Table 3: theta1 = 12.8 (95% CI 4.8); CL intercept of the additive linear renal-plus-non-renal model
    lvp     <- log(45.5); label("Peripheral volume of distribution V2 (L)")     # Leger 2004 Table 3: theta5 = 45.5 (95% CI 7.0)
    lq      <- log(49.2); label("Intercompartmental clearance Q (L/h)")         # Leger 2004 Table 3: theta6 = 49.2 (95% CI 16.9)
    lfdepot <- log(0.324); label("Oral bioavailability F (fraction)")           # Leger 2004 Table 3: theta7 = 32.4% (95% CI 3.9)
    lka     <- log(1.7);  label("First-order absorption rate constant Ka (1/h)") # Leger 2004 Table 3: theta9 = 1.7 (95% CI 0.6)
    ltlag   <- log(0.17); label("Oral absorption lag time (h)")                  # Leger 2004 Table 3: theta10 = 0.17 (95% CI 0.03)

    # Linear scaling slope on V1 (no allometric exponent: paper fits V1 = 0.58 * WT
    # as the final covariate form). lvc is log-transformed to keep the slope
    # positive; units are L/kg.
    lvc     <- log(0.58); label("Slope of V1 on body weight (L/kg)")             # Leger 2004 Table 3: theta4 = 0.58 (95% CI 0.13); central volume V1 = theta4 * body weight

    # Covariate-effect slopes on CL. Both are linear-scale because they appear
    # inside additive / one-minus-product covariate terms rather than as
    # exponentiated multipliers.
    e_crcl_cl   <- 2.1;  label("Slope of CL on Cockcroft-Gault CrCl (L/h per L/h of CrCl)") # Leger 2004 Table 3: theta2 = 2.1 (95% CI 1.0); additive linear CrCl-on-CL slope (CrCl in L/h)
    e_who_ps_cl <- 0.12; label("Fractional reduction of CL per unit WHO PS (1/unit)")        # Leger 2004 Table 3: theta3 = 0.12 (95% CI 0.09); multiplicative-linear PS effect (1 - 0.12 * WHO_PS)

    # Inter-individual variability (Leger 2004 Table 3 final-covariate-model
    # %CV column, cycle 1). Log-normal scale: omega^2 = log(CV^2 + 1).
    # Note: Leger 2004 also reports "interday variability" of CL (18%), V1 (49%),
    # and F (28%) as a separate interoccasion variability term; the standard
    # popPK extraction includes only the interindividual term and leaves IOV
    # out (it would require an OCC column in the event dataset). The
    # interindividual values below pin the published CV% values directly.
    etalcl     ~ 0.08075 # log(0.29^2 + 1); 29% CV on CL (Leger 2004 Table 3)
    etalvc     ~ 0.14157 # log(0.39^2 + 1); 39% CV on V1
    etalvp     ~ 0.19193 # log(0.46^2 + 1); 46% CV on V2
    etalq      ~ 0.49470 # log(0.80^2 + 1); 80% CV on Q
    etalfdepot ~ 0.04727 # log(0.22^2 + 1); 22% CV on F (oral only)
    etalka     ~ 0.34335 # log(0.64^2 + 1); 64% CV on Ka (oral only)
    etaltlag   ~ 0.01676 # log(0.13^2 + 1); 13% CV on Tlag (oral only)

    # Combined (additive + proportional) residual error. Leger 2004 Results
    # reports route-specific residual error values:
    #     IV   data: proportional 11% , additive 0.64 ug/L
    #     Oral data: proportional 17% , additive 0.09 ug/L
    # The structural difference between routes is attributable to different
    # HPLC lower limits of quantification at the two analytical sites
    # (Toulouse 0.5 ng/mL for IV, Rotterdam 0.1 ng/mL for oral). The encoded
    # model uses the oral-route residual error because (a) the paper's primary
    # contribution is oral topotecan PK and the limited-sampling strategy is
    # developed for oral, and (b) ~62% of the cohort (118 / 190) is oral. The
    # IV-route alternative values are documented in the vignette Assumptions
    # and deviations section so a user simulating IV-only data can swap them
    # in.
    propSd <- 0.17; label("Proportional residual error (fraction)") # Leger 2004 Results: 17% proportional residual on oral data
    addSd  <- 0.09; label("Additive residual error (ug/L)")          # Leger 2004 Results: 0.09 ug/L additive residual on oral data
  })

  model({
    # Convert raw Cockcroft-Gault CrCl from mL/min (covariate column) to L/h
    # for use in the structural CL formula. The conversion constant 0.06 =
    # 60 min/h / 1000 mL/L.
    crcl_Lh <- CRCL * 0.06

    # Typical-value clearance and individual log-normal IIV. Structural form
    # from Leger 2004 Table 3:
    #   CL = (intercept + slope * CRCL_Lh) * (1 - e_who_ps_cl * WHO_PS)
    cl <- (exp(lcl) + e_crcl_cl * crcl_Lh) * (1 - e_who_ps_cl * WHO_PS) * exp(etalcl)

    # Central volume V1 = slope * WT (linear scaling, no allometric exponent).
    vc <- exp(lvc + etalvc) * WT

    # Peripheral volume and intercompartmental clearance -- constant across
    # patients in the final covariate model.
    vp <- exp(lvp + etalvp)
    q  <- exp(lq  + etalq)

    # Absorption parameters (oral route).
    ka   <- exp(lka   + etalka)
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                 k12 * central - k21 * peripheral1

    # Apply oral bioavailability and lag-time to the depot compartment. IV
    # doses bypass the depot (cmt = central in the event table) and are not
    # affected by f(depot) or alag(depot).
    f(depot)    <- exp(lfdepot + etalfdepot)
    alag(depot) <- tlag

    # Dose in mg, volumes in L, concentration converted to ug/L via 1000.
    # central has units mg, so central / vc has units mg/L = ug/uL; multiply
    # by 1000 to match the paper's reported concentration units (ug/L).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
