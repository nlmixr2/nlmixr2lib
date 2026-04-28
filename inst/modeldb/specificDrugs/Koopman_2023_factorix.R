Koopman_2023_factorix <- function() {
  description <- "Two-compartment population PK model for recombinant factor IX-Fc fusion concentrate (rFIX-Fc, eftrenonacog alfa) in haemophilia B patients aged 2-71 years (Koopman 2023)"
  reference <- "Koopman SF, Goedhart TMHJ, Bukkems LH, et al. A new population pharmacokinetic model for recombinant factor IX-Fc fusion concentrate including young children with haemophilia B. Br J Clin Pharmacol. 2024;90(1):220-231. doi:10.1111/bcp.15881"
  vignette <- "Koopman_2023_factorix"
  units <- list(time = "hour", dosing = "IU", concentration = "IU/dL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL, Q, V1, V2 with reference weight 73 kg (Koopman 2023 Table 2 footnote and equations on p. 226). Exponents fixed at 0.75 for CL and Q, 1.00 for V1 and V2.",
      source_name        = "BW"
    ),
    AGE = list(
      description        = "Subject age in years",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation effect on CL centered at 15.8 years (median age of the development cohort). Effect form: CL multiplied by (1 - 0.0047 * (AGE - 15.8)) (Koopman 2023 equations on p. 226).",
      source_name        = "AGE"
    )
  )

  population <- list(
    n_subjects     = 37L,
    n_studies      = 2L,
    age_range      = "2-71 years",
    age_median     = "15.8 years (IQR 11-30)",
    weight_range   = "12-103 kg",
    weight_median  = "65.4 kg (IQR 33-77)",
    sex_female_pct = 0,
    race_ethnicity = "Not reported",
    disease_state  = "Severe (n = 35) or moderately severe (n = 2) haemophilia B; baseline FIX < 1 IU/dL in severe patients (median baseline FIX 1.0 IU/dL in non-severe).",
    dose_range     = "Median 36 IU/kg rFIX-Fc IV (range 10-132 IU/kg) prophylactically administered",
    regions        = "Netherlands (OPTI-CLOT TARGET, NTR7523) and United Kingdom (UK-EHL Outcomes Registry, NCT02938156)",
    paediatric_breakdown = "19 patients < 18 years (51%), 14 patients < 12 years (38%), 7 patients < 6 years (19%); paediatric median age 11 years (IQR 4-12, range 2-16); paediatric median weight 32.8 kg (IQR 16-32, range 12-52)",
    blood_samples  = "Median 5 FIX activity levels per PK profile (range 3-7); 287 measurements total, 3 below LLOQ excluded",
    notes          = "Real-world prospective + retrospective data pooled from OPTI-CLOT TARGET (NL) and UK-EHL Outcomes Registry; baseline demographics per Koopman 2023 Table 1. Hemophilia B is X-linked, so the cohort is essentially all male (sex_female_pct = 0)."
  )

  ini({
    # Structural parameters - typical values for the paper's reference patient
    # (BW = 73 kg, AGE = 15.8 years). Units: CL and Q in dL/h; V1 and V2 in dL.
    # Source: Koopman 2023 Table 2 ("New" column).
    lcl  <- log(1.41); label("Clearance for the reference 73 kg, 15.8-year-old patient (CL, dL/h)") # Koopman 2023 Table 2: CL = 1.41 dL/h
    lvc  <- log(73.1); label("Central volume of distribution for the reference 73 kg patient (V1, dL)") # Koopman 2023 Table 2: V1 = 73.1 dL
    lq   <- log(2.77); label("Intercompartmental clearance for the reference 73 kg patient (Q2, dL/h)") # Koopman 2023 Table 2: Q2 = 2.77 dL/h
    lvp  <- log(80.1); label("Peripheral volume of distribution for the reference 73 kg patient (V2, dL)") # Koopman 2023 Table 2: V2 = 80.1 dL

    # Covariate effect: linear-deviation effect of AGE on CL, centered at 15.8 years.
    # Form: CL multiplied by (1 - e_age_cl * (AGE - 15.8)). Sign convention: a positive
    # e_age_cl means CL decreases with increasing age above 15.8 years (Koopman 2023
    # Discussion: "typical clearance of a 73-kg patient would decrease from 1.89 dL/h
    # at age 20 years to 1.36 dL/h at 70 years").
    e_age_cl <- 0.0047; label("Linear age slope on CL, applied as (1 - e_age_cl * (AGE - 15.8)) (1/year)") # Koopman 2023 Table 2: age exponent on CL = 0.0047

    # Inter-individual variability. The paper footnote states
    # "IIV ... coefficient of variation calculated as: sqrt(variance) * 100%",
    # i.e., omega^2 = (CV/100)^2 directly (no log(CV^2 + 1) transform).
    # Correlation footnote: "covariance / (sqrt(variance1) * sqrt(variance2)) * 100%".
    # CL/V1 correlated block:
    #   omega^2_CL = 0.236^2 = 0.05570
    #   cov_CL_V1  = 0.44 * 0.236 * 0.316 = 0.03281
    #   omega^2_V1 = 0.316^2 = 0.09986
    etalcl + etalvc ~ c(0.05570,
                        0.03281, 0.09986)                   # Koopman 2023 Table 2: IIV CL = 23.6%, IIV V1 = 31.6%, corr CL:V1 = 44.0%
    # Independent IIV on V2: omega^2 = 0.412^2 = 0.16974
    etalvp ~ 0.16974                                        # Koopman 2023 Table 2: IIV V2 = 41.2%

    # Residual error: combined additive + proportional (Koopman 2023 Table 2)
    propSd <- 0.163; label("Proportional residual error (fraction)")        # Koopman 2023 Table 2: proportional error = 16.3%
    addSd  <- 1.04;  label("Additive residual error (IU/dL)")               # Koopman 2023 Table 2: additive error = 1.04 IU/dL
  })
  model({
    # Age effect on CL (linear-deviation, centered at 15.8 years).
    # Reference age 15.8 years matches the median age of the development cohort.
    age_cl <- 1 - e_age_cl * (AGE - 15.8)

    # Individual PK parameters with allometric weight scaling (reference 73 kg).
    # Allometric exponents fixed at theoretical values: 0.75 for CL/Q, 1.00 for V1/V2
    # (Koopman 2023 Table 2 and equations on p. 226). The paper fits IIV on CL, V1
    # and V2 only (no IIV on Q). Inter-occasion variability on CL was reported
    # (19.8%) but is omitted here because the static library model has no occasion
    # variable; it is documented in the vignette's Assumptions and deviations
    # section.
    cl <- exp(lcl + etalcl) * (WT / 73)^0.75 * age_cl
    vc <- exp(lvc + etalvc) * (WT / 73)^1.00
    q  <- exp(lq)           * (WT / 73)^0.75
    vp <- exp(lvp + etalvp) * (WT / 73)^1.00

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV model: dose enters central directly (rFIX-Fc is given IV).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # FIX activity: dose in IU and central volume in dL -> central / vc has units IU/dL.
    # The paper fits to baseline-corrected FIX activity (one-stage assay with endogenous
    # baseline and prior-product residual subtracted; Koopman 2023 Eqs. 1-3), so Cc here
    # represents the rFIX-Fc-attributable FIX activity above baseline.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
