Narwal_2013_sifalimumab <- function() {
  description <- "Two-compartment population PK model for sifalimumab (anti-IFN-alpha IgG1) in adult patients with systemic lupus erythematosus (Narwal 2013)"
  reference <- "Narwal R, Roskos LK, Robbie GJ. Population pharmacokinetics of sifalimumab, an investigational anti-interferon-alpha monoclonal antibody, in systemic lupus erythematosus. Clin Pharmacokinet. 2013;52(11):1017-1027. doi:10.1007/s40262-013-0085-2"
  vignette <- "Narwal_2013_sifalimumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL, V1 (central volume), and V2 (peripheral volume) with reference weight 75 kg; exponents 0.481, 0.489, and 0.646 respectively. Narwal 2013 Table 2 / Eqs. 3-5.",
      source_name        = "WT"
    ),
    BGENE21 = list(
      description        = "Baseline type I interferon gene signature from 21 IFN-inducible genes",
      units              = "(gene-signature score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL with reference value 32 and exponent 0.0558 (Narwal 2013 Eq. 3). Population median was 33 (range 0.63-87; Narwal 2013 Table 1).",
      source_name        = "BGENE21"
    ),
    COHDOSE = list(
      description        = "Randomized dose cohort, expressed in mg/kg",
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL with reference value 1 mg/kg and exponent 0.0542 (Narwal 2013 Eq. 3). The MI-CP152 study fixed each subject to one of 0.3, 1, 3, or 10 mg/kg Q14D across all infusions, so this is a subject-level covariate rather than a per-dose covariate. Narwal et al. acknowledge (Discussion, p. 1024) that the apparent dose effect may be a data artifact of the escalating cohort design since single-dose data in MI-CP126 were linear across 0.3-30 mg/kg.",
      source_name        = "DOSE"
    ),
    STEROID = list(
      description        = "Baseline systemic corticosteroid use (0 = no, 1 = yes)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no baseline corticosteroid use)",
      notes              = "Multiplicative effect on CL of the form (1 + 0.195 x STEROID), giving +19.5% CL for baseline-steroid users (Narwal 2013 Eq. 3, Table 2). Corticosteroids (e.g., methylprednisolone) were used throughout the study to control disease activity.",
      source_name        = "BSTEROID"
    )
  )

  population <- list(
    n_subjects     = 120L,
    n_studies      = 1L,
    study_id       = "MI-CP152 (NCT00482989)",
    age_range      = "18-71 years",
    age_median     = "43 years",
    weight_range   = "43.1-120 kg",
    weight_median  = "73 kg",
    sex_female_pct = 95,
    race_ethnicity = "Not tabulated in the main text",
    regions        = c("North America" = 71, "South America" = 29),
    disease_state  = "Moderately active adult systemic lupus erythematosus (SLE)",
    dose_range     = "0.3, 1, 3, or 10 mg/kg IV every 14 days over 30-60 min, up to 14 infusions",
    bsledai_median = "10 (range 2-34)",
    bgene21_median = "33 (range 0.63-87)",
    notes          = "Phase Ib MI-CP152 study (Narwal 2013 Section 3.1 / Table 1). A total of 2,370 serum concentrations from 120 evaluable subjects (one 10 mg/kg subject was excluded due to very low observed concentrations). Eight subjects had BGENE4 imputed to the population median because the 4-gene signature was unavailable; BGENE21 was available for all 120."
  )

  ini({
    # Structural parameters - typical values for a 75 kg subject with BGENE21 = 32,
    # DOSE = 1 mg/kg, and no baseline steroid use.
    lcl <- log(0.176);  label("Clearance CL for the reference subject (L/day)")                           # Narwal 2013 Table 2: theta1 = 0.176 L/day (176 mL/day)
    lvc <- log(2.90);   label("Central volume V1 for a 75 kg subject (L)")                                # Narwal 2013 Table 2: theta2 = 2.90 L
    lvp <- log(2.12);   label("Peripheral volume V2 for a 75 kg subject (L)")                             # Narwal 2013 Table 2: theta3 = 2.12 L
    lq  <- log(0.171);  label("Intercompartmental clearance Q (L/day)")                                   # Narwal 2013 Table 2: theta4 = 0.171 L/day

    # Covariate effects on CL (Eq. 3)
    e_wt_cl       <- 0.481;  label("Body-weight exponent on CL (unitless)")                               # Narwal 2013 Table 2: theta5
    e_bgene21_cl  <- 0.0558; label("BGENE21 exponent on CL (unitless)")                                   # Narwal 2013 Table 2: theta6
    e_cohdose_cl  <- 0.0542; label("Dose (mg/kg) exponent on CL (unitless)")                              # Narwal 2013 Table 2: theta7
    e_steroid_cl  <- 0.195;  label("Fractional change in CL for baseline steroid users (unitless)")       # Narwal 2013 Table 2: theta8

    # Covariate effects on volumes (Eqs. 4-5)
    e_wt_v1       <- 0.489;  label("Body-weight exponent on V1 (unitless)")                               # Narwal 2013 Table 2: theta9
    e_wt_v2       <- 0.646;  label("Body-weight exponent on V2 (unitless)")                               # Narwal 2013 Table 2: theta10

    # Inter-individual variability. CL, V1, V2 form a 3x3 log-normal block with
    # correlations CL-V1 = 0.557 and V1-V2 = 0.131 (Narwal 2013 Table 2).
    # The CL-V2 correlation is not reported and is taken as 0 (not included in
    # the NONMEM OMEGA block). Variances from %CV via omega^2 = log(1 + CV^2):
    #   omega^2_CL = log(1 + 0.28^2) = 0.075478
    #   omega^2_V1 = log(1 + 0.31^2) = 0.091758
    #   omega^2_V2 = log(1 + 0.58^2) = 0.289979
    # Covariances from r * sqrt(var1 * var2):
    #   cov(CL, V1) = 0.557 * sqrt(0.075478 * 0.091758) = 0.046354
    #   cov(V1, V2) = 0.131 * sqrt(0.091758 * 0.289979) = 0.021369
    etalcl + etalvc + etalvp ~ c(
      0.075478,
      0.046354, 0.091758,
      0.000000, 0.021369, 0.289979
    )
    # Q IIV is uncorrelated with the others (no correlation reported in Table 2).
    # omega^2_Q = log(1 + 0.71^2) = 0.408195
    etalq ~ 0.408195                                                                                     # Narwal 2013 Table 2: Q %CV = 71%

    # Residual error (proportional, 27.5% CV).
    propSd <- 0.275; label("Proportional residual error (fraction)")                                     # Narwal 2013 Table 2: residual error 27.5%
  })
  model({
    # Covariate-adjusted individual PK parameters (Eqs. 3-6).
    cl <- exp(lcl + etalcl) *
      (WT / 75)^e_wt_cl *
      (BGENE21 / 32)^e_bgene21_cl *
      (COHDOSE / 1)^e_cohdose_cl *
      (1 + e_steroid_cl * STEROID)
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_v1
    vp <- exp(lvp + etalvp) * (WT / 75)^e_wt_v2
    q  <- exp(lq  + etalq)

    # Two-compartment micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg and vc in L, so central/vc has units mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
