Zhou_1996_ganciclovir <- function() {
  description <- paste(
    "One-compartment population PK model for IV ganciclovir in newborns with acute",
    "symptomatic congenital CMV disease (Zhou 1996), with additive-linear",
    "renal-function scaling on clearance, additive-linear body-weight scaling on",
    "the central volume, and correlated between-subject variability on CL and Vc.",
    "Parameters transcribed from the Yang 2023 ganciclovir / valganciclovir",
    "population-PK model repository (Table 3), not from the primary publication;",
    "re-verify against Zhou 1996 when the primary is obtained.",
    sep = " "
  )
  reference <- paste(
    "Zhou XJ, Gruber W, Demmler G, Jacobs R, Reuman P, Adler S, Shelton M, Pass R,",
    "Britt B, Trang JM, et al. Population pharmacokinetics of ganciclovir in",
    "newborns with congenital cytomegalovirus infections. NIAID Collaborative",
    "Antiviral Study Group. Antimicrob Agents Chemother. 1996;40(9):2202-2205.",
    "doi:10.1128/AAC.40.9.2202.",
    "Parameters transcribed from Yang W, Mak W, Gwee A, Gu M, Wu Y, Shi Y, He Q,",
    "Xiang X, Han B, Zhu X. Establishment and Evaluation of a Parametric Population",
    "Pharmacokinetic Model Repository for Ganciclovir and Valganciclovir.",
    "Pharmaceutics. 2023;15(7):1801. doi:10.3390/pharmaceutics15071801 (Table 3).",
    sep = " "
  )
  vignette <- "Yang_2023_ganciclovir_model_repository"
  units    <- list(time = "hr", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Additive-linear effect on the central volume with a non-zero intercept:",
        "Vc = 0.627 + 0.437 * BW (L), so the slope 0.437 L/kg is a volume increment",
        "per kilogram and 0.627 L is the weight-independent intercept. Weight is not",
        "a covariate on clearance in this model. Yang 2023 Table 2 reports the",
        "newborn cohort weights as not reported (NR); Yang 2023 Table 1 uses 3 kg",
        "for its typical neonatal virtual patient.",
        sep = " "
      ),
      source_name        = "BW"
    ),
    CRCL = list(
      description        = "Approximated creatinine clearance from serum creatinine (ASCC)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Yang 2023 Table 3 footnote defines ASCC as 'approximated creatinine",
        "clearance from serum (mL/min/1.73 m^2)', i.e. a BSA-normalized",
        "creatinine-based estimate of renal function, stored under the canonical",
        "CRCL column. Additive-linear effect on clearance with a non-zero intercept:",
        "CL = 0.262 + 0.00271 * ASCC (L/h), so the slope is 0.00271 L/h per",
        "mL/min/1.73 m^2 and no centering or divisive normalization is applied",
        "(the same additive-linear pattern used in Georges_2009_ceftazidime.R).",
        "ASCC was the only covariate retained on CL (Yang 2023 Table 4).",
        "Yang 2023 Table 1 uses 30 umol/L serum creatinine for its typical neonatal",
        "virtual patient.",
        sep = " "
      ),
      source_name        = "ASCC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 27L,
    n_studies      = 1L,
    n_observations = 219L,
    age_median     = "Newborns (exact ages not reported in Yang 2023 Table 2).",
    weight_median  = "Not reported in Yang 2023 Table 2.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported.",
    disease_state  = "Newborns with acute symptomatic congenital CMV disease.",
    dose_range     = "A single IV ganciclovir dose of 4 or 6 mg/kg as a 1 h constant-rate infusion.",
    regions        = "United States (prospective; NIAID Collaborative Antiviral Study Group).",
    bioassay       = "HPLC, LLOQ 0.1 ug/mL.",
    notes          = paste(
      "Demographics and dosing from Yang 2023 Table 2. Sex, age detail and weights",
      "are recorded as not reported (NR) in the review's table. Sampling strategy",
      "not reported. Covariates tested: weight, ASCC and platelet count; only ASCC",
      "was retained on CL and weight on Vc (Yang 2023 Table 4). The primary study",
      "reported no model-evaluation results (published before 2000), but Yang 2023",
      "Section 3.2.1 states that after careful inspection this model's performance",
      "was comparable to the others in the repository. The primary study's stated",
      "purpose was to evaluate the effect of covariates on ganciclovir PK.",
      sep = " "
    )
  )

  ini({
    # Structural PK -- Yang 2023 Table 3, Zhou et al. (1996) row. Both CL and Vc
    # are additive-linear in a covariate with a non-zero intercept, so the
    # exponentiated intercepts below are the covariate-free components and the
    # e_* parameters are additive slopes (not power exponents). Clearance in L/h,
    # volume in L. IV-only model, so no absorption parameters.
    lcl <- log(0.262); label("Covariate-free intercept of clearance (L/h)")               # Yang 2023 Table 3 (Zhou 1996): CL = 0.262 + (0.00271 * ASCC)
    lvc <- log(0.627); label("Covariate-free intercept of central volume (L)")            # Yang 2023 Table 3 (Zhou 1996): Vc = 0.627 + (0.437 * BW)

    # Additive-linear covariate slopes.
    e_crcl_cl <- 0.00271; label("Additive slope of ASCC on CL (L/h per mL/min/1.73 m^2)")  # Yang 2023 Table 3 (Zhou 1996): + 0.00271 * ASCC
    e_wt_vc   <- 0.437  ; label("Additive slope of WT on Vc (L per kg)")                   # Yang 2023 Table 3 (Zhou 1996): + 0.437 * BW

    # Between-subject variability. Yang 2023 Methods: %CV = sqrt(omega^2) * 100%,
    # so the diagonal variances are (BSV% / 100)^2: CL 35.4% -> 0.354^2 = 0.125316,
    # Vc 30.1% -> 0.301^2 = 0.090601.
    #
    # Yang 2023 Table 3 additionally reports "COV = 28.5" between CL and Vc, with
    # the footnote "COV: covariance between CL and Vc". Read literally as a
    # covariance of 28.5 (or 0.285) the matrix is not positive definite -- the
    # maximum admissible covariance for these two variances is
    # 0.354 * 0.301 = 0.106554, so 0.285 would imply a correlation of 2.67.
    # The only interpretation consistent with a valid covariance matrix is that
    # 28.5 is a CORRELATION of 0.285, giving
    # cov = 0.285 * 0.354 * 0.301 = 0.03036789. That reading is used here and is
    # documented in the validation vignette Assumptions and deviations section.
    etalcl + etalvc ~ c(0.125316,
                        0.03036789, 0.090601)  # Yang 2023 Table 3 (Zhou 1996): BSV CL = 35.4%, BSV Vc = 30.1%, COV = 28.5 read as correlation 0.285

    # Residual unexplained variability: proportional.
    propSd <- 0.0846; label("Proportional residual error (fraction)")  # Yang 2023 Table 3 (Zhou 1996): 8.46% proportional error
  })

  model({
    # Additive-linear covariate models with exponential between-subject variability
    # applied to the assembled typical value (the Georges_2009_ceftazidime.R
    # pattern for additive-intercept clearance).
    cl <- (exp(lcl) + e_crcl_cl * CRCL) * exp(etalcl)
    vc <- (exp(lvc) + e_wt_vc   * WT)   * exp(etalvc)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, volume in L -> concentration in mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
