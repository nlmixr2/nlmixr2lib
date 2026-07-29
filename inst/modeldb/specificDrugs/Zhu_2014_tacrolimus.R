Zhu_2014_tacrolimus <- function() {
  description <- "Two-compartment population PK model for oral tacrolimus in Chinese adult liver transplant recipients (Zhu 2014), with first-order absorption, a power-form joint DOSE x POD covariate effect on apparent clearance, log-normal IIV on CL/F, V2/F, Q/F, V3/F, and ka, and proportional residual error. Bioavailability was not estimated; the structural disposition parameters are apparent values (CL/F, V/F, Q/F)."
  reference <- "Zhu L, Wang H, Sun X, Rao W, Qu W, Zhang Y, Sun L. The Population Pharmacokinetic Models of Tacrolimus in Chinese Adult Liver Transplantation Patients. Journal of Pharmaceutics. 2014;2014:713650. doi:10.1155/2014/713650"
  vignette <- "Zhu_2014_tacrolimus"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    DOSE = list(
      description        = "Total tacrolimus oral dose per day",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying within subject; reflects the daily-total tacrolimus dose at the time of the record (twice-daily oral capsules adjusted by therapeutic drug monitoring). Zhu 2014 Table 1: mean 5.31 mg/d, median 5 mg/d, range 1-10.5 mg/d. Enters the final covariate equation as a power term on CL/F (Eq. 3): CL/F = theta_CL/F * DOSE^theta_DOSE * POD^theta_POD; theta_DOSE = 0.371. Because the relationship is multiplicative and a power form, DOSE must be > 0 -- the dataset is restricted to on-treatment records by construction.",
      source_name        = "DOSE"
    ),
    POD = list(
      description        = "Postoperative days since liver transplantation",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying within subject; rises monotonically from the transplant date. Zhu 2014 Table 1: mean 20.71 days, median 14 days, range 2-85 days. Enters the final covariate equation as a power term on CL/F (Eq. 3): CL/F = theta_CL/F * DOSE^theta_DOSE * POD^theta_POD; theta_POD = 0.127. POD must be > 0 because POD^theta_POD is undefined at POD = 0; the source dataset begins at POD = 2 d post-transplant.",
      source_name        = "POD"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 47L,
    n_studies        = 1L,
    n_observations   = 435L,
    age_range        = "25-78 years",
    age_median       = "60 years",
    age_mean_sd      = "57.47 +/- 11.16 years",
    weight_range     = "not collected (excluded by Zhu 2014 due to extensive missing weight data among inpatients)",
    sex_female_pct   = 100 * 20 / 47,
    race_ethnicity   = c(Asian = 100),
    disease_state    = "Adult liver transplant recipients on triple immunosuppression (tacrolimus + mycophenolate mofetil + corticosteroids) during inpatient hospitalisation after orthotopic liver transplantation.",
    dose_range       = "Oral tacrolimus capsules (0.5 mg and 1 mg). Therapy initiated at 0.1-0.15 mg/kg twice daily and adjusted by therapeutic drug monitoring to a trough target of 10-15 ng/mL during the first three months posttransplant. Observed daily doses 1-10.5 mg/d (median 5 mg/d).",
    pod_range        = "2-85 days posttransplant (median 14 d, mean 20.71 d).",
    regions          = "China (Tianjin First Central Hospital, Beijing Friendship Hospital).",
    co_medication    = "Mycophenolate mofetil and corticosteroids as part of triple immunosuppressive regimen. Patients on fluconazole, diltiazem, or other strong CYP3A4 inhibitors/inducers were not specifically excluded but were rare in the dataset.",
    sampling_design  = "Steady-state inpatient sampling. Daily predose troughs after transplantation until concentrations stabilised; thereafter blood samples 3x weekly or more frequently. A subset had full profiles (predose, 0.3, 1, 1.5, 2, 4, 6, 8, 12 h postdose).",
    assay            = "Microparticle enzyme immunoassay (MEIA, IMx platform). LLOQ 1.5 ng/mL; linear 1.5-30 ng/mL; samples above LOQ diluted per manufacturer. The antitacrolimus monoclonal recognises parent drug plus three metabolites (M-II, M-III, M-V); cross-reactivity with other metabolites was below the minimum detectable sensitivity.",
    notes            = "Single-centre prospective TDM cohort (Tianjin First Hospital, China, 2008-2011) -- 27 male and 20 female adults. Body weight and transplant type (whole vs cut-down graft) were excluded as candidate covariates because of extensive missing values among the inpatient population. Eleven candidate covariates were screened; only DOSE and POD on CL/F survived forward inclusion (alpha = 0.01) and backward elimination (alpha = 0.005). DOSE on V2/F was tested but eliminated (Table 3, eliminate-DOSE-on-V2 row: dOFV +4.8, p > 0.01)."
  )

  ini({
    # Structural PK -- Zhu 2014 Table 4 final-model estimates. Time in hours;
    # apparent clearances (CL/F, Q/F) in L/h; apparent volumes (V2/F = central,
    # V3/F = peripheral) in L; ka in 1/h. Bioavailability was not estimated; all
    # disposition parameters are apparent values (CL/F, V/F, Q/F).
    #
    # The Table 4 "CL/F = 11.2 L/h" is the theta_CL/F coefficient evaluated at
    # DOSE = 1 mg/d and POD = 1 d (Eq. 3 reference point). The typical-value
    # CL/F at the cohort medians (DOSE = 5 mg/d, POD = 14 d) is
    # 11.2 * 5^0.371 * 14^0.127 = 11.2 * 1.835 * 1.401 ~ 28.8 L/h, close to the
    # basic-model CL/F of 30.2 L/h reported in Table 2.
    lka  <- log(0.723)  ; label("Absorption rate constant ka (1/h)")  # Zhu 2014 Table 4 ka = 0.723 1/h
    lcl  <- log(11.2)   ; label("Apparent oral clearance coefficient theta_CL/F at DOSE = 1 mg/d, POD = 1 d (L/h)")  # Zhu 2014 Table 4 CL/F = 11.2 L/h
    lvc  <- log(406)    ; label("Apparent central volume V2/F (L)")   # Zhu 2014 Table 4 V2 = 406 L
    lq   <- log(57.3)   ; label("Apparent inter-compartmental clearance Q/F (L/h)")  # Zhu 2014 Table 4 Q = 57.3 L/h
    lvp  <- log(503)    ; label("Apparent peripheral volume V3/F (L)") # Zhu 2014 Table 4 V3 = 503 L

    # Covariate effects on CL/F -- Zhu 2014 Eq. 3 final covariate model:
    # CL/F = theta_CL/F * DOSE^theta_DOSE * POD^theta_POD.
    e_dose_cl <- 0.371 ; label("Power exponent of DOSE (mg/d) on CL/F (unitless)") # Zhu 2014 Table 4 theta_DOSE = 0.371
    e_pod_cl  <- 0.127 ; label("Power exponent of POD (d) on CL/F (unitless)")     # Zhu 2014 Table 4 theta_POD = 0.127

    # Inter-individual variability -- Zhu 2014 Table 4 reports IIV as %CV from
    # a log-normal model (Eq. 1: theta_ij = theta * exp(eta_ij)). Convert to
    # log-scale variance via omega^2 = log(CV^2 + 1):
    #   CL/F  CV  16.2% -> log(0.162^2 + 1) = 0.02591
    #   V2/F  CV 163%   -> log(1.63^2 + 1)  = 1.29657
    #   Q/F   CV  19.7% -> log(0.197^2 + 1) = 0.03808
    #   V3/F  CV 199%   -> log(1.99^2 + 1)  = 1.60139
    #   ka    CV  74.3% -> log(0.743^2 + 1) = 0.43972
    # The paper reports diagonal IIV only (no correlation block); each eta is
    # therefore independent. The V2/F and V3/F IIVs (163% and 199% CV) are very
    # large with poor RSE (164% and 231%) -- see vignette Assumptions and
    # deviations for the identifiability caveat.
    etalcl ~ 0.02591  # Zhu 2014 Table 4 omega CL/F = 16.2%
    etalvc ~ 1.29657  # Zhu 2014 Table 4 omega V2/F = 163%
    etalq  ~ 0.03808  # Zhu 2014 Table 4 omega Q/F = 19.7%
    etalvp ~ 1.60139  # Zhu 2014 Table 4 omega V3/F = 199%
    etalka ~ 0.43972  # Zhu 2014 Table 4 omega ka = 74.3%

    # Residual unexplained variability -- Zhu 2014 Eq. 2:
    # Y = IPRED * (1 + epsilon) -- proportional error model on the linear
    # concentration scale; sigma = 26.54%.
    propSd <- 0.2654 ; label("Proportional residual error (fraction)") # Zhu 2014 Table 4 sigma = 26.54%
  })

  model({
    # Individual PK parameters with Zhu 2014 Eq. 3 covariate equation on CL/F.
    # DOSE (mg/d) and POD (d) are time-varying record-level covariates; both
    # must be strictly positive (the power form is undefined at zero). The
    # source dataset begins at POD >= 2 and tacrolimus dosing is non-zero by
    # cohort construction.
    cl <- exp(lcl + etalcl) * DOSE^e_dose_cl * POD^e_pod_cl
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)
    ka <- exp(lka + etalka)

    # Two-compartment oral PK with first-order absorption. Dose lands in
    # `depot`; bioavailability is not parameterised (CL/F and V/F are
    # apparent). Micro-constants are spelled out for clarity.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Internal concentration is mg/L (dose in mg, vc in L). Tacrolimus assay
    # reports concentrations in ng/mL (= ug/L). 1 mg/L = 1000 ng/mL, so
    # multiply by 1000 to map to the assay scale.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
