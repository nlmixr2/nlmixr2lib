Sun_2024_tacrolimus <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral immediate-release tacrolimus in adult renal transplant recipients at least one year post-transplant who were co-administered oral voriconazole. The absorption rate constant is fixed. Apparent clearance falls linearly as the measured voriconazole concentration rises, falls with red blood cell count and rises with platelet count; apparent volume of distribution falls linearly with voriconazole concentration. Combined additive and proportional residual error."
  reference <- paste(
    "Sun Z-H, Zhao Y-C, Li J-K, Peng F, Yu F, Zhang B-K, Yan M.",
    "Population dynamics analysis of the interaction between",
    "tacrolimus and voriconazole in renal transplant recipients.",
    "Front Pharmacol. 2024;15:1502097. doi:10.3389/fphar.2024.1502097.",
    "Parameter values are the final-model estimates in Table 3, cross-checked",
    "against the bootstrap summary in Table 4; covariate forms from the",
    "section 3.2 equations for V and Cl.",
    sep = " "
  )
  vignette <- "Sun_2024_tacrolimus"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  # Section 2.3.1: tacrolimus measured with the ARCHITECT Tacrolimus
  # chemiluminescent immunoassay, a whole-blood assay, so V/F and CL/F are
  # whole-blood apparent parameters.
  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    CONC_VORI_NGML = list(
      description = "Measured voriconazole concentration in the co-administered patient (the source's C_VRC / Cvrc); 0 when no voriconazole is detectable",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "TIME-VARYING within subject: voriconazole was measured alongside tacrolimus (section 2.3.1). UNIT CONVERSION: Sun 2024 reports Cvrc in ug/mL (Table 1 median 0.95 [IQR 0.00, 2.64]; assay linear range 0.35-11.26 ug/mL) whereas the canonical column is ng/mL, so the model divides by 1000 before applying the published coefficients. The effect is LINEAR and centred at 1.70 ug/mL on both CL/F and V/F, as printed in the section 3.2 equations: (1 + (Cvrc - 1.70) * theta). Because the term is linear, CL/F reaches zero at Cvrc = 1.70 + 1/0.12 = 10.03 ug/mL and V/F at 1.70 + 1/0.07 = 15.99 ug/mL; the model must not be used near or above those concentrations (the cohort's upper quartile is 2.64 ug/mL).",
      source_name = "Cvrc"
    ),
    RBC = list(
      description = "Red blood cell count",
      units = "10^12 cells/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Laboratory value collected with the tacrolimus samples; time-varying status not stated in the paper. Enters as the power term (RBC / 3.49)^-1.34 on CL/F, so CL/F FALLS as RBC rises (section 3.2 equation; the Abstract and Results say 'negative correlation'). 3.49 is the centring constant printed in the equation; the Table 1 cohort median is 3.44 (IQR 2.82-4.12), so 3.49 is presumably the cohort mean.",
      source_name = "RBC"
    ),
    PLT = list(
      description = "Platelet count",
      units = "10^9 cells/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Laboratory value collected with the tacrolimus samples; time-varying status not stated in the paper. Enters as the power term (PLT / 181.24)^0.46 on CL/F, so CL/F RISES with platelet count (section 3.2 equation). 181.24 is the centring constant printed in the equation; the Table 1 cohort median is 170 (IQR 129-220), so 181.24 is presumably the cohort mean.",
      source_name = "PLT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 51L,
    n_studies = 1L,
    age_range = "median 45 years (IQR 38-55); adults >= 18 years",
    weight_range = "median 59.6 kg (IQR 51.1-67.2)",
    sex_female_pct = 27.45,
    disease_state = "Adult renal transplant recipients at least one year post-transplant on an oral triple immunosuppressive regimen (immediate-release tacrolimus twice daily + mycophenolate mofetil + glucocorticoid) and receiving concomitant oral voriconazole, mostly for invasive fungal infection",
    dose_range = "Oral immediate-release tacrolimus, typically twice daily, clinician-titrated by therapeutic drug monitoring; doses not tabulated",
    regions = "Single centre, The Second Xiangya Hospital of Central South University, Changsha, Hunan, People's Republic of China",
    notes = "Prospective non-interventional cohort, April 2023 to September 2024 (ethics approval LYEC2024-K0106). 281 tacrolimus whole-blood concentrations (at least five per patient) from a sparse design with three randomised sampling groups (0, 0.5, 1 h; 2, 4, 6 h; 8 h and pre-dose). Observed tacrolimus concentration median 5.80 ng/mL (IQR 4.20-8.65); voriconazole concentration median 0.95 ug/mL (IQR 0.00-2.64); RBC median 3.44 x 10^12/L (IQR 2.82-4.12); PLT median 170 x 10^9/L (IQR 129-220); serum creatinine median 172 umol/L (IQR 120-277); albumin median 34.4 g/L; haematocrit median 31.8%. Race/ethnicity not reported (single Chinese centre). Table 1 reports 37 male (72.55%) and 14 female; section 3.1 says 38 males. Estimated in Phoenix NLME 8.1; 21 covariates screened by forward inclusion (p < 0.05) and backward elimination (p < 0.01); validated by 1500-sample bootstrap and VPC."
  )

  ini({
    # Structural parameters: Sun 2024 Table 3 (final model), repeated with
    # bootstrap summaries in Table 4. Ka is labelled 'Ka (Fix)' with Stderr 0
    # and a degenerate 3.09-3.09 confidence interval. V/F and CL/F are the
    # typical values at the centring covariates Cvrc = 1.70 ug/mL,
    # RBC = 3.49 x 10^12/L and PLT = 181.24 x 10^9/L.
    lka <- fixed(log(3.09)); label("Absorption rate constant Ka (1/h)")                                                            # Sun 2024 Table 3: Ka (Fix) = 3.09 1/h (Stderr 0.00, CI 3.09-3.09)
    lvc <- log(1635.34);     label("Apparent volume of distribution V/F at Cvrc = 1.70 ug/mL (L)")                                  # Sun 2024 Table 3: V/F = 1635.34 L (bootstrap median 1626.24, 95% CI 1313.12-2013.20)
    lcl <- log(5.73);        label("Apparent clearance CL/F at Cvrc = 1.70 ug/mL, RBC = 3.49 x 10^12/L, PLT = 181.24 x 10^9/L (L/h)")  # Sun 2024 Table 3: CL/F = 5.73 L/h (bootstrap median 5.74, 95% CI 4.65-7.07)

    # Covariate effects: Sun 2024 Table 3 rows dCldRBC, dCldPLT, dCldCvRc,
    # dVdCVRC, with the functional forms of the section 3.2 equations
    #   V  = 1635.34 * (1 + (Cvrc - 1.70) * (-0.07)) * exp(nV)
    #   Cl = 5.73 * (RBC/3.49)^(-1.34) * (PLT/181.24)^(0.46)
    #          * (1 + (Cvrc - 1.70) * (-0.12)) * exp(nCL)
    # with Cvrc in ug/mL.
    e_rbc_cl            <- -1.34; label("Power exponent of (RBC / 3.49 x 10^12/L) on CL/F (unitless)")                             # Sun 2024 Table 3: dCldRBC = -1.34 (Stderr 0.29; CI -1.92 to -0.77)
    e_plt_cl            <-  0.46; label("Power exponent of (PLT / 181.24 x 10^9/L) on CL/F (unitless)")                            # Sun 2024 Table 3: dCldPLT = 0.46 (Stderr 0.20; CI 0.07 to 0.85)
    e_conc_vori_ngml_cl <- -0.12; label("Linear coefficient of (voriconazole concentration - 1.70 ug/mL) on CL/F (per ug/mL)")    # Sun 2024 Table 3: dCldCvRc = -0.12 (Stderr 0.02; CI -0.16 to -0.08)
    e_conc_vori_ngml_vc <- -0.07; label("Linear coefficient of (voriconazole concentration - 1.70 ug/mL) on V/F (per ug/mL)")     # Sun 2024 Table 3: dVdCVRC = -0.07 (Stderr 0.02; CI -0.11 to -0.03)

    # Inter-individual variability: Sun 2024 Table 4 reports these as
    # 'omega 2 V' and 'omega 2 CL', i.e. variances of the log-normal random
    # effects (exp(nV), exp(nCL) in the section 3.2 equations); no covariance
    # is reported, so the block is diagonal.
    etalvc ~ 0.091  # Sun 2024 Table 4: omega^2 V = 0.091 (bootstrap median 0.091)
    etalcl ~ 0.039  # Sun 2024 Table 4: omega^2 CL = 0.039 (bootstrap median 0.063)

    # Residual error: Phoenix NLME 'Add-Multiplicative' model selected in
    # Table 2 (1a_Add_Multiplicative), whose variance is
    # stdev0^2 + (CMultStdev * C)^2 -- the nlmixr2 combined2 form.
    propSd <- 0.36; label("Proportional residual error (fraction)")                          # Sun 2024 Table 3: tvCMultStdev = 0.36 (Stderr 0.02); section 3.2: CMultStdev = 0.36
    addSd  <- 0.04; label("Additive residual error for tacrolimus whole blood (ng/mL)")       # Sun 2024 Table 3: stdev0 = 0.04 (Stderr 0.007); Table 4 sigma = 0.04
  })

  model({
    # The published coefficients are on the ug/mL scale the paper reports
    # voriconazole in; the canonical covariate column is ng/mL.
    cvrc <- CONC_VORI_NGML / 1000

    ka <- exp(lka)
    vc <- exp(lvc + etalvc) * (1 + (cvrc - 1.70) * e_conc_vori_ngml_vc)
    cl <- exp(lcl + etalcl) *
      (RBC / 3.49)^e_rbc_cl *
      (PLT / 181.24)^e_plt_cl *
      (1 + (cvrc - 1.70) * e_conc_vori_ngml_cl)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Doses are in mg and V/F in L, so central / vc is mg/L; whole-blood
    # tacrolimus is reported in ng/mL and 1 mg/L = 1000 ng/mL.
    Cc <- 1000 * central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
