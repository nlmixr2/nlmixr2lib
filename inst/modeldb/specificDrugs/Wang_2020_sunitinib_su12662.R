Wang_2020_sunitinib_su12662 <- function() {
  description <- "Two-compartment population PK model with first-order absorption and lag time for the active sunitinib metabolite SU012662 (N-desethyl sunitinib), fitted separately from the parent to SU012662 plasma concentrations after oral sunitinib in children and young adults (2-21 years) with refractory solid tumours. An assumed 21% conversion of the sunitinib dose to SU012662 enters as a fixed fraction of the dose; BSA is a power covariate on apparent clearance and apparent central volume (Wang 2020)."
  reference <- paste(
    "Wang E, DuBois SG, Wetmore C, Khosravan R.",
    "Population pharmacokinetics-pharmacodynamics of sunitinib in pediatric",
    "patients with solid tumors.",
    "Cancer Chemother Pharmacol. 2020;86(2):181-192.",
    "doi:10.1007/s00280-020-04106-z.",
    "The parent sunitinib model is modellib('Wang_2020_sunitinib').",
    sep = " "
  )
  vignette <- "Wang_2020_sunitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    BSA = list(
      description = "Baseline body surface area (DuBois and DuBois formula)",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline value; reference 1.47 m^2 (cohort median, Wang 2020 Table 2 and Table 3 footnote c). Power effects on CL/F and Vc/F.",
      source_name = "BSA"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested in a separate SCM run with BSA replaced by body weight; the BSA model had the lower OFV and was selected (Wang 2020 Results).",
      source_name = "WT"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "sunitinib (fraction converted to SU012662)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central_su12662 = list(analyte = "SU012662", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_su12662 = list(analyte = "SU012662", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 59L,
    n_studies = 2L,
    age_range = "2-21 years",
    weight_range = "16.2-100 kg (median 50.4 kg)",
    bsa_range = "0.66-2.14 m^2 (median 1.47 m^2)",
    sex_female_pct = 52.5,
    race_ethnicity = c(Asian = 5.1, NonAsian = 89.8, Unknown = 5.1),
    disease_state = "Children and young adults with refractory solid tumours, predominantly high-grade glioma, ependymoma, brain stem glioma, or sarcoma (Children's Oncology Group studies ADVL0612 and ACNS1021).",
    dose_range = "Sunitinib 15 or 20 mg/m^2 orally once daily on schedule 4/2 (4 weeks on, 2 weeks off).",
    regions = "United States and Canada (Children's Oncology Group).",
    notes = "Baseline demographics from Wang 2020 Table 2. 340 SU012662 plasma observations; LLOQ 1 ng/mL (ADVL0612) or 0.1 ng/mL (ACNS1021)."
  )

  ini({
    # Assumed conversion of sunitinib to SU012662 (Wang 2020 Results: 'a conversion
    # of 21% of sunitinib to SU012662 was assumed'), applied to the sunitinib dose.
    lfdepot <- fixed(log(0.21)); label("Assumed fraction of the sunitinib dose converted to SU012662 (unitless)") # Results 'Sunitinib and SU012662 base and final PK models' paragraph 1

    # Structural parameters (Wang 2020 Table 3, SU012662 final model, typical patient with BSA 1.47 m^2)
    lka <- log(0.28); label("Apparent first-order input rate constant ka for SU012662 (1/h)") # Table 3 SU012662 'ka' = 0.28 (RSE 36.8%)
    ltlag <- log(0.46); label("Input lag time tlag (h)") # Table 3 SU012662 'tlag' = 0.46 (RSE 76.9%)
    lcl_su12662 <- log(10.9); label("Apparent SU012662 clearance CL/F at BSA 1.47 m^2 (L/h)") # Table 3 SU012662 'CL/F' = 10.9 (RSE 7.5%)
    lvc_su12662 <- log(1030); label("Apparent SU012662 central volume Vc/F at BSA 1.47 m^2 (L)") # Table 3 SU012662 'Vc/F' = 1030 (RSE 15.3%)
    lvp_su12662 <- log(122); label("Apparent SU012662 peripheral volume Vp/F (L)") # Table 3 SU012662 'Vp/F' = 122 (RSE 72.5%)
    lq_su12662 <- log(17.8); label("Apparent SU012662 intercompartmental clearance Q/F (L/h)") # Table 3 SU012662 'Q/F' = 17.8 (RSE 110%)

    # Covariate effects
    e_bsa_cl_su12662 <- 0.843; label("Power exponent of BSA/1.47 on SU012662 CL/F (unitless)") # Results text CL/F = 10.9 * (BSA/1.47)^0.843; Table 3 'BSA on CL/F' = 0.84 (RSE 30.8%)
    e_bsa_vc_su12662 <- 1.72; label("Power exponent of BSA/1.47 on SU012662 Vc/F (unitless)") # Results text Vc/F = 1030 * (BSA/1.47)^1.72; Table 3 'BSA on Vc/F' = 1.72 (RSE 23.5%)

    # IIV: Table 3 reports omega as CV%; omega^2 = log(1 + CV^2)
    etalcl_su12662 ~ 0.2081 # Table 3 SU012662 'omega (CL/F)' = 48.1%
    etalvc_su12662 ~ 0.2223 # Table 3 SU012662 'omega (Vc/F)' = 49.9%
    etalka ~ 0.4215 # Table 3 SU012662 'omega (ka)' = 72.4%

    # Residual error
    propSd_su12662 <- 0.231; label("Proportional residual error on SU012662 (fraction)") # Table 3 SU012662 'sigma' = 23.1% (RSE 4.7%)
  })

  model({
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)
    fdepot <- exp(lfdepot)
    cl_su12662 <- exp(lcl_su12662 + etalcl_su12662) * (BSA / 1.47)^e_bsa_cl_su12662
    vc_su12662 <- exp(lvc_su12662 + etalvc_su12662) * (BSA / 1.47)^e_bsa_vc_su12662
    vp_su12662 <- exp(lvp_su12662)
    q_su12662 <- exp(lq_su12662)

    kel_su12662 <- cl_su12662 / vc_su12662
    k12_su12662 <- q_su12662 / vc_su12662
    k21_su12662 <- q_su12662 / vp_su12662

    d/dt(depot) <- -ka * depot
    d/dt(central_su12662) <- ka * depot - kel_su12662 * central_su12662 -
      k12_su12662 * central_su12662 + k21_su12662 * peripheral1_su12662
    d/dt(peripheral1_su12662) <- k12_su12662 * central_su12662 - k21_su12662 * peripheral1_su12662

    f(depot) <- fdepot
    alag(depot) <- tlag

    # mg / L * 1000 = ng/mL
    Cc_su12662 <- central_su12662 / vc_su12662 * 1000
    Cc_su12662 ~ prop(propSd_su12662)
  })
}
