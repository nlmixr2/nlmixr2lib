Danielak_2017_treosulfan <- function() {
  description <- "Two-compartment IV-infusion population PK model for treosulfan (TREO) in pediatric patients undergoing conditioning prior to hematopoietic stem cell transplantation (Danielak 2017). Allometric body-weight scaling normalised to a 70 kg adult typical value with exponents fixed at 0.75 on CL and 1 on V1 and V2; Q has no weight covariate. Correlated IIV on CL and V1 (Cl-V1 correlation 0.714); independent IIV on Q. Proportional residual error."
  reference   <- "Danielak D, Twardosz J, Kasprzyk A, Wachowiak J, Kalwak K, Glowka F. Population pharmacokinetics of treosulfan and development of a limited sampling strategy in children prior to hematopoietic stem cell transplantation. Eur J Clin Pharmacol. 2018 Jan;74(1):79-89. doi:10.1007/s00228-017-2344-x"
  vignette    <- "Danielak_2017_treosulfan"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Treated as baseline (single-day sampling on first day of therapy). Allometric scaling on CL (exponent 0.75 fixed), V1 (1 fixed), and V2 (1 fixed); not applied to Q (paper reports addition of WT on Q worsened the fit). Reference weight 70 kg adult typical (Danielak 2017 Eq. 9 + Table 2; cohort range 7.7-52 kg).",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 15L,
    n_studies      = 1L,
    age_range      = "0.4-15 years",
    age_mean       = "7.8 +/- 4.9 years",
    weight_range   = "7.7-52 kg",
    weight_mean    = "26.9 +/- 15.7 kg",
    bsa_range      = "0.25-1.63 m^2",
    bsa_mean       = "0.95 +/- 0.44 m^2",
    sex_female_pct = 20,
    disease_state  = "Pediatric patients with malignant or non-malignant disorders receiving treosulfan-based conditioning regimens prior to allogeneic hematopoietic stem cell transplantation. Diagnoses include ALL (4), AML (1), CML (1), neuroblastoma (2), Ewing sarcoma (2), adrenoleukodystrophy (2), Diamond-Blackfan anemia (1), severe congenital neutropenia (1), and Wiskott-Aldrich syndrome (1).",
    dose_range     = "Treosulfan IV infusion at 10, 12, or 14 g/m^2 daily as a 1 h or 2 h infusion (sampled on day 1 of therapy). Specific regimens: 10 g/m^2 over 1 h (n = 1); 12 g/m^2 over 1 h (n = 4); 12 g/m^2 over 2 h (n = 4); 14 g/m^2 over 2 h (n = 6). Body surface area calculated by the Mosteller method.",
    regions        = "Poland (Poznan and Wroclaw)",
    n_observations = "110 plasma treosulfan concentrations across 15 subjects (sparse sampling: 6 timepoints in 7 subjects; 10 timepoints in 8 subjects, all on day 1).",
    notes          = "Patients recruited 2007-2011 at the Department of Oncology, Hematology and Pediatric Transplantation (Poznan University of Medical Sciences) and the Department of Pediatric Hematology, Oncology and Bone Marrow Transplantation (Wroclaw Medical University). Creatinine clearance available for 8 of 15 patients (123 +/- 60 mL/min, range 71-239); not tested as a covariate due to missingness and normal renal function in the available data. Patient sex was tested and not significant (only 3 girls). Demographics from Danielak 2017 Table 1; model and covariate analysis from Eq. 9 and Table 2."
  )

  ini({
    # Structural parameters at adult reference 70 kg (Danielak 2017 Table 2,
    # final-model SAEM estimates, percent-RSE in parentheses). Reference
    # weight is the conventional 70 kg adult typical value used with the
    # paper's allometric exponents (Eq. 9 + Methods "Covariate selection").
    lcl <- log(14.7);  label("Clearance at 70 kg reference (CL, L/h)")                # Danielak 2017 Table 2: CL 14.7 L/h/70 kg (RSE 6.9 percent)
    lvc <- log(26.0);  label("Central volume at 70 kg reference (V1, L)")             # Danielak 2017 Table 2: V1 26.0 L/70 kg (RSE 14.0 percent)
    lq  <- log(2.25);  label("Intercompartmental clearance (Q, L/h)")                 # Danielak 2017 Table 2: Q 2.25 L/h (RSE 22.2 percent); no weight covariate
    lvp <- log(9.93);  label("Peripheral volume at 70 kg reference (V2, L)")          # Danielak 2017 Table 2: V2 9.93 L/70 kg (RSE 9.0 percent)

    # Allometric exponents on weight; paper estimated 0.804 (CL), 0.959 (V1),
    # 0.925 (V2) and then fixed them at 0.75 (CL) and 1 (V1, V2) per
    # Anderson and Holford (Danielak 2017 Results paragraph 1 and Table 2
    # entries marked "(fixed)"). Q has no weight covariate.
    e_wt_cl <- fixed(0.75);  label("Allometric exponent on CL (unitless, fixed)")     # Danielak 2017 Table 2: "beta Cl, weight 0.75 (fixed)"
    e_wt_vc <- fixed(1.0);   label("Allometric exponent on V1 (unitless, fixed)")     # Danielak 2017 Table 2: "beta V1, weight 1 (fixed)"
    e_wt_vp <- fixed(1.0);   label("Allometric exponent on V2 (unitless, fixed)")     # Danielak 2017 Table 2: "beta V2, weight 1 (fixed)"

    # Correlated IIV on CL and V1. Paper reports omega-style IIV in percent
    # (Danielak 2017 Table 2: omega_Cl 25.5%, omega_V1 51.4%, omega_Cl-V1
    # 71.4%). With Monolix's exponential / lognormal IIV model, the
    # diagonal omegas are SDs of eta on the log scale (omega = 0.255,
    # 0.514) and the off-diagonal is the correlation (rho = 0.714).
    # Internal variance / covariance:
    #   var(etalcl)  = 0.255^2          = 0.065025
    #   var(etalvc)  = 0.514^2          = 0.264196
    #   cov(etalcl, etalvc) = 0.714 * 0.255 * 0.514 = 0.093596
    etalcl + etalvc ~ c(0.065025,
                        0.093596, 0.264196)                                            # Danielak 2017 Table 2: IIV CL 25.5% (RSE 19.8); IIV V1 51.4% (RSE 20.0); Cl-V1 71.4% (RSE 20.7)

    # Independent IIV on Q. IIV on V2 was estimated but the standard error
    # exceeded 100 percent so the authors removed it from the final model
    # (Danielak 2017 Results paragraph 1).
    etalq ~ 0.148996                                                                   # Danielak 2017 Table 2: IIV Q 38.6% (RSE 52.3); var = 0.386^2

    # Proportional residual error. Paper retained the proportional-only
    # model (Danielak 2017 Methods "Structural and error model selection"
    # and Table 2).
    propSd <- 0.188;  label("Proportional residual error (fraction)")                  # Danielak 2017 Table 2: Residual proportional error 0.188 (RSE 9.01)
  })

  model({
    # Reference body weight for allometric scaling (Danielak 2017 Eq. 9
    # standardises WT to a 70 kg adult typical value).
    ref_wt <- 70

    # Size scaling on baseline body weight.
    size_cl <- (WT / ref_wt)^e_wt_cl
    size_vc <- (WT / ref_wt)^e_wt_vc
    size_vp <- (WT / ref_wt)^e_wt_vp

    # Individual PK parameters. Q has no covariate scaling (Danielak 2017
    # Results: "addition of bodyweight as a covariate for Q increased the
    # MOFV and worsened the model fit").
    cl <- exp(lcl + etalcl) * size_cl
    vc <- exp(lvc + etalvc) * size_vc
    q  <- exp(lq  + etalq)
    vp <- exp(lvp)          * size_vp

    # Micro-constants for the linear two-compartment ODE.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV-infusion ODE system. Dose enters `central`
    # directly via the rxode2 event-table infusion (rate/duration).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Observation. Dose in mg / V in L -> Cc in mg/L (matches the
    # plasma concentrations in Danielak 2017 Fig. 2 and the HPLC-MS/MS
    # quantitation range 0.0566 mg/L - 1.59 mg/L from the Methods).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
