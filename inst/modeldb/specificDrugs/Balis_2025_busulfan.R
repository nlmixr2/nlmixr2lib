Balis_2025_busulfan <- function() {
  description <- "One-compartment IV population PK model for busulfan in 328 infants and children (median age 5.9 years, range 0.11-28.1) receiving busulfan-containing hematopoietic stem cell transplant conditioning regimens with therapeutic drug monitoring (Balis 2025). Body surface area enters both clearance and volume as an estimated power effect normalized to the cohort median BSA of 0.81 m^2 (exponent 1.14 on CL, 1.33 on V), which the stepwise covariate search selected over body weight. The BSA-scaled model underpins the paper's proposed dosing method: 100 mg/m^2/day for patients with BSA >= 0.5 m^2 plus a BSA-banded dosing table for smaller infants. The residual error model was multiplicative but its magnitude was not reported, so the proportional residual standard deviation is encoded as fixed(0); see the vignette Errata."
  reference   <- "Balis FM, Rieger E, Bunin NJ, Gardiner J, Shaw LM, Olson TS, Milone MC (2025) New approach to busulfan dosing in infants and children based on a population pharmacokinetic analysis. Cancer Chemother Pharmacol 95:32. doi:10.1007/s00280-025-04757-w"
  vignette    <- "Balis_2025_busulfan"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area on the day prior to the first busulfan dose",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Selected by the Phoenix NLME Stepwise Covariate Search as a covariate on both CL and V, in preference to body weight (Balis 2025 Methods, 'Pharmacokinetic analysis'; Supplemental Table 6). Entered as a power effect normalized to the cohort median 0.81 m^2 (Methods: 'Covariates were normalized to median values of BSA (0.81 m2) and body weight (20.5 kg)'). Phoenix names the two exponents dCLdBSA and dVdBSA; here they are e_bsa_cl and e_bsa_vc. The paper does not state which BSA formula was used to compute the column, only that weight, height and BSA were abstracted from the medical record (Supplemental Table 1).",
      source_name        = "BSA"
    )
  )

  compartmentData <- list(
    central = list(analyte = "busulfan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 328,
    n_studies      = 1,
    age_range      = "0.11-28.1 years",
    age_median     = "5.9 years",
    weight_range   = "3.4-126 kg",
    weight_median  = "20.5 kg",
    bsa_range      = "0.21-2.42 m^2",
    bsa_median     = "0.81 m^2",
    sex_female_pct = 40.5,
    disease_state  = "Malignant (n = 172; 131 leukemia, 36 neuroblastoma, 5 lymphoma) and non-malignant (n = 156; 97 immunodeficiency, 47 non-malignant hematologic, 12 inborn error of metabolism) disease treated with a busulfan-containing hematopoietic stem cell transplant conditioning regimen.",
    dose_range     = "3.2 mg/kg/dose infused over 3 h daily x 4 doses, or 0.8 mg/kg/dose infused over 2 h q6h x 16 doses (12.8 mg/kg/course on both schedules). 38 children >10 kg and <48 months on the daily schedule received 4 mg/kg/dose; 44 patients on the q6h schedule received 1 mg/kg/dose. Doses were adjusted after day 1 to target an AUCinf of 14.8-24.6 mg*h/L (3600-6000 uM*min); the actual administered doses were used for the PK analysis.",
    regions        = "United States (single institution: Children's Hospital of Philadelphia, Feb 2007 to May 2023)",
    renal_function = "Laboratory values for kidney and liver function were within the reference interval for essentially all subjects and were therefore not tested as covariates (Balis 2025 Methods, 'Population').",
    notes          = "Baseline demographics from Balis 2025 Table 2 (195 males, 133 females; 3339 busulfan concentrations total, 2149 of which were measured after the first dose and used for the covariate model per Supplemental Figure 3). 174 patients were treated on the q6h schedule and 154 on the daily schedule. Retrospective single-institution therapeutic drug monitoring data set; model fitted in Phoenix NLME v8.3 with FOCE-ELS. The parameters encoded here are the final covariate model of Supplemental Table 6, fitted to the first-dose (day 1) data."
  )

  ini({
    # Structural parameters -- typical values at the reference BSA of 0.81 m^2
    lcl <- log(4.42); label("Clearance at BSA 0.81 m^2 (L/h)")                                # Balis 2025 Supplemental Table 6: tvCl = 4.42 L/h (RSE 1.33%; bootstrap 4.41). 4.42 / 0.81 = 5.46 L/(h*m^2), matching the 5.47 L/(h*m^2) of Table 3
    lvc <- log(14.4); label("Central volume of distribution at BSA 0.81 m^2 (L)")              # Balis 2025 Supplemental Table 6: tvV = 14.4 L (RSE 0.91%; bootstrap 14.3). 14.4 / 0.81 = 17.8 L/m^2, matching Table 3

    # Covariate effects -- power form on (BSA / 0.81 m^2), both estimated
    e_bsa_cl <- 1.14; label("Power exponent of (BSA / 0.81 m^2) on CL (unitless)")             # Balis 2025 Supplemental Table 6: dCLdBSA = 1.14 (RSE 2.20%; bootstrap 1.14)
    e_bsa_vc <- 1.33; label("Power exponent of (BSA / 0.81 m^2) on Vc (unitless)")             # Balis 2025 Supplemental Table 6: dVdBSA = 1.33 (RSE 1.46%; bootstrap 1.34)

    # IIV -- exponential on CL and V (Balis 2025 Methods: Vi = tvV * exp(etaV), CLi = tvCL * exp(etaCL))
    etalcl ~ 0.0572  # Balis 2025 Supplemental Table 6 caption: Omega matrix variance for hCL = 0.0572 (CV 24.3%; the earlier supplement revision reports the same quantity as 24.0% and Figure 2 reports a 23% C.V. for BSA-normalized CLi)
    etalvc ~ 0.0279  # Balis 2025 Supplemental Table 6 caption: Omega matrix variance for hV = 0.0279 (CV 16.8%; the earlier supplement revision reports the same quantity as 16.6%)

    # Residual error
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - magnitude not reported in source)")  # Balis 2025 Methods states only "A multiplicative random error model ... was used"; no sigma value appears in the paper or either supplement revision. See vignette Errata.
  })

  model({
    # Individual PK parameters. BSA enters CL and V as a power effect normalized
    # to the cohort median BSA of 0.81 m^2. The power form (rather than a linear
    # or exponential covariate form) is confirmed by the paper's own subgroup
    # fits: at the BSA implied by Table 3 for patients <=3 years (0.429 m^2) this
    # gives CL 2.14 L/h and V 6.19 L against the reported 2.12 L/h and 6.24 L,
    # and at the <=1 year BSA (0.342 m^2) it gives CL 1.65 L/h and V 4.58 L
    # against the reported 1.46 L/h and 4.55 L.
    cl <- exp(lcl + etalcl) * (BSA / 0.81)^e_bsa_cl
    vc <- exp(lvc + etalvc) * (BSA / 0.81)^e_bsa_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, volumes in L -> central/vc has units mg/L (equivalent to ug/mL).
    # The source reports concentrations in ng/mL (1 mg/L = 1000 ng/mL) and
    # exposure in the harmonized busulfan plasma exposure unit mg*h/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
