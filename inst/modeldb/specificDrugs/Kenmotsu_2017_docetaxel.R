Kenmotsu_2017_docetaxel <- function() {
  description <- "Two-compartment population PK model for total (bound + unbound) docetaxel after a 1-h intravenous infusion of 60 mg/m^2 in Japanese adults with non-small cell lung cancer, including patients aged 75 years and over (Kenmotsu 2017). Clearance is proportional to body surface area and decreases with serum alpha-1 acid glycoprotein (power exponent -0.495, reference 121 mg/dL); the central volume decreases with serum albumin (power exponent -0.550, reference 3.7 g/dL). The NONMEM TRANS3 parameterisation (CL, Vc, Q, Vss) is kept, with the peripheral volume derived as Vss - Vc. Log-normal IIV on CL, Vc and Vss; exponential (proportional) residual error. The companion exposure-toxicity model for the percent decrease in absolute neutrophil count at nadir is modellib('Kenmotsu_2017_docetaxel_anc')."
  reference <- "Kenmotsu H, Imamura CK, Ono A, Omori S, Nakashima K, Wakuda K, Taira T, Naito T, Murakami H, Takahashi T, Tanigawara Y. The effects of advanced age and serum alpha1-acid glycoprotein on docetaxel unbound exposure and dose-limiting toxicity in cancer patients. Br J Clin Pharmacol. 2017;83(11):2416-2425. doi:10.1111/bcp.13354. PMCID: PMC5651322. Correction: Br J Clin Pharmacol. 2018;84(10):2449. doi:10.1111/bcp.13694 (relabels the Figure 2 y-axes as clearance; no model value changed)."
  vignette <- "Kenmotsu_2017_docetaxel"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "docetaxel", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "docetaxel", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed per subject. Clearance is directly proportional to BSA with no reference value or estimated exponent (Table 2: CL = theta1 * BSA * (AAG/121)^theta2), so theta1 is a clearance per m^2. BSA formula not stated. Cohort medians 1.54 m^2 (range 1.37-1.86, age < 75 years) and 1.64 m^2 (range 1.19-2.09, age >= 75 years), Table 1.",
      source_name = "BSA"
    ),
    AAG = list(
      description = "Serum alpha-1 acid glycoprotein concentration at baseline",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed per subject (measured before docetaxel administration). The source reports AAG in mg/dL with reference value 121 mg/dL in the Table 2 clearance equation; in canonical g/L the reference is 1.21 g/L (mg/dL / 100), and the ratio AAG/1.21 is identical to the paper's AAG/121. Cohort medians 1.18 g/L (range 0.71-2.64, age < 75 years) and 1.28 g/L (range 0.50-2.49, age >= 75 years), Table 1.",
      source_name = "AAG"
    ),
    ALB = list(
      description = "Serum albumin concentration at baseline",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed per subject (measured before docetaxel administration). The source reports ALB in g/dL with reference value 3.7 g/dL in the Table 2 central-volume equation; in canonical g/L the reference is 37 g/L (g/dL * 10), and the ratio ALB/37 is identical to the paper's ALB/3.7. Cohort median 37 g/L in both age groups (ranges 22-48 and 27-44 g/L), Table 1.",
      source_name = "ALB"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 51L,
    n_studies = 1L,
    age_range = "34-84 years",
    age_median = "63 years (age < 75 cohort, n = 34); 77 years (age >= 75 cohort, n = 17)",
    bsa_range = "1.19-2.09 m^2",
    sex_female_pct = 15.7,
    race_ethnicity = c(Asian = 100),
    disease_state = "Non-small cell lung cancer, docetaxel monotherapy",
    dose_range = "60 mg/m^2 as a 1-h intravenous infusion every 3 weeks (first course sampled), after 8 mg dexamethasone",
    regions = "Japan (Shizuoka Cancer Center; UMIN000014701)",
    aag_range = "50-264 mg/dL",
    alb_range = "2.2-4.8 g/dL",
    notes = "Prospective single-centre study enrolling 51 Japanese patients (43 male, 8 female) between August 2014 and May 2015: 34 aged < 75 years and 17 aged >= 75 years (Table 1). Five serum samples per patient during the first course (pre-dose, end of infusion, 10-60 min, 2-5 h and 12-24 h after the end of infusion). Total docetaxel by UPLC-MS/MS, calibration range 10-4000 ng/mL."
  )

  ini({
    lcl <- log(14.5)
    label("Total docetaxel clearance per m^2 body surface area at AAG = 121 mg/dL, theta1 (L/h/m^2)") # Table 2, theta1 = 14.5 (SE 0.516)
    e_aag_cl <- -0.495
    label("Power exponent of AAG/121 mg/dL on CL, theta2 (unitless)") # Table 2, theta2 = -0.495 (SE 0.0816)
    lvc <- log(8.89)
    label("Central volume at ALB = 3.7 g/dL, theta3 (L)") # Table 2, theta3 = 8.89 (SE 0.403)
    e_alb_vc <- -0.550
    label("Power exponent of ALB/3.7 g/dL on Vc, theta4 (unitless)") # Table 2, theta4 = -0.550 (SE 0.185)
    lvss <- log(129)
    label("Volume of distribution at steady state, Vss (L)") # Table 2, Vss = 129 (SE 8.28)
    lq <- log(11.0)
    label("Intercompartmental clearance, Q (L/h)") # Table 2, Q = 11.0 (SE 0.542)

    etalcl ~ 0.0426 # Table 2, omega2 CL = 0.0426 (SE 0.00903)
    etalvc ~ 0.00885 # Table 2, omega2 Vc = 0.00885 (SE 0.00923)
    etalvss ~ 0.0176 # Table 2, omega2 Vss = 0.0176 (SE 0.00987)

    propSd <- 0.18974
    label("Proportional (exponential) residual error (fraction)") # Table 2, sigma2 = 0.0360 -> sqrt(0.0360) = 0.18974; Equation 2 C = Cpred * exp(eps)
  })

  model({
    # Table 2: CL = theta1 * BSA * (AAG / 121)^theta2 with AAG in mg/dL;
    # AAG / 1.21 with AAG in g/L is the same ratio.
    cl <- exp(lcl + etalcl) * BSA * (AAG / 1.21)^e_aag_cl
    # Table 2: Vc = theta3 * (ALB / 3.7)^theta4 with ALB in g/dL;
    # ALB / 37 with ALB in g/L is the same ratio.
    vc <- exp(lvc + etalvc) * (ALB / 37)^e_alb_vc
    vss <- exp(lvss + etalvss)
    q <- exp(lq)
    # NONMEM ADVAN3 TRANS3: the peripheral volume is Vss - Vc.
    vp <- vss - vc

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d / dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L = ug/mL
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
