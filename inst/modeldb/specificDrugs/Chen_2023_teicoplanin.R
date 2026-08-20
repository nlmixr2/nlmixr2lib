Chen_2023_teicoplanin <- function() {
  description <- "Two-compartment IV infusion population PK model for teicoplanin in critically ill adults with sepsis in the intensive care unit, with CKD-EPI estimated glomerular filtration rate as a power covariate on clearance (Chen 2023)"
  reference   <- "Chen CY, Xie M, Gong J, Yu N, Wei R, Lei LL, Zhao SM, Li RM, Dong X, Zhang XL, Zhou Y, Li SL, Cui YM. Population pharmacokinetic analysis and dosing regimen optimization of teicoplanin in critically ill patients with sepsis. Front Pharmacol. 2023;14:1132367. doi:10.3389/fphar.2023.1132367"
  vignette    <- "Chen_2023_teicoplanin"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "teicoplanin", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "teicoplanin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Body-surface-area-normalized glomerular filtration rate estimated with the CKD-EPI equation (Levey 2009), the only covariate retained in the final model (on clearance)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column GFR. Creatinine-based CKD-EPI estimate, already normalized to 1.73 m^2 body surface area, stored under the canonical CRCL column per inst/references/covariate-columns.md. Entered the final model as the power term (GFR/71.88)^0.437 on CL (Chen 2023 Eq. 6). The centering value 71.88 mL/min/1.73 m^2 is the study-population median as printed in Eq. 6; Table 1 reports the same median rounded to 71.9 (range 11.0-124). Renal-function distribution of the 59 subjects (Chen 2023 Results 3.1): 23.7% GFR >= 90, 33.9% GFR 60-90, 23.7% GFR 30-60, 8.5% GFR 15-30, 10.2% GFR < 15 mL/min/1.73 m^2. Patients receiving renal replacement therapy were excluded by protocol, so the covariate is not defined for dialysis records.",
      source_name        = "GFR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 59L,
    n_studies      = 1L,
    age_range      = "28.0-92.0 years (Chen 2023 Table 1)",
    age_median     = "72.0 years",
    weight_range   = "35.0-90.0 kg (Chen 2023 Table 1)",
    weight_median  = "65.0 kg",
    sex_female_pct = 37.3,
    race_ethnicity = "Not reported by category; single-centre Chinese ICU cohort described in the Conclusion as a 'prospective cohort of Chinese septic patients'",
    disease_state  = "Critically ill adults (age >= 18 years) admitted to the intensive care unit with sepsis diagnosed by the Sepsis-3.0 criteria (Singer 2016), with confirmed or suspected gram-positive coccal infection and an expected teicoplanin course of >= 4 days. Excluded: pregnant and lactating patients, infections at special sites (endocardium, bone, joints), and any treatment that may affect drug elimination such as renal replacement therapy.",
    dose_range     = "Teicoplanin 400 mg IV infusion every 12 h for 3 doses (loading), followed by 400 mg once daily, with subsequent adjustment per clinical response. Infusion duration 30-60 min. The published Monte Carlo simulations (Chen 2023 Table 3) explored loading doses of 10-15 mg/kg q12h for 3 or 5 doses and maintenance doses of 3-15 mg/kg every 24-72 h in a typical 65 kg patient.",
    regions        = "China (Peking University First Hospital ICU, Beijing)",
    renal_function = "GFR (CKD-EPI) median 71.9 mL/min/1.73 m^2, range 11.0-124 (Chen 2023 Table 1). Serum creatinine median 80.3 umol/L, range 28.4-429. Patients on renal replacement therapy were excluded.",
    hepatic_function = "ALT median 15.0 IU/L (range 4.00-993), AST median 26.0 IU/L (range 11.0-2390), total bilirubin median 18.7 umol/L (range 5.50-244), direct bilirubin median 6.80 umol/L (range 0.600-141) (Chen 2023 Table 1).",
    screened_covariates = "Tested-but-not-retained covariates (Chen 2023 Methods 2.4 and Results 3.2): sex, age, height, weight, BMI, ALT, AST, total bilirubin, direct bilirubin, total protein, albumin, white blood count, serum creatinine and blood urea nitrogen. GFR on CL was the only covariate retained, decreasing the objective function value by 17.95 and reducing omega_CL by 9.8%.",
    notes          = "Prospective open-label PPK study; 249 serum concentrations from 59 subjects (22 female, 37 male). Sampling: trough samples immediately before the 3rd, 4th, 5th and 6th doses in every patient, plus one extra sample per patient assigned by thirds (immediately after the 5th infusion, 1 h after the 5th dose, or 1 h before the 6th dose). Assay: HPLC with UV detection at 240 nm (Agilent 1100); LOQ 3.125 ug/mL, linear 3.125-100 ug/mL (r = 0.9994); intra- and inter-day precision RSD < 10%. Estimation: NONMEM 7.4, FOCE with eta-epsilon interaction, ADVAN3 TRANS4. Model qualified by goodness-of-fit, 1000-sample bootstrap (923 successful runs) and prediction-corrected VPC. Additional baseline laboratory values (Chen 2023 Table 1): total protein median 52.3 g/L (33.5-71.0), albumin median 29.5 g/L (15.5-39.4), white blood count median 11.5 x10^9/L (2.39-25.9), blood urea nitrogen median 10.5 mmol/L (1.09-81.6), BMI median 24.1 kg/m^2 (15.6-32.0), height median 1.67 m (1.33-1.80)."
  )

  ini({
    # Structural fixed-effect parameters from Chen 2023 Table 2 ("Estimate (RSE%)"
    # column of the final model) and the printed final-model equations 6-9.
    lcl       <- log(1.03); label("Clearance (L/h)")                        # Chen 2023 Table 2: CL = 1.03 L/h (RSE 16.6%); Eq. 6
    lvc       <- log(20.1); label("Central volume V1 (L)")                  # Chen 2023 Table 2: V1 = 20.1 L (RSE 12.9%); Eq. 7
    lq        <- log(3.12); label("Intercompartmental clearance Q (L/h)")   # Chen 2023 Table 2: Q = 3.12 L/h (RSE 10.9%); Eq. 8
    lvp       <- log(101);  label("Peripheral volume V2 (L)")               # Chen 2023 Table 2: V2 = 101 L (RSE 12.7%); Eq. 9

    # Covariate effect: power model on BSA-normalized CKD-EPI GFR, centered at the
    # study-population median. Chen 2023 Eq. 6 prints
    #   CL (L/h) = 1.03 * (GFR / 71.88)^0.437 * e^0.29
    # so the centering constant 71.88 mL/min/1.73 m^2 comes from the equation
    # (Table 1 reports the same median rounded to 71.9).
    e_crcl_cl <- 0.437; label("Power exponent on (CRCL / 71.88) for CL (unitless)") # Chen 2023 Table 2: theta_CL_GFR = 0.437 (RSE 23.8%); Eq. 6

    # Between-subject variability. Chen 2023 Table 2 reports interindividual
    # variability as a percentage, and Eqs. 6-9 print the corresponding exponential
    # IIV terms as e^0.29 (CL), e^0.37 (V1), e^0.29 (Q) and e^0.10 (V2). Squaring the
    # Table 2 percentages reproduces those printed exponents exactly
    # (0.539^2 = 0.290, 0.607^2 = 0.368, 0.543^2 = 0.295, 0.321^2 = 0.103), which
    # establishes that the published percentages are omega (the SD on the log scale)
    # x 100, i.e. the usual NONMEM approximate-CV convention. The more precise
    # three-significant-figure Table 2 values are used here rather than the two-decimal
    # exponents printed in the equations. No off-diagonal covariances were published.
    etalcl ~ 0.290  # Chen 2023 Table 2: IIV on CL = 53.9% (RSE 10.1%) -> omega^2 = 0.539^2 = 0.290; Eq. 6 prints e^0.29
    etalvc ~ 0.368  # Chen 2023 Table 2: IIV on V1 = 60.7% (RSE 23.9%) -> omega^2 = 0.607^2 = 0.368; Eq. 7 prints e^0.37
    etalq  ~ 0.295  # Chen 2023 Table 2: IIV on Q  = 54.3% (RSE 15.5%) -> omega^2 = 0.543^2 = 0.295; Eq. 8 prints e^0.29
    etalvp ~ 0.103  # Chen 2023 Table 2: IIV on V2 = 32.1% (RSE 16.0%) -> omega^2 = 0.321^2 = 0.103; Eq. 9 prints e^0.10

    # Proportional residual error only (Chen 2023 Results 3.2: "a two-compartment model
    # with proportional residual best described the data"; Eq. 5 C_ij = C_pred,ij * (1 + eps_prop,ij)).
    propSd <- 0.174; label("Proportional residual error (fraction)")        # Chen 2023 Table 2: proportional error = 17.4%, encoded here as a fraction
  })
  model({
    # Individual PK parameters. GFR enters as a power term on CL only; V1, Q and V2
    # carry log-normal IIV with no covariate effects (Chen 2023 Eqs. 6-9).
    cl <- exp(lcl + etalcl) * (CRCL / 71.88)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # Micro-constants for the explicit two-compartment ODE (NONMEM ADVAN3 TRANS4).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # IV infusion into the central compartment (30-60 min infusions in the study);
    # first-order distribution to peripheral1 and first-order elimination from central.
    # Dose in mg and volumes in L give central / vc in mg/L (equivalently ug/mL).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
