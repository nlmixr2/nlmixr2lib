Wang_2020_caspofungin <- function() {
  description <- "Two-compartment population PK model for intravenous caspofungin in adult lung transplant recipients during the early postoperative period, with and without veno-venous extracorporeal membrane oxygenation (ECMO) (Wang 2020). Operative time enters as a power covariate (centred at 5 h) on clearance and central volume, male sex as an additive shift in central volume, and the sequential organ failure assessment (SOFA) score as a power covariate on intercompartmental clearance. ECMO was screened and not retained. Exponential IIV on CL, Vc and Vp; additive residual error. NONMEM 7.2; 19 patients, 31 profiles, 271 samples."
  reference <- "Wang Q, Zhang Z, Liu D, Chen W, Cui G, Li P, Zhang X, Li M, Zhan Q, Wang C. Population pharmacokinetics of caspofungin among extracorporeal membrane oxygenation patients during the postoperative period of lung transplantation. Antimicrob Agents Chemother. 2020;64(11):e00687-20. doi:10.1128/AAC.00687-20. PMC7577146."
  vignette <- "Wang_2020_caspofungin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Caspofungin was given as an intravenous infusion and TOTAL plasma
  # caspofungin was assayed by UPLC-MS/MS with caspofungin acetate-d4 as the
  # internal standard (Methods, 'Drug assay').
  compartmentData <- list(
    central = list(analyte = "caspofungin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "caspofungin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    T_SURG = list(
      description = "Duration of the lung transplantation operation (operative time)",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject. Enters as (T_SURG / 5)^e on clearance and on",
        "central volume, exactly as the Results equations print it: the",
        "divisor 5 h is the paper's centring value (Methods: continuous",
        "covariates are divided by their median). Group medians (Table 1):",
        "ECMO group 4.0 h (IQR 3.5-5.0), non-ECMO control group B 5.5 h",
        "(IQR 5.0-5.8). Must be strictly positive (power term)."
      ),
      source_name = "OPT"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "The paper's SEX column is a MALE indicator: the final equation is",
        "Vc = (VcTV + SEX * theta_SEX_Vc) * ..., theta = +0.62 L, and the",
        "abstract/Discussion state that 'the male gender' is associated with an",
        "increased distribution volume, so SEX = 1 for men and 0 for women.",
        "Applied in model() as (1 - SEXF) to keep the paper's positive",
        "coefficient; the typical Vc value 2.21 L is therefore the FEMALE",
        "value. Female 3/12 (25.0%) in the ECMO group and 2/7 (28.6%) in",
        "control group B (Table 1)."
      ),
      source_name = "SEX"
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score on the pharmacokinetic sampling day",
      units = "(score)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Recorded on each sampling day (Methods, 'ECMO apparatus and data",
        "collection'), so it can differ between the ECMO and post-weaning",
        "occasions of the same patient. Enters as (SOFA / 7)^1.98 on the",
        "intercompartmental clearance. The Q equation is NOT printed in the",
        "paper (only the CL and Vc equations are); the power form follows the",
        "Methods' general continuous-covariate equation (Cov_i / Cov_median)^",
        "theta, and the centring value 7 is the median SOFA score of the ECMO",
        "group and of control group B (Table 2: 7 (IQR 6-10) and 7 (6-8);",
        "control group A, the ECMO patients after weaning, 6 (5-9)). Must be",
        "strictly positive (power term)."
      ),
      source_name = "SOFA"
    )
  )

  # Screened but not retained (Methods, 'Establishment of a covariate
  # model'). Documentation only; none is referenced in model().
  covariatesDataExcluded <- list(
    ECMO_STATUS = list(
      description = "ECMO treatment-status indicator (1 = on ECMO at the sampling occasion)",
      units = "(binary)",
      type = "binary",
      notes = "Screened ('presence of ECMO'), not retained -- the paper's headline negative finding. All ECMO patients were on veno-venous ECMO; median ECMO duration 47.8 h and blood flow 2.9 L/min."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened, not retained. Median 65 (ECMO) and 59 (control B) years (Table 1)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened, not retained. Median 64 (ECMO) and 65 (control B) kg (Table 1)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened, not retained (Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened, not retained (Table 2)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened, not retained (Table 2)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened, not retained. Table 2 prints the unit as g/dl but the values (40-42.5) are g/L; no patient was hypoalbuminaemic."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened, not retained (Table 2)."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Screened, not retained. Raw mL/min, not BSA-normalised; median 81.7-93.0 mL/min by group (Table 2)."
    ),
    PCT = list(
      description = "Procalcitonin",
      units = "ug/L",
      type = "continuous",
      notes = "Screened, not retained (Table 2)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 19L,
    n_studies = 1L,
    n_samples = 271L,
    n_profiles = "31 (12 ECMO-group profiles on ECMO, 12 self-control profiles of the same patients after ECMO weaning (control group A), 7 non-ECMO profiles (control group B))",
    age_range = "Median 65 years (IQR 60-67) in the ECMO group and 59 (56-62) in control group B (Table 1)",
    weight_range = "Median 64 kg (IQR 59-69.3) in the ECMO group and 65 (53-65) in control group B (Table 1)",
    sex_female_pct = 26.3,
    race_ethnicity = "Not reported; single centre in Beijing, China",
    disease_state = "Adult lung transplant recipients (single- or double-lung) in the first days after surgery, receiving caspofungin as universal antifungal prophylaxis. Twelve required veno-venous ECMO for more than 24 h after surgery and were weaned; seven never needed ECMO. Patients on renal replacement therapy were excluded; none was hypoalbuminaemic and all had a negative 24-h fluid balance.",
    dose_range = "50 mg every 24 h by intravenous infusion, no loading dose, first dose immediately after surgery. The infusion duration is not stated.",
    regions = "China (China-Japan Friendship Hospital, Beijing; October 2017 - March 2018; NCT03766282)",
    notes = paste(
      "Sampling: pre-dose and 0.5, 1, 2, 4, 8, 12 and 24 h after the start of",
      "the infusion, at the 2nd dose (ECMO group, control group B) and at the",
      "first dose after ECMO weaning (control group A). One-, two- and",
      "three-compartment models were tested; two compartments were selected.",
      "Bootstrap (n = 1000) medians in Table 4 (CL 0.26 L/h, Vc 2.61 L, Q",
      "1.00 L/h) lie ABOVE the point estimates, and the point estimates for",
      "CL, Vc and Q fall outside their own bootstrap 95% CIs; see the vignette."
    )
  )

  ini({
    # Table 4 'Estimate' column (final model). The typical values apply at
    # T_SURG = 5 h, SOFA = 7 and female sex (SEX = 0).
    lcl <- log(0.21)
    label("Clearance at T_SURG = 5 h (L/h)") # Table 4: CL = 0.21 L/h (RSE 8.0%)
    lvc <- log(2.21)
    label("Central volume, female, at T_SURG = 5 h (L)") # Table 4: Vc = 2.21 L (RSE 5.0%)
    lvp <- log(2.87)
    label("Peripheral volume (L)") # Table 4: Vp = 2.87 L (RSE 16.0%)
    lq <- log(0.84)
    label("Intercompartmental clearance at SOFA = 7 (L/h)") # Table 4: Q = 0.84 L/h (RSE 11.0%)

    e_t_surg_cl <- 1.30
    label("Power exponent of (T_SURG / 5) on CL (unitless)") # Table 4: Theta OPT on CL = 1.30 (RSE 15.0%)
    e_t_surg_vc <- 0.93
    label("Power exponent of (T_SURG / 5) on Vc (unitless)") # Table 4: Theta OPT on VC = 0.93 (RSE 14.0%)
    e_sex_vc <- 0.62
    label("Additive increase in typical Vc for men (L)") # Table 4: Theta SEX on VC = 0.62 (RSE 23.0%); additive per the Results Vc equation
    e_sofa_q <- 1.98
    label("Power exponent of (SOFA / 7) on Q (unitless)") # Table 4: Theta SOFA on Q = 1.98 (RSE 21.0%)

    # IIV: exponential model (Methods). Table 4 values read as omega^2
    # (variance of eta); see vignette Assumptions.
    etalcl ~ 0.04 # Table 4: IIV CL = 0.04 (RSE 21.3%)
    etalvc ~ 0.01 # Table 4: IIV Vc = 0.01 (RSE 15.5%)
    etalvp ~ 0.23 # Table 4: IIV Vp = 0.23 (RSE 21.2%)

    addSd <- 0.73
    label("Additive residual error (mg/L)") # Table 4: Residual error, Additive = 0.73 mg/L
  })

  model({
    # Results, final-model equations:
    #   CL = CL_TV * (OPT / 5)^theta_OPT_CL * exp(eta_CL)
    #   Vc = (Vc_TV + SEX * theta_SEX_Vc) * (OPT / 5)^theta_OPT_Vc * exp(eta_Vc)
    # with SEX = 1 for men (= 1 - SEXF). The Q equation is not printed; the
    # Methods' power form centred at the median SOFA (7) is used.
    cl <- exp(lcl + etalcl) * (T_SURG / 5)^e_t_surg_cl
    vc <- (exp(lvc) + e_sex_vc * (1 - SEXF)) * (T_SURG / 5)^e_t_surg_vc * exp(etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq) * (SOFA / 7)^e_sofa_q

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
