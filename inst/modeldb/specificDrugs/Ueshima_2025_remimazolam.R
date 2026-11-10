Ueshima_2025_remimazolam <- function() {
  description <- "Joint population PK/PD model for remimazolam, its inactive metabolite CNS 7054, and the bispectral index (BIS) in older Japanese adults during the first minutes of remimazolam-only induction of general anesthesia (Ueshima 2025). Remimazolam is described by a three-compartment model; the whole of remimazolam clearance feeds a single transit compartment that delays the appearance of CNS 7054, which is described by a two-compartment model. An effect compartment driven by the remimazolam plasma concentration feeds a maximum-inhibition (Emax) model on BIS with a fixed baseline of 100. Body weight scales remimazolam clearance allometrically (exponent fixed at 0.75, reference 60.4 kg) and the literature-fixed CNS 7054 peripheral volume and intercompartmental clearance linearly per 75 kg; no covariate was retained on the PD."
  reference <- "Ueshima S, Nagahara S, Nakaya K, Kawamura M, Hira D, Matsuura Y, Toyama H, Kawamae K, Kakumoto M. Pharmacokinetic and pharmacodynamic characterization of remimazolam in older Japanese adults who underwent general anesthesia at early stage of infusion. Sci Rep. 2025;15:24979. doi:10.1038/s41598-025-10015-z"
  vignette <- "Ueshima_2025_remimazolam"
  units <- list(time = "min", dosing = "mg", concentration = "ng/mL")

  # What each ODE state holds, verified against Ueshima 2025 Figure 1 (model
  # schematic): V1 remimazolam central with V2 and V3 peripheral, CL from V1
  # into "Transit", kT from Transit into V1M (CNS 7054 central) with V2M
  # peripheral, and the effect compartment fed from V1 by kE0.
  compartmentData <- list(
    central = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    transit1 = list(
      analyte = "remimazolam-derived CNS 7054 precursor",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central_cns7054 = list(analyte = "CNS 7054", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_cns7054 = list(analyte = "CNS 7054", units = "mg", specimen = "plasma", verified = TRUE),
    effect = list(analyte = "remimazolam", units = "ng/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight at baseline.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Two weight terms. (1) Remimazolam CL = theta_CL * (BW/60.4)^0.75 (Ueshima 2025 Table 2 row equation); 60.4 kg is the cohort median (Table 1) and the 0.75 exponent is the fixed allometric value of Methods Equation 3. (2) CNS 7054 V2M = theta_V2M * (BW/75) and Q2M = theta_Q2M * (BW/75) (Table 3 row equations, units 'L/75 kg' and 'L/min/75 kg'), carried with their fixed literature values from Schuttler 2020. No other parameter is weight-scaled in the printed final-model equations.",
      source_name = "BW"
    )
  )

  # Screened by stepwise forward inclusion / backward elimination on the PK
  # and PD parameters (Methods Equations 3-4) and not retained in either
  # final model (Results, 'Population pharmacokinetic modeling' and
  # 'Population pharmacodynamic modeling'). Cohort summaries are Table 1.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline.",
      units = "years",
      type = "continuous",
      notes = "Screened as a median-normalised power function on PK and PD parameters. On EC50 it lowered the OFV by 4.74 in forward inclusion but did not survive backward elimination. Cohort median 71.7 years, range 42.2-88.6."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      notes = "Screened as a dichotomous covariate (Methods Equation 4); not retained. Cohort 12 male / 8 female."
    ),
    BMI = list(
      description = "Body mass index at baseline.",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened; not retained. Cohort median 23.1 kg/m^2, range 18.1-28.0."
    ),
    CREAT = list(
      description = "Serum creatinine.",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened; not retained. Cohort median 0.81 mg/dL, range 0.50-1.28."
    ),
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation.",
      units = "mL/min",
      type = "continuous",
      notes = "Screened; not retained. Raw Cockcroft-Gault mL/min, NOT BSA-normalised (the canonical CRCL unit is mL/min/1.73 m^2). Cohort median 72.6 mL/min, range 41.6-115.6."
    ),
    AST = list(
      description = "Aspartate aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained. Cohort median 21 IU/L, range 15-48."
    ),
    ALT = list(
      description = "Alanine aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained. Cohort median 21 IU/L, range 8-57."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20,
    n_studies = 1,
    age_range = "42.2-88.6 years",
    age_median = "71.7 years",
    weight_range = "44.5-77.1 kg",
    weight_median = "60.4 kg",
    bmi_range = "18.1-28.0 kg/m2",
    sex_female_pct = 40,
    race_ethnicity = c(Asian = 100),
    disease_state = "Adult surgical patients undergoing general anesthesia induced with remimazolam alone (no remifentanil or fentanyl before the last sample); 14 of 20 (70%) aged 65 years or older. Severe liver cirrhosis (Child-Pugh C) excluded.",
    dose_range = "Remimazolam IV infusion at 12 mg/kg/h (0.2 mg/kg/min) until loss of consciousness (median 1.71 min, range 1.20-3.00), then 1.0 mg/kg/h (0.017 mg/kg/min).",
    regions = "Japan (Yamagata University Hospital, April 2021 - October 2022)",
    n_observations = "60 remimazolam and 60 CNS 7054 plasma concentrations and 60 BIS values, sampled at 1, 3 and 5 min after loss of consciousness.",
    notes = "Demographics from Ueshima 2025 Table 1. Sequential PK then PD analysis in NONMEM 7.3 (FOCE-I); the PD model was fitted on individual Bayesian PK parameters."
  )

  ini({
    # -----------------------------------------------------------------
    # Remimazolam (parent) disposition. Ueshima 2025 Table 2.
    # -----------------------------------------------------------------
    lcl <- log(1.38); label("Remimazolam clearance at 60.4 kg, CL (L/min)") # Table 2: theta_CL = 1.38 L/min (95% CI 1.12-1.62; RSE 9.49%)
    lvc <- log(6.38); label("Remimazolam central volume, V1 (L)") # Table 2: theta_V1 = 6.38 L (95% CI 4.28-8.48; RSE 16.8%)
    lvp <- log(4.11); label("Remimazolam peripheral volume, V2 (L)") # Table 2: theta_V2 = 4.11 L (95% CI 2.40-5.82; RSE 21.3%)
    lvp2 <- fixed(log(15.5)); label("Remimazolam second peripheral volume, V3 (L)") # Table 2: theta_V3 = 15.5 FIX, literature value (Schuttler 2020)
    lq <- fixed(log(1.04)); label("Remimazolam intercompartmental clearance V1-V2, Q2 (L/min)") # Table 2: theta_Q2 = 1.04 FIX, literature value (Schuttler 2020)
    lq2 <- fixed(log(0.19)); label("Remimazolam intercompartmental clearance V1-V3, Q3 (L/min)") # Table 2: theta_Q3 = 0.19 FIX, literature value (Schuttler 2020)

    # Allometric exponent on CL, Methods Equation 3 ('set to 0.75' for
    # clearance) and the Table 2 row equation CL = theta_CL * (BW/60.4)^0.75.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on remimazolam CL (unitless)") # Table 2 row equation: CL = theta_CL * (BW/60.4)^0.75

    # -----------------------------------------------------------------
    # CNS 7054 formation and disposition. Ueshima 2025 Table 3.
    # -----------------------------------------------------------------
    lktr <- fixed(log(0.024)); label("Transit rate constant for CNS 7054 formation, kT (1/min)") # Table 3: theta_kT = 0.024 FIX, literature value (Schuttler 2020)
    lcl_cns7054 <- log(0.123); label("CNS 7054 clearance, CLM (L/min)") # Table 3: theta_CLM = 0.123 L/min (95% CI 0.0624-0.184; RSE 25.1%)
    lvc_cns7054 <- log(0.423); label("CNS 7054 central volume, V1M (L)") # Table 3: theta_V1M = 0.423 L (95% CI 0.202-0.644; RSE 26.7%)
    lvp_cns7054 <- fixed(log(7.4)); label("CNS 7054 peripheral volume per 75 kg, V2M (L)") # Table 3: theta_V2M = 7.4 L/75 kg FIX, literature value (Schuttler 2020)
    lq_cns7054 <- fixed(log(0.13)); label("CNS 7054 intercompartmental clearance per 75 kg, Q2M (L/min)") # Table 3: theta_Q2M = 0.13 L/min/75 kg FIX, literature value (Schuttler 2020)

    # -----------------------------------------------------------------
    # Pharmacodynamics. Ueshima 2025 Table 4 and Methods Equation 5:
    #   BIS = 100 - Emax * Ce / (EC50 + Ce)
    # -----------------------------------------------------------------
    lec50 <- log(79.6); label("Effect-site remimazolam concentration at half-maximal BIS reduction, EC50 (ng/mL)") # Table 4: theta_EC50 = 79.6 ng/mL (95% CI 5.7-153; RSE 47.4%)
    limax <- log(54.0); label("Maximum reduction in BIS from the baseline of 100, Emax (BIS units)") # Table 4: theta_Emax = 54.0 (95% CI 47.4-60.6; RSE 6.22%)
    lke0 <- fixed(log(0.135)); label("Effect-compartment equilibration rate constant, kE0 (1/min)") # Table 4: theta_ke0 = 0.135 FIX, literature value (Zhou 2020)

    # -----------------------------------------------------------------
    # Inter-individual variability, exponential model (Methods Equation 2).
    # The tabulated percentages are read as CVs and converted with
    # omega^2 = log(CV^2 + 1); see the vignette 'Assumptions and deviations'.
    # -----------------------------------------------------------------
    etalcl ~ 0.0262 # Table 2: eta_CL = 16.3% (95% CI 7.95-24.6; RSE 26.1%); log(0.163^2 + 1)
    etalvc ~ 0.5923 # Table 2: eta_V1 = 89.9% (95% CI 39.9-140; RSE 28.4%); log(0.899^2 + 1)
    etalcl_cns7054 ~ 0.3680 # Table 3: eta_CLM = 66.7% (95% CI 34.6-98.8; RSE 24.6%); log(0.667^2 + 1)
    etalvc_cns7054 ~ 0.2865 # Table 3: eta_V1M = 57.6% (95% CI 36.6-78.6; RSE 18.6%); log(0.576^2 + 1)
    etalec50 ~ 0.6791 # Table 4: eta_EC50 = 98.6% (95% CI 29.8-167; RSE 35.6%); log(0.986^2 + 1); Results prose says 35.1%, see vignette

    # -----------------------------------------------------------------
    # Residual variability, exponential model (Methods Equation 1), encoded
    # as a log-normal error with the tabulated percentage as the log-scale SD.
    # -----------------------------------------------------------------
    expSd <- 0.0879; label("Exponential residual error on remimazolam (log-scale SD)") # Table 2: epsilon_R = 8.79% (95% CI 6.24-11.3; RSE 14.8%)
    expSd_cns7054 <- 0.0918; label("Exponential residual error on CNS 7054 (log-scale SD)") # Table 3: epsilon_C = 9.18% (95% CI 1.28-17.1; RSE 43.9%)
    expSd_BIS <- 0.138; label("Exponential residual error on BIS (log-scale SD)") # Table 4: epsilon_B = 13.8% (95% CI 10.7-16.9; RSE 11.6%)
  })

  model({
    # Reference weights printed in the final-model row equations: the cohort
    # median for remimazolam CL (Table 2) and the 75 kg standard of the
    # literature-fixed CNS 7054 peripheral parameters (Table 3).
    ref_wt_cl <- 60.4
    ref_wt_cns7054 <- 75

    # Amounts in mg and volumes in L give mg/L; the paper reports
    # concentrations and EC50 in ng/mL.
    mgL_to_ngmL <- 1000

    # Baseline BIS, the constant of Methods Equation 5 (not estimated).
    bis_base <- 100

    # Remimazolam (parent) individual parameters.
    cl <- exp(lcl + etalcl) * (WT / ref_wt_cl)^e_wt_cl
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    vp2 <- exp(lvp2)
    q <- exp(lq)
    q2 <- exp(lq2)
    ktr <- exp(lktr)

    # CNS 7054 individual parameters; V2M and Q2M scale linearly with weight
    # per the Table 3 row equations.
    cl_cns7054 <- exp(lcl_cns7054 + etalcl_cns7054)
    vc_cns7054 <- exp(lvc_cns7054 + etalvc_cns7054)
    vp_cns7054 <- exp(lvp_cns7054) * (WT / ref_wt_cns7054)
    q_cns7054 <- exp(lq_cns7054) * (WT / ref_wt_cns7054)

    # Pharmacodynamic individual parameters.
    ec50 <- exp(lec50 + etalec50)
    imax <- exp(limax)
    ke0 <- exp(lke0)

    # Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2
    kelm <- cl_cns7054 / vc_cns7054
    k56 <- q_cns7054 / vc_cns7054
    k65 <- q_cns7054 / vp_cns7054

    # Three-compartment remimazolam disposition with IV input to V1.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # All remimazolam clearance passes through one transit compartment into
    # CNS 7054 central (Figure 1). No molecular-weight conversion is stated,
    # so the transfer is 1:1 in mass.
    d/dt(transit1) <- kel * central - ktr * transit1
    d/dt(central_cns7054) <- ktr * transit1 - kelm * central_cns7054 -
      k56 * central_cns7054 + k65 * peripheral1_cns7054
    d/dt(peripheral1_cns7054) <- k56 * central_cns7054 - k65 * peripheral1_cns7054

    # Plasma concentrations (ng/mL).
    Cc <- central / vc * mgL_to_ngmL
    Cc_cns7054 <- central_cns7054 / vc_cns7054 * mgL_to_ngmL

    # Effect compartment (Figure 1, kE0 in and out); the state is the
    # effect-site concentration in ng/mL.
    d/dt(effect) <- ke0 * (Cc - effect)

    # Maximum-inhibition model on BIS, Methods Equation 5.
    BIS <- bis_base - imax * effect / (ec50 + effect)

    Cc ~ lnorm(expSd)
    Cc_cns7054 ~ lnorm(expSd_cns7054)
    BIS ~ lnorm(expSd_BIS)
  })
}
