Liang_2023_polymyxinB <- function() {
  description <- "Two-compartment intravenous population PK model for polymyxin B in critically ill adults not receiving CRRT or ECMO (Liang 2023). Fitted by nonparametric adaptive grid (NPAG) in Pmetrics. Serum albumin is a power covariate on CL and age is a power covariate on Vc, both normalized to the cohort median (ALB 31.45 g/L, age 68 years) with exponents fixed at -0.95 and +0.95 respectively. Combined additive plus proportional residual error derived from the Pmetrics assay-error polynomial scaled by the final gamma."
  reference <- paste(
    "Liang D, Liang Z, Deng G, Cen A, Luo D, Zhang C, Ni S.",
    "Population pharmacokinetic analysis and dosing optimization of",
    "polymyxin B in critically ill patients.",
    "Front Pharmacol. 2023;14:1122310.",
    "doi:10.3389/fphar.2023.1122310. PMCID PMC10090446.",
    sep = " "
  )
  vignette <- "Liang_2023_polymyxinB"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "polymyxinB", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "polymyxinB", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    ALB = list(
      description        = "Serum albumin concentration. Power covariate on CL, normalized to the cohort median of 31.45 g/L (Table 1). Higher albumin gives lower total-drug clearance (exponent -0.95), consistent with the Discussion: 'the polymyxin B CL lowers as the level of ALB rises'.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value in the source analysis. Cohort median 31.45 g/L, range 23.1 - 41.5 g/L (Table 1). The paper reports ALB in g/L (SI), which matches the canonical register units, so no unit conversion is applied inside model(). Retained in the final model (model 6, Table 2) with a -2LL drop > 6.63 and R^2_pop improved from 0.485 to 0.573.",
      source_name        = "ALB"
    ),
    AGE = list(
      description        = "Age. Power covariate on Vc, normalized to the cohort median of 68 years (Table 1). Older patients have a larger central volume (exponent +0.95).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed. Cohort median 68 years, range 31 - 94 years (Table 1); the Monte Carlo simulations used the 5th percentile (34 y), median (68 y), and 95th percentile (93 y). Retained in the final model (model 6, Table 2); the Discussion notes this is the first report of age as a determinant of polymyxin B volume of distribution.",
      source_name        = "age"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened against volume of distribution by linear regression but excluded: 'CCR and WT were excluded because CCR had a smaller correlation coefficient with CL (R^2 = 0.01) and WT with V (R^2 = 0.01)' (Results 3.2). Cohort median 60 kg, range 50 - 80 kg (Table 1). Also used internally to compute CCR via Cockcroft-Gault (Table 1 footnote b)."
    ),
    CRCL = list(
      description = "Creatinine clearance estimated by the Cockcroft-Gault equation",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened against CL by linear regression but excluded (R^2 = 0.01; Results 3.2). Cohort median 68.29 mL/min, range 34.05 - 192.59 mL/min (Table 1). Computed as CCR (mL/min) = [(140 - age) x WT(kg)] / (72 x Scr(umol/L) / 88.4), multiplied by 0.85 for women (Table 1 footnote b). Note the source reports plain Cockcroft-Gault mL/min, not the BSA-normalized mL/min/1.73 m^2 of the canonical CRCL register entry."
    ),
    SCR = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Collected as a baseline demographic (median 75 umol/L, range 33 - 125; Table 1) and used only to derive CCR; not screened directly as a model covariate."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Recorded in Table 1 (5 of 22 female, 22.7%) and used only in the Cockcroft-Gault CCR adjustment factor of 0.85 for women; not screened as a PK model covariate."
    ),
    SCORE_APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation II score",
      units       = "points",
      type        = "continuous",
      notes       = "Screened as a power covariate on Vc in model 5 of Table 2 -- V0*(APACHE-II/21.5)^(-0.95) -- with a correlation against Vc of R^2 = 0.17, the same as age. Not retained: model 5 gave -2LL 143.7 versus 139.9 for the age model (model 4), so age was carried into the final model 6. Cohort median 21.5 points, range 15 - 46 (Table 1). SCORE_APACHE_II is not a ratified entry in inst/references/covariate-columns.md; it is used here only as a documentation label in covariatesDataExcluded, following the Mori_2018_zoledronicAcid.R precedent for screened-but-unretained covariates. It is deliberately NOT mapped to the existing SAPS_II canonical, which is a different severity score with different scoring rules."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 22,
    n_studies      = 1,
    n_observations = 64,
    age_range      = "31-94 years",
    age_median     = "68 years",
    weight_range   = "50-80 kg",
    weight_median  = "60 kg",
    sex_female_pct = 22.7,
    disease_state  = "critically ill adults (>= 18 years) with multidrug-resistant gram-negative bacterial infection; 50.0% pulmonary infection, 45.5% pulmonary infection with sepsis/septicemia, 4.5% intracranial infection. Pathogens: Acinetobacter baumannii 50.0%, Pseudomonas aeruginosa 13.6%, Klebsiella pneumoniae 9.1%, co-infection with 2 or more organisms 27.3%",
    renal_function = "baseline serum creatinine median 75 umol/L (33-125); baseline Cockcroft-Gault creatinine clearance median 68.29 mL/min (34.05-192.59)",
    albumin        = "baseline serum albumin median 31.45 g/L (23.1-41.5)",
    severity       = "APACHE-II score median 21.5 points (15-46)",
    dose_range     = "intravenous infusion over 1-2 h; loading dose 100-150 mg, maintenance dose 50-75 mg q12h",
    regions        = "China (single center, Guangzhou)",
    notes          = "Baseline demographics from Table 1. Patients receiving continuous renal replacement therapy (n = 10) or extracorporeal membrane oxygenation (n = 2) at the time of blood collection were excluded, so the model must not be extrapolated to CRRT or ECMO patients (Discussion, study limitations). Sampling was at steady state, >= 48 h after the start of treatment: 24 trough samples (within 1 h before a dose), 23 peak samples (within 1 h after the end of infusion), and 17 samples 6-8 h after infusion. Total (not unbound) polymyxin B was measured as the sum of polymyxin B1, B2, and B1-Ile by LC-MS/MS."
  )

  ini({
    # Structural parameters. Values are the Table 3 MEDIAN column of the final
    # model 6 nonparametric parameter distribution, which is the quantity that
    # maps onto exp(l<param>) under the log-normal IIV encoding used here. The
    # Table 3 mean +/- SD values quoted in the Abstract and Discussion are given
    # alongside each line for reference. See the vignette "Assumptions and
    # deviations" section for the full rationale.
    lcl <- log(1.22)   ; label("Clearance (L/h)")                     # Table 3, median (mean 1.24, SD 0.38)
    lvc <- log(13.52)  ; label("Central volume of distribution (L)")  # Table 3, median (mean 16.64, SD 12.74)
    lq  <- log(2.34)   ; label("Intercompartmental clearance (L/h)")  # Table 3, median (mean 3.04, SD 2.27)
    lvp <- log(76.53)  ; label("Peripheral volume of distribution (L)") # Table 3, median (mean 66.20, SD 36.25)

    # Covariate effects, both power form normalized to the cohort median.
    # Wrapped in fixed() because the exponents are printed as literal constants
    # in the Table 2 model definitions and no estimate, SE, RSE, or CI is
    # reported for either. Confirmed by the Table 2 information-criterion
    # arithmetic: AIC - (-2LL) is 11.0-11.1 for every two-compartment model
    # (models 2-6) whether zero, one, or two covariate relationships are
    # present, so adding a covariate did not add an estimated parameter. Only
    # CL, Q, Vc, and Vp carry estimates in Table 3.
    e_alb_cl <- fixed(-0.95) ; label("Exponent of (ALB/31.45) on CL (unitless)")  # Table 2, model 6: CL0*(ALB/31.45)^(-0.95)
    e_age_vc <- fixed(0.95)  ; label("Exponent of (AGE/68) on Vc (unitless)")     # Table 2, model 6: V0*(age/68)^0.95

    # Between-subject variability. Table 3 reports BSV as %CV of the
    # nonparametric parameter distribution (equal to SD/mean of the Table 3
    # Mean and SD columns). Converted to log-normal variance via
    # omega^2 = log(CV^2 + 1). No parameter correlations are reported, so the
    # etas are diagonal.
    etalcl ~ 0.088837  # Table 3, BSV 30.48 %CV
    etalvc ~ 0.461209  # Table 3, BSV 76.55 %CV
    etalq  ~ 0.444176  # Table 3, BSV 74.78 %CV
    etalvp ~ 0.262261  # Table 3, BSV 54.76 %CV

    # Residual unexplained variability. Pmetrics weights each observation by
    # the reciprocal of SD x gamma, where the assay-error polynomial is
    # SD = C0 + C1*Cobs + C2*Cobs^2 + C3*Cobs^3 with C0 = 0.1, C1 = 0.15,
    # C2 = C3 = 0 (Methods 2.4), and the final-cycle gamma was 0.72
    # (Results 3.2). The total residual SD is therefore
    # 0.72 * (0.1 + 0.15*C) = 0.072 + 0.108*C mg/L, i.e. a combined additive
    # plus proportional error. Both terms are fixed() because C0, C1, and
    # gamma are reported without uncertainty and C0/C1 were set a priori
    # ("scientifically set to 0.1, 0.15, 0, and 0").
    addSd  <- fixed(0.072) ; label("Additive residual error (mg/L)")       # Methods 2.4 C0 = 0.1 x Results 3.2 gamma = 0.72
    propSd <- fixed(0.108) ; label("Proportional residual error (fraction)") # Methods 2.4 C1 = 0.15 x Results 3.2 gamma = 0.72
  })

  model({
    # Individual parameters. Power covariate models normalized to the cohort
    # medians reported in Table 1 (ALB 31.45 g/L, age 68 years).
    cl <- exp(lcl + etalcl) * (ALB / 31.45)^e_alb_cl
    vc <- exp(lvc + etalvc) * (AGE / 68)^e_age_vc
    q  <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ODE system. Reproduces Eqs 1 and 2 of the source:
    #   dX1/dt = RateIV + (Q/Vp)*X2 - (CL/Vc + Q/Vc)*X1
    #   dX2/dt = (Q/Vc)*X1 - (Q/Vp)*X2
    # RateIV (the intravenous infusion rate) enters via the event table.
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Observation and error
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
