Xiang_2025_tacrolimus_fpg <- function() {
  description <- "Sequential exposure-response (PK/linear) model relating oral tacrolimus trough concentration to fasting plasma glucose in adult renal transplant recipients, as a quantitative model of post-transplantation diabetes mellitus risk. Fasting plasma glucose rises linearly from its pre-transplantation baseline with the predicted tacrolimus concentration; the baseline rises as a power function of age. The tacrolimus PK is the one-compartment model of the same paper, carried here as fixed parameters because the PK/PD model was fitted sequentially on the population PK step; see Xiang_2025_tacrolimus. The companion renal-function model in the same paper is Xiang_2025_tacrolimus_egfr."
  reference <- paste(
    "Xiang Q, Yang Y, Li G, Chen S, Yang Y, Liu L, Yu X.",
    "Population Pharmacokinetic/Pharmacodynamic Modeling of Tacrolimus in Renal",
    "Transplant Recipients: Impact of CYP3A5 Genotype and Wuzhi Capsule",
    "Co-Medication. Drug Des Devel Ther. 2025;19:8375-8389.",
    "doi:10.2147/DDDT.S542786.",
    "Parameter values are the final-model estimates in Table 2 ('C0-FPG model'",
    "block, plus the 'PPK model' block for the fixed PK); the structural",
    "equations are Eq. 1, 2, 5, 6 and 7.",
    sep = " "
  )
  vignette <- "Xiang_2025_tacrolimus"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  # Issue #482: the ODE states are the tacrolimus PK carried over from the
  # population PK step. Fasting plasma glucose is an algebraic direct effect of
  # the predicted concentration and is not a state.
  compartmentData <- list(
    depot   = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    POD = list(
      description        = "Post-operative day: days elapsed since renal transplantation",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = "TIME-VARYING within subject. Enters the fixed PK layer as the power scaling (POD / 34)^0.109 on CL/F; 34 days is the median post-operative day of the index cohort (Xiang 2025 Table 1). The power form is undefined at POD = 0, so supply POD >= 1. It has no direct effect on fasting plasma glucose; it acts only through the predicted tacrolimus concentration.",
      source_name        = "POD"
    ),
    CONMED_WUZHI = list(
      description        = "Concomitant Wuzhi capsule (Schisandra sphenanthera extract) indicator (1 = receiving Wuzhi capsule, 0 = not)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant Wuzhi capsule)",
      notes              = "Source column WZ, defined below Xiang 2025 Eq. 6 as 'WZ = 1 if WZ is present; otherwise = 0'. Enters the fixed PK layer as exp(-0.211 * CONMED_WUZHI) on CL/F. It has no direct effect on fasting plasma glucose; it acts only through the predicted tacrolimus concentration, which is why the Wuzhi capsule raises simulated FPG.",
      source_name        = "WZ"
    ),
    CYP3A5_EXPR = list(
      description        = "CYP3A5 expresser status (1 = carries at least one CYP3A5*1 allele, genotype *1/*1 or *1/*3; 0 = CYP3A5*3/*3 non-expresser)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5*3/*3 non-expresser)",
      notes              = "VALUE INVERSION relative to the source: Xiang 2025 codes 'Genotype = 1 if the genotype is CYP3A5*3/*3' below Eq. 6, so CYP3A5_EXPR = 1 - Genotype and the model applies the published coefficient to (1 - CYP3A5_EXPR). Enters the fixed PK layer only; it has no direct effect on fasting plasma glucose.",
      source_name        = "Genotype"
    ),
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed. Enters as the power scaling (AGE / 38)^0.232 on the pre-transplantation baseline fasting plasma glucose (Xiang 2025 Eq. 7); 38 years is the median age of the index cohort. This was the only covariate retained on FPG0 (dOFV = 18.765, df = 1, p < 0.001). Amlodipine co-medication was significant in forward inclusion (dOFV = -8.613) but failed backward elimination and is recorded in covariatesDataExcluded. Age does not modify the tacrolimus effect on FPG, only its baseline.",
      source_name        = "Age"
    )
  )

  # Screened in the paper's covariate analysis but NOT retained in the final
  # model, so documented rather than used (Xiang 2025 Discussion and
  # Supplementary Table 3): amlodipine reached the forward-inclusion threshold
  # on FPG0 but did not meet the backward-elimination criterion (p < 0.001).
  covariatesDataExcluded <- list(
    CONMED_AMLODIPINE = list(
      description        = "Concomitant amlodipine indicator (1 = receiving amlodipine, 0 = not)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant amlodipine)",
      notes              = "Screened on baseline fasting plasma glucose and significant in forward inclusion only (dOFV = -8.613, p < 0.01; Xiang 2025 Discussion and Supplementary Table 3). It did not meet the backward-elimination criterion (p < 0.001) and is absent from the final model, and no point estimate is published for it. 73% of the index cohort received amlodipine."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 100L,
    n_studies      = 1L,
    age_range      = "19-65 years (median 38)",
    age_median     = "38 years",
    weight_range   = "31-99 kg (median 60)",
    weight_median  = "60 kg",
    sex_female_pct = 36.0,
    race_ethnicity = c(Han = 92.0, Other = 8.0),
    disease_state  = "Adult renal transplant recipients on a tacrolimus-based immunosuppressive regimen; baseline fasting plasma glucose 1.45-20.79 mmol/L (mean 6.01, median 5.66)",
    dose_range     = "Oral tacrolimus started within 24 h of transplantation at 0.05-0.15 mg/kg/day divided q12h, then dose-adjusted to the trough TDM target; observed doses 0.5-5.5 mg per administration (median 2.5)",
    regions        = "Single centre, Chongqing, People's Republic of China",
    notes          = "Retrospective cohort of 126 renal transplant recipients randomly split 4:1 into an index group (n = 100) used to build the model and an external validation group (n = 26). 2055 fasting plasma glucose observations. Baseline demographics from Xiang 2025 Table 1 (index dataset column). Estimated in Phoenix NLME 8.3.5, sequentially on the population PK step. A fasting plasma glucose above 7.0 mmol/L is the paper's threshold for likely post-transplantation diabetes mellitus."
  )

  ini({
    # -------------------------------------------------------------------------
    # Tacrolimus PK, fixed from the population PK step of the same paper
    # (Xiang 2025 Eq. 5 and 6, Table 2 'PPK model' block). Methods: "The PK/PD
    # models were sequentially built following the population PK model", and
    # Table 2 does not re-report PK parameters in the C0-FPG block, so none of
    # these were re-estimated here. See modellib('Xiang_2025_tacrolimus') for
    # the estimated form.
    # -------------------------------------------------------------------------
    lka <- fixed(log(3.09)); label("Absorption rate constant Ka (1/h)")                 # Xiang 2025 Table 2 PPK model: Ka = 3.09 h^-1 FIX
    lvc <- fixed(log(2248)); label("Apparent volume of distribution V/F (L)")           # Xiang 2025 Eq. 5 and Table 2: V/F = 2248 L
    lcl <- fixed(log(33.7)); label("Apparent clearance CL/F at POD 34 d, CYP3A5 expresser, no Wuzhi capsule (L/h)")  # Xiang 2025 Eq. 6 and Table 2: CL/F = 33.7 L/h

    e_pod_cl          <- fixed(0.109);  label("Power exponent of (POD / 34 days) on CL/F (unitless)")                                  # Xiang 2025 Table 2: POD on CL/F = 0.109
    e_conmed_wuzhi_cl <- fixed(-0.211); label("Exponential effect of concomitant Wuzhi capsule on CL/F (unitless)")                    # Xiang 2025 Table 2: WZ on CL/F = -0.211
    e_cyp3a5_expr_cl  <- fixed(-0.381); label("Exponential effect of CYP3A5 non-expresser status on CL/F, applied to (1 - CYP3A5_EXPR) (unitless)")  # Xiang 2025 Table 2: CYP3A5*3/*3 on CL/F = -0.381

    # -------------------------------------------------------------------------
    # Fasting plasma glucose. Xiang 2025 Eq. 2:  FPG = FPG0 + beta * Cp,t
    # and Eq. 7:  FPG0 = 4.2 * (Age/38)^0.232 * e^(etaFPG0)
    #
    # beta was FIXED rather than estimated. Results: "In the C0-FPG model, beta
    # was fixed at 0.264"; Discussion: "the cumulative incidence rate of PTDM
    # within six months following renal transplantation (0.264) was set as the
    # slope beta" (the authors' reference 39). None of the direct-effect,
    # indirect-effect or linear models could estimate all PD parameters
    # precisely at once, so beta was fixed to avoid overparameterisation, and a
    # +/-30% sensitivity analysis changed effect values by < 15%
    # (Supplementary Table 5).
    # -------------------------------------------------------------------------
    lrbase_fpg      <- log(4.2);          label("Baseline fasting plasma glucose FPG0 before transplantation, at age 38 years (mmol/L)")  # Xiang 2025 Eq. 7 and Table 2: FPG0 = 4.2 mmol/L (RSE 2.3%; bootstrap median 4.2, 95% CI 4.1-4.3)
    lslope          <- fixed(log(0.264)); label("Slope of the linear tacrolimus effect on fasting plasma glucose (mmol/L per ng/mL)")     # Xiang 2025 Table 2: beta = 0.264 FIX; set from the published 6-month cumulative PTDM incidence, not estimated
    e_age_rbase_fpg <- 0.232;             label("Power exponent of (AGE / 38 years) on FPG0 (unitless)")                                  # Xiang 2025 Eq. 7 and Table 2: Age on FPG0 = 0.232 (RSE 23.9%; bootstrap 95% CI 0.131-0.345)

    # -------------------------------------------------------------------------
    # Inter-individual variability. The PK IIVs are carried over fixed from the
    # population PK step so that exposure variability propagates into the
    # simulated FPG; the FPG0 IIV is the quantity estimated here. Table 2
    # reports each IIV as a percentage, read as a coefficient of variation of
    # the log-normal parameter distribution and converted with
    # omega^2 = log(CV^2 + 1); see the vignette's Assumptions and deviations.
    #
    # Xiang 2025 carries NO inter-individual variability on the drug-effect
    # parameter beta, because beta is fixed. The Discussion names this as a
    # limitation: "the parameters associated with changes (beta for FPG and
    # Imax/IC50 for eGFR) are fixed, and there is no way to account for
    # individual variability in PD".
    # -------------------------------------------------------------------------
    etalcl ~ fixed(0.1009994)  # Xiang 2025 Table 2: IIV CL/F = 32.6% -> log(1 + 0.326^2); carried over from the population PK step
    etalvc ~ fixed(0.5794162)  # Xiang 2025 Table 2: IIV V/F  = 88.6% -> log(1 + 0.886^2); carried over from the population PK step

    etalrbase_fpg ~ 0.0284903  # Xiang 2025 Table 2: IIV FPG0 = 17.0% (RSE 19.9%; bootstrap 95% CI 13.6-19.5) -> log(1 + 0.170^2)

    # -------------------------------------------------------------------------
    # Residual error. Results: "The proportional error model was found to best
    # describe the residuals of the PK/PD models."
    # -------------------------------------------------------------------------
    propSd <- 0.231; label("Proportional residual error for fasting plasma glucose (fraction)")  # Xiang 2025 Table 2: Residual error Proportional = 23.1% (RSE 1.2%; bootstrap 95% CI 21.7-24.6)
  })

  model({
    # -----------------------------------------------------------------------
    # Tacrolimus PK layer, identical to Xiang_2025_tacrolimus (Eq. 5 and 6).
    # -----------------------------------------------------------------------
    ka <- exp(lka)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl) *
      (POD / 34)^e_pod_cl *
      exp(e_conmed_wuzhi_cl * CONMED_WUZHI) *
      exp(e_cyp3a5_expr_cl * (1 - CYP3A5_EXPR))

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in mg over V/F in L gives mg/L; 1 mg/L = 1000 ng/mL.
    Cc <- 1000 * central / vc

    # -----------------------------------------------------------------------
    # Individual baseline fasting plasma glucose, Xiang 2025 Eq. 7.
    # -----------------------------------------------------------------------
    rbase_fpg <- exp(lrbase_fpg + etalrbase_fpg) * (AGE / 38)^e_age_rbase_fpg

    slope <- exp(lslope)

    # -----------------------------------------------------------------------
    # Direct linear drug effect, Xiang 2025 Eq. 2: FPG = FPG0 + beta * Cp,t.
    # Cp,t is defined in the paper as "the trough concentration of tacrolimus
    # at time t"; every fasting-plasma-glucose observation in the source data
    # was paired with a trough sample. The driver here is the model-predicted
    # concentration at the observation time, so this reproduces the fitted
    # relationship exactly at trough times and interpolates it within the
    # dosing interval. Evaluate simulated FPG at trough times to stay on the
    # relationship the authors fitted; see the vignette.
    # -----------------------------------------------------------------------
    fpg <- rbase_fpg + slope * Cc

    fpg ~ prop(propSd)
  })
}
