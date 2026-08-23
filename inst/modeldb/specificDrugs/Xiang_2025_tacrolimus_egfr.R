Xiang_2025_tacrolimus_egfr <- function() {
  description <- "Sequential exposure-response (Imax) model relating oral tacrolimus trough concentration to the estimated glomerular filtration rate in adult renal transplant recipients, as a quantitative model of calcineurin-inhibitor nephrotoxicity. The eGFR falls from its post-transplant stable baseline by a saturable Imax function of the predicted tacrolimus concentration; the baseline falls with age and rises with hemoglobin. Both Imax and IC50 were fixed at clinically anchored values because they could not be estimated interpretably. The tacrolimus PK is the one-compartment model of the same paper, carried here as fixed parameters because the PK/PD model was fitted sequentially on the population PK step; see Xiang_2025_tacrolimus. The companion glucose model in the same paper is Xiang_2025_tacrolimus_fpg."
  reference <- paste(
    "Xiang Q, Yang Y, Li G, Chen S, Yang Y, Liu L, Yu X.",
    "Population Pharmacokinetic/Pharmacodynamic Modeling of Tacrolimus in Renal",
    "Transplant Recipients: Impact of CYP3A5 Genotype and Wuzhi Capsule",
    "Co-Medication. Drug Des Devel Ther. 2025;19:8375-8389.",
    "doi:10.2147/DDDT.S542786.",
    "Parameter values are the final-model estimates in Table 2 ('C0-eGFR model'",
    "block, plus the 'PPK model' block for the fixed PK); the structural",
    "equations are Eq. 1, 3, 5, 6 and 8.",
    sep = " "
  )
  vignette <- "Xiang_2025_tacrolimus"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  # Issue #482: the ODE states are the tacrolimus PK carried over from the
  # population PK step. The estimated glomerular filtration rate is an
  # algebraic direct effect of the predicted concentration and is not a state.
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
      notes              = "TIME-VARYING within subject. Enters the fixed PK layer as the power scaling (POD / 34)^0.109 on CL/F; 34 days is the median post-operative day of the index cohort (Xiang 2025 Table 1). The power form is undefined at POD = 0, so supply POD >= 1. It has no direct effect on eGFR; it acts only through the predicted tacrolimus concentration. Note the eGFR model's own time origin is different from the PK model's: eGFR0 is anchored at the day renal function became stable, taken as post-operative day 7 (Xiang 2025 Methods, 'Simulations'). POD was significant on eGFR0 in forward inclusion but failed backward elimination; see covariatesDataExcluded.",
      source_name        = "POD"
    ),
    CONMED_WUZHI = list(
      description        = "Concomitant Wuzhi capsule (Schisandra sphenanthera extract) indicator (1 = receiving Wuzhi capsule, 0 = not)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant Wuzhi capsule)",
      notes              = "Source column WZ, defined below Xiang 2025 Eq. 6 as 'WZ = 1 if WZ is present; otherwise = 0'. Enters the fixed PK layer as exp(-0.211 * CONMED_WUZHI) on CL/F. It has no direct effect on eGFR; it acts only through the predicted tacrolimus concentration, which is the paper's mechanism for the headline result that CYP3A5*3/*3 patients taking the Wuzhi capsule can lose more than 20% of baseline eGFR even at 2 mg q12h.",
      source_name        = "WZ"
    ),
    CYP3A5_EXPR = list(
      description        = "CYP3A5 expresser status (1 = carries at least one CYP3A5*1 allele, genotype *1/*1 or *1/*3; 0 = CYP3A5*3/*3 non-expresser)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5*3/*3 non-expresser)",
      notes              = "VALUE INVERSION relative to the source: Xiang 2025 codes 'Genotype = 1 if the genotype is CYP3A5*3/*3' below Eq. 6, so CYP3A5_EXPR = 1 - Genotype and the model applies the published coefficient to (1 - CYP3A5_EXPR). Enters the fixed PK layer only; it has no direct effect on eGFR.",
      source_name        = "Genotype"
    ),
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed. Enters as the power scaling (AGE / 38)^-0.403 on the stable post-transplant baseline eGFR (Xiang 2025 Eq. 8); 38 years is the median age of the index cohort. Retained with dOFV = 39.013, df = 1, p < 0.001. Age does not modify the tacrolimus effect on eGFR, only its baseline.",
      source_name        = "Age"
    ),
    HGB = list(
      description        = "Blood hemoglobin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed in this model. Enters as the power scaling (HGB / 111)^0.254 on the stable post-transplant baseline eGFR (Xiang 2025 Eq. 8); 111 g/L is the median hemoglobin of the index cohort (range 50-178 g/L, Table 1). Units are SI g/L, not g/dL. Retained with dOFV = 17.784, df = 1, p < 0.001. The paper reads hemoglobin as an indirect marker of renal function, because erythropoietin secretion falls as the graft's function declines (Xiang 2025 Discussion).",
      source_name        = "HB"
    )
  )

  # Screened in the paper's covariate analysis but NOT retained in the final
  # model, so documented rather than used (Xiang 2025 Discussion and
  # Supplementary Table 4): both reached the forward-inclusion threshold on
  # eGFR0 but neither met the backward-elimination criterion (p < 0.001), and
  # no point estimate is published for either.
  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on baseline eGFR and significant in forward inclusion only (dOFV = -8.451, p < 0.01; Xiang 2025 Discussion and Supplementary Table 4). It did not meet the backward-elimination criterion (p < 0.001) and is absent from the final model, with no published point estimate. Index-cohort weight 31-99 kg (median 60)."
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
    disease_state  = "Adult renal transplant recipients on a tacrolimus-based immunosuppressive regimen; observed eGFR 3.1-171.3 mL/min/1.73 m^2 (mean 55.18, median 54.5) and hemoglobin 50-178 g/L (median 111)",
    dose_range     = "Oral tacrolimus started within 24 h of transplantation at 0.05-0.15 mg/kg/day divided q12h, then dose-adjusted to the trough TDM target; observed doses 0.5-5.5 mg per administration (median 2.5)",
    regions        = "Single centre, Chongqing, People's Republic of China",
    notes          = "Retrospective cohort of 126 renal transplant recipients randomly split 4:1 into an index group (n = 100) used to build the model and an external validation group (n = 26). 2334 eGFR observations. Baseline demographics from Xiang 2025 Table 1 (index dataset column). Estimated in Phoenix NLME 8.3.5, sequentially on the population PK step. The paper's renal-injury endpoint is a relative fall in eGFR from baseline of more than 20%, computed as ReGFR = (eGFR0 - eGFR) / eGFR0 * 100% (Eq. 4)."
  )

  ini({
    # -------------------------------------------------------------------------
    # Tacrolimus PK, fixed from the population PK step of the same paper
    # (Xiang 2025 Eq. 5 and 6, Table 2 'PPK model' block). Methods: "The PK/PD
    # models were sequentially built following the population PK model", and
    # Table 2 does not re-report PK parameters in the C0-eGFR block, so none of
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
    # Estimated glomerular filtration rate. Xiang 2025 Eq. 3:
    #   eGFR = eGFR0 - Imax * Cp,t / (IC50 + Cp,t)
    # and Eq. 8:
    #   eGFR0 = 52.3 * (Age/38)^-0.403 * (HB/111)^0.254 * e^(etaeGFR0)
    #
    # Imax and IC50 were both FIXED rather than estimated. Discussion:
    # "Different effect models were explored to link C0 and eGFR, but it was
    # discovered that the parameter estimates of IC50 and Imax cannot be
    # reasonably interpreted in clinical practice." IC50 was then anchored at
    # 10 ng/mL because "the incidence of nephrotoxicity will greatly increase
    # when the C0 of tacrolimus is higher than 10 ng/mL" (the authors'
    # reference 43), and Imax at 30 mL/min/1.73 m^2 because an eGFR below that
    # value is stage 4 chronic kidney disease (reference 44). A +/-30%
    # sensitivity analysis changed effect values by < 15% (Supplementary
    # Table 5).
    #
    # Imax is an ABSOLUTE maximum reduction in mL/min/1.73 m^2, not a fraction
    # of baseline, so Eq. 3 subtracts rather than multiplies.
    # -------------------------------------------------------------------------
    lrbase_egfr      <- log(52.3);      label("Stable post-transplant baseline eGFR0, at age 38 years and hemoglobin 111 g/L (mL/min/1.73 m^2)")  # Xiang 2025 Eq. 8 and Table 2: eGFR0 = 52.3 mL/min/1.73 m^2 (RSE 2.8%; bootstrap median 52.0, 95% CI 42.6-61.2)
    limax            <- fixed(log(30)); label("Maximum tacrolimus-attributable reduction in eGFR, Imax (mL/min/1.73 m^2)")                         # Xiang 2025 Table 2: Imax = 30 mL/min/1.73 m^2 FIX; anchored at the stage-4 chronic-kidney-disease threshold, not estimated
    lic50            <- fixed(log(10)); label("Tacrolimus concentration producing half of Imax, IC50 (ng/mL)")                                     # Xiang 2025 Table 2: IC50 = 10 ng/mL FIX; anchored at the published trough above which nephrotoxicity incidence rises steeply, not estimated
    e_age_rbase_egfr <- -0.403;         label("Power exponent of (AGE / 38 years) on eGFR0 (unitless)")                                            # Xiang 2025 Eq. 8 and Table 2: Age on eGFR0 = -0.403 (RSE 14.8%; bootstrap 95% CI -0.525 to -0.254)
    e_hgb_rbase_egfr <- 0.254;          label("Power exponent of (HGB / 111 g/L) on eGFR0 (unitless)")                                             # Xiang 2025 Eq. 8 and Table 2: HB on eGFR0 = 0.254 (RSE 4.9%; bootstrap 95% CI 0.132-0.417)

    # -------------------------------------------------------------------------
    # Inter-individual variability. The PK IIVs are carried over fixed from the
    # population PK step so that exposure variability propagates into the
    # simulated eGFR; the eGFR0 IIV is the quantity estimated here. Table 2
    # reports each IIV as a percentage, read as a coefficient of variation of
    # the log-normal parameter distribution and converted with
    # omega^2 = log(CV^2 + 1); see the vignette's Assumptions and deviations.
    #
    # Xiang 2025 carries NO inter-individual variability on the drug-effect
    # parameters Imax and IC50, because both are fixed. The Discussion names
    # this as a limitation: fixing them "made it difficult to estimate the
    # inter-individual variation of PD indicators".
    # -------------------------------------------------------------------------
    etalcl ~ fixed(0.1009994)  # Xiang 2025 Table 2: IIV CL/F = 32.6% -> log(1 + 0.326^2); carried over from the population PK step
    etalvc ~ fixed(0.5794162)  # Xiang 2025 Table 2: IIV V/F  = 88.6% -> log(1 + 0.886^2); carried over from the population PK step

    etalrbase_egfr ~ 0.0476857  # Xiang 2025 Table 2: IIV eGFR0 = 22.1% (RSE 13.8%; bootstrap 95% CI 18.1-25.6) -> log(1 + 0.221^2)

    # -------------------------------------------------------------------------
    # Residual error. Results: "The proportional error model was found to best
    # describe the residuals of the PK/PD models."
    # -------------------------------------------------------------------------
    propSd <- 0.172; label("Proportional residual error for eGFR (fraction)")  # Xiang 2025 Table 2: Residual error Proportional = 17.2% (RSE 0.6%; bootstrap 95% CI 15.6-18.7)
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
    # Individual stable post-transplant baseline eGFR, Xiang 2025 Eq. 8.
    # -----------------------------------------------------------------------
    rbase_egfr <- exp(lrbase_egfr + etalrbase_egfr) *
      (AGE / 38)^e_age_rbase_egfr *
      (HGB / 111)^e_hgb_rbase_egfr

    imax <- exp(limax)
    ic50 <- exp(lic50)

    # -----------------------------------------------------------------------
    # Direct Imax drug effect, Xiang 2025 Eq. 3:
    #   eGFR = eGFR0 - Imax * Cp,t / (IC50 + Cp,t)
    # Cp,t is defined in the paper as "the trough concentration of tacrolimus
    # at time t"; every eGFR observation in the source data was paired with a
    # trough sample. The driver here is the model-predicted concentration at
    # the observation time, so this reproduces the fitted relationship exactly
    # at trough times and interpolates it within the dosing interval. Evaluate
    # simulated eGFR at trough times to stay on the relationship the authors
    # fitted; see the vignette.
    # -----------------------------------------------------------------------
    egfr <- rbase_egfr - imax * Cc / (ic50 + Cc)

    egfr ~ prop(propSd)
  })
}
