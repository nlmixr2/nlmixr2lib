Frobel_2013_ciclosporin <- function() {
  description <- "Parametric time-to-event (TTE) model for the first acute rejection (AR) event after paediatric kidney transplantation in patients receiving oral ciclosporin A (Neoral microemulsion). The baseline hazard is a five-interval step-function exponential with break-points at 5, 8, 25, and 100 days after transplantation. The final model carries no covariates: 15 candidate covariates (including ciclosporin AUC, baseline AUC, demographics, donor characteristics, HLA mismatches, dialysis time, basiliximab induction) were screened by univariate testing, stepwise covariate modelling, cross-validated SCM, and bootstrap-SCM, and none reached statistical significance or clinical relevance. The model output `sur` is the probability of remaining acute-rejection-free at time t; `hazard` and `cumhaz` are exposed as derived outputs."
  reference <- paste(
    "Frobel A-K, Karlsson MO, Backman JT, Hoppu K, Qvist E, Seikku P,",
    "Jalanko H, Holmberg C, Keizer RJ, Fanta S, Jonsson S.",
    "A time-to-event model for acute rejections in paediatric renal",
    "transplant recipients treated with ciclosporin A.",
    "Br J Clin Pharmacol. 2013;76(4):603-615.",
    "doi:10.1111/bcp.12121.",
    sep = " "
  )
  vignette <- "Frobel_2013_ciclosporin"
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; the time-to-event model has no PK/PD parameters)",
    concentration = "probability (the model output `sur` is a survival probability, not a drug concentration)"
  )

  covariateData <- list()

  # All 15 covariates that the authors screened but did not retain in the
  # final model. Kept here for provenance only. checkModelConventions()
  # treats covariatesDataExcluded as documentation: entries must not be
  # referenced in model().
  covariatesDataExcluded <- list(
    AUC = list(
      description        = "Daily ciclosporin A exposure (h*mg/L). Computed from the individual empirical-Bayes estimates of bioavailability and total clearance from the Fanta et al. popPK model and the daily dose (Methods Equation 3).",
      units              = "h*mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in univariate testing (delta-OFV = 1.15, P = 0.2827, Table 4); not retained. Time-varying via dose changes and EBE-driven CL fluctuations.",
      source_name        = "AUC"
    ),
    AUC_BASELINE = list(
      description        = "Baseline (at-transplantation) daily ciclosporin A exposure.",
      units              = "h*mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in univariate testing (delta-OFV = 2.01, P = 0.1563, Table 4); preferred over time-varying AUC in further screening but still not retained.",
      source_name        = "Baseline AUC"
    ),
    AGE = list(
      description        = "Patient age (years).",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Screened delta-OFV = 1.87, P = 0.1716 (Table 4); not retained.",
      source_name        = "Age"
    ),
    AGE_BASELINE = list(
      description        = "Patient age at transplantation (years).",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened delta-OFV = 2.21, P = 0.1369 (Table 4); not retained.",
      source_name        = "Baseline age"
    ),
    WT = list(
      description        = "Patient body weight (kg).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Screened delta-OFV = 2.70, P = 0.1006 (Table 4); not retained.",
      source_name        = "Bodyweight"
    ),
    WT_BASELINE = list(
      description        = "Patient body weight at transplantation (kg).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened delta-OFV = 3.85, P = 0.0497 (Table 4); selected in the forward SCM but eliminated in the backward step at P < 0.01. The bootstrap-SCM inclusion frequency was 43% (the fourth-most-included covariate).",
      source_name        = "Baseline bodyweight"
    ),
    SEXF = list(
      description        = "Patient sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male, the most frequent category in the cohort: 57/87)",
      notes              = "Screened delta-OFV = 2.99, P = 0.0836 (Table 4); selected in the forward SCM but eliminated in the backward step at P < 0.01. The bootstrap-SCM inclusion frequency was 46% (tied for second-most-included covariate). Encoded against the canonical SEXF = 1 (female) convention.",
      source_name        = "Sex"
    ),
    DIALYSIS_TIME = list(
      description        = "Time on dialysis before transplantation (years).",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened delta-OFV = 4.80, P = 0.0284 (Table 4); selected first in the forward SCM but eliminated in the backward step at P < 0.01. The bootstrap-SCM inclusion frequency was 58% (the most-included covariate).",
      source_name        = "Dialysis time"
    ),
    BASILIXIMAB = list(
      description        = "Basiliximab induction at transplantation indicator (1 = yes, 0 = no).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (yes, the most frequent category in the cohort: 53/87 after 1999)",
      notes              = "Screened delta-OFV = 3.68, P = 0.0550 (Table 4); not selected in the SCM forward step. The bootstrap-SCM inclusion frequency was 46% (tied for second-most-included covariate).",
      source_name        = "Basiliximab induction"
    ),
    DIAGNOSIS_CNF = list(
      description        = "Diagnosis dichotomous indicator: 1 = congenital nephrotic syndrome (CNF), 0 = any other diagnosis (posterior urethral valve, polycystic kidney disease, nephronophthisis, other).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any non-CNF diagnosis)",
      notes              = "Screened delta-OFV = 3.76, P = 0.0524 (Table 4); dichotomous form preferred over the polychotomous form (delta-OFV = 4.97, P = 0.2900). Not retained.",
      source_name        = "Diagnosis (dichotomous)"
    ),
    CMV_DONOR = list(
      description        = "Donor CMV positive indicator (1 = yes, 0 = no).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (yes, the most frequent category: 61/87)",
      notes              = "Screened delta-OFV = 0.55, P = 0.4592 (Table 4); not retained.",
      source_name        = "CMV donor"
    ),
    CMV_RECIPIENT = list(
      description        = "Recipient CMV positive indicator (1 = yes, 0 = no).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no, the most frequent category: 57/87)",
      notes              = "Screened delta-OFV = 0.08, P = 0.7734 (Table 4); not retained.",
      source_name        = "CMV recipient"
    ),
    DONOR_AGE = list(
      description        = "Donor age (years).",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened delta-OFV = 0.17, P = 0.6765 (Table 4); not retained.",
      source_name        = "Donor age"
    ),
    DONOR_SEXF = list(
      description        = "Donor sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male, the most frequent category: 53/87)",
      notes              = "Screened delta-OFV = 0.04, P = 0.8345 (Table 4); not retained.",
      source_name        = "Donor sex"
    ),
    DONOR_TYPE = list(
      description        = "Donor type indicator: 1 = living related donor, 0 = cadaver donor.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (cadaver donor, the most frequent category: 71/87)",
      notes              = "Screened delta-OFV = 0.98, P = 0.3230 (Table 4); not retained.",
      source_name        = "Donor type"
    ),
    COLD_ISCHAEMIA_TIME = list(
      description        = "Renal transplant cold ischaemia time (hours).",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened delta-OFV = 0.92, P = 0.3381 (Table 4); not retained.",
      source_name        = "Cold ischaemia time"
    ),
    HLA_AB_MM_ANY = list(
      description        = "HLA-AB mismatch dichotomous indicator: 1 = any HLA-AB mismatch, 0 = zero HLA-AB mismatches.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (any mismatch, the most frequent category)",
      notes              = "Dichotomous form (delta-OFV = 0.59, P = 0.4420) and polychotomous form (delta-OFV = 3.38, P = 0.3363, 3 df) screened (Table 4); not retained.",
      source_name        = "HLA-AB mismatches"
    ),
    HLA_DR_MM_ANY = list(
      description        = "HLA-DR mismatch dichotomous indicator: 1 = any HLA-DR mismatch, 0 = zero HLA-DR mismatches.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (zero mismatches, used as reference per the paper's univariate framing)",
      notes              = "Dichotomous form (delta-OFV = 1.46, P = 0.2262) and polychotomous form (delta-OFV = 2.34, P = 0.3097, 2 df) screened (Table 4); the dichotomous form was preferred for the SCM. Not retained.",
      source_name        = "HLA-DR mismatches"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 87L,
    n_studies       = 1L,
    age_range       = "0.67-19.78 years (median 7.13 across all observations; baseline age range 0.67-18.17 years, median 4.51 years at transplantation)",
    weight_range    = "8.5-91.5 kg (median 21.6 across all observations; baseline weight range 8.8-68.7 kg, median 19.1 kg at transplantation)",
    sex_female_pct  = 34.5,
    race_ethnicity  = "Caucasian (homogeneous; single-centre Helsinki cohort)",
    disease_state   = "Paediatric kidney transplant recipients. Underlying diagnoses: congenital nephrotic syndrome of the Finnish type (CNF, NPHS1) 33%, posterior urethral valve 10%, nephronophthisis 8%, polycystic kidney disease 7%, other 40%. Median time on dialysis pre-transplant 0.96 years (range 0.01-3.86).",
    dose_range      = "Oral ciclosporin A (Neoral microemulsion); target trough 300 ug/L immediately post-transplant, reduced to 100 ug/L after 6 months. Median daily dose 185 mg (range 15-1100); median weight-normalised daily dose 9.57 mg/kg/day (range 0.76-36.92). Median daily AUC 14.40 h*mg/L (range 0.70-42.63).",
    regions         = "Finland (Children's Hospital, Helsinki; single-centre)",
    observation_window = "Median 3 years post-transplant (range 31 days to 14 years); the longest follow-up was 5111 days (~14 years) in one patient. Observation ended either at age ~18-19 years (transition to adult care) or at the end of data extraction in April 2006.",
    co_medication   = "Triple immunosuppression (methylprednisolone + azathioprine + ciclosporin A) before September 1999; from September 1999 onwards basiliximab induction (10 mg if WT < 30 kg, 20 mg if WT > 30 kg, intraoperatively and on day 4) was added, and azathioprine could be replaced by mycophenolate mofetil. Methylprednisolone was tapered from 1 mg/kg/day to 0.25 mg/kg/day at 3 weeks and 0.37 mg/kg every other day at 3 months.",
    notes           = "Retrospective single-centre analysis (1995-2006) of consecutive paediatric kidney transplant recipients at the only Finnish paediatric renal transplantation centre. 89 patients were eligible; 2 were excluded for incomplete records, giving N = 87. 54/87 (62%) experienced a first acute rejection (AR) during the observation window; the remaining 33/87 (38%) were treated as right-censored. AR was diagnosed by routine fine-needle aspiration biopsy (taken on day 5 post-transplant and at least twice weekly until discharge, plus whenever AR was clinically suspected). The PK was previously described by Fanta et al. (2007 / 2008) and supplied this analysis with individual empirical-Bayes estimates of CL and F used to compute AUC; the present paper develops only the TTE model on top of those PK posteriors."
  )

  ini({
    # Hazards on the natural log scale. Estimates are from Table 3
    # (NONMEM final estimates, $ESTIMATION METHOD=0 LIKE, parametric
    # survival in NONMEM 7.2.0). The hazard units are 1/day; the log
    # form is taken inside ini() so the natural typical-value
    # parameterisation in model() is exp(lh<n>).
    lh1 <- log(0.00465); label("Log baseline hazard for t <= 5 days post-transplant (1/day)")  # Table 3, row 1
    lh2 <- log(0.05780); label("Log baseline hazard for 5 < t <= 8 days post-transplant (1/day)")  # Table 3, row 2
    lh3 <- log(0.01870); label("Log baseline hazard for 8 < t <= 25 days post-transplant (1/day)")  # Table 3, row 3
    lh4 <- log(0.00470); label("Log baseline hazard for 25 < t <= 100 days post-transplant (1/day)")  # Table 3, row 4
    lh5 <- log(0.00013); label("Log baseline hazard for t > 100 days post-transplant (1/day)")  # Table 3, row 5

    # No IIV. The paper explicitly states that "as only one observation
    # was available per individual, random effects on the baseline
    # hazard could not be estimated, i.e. the same baseline hazard was
    # assumed for all subjects" (Methods, "Development of the base
    # model" section). No eta* parameters are added.

    # No residual error. The NONMEM run uses $ESTIMATION METHOD=0 LIKE
    # so the likelihood is the survival/event density itself, not an
    # observation-error model. This nlmixr2 translation is intended for
    # forward simulation; `hazard`, `cumhaz`, and `sur` are exposed as
    # derived outputs.
  })
  model({
    # Back-transformed piecewise hazards (1/day).
    h1 <- exp(lh1)
    h2 <- exp(lh2)
    h3 <- exp(lh3)
    h4 <- exp(lh4)
    h5 <- exp(lh5)

    # Piecewise-exponential step-function baseline hazard, with
    # break-points identified by tentatively removing steps of an
    # initial 15-step curve and keeping reductions that did not
    # significantly increase OFV (Methods, "Development of the base
    # model"; Figure 1). The final time cut-offs in days are 5, 8,
    # 25, and 100.
    hazard <- ifelse(t <= 5, h1,
              ifelse(t <= 8, h2,
              ifelse(t <= 25, h3,
              ifelse(t <= 100, h4, h5))))

    # Cumulative hazard and survival probability.
    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    sur <- exp(-cumhaz)
  })
}
