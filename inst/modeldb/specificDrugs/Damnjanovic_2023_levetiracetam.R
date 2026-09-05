Damnjanovic_2023_levetiracetam <- function() {
  description <- "One-compartment population PK model with first-order absorption and elimination for levetiracetam (LEV) in Serbian children aged 2-18 years on dual antiepileptic therapy (Damnjanovic 2023 Table 2A). Body weight is the only retained covariate and enters both apparent volume and apparent clearance as fixed-exponent allometric terms (1 for V/F, 0.75 for CL/F) referenced to the 37.1 kg cohort mean. Ka was FIXED at 2.6 1/h from the literature because the trough-only sampling carried no absorption information. Inter-individual variability on V/F and CL/F is a correlated block (r = 0.86) and residual error is additive. Fit in Monolix 2021R2 to a single steady-state trough per patient."
  reference   <- "Damnjanovic I, Tsyplakova N, Stefanovic N, Tosic T, Catic-Djordjevic A, Karalis V. Joint use of population pharmacokinetics and machine learning for optimizing antiepileptic treatment in pediatric population. Ther Adv Drug Saf. 2023;14:20420986231181337. doi:10.1177/20420986231181337. PMCID PMC10288421. Parameters from Table 2(a); cohort demographics from Table 1."
  vignette    <- "Damnjanovic_2023_pediatric_antiepileptics"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot   = list(analyte = "levetiracetam", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "levetiracetam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference weight 37.1 kg = the cohort mean body weight (Damnjanovic 2023 Table 1A, 'Body weight (kg) / Mean' row; the Results narrative calls the same number the median). Enters as (WT/37.1)^1 on V/F and (WT/37.1)^0.75 on CL/F, both exponents FIXED (Methods: 'Using fixed exponents, allometric scaling was applied to the apparent volumes of distribution (V) and clearance (Cl) parameters (1 for V and 0.75 for Cl)'). The paper's covariate is named beta_V_logBW / beta_Cl_logBW, i.e. the coefficient of a log-transformed weight; a log(WT) term without centring would put V/F near 900 L, so the transformation must be log(WT/WTref) and the Methods statement that continuous covariates were 'centered on their mean value' fixes WTref at the cohort mean. Baseline weight; the paper reports one trough per patient so no time-varying weight schedule exists.",
      source_name        = "BW"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 29,
    n_studies       = 1,
    n_observations  = 29,
    age_range       = "2-18 years (Damnjanovic 2023 Methods, inclusion criteria).",
    age_median      = "Median 11 years; mean 10.9 years; interquartile range 7 years (Table 1A).",
    weight_range    = "Not reported as a range; median 35 kg, mean 37.1 kg, interquartile range 21 kg (Table 1A).",
    weight_median   = "35 kg (Table 1A).",
    sex_female_pct  = 56.3,
    race_ethnicity  = c(NotReported = 100),
    disease_state   = "Children with diagnosed epilepsy (ICD-10 G40) on dual antiepileptic therapy. Poor renal or hepatic function and other serious disease states were exclusion criteria (Methods).",
    dose_range      = "Not reported. Doses were the patients' own prescribed maintenance regimens; the paper states they were 'in accordance with leading recommendations' (Discussion).",
    regions         = "Serbia (single centre: Clinic of Pediatric Internal Medicine, University Clinical Center of Nis).",
    co_medication   = "All patients were on dual antiepileptic therapy. Levetiracetam was given with valproic acid (VA/LEV, n = 20) or with lamotrigine (LTG/LEV, n = 9). Coadministration of VA or LTG was tested on LEV clearance and was NOT statistically significant (Discussion paragraph 4), so no comedication covariate appears in this model.",
    notes           = "The whole study enrolled 71 children across three dual-therapy regimens (VA/LTG n = 42, VA/LEV n = 20, LTG/LEV n = 9; Table 1B). The LEV sub-cohort is the 29 patients whose regimen contained levetiracetam. That denominator is confirmed by the paper's own reporting: 68.97% of LEV concentrations were within the 12-46 mg/L reference range = 20/29, and 27.59% were below it = 8/29 (Results paragraph 1). Exactly one steady-state trough concentration was drawn per patient, immediately before the next morning dose, so all parameters were identified from a single sample per subject via a stepwise fix-and-release estimation strategy (Methods, 'Population pharmacokinetics')."
  )

  ini({
    # Structural parameters - final-model estimates from Damnjanovic 2023 Table 2(a).

    # Ka FIXED at the literature value: "Given that the collected LEV steady-state
    # trough plasma concentration data did not provide information about the extent
    # and rate of absorption processes, the absorption rate (Ka) was fixed at 2.6 h-1
    # according to the values in the existing literature." (Results, LEV paragraph 1).
    # Table 2(a) accordingly prints no SE / RSE for Ka.
    lka <- fixed(log(2.6));   label("Absorption rate constant (Ka, 1/h)")                  # Table 2(a): Ka = 2.6 1/h, no SE/RSE reported
    lvc <- log(25.01);        label("Apparent central volume of distribution (V/F, L)")    # Table 2(a): V = 25.01 L (SE 5.65, RSE 22.6%)
    lcl <- log(1.51);         label("Apparent oral clearance (CL/F, L/h)")                 # Table 2(a): Cl = 1.51 L/h (SE 0.27, RSE 18.1%)

    # Allometric exponents on body weight, both FIXED (Table 2(a) prints "-" for
    # SE and RSE on both rows; Methods states the exponents were fixed at 1 for V
    # and 0.75 for Cl). The p < 0.001 in Table 2(a) is the Wald test on the
    # covariate's inclusion, not on the exponent value.
    e_wt_vc <- fixed(1);      label("Allometric exponent: WT on V/F (unitless)")           # Table 2(a): beta_V_logBW = 1, p < 0.001
    e_wt_cl <- fixed(0.75);   label("Allometric exponent: WT on CL/F (unitless)")          # Table 2(a): beta_Cl_logBW = 0.75, p < 0.001

    # IIV. Monolix reports omega as the STANDARD DEVIATION of the random effect on
    # the log scale for lognormally-distributed parameters; nlmixr2's ini() takes
    # VARIANCES, so each omega is squared here:
    #   omega_V  = 0.84  -> var = 0.7056
    #   omega_Cl = 0.59  -> var = 0.3481
    # The off-diagonal is the covariance implied by the reported correlation:
    #   cov = corr_V_Cl * omega_V * omega_Cl = 0.86 * 0.84 * 0.59 = 0.426216
    # The resulting 2x2 block has determinance 0.7056*0.3481 - 0.426216^2 = 0.0640 > 0,
    # so it is positive definite and needs no nudge (failure pattern 1).
    etalvc + etalcl ~ c(0.7056,
                        0.426216, 0.3481)                                                 # Table 2(a): omega_V = 0.84 (RSE 19.1%), omega_Cl = 0.59 (RSE 17.2%), corr_V_Cl = 0.86 (RSE 21.4%)

    # Residual error. "The constant error model produced the best residual
    # variability performance of any residual error model studied." (Results, LEV
    # paragraph 1). Monolix's constant error model is y = f + a * eps, i.e. purely
    # additive with SD = a on the concentration scale.
    addSd <- 3.82;            label("Additive residual error SD (mg/L)")                   # Table 2(a): a = 3.82 (SE 0.96, RSE 25.0%)
  })

  model({
    # No IIV on Ka: Table 2(a) reports no omega_Ka, consistent with Ka being fixed.
    ka <- exp(lka)

    # Fixed-exponent allometry on the 37.1 kg cohort mean weight.
    vc <- exp(lvc + etalvc) * (WT / 37.1)^e_wt_vc
    cl <- exp(lcl + etalcl) * (WT / 37.1)^e_wt_cl

    kel <- cl / vc

    # One compartment, first-order oral absorption and first-order elimination.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in mg, V/F in L -> mg/L, the unit of the paper's reference range
    # (12-46 mg/L, Discussion paragraph 3).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
