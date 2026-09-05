Damnjanovic_2023_lamotrigine <- function() {
  description <- "One-compartment population PK model with first-order absorption and elimination for lamotrigine (LTG) in Serbian children aged 2-18 years on dual antiepileptic therapy (Damnjanovic 2023 Table 2B). Body weight enters apparent volume as a power term with an ESTIMATED exponent of 2.83 referenced to the 37.1 kg cohort mean; apparent clearance carries the patient's own total daily lamotrigine dose as an exponential-linear term and a binary valproate-comedication term that reduces CL/F by 46% (exp(-0.61)), reproducing valproate's known inhibition of lamotrigine glucuronidation. Ka was FIXED at 1.57 1/h from the literature because the trough-only sampling carried no absorption information. Residual error is proportional. Fit in Monolix 2021R2 to a single steady-state trough per patient."
  reference   <- "Damnjanovic I, Tsyplakova N, Stefanovic N, Tosic T, Catic-Djordjevic A, Karalis V. Joint use of population pharmacokinetics and machine learning for optimizing antiepileptic treatment in pediatric population. Ther Adv Drug Saf. 2023;14:20420986231181337. doi:10.1177/20420986231181337. PMCID PMC10288421. Parameters from Table 2(b); cohort demographics from Table 1."
  vignette    <- "Damnjanovic_2023_pediatric_antiepileptics"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot   = list(analyte = "lamotrigine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lamotrigine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference weight 37.1 kg = the cohort mean body weight (Damnjanovic 2023 Table 1A, 'Body weight (kg) / Mean' row). Enters V/F only, as (WT/37.1)^2.83; the Results narrative describes it as 'BW on apparent V (using a weight-centered individual model)'. Unlike the levetiracetam and valproic-acid models, this exponent was ESTIMATED rather than fixed at the allometric 1 (SE 0.52, RSE 18.4%). An exponent of 2.83 is far above any physiologically defensible allometric value for a distribution volume; it is transcribed as published and the consequence for the simulated weight range is quantified in the vignette Errata. Unlike the other two models, lamotrigine CL/F carries NO weight term at all in Table 2(b).",
      source_name        = "BW"
    ),
    DOSE_LTG_MGD = list(
      description        = "Patient's own total daily lamotrigine dose",
      units              = "mg/d",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters CL/F as the exponential-linear term exp(0.0056 * DOSE_LTG_MGD), i.e. a 0.56% rise in apparent clearance per mg/day of lamotrigine. Per-record covariate carrying the current daily dose LEVEL, summed across the day for a twice-daily regimen; distinct from the rxode2 event-table amt column, which carries each individual administration. The paper reports no cohort mean daily dose anywhere in the article or supplement, so the covariate cannot be mean-centred and is applied untransformed - consistent with Methods, which offered continuous covariates 'either untransformed or centered on their mean value', and with the fact that only the untransformed reading reproduces the pediatric lamotrigine clearances of 0.7-1.5 L/h tabulated in the paper's own Supplemental Table S3. See vignette Errata.",
      source_name        = "DailyDose"
    ),
    CONMED_VPA = list(
      description        = "Concomitant valproate (valproic acid) therapy",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = lamotrigine given without valproate (in this cohort, the LTG/LEV regimen)",
      notes              = "Table 2(b) calls this covariate 'Regimen' and the table footnote defines it as 'the relationship between Cl and therapeutic regimen (whether existence of valproic acid)'. Coded 1 when the child's dual therapy contains valproate (VA/LTG, n = 42) and 0 when it does not (LTG/LEV, n = 9). The sign is anchored by the Discussion: 'It was found that the coadministration of VA, which is an inhibitor, decreases LTG clearance', so the negative coefficient must attach to the valproate-present stratum and the valproate-free stratum is the reference. exp(-0.61) = 0.543, a 45.7% reduction in CL/F.",
      source_name        = "Regimen"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 51,
    n_studies       = 1,
    n_observations  = 51,
    age_range       = "2-18 years (Damnjanovic 2023 Methods, inclusion criteria).",
    age_median      = "Median 11 years; mean 10.9 years; interquartile range 7 years (Table 1A).",
    weight_range    = "Not reported as a range; median 35 kg, mean 37.1 kg, interquartile range 21 kg (Table 1A).",
    weight_median   = "35 kg (Table 1A).",
    sex_female_pct  = 56.3,
    race_ethnicity  = c(NotReported = 100),
    disease_state   = "Children with diagnosed epilepsy (ICD-10 G40) on dual antiepileptic therapy. Poor renal or hepatic function and other serious disease states were exclusion criteria (Methods).",
    dose_range      = "Not reported. Doses were the patients' own prescribed maintenance regimens. The daily dose is a model covariate but its cohort distribution is not tabulated anywhere in the article or its supplement.",
    regions         = "Serbia (single centre: Clinic of Pediatric Internal Medicine, University Clinical Center of Nis).",
    co_medication   = "All patients were on dual antiepileptic therapy. Lamotrigine was given with valproic acid (VA/LTG, n = 42) or with levetiracetam (LTG/LEV, n = 9). Valproate coadministration is a retained covariate on CL/F.",
    notes           = "The whole study enrolled 71 children across three dual-therapy regimens (VA/LTG n = 42, VA/LEV n = 20, LTG/LEV n = 9; Table 1B). The LTG sub-cohort is the 51 patients whose regimen contained lamotrigine. That denominator is confirmed by the paper's own reporting: 86.27% of LTG concentrations were within the 3-15 mg/L reference range = 44/51 (Results paragraph 1). Exactly one steady-state trough concentration was drawn per patient, immediately before the next morning dose, so all parameters were identified from a single sample per subject via a stepwise fix-and-release estimation strategy (Methods, 'Population pharmacokinetics')."
  )

  ini({
    # Structural parameters - final-model estimates from Damnjanovic 2023 Table 2(b).

    # Ka FIXED at the literature value: "Given that the gathered LTG steady-state
    # trough plasma concentration data did not give information regarding the extent
    # and rate of absorption processes, the absorption rate was set at 1.57 h-1 based
    # on literature values." (Results, LTG paragraph). Table 2(b) prints no SE / RSE.
    lka <- fixed(log(1.57));  label("Absorption rate constant (Ka, 1/h)")                  # Table 2(b): Ka = 1.57 1/h, no SE/RSE reported
    lvc <- log(5.15);         label("Apparent central volume of distribution (V/F, L)")    # Table 2(b): V = 5.15 L (SE 1.18, RSE 22.9%)
    lcl <- log(0.15);         label("Apparent oral clearance (CL/F, L/h)")                 # Table 2(b): Cl = 0.15 L/h (SE 0.02, RSE 13.3%)

    # Weight on V/F. This exponent was ESTIMATED (SE and RSE are both printed),
    # unlike the fixed 1 / 0.75 pair used for levetiracetam and valproic acid.
    e_wt_vc <- 2.83;          label("Power exponent: WT on V/F, referenced to 37.1 kg (unitless)")   # Table 2(b): beta_V_logBW = 2.83 (SE 0.52, RSE 18.4%, p < 0.001)

    # Lamotrigine's own total daily dose on CL/F, exponential-linear and
    # untransformed (see the covariateData note and the vignette Errata for why
    # the covariate is not mean-centred: the paper never reports a mean daily dose).
    e_dose_ltg_cl <- 0.0056;  label("Effect of total daily lamotrigine dose on CL/F (per mg/day, exponential-linear)")  # Table 2(b): beta_Cl_DailyDose = 0.0056 (SE 0.0007, RSE 13.0%, p < 0.001)

    # Valproate comedication on CL/F. Table 2(b) prints the SE for this row as
    # "-0.13", which is impossible - a standard error cannot be negative. The
    # magnitude is confirmed by the printed RSE: 0.61 * 0.222 = 0.135, so the
    # leading minus is a typesetting carry-over from the estimate above it. Only
    # the point estimate is used by the model, so this does not change any value.
    e_conmed_vpa_cl <- -0.61; label("Effect of concomitant valproate on CL/F (log-scale shift)")     # Table 2(b): beta_Cl_Regimen = -0.61 (SE printed as -0.13, RSE 22.2%, p < 0.001)

    # IIV. Monolix reports omega as the STANDARD DEVIATION of the log-scale random
    # effect; nlmixr2's ini() takes VARIANCES, so each omega is squared:
    #   omega_V  = 0.32 -> var = 0.1024
    #   omega_Cl = 0.28 -> var = 0.0784
    # Table 2(b) reports no V-Cl correlation for lamotrigine (unlike Table 2(a) for
    # levetiracetam), so the two etas are independent here.
    etalvc ~ 0.1024                                                                       # Table 2(b): omega_V = 0.32 (SE 0.07, RSE 23.4%)
    etalcl ~ 0.0784                                                                       # Table 2(b): omega_Cl = 0.28 (SE 0.07, RSE 22.7%)

    # Residual error. "The residual variability was estimated using a proportional
    # error model." (Results, LTG paragraph). Monolix's proportional error model is
    # y = f + b * f * eps, so b is the proportional SD as a fraction.
    propSd <- 0.15;           label("Proportional residual error (fraction)")              # Table 2(b): b = 0.15 (SE 0.03, RSE 19.8%)
  })

  model({
    # No IIV on Ka: Table 2(b) reports no omega_Ka, consistent with Ka being fixed.
    ka <- exp(lka)

    # Weight on V/F only. Table 2(b) carries NO weight term on CL/F.
    vc <- exp(lvc + etalvc) * (WT / 37.1)^e_wt_vc

    # CL/F covariate model on the log-parameter scale, as Monolix parameterises a
    # lognormal parameter: log(Cl_i) = log(Cl_pop) + beta_dose * DailyDose
    #                                 + beta_regimen * VPA + eta_i
    cl <- exp(lcl + etalcl) *
      exp(e_dose_ltg_cl * DOSE_LTG_MGD) *
      exp(e_conmed_vpa_cl * CONMED_VPA)

    kel <- cl / vc

    # One compartment, first-order oral absorption and first-order elimination.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in mg, V/F in L -> mg/L, the unit of the paper's reference range
    # (3-15 mg/L, Discussion paragraph 5; the article prints "3-15 mg/ml", a
    # typographic slip for mg/L - see the vignette Errata).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
