Damnjanovic_2023_valproic_acid <- function() {
  description <- "One-compartment population PK model with first-order absorption and elimination for total plasma valproic acid (VA) in Serbian children aged 2-18 years on dual antiepileptic therapy (Damnjanovic 2023 Table 2C). Apparent volume carries a fixed allometric weight exponent of 1 referenced to the 37.1 kg cohort mean plus an estimated age effect centred on the 10.9-year cohort mean; apparent clearance rises with the patient's own total daily valproate dose through an exponential-linear term. Ka was FIXED at 1.68 1/h from the literature because the trough-only sampling carried no absorption information. Residual error is proportional. Fit in Monolix 2021R2 to a single steady-state trough per patient."
  reference   <- "Damnjanovic I, Tsyplakova N, Stefanovic N, Tosic T, Catic-Djordjevic A, Karalis V. Joint use of population pharmacokinetics and machine learning for optimizing antiepileptic treatment in pediatric population. Ther Adv Drug Saf. 2023;14:20420986231181337. doi:10.1177/20420986231181337. PMCID PMC10288421. Parameters from Table 2(c); cohort demographics from Table 1."
  vignette    <- "Damnjanovic_2023_pediatric_antiepileptics"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot   = list(analyte = "valproic acid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "valproic acid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference weight 37.1 kg = the cohort mean body weight (Damnjanovic 2023 Table 1A, 'Body weight (kg) / Mean' row). Enters V/F only, as (WT/37.1)^1, with the exponent FIXED (Table 2(c) prints '-' for SE and RSE; Methods: 'Using fixed exponents, allometric scaling was applied to the apparent volumes of distribution (V) and clearance (Cl) parameters (1 for V and 0.75 for Cl)'). The Results narrative for valproic acid says 'BW on apparent V (with an allometric exponent of 0.75)', which contradicts both Table 2(c) (beta_V_logBW = 1) and the Methods statement that 0.75 is the CLEARANCE exponent; the table value of 1 is used. Table 2(c) carries no weight term on CL/F at all.",
      source_name        = "BW"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters V/F as exp(0.07 * (AGE - 10.9)), centred on the 10.9-year cohort mean age (Table 1A, 'Age (years) / Mean' row), following the Methods statement that continuous covariates were 'centered on their mean value'. Centring is what makes V/F = 15.61 L the typical value for the typical child, which is the reading the authors themselves assert when they write that their apparent Cl and V estimates 'were very close to those reported in other PopPK studies' (Discussion paragraph 8) - the paediatric valproate volumes tabulated in Supplemental Table S4 span 2.88-22.12 L, bracketing 15.61 L, whereas the uncentred reading would give 33.5 L, above every one of them. The Results text calls the effect positive (+0.07, Table 2(c)); the PCA paragraph later refers to 'the negative value of beta_V_Age (i.e. equal to -0.07, Table 2)'. Table 2(c) and the Results text agree on the positive sign and are used. See vignette Errata.",
      source_name        = "Age"
    ),
    DOSE_VPA_MGD = list(
      description        = "Patient's own total daily valproic acid dose",
      units              = "mg/d",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters CL/F as the exponential-linear term exp(0.0012 * DOSE_VPA_MGD), i.e. a 0.12% rise in apparent clearance per mg/day of valproate. Table 2(c) names the coefficient beta_Cl_VA; the Results narrative identifies the covariate as 'daily dose on Cl' for valproic acid, and the Discussion confirms that coadministration of the other two antiepileptics was NOT significant on valproate clearance, so beta_Cl_VA is valproate's own daily dose rather than a comedication flag. Per-record covariate carrying the current daily dose LEVEL, summed across the day; distinct from the rxode2 event-table amt column, which carries each individual administration. The paper reports no cohort mean daily dose anywhere in the article or supplement, so the covariate cannot be mean-centred and is applied untransformed - the untransformed reading also puts CL/F at 0.2-0.7 L/h across plausible paediatric valproate doses, inside the 0.047-0.854 L/h range tabulated in the paper's own Supplemental Table S4, whereas the centred reading would fix it at 0.12 L/h. See vignette Errata.",
      source_name        = "VA"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 62,
    n_studies       = 1,
    n_observations  = 62,
    age_range       = "2-18 years (Damnjanovic 2023 Methods, inclusion criteria).",
    age_median      = "Median 11 years; mean 10.9 years; interquartile range 7 years (Table 1A).",
    weight_range    = "Not reported as a range; median 35 kg, mean 37.1 kg, interquartile range 21 kg (Table 1A).",
    weight_median   = "35 kg (Table 1A).",
    sex_female_pct  = 56.3,
    race_ethnicity  = c(NotReported = 100),
    disease_state   = "Children with diagnosed epilepsy (ICD-10 G40) on dual antiepileptic therapy. Poor renal or hepatic function and other serious disease states were exclusion criteria (Methods).",
    dose_range      = "Not reported. Doses were the patients' own prescribed maintenance regimens. The daily dose is a model covariate but its cohort distribution is not tabulated anywhere in the article or its supplement.",
    regions         = "Serbia (single centre: Clinic of Pediatric Internal Medicine, University Clinical Center of Nis).",
    co_medication   = "All patients were on dual antiepileptic therapy. Valproic acid was given with lamotrigine (VA/LTG, n = 42) or with levetiracetam (VA/LEV, n = 20). 'Any influence from the coadministration of the other two drugs, LEV and LTG, was not found to be statistically significant' (Discussion paragraph 6), so no comedication covariate appears in this model.",
    notes           = "The whole study enrolled 71 children across three dual-therapy regimens (VA/LTG n = 42, VA/LEV n = 20, LTG/LEV n = 9; Table 1B). The VA sub-cohort is the 62 patients whose regimen contained valproic acid. That denominator is confirmed by the paper's own reporting: 93.55% of VA concentrations were within the reference range = 58/62 (Results paragraph 1). Exactly one steady-state trough concentration was drawn per patient, immediately before the next morning dose, so all parameters were identified from a single sample per subject via a stepwise fix-and-release estimation strategy (Methods, 'Population pharmacokinetics'). Valproate was the only one of the three drugs under routine therapeutic drug monitoring at the centre."
  )

  ini({
    # Structural parameters - final-model estimates from Damnjanovic 2023 Table 2(c).

    # Ka FIXED at the literature value: "Given that the VA steady-state trough
    # plasma concentration data collected did not provide information about the
    # amount and pace of absorption processes, Ka was set at 1.68 h-1 based on the
    # literature Ka values." (Results, VA paragraph). Table 2(c) prints no SE / RSE.
    lka <- fixed(log(1.68));  label("Absorption rate constant (1/h)")                  # Table 2(c): Ka = 1.68 1/h, no SE/RSE reported
    lvc <- log(15.61);        label("Apparent central volume of distribution (L)")    # Table 2(c): V = 15.61 L (SE 4.84, RSE 31.0%)
    lcl <- log(0.12);         label("Apparent oral clearance (L/h)")                 # Table 2(c): Cl = 0.12 L/h (SE 0.013, RSE 10.4%)

    # Allometric exponent on body weight, FIXED (Table 2(c) prints "-" for SE and
    # RSE). The p < 0.001 is the Wald test on the covariate's inclusion.
    e_wt_vc <- fixed(1);      label("Allometric exponent: WT on V/F (unitless)")           # Table 2(c): beta_V_logBW = 1, p < 0.001

    # Age on V/F, estimated, centred on the 10.9-year cohort mean.
    e_age_vc <- 0.07;         label("Effect of age on V/F, centred at 10.9 years (per year, exponential-linear)")  # Table 2(c): beta_V_Age = 0.07 (SE 0.016, RSE 23.3%, p = 0.032)

    # Valproate's own total daily dose on CL/F, exponential-linear and untransformed
    # (see the covariateData note and the vignette Errata for why the covariate is
    # not mean-centred: the paper never reports a mean daily dose).
    e_dose_vpa_cl <- 0.0012;  label("Effect of total daily valproate dose on CL/F (per mg/day, exponential-linear)")  # Table 2(c): beta_Cl_VA = 0.0012 (SE 0.0001, RSE 9.7%, p < 0.001)

    # IIV. Monolix reports omega as the STANDARD DEVIATION of the log-scale random
    # effect; nlmixr2's ini() takes VARIANCES, so each omega is squared:
    #   omega_V  = 0.33  -> var = 0.1089
    #   omega_Cl = 0.089 -> var = 0.007921
    # Table 2(c) reports no V-Cl correlation for valproic acid, so the two etas are
    # independent here.
    etalvc ~ 0.1089                                                                       # Table 2(c): omega_V = 0.33 (SE 0.081, RSE 24.5%)
    etalcl ~ 0.007921                                                                     # Table 2(c): omega_Cl = 0.089 (SE 0.024, RSE 26.7%)

    # Residual error. "The residual variability was estimated using a proportional
    # error model." (Results, VA paragraph). Monolix's proportional error model is
    # y = f + b * f * eps, so b is the proportional SD as a fraction.
    propSd <- 0.14;           label("Proportional residual error (fraction)")              # Table 2(c): b = 0.14 (SE 0.037, RSE 27.1%)
  })

  model({
    # No IIV on Ka: Table 2(c) reports no omega_Ka, consistent with Ka being fixed.
    ka <- exp(lka)

    # V/F: fixed allometric weight term times a mean-centred exponential age term.
    vc <- exp(lvc + etalvc) *
      (WT / 37.1)^e_wt_vc *
      exp(e_age_vc * (AGE - 10.9))

    # CL/F covariate model on the log-parameter scale, as Monolix parameterises a
    # lognormal parameter: log(Cl_i) = log(Cl_pop) + beta_dose * DailyDose + eta_i
    cl <- exp(lcl + etalcl) * exp(e_dose_vpa_cl * DOSE_VPA_MGD)

    kel <- cl / vc

    # One compartment, first-order oral absorption and first-order elimination.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in mg, V/F in L -> mg/L, the unit of the paper's therapeutic reference
    # range for total valproate.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
