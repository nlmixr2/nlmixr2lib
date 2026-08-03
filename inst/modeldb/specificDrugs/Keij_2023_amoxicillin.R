Keij_2023_amoxicillin <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption for oral",
    "and intravenous amoxicillin in preterm and term neonates treated for",
    "possible serious bacterial infection, pooled from three studies (RAIN,",
    "Maastricht/Pullen, and SATT). Clearance carries fixed allometric weight",
    "scaling (exponent 0.75, reference 70 kg) plus two estimated maturation",
    "power terms on postnatal age (reference 6.8 days) and gestational age",
    "(reference 35.8 weeks); central volume carries fixed linear weight",
    "scaling (exponent 1.00, reference 70 kg). Oral bioavailability is 87.3%",
    "and absorption is slow (Ka 0.085 1/h), so oral profiles are",
    "absorption-rate-limited. Interindividual variability is on clearance",
    "only; residual error is combined additive plus proportional",
    "(Keij 2023)."
  )
  reference <- paste(
    "Keij FM, Schouwenburg S, Kornelisse RF, Preijers T, Mir F, Degraeuwe P,",
    "Stolk LM, van Driel A, Kenter S, van der Sluijs J, Heidema J,",
    "den Butter PCP, Reiss IKM, Allegaert K, Tramper-Stranders GA,",
    "Koch BCP, Flint RB. Oral and Intravenous Amoxicillin Dosing",
    "Recommendations in Neonates: A Pooled Population Pharmacokinetic Study.",
    "Clin Infect Dis. 2023;77(11):1595-1603. doi:10.1093/cid/ciad432"
  )
  vignette <- "Keij_2023_amoxicillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Current (not birth) body weight of the neonate.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column BW ('current bodyweight'). Enters both clearance and",
        "central volume as an allometric power term normalized to 70 kg",
        "(Keij 2023 Table 2 equations and Table 2 footnote: 'Current",
        "bodyweight is scaled to 70 kg'). Both exponents were FIXED, not",
        "estimated: 0.75 on CL and 1.00 on V. Supplementary Appendix B",
        "('Model Development') states that allometric scaling was tested",
        "several ways and that the 'fixed power' approach was the one",
        "included in the final model. The 70 kg reference is far outside the",
        "observed range (0.5-5.0 kg, median 2.6 kg), so exp(lcl) = 3.22 L/h",
        "and exp(lvc) = 43 L are extrapolated adult-size anchors rather than",
        "neonatal typical values; at the cohort-median 2.6 kg they",
        "correspond to CL = 0.27 L/h and V = 1.60 L. Time-varying in",
        "principle; the source datasets recorded weight during therapy.",
        "The neonatal cohort range is 0.5-5.0 kg, so simulations far above",
        "that range are extrapolation."
      ),
      source_name        = "BW"
    ),
    PNA = list(
      description        = "Postnatal (chronological) age since birth.",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Keij 2023 reports postnatal age in DAYS (Table 2 clearance equation",
        "uses (PNA / 6.8) with PNA in days; the 6.8-day reference is the",
        "pooled dataset median per the Table 2 footnote 'Both postnatal age",
        "and gestational age are scaled to dataset median'). The canonical",
        "nlmixr2lib PNA column carries MONTHS, so model() reparameterises",
        "the reference as 6.8 / 30.4375 = 0.2234 months; numerator and",
        "denominator carry the same units factor, so the ratio and hence the",
        "exponent are unchanged. This follows the Zhao_2018_omeprazole",
        "precedent. Users supply PNA in months.",
        "IMPORTANT: the power form (PNA / 6.8)^0.357 is singular at PNA = 0",
        "-- clearance goes to zero on the day of birth. The authors' own",
        "simulation cohorts include PNA = 0 patients (Supplementary Table",
        "S2: 'Lowest 1: GA 25, PNA 0, weight 700g'), so the published model",
        "cannot be evaluated at PNA = 0 as written. Use PNA >= 1 day",
        "(0.0329 months) when simulating. Studied range was PNA 0-59 days",
        "(pooled median 1 day, IQR 0-4; the SATT cohort extends to 59 days)."
      ),
      source_name        = "PNA"
    ),
    GA = list(
      description        = "Gestational age at birth. Time-fixed per subject.",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters clearance as (GA / 35.8)^2.37, where 35.8 weeks is the",
        "pooled dataset median gestational age quoted in the Keij 2023",
        "abstract and used as the centering value in the Table 2 clearance",
        "equation (Table 2 footnote: 'Both postnatal age and gestational age",
        "are scaled to dataset median'). Note that Table 1 reports a median",
        "GA of 37.4 weeks [IQR 31.7-39.86] for the pooled cohort; the",
        "equation's centering constant is nevertheless printed as 35.8, and",
        "the equation is authoritative. Observed range 24.9-42.4 weeks. The",
        "steep exponent (2.37) means clearance is highly sensitive to GA:",
        "a 25-week neonate has roughly one fifth the size-normalized",
        "clearance of a 41-week neonate."
      ),
      source_name        = "GA"
    )
  )

  # Covariates that Keij 2023 screened during covariate analysis but did NOT
  # retain in the final model. Documentation only -- no published point
  # estimates exist for these, so they carry no effect in model().
  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Female sex indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Tested as a categorical covariate (proportional model) during",
        "stepwise forward inclusion / backward elimination and not retained",
        "in the final model (Keij 2023 Results, 'Final Population",
        "Pharmacokinetic Model'; Supplementary Figure S7 shows sex vs",
        "ETACL boxplots). Cohort was 42.9% female (Table 1). No point",
        "estimate is published, so no effect is encoded."
      ),
      source_name        = "SEX"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age plus postnatal age).",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tested as a continuous covariate and not retained: the authors",
        "found that GA and PNA entered SEPARATELY described maturation of",
        "amoxicillin clearance better than the combined postmenstrual-age",
        "term (Keij 2023 Results: 'GA and PNA were found to best describe",
        "maturation of amoxicillin clearance'). Pooled median PMA was 38.29",
        "weeks [IQR 31.9-40.71] (Table 1). No point estimate is published,",
        "so no effect is encoded."
      ),
      source_name        = "PMA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 261,
    n_studies      = 3,
    n_observations = 938,
    age_range      = "postnatal age 0-59 days (pooled median 1 day, IQR 0-4)",
    ga_range       = "gestational age 24.9-42.4 weeks (median 35.8 weeks per the abstract; Table 1 reports median 37.4 weeks, IQR 31.7-39.86)",
    pma_range      = "postmenstrual age median 38.29 weeks (IQR 31.9-40.71)",
    weight_range   = "0.5-5.0 kg",
    weight_median  = "2.6 kg (IQR 1.6-3.5)",
    sex_female_pct = 42.9,
    disease_state  = "preterm and term neonates treated for possible serious bacterial infection (pSBI)",
    dose_range     = "intravenous median 50 mg/kg/dose (range 9.4-112.9); oral median 78.9 mg/kg/dose (range 23.4-100)",
    routes         = "oral (79 neonates, 123 concentrations) and intravenous (182 neonates, 815 concentrations); 7 patients switched IV to oral",
    regions        = "The Netherlands (RAIN and Maastricht cohorts) and Pakistan (SATT cohort, Karachi)",
    notes          = paste(
      "Pooled analysis of three datasets (Keij 2023 Table 1): (1) RAIN",
      "(Reduction of intravenous Antibiotics In Neonates), n = 39, oral,",
      "PMA >= 35 weeks and weight >= 2 kg, all co-administered",
      "amoxicillin + clavulanic acid 4:1; (2) 'Maastricht' (Pullen et al.),",
      "n = 182, intravenous, PNA <= 9 days, co-medication indomethacin;",
      "(3) SATT (Simplified Antibiotic Therapy Trial, Karachi, Pakistan),",
      "n = 40, oral, PNA 0-59 days, co-medication unknown.",
      "Race/ethnicity was not reported. Modelling was in NONMEM 7.4.",
      "Covariates screened and NOT retained: sex, postmenstrual age,",
      "study centre (which indirectly encodes the clavulanic-acid and",
      "indomethacin co-medication differences between cohorts), and",
      "small-for-gestational-age status; see Supplementary Figures S6-S7.",
      "Only gestational age and postnatal age were retained on clearance",
      "(dOFV -100 and -123 respectively, both P < .001).",
      "Free (unbound) concentrations were NOT measured and total",
      "concentrations were NOT corrected for protein binding; the authors",
      "note neonatal amoxicillin protein binding is only 10-14%, so the",
      "paper's %fT>MIC target attainment is computed on total",
      "concentrations. 13 concentrations above the ULOQ and 5 below the",
      "LLOQ were retained in the dataset rather than excluded."
    )
  )

  ini({
    # --- Structural parameters (reference: WT 70 kg, PNA 6.8 d, GA 35.8 wk) ---
    lka <- log(0.085)
    label("Absorption rate constant (1/h)")
    # Keij 2023 Table 2, 'Ka (h-1)': 0.085 (RSE 25%); bootstrap 0.091 (0.06-0.11)

    lcl <- log(3.22)
    label("Typical clearance TVCL at WT 70 kg, PNA 6.8 d, GA 35.8 wk (L/h)")
    # Keij 2023 Table 2, 'TVCL (L/h)': 3.22 (RSE 3%); bootstrap 3.22 (3.09-3.35)

    lvc <- log(43)
    label("Typical central volume TVV1 at WT 70 kg (L)")
    # Keij 2023 Table 2, 'TVV1 (L)': 43 (RSE 2%); bootstrap 43 (41.61-44.32)

    lfdepot <- log(0.873)
    label("Oral bioavailability (fraction)")
    # Keij 2023 Table 2, 'F (%)': 0.873 (RSE 16%); bootstrap 0.832 (0.72-1.02).
    # Abstract and Discussion both quote this as 87% oral bioavailability.

    # --- Allometric exponents: FIXED, not estimated ---
    # Supplementary Appendix B, 'Model Development': allometric scaling was
    # tested several ways and the 'fixed power (included in the final model)'
    # approach was retained. Table 2 prints both exponents inline in the
    # equations with no RSE and no bootstrap CI, confirming they were fixed.
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on clearance (unitless)")
    # Keij 2023 Table 2 equation: CL_A = TVCL * (BW/70)^0.75 * ...

    e_wt_vc <- fixed(1.00)
    label("Allometric exponent of body weight on central volume (unitless)")
    # Keij 2023 Table 2 equation: V_A = TVV1 * (BW/70)^1.00

    # --- Estimated maturation covariate effects on clearance ---
    e_pna_cl <- 0.357
    label("Postnatal-age power exponent on clearance (unitless; reference 6.8 days)")
    # Keij 2023 Table 2, 'theta_PNA': 0.357 (RSE 8%); bootstrap 0.359 (0.30-0.41)

    e_ga_cl <- 2.37
    label("Gestational-age power exponent on clearance (unitless; reference 35.8 weeks)")
    # Keij 2023 Table 2, 'theta_GA': 2.37 (RSE 6%); bootstrap 2.37 (2.13-2.61)

    # --- Interindividual variability (clearance only) ---
    etalcl ~ 0.0713
    # Keij 2023 Table 2, 'Interindividual variability / Clearance (%)': 26.7%
    # (shrinkage 21%). The bootstrap column for this row reads
    # '26.7 (.05-.09)', i.e. the CI is reported on the VARIANCE scale while
    # the point estimate is reported as a CV percentage. omega^2 = 0.267^2 =
    # 0.0713, which falls inside the bootstrap interval (0.05-0.09) and so
    # confirms the CV-percent reading. (The strict log-normal conversion
    # log(CV^2 + 1) = 0.0689 also falls in that interval; the difference is
    # immaterial at this CV and the direct-square form is used here.)
    # No IIV was reported on Ka, V, or F, so none is encoded.

    # --- Residual unexplained variability (combined / 'mixed' error model) ---
    propSd <- 0.132
    label("Proportional residual error (fraction)")
    # Keij 2023 Table 2, 'Proportional error': 0.132 (RSE 10%);
    # bootstrap 0.133 (0.11-0.15). Reported as an SD on the fractional scale.

    addSd <- 4.48
    label("Additive residual error (mg/L)")
    # Keij 2023 Table 2, 'Additive error (mg/L)': 4.48 (RSE 12%);
    # bootstrap 4.38 (3.59-5.37). The mg/L unit in the row label confirms
    # this is an SD in concentration units rather than a variance.
    # Results text: 'A mixed error model was used to describe the residual
    # variability.'
  })

  model({
    # 1. Derived covariate terms
    #    Keij 2023 Table 2 (and Supplementary Appendix B):
    #      CL_A (L/h) = TVCL * (BW/70)^0.75 * (PNA/6.8)^theta_PNA * (GA/35.8)^theta_GA
    #      V_A  (L)   = TVV1 * (BW/70)^1.00
    #    PNA reference is 6.8 DAYS in the paper; the canonical PNA column is
    #    in MONTHS, so the reference is converted once here. The ratio is
    #    unit-invariant, so the exponent transfers unchanged.
    pna_ref_months <- 6.8 / 30.4375
    f_pna <- (PNA / pna_ref_months)^e_pna_cl
    f_ga <- (GA / 35.8)^e_ga_cl

    # 2. Individual parameters
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * f_pna * f_ga
    vc <- exp(lvc) * (WT / 70)^e_wt_vc

    # 3. Micro-constants
    kel <- cl / vc

    # 4. ODE system (one compartment with a first-order absorption depot)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 5. Bioavailability -- applies to oral (depot) doses only. Intravenous
    #    doses are administered directly into `central` and therefore carry
    #    F = 1, matching the study design (Keij 2023 Table 1: RAIN and SATT
    #    oral, Maastricht intravenous).
    f(depot) <- exp(lfdepot)

    # 6. Observation and error
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
