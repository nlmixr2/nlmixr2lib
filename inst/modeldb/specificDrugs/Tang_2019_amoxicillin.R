Tang_2019_amoxicillin <- function() {
  description <- "Two-compartment population PK model with first-order elimination for intravenous amoxicillin in Chinese neonates and young infants (Tang 2019). Current weight enters as a fixed allometric power on both volumes (exponent 1) and on CL and Q (exponent 0.75); CL is further modulated by a maturation factor F_age that is the product of two power functions of gestational age and postnatal age. Interindividual variability is estimated on the peripheral volume V2 and on CL only; residual variability follows an exponential model (proportional in linear space)."
  reference <- paste(
    "Tang B-H, Wu Y-E, Kou C, Qi Y-J, Qi H, Xu H-Y, Leroux S, Huang X, Zhou Y,",
    "Zheng Y, Jacqz-Aigrain E, Shen A-D, Zhao W.",
    "Population Pharmacokinetics and Dosing Optimization of Amoxicillin in Neonates",
    "and Young Infants. Antimicrob Agents Chemother. 2019;63(2):e02336-18.",
    "doi:10.1128/AAC.02336-18.",
    sep = " "
  )
  vignette <- "Tang_2019_amoxicillin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Current body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying current weight on the day of the PK sample. Tang 2019 reports the column in grams (median 3,210 g; range 1,060 to 4,580 g) and uses 3,210 g as the allometric reference. The canonical WT column is in kg, so the reference is rescaled to 3.21 kg and the (CW / 3210)^a relationships from Table 2 become (WT / 3.21)^a inside model() with no change in the exponents.",
      source_name        = "CW"
    ),
    GA = list(
      description        = "Gestational age at birth",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Tang 2019 cohort median 38.14 weeks (range 28.3 to 41.4 weeks) is used as the reference inside the F_age maturation factor for CL.",
      source_name        = "GA"
    ),
    PNA = list(
      description        = "Postnatal age (chronological time since birth)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Tang 2019 reports PNA in DAYS (cohort median 7 days, range 1 to 37 days) and parameterises F_age as (PNA_days / 7)^0.28. The canonical PNA column is in MONTHS, so the reference is rescaled to 7 / 30.4375 = 0.2300 months and the model code uses (PNA_months / 0.2300)^0.28; numerator and denominator carry the same unit so the dynamic relationship and the estimated exponent are unchanged. Same rescaling precedent as Zhao 2018 omeprazole.",
      source_name        = "PNA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 187L,
    n_studies      = 1L,
    n_observations = "224 amoxicillin plasma concentrations (Tang 2019 Results 'Model building'); 18 below the LLOQ of 0.5 ug/mL were imputed at LLOQ/2 = 0.25 ug/mL by the source authors.",
    age_range      = "Neonates and young infants; postmenstrual age 28.4 to 46.3 weeks (Tang 2019 Table 1)",
    age_median     = "PMA median 39.0 weeks; PNA median 7 days (Table 1)",
    weight_range   = "Current weight 1.06 to 4.58 kg (Tang 2019 Table 1)",
    weight_median  = "Current weight 3.21 kg (Table 1)",
    sex_female_pct = 47.6,  # 89 female / 187 total per Table 1 footnote a
    race_ethnicity = c(Asian = 100),
    disease_state  = "Chinese neonates and young infants admitted to neonatal intensive care units who received amoxicillin as part of regular antimicrobial treatment; postmenstrual age < 48 weeks (Tang 2019 Methods 'Study design').",
    dose_range     = "Intravenous amoxicillin 25 mg/kg twice daily administered as a bolus over 5 minutes or as an infusion over 30 minutes; unit-dose median 23.7 mg/kg/dose (range 15.7 to 35.2 mg/kg/dose) per Table 1.",
    regions        = "China (Beijing Obstetrics and Gynecology Hospital, Beijing Children's Hospital, Shandong Provincial Qianfoshan Hospital).",
    notes          = "Demographics from Tang 2019 Table 1 and Methods 'Study design' / 'Analytical method for amoxicillin'. Opportunistic sampling design with HPLC quantification (LLOQ 0.5 ug/mL, intra- and interday CV 3.05% and 4.30%). Free fraction not measured; serum protein binding considered negligible (~ 10%). NONMEM v7.2 with FOCE-I; PsN v2.30 for bootstrap. The 'Final model' column of Table 2 is used here (rather than the 'Total model' column of Table 3 that pools the original 187-subject dataset with 48 external-validation subjects); the parameter values differ slightly between the two columns and Tang 2019 reports the 187-subject final-model estimates as the primary analysis."
  )

  ini({
    # Final-model estimates from Tang 2019 Table 2 ('Final estimate' column).
    # Bootstrap medians and 5th-95th-percentile intervals are reported
    # alongside in Table 2 and agree closely with the final estimates,
    # indicating model stability. Volumes and clearances are in litres and
    # L/h respectively; weight is normalised against the cohort median
    # CW = 3,210 g = 3.21 kg. The 'Total model' values in Table 3 (which
    # pool the 187-subject development cohort with the 48-subject external
    # validation cohort) are NOT used here -- the primary published model
    # is the development-cohort final model from Table 2.
    #
    # Source-paper definitions reproduced here for traceability:
    #   V1 = theta1 * (CW / 3210)
    #   V2 = theta2 * (CW / 3210)
    #   Q  = theta3 * (CW / 3210)^0.75
    #   CL = theta4 * (CW / 3210)^0.75 * F_age
    #   F_age = (GA / 38.14)^theta5 * (PNA / 7)^theta6
    # (CW in grams, GA in weeks, PNA in days; this skill reparameterises
    # CW to kg and PNA to months -- see covariateData notes.)

    # Structural PK parameters (reference 3.21 kg, GA = 38.14 weeks, PNA = 7 days).
    lvc <- log(1.48);  label("Central volume of distribution V1 at 3.21 kg reference weight (L)")            # Tang 2019 Table 2: theta1 = 1.48 L (RSE 7.80%)
    lvp <- log(2.42);  label("Peripheral volume of distribution V2 at 3.21 kg reference weight (L)")         # Tang 2019 Table 2: theta2 = 2.42 L (RSE 28.30%)
    lq  <- log(0.17);  label("Inter-compartmental clearance Q at 3.21 kg reference weight (L/h)")            # Tang 2019 Table 2: theta3 = 0.17 L/h (RSE 23.70%)
    lcl <- log(0.81);  label("Clearance CL at 3.21 kg reference weight, GA = 38.14 wk, PNA = 7 d (L/h)")     # Tang 2019 Table 2: theta4 = 0.81 L/h (RSE 7.50%)

    # Allometric weight exponents -- held fixed at canonical 1.0 (volumes)
    # and 0.75 (clearances) per Tang 2019 Results 'Covariate analysis':
    # "The current weight (CW) was a priori incorporated into the basic
    # model by the allometric size approach (allometric coefficients were
    # 1 for V1 and V2 and 0.75 for CL and Q)." No RSE is reported for
    # these exponents, consistent with fixed allometric scaling.
    e_wt_vc_vp <- fixed(1.00); label("Allometric (WT) exponent shared across V1 and V2 (unitless, fixed)")   # Tang 2019 Results 'Covariate analysis' (a priori fixed)
    e_wt_cl_q  <- fixed(0.75); label("Allometric (WT) exponent shared across CL and Q (unitless, fixed)")    # Tang 2019 Results 'Covariate analysis' (a priori fixed)

    # F_age maturation exponents on CL (estimated).
    e_ga_cl  <- 4.19;  label("Gestational-age power exponent on CL (unitless; GA reference 38.14 weeks)")    # Tang 2019 Table 2: theta5 = 4.19 (RSE 12.60%)
    e_pna_cl <- 0.28;  label("Postnatal-age power exponent on CL (unitless; PNA reference 7 days = 0.2300 months)") # Tang 2019 Table 2: theta6 = 0.28 (RSE 20.30%)

    # Inter-individual variability (exponential model in Tang 2019;
    # log-normal variance via omega^2 = log(1 + CV^2)).
    # V2: CV 80.00% -> log(1 + 0.80^2) = log(1.64) = 0.49470
    etalvp ~ 0.49470  # Tang 2019 Table 2: IIV(V2) = 80.00% CV (RSE 55.00%)
    # CL: CV 40.00% -> log(1 + 0.40^2) = log(1.16) = 0.14842
    etalcl ~ 0.14842  # Tang 2019 Table 2: IIV(CL) = 40.00% CV (RSE 33.50%)

    # Residual error. Tang 2019 Results 'Model building': "Residual
    # variability was best described by an exponential model." A NONMEM
    # exponential residual on log-transformed observations is equivalent to
    # a proportional residual on the linear scale; the reported 35.00% CV
    # maps directly to propSd = 0.35.
    propSd <- 0.35;  label("Proportional residual error (fraction)")  # Tang 2019 Table 2: residual = 35.00% CV (RSE 13.60%)
  })

  model({
    # Maturation factor F_age on CL.
    #
    # Tang 2019 source form (Table 2 / Table 3 footnote a):
    #   F_age = (GA_weeks / 38.14)^4.19 * (PNA_days / 7)^0.28
    # Reparameterised to canonical PNA (months) without changing the
    # estimated exponent: 7 days / (30.4375 days / month) = 0.2300 months,
    # so (PNA_months / 0.2300)^0.28 is equivalent to (PNA_days / 7)^0.28.
    # GA is already in canonical units (weeks).
    pna_ref_months <- 7 / 30.4375  # 0.2300 months
    f_age <- (GA / 38.14)^e_ga_cl * (PNA / pna_ref_months)^e_pna_cl

    # Individual PK parameters. Volumes scale with WT^1 and clearances
    # with WT^0.75 relative to the 3.21 kg cohort-median reference; CL
    # carries the F_age maturation factor in addition to the WT scaling.
    vc <- exp(lvc)          * (WT / 3.21)^e_wt_vc_vp
    vp <- exp(lvp + etalvp) * (WT / 3.21)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 3.21)^e_wt_cl_q
    cl <- exp(lcl + etalcl) * (WT / 3.21)^e_wt_cl_q * f_age

    # Micro-constants for the 2-compartment IV model.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system. Amoxicillin is administered intravenously (5-min bolus
    # or 30-min infusion); dose goes directly into the central compartment
    # with no depot / first-order absorption.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Plasma amoxicillin concentration in mg/L (= ug/mL, the unit used
    # throughout Tang 2019). Dose in mg / V in L = mg/L; the cohort
    # median V1 = 1.48 L gives Cmax ~ 50 ug/mL after a 75 mg dose, in
    # line with the observed range up to 73.6 ug/mL in Tang 2019 Figure 1.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
