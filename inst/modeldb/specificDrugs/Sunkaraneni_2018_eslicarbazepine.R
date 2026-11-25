Sunkaraneni_2018_eslicarbazepine <- function() {
  description <- paste0(
    "One-compartment population PK model for ESLICARBAZEPINE, the active ",
    "metabolite of the once-daily antiepileptic prodrug eslicarbazepine ",
    "acetate (ESL), in pediatric patients aged 2-18 years with ",
    "partial-onset seizures (Sunkaraneni 2018, n = 146, 857 plasma ",
    "concentrations from studies BIA-2093-202 and BIA-2093-305). ",
    "Formulation-specific first-order absorption (tablet versus oral ",
    "suspension, both rate constants fixed from the rich study 202 data) ",
    "and first-order elimination. Apparent clearance and volume are ",
    "referenced to the ESL dose and scaled allometrically to a 33 kg ",
    "subject with fixed exponents 0.75 and 1. Concomitant levetiracetam ",
    "lowers clearance by 17.6%, phenobarbital-like inducers (phenobarbital, ",
    "phenytoin, primidone) raise it by 62.6%, and concomitant carbamazepine ",
    "lowers relative bioavailability to 0.679. Combined additive and ",
    "proportional residual error on the ng/mL scale."
  )
  reference <- paste(
    "Sunkaraneni S, Ludwig E, Fiedler-Kelly J, Hopkins S, Galluppi G,",
    "Blum D. Modeling and simulations to support dose selection for",
    "eslicarbazepine acetate therapy in pediatric patients with",
    "partial-onset seizures. J Pharmacokinet Pharmacodyn.",
    "2018;45(4):649-658. doi:10.1007/s10928-018-9596-7.",
    sep = " "
  )
  vignette <- "Sunkaraneni_2018_eslicarbazepine"

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of apparent clearance with fixed exponent 0.75",
        "and of apparent volume with fixed exponent 1, both centred at",
        "33.0 kg (the printed typical-value equations under Table 3).",
        "Observed range about 10-79 kg (Table 2 weight categories); mean",
        "34.7 kg, SD 16.0. Freely estimated exponents (0.648 on CL, 1.19",
        "on V) were tested and rejected in favour of the theoretical",
        "values (Discussion)."
      ),
      source_name = "WTKG"
    ),
    CONMED_LEV = list(
      description = "Concomitant levetiracetam; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no levetiracetam)",
      notes = paste(
        "Proportional shift -0.176 on apparent clearance,",
        "CL * (1 - 0.176 * LEVEF). The authors caution that this small",
        "effect, absent in adults, may be an artefact of sample size and",
        "not a true difference (Discussion)."
      ),
      source_name = "LEVEF"
    ),
    CONMED_PB = list(
      description = "Concomitant phenobarbital; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no phenobarbital)",
      notes = paste(
        "One of the three constituents of the paper's composite",
        "phenobarbital-like AED flag PHENLF ('phenobarbital, phenytoin, and",
        "primidone'). The flag is rebuilt in model() as the logical OR",
        "1 - (1 - CONMED_PB) * (1 - CONMED_PHT) * (1 - CONMED_PRM), the same",
        "decomposition as Gidal_2018_eslicarbazepine.R. Carbamazepine is",
        "NOT part of the flag; it has its own bioavailability term."
      ),
      source_name = "PHENLF (component)"
    ),
    CONMED_PHT = list(
      description = "Concomitant phenytoin; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no phenytoin)",
      notes = paste(
        "Second constituent of the composite PHENLF flag; see the",
        "CONMED_PB notes."
      ),
      source_name = "PHENLF (component)"
    ),
    CONMED_PRM = list(
      description = "Concomitant primidone; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no primidone)",
      notes = paste(
        "Third constituent of the composite PHENLF flag; see the",
        "CONMED_PB notes."
      ),
      source_name = "PHENLF (component)"
    ),
    CONMED_CBZ = list(
      description = "Concomitant carbamazepine; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no carbamazepine)",
      notes = paste(
        "Acts on relative bioavailability only: F1 = 0.679 during",
        "carbamazepine use, 1 otherwise (Table 3), i.e. 32.1% lower",
        "bioavailability and AUC0-24,ss (Results). Unlike the adult model",
        "of Gidal 2018 there is no carbamazepine dose term."
      ),
      source_name = "F1 (carbamazepine use)"
    ),
    FORM_TABLET = list(
      description = paste(
        "ESL formulation; 1 = oral tablet, 0 = oral suspension."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (oral suspension)",
      notes = paste(
        "Selects between the two fixed first-order absorption rate",
        "constants, KAT = 0.895 1/h for the tablet and KAO = 4.18 1/h for",
        "the oral suspension (Table 3). The comparator oral liquid is a",
        "SUSPENSION. Formulation has no effect on bioavailability; the",
        "formulations are bioequivalent (Discussion)."
      ),
      source_name = "formulation (KAT / KAO)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = "Screened on CL and V in the stepwise covariate search (Methods); not retained."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate.",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Screened on CL and V (Methods, reference [16] for the estimating equation); not retained."
    ),
    HT = list(
      description = "Height.",
      units = "cm",
      type = "continuous",
      notes = "Screened on CL and V (Methods); not retained."
    ),
    RACE_ASIAN = list(
      description = "Asian race; 1 = Asian, 0 = otherwise (all others were White).",
      units = "(binary)",
      type = "binary",
      notes = "Race screened on CL and V (Methods); not retained. 6.8% Asian (Table 2)."
    ),
    SEXF = list(
      description = "Sex; 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL and V (Methods); not retained, unlike the adult model where sex acts on V/F."
    ),
    CONMED_VPA = list(
      description = "Concomitant valproic acid; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      notes = "Screened (used by > 10% of the analysis population, Methods); not retained."
    ),
    CONMED_LAMOTRIGINE = list(
      description = "Concomitant lamotrigine; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      notes = "Screened (used by > 10% of the analysis population, Methods); not retained."
    ),
    CONMED_TPM = list(
      description = "Concomitant topiramate; 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      notes = "Screened (used by > 10% of the analysis population, Methods); not retained."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "eslicarbazepine",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "eslicarbazepine",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 146L,
    n_studies = 2L,
    n_observations = "857 plasma eslicarbazepine concentrations",
    age_range = "2-18 years",
    age_median = "mean 10.0 years, SD 4.0 (Table 2)",
    weight_range = "about 10-79 kg (Table 2 weight categories)",
    weight_median = "mean 34.7 kg, SD 16.0 (Table 2); reference 33.0 kg",
    sex_female_pct = 53.4,
    race_ethnicity = c(White = 93.2, Asian = 6.8),
    disease_state = "pediatric patients with partial-onset seizures on adjunctive ESL therapy",
    dose_range = paste(
      "ESL orally once daily, 5-30 mg/kg/day; study 202 (n = 26 in the",
      "analysis) 5, 15 and 30 mg/kg/day in consecutive 4-week periods",
      "(maximum 1800 mg/day) with rich sampling 0.5-24 h post-dose; study",
      "305 (n = 120) titration 10 and 20 mg/kg/day and maintenance 10, 20",
      "or 30 mg/kg/day (maximum 1200 mg/day) with about 3 sparse samples",
      "per subject"
    ),
    formulations = "oral tablet and oral suspension",
    co_medication = paste(
      "concomitant AEDs used by > 10% of the population were screened:",
      "carbamazepine, valproic acid, lamotrigine, topiramate,",
      "levetiracetam and phenobarbital-like inducers (phenobarbital,",
      "primidone, phenytoin)"
    ),
    notes = paste(
      "Estimation in NONMEM 7.1.2. The structural model is the adult",
      "one-compartment adjunctive-therapy model refined for children.",
      "Model evaluation by prediction-corrected VPC (Fig. 2) and",
      "non-parametric bootstrap (Table 3)."
    )
  )

  ini({
    # Table 3 and the typical-value equations printed beneath it:
    #   CL_i = 1.69 * (WTKG_i / 33.0)^0.75 * (1 - 0.176 * LEVEF_i)
    #                                      * (1 + 0.626 * PHENLF_i)
    #   V_i  = 32.8 * WTKG_i / 33.0
    lka_tablet <- fixed(log(0.895)) ; label("First-order absorption rate constant, tablet (1/h)")  # Table 3, 'KAT' 0.895 1/h, Fixed (from study 202, Discussion)
    lka_susp <- fixed(log(4.18)) ; label("First-order absorption rate constant, oral suspension (1/h)")  # Table 3, 'KAO' 4.18 1/h, Fixed (from study 202, Discussion)
    lcl <- log(1.69) ; label("Apparent clearance for a 33 kg subject without levetiracetam or phenobarbital-like AEDs (L/h)")  # Table 3, 'CL' 1.69 L/h, 2.92 %SEM; CL equation
    lvc <- log(32.8) ; label("Apparent volume of distribution for a 33 kg subject (L)")  # Table 3, 'V' 32.8 L, 4.78 %SEM; V equation

    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of (WT / 33 kg) on clearance (unitless)")  # CL equation and Discussion, exponent fixed at 0.75
    e_wt_vc <- fixed(1) ; label("Allometric exponent of (WT / 33 kg) on volume (unitless)")  # V equation and Discussion, exponent fixed at 1

    e_conmed_lev_cl <- -0.176 ; label("Proportional shift in clearance with concomitant levetiracetam (fraction)")  # Table 3, 'CL: proportional shift for levetiracetam use' -0.176, 25.6 %SEM
    e_pblike_cl <- 0.626 ; label("Proportional shift in clearance with phenobarbital-like AEDs (fraction)")  # Table 3, 'CL: proportional shift for use of phenobarbital-like AEDs' 0.626, 18.8 %SEM
    e_conmed_cbz_fdepot <- 0.679 ; label("Relative bioavailability during concomitant carbamazepine use (fraction)")  # Table 3, 'F1: relative bioavailability during carbamazepine use' 0.679, 6.76 %SEM

    # IIV printed as %CV; variance on the log scale omega^2 = log(CV^2 + 1).
    # One KA eta shared by both formulations ('variance assumed to be the
    # same for the oral suspension and tablet formulations', Results).
    etalka ~ fixed(0.531947)  # Table 3, KAT / KAO IIV 83.8 %CV, no SEM: log(0.838^2 + 1) = 0.531947
    etalcl ~ 0.060625  # Table 3, CL IIV 25.0 %CV, 15.5 %SEM: log(0.25^2 + 1) = 0.060625
    etalvc ~ 0.017274  # Table 3, V IIV 13.2 %CV, 65.1 %SEM: log(0.132^2 + 1) = 0.017274

    # Residual: footnote a gives %CV = sqrt(F^2 * 0.0543 + 107,000) / F * 100,
    # so both printed values are variances on the ng/mL scale.
    propSd <- 0.233024 ; label("Proportional residual error (fraction)")  # Table 3, 'RV CCV component' 0.0543 (variance), 11.6 %SEM: sqrt(0.0543) = 0.233024
    addSd <- 327.1085 ; label("Additive residual error (ng/mL)")  # Table 3, 'RV additive component' 107,000 (variance, ng^2/mL^2), 53.2 %SEM: sqrt(107000) = 327.1085
  })

  model({
    # Composite phenobarbital-like AED flag (PHENLF): logical OR of the three
    # constituent comedications.
    pblike <- 1 - (1 - CONMED_PB) * (1 - CONMED_PHT) * (1 - CONMED_PRM)

    ka <- exp(lka_tablet * FORM_TABLET + lka_susp * (1 - FORM_TABLET) + etalka)
    cl <- exp(lcl + etalcl) * (WT / 33)^e_wt_cl *
      (1 + e_conmed_lev_cl * CONMED_LEV) *
      (1 + e_pblike_cl * pblike)
    vc <- exp(lvc + etalvc) * (WT / 33)^e_wt_vc

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # F1 = 0.679 with carbamazepine, 1 without.
    f(depot) <- 1 + (e_conmed_cbz_fdepot - 1) * CONMED_CBZ

    # Dose in mg and vc in L give mg/L (ug/mL); x1000 gives the ng/mL scale of
    # the residual-error model (footnote a, F range 100-50,000).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
