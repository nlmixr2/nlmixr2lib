Kohl_2010_posaconazole <- function() {
  description <- "One-compartment population PK model for prophylactic oral posaconazole in adult allogeneic stem cell transplant recipients with hematological malignancies (Kohl 2010); ka fixed, age and concurrent diarrhea as covariates."
  reference <- "Kohl V, Muller C, Cornely OA, Abduljalil K, Fuhr U, Vehreschild JJ, Scheid C, Hallek M, Ruping MJGT. Factors Influencing Pharmacokinetics of Prophylactic Posaconazole in Patients Undergoing Allogeneic Stem Cell Transplantation. Antimicrob Agents Chemother. 2010;54(1):207-212. doi:10.1128/AAC.01027-09"
  vignette <- "Kohl_2010_posaconazole"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject; enters apparent V/F via a linear additive deviation from the cohort median (49 years).",
      source_name        = "Age"
    ),
    DIARRHEA = list(
      description        = "Concurrent clinical diarrhea indicator (1 = diarrhea, 0 = no diarrhea)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no diarrhea)",
      notes              = "Time-fixed per subject in Kohl 2010 (each patient carries a single yes/no flag over the TDM window). Enters CL/F and V/F via the shared power-form multiplier e_diarrhea_cl_vc^DIARRHEA; the equal-magnitude effect on both is algebraically equivalent to a 1.7-fold reduction in oral bioavailability F (F_with/F_without ~= 1/1.69 ~= 0.59 per Kohl 2010 Discussion).",
      source_name        = "Diarrhea"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 32,
    n_studies      = 1,
    age_range      = "17-66 years",
    age_median     = "49.5 years",
    weight_range   = "49-115 kg",
    weight_median  = "68.5 kg",
    sex_female_pct = 50.0,
    race_ethnicity = c(Caucasian = 93.8, Asian = 6.2),
    disease_state  = "Adult allogeneic hematopoietic stem cell transplant (SCT) recipients with hematological malignancies (acute myelogenous leukemia 53.1%, lymphoma 18.8%, chronic lymphocytic leukemia 12.5%, chronic myelogenous leukemia 6.2%, acute lymphocytic leukemia / idiopathic thrombocytopenia / plasma cell leukemia 3.1% each), on prophylactic oral posaconazole.",
    dose_range     = "200 mg oral suspension three times daily (standard prophylactic dose; some dose changes during therapy enabled separate identifiability of CL/F and V/F)",
    regions        = "Germany (single center: University Hospital of Cologne)",
    notes          = "Demographics per Kohl 2010 Table 1. 22/32 (68.6%) patients had diarrhea over the TDM observation window; 16/32 (50.0%) had concomitant fever. Concomitant medications: cyclosporine 81.3%, pantoprazole 81.3%, ranitidine 50.0%, tacrolimus 25.0%. 149 trough serum posaconazole concentrations across the cohort (median 5 samples / patient, range 1-12). Patients unable to take oral medication (severe oral mucositis or vomiting) were switched to intravenous antifungals and are not part of the analysis."
  )

  ini({
    # Structural parameters -- log-transformed where positive.
    # Reference subject: 49 years old (cohort median), no diarrhea.
    lka       <- fixed(log(0.4));   label("Absorption rate constant (ka, 1/h)")            # Kohl 2010 Methods assumption (iii): fixed at 0.4 /h from reference 5
    lcl       <- log(67.0);          label("Apparent clearance at reference (CL/F, L/h)")  # Kohl 2010 Table 4 final model, no diarrhea
    lvc       <- log(2250);          label("Apparent volume of distribution at reference age 49 yr (V/F, L)")  # Kohl 2010 Table 4 final model, no diarrhea, age 49
    e_age_vc          <- -123;       label("Linear-additive age effect on V/F (L per year above 49)")          # Kohl 2010 Table 4 final model: V/F decreases 123 L per year of age above the cohort median 49 yr
    e_diarrhea_cl_vc  <- 1.69;       label("Shared diarrhea multiplier on CL/F and V/F (theta_Di)")            # Kohl 2010 Table 4 final model: 113.2/67.0 = 3802.5/2250 = 1.69 for both apparent parameters (Discussion calls this a 1.7-fold reduction in oral bioavailability F; equivalently F_with/F_without ~= 1/1.69 ~= 0.59); applied as theta_Di^Diarrhea per Kohl 2010 Table 3 model 2

    # Inter-individual variability on CL/F (exponential).
    # CV = 26.9% per Kohl 2010 Table 4; omega^2 = log(CV^2 + 1) = log(0.269^2 + 1) ~= 0.06985.
    etalcl    ~ 0.06985                                                                                       # Kohl 2010 Table 4 final model: IIV CL/F = 26.9% CV

    # Residual error: proportional, CV = 42% per Kohl 2010 Table 4.
    propSd    <- 0.42;               label("Proportional residual error (fraction)")                          # Kohl 2010 Table 4 final model: residual variability CV 42.0%
  })

  model({
    # Power-form multiplicative diarrhea effect shared by CL/F and V/F (theta_Di^Diarrhea per Kohl 2010 Table 3 / 4).
    diarrhea_eff <- e_diarrhea_cl_vc^DIARRHEA

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * diarrhea_eff
    vc <- (exp(lvc) + e_age_vc * (AGE - 49)) * diarrhea_eff

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL = 1000 ug/L; Kohl 2010 reports concentrations
    # in ug/L (e.g. mean 411 ug/L) so multiply central / vc by 1000 to land in ug/L.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
