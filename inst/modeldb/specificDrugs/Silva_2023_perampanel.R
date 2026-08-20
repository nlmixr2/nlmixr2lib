Silva_2023_perampanel <- function() {
  description <- "One-compartment popPK model of perampanel in Portuguese adults with refractory epilepsy on chronic therapeutic drug monitoring (Silva 2023). NONMEM FOCE-I fit of 72 steady-state plasma concentrations from 44 patients; concentrations were obtained during the elimination phase only (9.7-24 h post-dose) so first-order absorption could not be estimated, and dosing was modelled as an IV bolus into the central compartment (perampanel oral bioavailability is ~100% per the summary of product characteristics). Concomitant enzyme-inducing antiepileptic drugs (carbamazepine, oxcarbazepine, phenobarbital, or phenytoin) increase apparent clearance 2.76-fold; body-mass index enters the volume of distribution as a power function centred at the sample median BMI of 25.1 kg/m^2. Inter-patient variability is exponential on clearance only; residual error is proportional."
  reference <- "Silva R, Colom H, Bicker J, Almeida A, Silva A, Sales F, Santana I, Falcao A, Fortuna A. Population Pharmacokinetic Analysis of Perampanel in Portuguese Patients Diagnosed with Refractory Epilepsy. Pharmaceutics. 2023;15(6):1704. doi:10.3390/pharmaceutics15061704"
  vignette <- "Silva_2023_perampanel"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central = list(analyte = "perampanel", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    BMI = list(
      description        = "Body mass index at baseline",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed. Enters the volume of distribution as (BMI/25.1)^e_bmi_vc where 25.1 kg/m^2 is the sample median reported by Silva 2023 (Table 1 median 25.14, rounded to 25.1 in the model equation on page 6 and Table S5 model 3). Studied range 15.76-36.20 kg/m^2.",
      source_name        = "BMI"
    ),
    CONMED_EIAED = list(
      description        = "Concomitant enzyme-inducing antiepileptic drug (any of carbamazepine, oxcarbazepine, phenobarbital, or phenytoin)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant EIAED)",
      notes              = "Silva 2023 Section 2.4.2 defines the EIAED bucket as carbamazepine, oxcarbazepine, phenobarbital, and phenytoin. 20 of 44 patients (45.5%) were co-prescribed at least one EIAED (Section 3, page 5). Multiplicative effect on apparent clearance (cl *= e_eiaed_cl^CONMED_EIAED) with e_eiaed_cl = 2.76 (Table 2 final model).",
      source_name        = "IND"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 44L,
    n_studies      = 1L,
    age_range      = "19-76 years",
    age_median     = "33 years",
    weight_range   = "45-99 kg",
    weight_median  = "77.5 kg",
    sex_female_pct = 38.6,
    disease_state  = "Adults with refractory epilepsy (two or more appropriate and tolerated antiepileptic drug schedules failed to achieve sustained seizure freedom); most patients on polytherapy with one or more concomitant AEDs alongside perampanel",
    dose_range     = "2-10 mg once-daily oral perampanel administered at bedtime (Silva 2023 Table 1; 12 mg extrapolated in model-based simulations)",
    regions        = "Portugal (Refractory Epilepsy Centre, Coimbra Hospital and University Centre, EPE)",
    notes          = "Retrospective single-centre observational study of TDM records collected April 2019-December 2022 (Silva 2023 Section 2.1). 72 steady-state plasma concentrations sampled 9.7-14 h post-dose (n=42) or 20.5-24 h post-dose (n=30). Concomitant AED landscape: levetiracetam 63.6%, carbamazepine 34.1%, clobazam 25.0%, eslicarbazepine acetate 22.7%, other AEDs at lower frequencies (Silva 2023 Table 1); EIAED co-prescription 45.5%. Median BMI 25.14 kg/m^2 (range 15.76-36.20)."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK parameters. Silva 2023 Table 2 final-model column.
    # A one-compartment model with first-order elimination was chosen;
    # dosing enters `central` directly (IV-bolus parameterisation)
    # because the elimination-phase sampling (9.7-24 h post-dose) does
    # not support estimation of a first-order absorption rate constant
    # (Silva 2023 Section 2.4.1). Perampanel oral bioavailability is
    # ~100% (Discussion page 3), so an IV bolus into central produces
    # the same steady-state trough profile as an oral dose with rapid
    # absorption and F = 1.
    # ------------------------------------------------------------------
    lcl <- log(0.419); label("Perampanel apparent clearance CL/F at no-EIAED reference (L/h)")             # Silva 2023 Table 2 final-model TVCL = 0.419
    lvc <- log(29.50); label("Perampanel apparent volume of distribution V/F at BMI 25.1 kg/m^2 (L)")     # Silva 2023 Table 2 final-model TVV = 29.50

    # ------------------------------------------------------------------
    # Covariate effects. Silva 2023 Table 2 final-model column and
    # Table S5 model 3 (the final covariate model). The EIAED effect
    # is a multiplier on typical CL when at least one EIAED is
    # co-administered; the BMI effect is a power function on typical
    # V/F centred at the sample median (25.1 kg/m^2) with a
    # positive-direction exponent of 2.12.
    # ------------------------------------------------------------------
    e_eiaed_cl <- 2.76; label("Multiplicative effect of concomitant EIAED on perampanel clearance (unitless)")  # Silva 2023 Table 2 final-model IND_CL = 2.76
    e_bmi_vc   <- 2.12; label("Power-function exponent of BMI (centred at 25.1 kg/m^2) on V/F (unitless)")     # Silva 2023 Table 2 final-model BMI_V = 2.12

    # ------------------------------------------------------------------
    # Inter-individual variability. Silva 2023 modelled IPV
    # exponentially on CL only (Section 3.1); IPV on V was not
    # supported by the data (Silva 2023 Section 2.4.1). Reported as
    # %CV = 30.82 in Table 2 (final-model IPV_CL); converted to the
    # log-normal variance via omega^2 = log(1 + CV^2) = log(1 + 0.3082^2)
    # = 0.09067.
    # ------------------------------------------------------------------
    etalcl ~ 0.09067  # CV = 30.82% -> log(1 + 0.3082^2); Silva 2023 Table 2 final-model IPV_CL

    # ------------------------------------------------------------------
    # Residual error. Silva 2023 modelled the residual error as
    # proportional (Section 3.1); reported as 6.44% in Table 2
    # final-model RE_proportional.
    # ------------------------------------------------------------------
    propSd <- 0.0644; label("Proportional residual error on perampanel plasma concentration (fraction)")  # Silva 2023 Table 2 final-model RE_proportional = 6.44%
  })

  model({
    # 1. Individual PK parameters
    #    CL_i = TVCL * e_eiaed_cl^CONMED_EIAED * exp(eta_lcl)
    #    V_i  = TVV  * (BMI / 25.1)^e_bmi_vc
    #    Silva 2023 Equations (1) and (2) / Table S5 model 3.
    cl <- exp(lcl + etalcl) * e_eiaed_cl^CONMED_EIAED
    vc <- exp(lvc) * (BMI / 25.1)^e_bmi_vc

    # 2. Micro-constant
    kel <- cl / vc

    # 3. One-compartment ODE, IV bolus into central
    d/dt(central) <- -kel * central

    # 4. Observation and error
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
