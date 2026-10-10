Abouellil_2023_remdesivir <- function() {
  description <- paste(
    "Six-compartment population PK model for intravenous remdesivir and",
    "its two plasma metabolites GS-704277 and GS-441524 in healthy adults",
    "(Abouellil 2023). Each analyte has a central and a peripheral",
    "compartment. Remdesivir is converted to GS-704277 from both of its",
    "compartments (central to central, and peripheral to peripheral);",
    "GS-704277 is converted to GS-441524 from its central compartment only;",
    "each analyte is also eliminated from its central compartment. States",
    "are molar amounts (umol), so doses must be supplied in umol of",
    "remdesivir. Fitted in Monolix to digitised MEAN concentration profiles",
    "of the Gilead single-ascending-dose study, so the random effects are",
    "between-dose-cohort, not between-subject, variability."
  )
  reference <- "Abouellil A, Bilal M, Taubert M, Fuhr U. A population pharmacokinetic model of remdesivir and its major metabolites based on published mean values from healthy subjects. Naunyn Schmiedebergs Arch Pharmacol. 2023;396(1):73-82. doi:10.1007/s00210-022-02292-6"
  vignette <- "Abouellil_2023_remdesivir"

  # Abouellil 2023 Figures 2 and 4 plot every analyte in nmol/L, and the
  # parent-to-metabolite fluxes in Table 1 carry no molecular-weight
  # factor, so the model runs in molar space. Doses are umol of
  # remdesivir (umol = mg / 602.58 * 1000); umol / L * 1000 gives nmol/L.
  units <- list(time = "h", dosing = "umol", concentration = "nmol/L")

  compartmentData <- list(
    central = list(analyte = "remdesivir", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "remdesivir", units = "umol", specimen = "tissue", verified = TRUE),
    central_gs704277 = list(analyte = "GS-704277", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1_gs704277 = list(analyte = "GS-704277", units = "umol", specimen = "tissue", verified = TRUE),
    central_gs441524 = list(analyte = "GS-441524", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1_gs441524 = list(analyte = "GS-441524", units = "umol", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "18-55 years",
    bmi_range = "18-30 kg/m^2",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not tabulated; the Discussion describes the source cohort as having 'a focus on a Hispanic population'.",
    disease_state = "Healthy volunteers (male, and non-pregnant non-lactating female).",
    dose_range = "Single 2-h intravenous infusions of 3, 10, 30, 75, 150 and 225 mg remdesivir (the six cohorts fitted). The source phase I programme also gave once-daily 1-h infusions for 7 and 14 days, which were not part of the fit (Table 2 caption).",
    regions = "United States (Gilead Sciences phase I programme).",
    notes = "Abouellil 2023 had no individual data: mean concentration-time profiles of the Gilead phase I single-ascending-dose study (Humeniuk et al. 2020, Antimicrob Agents Chemother) were digitised with GetData Graph Digitizer and fitted in Monolix 2019R2. Each dose cohort contributes one mean profile per analyte, so the random effects describe between-cohort variability of the mean profile (Methods), not between-subject variability. The number of subjects per cohort is not reported in Abouellil 2023."
  )

  ini({
    # Structural parameters: Abouellil 2023 Table 2 (population estimates
    # of the fixed effects; SD of the random effects in parentheses).
    # Table 1 and parts of the text print the intermediate as 'GS-774277';
    # Table 2, Figure 1 and the Abstract print GS-704277, which is the
    # correct Gilead code.

    # Remdesivir
    lcl <- log(18.1)
    label("Remdesivir clearance from the central compartment (L/h)") # Table 2, Remdesivir column, 'Total body clearance' = 18.1 L/h
    lvc <- log(4.89)
    label("Remdesivir central volume of distribution (L)") # Table 2, Remdesivir column, 'Central compartment volume of distribution' = 4.89 L
    lvp <- log(46.5)
    label("Remdesivir peripheral volume of distribution (L)") # Table 2, Remdesivir column, 'Peripheral compartment volume of distribution' = 46.5 L
    lq <- log(13.2)
    label("Remdesivir intercompartmental clearance (L/h)") # Table 2, Remdesivir column, 'Inter-compartmental clearance' = 13.2 L/h

    # Remdesivir -> GS-704277 formation, from both remdesivir compartments
    lcl_form_gs704277_central <- log(16.9)
    label("GS-704277 formation clearance from central remdesivir (L/h)") # Table 2, GS-704277 column, 'Central formation clearance' = 16.9 L/h
    lcl_form_gs704277_peripheral1 <- log(18.9)
    label("GS-704277 formation clearance from peripheral remdesivir (L/h)") # Table 2, GS-704277 column, 'Peripheral formation clearance' = 18.9 L/h

    # GS-704277
    lcl_gs704277 <- log(36.9)
    label("GS-704277 clearance from its central compartment (L/h)") # Table 2, GS-704277 column, 'Total body clearance' = 36.9 L/h
    lvc_gs704277 <- log(96.4)
    label("GS-704277 central volume of distribution (L)") # Table 2, GS-704277 column, 'Central compartment volume of distribution' = 96.4 L
    lvp_gs704277 <- log(8.64)
    label("GS-704277 peripheral volume of distribution (L)") # Table 2, GS-704277 column, 'Peripheral compartment volume of distribution' = 8.64 L
    lq_gs704277 <- log(0.12)
    label("GS-704277 intercompartmental clearance (L/h)") # Table 2, GS-704277 column, 'Inter-compartmental clearance' = 0.12 L/h

    # GS-704277 -> GS-441524 formation, central compartment only
    lcl_form_gs441524 <- log(50.5)
    label("GS-441524 formation clearance from central GS-704277 (L/h)") # Table 2, GS-441524 column, 'Central formation clearance' = 50.5 L/h

    # GS-441524
    lcl_gs441524 <- log(4.74)
    label("GS-441524 clearance from its central compartment (L/h)") # Table 2, GS-441524 column, 'Total body clearance' = 4.74 L/h
    lvc_gs441524 <- log(26.2)
    label("GS-441524 central volume of distribution (L)") # Table 2, GS-441524 column, 'Central compartment volume of distribution' = 26.2 L
    lvp_gs441524 <- log(66.2)
    label("GS-441524 peripheral volume of distribution (L)") # Table 2, GS-441524 column, 'Peripheral compartment volume of distribution' = 66.2 L
    lq_gs441524 <- log(55)
    label("GS-441524 intercompartmental clearance (L/h)") # Table 2, GS-441524 column, 'Inter-compartmental clearance' = 55 L/h

    # Between-dose-cohort variability. Table 2 prints the Monolix SD of
    # each log-normal random effect (caption: 'SD of the random effects');
    # the variances below are SD^2. The seven parameters carrying a random
    # effect are those listed in Results ('The final model had exponential
    # inter-cohort variabilities for ...'); all other parameters carry none.
    etalcl ~ 0.1521 # Table 2, Remdesivir 'Total body clearance' SD 0.39; 0.39^2
    etalcl_form_gs704277_central ~ 0.0625 # Table 2, GS-704277 'Central formation clearance' SD 0.25; 0.25^2
    etalcl_form_gs704277_peripheral1 ~ 0.2809 # Table 2, GS-704277 'Peripheral formation clearance' SD 0.53; 0.53^2
    etalcl_gs704277 ~ 0.0961 # Table 2, GS-704277 'Total body clearance' SD 0.31; 0.31^2
    etalcl_form_gs441524 ~ 0.0729 # Table 2, GS-441524 'Central formation clearance' SD 0.27; 0.27^2
    etalvc_gs441524 ~ 0.5041 # Table 2, GS-441524 'Central compartment volume of distribution' SD 0.71; 0.71^2
    etalvp_gs441524 ~ 0.0576 # Table 2, GS-441524 'Peripheral compartment volume of distribution' SD 0.24; 0.24^2

    # Residual error. Results: 'the error models which matched the data
    # best were proportional for remdesivir, combined for GS-704277, and
    # proportional for GS-441524'. No residual-error magnitude is reported
    # in the paper or its supplement, so each term is held at zero to keep
    # the selected structure visible; simulations are free of residual
    # noise.
    propSd <- fixed(0)
    label("Proportional residual SD, remdesivir (fraction); magnitude not reported") # Results: proportional error model selected, value not reported
    addSd_gs704277 <- fixed(0)
    label("Additive residual SD, GS-704277 (nmol/L); magnitude not reported") # Results: combined error model selected, value not reported
    propSd_gs704277 <- fixed(0)
    label("Proportional residual SD, GS-704277 (fraction); magnitude not reported") # Results: combined error model selected, value not reported
    propSd_gs441524 <- fixed(0)
    label("Proportional residual SD, GS-441524 (fraction); magnitude not reported") # Results: proportional error model selected, value not reported
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)

    cl_form_gs704277_central <- exp(lcl_form_gs704277_central + etalcl_form_gs704277_central)
    cl_form_gs704277_peripheral1 <- exp(lcl_form_gs704277_peripheral1 + etalcl_form_gs704277_peripheral1)

    cl_gs704277 <- exp(lcl_gs704277 + etalcl_gs704277)
    vc_gs704277 <- exp(lvc_gs704277)
    vp_gs704277 <- exp(lvp_gs704277)
    q_gs704277 <- exp(lq_gs704277)

    cl_form_gs441524 <- exp(lcl_form_gs441524 + etalcl_form_gs441524)

    cl_gs441524 <- exp(lcl_gs441524)
    vc_gs441524 <- exp(lvc_gs441524 + etalvc_gs441524)
    vp_gs441524 <- exp(lvp_gs441524 + etalvp_gs441524)
    q_gs441524 <- exp(lq_gs441524)

    # Concentrations (umol/L) as defined at the top of Abouellil 2023 Table 1
    c_central <- central / vc
    c_peripheral1 <- peripheral1 / vp
    c_central_gs704277 <- central_gs704277 / vc_gs704277
    c_peripheral1_gs704277 <- peripheral1_gs704277 / vp_gs704277
    c_central_gs441524 <- central_gs441524 / vc_gs441524
    c_peripheral1_gs441524 <- peripheral1_gs441524 / vp_gs441524

    # Abouellil 2023 Table 1, one equation per state. Flows are
    # clearance x concentration (umol/h). The 1:1 molar transfer between
    # moieties means each metabolite clearance and volume is apparent with
    # respect to the (unknown) fraction of the upstream flux that appears
    # in plasma as that metabolite.
    d/dt(central) <- q * c_peripheral1 - q * c_central - cl * c_central -
      cl_form_gs704277_central * c_central
    d/dt(peripheral1) <- q * c_central - q * c_peripheral1 -
      cl_form_gs704277_peripheral1 * c_peripheral1
    d/dt(central_gs704277) <- cl_form_gs704277_central * c_central +
      q_gs704277 * c_peripheral1_gs704277 - q_gs704277 * c_central_gs704277 -
      cl_gs704277 * c_central_gs704277 - cl_form_gs441524 * c_central_gs704277
    d/dt(peripheral1_gs704277) <- q_gs704277 * c_central_gs704277 +
      cl_form_gs704277_peripheral1 * c_peripheral1 -
      q_gs704277 * c_peripheral1_gs704277
    d/dt(central_gs441524) <- cl_form_gs441524 * c_central_gs704277 +
      q_gs441524 * c_peripheral1_gs441524 - cl_gs441524 * c_central_gs441524 -
      q_gs441524 * c_central_gs441524
    d/dt(peripheral1_gs441524) <- q_gs441524 * c_central_gs441524 -
      q_gs441524 * c_peripheral1_gs441524

    # Plasma concentrations in nmol/L, the scale of Figures 2 and 4
    Cc <- 1000 * c_central
    Cc_gs704277 <- 1000 * c_central_gs704277
    Cc_gs441524 <- 1000 * c_central_gs441524

    Cc ~ prop(propSd)
    Cc_gs704277 ~ add(addSd_gs704277) + prop(propSd_gs704277)
    Cc_gs441524 ~ prop(propSd_gs441524)
  })
}
