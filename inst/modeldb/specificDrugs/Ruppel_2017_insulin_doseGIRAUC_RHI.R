Ruppel_2017_insulin_doseGIRAUC_RHI <- function() {
  description <- "Population dose-response (dose-GIRAUC) model for subcutaneous regular human insulin (RHI) in healthy volunteers (Ruppel 2017, RHI-only fit). Algebraic sigmoid Emax model, without ODEs, of the 0-20 h area under the insulin-induced glucose infusion rate (GIR above baseline) versus RHI dose, fitted to GIR profiles simulated from the companion PK-GIR model (Ruppel_2017_insulin) for doses up to 80 IU. The maximum GIRAUC is not estimated: it is each subject's GIRmax times the 20 h integration window, carried here as a fixed typical value with the PK-GIR between-subject variability of GIRmax. Inter-individual variability on ED50 and the Hill coefficient; proportional residual error."
  reference <- paste(
    "Ruppel D, Dahmen R, Boss A, Jager R, Grant M, Baughman R, Klabunde T.",
    "A Population Dose-Response Model for Inhaled Technosphere Insulin",
    "Administered to Healthy Subjects.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(6):365-372.",
    "doi:10.1002/psp4.12189.",
    sep = " "
  )
  vignette <- "Ruppel_2017_insulin"
  units <- list(
    time = "n/a (static dose-response model; no time dimension, GIRAUC is a 0-20 h integral)",
    dosing = "IU (carried in the DOSE covariate; no dose events)",
    concentration = "n/a (output GIRAUC in g/kg: grams of glucose per kg body weight)",
    GIRAUC = "g/kg"
  )

  covariateData <- list(
    DOSE = list(
      description = "Single subcutaneous regular human insulin dose",
      units = "IU",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose in international units (IU). The model has no time dimension: each record is one dose level, and GIRAUC is the 0-20 h integral of the insulin-induced GIR for that dose (Equation 3). Simulated doses in the source ranged up to 80 IU.",
      source_name = "DOSE"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 31L,
    n_studies = 1L,
    age_range = "18-53 years",
    weight_range = "49-100 kg",
    sex_female_pct = 13 / 31 * 100,
    disease_state = "Healthy volunteers",
    dose_range = "Simulated single s.c. RHI doses up to 80 IU (clinical study tested 15 IU)",
    regions = "United States (Profil Institute for Clinical Research, Chula Vista, CA)",
    notes = "Fitted to GIRAUC values computed from GIR(t) profiles simulated with the individual post-hoc PK-GIR parameters of the 31 analysed subjects of study NCT01490762 (Results, 'Dose-response model'); not fitted to observed data."
  )

  ini({
    lemax <- fixed(log(14.287)); label("Typical maximum 0-20 h GIRAUC, GIRmax x 1200 min (g/kg)") # Equation 4 and Methods: GIRAUC,max = 11.9 mg/kg/min * 1200 min = 14.3 g/kg; 11.9059 * 1.2 from the PK-GIR stream's fixed GIRmax
    led50 <- log(43); label("RHI dose producing half of the maximum GIRAUC, ED50,RHI (IU)") # Table 2, RHI column: ED50,RHI = 43 (36-51) IU
    lhill <- log(1.15); label("Hill coefficient of the dose-GIRAUC model, gamma (unitless)") # Table 2, RHI column: gamma = 1.15 (1.08-1.23)

    # Table 2 omegas printed as 100 * sqrt(omega^2), the convention of the
    # companion PK and PK-GIR tables.
    etalemax ~ fixed(0.0961) # GIRAUC,max is an individual covariate (GIRmax,i * 20 h; Methods); its spread is the PK-GIR Table 1 omega Emax = 31 %
    etaled50 ~ 0.2401 # Table 2, RHI column: omega ED50,RHI = 49 %
    etalhill ~ 0.0289 # Table 2, RHI column: omega gamma = 17 %

    propSd <- 0.04; label("Proportional residual error on GIRAUC (fraction)") # Table 2, RHI column: sigma p = 4 %
  })
  model({
    emax <- exp(lemax + etalemax)
    ed50 <- exp(led50 + etaled50)
    hill <- exp(lhill + etalhill)

    GIRAUC <- emax * DOSE^hill / (ed50^hill + DOSE^hill)
    GIRAUC ~ prop(propSd)
  })
}
