Ruppel_2017_insulin_doseGIRAUC_TIRHI <- function() {
  description <- "Population dose-response (dose-GIRAUC) model for inhaled Technosphere insulin (TI) and subcutaneous regular human insulin (RHI) fitted jointly with a common Hill coefficient, in healthy volunteers (Ruppel 2017, common TI/RHI fit). Algebraic sigmoid Emax model, without ODEs, of the 0-20 h area under the insulin-induced glucose infusion rate (GIR above baseline) versus dose, with a formulation-specific ED50; the ED50 ratio (about 5.1) is the paper's TI:RHI equivalent-dose conversion factor. Fitted to GIR profiles simulated from the companion PK-GIR model (Ruppel_2017_insulin). The maximum GIRAUC is each subject's GIRmax times 20 h, carried as a fixed typical value with the PK-GIR between-subject variability of GIRmax. Inter-individual variability on both ED50s and the Hill coefficient; proportional residual error."
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
    dosing = "U for Technosphere insulin or IU for regular human insulin (carried in the DOSE covariate; no dose events)",
    concentration = "n/a (output GIRAUC in g/kg: grams of glucose per kg body weight)",
    GIRAUC = "g/kg"
  )

  covariateData <- list(
    DOSE = list(
      description = "Single insulin dose: inhaled Technosphere insulin in labelled TI units (U) or subcutaneous regular human insulin in international units (IU)",
      units = "U (TI) or IU (RHI)",
      type = "continuous",
      reference_category = NULL,
      notes = "The model has no time dimension: each record is one dose level, and GIRAUC is the 0-20 h integral of the insulin-induced GIR for that dose (Equation 3). Simulated doses in the source ranged up to 240 U TI and 80 IU RHI.",
      source_name = "DOSE"
    ),
    FORM_INSULIN_TI = list(
      description = "Insulin formulation indicator: 1 = inhaled Technosphere insulin (TI), 0 = subcutaneous regular human insulin (RHI)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Selects the formulation-specific ED50 and its random effect, reproducing the dose-GIRAUC control stream's IF(TYPE==1) (TI) / IF(TYPE==2) (RHI) branches.",
      source_name = "TYPE (1 = TI, 2 = RHI)"
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
    dose_range = "Simulated single doses up to 240 U TI and up to 80 IU s.c. RHI (clinical study tested 10-80 U TI and 15 IU RHI)",
    regions = "United States (Profil Institute for Clinical Research, Chula Vista, CA)",
    notes = "Fitted to GIRAUC values computed from GIR(t) profiles simulated with the individual post-hoc PK-GIR parameters of the 31 analysed subjects of study NCT01490762 (Results, 'Dose-response model'); not fitted to observed data."
  )

  ini({
    lemax <- fixed(log(14.287)); label("Typical maximum 0-20 h GIRAUC, GIRmax x 1200 min (g/kg)") # Equation 4 and Methods: GIRAUC,max = 11.9 mg/kg/min * 1200 min = 14.3 g/kg; 11.9059 * 1.2 from the PK-GIR stream's fixed GIRmax
    led50_lung <- log(245); label("Inhaled TI dose producing half of the maximum GIRAUC, ED50,TI (U)") # Table 2, TI/RHI column: ED50,TI = 245 (201-290) U
    led50_sc <- log(48); label("s.c. RHI dose producing half of the maximum GIRAUC, ED50,RHI (IU)") # Table 2, TI/RHI column: ED50,RHI = 48 (39-57) IU
    lhill <- log(1.08); label("Common Hill coefficient of the dose-GIRAUC model, gamma (unitless)") # Table 2, TI/RHI column: gamma = 1.08 (0.99-1.16)

    # Table 2 omegas printed as 100 * sqrt(omega^2), the convention of the
    # companion PK and PK-GIR tables.
    etalemax ~ fixed(0.0961) # GIRAUC,max is an individual covariate (GIRmax,i * 20 h; Methods); its spread is the PK-GIR Table 1 omega Emax = 31 %
    etaled50_lung ~ 0.16 # Table 2, TI/RHI column: omega ED50,TI = 40 %
    etaled50_sc ~ 0.2025 # Table 2, TI/RHI column: omega ED50,RHI = 45 %
    etalhill ~ 0.0289 # Table 2, TI/RHI column: omega gamma = 17 %

    propSd <- 0.21; label("Proportional residual error on GIRAUC (fraction)") # Table 2, TI/RHI column: sigma p = 21 %
  })
  model({
    emax <- exp(lemax + etalemax)
    ed50 <- FORM_INSULIN_TI * exp(led50_lung + etaled50_lung) +
      (1 - FORM_INSULIN_TI) * exp(led50_sc + etaled50_sc)
    hill <- exp(lhill + etalhill)

    GIRAUC <- emax * DOSE^hill / (ed50^hill + DOSE^hill)
    GIRAUC ~ prop(propSd)
  })
}
