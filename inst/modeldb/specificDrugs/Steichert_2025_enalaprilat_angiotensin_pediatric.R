Steichert_2025_enalaprilat_angiotensin_pediatric <- function() {
  description <- "Population PD model of the plasma angiotensin II / angiotensin I ratio versus serum enalaprilat concentration in ACE-inhibitor-naive children with heart failure treated with enalapril orodispersible minitablets (Steichert 2025, Pharmaceutics; LENA studies). Direct full-inhibition Imax model without a sigmoidicity factor, driven by the enalaprilat concentration supplied as the time-varying covariate CP_ENALAPRILAT_NGML (the paper fitted the simultaneously measured serum concentrations; no PK model is part of this fit). IIV on IC50 only and proportional residual error. To simulate from a dose, generate CP_ENALAPRILAT_NGML from the paediatric enalapril/enalaprilat PK model of the same group (Steichert_2025_enalapril_enalaprilat_pediatric), as the paper's Conclusions suggest."
  reference <- "Steichert M, Cawello W, Burckhardt BB, Suessenbach FK, Laeer S; LENA Consortium. Angiotensin II/Angiotensin I Ratio as a New Pharmacodynamic Parameter for Population Modelling in Healthy Adults and Children with Heart Failure Treated with Enalapril. Pharmaceutics. 2025;17(10):1345. doi:10.3390/pharmaceutics17101345"
  vignette <- "Steichert_2025_enalaprilat_angiotensin"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  covariateData <- list(
    CP_ENALAPRILAT_NGML = list(
      description = "Serum enalaprilat concentration at the time of the angiotensin sample, the direct driver of the Imax inhibition of the angiotensin II / angiotensin I ratio.",
      units = "ng/mL (= ug/L)",
      type = "continuous",
      reference_category = NULL,
      notes = "Measured simultaneously with each angiotensin sample in the source data (Section 2.7.1); predose values below the LLOQ were set to 0 (Section 3.2). Observed range 0-18.3 ug/L. Set to 0 before the first dose.",
      source_name = "serum enalaprilat concentration"
    )
  )

  covariatesDataExcluded <- list()

  population <- list(
    species = "human",
    n_subjects = 27,
    n_studies = 3,
    age_range = "0.07-2.24 years",
    age_median = "0.36 years",
    weight_range = "3.2-13.0 kg",
    weight_median = "4.8 kg",
    sex_female_pct = 55.6,
    disease_state = "ACE-inhibitor-naive children with heart failure due to congenital heart disease (88.9%) or dilated cardiomyopathy (11.1%); modified Ross score 0-9 (median 4).",
    dose_range = "Enalapril maleate orodispersible minitablets (0.25 mg and 1 mg); median first dose 0.06 mg/kg (range 0.03-0.08 mg/kg); on repeated dosing mostly twice daily, median daily dose 0.11 mg/kg (range 0.06-0.27 mg/kg).",
    regions = "Austria, Germany, Hungary, the Netherlands, Serbia (LENA PK-bridging studies EudraCT 2015-002335-17 and 2015-002396-18 and follow-up study 2015-002397-21, 2016-2018).",
    n_observations = "54 paired serum enalaprilat concentrations and angiotensin II / angiotensin I ratios: 16 before the first dose, 12 after the first dose, 26 after repeated doses (3.9 days to 11.8 months). Angiotensin II below the LLOQ was replaced by LLOQ/2 in 55.6% of ratios.",
    notes = "Demographics from Table 1 (per-observation summaries) and Sections 2.2 and 3.2. Angiotensin I and II by multiplex LC-HRMS."
  )

  ini({
    le0 <- log(0.19); label("Baseline angiotensin II / angiotensin I ratio E0 (unitless)") # Table 3: E0 = 0.19 (RSE 8.2%)
    lic50 <- log(1.19); label("Enalaprilat concentration giving half-maximal inhibition IC50 (ug/L)") # Table 3: IC50 = 1.19 ug/L (RSE 17.9%)

    # Monolix log-normal IIV reported as CV%; omega^2 = log(1 + 0.5992^2).
    etalic50 ~ 0.30678 # Table 3: IIV IC50 = 59.92 CV%

    propSd <- 0.37; label("Proportional residual error on the angiotensin II / angiotensin I ratio (fraction)") # Table 3: proportional error = 0.37
  })

  model({
    e0 <- exp(le0)
    ic50 <- exp(lic50 + etalic50)

    # Full-inhibition Imax model without a sigmoidicity factor (Section 3.6.1),
    # Monolix library form E = E0 * (1 - C / (C + IC50)).
    ang2_ang1_ratio <- e0 * (1 - CP_ENALAPRILAT_NGML / (CP_ENALAPRILAT_NGML + ic50))

    ang2_ang1_ratio ~ prop(propSd)
  })
}
