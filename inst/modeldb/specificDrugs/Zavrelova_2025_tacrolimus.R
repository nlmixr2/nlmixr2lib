Zavrelova_2025_tacrolimus <- function() {
  description <- "One-compartment population PK model for continuous-infusion intravenous tacrolimus (whole blood) in adults after allogeneic haematopoietic stem cell transplantation (Zavrelova 2025); clearance declines exponentially with age and is 23% lower in a repeat transplantation, and volume of distribution carries no covariate"
  reference <- "Zavrelova A, Zibridova K, Radocha J, Cermakova E, Rozsivalova P, Zak P, Visek B, Lanska M, Stevkova J, Merdita S, Slanar O, Sima M. Tacrolimus exposure during the three-month period following allogeneic stem cell transplantation predicts overall survival. Front Pharmacol. 2025;16:1517083. doi:10.3389/fphar.2025.1517083"
  vignette <- "Zavrelova_2025_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Tacrolimus was measured in EDTA whole blood by CMIA on an Abbott ARCHITECT
  # i1000SR (Methods 2.1).
  compartmentData <- list(
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Entered UNCENTRED on the log scale, exactly as printed in the source: CL = 7.94 * exp(-0.0085 * age) (Abstract; Results 3.5). CL_pop is therefore the extrapolated age-0 intercept, not the value at a reference age. Results 3.5 quotes the worked example of a 57-year-old (stated as the population median) first-transplant patient, CL = 4.89 L/h; Table 1 gives median 58 years (range 20-72) for the 86 first-transplant patients.",
      source_name = "age"
    ),
    HCT_PRIOR = list(
      description = "Repeat (second or later) allogeneic stem cell transplantation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (first allogeneic transplantation)",
      notes = "1 when the tacrolimus course belongs to a repeat allogeneic stem cell transplantation, i.e. the patient had at least one prior allogeneic HCT before the index transplant; 0 for a first transplantation. 14 of the 100 transplantations in the PK dataset were repeats (Results 3.5). Gates the exp(-0.26) = 0.77 fold change in CL (Table 4; Results 3.5, 'reduced by 23%').",
      source_name = "repeated Tx"
    )
  )

  # Tested as covariates (Methods 2.3) but not retained in the final model --
  # 'tacrolimus Vd did not correlate with any of the tested variables' and only
  # age and repeated transplantation entered CL (Results 3.5). The categorical
  # covariates diagnosis, remission before transplantation, CMV reactivation,
  # preparative (conditioning) regimen and sequential chemotherapy were also
  # screened; they have no canonical column and are listed in population$notes.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Median 82 kg (range 49-133; Table 1, first-transplant patients). Screened, not retained."
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      notes = "Median 172.5 cm (range 152-200; Table 1). Screened, not retained."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Median 26.0 kg/m^2 (range 16.8-42.7; Table 1). Screened, not retained."
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Median 1.96 m^2 (range 1.49-2.49; Table 1). Screened, not retained."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Tested as a continuous covariate (Methods 2.3); not summarised in Table 1. Not retained."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "36 of 86 first-transplant patients female (41.9%; Table 1). Tested as a categorical covariate, not retained."
    ),
    CONMED_AZOLE = list(
      description = "Concomitant azole antifungal",
      units = "(binary)",
      type = "binary",
      notes = "Co-medication with azole antifungals was tested in a dose-independent manner (Methods 2.3) and not retained. Every patient received an azole (voriconazole 28, posaconazole 20, fluconazole 52; Results 3.1), so there was no azole-free reference group; the Discussion reports median individual CL of 4.50, 4.77 and 5.12 L/h in the voriconazole, posaconazole and fluconazole subgroups."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 100,
    n_studies = 1,
    n_observations = 349,
    age_range = "20-72 years (first-transplant patients, Table 1)",
    age_median = "58 years (Table 1); 57 years is quoted as the median in the Results 3.5 worked example",
    weight_range = "49-133 kg",
    weight_median = "82 kg",
    sex_female_pct = 41.9,
    race_ethnicity = "White (Discussion: 'there were only white patients in our study population').",
    disease_state = "Adults (age >= 18 years) undergoing allogeneic peripheral blood stem cell transplantation, almost all for haematological malignancy: AML 50.0%, MDS 14.0%, ALL 9.3%, CLL + lymphoma 12.8%, PCL + MM 3.5%, myeloproliferative disease, histiocytosis and HLH (Table 1). Disease risk high in 59.3%; myeloablative conditioning in 59.3%.",
    dose_range = "Tacrolimus started at 0.03 mg/kg as a continuous intravenous infusion on the day before transplantation and switched to oral administration around day 11; doses adjusted to a 5-15 ng/mL trough target and when voriconazole or posaconazole was co-administered (Methods 2.1).",
    regions = "Czech Republic (University Hospital Hradec Kralove).",
    co_medication = "Mycophenolate mofetil, rabbit anti-thymocyte globulin and an azole antifungal (voriconazole, posaconazole or fluconazole) in every patient; no glucocorticoids.",
    notes = "Retrospective single-centre TDM analysis of transplantations performed 2016-2018. The PK model was fitted only to concentrations drawn during the initial continuous intravenous infusion phase (Methods 2.3), 349 whole-blood concentrations from 100 transplantations (86 first + 14 repeat; 2-7 per case). Demographics in Table 1 describe the 86 first-transplant patients used for the survival analysis. Categorical covariates screened without a canonical column: diagnosis, remission before transplantation, CMV reactivation, preparative regimen and sequential chemotherapy (Methods 2.3). The paper's Cox survival analysis (Tables 2-3) is a statistical regression on dichotomised exposure and is not part of this model."
  )

  ini({
    # Structural typical values. Age enters UNCENTRED (see covariateData$AGE),
    # so lcl is the extrapolated age-0 intercept of CL = 7.94 * exp(-0.0085 * age).
    lvc <- log(127.1); label("Volume of distribution (L)")                                   # Table 4 (Vd_pop = 127.1 L, RSE 7.9%)
    lcl <- log(7.94);  label("Clearance, age-0 intercept of the uncentred age model (L/h)")  # Table 4 (CL_pop = 7.94 L/h, RSE 15.3%)

    # Covariate effects on the natural-log scale of CL.
    e_age_cl       <- -0.0085; label("Effect of age on log CL (1/year)")                          # Table 4 (beta_CL_age = -0.0085, RSE 32.2%)
    e_hctprior_cl  <- -0.26;   label("Effect of repeat transplantation on log CL (unitless)")     # Table 4 (beta_CL_repeated Tx = -0.26, RSE 37.2%); exp(-0.26) = 0.77

    # IIV. Table 4 reports Omega as the standard deviation of the random effect
    # on the log scale, so the variance is Omega^2.
    etalvc ~ 0.45^2  # Table 4 (Omega_Vd = 0.45, RSE 18.6%)
    etalcl ~ 0.29^2  # Table 4 (Omega_CL = 0.29, RSE 8.6%)

    # Combined residual error (Results 3.5). Table 4 prints the constant term
    # without units as 0.0018; read as mg/L (= 1.8 ng/mL) because dose (mg) and
    # Vd (L) put model predictions in mg/L. See the vignette Assumptions section.
    addSd  <- 1.8;  label("Additive residual error (ng/mL)")          # Table 4 (Constant = 0.0018 mg/L, RSE 15.6%) x 1000
    propSd <- 0.17; label("Proportional residual error (fraction)")   # Table 4 (Proportional = 0.17, RSE 11.2%)
  })
  model({
    # Individual parameters. Source equations (Results 3.5):
    #   Vd = Vd_pop
    #   CL = CL_pop * exp(beta_CL_age * age) * exp(beta_CL_repeated Tx [if yes])
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + e_age_cl * AGE + e_hctprior_cl * HCT_PRIOR + etalcl)

    kel <- cl / vc

    # One compartment with linear elimination; the model was fitted to the
    # continuous intravenous infusion phase only, so there is no absorption.
    d/dt(central) <- -kel * central

    # Dose in mg, volume in L -> mg/L; x 1000 -> ng/mL (= ug/L)
    Cc <- 1000 * central / vc
    # Monolix 'combined' error, whose default form (combined1) sums the SDs:
    # y = f + (a + b * f) * e.
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
