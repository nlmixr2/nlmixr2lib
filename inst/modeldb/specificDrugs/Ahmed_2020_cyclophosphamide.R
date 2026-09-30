Ahmed_2020_cyclophosphamide <- function() {
  description <- "One-compartment population PK model for intravenous cyclophosphamide in Ethiopian women with breast cancer (Ahmed 2020), with the cyclophosphamide dosage regimen (500 vs 600 mg/m^2 based) on clearance and volume and body surface area on volume, linked to a linear direct-response (empiric) model of absolute neutrophil count driven by the cumulative cyclophosphamide AUC."
  reference <- "Ahmed JH, Makonnen E, Bisaso RK, Mukonzo JK, Fotoohi A, Aseffa A, Howe R, Hassan M, Aklillu E. Population Pharmacokinetic, Pharmacogenetic, and Pharmacodynamic Analysis of Cyclophosphamide in Ethiopian Breast Cancer Patients. Front Pharmacol. 2020;11:406. doi:10.3389/fphar.2020.00406"
  vignette <- "Ahmed_2020_cyclophosphamide"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L",
    ANC = "cells/mm^3"
  )

  # `auc` is a bookkeeping integrator of the plasma concentration; it holds
  # no drug and exists only to supply the AUC regressor of the empiric
  # neutrophil model (Ahmed 2020 Equation 5), following the auc_<scope>
  # pattern of Bender_2024_mosunetuzumab.R.
  paper_specific_compartments <- c("auc")

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power model on the apparent volume only, normalised to 1.58 m^2: V = V_pop * (BSA/1.58)^0.861 (Ahmed 2020 Table 4). Methods Equation 1 describes the normalising value as the cohort median; the cohort mean was 1.59 +/- 0.20 m^2 (Table 1). Baseline (time-fixed) value; the per-patient dose was also calculated from it.",
      source_name = "BSA"
    ),
    DOSE_HIGH = list(
      description = "Cyclophosphamide 600 mg/m^2 based regimen indicator (1 = 600 mg/m^2 based AC or AC-T regimen, 0 = 500 mg/m^2 based FAC regimen)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 in the register convention; the source reference level is the 600 mg/m^2 regimen (DOSE_HIGH = 1)",
      notes = "Ahmed 2020 stratified patients into a 600 mg/m^2 cyclophosphamide based regimen (AC: doxorubicin 50 mg/m^2 + cyclophosphamide 600 mg/m^2; AC-T: doxorubicin 60 mg/m^2 + cyclophosphamide 600 mg/m^2 followed by paclitaxel; 161 patients, 60.3%) and a 500 mg/m^2 based regimen (FAC: 5-fluorouracil 500 mg/m^2 + doxorubicin 50 mg/m^2 + cyclophosphamide 500 mg/m^2; 106 patients, 39.7%) (Methods 'Patients'; Table 1). Table 4 takes the 600 mg/m^2 regimen as the reference and applies (1 + THETA) to the 500 mg/m^2 regimen, so the model applies the coefficient to the complement (1 - DOSE_HIGH), as in Jonsson_2016_tanezumab.R. Only two dose levels exist in the data, so the indicator names a level rather than binning a continuous dose. The regimens also differ in co-administered drugs (5-fluorouracil is given only in FAC), so the effect is a regimen effect and must not be read as a pure cyclophosphamide dose nonlinearity.",
      source_name = "CPA regimen (500 vs 600 mg/m^2 based)"
    )
  )

  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Tested in the stepwise covariate model but not retained (Ahmed 2020 Methods and Table 4). The reported BMI-group difference in volume (33.5 L underweight vs 48.1 L overweight vs 51.9 L obese) is a post hoc ANOVA on empirical Bayes estimates, not a model term."
    ),
    SNP_CYP3A5_RS776746 = list(
      description = "CYP3A5*3 genotype (rs776746)",
      units = "(count of variant alleles or carrier indicator)",
      type = "categorical",
      notes = "Tested in the stepwise covariate model but not retained; the reported lower elimination rate constant in CYP3A5*3/*6 carriers is a post hoc t-test on empirical Bayes estimates (Ahmed 2020 Results, Figures 4 and 5)."
    ),
    SNP_CYP2C9_RS1799853 = list(
      description = "CYP2C9*2 genotype (rs1799853), analysed together with CYP2C9*3 (rs1057910)",
      units = "(count of variant alleles or carrier indicator)",
      type = "categorical",
      notes = "Tested in the stepwise covariate model but not retained; the reported higher clearance in CYP2C9*2/*3 carriers is a post hoc comparison of empirical Bayes estimates (Ahmed 2020 Results, Figures 4 and 6). CYP2B6, CYP2C19, CYP2J2, POR and ABCB1 genotypes, cardiovascular comorbidity, HIV status, AST, ALT, ALP, serum creatinine and BUN were also screened and not retained."
    )
  )

  compartmentData <- list(
    central = list(analyte = "cyclophosphamide", units = "mg", specimen = "plasma", verified = TRUE),
    auc = list(analyte = "cyclophosphamide", units = "mg*h/L", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 267,
    n_studies = 1,
    age_range = "IQR 33-48 years",
    age_median = "38 years",
    bsa_mean = "1.59 +/- 0.20 m^2 (mean +/- SD)",
    bmi_mean = "23.78 +/- 4.8 kg/m^2 (mean +/- SD)",
    sex_female_pct = 100,
    race_ethnicity = c(Black = 100),
    disease_state = "breast cancer, first cycle of cyclophosphamide-containing chemotherapy (AC, AC-T or FAC)",
    dose_range = "cyclophosphamide 500 or 600 mg/m^2 as a 30-min IV infusion; absolute dose 600-1,150 mg (median 870 mg; 777.5 mg in the 500 mg/m^2 group and 930 mg in the 600 mg/m^2 group)",
    regions = "Ethiopia (Tikur Anbessa Specialized Hospital, Addis Ababa)",
    baseline_anc = "median 3645.5 cells/mm^3 (IQR 2645-4765)",
    n_observations = "532 plasma cyclophosphamide concentrations (average 2 per patient); ANC at baseline and day 20",
    notes = "Baseline characteristics from Ahmed 2020 Table 1; genotype frequencies in Table 2. Table 1 labels neutrophil counts '10^3 cells/mm^3' but the values (median 3645.5) and the inclusion criterion (ANC >= 1,500/mm^3) are in cells/mm^3."
  )

  ini({
    # ---- Pharmacokinetics (Ahmed 2020 Table 3, final model) ----
    # Typical values refer to a patient on the 600 mg/m^2 based regimen with
    # BSA = 1.58 m^2.
    lcl <- log(5.41); label("Clearance for the 600 mg/m^2 regimen (L/h)") # Ahmed 2020 Table 3 final model CL = 5.41 L/h (%SE 8.4)
    lvc <- log(46.5); label("Volume of distribution for the 600 mg/m^2 regimen at BSA 1.58 m^2 (L)") # Ahmed 2020 Table 3 final model VD = 46.5 L (%SE 13.2)

    # ---- Covariate effects (Ahmed 2020 Table 4) ----
    e_dose_high_cl <- -0.323; label("Fractional change in CL for the 500 mg/m^2 vs 600 mg/m^2 regimen (unitless)") # Ahmed 2020 Table 4 CL 'CPA dose 500 mg/m2 (1 + THETA1)' = -0.323
    e_dose_high_vc <- -0.371; label("Fractional change in V for the 500 mg/m^2 vs 600 mg/m^2 regimen (unitless)") # Ahmed 2020 Table 4 VD 'CPA dose 500 mg/m2 (1 + THETA1)' = -0.371
    e_bsa_vc <- 0.861; label("Power exponent for BSA on V (unitless)") # Ahmed 2020 Table 4 VD '(BSA/1.58)**THETA2' = 0.861

    # ---- Inter-individual variability (Ahmed 2020 Table 3, final model) ----
    # %CV = sqrt(exp(omega^2) - 1) (Methods), so omega^2 = log(CV^2 + 1).
    etalcl ~ 0.19500 # Ahmed 2020 Table 3 final model BSV CL = 46.4%; log(0.464^2 + 1) = 0.19500
    etalvc ~ 0.12123 # Ahmed 2020 Table 3 final model BSV VD = 35.9%; log(0.359^2 + 1) = 0.12123

    # ---- PK residual error (Ahmed 2020 Table 3) ----
    addSd <- 1.54; label("Additive residual error on cyclophosphamide concentration (mg/L)") # Ahmed 2020 Table 3 'Additive error 1' = 1.54 mg/L (%SE 24)

    # ---- Pharmacodynamics: linear direct response on ANC (Ahmed 2020 Results) ----
    lrbase_anc <- log(3450); label("Baseline absolute neutrophil count ANC0 (cells/mm^3)") # Ahmed 2020 Results 'Pharmacodynamic Modeling': ANCo = 3,450
    slope_anc <- -1.42; label("Linear change in ANC per unit cyclophosphamide AUC ((cells/mm^3)/(umol*h/L))") # Ahmed 2020 Results 'Pharmacodynamic Modeling': slope = -1.42

    # ---- PD residual error ----
    # Ahmed 2020 Results: a combined proportional and additive error model was
    # selected for ANC, but neither magnitude is reported, so both are held
    # at zero.
    propSd_ANC <- fixed(0); label("Proportional residual error on ANC (fraction); magnitude not reported") # Ahmed 2020 Results 'Pharmacodynamic Modeling' (combined error selected, values not reported)
    addSd_ANC <- fixed(0); label("Additive residual error on ANC (cells/mm^3); magnitude not reported") # Ahmed 2020 Results 'Pharmacodynamic Modeling' (combined error selected, values not reported)
  })

  model({
    # 1. Covariate multipliers. The source reference level is the 600 mg/m^2
    #    regimen (DOSE_HIGH = 1); the (1 + THETA) factor applies to the
    #    500 mg/m^2 regimen (Ahmed 2020 Table 4).
    dose_cl <- 1 + e_dose_high_cl * (1 - DOSE_HIGH)
    dose_vc <- 1 + e_dose_high_vc * (1 - DOSE_HIGH)

    # 2. Individual PK parameters (Ahmed 2020 Methods Equations 1 and 2;
    #    Table 4)
    cl <- exp(lcl + etalcl) * dose_cl
    vc <- exp(lvc + etalvc) * (BSA / 1.58)^e_bsa_vc * dose_vc

    # 3. Micro-constant
    kel <- cl / vc

    # 4. One-compartment disposition after a 30-min IV infusion into central
    d/dt(central) <- -kel * central
    Cc <- central / vc

    # Cumulative plasma AUC (mg*h/L) from the time of dosing
    d/dt(auc) <- Cc

    # 5. Empiric (linear direct response) neutrophil model, Ahmed 2020
    #    Equation 5. The AUC regressor is in umol*h/L, the unit Table 3 reports
    #    for AUC0-day20, obtained from mg*h/L with the cyclophosphamide
    #    molecular weight 261.09 g/mol. The slope is the reported -1.42; it is
    #    added (not subtracted, as Equation 5 is printed) because the Results
    #    state a one-unit AUC increase decreases the count by 1.42 and every
    #    Figure 7 prediction lies at or below ANC0.
    auc_umol <- auc * 1000 / 261.09
    rbase_anc <- exp(lrbase_anc)
    ANC <- rbase_anc + slope_anc * auc_umol

    Cc ~ add(addSd)
    ANC ~ add(addSd_ANC) + prop(propSd_ANC)
  })
}
