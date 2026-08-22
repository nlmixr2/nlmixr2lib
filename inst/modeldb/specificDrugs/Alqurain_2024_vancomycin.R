Alqurain_2024_vancomycin <- function() {
  description <- "One-compartment IV population PK model for vancomycin in non-critical-care adult (>= 40 years) medical-ward inpatients in the Eastern Province of Saudi Arabia, built from routine therapeutic-drug-monitoring trough concentrations. Clearance is an uncentered exponential (log-linear) function of C-reactive protein and Cockcroft-Gault creatinine clearance, so exp(lcl) is the clearance at CRP = 0 and CRCL = 0 rather than a typical-patient clearance; volume of distribution carries no covariates."
  reference   <- "Alqurain AA, Alrashidi LN, Aloraifej SK, Alkhalifah M, Alsayed HA, Abohelaika S, Alshabeeb MA, Aldhafeeri AS, Almuslim M, Bumozah TN, Alomar MJ, Alshehab AA, Alamer AA, Al-Matouq J, Bidasee KR, Alomar FA. Factors Affecting Vancomycin Trough Concentration; a Population Pharmacokinetic Model in Non-Critical Care Saudi Patients. Drug Des Devel Ther. 2024;18:6185-6198. doi:10.2147/DDDT.S496512"
  vignette    <- "Alqurain_2024_vancomycin"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the source: the single disposition
  # compartment holds vancomycin and is sampled as plasma trough
  # concentration by the EMIT immunoassay (Methods, "Vancomycin Assay").
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRP = list(
      description        = "C-reactive protein concentration measured at the time the vancomycin trough sample was drawn",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Standard (not high-sensitivity) assay; NOT a baseline value -- Alqurain 2024 Discussion states 'the CRP values were",
        "measured at the time of vancomycin trough level collection, and they are not reflective to the actual inflammatory",
        "state when initiating vancomycin therapy.' Cohort mean 80.6 +/- 32.5 (Table 1); the Simulx simulation scenarios use",
        "20 (low), 80 (normal) and 180 (high) (Methods, 'Simulation of Vancomycin Exposure'; Figure 6 panel labels).",
        "UNITS CAVEAT: the source labels this column 'mg/dL' in Table 1 and in the Figure 6 caption, but the numeric range",
        "(10-180, mean 80.6) is two orders of magnitude above any plausible mg/dL CRP and is exactly the range routinely",
        "reported in mg/L for a serious-infection cohort. The canonical CRP register unit mg/L is therefore used here and the",
        "source numbers are carried across unchanged, so e_crp_cl applies per mg/L. Enters CL uncentered and untransformed:",
        "CL = exp(lcl) * exp(e_crp_cl * CRP + e_crcl_cl * CRCL). The uncentered form is not an inference -- it is the only",
        "form that reproduces Supplementary Table 9 (median simulated CL falls from 3 L/h at CRP 20 to 0.001 L/h at CRP 180,",
        "a ~6000-fold span that a log-transformed or mean-centered covariate cannot generate) and Figure 6.",
        sep = " "
      ),
      source_name        = "CRP"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Alqurain 2024 Methods, 'Study Design, Subjects, and Sample Collection': 'CrCl was calculated using the",
        "Cockcroft-Gault formula.' The paper applies no BSA normalization, so the column is raw mL/min; the canonical CRCL",
        "register entry explicitly accepts raw Cockcroft-Gault mL/min with the assay form documented per model (precedents:",
        "Georges_2009_ceftazidime.R, Chen_2023_nemonoxacin.R, Delattre_2010_amikacin.R). Cohort mean 61.1 +/- 48.2 mL/min",
        "with 35% of patients at <= 30 mL/min (Table 1); the Simulx simulation scenarios use 90 (normal), 60 (moderately",
        "low) and 30 (severely low) mL/min. Enters CL uncentered and untransformed alongside CRP. Patients with chronic",
        "kidney disease were an exclusion criterion, so the low-CrCl tail is acute rather than chronic renal impairment.",
        sep = " "
      ),
      source_name        = "CrCl"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 124L,
    n_observations = 172L,
    n_studies      = 1L,
    age_range      = "40-99 years (inclusion criterion >= 40 years); median 79 years (IQR 50-86). Age bands (Table 1): 40-49 8 (7%), 50-59 8 (7%), 60-69 42 (33%), 70-79 31 (25%), 80-89 26 (21%), 90-99 9 (7%).",
    age_median     = "79 years",
    weight_range   = "Mean 72 +/- 20 kg (Table 1)",
    sex_female_pct = 41,
    race_ethnicity = "Not reported; single-country cohort recruited in the Eastern Province of Saudi Arabia",
    disease_state  = "Adults admitted to general medical wards and initiated on systemic vancomycin. Patients in the emergency department, intensive care unit or surgical wards, patients with chronic kidney disease, and pregnant women were excluded, so this is explicitly a non-critical-care population. Most patients nonetheless had some degree of (acute) renal impairment: mean CrCl 61.1 +/- 48.2 mL/min, 35% at <= 30 mL/min. Mean serum albumin 29.5 +/- 8 g/L; mean serum creatinine 183 +/- 171 (units reported as mg/dL in Table 1, but the magnitude is consistent with umol/L).",
    renal_function = "Mean Cockcroft-Gault CrCl 61.1 +/- 48.2 mL/min; > 120 mL/min 17 (14%), 91-120 21 (17%), 61-90 14 (11%), 31-60 28 (23%), <= 30 mL/min 44 (35%) (Table 1)",
    dose_range     = "500-1750 mg per intravenous dose (Table 1: 500 mg 25%, 750 mg 15%, 850 mg 2%, 1000 mg 40%, 1200 mg 1%, 1250 mg 9%, 1500 mg 6%, 1750 mg 2%). Doses were standardized to 500 mg for figure presentation only; actual doses were used in the modelling.",
    regions        = "Saudi Arabia (Eastern Province): Al-Mana General Hospital Al-Khobar, Qatif Central Hospital, and Dammam Medical Complex",
    notes          = "Multicentre retrospective study of routine therapeutic-drug-monitoring records collected 1 January - 31 December 2022; one visit per patient. 172 trough concentrations from 124 patients. Baseline demographics are Table 1 of Alqurain 2024. Concentrations were measured by enzyme-multiplied immunoassay (EMIT) with a 2-50 mg/L calibration range, so only trough samples inform the fit -- the authors note this prevented resolution of a distribution phase and is why a one-compartment model was retained."
  )

  ini({
    # All estimates are the "Final Model / Population" column of Alqurain 2024
    # Table 2. Note that the Results prose ("The estimated V was 224.4 L and
    # the estimated CL was 3.52 L/h") mixes the BASE-model V (224.37 L) with
    # the FINAL-model CL; Table 2's final column and Supplementary Table 9
    # (median simulated V 192-198 L across all nine scenario groups) both give
    # the final V as 193.65 L, which is the value used here.
    #
    # The model was fitted in Monolix 2023 with log-normally distributed
    # individual parameters (Methods equ. 1), so the parameters below are the
    # log of the reported population values.

    lvc <- log(193.65); label("Volume of distribution (L)")  # Alqurain 2024 Table 2, Final Model column, "V (L)" = 193.65 (RSE 7.41%)

    # Uncentered exponential covariate model on CL: exp(lcl) is CL at
    # CRP = 0 mg/L AND CRCL = 0 mL/min, not the clearance of a typical
    # patient. See the model() block for the algebraic form and the
    # covariateData notes for the evidence that the form is uncentered.
    lcl <- log(3.52);   label("Clearance at CRP = 0 mg/L and CRCL = 0 mL/min (L/h)")  # Alqurain 2024 Table 2, Final Model column, "CL (L/h)" = 3.52 (RSE 62.5%)

    # Covariate effects (Alqurain 2024 Table 2, Final Model column). Both are
    # slopes on the natural-log scale of CL, per raw unit of the covariate.
    e_crp_cl  <- -0.05;  label("Exponential effect of CRP on CL (per mg/L)")        # Alqurain 2024 Table 2, Final Model column, "CRP effect on CL" = -0.05 (RSE 28%)
    e_crcl_cl <- 0.0088; label("Exponential effect of CRCL on CL (per mL/min)")     # Alqurain 2024 Table 2, Final Model column, "CrCl effect on CL" = 0.0088 (RSE 48.4%)

    # Inter-individual variability. Table 2's footnote states "Omega,
    # Inter-individual variability presented as standard deviation", and
    # Monolix reports omega as an SD, so the reported 0.33 / 0.81 are SDs of
    # the log-normal etas and the variances entered here are their squares.
    # These reproduce the paper's quoted BSV percentages directly: the
    # abstract and Results report BSV for CL falling "from 173% to 81%",
    # i.e. 100 * omega_CL for the base (1.73) and final (0.81) models.
    etalvc ~ 0.1089  # Alqurain 2024 Table 2, Final Model row Omega V  = 0.33 (SD) -> variance 0.33^2 = 0.1089 (RSE 17.2%)
    etalcl ~ 0.6561  # Alqurain 2024 Table 2, Final Model row Omega Cl = 0.81 (SD) -> variance 0.81^2 = 0.6561 (RSE 44%)

    # Residual error: proportional only. Monolix's proportional error model is
    # y = f + b * f * e with e ~ N(0, 1), so the reported b IS the linear-scale
    # proportional SD and maps directly onto propSd.
    propSd <- 0.30; label("Proportional residual error (fraction)")  # Alqurain 2024 Table 2, Final Model column, "Residual error, b" = 0.3 (RSE 13.7%)
  })

  model({
    # Individual PK parameters. The covariate model is uncentered and
    # untransformed on both covariates:
    #
    #   CL_i = 3.52 * exp(-0.05 * CRP + 0.0088 * CRCL) * exp(eta_CL)
    #   V_i  = 193.65 * exp(eta_V)
    #
    # This form is confirmed by the paper's own Monte-Carlo output rather than
    # assumed: substituting the nine (CRP, CrCl) scenario pairs of
    # Supplementary Table 9 reproduces the tabulated median clearances across
    # their full ~6000-fold span (3, 0.2, 0.001, 2.2, 0.2, 0.0007, 1.6, 0.1,
    # 0.0005 L/h), and the resulting concentration-time profiles reproduce
    # Figure 6 (every panel starts at 1000/193.65 = 5.2 mg/L after the first
    # 1 g dose and the high-CRP panels plateau near 4 x 5.2 = 20.7 mg/L
    # because clearance is then effectively zero over 48 h). A
    # median-centered or log-transformed covariate form reproduces neither.
    cl <- exp(lcl + etalcl) * exp(e_crp_cl * CRP + e_crcl_cl * CRCL)
    vc <- exp(lvc + etalvc)

    # One-compartment IV disposition with linear elimination (Results,
    # "Population Pharmacokinetic Parameters of Vancomycin": "vancomycin PK
    # was best described by infusion administration, no-delay absorption,
    # one-compartment distribution, and linear elimination").
    kel <- cl / vc

    # Doses are administered into the central compartment as intravenous
    # infusions; the infusion duration is a property of the event table
    # (rate/dur) and is not reported by the source, which had trough-only
    # sampling. Supplementary Table 1 shows the bolus and infusion variants
    # returning identical -2LL, AIC, BIC and parameter estimates, so the
    # infusion duration is not identifiable from these data.
    d/dt(central) <- -kel * central

    # Vancomycin plasma concentration in mg/L (doses in mg, volume in L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
