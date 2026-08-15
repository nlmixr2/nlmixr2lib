Roberts_2025_remdesivir <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for intravenous remdesivir",
    "and its circulating nucleoside metabolite GS-441524 in hospitalised",
    "adults with COVID-19 (Roberts 2025). Each compound is described by a",
    "one-compartment model with first-order elimination, coupled in series,",
    "so the system as a whole is the two-compartment model the paper",
    "describes ('one compartment for each compound'). The entire model is",
    "written in molar units: the source converted both analytes from ng/mL",
    "to uM before fitting so that the parent-to-metabolite flux carries 1:1",
    "stoichiometry despite the molecular-weight difference (remdesivir",
    "602.585 g/mol, GS-441524 291.26 g/mol). Ten percent of remdesivir",
    "clearance was fixed to renal excretion of unchanged drug, so the",
    "remaining 90% (fm) is the flux that forms GS-441524. The fraction of",
    "that flux which actually appears as measurable plasma GS-441524 is not",
    "identifiable from plasma data alone, so the metabolite clearance and",
    "volume are apparent parameters (source CL/fm and V/fm) and the",
    "metabolite state carries an apparent amount; metabolite concentrations",
    "are nevertheless predicted correctly because the apparent volume and",
    "the apparent amount are scaled by the same unknown factor. Retained",
    "covariates are eGFR on GS-441524 clearance (power, exponent 1.12,",
    "reference 80 mL/min/1.73 m^2) and age on GS-441524 volume (power,",
    "exponent -1.15, reference 68.5 years); no covariate was retained on",
    "either remdesivir parameter. Between-subject variability was estimated",
    "on the clearance and volume of both compounds. Residual error is",
    "combined proportional-plus-additive for remdesivir and proportional",
    "for GS-441524.")
  reference <- "Roberts DM, Liu X, Parker SL, Burke A, Peek J, Carland JE, et al. Population Pharmacokinetic Modelling of Remdesivir and Its Metabolite GS-441524 in Hospitalised Patients with COVID-19. Clin Pharmacokinet. 2025;64(5):723-735. doi:10.1007/s40262-025-01496-2"
  vignette <- "Roberts_2025_remdesivir"
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # Both states are molar amounts (umol) in plasma. `central_gs441524`
  # holds an APPARENT amount: the true GS-441524 amount divided by the
  # unidentifiable metabolite fraction that the source folds into its
  # CL/fm and V/fm parameters. Dividing that apparent amount by the
  # apparent volume `vc_gs441524` returns the true plasma concentration,
  # because the same unknown factor scales numerator and denominator.
  # Verified against the MLXTRAN source listing in the Supplementary
  # Material ("MONOLIX model file", states A1 and A2).
  compartmentData <- list(
    central          = list(analyte = "remdesivir", units = "umol", specimen = "plasma", verified = TRUE),
    central_gs441524 = list(analyte = "GS-441524", units = "umol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Estimated glomerular filtration rate calculated from serum creatinine with the CKD-EPI equation, BSA-normalised.",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Roberts 2025 Methods 2.1: 'Estimated glomerular filtration rate (eGFR) was calculated on the basis of serum creatinine concentrations using the chronic kidney disease epidemiology collaboration formula (CKD-EPI) formula', with caution applied in patients with significant acute kidney injury. Patients with eGFR < 30 mL/min/1.73 m^2 were excluded, so the model is not informed below that value. Enters as a power effect on GS-441524 apparent clearance only, normalised to 80 mL/min/1.73 m^2 (the whole-cohort median): (CL/fm)_i = 15.9 * (eGFR_i / 80)^1.12. Baseline (not time-varying) eGFR was used.",
      source_name        = "eGFR"
    ),
    AGE = list(
      description        = "Age at enrolment.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Roberts 2025 Results 3.2 covariate equation: (V/fm)_i = 429 * (Age_i / 68.5)^-1.15. Note the normalising constant is 68.5 years, which is NOT the median age quoted elsewhere in the paper (70 years for the whole cohort, 69 years for the model-building subset); 68.5 is the value printed in the covariate equation and is what reproduces the reported typical volume, so it is used here verbatim. Cohort age range 25-97 years.",
      source_name        = "Age"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 33L,
    n_studies      = 1L,
    n_observations = "Model-building dataset: 49 remdesivir plasma concentrations from 12 patients and 153 GS-441524 plasma concentrations from 25 patients (Roberts 2025 Results 3.2). The 8 remaining patients formed an external-validation set used only to assess predictive performance (GS-441524 concentrations only).",
    age_range      = "25-97 years (whole cohort; median 70, IQR 60-77). Model-building subset median 69 years (IQR 60-73; range 25-97); external-validation subset median 78 years (IQR 64-80; range 25-86).",
    age_median     = "70 years",
    weight_range   = "44-132 kg (model-building subset; median 92 kg, IQR 65-104). Weight was not recorded for the external-validation subset.",
    weight_median  = "92 kg",
    sex_female_pct = 36.4,
    disease_state  = "Adults (>18 years) admitted to hospital with SARS-CoV-2 infection confirmed by nucleic acid amplification test within 14 days of symptom onset, prescribed remdesivir for treatment of COVID-19. Worst disease severity in model-building survivors: mild 20%, moderate 16%, severe 40%, critical 20%. 72% received oxygen therapy, 30% were admitted to intensive care and 12% received ventilatory support. No patient had pre-existing chronic liver disease or received kidney replacement or vasopressor therapy.",
    renal_function = "Median baseline eGFR 80 mL/min/1.73 m^2 (range 33-124) in the whole cohort; model-building subset median 72 (IQR 59-87; range 34-124). Chronic kidney disease (eGFR 30-60) in 40% of the model-building subset. Patients with eGFR < 30 mL/min/1.73 m^2 were excluded from the study, so the model is uninformed below that value.",
    co_medication  = "Corticosteroids (mostly dexamethasone) 88%, baricitinib 48%, tocilizumab 8%, sotrovimab 8%, antibiotics 24% (model-building subset; Table 1).",
    dose_range     = "Licensed regimen only: remdesivir 200 mg intravenously on day 1 followed by 100 mg once daily from day 2 up to day 5 or day 10, each administered as a 60-minute infusion. No alternative regimen was studied clinically; the alternative regimens explored in the paper are simulation only.",
    regions        = "Four hospitals in Australia (St Vincent's Hospital Sydney, The University of Queensland, Royal Brisbane and Women's Hospital and one further site), July 2021 to August 2022.",
    notes          = "Prospective, open-label, multi-centre, observational study. Covariates tested but not retained: gender, height, body weight, BMI, SOFA score, serum albumin, bilirubin, ALT, AST, ALP and GGT (Roberts 2025 Methods 2.4). Baseline demographics are in Roberts 2025 Table 1. Estimation was by SAEM in Monolix 2021R2; the final model was checked by a 1000-run bootstrap (Rsmlx 4.0.2) and externally validated against the held-out 8 patients (median prediction error -15.2%, mean -19.5%, RMSE 30.5%)."
  )

  ini({
    # ------------------------------------------------------------------
    # REMDESIVIR (PARENT) STRUCTURAL PARAMETERS
    # ------------------------------------------------------------------
    # Roberts 2025 Table 2, "Final model / Mean" column. Remdesivir is
    # given as an intravenous infusion, so CL and V are true (not
    # apparent) parameters. Volumes are in L and clearances in L/h;
    # amounts are in umol, so central/vc returns uM directly and no
    # unit multiplier is needed.

    lcl <- log(105)
    label("Remdesivir clearance CL (L/h)")                                       # Roberts 2025 Table 2: remdesivir CL = 105 L/h (RSE 17.8%; bootstrap median 100, 95% CI 69.3-154)

    lvc <- log(121)
    label("Remdesivir central volume of distribution V (L)")                     # Roberts 2025 Table 2: remdesivir V = 121 L (RSE 23.7%; bootstrap median 111, 95% CI 45.9-195.5)

    # Fraction of remdesivir clearance routed down the metabolic pathway
    # that forms GS-441524. Fixed, not estimated: the source fixed renal
    # clearance of unchanged remdesivir at 10% of total clearance on the
    # basis of external mass-balance data, leaving 90% metabolic. The
    # supplementary MLXTRAN listing codes this literally as
    # "CLr = 0.1*CL" and "CLtr = 0.9*CL", with only CLtr feeding the
    # metabolite state.
    #
    # CAUTION on symbols: this `fm` is NOT the same quantity as the
    # "fm" appearing in the source's CL/fm and V/fm parameter names.
    # This one is the metabolic (non-renal) fraction of remdesivir
    # clearance and is fixed at 0.9. The source's fm is the further,
    # unidentifiable fraction of metabolised remdesivir that appears in
    # plasma as measurable GS-441524; it is never estimated and is
    # absorbed into the apparent metabolite clearance and volume below.
    fm <- fixed(0.9)
    label("Fraction of remdesivir clearance forming GS-441524 (unitless)")  # Roberts 2025 Methods 2.3: "As 10% of the administered remdesivir is excreted unchanged in the urine [19], renal clearance of remdesivir was fixed to 10% of the total remdesivir clearance." Supplementary Material MLXTRAN listing: CLr = 0.1*CL; CLtr = 0.9*CL

    # ------------------------------------------------------------------
    # GS-441524 (METABOLITE) STRUCTURAL PARAMETERS
    # ------------------------------------------------------------------
    # Apparent parameters: the source reports CL/fm and V/fm because the
    # metabolite fraction is not identifiable from plasma data alone.
    # Both are carried here exactly as reported; predicted metabolite
    # concentrations are unaffected by the unknown scaling (see the
    # compartmentData note above). Values are the typical values at the
    # covariate reference points (eGFR 80 mL/min/1.73 m^2, age 68.5 y).

    lcl_gs441524 <- log(15.9)
    label("GS-441524 apparent clearance CL/fm (L/h) at eGFR 80 mL/min/1.73 m^2")  # Roberts 2025 Table 2: GS-441524 CL/fm = 15.9 L/h (RSE 8.39%; bootstrap median 16.0, 95% CI 13.2-19.6). Reference eGFR from the Results 3.2 covariate equation.

    lvc_gs441524 <- log(429)
    label("GS-441524 apparent volume of distribution V/fm (L) at age 68.5 years")  # Roberts 2025 Table 2: GS-441524 V/fm = 429 L (RSE 11.4%; bootstrap median 429, 95% CI 353-512). Reference age from the Results 3.2 covariate equation.

    # ------------------------------------------------------------------
    # COVARIATE EFFECTS
    # ------------------------------------------------------------------
    # Both retained covariates act on GS-441524 only, both as power
    # models, per the two equations printed in Roberts 2025 Results 3.2:
    #   (Cl/fm)_i = 15.9 * (eGFR_i / 80)^1.12
    #   (V/fm)_i  = 429  * (Age_i  / 68.5)^-1.15
    # No covariate was retained on remdesivir CL or V.

    e_crcl_cl_gs441524 <- 1.12
    label("Power exponent for eGFR on GS-441524 apparent clearance (unitless)")   # Roberts 2025 Table 2 "eGFR effect on CL" = 1.12 (bootstrap median 1.23, 95% CI 0.78-1.89) and the Results 3.2 equation exponent

    e_age_vc_gs441524 <- -1.15
    label("Power exponent for age on GS-441524 apparent volume (unitless)")       # Roberts 2025 Table 2 "Age effect on V" = -1.15 (bootstrap median -1.18, 95% CI -2.82 to -0.51) and the Results 3.2 equation exponent

    # ------------------------------------------------------------------
    # BETWEEN-SUBJECT VARIABILITY
    # ------------------------------------------------------------------
    # Roberts 2025 Methods 2.3: "All individual parameters were assumed
    # to be log-normally distributed. The between-subject variability
    # (BSV or omega) was described using an exponential model." Table 2
    # heads the block "Between subject variability (%)", so the reported
    # values are read as %CV of a log-normal and converted with the
    # standard identity omega^2 = log(1 + CV^2):
    #   remdesivir CL   53.1% -> log(1 + 0.531^2) = 0.24839
    #   remdesivir V    61.4% -> log(1 + 0.614^2) = 0.31990
    #   GS-441524 CL/fm 36.8% -> log(1 + 0.368^2) = 0.12701
    #   GS-441524 V/fm  46.0% -> log(1 + 0.460^2) = 0.19194
    # See the vignette Errata: the alternative reading (the values are
    # omega itself, expressed as a percentage) cannot be excluded from
    # the paper's wording, but changes only the width of the simulated
    # prediction intervals, never a typical-value prediction. Monolix
    # reported no BSV correlations, so the matrix is diagonal.

    etalcl ~ 0.24839                                                             # Roberts 2025 Table 2: BSV remdesivir CL = 53.1% CV (RSE 29.5%, shrinkage 10.0%) -> omega^2 = log(1 + 0.531^2)
    etalvc ~ 0.31990                                                             # Roberts 2025 Table 2: BSV remdesivir V = 61.4% CV (RSE 31.2%, shrinkage 12.8%) -> omega^2 = log(1 + 0.614^2)
    etalcl_gs441524 ~ 0.12701                                                    # Roberts 2025 Table 2: BSV GS-441524 CL/fm = 36.8% CV (RSE 17.7%, shrinkage 6.91%) -> omega^2 = log(1 + 0.368^2)
    etalvc_gs441524 ~ 0.19194                                                    # Roberts 2025 Table 2: BSV GS-441524 V/fm = 46.0% CV (RSE 22.1%, shrinkage -5.52%) -> omega^2 = log(1 + 0.460^2)

    # ------------------------------------------------------------------
    # RESIDUAL ERROR
    # ------------------------------------------------------------------
    # Roberts 2025 Results 3.2: "The combined (proportional plus
    # additive) and proportional error model were used to describe
    # residual variability for remdesivir and GS-441524, respectively."
    # The additive term is on the molar concentration scale the model
    # was fitted on, so its units are uM.

    addSd <- 0.014
    label("Remdesivir additive residual SD (uM)")                                # Roberts 2025 Table 2: remdesivir additive residual = 0.014 uM (RSE 38.1%; bootstrap median 0.016, 95% CI 0.0008-0.039)

    propSd <- 0.62
    label("Remdesivir proportional residual SD (fraction)")                      # Roberts 2025 Table 2: remdesivir proportional error = 0.62 (RSE 17.2%; bootstrap median 0.61, 95% CI 0.44-0.74)

    propSd_gs441524 <- 0.16
    label("GS-441524 proportional residual SD (fraction)")                       # Roberts 2025 Table 2: GS-441524 proportional error = 0.16 (RSE 7.06%; bootstrap median 0.16, 95% CI 0.11-0.20)
  })

  model({
    # ------------------------------------------------------------------
    # Individual parameters
    # ------------------------------------------------------------------
    # Covariate effects are applied as multiplicative power terms
    # exactly as printed in Roberts 2025 Results 3.2. No covariate acts
    # on remdesivir.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    cl_gs441524 <- exp(lcl_gs441524 + etalcl_gs441524) *
      (CRCL / 80)^e_crcl_cl_gs441524
    vc_gs441524 <- exp(lvc_gs441524 + etalvc_gs441524) *
      (AGE / 68.5)^e_age_vc_gs441524

    # ------------------------------------------------------------------
    # Micro-constants
    # ------------------------------------------------------------------
    # Supplementary MLXTRAN listing, transcribed one-for-one:
    #   CLr  = 0.1*CL ; k10 = CLr/V1    (renal loss of unchanged drug)
    #   CLtr = 0.9*CL ; k12 = CLtr/V1   (conversion to GS-441524)
    #   k20  = CLmi/V2                  (GS-441524 elimination)
    # k10 + k12 is the total remdesivir elimination rate constant, so
    # the parent decays with CL/V regardless of how the split is made;
    # the split determines only how much flux reaches the metabolite.
    kel          <- cl / vc
    kform        <- fm * kel
    kel_gs441524 <- cl_gs441524 / vc_gs441524

    # ------------------------------------------------------------------
    # ODE system
    # ------------------------------------------------------------------
    # Remdesivir is dosed into `central` as an intravenous infusion.
    # Both states are in umol and the source fitted in molar units, so
    # the parent-to-metabolite flux is 1:1 with no molecular-weight
    # factor. Doses expressed in mg must be converted before use:
    # umol = mg / 602.585 * 1000.
    d/dt(central)          <- -kel * central
    d/dt(central_gs441524) <-  kform * central - kel_gs441524 * central_gs441524

    # ------------------------------------------------------------------
    # Observations
    # ------------------------------------------------------------------
    # Amounts in umol divided by volumes in L give uM, the scale the
    # source fitted and reports on.
    Cc          <- central          / vc
    Cc_gs441524 <- central_gs441524 / vc_gs441524

    Cc          ~ prop(propSd) + add(addSd)
    Cc_gs441524 ~ prop(propSd_gs441524)
  })
}
