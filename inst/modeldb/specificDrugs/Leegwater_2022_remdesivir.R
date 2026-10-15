Leegwater_2022_remdesivir <- function() {
  description <- paste(
    "Integrated parent-metabolite population PK model for intravenous",
    "remdesivir and its circulating nucleoside metabolite GS-441524 in",
    "non-critically ill hospitalized adults with COVID-19 and hypoxemia",
    "(Leegwater 2022). One compartment for each analyte, coupled in series",
    "(NONMEM ADVAN5). The estimated remdesivir clearance is the metabolic",
    "arm that forms GS-441524; renal excretion of unchanged remdesivir is",
    "fixed at 10% of that metabolic clearance, so total remdesivir",
    "clearance is 1.1 times the metabolic clearance and the two arms share",
    "one random effect. BSA-normalized eGFR (CKD-EPI) enters GS-441524",
    "clearance as a power function referenced to 94 mL/min/1.73 m^2. The",
    "parent-to-metabolite transfer is mass for mass with no molecular-weight",
    "correction, as in the source control stream, so the GS-441524 volume",
    "and clearance are apparent values that absorb the molar-mass ratio and",
    "any unformed fraction. Residual error is additive on the linear",
    "concentration scale for each analyte.",
    sep = " "
  )
  reference <- paste(
    "Leegwater E, Moes DJAR, Bosma LBE, Ottens TH, van der Meer IM,",
    "van Nieuwkoop C, Wilms EB. Population Pharmacokinetics of Remdesivir",
    "and GS-441524 in Hospitalized COVID-19 Patients. Antimicrob Agents",
    "Chemother. 2022;66(6):e00254-22. doi:10.1128/aac.00254-22.",
    "Final estimates from Table 2; the compartment structure, the",
    "0.1 x metabolic-clearance renal arm, the eGFR power term, the",
    "OMEGA variances and the additive residual-error form from the NONMEM",
    "control stream in Supplement 2 ('FINAL model RDV+ GS-441524', ADVAN5",
    "TRANS1).",
    sep = " "
  )
  vignette <- "Leegwater_2022_remdesivir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description = "Estimated glomerular filtration rate calculated with the CKD-EPI equation and reported BSA-normalized",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters GS-441524 clearance only, as a power function normalized to",
        "94 mL/min/1.73 m^2, the cohort median (Table 1; control stream",
        "'TVCLm = THETA(3)*((GFR/94)**THETA(7))'). Baseline value; the",
        "study sampled only the first dosing day. Cohort range 8-119",
        "mL/min/1.73 m^2 (Table 1), so the model is informed by very few",
        "patients with severe renal impairment. Results: adding eGFR",
        "lowered the objective function by 20 points and explained 66% of",
        "the IIV in GS-441524 clearance.",
        sep = " "
      ),
      source_name = "GFR"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested as a covariate (Methods, 'Pharmacokinetic analysis') but not retained in the final model."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested (including allometric scaling to 70 kg) but not retained in the final model."
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Tested but not retained in the final model."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Tested but not retained in the final model."
    ),
    CRP = list(
      description = "C-reactive protein",
      units = "mg/L",
      type = "continuous",
      notes = "Tested but not retained in the final model."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested but not retained in the final model."
    )
  )

  compartmentData <- list(
    central = list(analyte = "remdesivir", units = "mg", specimen = "plasma", verified = TRUE),
    central_gs441524 = list(
      analyte = "GS-441524",
      units = "mg (remdesivir mass equivalents)",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 17L,
    n_studies = 1L,
    n_observations = "84 blood samples, each assayed for remdesivir and GS-441524, all drawn on the first day of therapy. 34% of remdesivir concentrations were below the limit of quantification (4 ug/L) and 25% below the limit of detection (1 ug/L); none of the GS-441524 concentrations were below the limit of quantification (12 ug/L).",
    age_range = "31-74 years",
    age_median = "55 years",
    weight_range = "65-122 kg",
    weight_median = "92 kg",
    sex_female_pct = 5.9,
    race_ethnicity = "Not reported.",
    disease_state = paste(
      "Adults hospitalized on the general ward with RT-PCR-confirmed",
      "COVID-19 who needed supplemental oxygen (WHO ordinal scale 5 for all",
      "patients; median oxygen requirement 4 L/min, range 1-15). All",
      "received concomitant dexamethasone. Comorbidities: diabetes",
      "mellitus 35%, cardiovascular disease 29%. Three patients (17.6%)",
      "were later admitted to the ICU and one died.",
      sep = " "
    ),
    renal_function = "eGFR (CKD-EPI) median 94 mL/min/1.73 m^2 (range 8-119); serum creatinine median 75 umol/L (range 46-573).",
    hepatic_function = "ALT median 36 U/L (range 20-150); total bilirubin median 8 umol/L (range 2-18); albumin median 37 g/L (range 30-47).",
    dose_range = "Remdesivir 200 mg intravenously on day 1 followed by 100 mg once daily on days 2-5, each infused over 1 to 2 h.",
    regions = "The Netherlands (Haga Teaching Hospital, The Hague), January to July 2021.",
    notes = paste(
      "Prospective observational PK study. Six samples were scheduled at",
      "0.5, 1.5, 2.5, 6, 12 and 23 h after the end of the first infusion.",
      "Baseline demographics from Table 1. Median BMI 30.86 kg/m^2 (range",
      "21.72-41.21) and median BSA 2.11 m^2 (range 1.77-2.52). Below-LOQ",
      "remdesivir data were modelled with the all-data method of Keizer et",
      "al. Estimation by FOCE-I in NONMEM 7.4; parameter precision from a",
      "1,000-sample nonparametric bootstrap.",
      sep = " "
    )
  )

  ini({
    # Remdesivir (parent). Table 2 'Final model / Mean value' column,
    # identical to $THETA(1)-$THETA(2) of the Supplement 2 control stream.
    # The control stream's own comment calls THETA(1) the 'non metabolism'
    # clearance, but the code uses it as CLtoM (K12 = CLtoM/V1, the flux
    # into GS-441524) and Table 2 labels it 'Metabolic CL'; the code is
    # what was fitted.
    lcl_met <- log(207); label("Remdesivir metabolic (GS-441524-forming) clearance (L/h)")                           # Table 2 'Metabolic CL (L/h)' 207, RSE 13%, bootstrap median 209 (95% CI 152-279); Supplement 2 $THETA(1) 207, CLtoM = TVCL*EXP(ETA(1))
    lvc <- log(157); label("Remdesivir volume of distribution (L)")                                                     # Table 2 'V (L)' 157, RSE 19%, bootstrap median 164 (95% CI 98.9-259); Supplement 2 $THETA(2) 157

    # Renal (unchanged) remdesivir clearance is not a separate THETA: the
    # control stream sets CL = 0.1 * CLtoM, so the metabolic arm is 10
    # times the renal arm and Table 2 reports the renal arm as
    # 'Renal CL 20.7 fixed' (= 0.1 * 207). Results: fixed because about
    # 10% of administered remdesivir is excreted unchanged in urine.
    clrat_gs441524 <- fixed(10); label("Ratio of GS-441524-forming clearance to unchanged renal remdesivir clearance (unitless)") # Supplement 2 'CL = 0.1 * CLtoM' -> ratio 1/0.1 = 10; Table 2 'Renal CL (L/h) 20.7 fixed' = 207/10

    # GS-441524 (metabolite). Table 2 and $THETA(3), $THETA(4), $THETA(7).
    lcl_gs441524 <- log(27.6); label("GS-441524 clearance at eGFR 94 mL/min/1.73 m^2 (L/h)")                        # Table 2 GS-441524 'CL (L/h)' 27.6, RSE 17%, bootstrap median 28.1 (95% CI 20.7-38.9); Supplement 2 $THETA(3) 27.6
    lvc_gs441524 <- log(1060); label("GS-441524 volume of distribution (L)")                                            # Table 2 GS-441524 'V (L)' 1,060, RSE 11%, bootstrap median 1,062 (95% CI 834-1,270); Supplement 2 $THETA(4) 1060
    e_crcl_cl_gs441524 <- 1.76; label("Power exponent on (CRCL/94) for GS-441524 clearance (unitless)")              # Table 2 'eGFR on CL' 1.76, RSE 22%, bootstrap median 1.68 (95% CI 0.31-2.91); Supplement 2 $THETA(7) 1.76, (GFR/94)**THETA(7)

    # Between-subject variability. Every parameter enters as TV*EXP(ETA),
    # and the Table 2 percentages are sqrt(omega^2) of the $OMEGA block:
    # sqrt(0.151) = 38.9%, sqrt(0.229) = 47.9%, sqrt(0.225) = 47.4%,
    # sqrt(0.184) = 42.9%. ini() takes omega^2 directly from the control
    # stream. The block is diagonal. The remdesivir eta sits on the
    # metabolic clearance and so also scales the renal arm.
    etalcl_met ~ 0.151                                                                                                  # Supplement 2 $OMEGA(1) 0.151; Table 2 IIV 'Remdesivir nonrenal CL' 38.9%, RSE 22%, shrinkage 12%
    etalvc ~ 0.229                                                                                                      # Supplement 2 $OMEGA(2) 0.229; Table 2 IIV 'Remdesivir V' 47.9%, RSE 23%, shrinkage 24%
    etalcl_gs441524 ~ 0.225                                                                                             # Supplement 2 $OMEGA(3) 0.225; Table 2 IIV 'GS-441524 CL' 47.4%, RSE 29%, shrinkage 28%
    etalvc_gs441524 ~ 0.184                                                                                             # Supplement 2 $OMEGA(4) 0.184; Table 2 IIV 'GS-441524 V' 42.9%, RSE 17%, shrinkage 0%

    # Residual error. Supplement 2 $ERROR: Y = IPRED + W*EPS(1) with
    # $SIGMA 1 FIX and W = THETA(5) for remdesivir, THETA(6) for
    # GS-441524 (CMT 2), i.e. a pure additive SD on the linear mg/L scale.
    addSd <- 0.0294; label("Remdesivir additive residual SD (mg/L)")                                                   # Table 2 residual variability 'Remdesivir' 0.0294, RSE 12%, bootstrap median 0.0269 (95% CI 0.0091-0.0468); Supplement 2 $THETA(5)
    addSd_gs441524 <- 0.0140; label("GS-441524 additive residual SD (mg/L)")                                           # Table 2 residual variability 'GS-441524' 0.0140, RSE 11%, bootstrap median 0.0136 (95% CI 0.0083-0.0189); Supplement 2 $THETA(6)
  })

  model({
    # Individual parameters (Supplement 2 $PK).
    cl_met <- exp(lcl_met + etalcl_met) # CLtoM
    cl_nonmet <- cl_met / clrat_gs441524 # CL = 0.1 * CLtoM
    vc <- exp(lvc + etalvc) # V1
    cl_gs441524 <- exp(lcl_gs441524 + etalcl_gs441524) * (CRCL / 94)^e_crcl_cl_gs441524 # CLm
    vc_gs441524 <- exp(lvc_gs441524 + etalvc_gs441524) # V2

    # ADVAN5 rate constants: K10 = CL/V1, K12 = CLtoM/V1, K20 = CLm/V2.
    kel <- cl_nonmet / vc
    kform <- cl_met / vc
    kel_gs441524 <- cl_gs441524 / vc_gs441524

    # Remdesivir is infused into `central`. The transfer to GS-441524 is
    # amount for amount, exactly as ADVAN5's K12 is: the source applies no
    # molar-mass correction between remdesivir (602.6 g/mol) and GS-441524
    # (291.3 g/mol), so `central_gs441524` carries remdesivir-mass
    # equivalents and the apparent volume `vc_gs441524` absorbs the ratio.
    # Concentrations come out in the units the source fitted (mg/L).
    d/dt(central) <- -(kel + kform) * central
    d/dt(central_gs441524) <- kform * central - kel_gs441524 * central_gs441524

    Cc <- central / vc
    Cc_gs441524 <- central_gs441524 / vc_gs441524

    Cc ~ add(addSd)
    Cc_gs441524 ~ add(addSd_gs441524)
  })
}
