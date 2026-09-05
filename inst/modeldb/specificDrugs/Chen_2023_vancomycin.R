Chen_2023_vancomycin <- function() {
  description <- "One-compartment IV population PK model for vancomycin in Chinese paediatric patients (0.01-17.2 years) with varying renal function (Chen 2023). Clearance carries fixed allometric body weight (exponent 0.75, reference 70 kg), a sigmoid postmenstrual-age maturation function (TM50 37.0 weeks), a near-proportional power effect of BUN-and-creatinine-based estimated GFR capped at 120 mL/min/1.73 m^2 (exponent 1.01, reference 109 mL/min/1.73 m^2), and a 23.9% reduction in patients who underwent cardiothoracic surgery and were receiving vasoactive agents. Central volume scales linearly with body weight (exponent fixed to 1, reference 70 kg) and increases as serum albumin falls (power exponent 0.279 on the 36.2 g/L over albumin ratio). Residual variability is a combined exponential-plus-additive model, encoded here with the exponential arm replaced by its first-order-equivalent proportional arm because rxode2 cannot combine a log-normal and an additive residual."
  reference <- "Chen J, Huang X, Yu L, Li J, Yang R, Li L, Zhou J, Yao H, Bu S. Vancomycin population pharmacokinetics analysis in Chinese paediatric patients with varying degrees of renal function and ages: development of new practical dosing recommendations. J Antimicrob Chemother. 2023;78(8):2037-2051. doi:10.1093/jac/dkad202"
  vignette <- "Chen_2023_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Vancomycin was given by intermittent IV infusion
  # (Chen 2023 Methods, 'Patients and setting'), so the dose enters `central`
  # directly and there is no depot state. `specimen` is verified: the same
  # paragraph states that SERUM vancomycin concentrations were analysed by
  # HPLC.
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 70 kg. Chen 2023 Table 1 model-group median 11.5 kg (IQR 6.5-20.0, range 1.4-69.0). The allometric exponents are fixed, not estimated: 0.75 on CL (Chen 2023 Equations 1 and 5) and 1 on V (Equation 2; the Methods paragraph after Equation 10 states 'For maturation models III-V, allometric exponents for volume of distribution (V) were fixed to 1').",
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age, computed as gestational age plus postnatal age",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "WEEKS, not the register-default months -- Chen 2023 states TM50 = 37.0 weeks (Table 3) and the sigmoid maturation function of Equation 6 is only meaningful on the week scale, so this model declares and uses weeks throughout per the PAGE register entry's Notes. Chen 2023 Methods: 'Postmenstrual age (PMA) was calculated by adding GA to PNA. For patients > 1 year old, GA was imputed as 40 weeks.' Table 1 model-group median 134 weeks (IQR 66-344, range 34-938); only 3 subjects (0.4%) were below 37 weeks, and the paper's own Monte Carlo simulations excluded them, so the maturation term is poorly informed below term.",
      source_name        = "PMA"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate from the BUN-and-creatinine-based paediatric equation, BSA-normalized: eGFR (mL/min/1.73 m^2) = 40.7 * (height (m) / SCr (mg/dL))^0.64 * (30 / BUN (mg/dL))^0.202",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Chen 2023 Equation 11, restated in the Table 3 'The final model' block. Supply the UNCAPPED value: the model applies the 120 mL/min/1.73 m^2 ceiling internally, because the cap is a fitted feature of the covariate model rather than a data-preparation step (Table 2 Model D, dOFV -856.2 against the maturation model, versus -588.8 for the same equation uncapped as Model C). Reference value 109 mL/min/1.73 m^2 is the Table 1 model-group median for Equation 11 (IQR 87-141, range 12-379). Note the three unit conversions a user must make from the Chen 2023 Table 1 reporting units: height cm -> m, SCr umol/L -> mg/dL (divide by 88.4), BUN mmol/L -> mg/dL (multiply by 2.8). Stored under canonical CRCL, which covers BSA-normalized creatinine-based GFR estimates; the assay form here is the KDIGO-endorsed paediatric equation combining serum creatinine and blood urea nitrogen, which Chen 2023 found outperformed the bedside Schwartz, classic Schwartz and Shull equations (Table 2 Models C-J). Only 7 subjects (1.0%) had eGFR below 30 and the paper's simulations excluded them.",
      source_name        = "eGFR"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 36.2 g/L, the constant printed in the Chen 2023 Table 3 final-model V equation; Table 1 reports a model-group median of 36.3 g/L (IQR 31.5-40.5, range 18.6-57.9), so the reference is the rounded cohort median. The effect enters as (36.2 / ALB)^0.279, i.e. albumin is in the DENOMINATOR, so hypoalbuminaemia raises the central volume. Chen 2023 already reports albumin in SI g/L, matching the canonical unit, so no conversion is applied.",
      source_name        = "ALB"
    ),
    SURG_CARD_VASOACTIVE = list(
      description        = "Cardiothoracic surgery with concomitant vasoactive support",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no cardiothoracic surgery, or cardiothoracic surgery without concomitant vasoactive agents)",
      notes              = "Chen 2023 Methods, 'Patients and setting': 'We defined CTS patients as those who underwent CTS and received concomitant vasoactive agents (including milrinone, dopamine, epinephrine or norepinephrine) during vancomycin therapy.' Both halves of that conjunction are required for a 1; the surgery alone does not qualify. Table 1: 68 of 673 model-group subjects (10.1%) were CTS patients, and none of the 53 external-validation subjects were. Time-fixed per subject in the source analysis. The Chen 2023 Discussion attributes the reduced clearance to low cardiac output state and cardiopulmonary bypass raising acute-kidney-injury risk.",
      source_name        = "CTS"
    )
  )

  # Screened in the Chen 2023 covariate search but NOT retained in the final
  # model, so they are documentation only and are not referenced in model().
  # SCr and BUN additionally enter the retained model indirectly, as the two
  # laboratory inputs to the CRCL column above.
  covariatesDataExcluded <- list(
    HT = list(
      description        = "Body height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Not a covariate of the final model in its own right; it is the third input (converted to metres) to the Equation 11 eGFR carried by CRCL. Chen 2023 Table 1 model-group median 83 cm (IQR 65-113, range 45-178).",
      source_name        = "Height"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested as a standalone power covariate on CL (Chen 2023 Table 2 Model A, dOFV -441.7) and rejected in favour of the Equation 11 eGFR (Model D, dOFV -856.2). It remains an input to the CRCL column, where the equation requires mg/dL (divide the Table 1 umol/L values by 88.4). Table 1 model-group median 23.6 umol/L (IQR 17.8-32.6, range 5.0-274.3).",
      source_name        = "SCr"
    ),
    BUN = list(
      description        = "Blood urea nitrogen",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested as a standalone power covariate on CL (Chen 2023 Table 2 Model B, dOFV -315.7) and rejected in favour of the Equation 11 eGFR. It remains an input to the CRCL column, where the equation requires mg/dL (multiply the Table 1 mmol/L values by 2.8). Table 1 model-group median 3.5 mmol/L (IQR 2.4-5.0, range 0.7-29).",
      source_name        = "BUN"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened as a potential covariate (Chen 2023 Methods, final paragraph of 'PPK modelling') and not retained. Table 1 model group 397 of 673 male (59.0%), i.e. 41.0% female.",
      source_name        = "sex"
    ),
    HCT = list(
      description        = "Haematocrit",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened and not retained. Chen 2023 Table 1 model-group median 29.4% (IQR 26.5-32.9, range 13.1-49.5).",
      source_name        = "HCT"
    ),
    TBIL = list(
      description        = "Total bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened and not retained. Chen 2023 Table 1 model-group median 10.4 umol/L (IQR 6.0-18.4, range 1.0-284.3).",
      source_name        = "TBIL"
    ),
    ALT = list(
      description        = "Alanine aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened and not retained. Chen 2023 Table 1 model-group median 30 U/L (IQR 20-49, range 4-4120).",
      source_name        = "ALT"
    ),
    AST = list(
      description        = "Aspartate aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened and not retained. Chen 2023 Table 1 model-group median 45 U/L (IQR 30-83, range 5-7713).",
      source_name        = "AST"
    ),
    TPRO = list(
      description        = "Total serum protein",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened and not retained; serum albumin was the protein-binding covariate that survived onto V. Chen 2023 Table 1 model-group median 62.0 g/L (IQR 55.8-69.0, range 36.0-108.9).",
      source_name        = "TP"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 673L,
    n_studies        = 1L,
    n_sites          = 1L,
    n_concentrations = 1547L,
    age_range        = "0.01-17.23 years postnatal (eligibility was under 18 years); postmenstrual age 34-938 weeks",
    age_median       = "1.81 years postnatal (IQR 0.51-5.83); postmenstrual age 134 weeks (IQR 66-344)",
    weight_range     = "1.4-69.0 kg",
    weight_median    = "11.5 kg (IQR 6.5-20.0)",
    sex_female_pct   = 41.0,
    race_ethnicity   = c(Chinese = 100),
    disease_state    = "Hospitalized Chinese children treated with intravenous vancomycin for at least 3 days for a serious suspected or proven Gram-positive infection, with at least one vancomycin serum concentration available. 82.0% were admitted to an ICU and 10.1% had undergone cardiothoracic surgery with concomitant vasoactive agents. Leading indications were pneumonia (40.6%), primary bloodstream infection (40.0%) and CNS infection (26.0%). Excluded: ongoing renal replacement therapy and extracorporeal membrane oxygenation.",
    dose_range       = "10-20 mg/kg per dose every 6, 8 or 12 h by intermittent IV infusion over 1 or 3 h; daily dose median 40.0 mg/kg (IQR 39.2-57.8, range 19.7-70.0)",
    regions          = "China (Xinhua Hospital, Shanghai Jiao Tong University School of Medicine, Shanghai)",
    renal_function   = "eGFR by Chen 2023 Equation 11: median 109 mL/min/1.73 m^2 (IQR 87-141, range 12-379). Strata: <30 in 1.0%, 30-60 in 6.7%, 60-90 in 20.4%, 90-120 in 33.4%, >120 in 38.5%",
    notes            = "Retrospective single-centre therapeutic-drug-monitoring cohort collected 2013-06 to 2022-06 (Chen 2023 Table 1). 1547 concentrations comprising 779 troughs and 768 peaks; five values below the 0.3 mg/L detection limit were excluded. Assayed by HPLC over a linear range of 1.5-100 mg/L. Fitted in NONMEM 7.5 with FOCE-I; covariate selection by forward inclusion at p < 0.01 and backward elimination at p < 0.001. Reported eta shrinkage 9% on CL and 31% on V, and 28% epsilon shrinkage on both residual-error terms. A separate 53-subject, 77-concentration external-validation cohort was recruited at Shanghai Children's Medical Center between 2021-09 and 2022-06; it contained no cardiothoracic-surgery patients and only one peak concentration, and was not used for estimation."
  )

  ini({
    # Structural parameters (Chen 2023 Table 3). The reference subject weighs
    # 70 kg, is fully mature, has an eGFR of 109 mL/min/1.73 m^2, an albumin of
    # 36.2 g/L and has not undergone cardiothoracic surgery.
    lcl <- log(7.75); label("Clearance standardized to WT=70 kg, mature, eGFR=109 mL/min/1.73 m^2, non-CTS (L/h)")  # Chen 2023 Table 3: TV(CL) = 7.75 L/h (RSE 2.3%; bootstrap 95% CI 7.43-8.13)
    lvc <- log(36.2); label("Central volume standardized to WT=70 kg and ALB=36.2 g/L (L)")                          # Chen 2023 Table 3: TV(V) = 36.2 L (RSE 1.7%; bootstrap 95% CI 35.0-37.4)

    # Allometric exponents. Chen 2023 prints 0.75 inside Equations 1 and 5
    # without an RSE, and the Methods paragraph following Equation 10 states
    # that the volume exponent was fixed to 1 for maturation models III-V.
    # Model III is the one carried forward into the final model, so both
    # exponents are fixed rather than estimated.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT/70) for CL (unitless)")  # Chen 2023 Equations 1 and 5: (WT/70)^0.75, printed without uncertainty
    e_wt_vc <- fixed(1);    label("Allometric exponent on (WT/70) for V (unitless)")   # Chen 2023 Equation 2 and the Methods statement 'allometric exponents for volume of distribution (V) were fixed to 1'

    # Sigmoid postmenstrual-age maturation on CL (Chen 2023 Equation 6). Both
    # terms are estimated and are reported with an RSE in Table 3.
    pma_tm50 <- 37.0;  label("Postmenstrual age at 50% mature CL (TM50, weeks)")                    # Chen 2023 Table 3: TM50 = 37.0 weeks (RSE 4.6%; bootstrap 95% CI 33.5-40.2)
    pma_hill <- -1.63; label("Hill coefficient of the sigmoid CL maturation function (unitless)")   # Chen 2023 Table 3: Hill = -1.63 (RSE 13%; bootstrap 95% CI -2.11 to -1.24). Negative because Equation 6 writes MF = 1 / (1 + (PMA/TM50)^Hill)

    # Covariate effects on CL.
    e_crcl_cl <- 1.01;  label("Power exponent on (capped eGFR / 109) for CL (unitless)")           # Chen 2023 Table 3: theta_eGFR = 1.01 (RSE 3.4%; bootstrap 95% CI 0.94-1.08)
    e_cts_cl  <- 0.761; label("Multiplicative CL factor for cardiothoracic surgery with vasoactive support (unitless)")  # Chen 2023 Table 3: theta_1 = 0.761 (RSE 3.6%; bootstrap 95% CI 0.709-0.817), i.e. a 23.9% CL reduction, the 'roughly 24% reduction' of the Discussion

    # Covariate effect on V.
    e_alb_vc <- 0.279; label("Power exponent on (36.2 / ALB) for V (unitless)")  # Chen 2023 Table 3: theta_ALB = 0.279 (RSE 26.7%; bootstrap 95% CI 0.135-0.423)

    # Correlated between-subject variability on CL and V. Chen 2023 Table 3
    # reports the two diagonals as %CV and the off-diagonal as a covariance on
    # the OMEGA scale in the same block, so the diagonals are read on the same
    # scale as the covariance: omega^2 = (%CV / 100)^2. See the vignette
    # 'Assumptions and deviations' for the alternative log-normal reading and
    # why it changes nothing material.
    #   CL: 0.275^2 = 0.075625     V: 0.272^2 = 0.073984
    #   implied correlation 0.0536 / sqrt(0.075625 * 0.073984) = 0.717
    etalcl + etalvc ~ c(0.075625,
                        0.053600, 0.073984)  # Chen 2023 Table 3: IIV CL 27.5 %CV (RSE 6.0%), Cov(CL-V) 0.0536 (RSE 10.8%), IIV V 27.2 %CV (RSE 16.2%)

    # Residual variability. Chen 2023 Results, 'PPK modelling': 'The residual
    # variability was best described by a combined additive and exponential
    # model', i.e. Y = F * exp(eps1) + eps2.
    #
    # rxode2 (tested on 5.1.7) cannot express the EXPONENTIAL arm combined with
    # an additive arm: `Cc ~ lnorm(expSd) + add(addSd)` PARSES AND COMPILES, and
    # then fails at rxSolve() with "cannot find additive standard deviation for
    # 'Cc'". Verify at the solve, not the build -- a build-only check succeeds
    # and makes this note look wrong. The exponential arm is
    # therefore encoded as its first-order equivalent, a PROPORTIONAL arm, which
    # is what F * exp(eps1) reduces to for a small eps1: F * exp(eps1) ~=
    # F * (1 + eps1). The approximation costs 0.8% on the residual magnitude --
    # the paper's log-scale SD of 0.177 implies a linear-space CV of
    # sqrt(exp(0.177^2) - 1) = 0.1784, against the 0.177 encoded here. The
    # additive arm is carried unchanged. See the vignette 'Assumptions and
    # deviations' section.
    propSd <- 0.177; label("Proportional residual error (fraction)")  # Chen 2023 Table 3: exponential residual 17.7 %CV (RSE 10.6%; bootstrap 95% CI 15.7-19.6), encoded as the first-order-equivalent proportional arm
    addSd  <- 0.154; label("Additive residual SD (ug/mL)")            # Chen 2023 Table 3: additive residual 0.154 mg/L (RSE 44.5%; bootstrap 95% CI 0.026-0.291); 1 mg/L == 1 ug/mL
  })
  model({
    # 1. Derived covariate terms.

    # Renal-function ceiling. Chen 2023 Table 2 Model D caps the Equation 11
    # eGFR at 120 mL/min/1.73 m^2, and the Discussion explains why: Figure 1
    # shows vancomycin CL rising with eGFR only up to about 120, after which
    # the relationship flattens. The cap is part of the fitted covariate model,
    # so it is applied here rather than being left to the user's data
    # preparation.
    crcl_capped <- min(CRCL, 120)

    # Sigmoid postmenstrual-age maturation, Chen 2023 Equation 6:
    #   MF = 1 / (1 + (PMA / TM50)^Hill)
    # with Hill negative, so MF rises from 0 towards 1 as PMA increases.
    maturation_cl <- 1 / (1 + (PAGE / pma_tm50)^pma_hill)

    # 2. Individual PK parameters, Chen 2023 Table 3 'The final model':
    #   CL = TV(CL) * (WT/70)^0.75 * MF * (eGFR/109)^theta_eGFR * theta_1^CTS
    #   V  = TV(V)  * (WT/70)      * (36.2/ALB)^theta_ALB
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * maturation_cl *
      (crcl_capped / 109)^e_crcl_cl * e_cts_cl^SURG_CARD_VASOACTIVE
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * (36.2 / ALB)^e_alb_vc

    # 3. Micro-constants.
    kel <- cl / vc

    # 4. ODE system. Vancomycin is given by intermittent IV infusion, so the
    # dose enters `central` directly.
    d/dt(central) <- -kel * central

    # 5. Observation and error. Dose in mg and volume in L give mg/L, which is
    # the ug/mL of the declared concentration units.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
