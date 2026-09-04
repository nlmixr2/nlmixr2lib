Blackman_2026_methotrexate <- function() {
  description <- "Three-compartment population PK model for high-dose intravenous methotrexate in adult patients with leukemia, lymphoma, or sarcoma treated at Vanderbilt University Medical Center (Blackman 2026), with body surface area on all six disposition parameters (an estimated exponent on clearance and an exponent fixed at 1 on the central and peripheral volumes and both inter-compartmental clearances), serum creatinine and female sex as additional covariates on clearance, inter-individual variability on all six PK parameters, four-occasion inter-occasion variability on clearance, and proportional residual error. Parameter values are taken from the publication's Table 2 (Final Model column)."
  reference <- paste(
    "Blackman MH, Yelvington B, Beck C, Cortez M, Choi L. (2026).",
    "Development and validation of high-dose methotrexate population",
    "pharmacokinetic models to inform clinical decisions on dosing.",
    "Eur J Clin Pharmacol 82:156.",
    "doi:10.1007/s00228-026-04080-0.",
    sep = " "
  )
  vignette <- "Blackman_2026_methotrexate"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate (Methods 'Data and processing': 'Weight, height, and BSA varied over time and were also included as time-varying covariates in the popPK modeling'). Normalized to the cohort median of 1.97 m^2, which matches the Table 1 median BSA of 1.97 m^2 in both the training and test datasets. Enters every one of the six disposition parameters: the exponent on CL (`e_bsa_cl`) is estimated at 0.61, while the exponents on V1, Q2, V2, Q3, and V3 are fixed at 1 per Methods 'Population PK analysis': 'The effect of BSA on clearance was estimated, whereas its effect on all other PK parameters was fixed with an exponent of 1.' Cohort BSA range 1.32-3.01 m^2 (Table 1).",
      source_name        = "BSA"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate measured concurrently with the methotrexate concentrations; when no concurrent value existed the closest available value was carried (Methods 'Data and processing': 'SCR concentrations, typically measured concurrently with MTX concentrations, were used as a time-varying covariate. If an SCR measurement was not available at the time of MTX sampling, the closest available value was used'). Source values were reported in mg/dL and converted to umol/L by the authors; this model expects umol/L. Normalized to 68.08 umol/L in the Table 2 Final Model equation, close to but not identical with the Table 1 training-set median of 68.5 umol/L (overall median 68.1 umol/L) -- the model's centering constant is the per-record median rather than the per-subject median. Cohort range 44.2-132 umol/L (Table 1). Enters clearance only, as a power term with the estimated exponent `e_creat_cl` = -0.56.",
      source_name        = "SCR"
    ),
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Non-time-varying covariate (Methods 'Data and processing': 'age, which changed minimally over the data collection period, and sex were treated as non-time-varying covariates'). The Table 2 Final Model equation carries the female indicator directly as exp(theta9 * I(female)), so male is the reference category and no value inversion is needed. 60.0% female in the training dataset used to fit this model, 55.8% female overall (Table 1). See the model file comment on `e_sexf_cl` and the vignette Errata for the internal inconsistency between the near-zero point estimate and the reported drop in objective function value.",
      source_name        = "Sex"
    ),
    OCC = list(
      description        = "Integer-valued occasion indicator for inter-occasion-variability multiplexing.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Values 1, 2, 3, 4 identify the high-dose methotrexate treatment cycle within a subject. The paper fits a single inter-occasion variance on clearance shared across all occasions (Table 2 row 'delta CL (%CV)'), and does not state a fixed occasion count; patients contributed a median of 3 cycles (range 1-12) and a median of 4 dosing events (range 1-13) per Table 1. Four occasions are encoded here to span the cohort median, following the registered idiom in Jonsson_2011_ethambutol.R; because the variance is common to every occasion, extending the chain to more occasions is a mechanical copy of the `fix()` lines. Decomposed inside `model()` into binary indicators `oc1` .. `oc4` that multiplex the four IOV etas on log-CL. For single-occasion records pass OCC = 1 so the first IOV eta applies.",
      source_name        = "OCC"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "methotrexate", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "methotrexate", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "methotrexate", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 145L,
    n_studies      = 1L,
    age_range      = "18-84 years (median 56)",
    age_median     = "56 years",
    weight_range   = "42.4-141.0 kg (median 82.4)",
    weight_median  = "82.4 kg",
    bsa_range      = "1.32-2.66 m^2 (median 1.97)",
    bsa_median     = "1.97 m^2",
    sex_female_pct = 60.0,
    renal_function = "Serum creatinine 44.2-132 umol/L (median 68.5). Patients who received glucarpidase were excluded because it alters methotrexate clearance.",
    disease_state  = "Adults receiving a high-dose methotrexate-containing regimen for lymphoma (51.0%), leukemia (40.7%), or sarcoma (6.9%); disease missing for 1.4%.",
    dose_range     = "Intravenous methotrexate 200 mg/m^2 or higher; excluding loading doses the training-set dose was a median of 3.5 g/m^2 (range 0.2-12.1). Infusion duration median 4.7 h (range 2.0-27.9). Dosing intervals typically 2-4 weeks (median 3.4). Infusions were given as short-term single doses (59%), long-term single doses (19%), or a loading dose plus maintenance dose (22%).",
    regions        = "United States (single center: Vanderbilt University Medical Center), November 2017 through December 2022.",
    notes          = "Electronic-health-record cohort of 208 adults split randomly 70:30 into training (N = 145) and test (N = 63) datasets; the model in this file was fit to the TRAINING dataset only, so `n_subjects` is 145 rather than the 208-patient full cohort. Of 2,448 measured methotrexate concentrations, 53 (2.2%) were excluded in processing, leaving 2,395; concentrations beyond 96 h post-dose were excluded as sparse and typically below the quantification limit. Training-set patients contributed a median of 9 drug levels (range 2-52) and 4 dosing events (range 1-13) across a median of 3 cycles (range 1-12). Leucovorin was given to all patients starting 24 h after the methotrexate dose per institutional protocol; it is not represented in the model. Estimation was by SAEM in Monolix 2024R. Demographics from Table 1; final parameter estimates from Table 2 (Final Model column)."
  )

  ini({
    # Structural PK parameters -- Blackman 2026 Table 2, 'Final Model' column
    # (OFV = -2512.1, AIC = -2478.1, BIC = -2407.3). Each value is the typical
    # value at the reference covariates BSA = 1.97 m^2, SCR = 68.08 umol/L,
    # male sex. NONMEM/Monolix V1/V2/V3 map onto the nlmixr2lib canonical
    # volumes vc/vp/vp2 and Q2/Q3 onto q/q2.
    lcl  <- log(9.41)  ; label("Clearance CL at the reference covariates (L/h)")                                  # Table 2 Final Model theta1 = 9.41 (SE 0.37) [8.72,10.16]
    lvc  <- log(30.67) ; label("Central volume of distribution V1 at BSA = 1.97 m^2 (L)")                          # Table 2 Final Model theta2 = 30.67 (SE 1.42) [28.01,33.57]
    lq   <- log(0.87)  ; label("Inter-compartmental clearance Q2 to peripheral1 at BSA = 1.97 m^2 (L/h)")          # Table 2 Final Model theta3 = 0.87 (SE 0.07) [0.73,1.03]
    lvp  <- log(5.62)  ; label("First peripheral volume of distribution V2 at BSA = 1.97 m^2 (L)")                 # Table 2 Final Model theta4 = 5.62 (SE 0.29) [5.09,6.21]
    lq2  <- log(0.14)  ; label("Inter-compartmental clearance Q3 to peripheral2 at BSA = 1.97 m^2 (L/h)")          # Table 2 Final Model theta5 = 0.14 (SE 0.01) [0.11,0.17]
    lvp2 <- log(10.60) ; label("Second peripheral volume of distribution V3 at BSA = 1.97 m^2 (L)")                # Table 2 Final Model theta6 = 10.60 (SE 2.43) [6.88,16.34]

    # Covariate effects on clearance. The Table 2 Final Model closed form is
    #   CL = theta1 * (BSA/1.97)^theta7 * (SCR/68.08)^theta8 * exp(theta9 * I(female))
    # restated in the Results 'Population PK model' display equation as
    #   CL_ij = theta1 * (BSA_ij/1.97)^theta7 * (SCR_ij/68.08)^theta8
    #           * exp(theta9 * I(female)_i) * exp(eta_iCL + delta_iCL).
    e_bsa_cl   <-  0.61 ; label("Power exponent on (BSA / 1.97 m^2) for CL (unitless)")                            # Table 2 Final Model theta7 = 0.61 (SE 0.14) [0.39,0.95]
    e_creat_cl <- -0.56 ; label("Power exponent on (CREAT / 68.08 umol/L) for CL (unitless)")                      # Table 2 Final Model theta8 = -0.56 (SE 0.08); CI printed as [-0.42,-0.74], i.e. upper bound first
    # NOTE -- paper-internal inconsistency, transcribed as printed and NOT
    # adjusted. theta9 = 0.01 with SE 0.05 gives female CL only 1.0% above male
    # and a 95% CI [-0.08,0.11] that comfortably spans zero, yet the paper
    # retains sex on the strength of a 25.7-unit OFV drop (Results: adding sex
    # to the SCR model took the OFV reduction from 119.19 to 144.91; Discussion:
    # 'retention of sex in the final model results in a substantial reduction in
    # OFV (25.7 units)'). The estimate and its SE are mutually consistent
    # (0.01 +/- 1.96*0.05 = [-0.088,0.108], matching the printed CI), so the
    # printed point estimate is used. See the vignette Errata section.
    e_sexf_cl  <-  0.01 ; label("Exponential effect of female sex on CL (unitless, male = reference)")             # Table 2 Final Model theta9 = 0.01 (SE 0.05) [-0.08,0.11]

    # BSA exponents on the remaining five disposition parameters, fixed at 1 by
    # the authors rather than estimated. Methods 'Population PK analysis': 'The
    # effect of BSA on clearance was estimated, whereas its effect on all other
    # PK parameters was fixed with an exponent of 1. This approach was selected
    # due to the lack of a clear physiological rationale for alternative BSA
    # scaling and to preserve model parsimony.' Table 2 prints the exponent
    # literally as a superscript 1 on each of the V1/Q2/V2/Q3/V3 rows.
    e_bsa_vc  <- fixed(1) ; label("Power exponent on (BSA / 1.97 m^2) for V1 (unitless)")                          # Table 2 Final Model 'V1 = theta2 x (BSA/1.97)^1'
    e_bsa_q   <- fixed(1) ; label("Power exponent on (BSA / 1.97 m^2) for Q2 (unitless)")                          # Table 2 Final Model 'Q2 = theta3 x (BSA/1.97)^1'
    e_bsa_vp  <- fixed(1) ; label("Power exponent on (BSA / 1.97 m^2) for V2 (unitless)")                          # Table 2 Final Model 'V2 = theta4 x (BSA/1.97)^1'
    e_bsa_q2  <- fixed(1) ; label("Power exponent on (BSA / 1.97 m^2) for Q3 (unitless)")                          # Table 2 Final Model 'Q3 = theta5 x (BSA/1.97)^1'
    e_bsa_vp2 <- fixed(1) ; label("Power exponent on (BSA / 1.97 m^2) for V3 (unitless)")                          # Table 2 Final Model 'V3 = theta6 x (BSA/1.97)^1'

    # Inter-individual variability. Table 2 heads these rows 'omega<param>
    # (%CV)' and the Table 2 footnote states the etas 'follow a lognormal
    # distribution with means zero and covariance matrix with variances along
    # the diagonal'. The tabulated numbers are the lognormal SDs on the log
    # scale, NOT variances: the Results text quotes the same six numbers as
    # percentages ('CL = 9.41 L/hour (24%) ... Q3 = 0.14 L/hour (47%)'), and
    # only the SD reading reproduces them -- for omegaQ3, SD 0.47 gives
    # CV = sqrt(exp(0.47^2) - 1) = 49.7% ~ 47%, whereas reading 0.47 as a
    # variance would give CV = sqrt(exp(0.47) - 1) = 77.5%. The reported SEs
    # agree: omegaCL = 0.24 (SE 0.02) is an 8.3% RSE, below the sqrt(2/145) =
    # 11.8% floor a variance estimate on 145 subjects would carry. nlmixr2
    # takes variances, so each entry below is the tabulated SD squared.
    # Correlated random effects were tested and rejected (Results: 'Allowing
    # for correlation on the random effects did not improve the model'), so
    # these are independent diagonal entries rather than a block.
    etalcl  ~ 0.0576  # Table 2 Final Model omegaCL  = 0.24 (SE 0.02) [0.21,0.29]; variance = 0.24^2 = 0.0576
    etalvc  ~ 0.0841  # Table 2 Final Model omegaV1  = 0.29 (SE 0.03) [0.23,0.36]; variance = 0.29^2 = 0.0841
    etalq   ~ 0.1444  # Table 2 Final Model omegaQ2  = 0.38 (SE 0.12) [0.21,0.69]; variance = 0.38^2 = 0.1444
    etalvp  ~ 0.0900  # Table 2 Final Model omegaV2  = 0.30 (SE 0.05) [0.21,0.42]; variance = 0.30^2 = 0.0900
    etalq2  ~ 0.2209  # Table 2 Final Model omegaQ3  = 0.47 (SE 0.05) [0.39,0.56]; variance = 0.47^2 = 0.2209
    etalvp2 ~ 0.1156  # Table 2 Final Model omegaV3  = 0.34 (SE 0.28) [0.10,1.15]; variance = 0.34^2 = 0.1156

    # Inter-occasion variability on log-CL, the delta_iCL term of the Results
    # display equation. One variance is reported and it is shared by every
    # occasion; nlmixr2 has no NONMEM `$OMEGA BLOCK(1) SAME` shortcut, so the
    # first occasion carries the estimated variance and later occasions fix it
    # to the same value (the registered idiom -- see Jonsson_2011_ethambutol.R,
    # Aregbe_2012_alvespimycin.R, Xie_2019_agomelatine.R).
    etaiov_cl_1 ~ 0.0196       # Table 2 Final Model deltaCL = 0.14 (SE 0.01) [0.13,0.15]; variance = 0.14^2 = 0.0196 (estimated)
    etaiov_cl_2 ~ fix(0.0196)  # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_3 ~ fix(0.0196)  # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_4 ~ fix(0.0196)  # SAME-equivalent: equal to the occasion-1 IOV variance

    # Residual error. The base model was selected with proportional error
    # (Results: 'The selected base model was a three-compartment model with
    # proportional error') and Table 2 reports a single sigma_proportional row
    # for all three models. Reported on the linear (mg/L) scale.
    propSd <- 0.25 ; label("Proportional residual error (fraction)")                                               # Table 2 Final Model sigma_proportional = 0.25 (SE 0.01) [0.24,0.26]
  })

  model({
    # 1. Decompose the integer occasion column into binary indicators to
    #    multiplex the four inter-occasion-variability etas on log-CL. For
    #    single-occasion data pass OCC = 1 so the first IOV eta applies.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)

    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4

    # 2. Individual PK parameters. Clearance carries the estimated BSA
    #    exponent, the serum-creatinine power term, the exponential female-sex
    #    term, its own eta and the occasion-specific IOV eta; the five other
    #    disposition parameters carry BSA with the exponent fixed at 1 plus
    #    their own etas. Reference covariates: BSA 1.97 m^2, SCR 68.08 umol/L,
    #    male sex.
    cl  <- exp(lcl  + etalcl + iov_cl) * (BSA / 1.97)^e_bsa_cl *
      (CREAT / 68.08)^e_creat_cl * exp(e_sexf_cl * SEXF)
    vc  <- exp(lvc  + etalvc)  * (BSA / 1.97)^e_bsa_vc
    q   <- exp(lq   + etalq)   * (BSA / 1.97)^e_bsa_q
    vp  <- exp(lvp  + etalvp)  * (BSA / 1.97)^e_bsa_vp
    q2  <- exp(lq2  + etalq2)  * (BSA / 1.97)^e_bsa_q2
    vp2 <- exp(lvp2 + etalvp2) * (BSA / 1.97)^e_bsa_vp2

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 4. Three-compartment intravenous disposition. Methotrexate is given as an
    #    intravenous infusion, so the dose enters `central` directly and there
    #    is no absorption compartment.
    d/dt(central)     <- -kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # 5. Observation. Dose units mg, vc units L, so Cc is in mg/L -- the unit
    #    the paper reports drug levels in (Table 1 row 'Drug Level (mg/L)').
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
