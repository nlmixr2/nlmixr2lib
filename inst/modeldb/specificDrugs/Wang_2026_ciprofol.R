Wang_2026_ciprofol <- function() {
  description <- "Three-compartment intravenous population PK model for ciprofol (HSK3486) after a single 0.6 mg/kg bolus given over 30 s in Chinese pediatric surgical patients (Wang 2026; 25 children aged 1-9 years, ASA physical status I-II, scheduled for elective urologic surgery; 317 arterial plasma samples). A three-compartment model was significantly better than a two-compartment model (dOFV = 65.8, p < 0.001). All disposition parameters are reported per kilogram of body weight (CL 31.2 mL/min/kg, V1 506 mL/kg, Q2 28.2 mL/min/kg, V2 231 mL/kg, Q3 19.8 mL/min/kg, V3 1360 mL/kg), i.e. body weight enters every parameter with a linear (exponent 1) scaling; standard allometric scaling and age-dependent maturation functions on clearance were tested and did not improve the fit, and no further effect of weight, age, sex or BMI was detectable after per-kilogram normalisation. Blood urea nitrogen was the single retained covariate, acting on the central volume as a power of the ratio to the cohort median with an estimated exponent of -0.821, so that V1 falls from 0.770 L/kg at BUN 3 mmol/L to 0.384 L/kg at BUN 7 mmol/L; the authors judged the resulting exposure change clinically insignificant. Log-normal inter-individual variability was retained on CL, V1 (correlated, r = -0.821) and Q3 only, because the IIV estimates for V2, V3 and Q2 were close to zero. Residual error is combined proportional plus additive."
  reference <- paste(
    "Wang S, Li Y, Hu Z, Du L, Wang Y, Jiang X, Li L, Shangguan W. (2026).",
    "Population pharmacokinetics of a single bolus of ciprofol in Chinese",
    "pediatric patients. BMC Anesthesiology 26(1).",
    "doi:10.1186/s12871-026-03647-9. PMCID: PMC12930602.",
    "Chinese Clinical Trial Registry ChiCTR2200058405.",
    sep = " "
  )
  vignette <- "Wang_2026_ciprofol"
  units    <- list(time = "min", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: the paper measures ciprofol in
  # plasma obtained from arterial blood samples (Methods, 'Blood sampling
  # and sample handling' / 'Determination of ciprofol concentration'), and
  # the disposition compartments hold unchanged parent ciprofol.
  compartmentData <- list(
    central     = list(analyte = "ciprofol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ciprofol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "ciprofol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline; time-fixed over the 180 min sampling window).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Structural, not a screened covariate effect: Wang 2026 reports every disposition parameter per kilogram (Table 3 rows are 'CL, ml/min/kg', 'V1, ml/kg', etc.), so body weight multiplies CL, V1, Q2, V2, Q3 and V3 with a linear exponent of 1 and no reference weight. This is a reparameterisation rather than an estimated covariate relationship -- concentrations are invariant to WT provided the dose is also given per kilogram. The paper explicitly tested standard allometric scaling and established age-dependent maturation functions for clearance as structural covariates and did not retain them (Methods, 'Model selection and covariate analysis'), and found no further significant weight effect after per-kilogram dosing (Results, 'Population covariant analysis'). Source range 9.8-37 kg, mean 18.5 kg (Table 1).",
      source_name        = "Weight"
    ),
    BUN = list(
      description        = "Baseline blood urea nitrogen, the single covariate retained in the final model.",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Acts on the central volume as a power of the ratio to the cohort median (Methods Eq. 5 continuous-covariate power model; Results 'Population covariant analysis'): V1 = V1_TV * (BUN / 5)^(-0.821). Adding BUN on V1 decreased the objective function by 9.55 units (p < 0.005) and reduced the V1 inter-individual variability from 41.6% to 35.1%. The paper does not print the median BUN used as the normalising constant -- Table 1 reports the cohort MEAN as 5.2 mmol/L -- but it is recovered exactly from the two V1 values the paper does print: V1 = 0.770 L/kg at BUN 3 mmol/L and V1 = 0.384 L/kg at BUN 7 mmol/L (Results, 'Population covariant analysis'). Solving 506 * (3/ref)^(-0.821) = 770 gives ref = 5.003, and 506 * (7/ref)^(-0.821) = 384 gives ref = 5.002; Table 4 independently simulates at BUN = 3, 5 and 7 mmol/L with 5 as the central level. The normalising constant is therefore 5 mmol/L. Note that the printed exponent -0.821 coincides numerically with the printed CL-V1 IIV correlation of -0.821; both were verified independently (the exponent by the V1 ratio 0.770/0.384 = 2.005 = (3/7)^(-0.821), the correlation by the Results text 'a correlation coefficient of -0.821 for the IIV of CL and V1'). Source range 2.9-7.6 mmol/L, mean 5.2 mmol/L (Table 1).",
      source_name        = "BUN"
    )
  )

  # Covariates the paper screened but did not retain in the final model
  # (Methods, 'Model selection and covariate analysis'; Results, 'Population
  # covariant analysis': 'After inclusion of BUN as a covariate, there were
  # no further significant effects of weight, age, sex, BMI, or other
  # laboratory parameters detected.'). Documentation only -- these are not
  # referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at surgery.", units = "years", type = "continuous",
      notes = "Screened as a candidate covariate and as the driver of established age-dependent maturation functions on clearance; neither was retained. Cohort 1.0-9.0 years, mean 4.2 years, enrolled in three balanced age strata (toddlers 1-2 y, preschoolers 3-5 y, school-age 6-9 y; Table 1)."
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.", units = "(binary)", type = "binary",
      notes = "Screened, not retained. 14 male / 11 female (Table 1). Source column 'Male (Female)' counts; SEXF = 1 - SEXM."
    ),
    BMI = list(
      description = "Body mass index.", units = "kg/m^2", type = "continuous",
      notes = "Screened, not retained. Enrolment required BMI between the 25th and 75th percentile for age and sex, so the cohort is deliberately non-obese; cohort 13.2-20.5 kg/m^2, mean 16.7 (Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase.", units = "U/L", type = "continuous",
      notes = "Screened, not retained. Cohort 9-27 U/L, mean 14.6 (Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase.", units = "U/L", type = "continuous",
      notes = "Screened, not retained. Cohort 19-50 U/L, mean 30.9 (Table 1)."
    ),
    TBILI = list(
      description = "Total bilirubin.", units = "umol/L", type = "continuous",
      notes = "Screened, not retained. Cohort 3.8-23.6 umol/L, mean 8.4 (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin.", units = "g/L", type = "continuous",
      notes = "Screened, not retained. Cohort 40.3-49.2 g/L, mean 45.5 (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine.", units = "umol/L", type = "continuous",
      notes = "Screened, not retained. Cohort 18.5-58.3 umol/L, mean 32.6 (Table 1). BUN was the only renal-function marker retained."
    ),
    HGB = list(
      description = "Hemoglobin.", units = "g/L", type = "continuous",
      notes = "Screened, not retained. Cohort 110-145 g/L, mean 126.7 (Table 1)."
    ),
    TPRO = list(
      description = "Total serum protein.", units = "g/L", type = "continuous",
      notes = "Screened, not retained. Cohort 57.5-75.4 g/L, mean 68.4 (Table 1)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 25L,
    n_studies      = 1L,
    age_range      = "1.0-9.0 years",
    age_median     = "4.2 years (cohort mean; Table 1 reports mean (min, max), not median)",
    weight_range   = "9.8-37 kg",
    weight_median  = "18.5 kg (cohort mean; Table 1 reports mean (min, max), not median)",
    sex_female_pct = 44,
    disease_state  = "Healthy, non-obese children with American Society of Anesthesiologists physical status I or II scheduled for elective urologic surgery of anticipated duration > 2 h. Exclusions included predicted difficult airway, acute upper respiratory infection, uncontrolled asthma, cardiomyopathy, hepatic disease, renal insufficiency, and sedative-hypnotic use within 7 days before surgery.",
    dose_range     = "Single 0.6 mg/kg intravenous bolus of ciprofol administered over 30 s at induction of general anaesthesia.",
    regions        = "China (Second Affiliated Hospital and Yuying Children's Hospital of Wenzhou Medical University, Wenzhou, Zhejiang), January-August 2023.",
    notes          = "27 children were enrolled in three balanced age strata (toddlers 1-2 years n = 9, preschoolers 3-5 years n = 9, school-age children 6-9 years n = 9); one withdrew when the operation was shortened to under 2 h and one was excluded for a blocked arterial catheter, leaving 25 children (14 male, 11 female; final strata 9 / 9 / 7). Thirteen arterial samples were planned per patient (pre-dose and 2, 4, 6, 8, 10, 20, 30, 45, 60, 90, 120 and 180 min post-injection); 8 of 325 samples were not collected, mostly the 180 min late-elimination sample, giving 317 samples for the analysis. Ciprofol was quantified by UPLC-APCI-MS/MS over 5-20000 ng/mL with an LLOQ of 5 ng/mL; no concentration fell below the limit of quantitation. All patients received midazolam 0.1-0.2 mg/kg before induction and fentanyl 2.0 ug/kg with ciprofol, followed by cisatracurium and sevoflurane maintenance. Fitted in NONMEM 7.4; the final model was evaluated by goodness-of-fit plots, a 1000-replicate bootstrap and a visual predictive check. Trial registration ChiCTR2200058405."
  )

  ini({
    # Final model estimates (Wang 2026 Table 3). Time in minutes, dose in mg,
    # central-compartment concentration in mg/L (= ug/mL). The paper reports
    # every disposition parameter per kilogram of body weight and in mL, so
    # each value below is the printed number divided by 1000 to convert
    # mL -> L; body weight is reapplied in model() with a linear exponent.
    #
    # The paper's three-compartment parameterisation maps to the nlmixr2lib
    # canonical names as:
    #   V1 (paper, central)          -> vc
    #   V2 (paper, fast peripheral)  -> vp
    #   V3 (paper, slow peripheral)  -> vp2
    #   Q2 (paper, central <-> fast) -> q
    #   Q3 (paper, central <-> slow) -> q2
    lcl  <- log(0.0312); label("Clearance per kilogram body weight (L/min/kg)")                          # Wang 2026 Table 3: CL = 31.2 ml/min/kg (RSE 5.0%; bootstrap median 31.2, 2.5th-97.5th 28.2-34.2)
    lvc  <- log(0.506);  label("Central volume of distribution V1 per kilogram body weight (L/kg)")      # Wang 2026 Table 3: V1 = 506 ml/kg (RSE 8.1%; bootstrap median 492, 377-576)
    lq   <- log(0.0282); label("Fast inter-compartmental clearance Q2 per kilogram body weight (L/min/kg)") # Wang 2026 Table 3: Q2 = 28.2 ml/min/kg (RSE 19.9%; bootstrap median 29.1, 20.7-68.6)
    lvp  <- log(0.231);  label("Fast peripheral volume of distribution V2 per kilogram body weight (L/kg)") # Wang 2026 Table 3: V2 = 231 ml/kg (RSE 14.7%; bootstrap median 237, 158-302)
    lq2  <- log(0.0198); label("Slow inter-compartmental clearance Q3 per kilogram body weight (L/min/kg)") # Wang 2026 Table 3: Q3 = 19.8 ml/min/kg (RSE 10.8%; bootstrap median 19.8, 16.0-24.7)
    lvp2 <- log(1.36);   label("Slow peripheral volume of distribution V3 per kilogram body weight (L/kg)") # Wang 2026 Table 3: V3 = 1360 ml/kg (RSE 12.4%; bootstrap median 1360, 1075-1699)

    # Blood urea nitrogen effect on the central volume, as a power of the
    # ratio to the cohort median (Wang 2026 Methods Eq. 5 continuous-covariate
    # power model). The normalising median of 5 mmol/L is not printed in the
    # paper and is recovered exactly from the two V1 values that are -- see
    # covariateData$BUN$notes for the algebra -- and is corroborated by the
    # Table 4 simulation levels of 3, 5 and 7 mmol/L.
    e_bun_vc <- -0.821;  label("Exponent of (BUN / 5 mmol/L) on the central volume V1 (unitless)")       # Wang 2026 Table 3: 'Covariate thetaBUN on V1' = -0.821 (RSE 28.6%; bootstrap median -0.834, -1.357 to -0.242)

    # Inter-individual variability (log-normal eta on the log-scale
    # parameters). Wang 2026 Table 3 reports IIV as a percentage under the
    # heading 'Inter-individual variability'; the same percent convention is
    # used for the residual-error rows in that table, where a percentage can
    # only be the standard deviation of the error term. The printed values
    # are therefore read as omega (the log-scale SD) expressed in percent,
    # and the variances below are omega^2. At these magnitudes the
    # alternative reading (percent = approximate CV of the log-normal) moves
    # omega by less than 2% relative, so the choice is not load-bearing.
    #
    # IIV was retained only on CL, V1 and Q3: the Results text states that the
    # IIV estimates for V2, V3 and Q2 were close to 0 and those parameters
    # were 'not included with IIV in the final model'.
    #
    # CL and V1 carry a correlated block. Table 3 row 'Corr_CL & V1' = -0.821,
    # confirmed by the Results text 'a decrease in OFV by 15.039, with a
    # correlation coefficient of -0.821 for the IIV of CL and V1'. The
    # covariance is r * omega_CL * omega_V1 = -0.821 * 0.263 * 0.351.
    etalcl + etalvc ~ c(0.069169, -0.075789, 0.123201)   # Wang 2026 Table 3: eta(CL) = 26.3% -> 0.263^2 = 0.069169; Corr_CL & V1 = -0.821 -> cov = -0.821 * 0.263 * 0.351 = -0.075789; eta(V1) = 35.1% -> 0.351^2 = 0.123201
    etalq2 ~ 0.146689                                     # Wang 2026 Table 3: eta(Q3) = 38.3% -> 0.383^2 = 0.146689 (RSE 15.5%; bootstrap median 36.6%, 25.6-49.0%; shrinkage 7.8%)

    # Combined proportional plus additive residual error (Wang 2026 Methods
    # Eq. 4 'combined additive and proportional model'; Table 3 'Residual
    # variability' rows). The additive term is printed under a '%' column
    # header carried down from the proportional row above it, but an additive
    # residual must carry concentration units: 4.83 in the paper's ug/L
    # reporting scale, which is 0.00483 mg/L here and sits essentially at the
    # assay LLOQ of 5 ng/mL (= 5 ug/L), the magnitude an additive term takes
    # when it is set by assay noise at the bottom of the calibration range.
    propSd <- 0.135;    label("Proportional residual error (fraction)")      # Wang 2026 Table 3: epsilon_prop = 13.5% (RSE 9.6%; bootstrap median 13.3%, 10.5-15.8%; shrinkage 11.2%)
    addSd  <- 0.00483;  label("Additive residual error (mg/L)")              # Wang 2026 Table 3: epsilon_add = 4.83 ug/L (RSE 28.8%; bootstrap median 4.63, 1.94-7.18) -> 0.00483 mg/L
  })

  model({
    # Individual PK parameters. Every disposition parameter is reported per
    # kilogram of body weight (Wang 2026 Table 3), so body weight enters each
    # one linearly with no reference weight; standard allometric scaling and
    # age-dependent maturation on clearance were tested and not retained
    # (Methods, 'Model selection and covariate analysis'). BUN acts on the
    # central volume only (Results, 'Population covariant analysis'); CL and
    # V3 were independent of every covariate tested.
    cl  <- exp(lcl  + etalcl) * WT
    vc  <- exp(lvc  + etalvc) * WT * (BUN / 5)^e_bun_vc
    q   <- exp(lq)            * WT
    vp  <- exp(lvp)           * WT
    q2  <- exp(lq2  + etalq2) * WT
    vp2 <- exp(lvp2)          * WT

    # Three-compartment disposition with first-order elimination from the
    # central compartment (Wang 2026 Results, 'Population PK model': a
    # three-compartment model was significantly better than a two-compartment
    # model, dOFV = 65.8, p < 0.001). Ciprofol is given intravenously, so the
    # dose enters `central` directly via the cmt column of the user data set;
    # the study administered it as a 0.6 mg/kg bolus over 30 s, i.e. a short
    # zero-order infusion.
    d/dt(central)     <-  q  / vp  * peripheral1 + q2 / vp2 * peripheral2 -
                          (cl + q + q2) / vc * central
    d/dt(peripheral1) <-  q  / vc  * central     - q  / vp  * peripheral1
    d/dt(peripheral2) <-  q2 / vc  * central     - q2 / vp2 * peripheral2

    # Ciprofol plasma concentration in the central compartment. Dose units mg,
    # vc units L -> Cc units mg/L (= ug/mL). The paper reports concentrations
    # in ug/L, which is 1000 * Cc.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
