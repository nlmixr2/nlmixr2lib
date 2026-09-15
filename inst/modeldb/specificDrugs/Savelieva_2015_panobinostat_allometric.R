Savelieva_2015_panobinostat_allometric <- function() {
  description <- "Three-compartment population PK model with lagged first-order absorption and linear elimination for the pan-deacetylase inhibitor panobinostat in adults with advanced hematologic and solid tumors (Savelieva 2015, SECOND final model). Developed on the same pooled intravenous and oral data set as Savelieva_2015_panobinostat and reported alongside it, this re-parameterization replaces the body-surface-area effects with fixed allometric weight scaling on every clearance and volume, replaces the distribution rate constants with intercompartmental clearances and peripheral volumes that each carry their own interindividual variability and age effect, and adds a formulation-dependent absorption lag. The authors note that this model did not satisfy NONMEM's default convergence criterion, although it fits substantially better by AIC and BIC."
  reference <- "Savelieva M, Woo MM, Schran H, Mu S, Nedelman J, Capdeville R. Population pharmacokinetics of intravenous and oral panobinostat in patients with hematologic and solid tumors. Eur J Clin Pharmacol. 2015;71(6):663-672. doi:10.1007/s00228-015-1846-7. Parameter estimates from Supplementary Table S3b; model code from Supplementary Table S3a."
  vignette <- "Savelieva_2015_panobinostat"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "panobinostat", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "panobinostat", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "panobinostat", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "panobinostat", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on every clearance and every volume, centered on 70 kg (Savelieva 2015 Supplementary Table S3a: '(WT0/70)**THETA(10)' on CL, '(WT0/70)**THETA(11)' on V2, '(WT0/70)**0.75' on Q3 and Q4, '(WT0/70)' on V3 and V4). Note that 70 kg is a standard reference weight and NOT the cohort median, which is 76.4 kg (range 41-196.4 kg); Table 3 of the paper accordingly tabulates its typical-value predictions at 76.4 kg rather than at 70 kg. The exponents are held constant at the canonical allometric values, 0.75 for clearances and 1 for volumes: the Results state that 'all clearances were assumed proportional to weight^0.75 and all volumes to weight^1', and Supplementary Table S3b reports Thetas 10 and 11 as exactly 0.750 and 1.000 with no standard error, no percent standard error and no bootstrap interval. Weight replaces the body surface area used by the first final model.",
      source_name        = "WT0"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effects centered on 61 years, the population median, applied to ALL six disposition parameters in this model - CL, V2, Q3, V3, Q4 and V4 (Savelieva 2015 Supplementary Table S3a). Extending the age effect to the intercompartmental clearances and peripheral volumes is exactly what distinguishes model 4 from model 3 in Table 2. Observed range 16-88 years, quartiles 51 and 70 years. The V2 exponent is -0.005 with a percent standard error of 213 and a bootstrap interval spanning zero, i.e. indistinguishable from no effect; it is retained here because the source retained it.",
      source_name        = "AGE0"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Caucasian, when RACE_BLACK and RACE_OTHER are also 0)",
      notes              = "Multiplicative factor applied as THETA^indicator on CL and on V2 only, not on the peripheral parameters (Savelieva 2015 Supplementary Table S3a: 'IF (RACE .EQ. 3) AS=1' then 'THETA(14)**(AS)'). The canonical 1 = Asian orientation matches the source coding. The three race indicators are mutually exclusive; all three equal to 0 selects the Caucasian reference. 27 of the 581 patients were Asian (Table 1).",
      source_name        = "AS (derived from RACE == 3)"
    ),
    RACE_BLACK = list(
      description        = "Black race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Caucasian, when RACE_ASIAN and RACE_OTHER are also 0)",
      notes              = "Multiplicative factor applied as THETA^indicator on CL and on V2 (Savelieva 2015 Supplementary Table S3a: 'IF (RACE .EQ. 2) BL=1' then 'THETA(16)**(BL)'). 34 of the 581 patients were Black (Table 1). The V2 factor of 1.817 is markedly larger than the 1.241 of the first final model and is imprecise (percent standard error 37), which is one of several signs that the race effects on V2 are less well determined in this parameterization.",
      source_name        = "BL (derived from RACE == 2)"
    ),
    RACE_OTHER = list(
      description        = "Race category 'other' indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Caucasian, when RACE_ASIAN and RACE_BLACK are also 0)",
      notes              = "Multiplicative factor applied as THETA^indicator on CL and on V2 (Savelieva 2015 Supplementary Table S3a: 'IF (RACE .EQ. 88) OT=1' then 'THETA(18)**(OT)'). 24 of the 581 patients fell in this category (Table 1). As in the first final model this group has the lowest clearance factor, 0.665, and correspondingly the highest typical exposure of the four race categories in Table 3.",
      source_name        = "OT (derived from RACE == 88)"
    ),
    FORM_PANO_CSF = list(
      description        = "Panobinostat clinical service formulation versus final market image indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (final market image, FMI, the formulation intended for commercialization)",
      notes              = "Selects BOTH the absorption rate constant and the absorption lag time between two estimated values each (Savelieva 2015 Supplementary Table S3a: 'KA=THETA(7)*(1-FORM)+THETA(8)*FORM' and 'ALAG1=THETA(20)*(1-FORM)+THETA(21)*FORM'). 1 = clinical service formulation (CSF), used in studies B2101, B2102 and B1101; 0 = final market image (FMI), used in every other oral study. 106 of the 494 orally dosed patients received the CSF. The formulation-dependent lag is the single addition that takes model 1 to model 2 in Table 2 and is credited in the Results with the better absorption-phase fit visible in the predictive checks. Formulation does not affect bioavailability, which uses one THETA(9) for both formulations. This covariate is irrelevant to intravenous dosing, which bypasses the depot entirely.",
      source_name        = "FORM"
    )
  )

  # Screened by the authors but NOT retained in either final model, plus the
  # body-size covariate this model replaced. Documented so the provenance of the
  # covariate search survives, without raising a declared-but-unreferenced
  # convention warning.
  covariatesDataExcluded <- list(
    BSA = list(
      description = "Baseline body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Explicitly REMOVED when this second final model was built: the Results state that 'BSA was removed from the model, and all clearances were assumed proportional to weight^0.75 and all volumes to weight^1'. BSA IS the body-size covariate of the first final model - see Savelieva_2015_panobinostat. Median 1.9 m^2, quartiles 1.8 and 2.1 m^2; computed by the authors with the Gehan-George formula BSA = 234.94 * (Weight^0.515 * Height^0.422) / 10000."
    ),
    HT = list(
      description = "Baseline height",
      units       = "cm",
      type        = "continuous",
      notes       = "Collected at screening and used only as an input to the Gehan-George BSA formula, never entered as a covariate in its own right; median 170 cm, range 143-198 cm. Missing for 35 patients, who were assigned the population median."
    ),
    BMI = list(
      description = "Baseline body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened as a body-size covariate on CL and V2 and not retained (Savelieva 2015 Methods, Analysis of the effects of covariates). Missing for 35 patients, who were assigned the population median."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL and V2 and not retained; the paper reports that covariate analysis 'showed no impact on panobinostat clearance and volume by patients' sex'. 219 of 581 patients were female (Table 1)."
    ),
    CRCL = list(
      description = "Baseline creatinine clearance",
      units       = "L/h",
      type        = "continuous",
      notes       = "Screened on CL and V2 and not retained; kidney function had no statistically significant effect. Computed by the authors with the Cockcroft-Gault formula from serum creatinine in micromoles per litre. Missing for 13 patients, who were assigned the population median."
    ),
    HEPIMP_MILD = list(
      description = "Mild hepatic impairment indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Part of the four-level liver-status covariate (normal 483, mild 91, moderate 6, severe 1; Table 1). Screened on CL and V2 and not retained. The Discussion attributes the null result to trial eligibility criteria that generally required baseline bilirubin at or below 1.5x ULN and AST/ALT at or below 2x ULN, and notes that dedicated organ-impairment studies DID find a significant exposure increase."
    ),
    HEPIMP_MOD = list(
      description = "Moderate hepatic impairment indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Part of the same four-level liver-status covariate; only 6 of 581 patients. Not retained."
    ),
    HEPIMP_SEV = list(
      description = "Severe hepatic impairment indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Part of the same four-level liver-status covariate; a single patient of 581. Not retained."
    ),
    TUMTP_OTHER = list(
      description = "Tumor-type indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tumor type was screened on CL and V2 across the pooled hematologic and solid-tumor population and not retained. Placeholder entry standing for the paper's tumor-type covariate as a whole; the source does not publish the individual category codes it tested."
    ),
    CONMED_AZOLE = list(
      description = "Strong CYP3A4/5 inhibitor comedication indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "One of five comedication groups screened on CL and not retained: drugs known to prolong QT, CYP2D6 substrates, strong CYP3A4/5 inhibitors, CYP3A4 inducers, and sensitive CYP3A4 substrates. Placeholder entry standing for that comedication screen. The Discussion attributes the null result to protocol guidance that recommended avoiding strong CYP3A4 inhibitors, and notes that a dedicated ketoconazole interaction study DID show a significant exposure increase."
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 581,
    n_studies       = 14,
    n_observations  = 7834,
    age_range       = "16-88 years",
    age_median      = "61 years (quartiles 51 and 70 years)",
    weight_range    = "41-196.4 kg",
    weight_median   = "76.4 kg",
    height_range    = "143-198 cm",
    height_median   = "170 cm",
    bsa_median      = "1.9 m^2 (quartiles 1.8 and 2.1 m^2)",
    sex_female_pct  = 37.7,
    race_ethnicity  = c(Caucasian = 85.4, Black = 5.9, Asian = 4.6, Other = 4.1),
    disease_state   = "Advanced hematologic and solid tumors, including cutaneous T-cell lymphoma, chronic myeloid leukemia, multiple myeloma, Hodgkin lymphoma, non-Hodgkin lymphoma and advanced solid tumors",
    hepatic_function = "Liver status graded on total bilirubin and AST against the upper limit of normal: normal 483, mild 91, moderate 6, severe 1 (Table 1)",
    dose_range      = "Intravenous 1.2-20 mg/m^2 daily under various intermittent regimens (studies A2101 and A2102, 87 patients); oral 10-80 mg/day in the phase 1 dose-escalation studies and 20-45 mg in the phase 2 and clinical pharmacology studies, most commonly 20 mg on days 1, 3 and 5 of each week (494 patients)",
    formulation     = "Clinical service formulation (CSF) in oral studies B2101, B2102 and B1101 (106 patients); final market image (FMI) in every other oral study (388 patients)",
    regions         = "International; study B1101 enrolled 13 Japanese patients and B1201 was conducted in Japan",
    notes           = "Same pooled data set as Savelieva_2015_panobinostat: 14 open-label phase 1 and phase 2 studies listed in Supplementary Table S1. This model is row 4 of the Table 2 model-development sequence and was built in response to reviewer suggestions. Its objective function is 32591 against 33758 for the first final model, and both AIC and BIC improve, but the paper records that this model did not satisfy NONMEM's default convergence criterion, so its standard errors should be read with that caveat."
  )

  ini({
    # Structural parameters - typical values for the reference patient: a
    # Caucasian of body weight 70 kg and age 61 years, i.e. every covariate term
    # equal to 1. Note that the reference weight is the standard 70 kg and NOT
    # the cohort median of 76.4 kg, so these typical values are not directly
    # comparable to the Table 3 predictions, which are tabulated at 76.4 kg.
    # Values are from Savelieva 2015 Supplementary Table S3b.
    lcl  <- log(28.833);  label("Clearance (L/h)")                                         # Suppl Table S3b, Theta 1 'CL (L/h)' = 28.833, PctSE 5.72, bootstrap Q10-Q90 27.046-31.184
    lvc  <- log(30.862);  label("Central volume of distribution V2 (L)")                   # Suppl Table S3b, Theta 2 'V2 (L)' = 30.862, PctSE 9.39, bootstrap Q10-Q90 26.869-33.765
    lq   <- log(32.751);  label("Intercompartmental clearance to peripheral1, Q3 (L/h)")   # Suppl Table S3b, Theta 3 'Q3 (L/h)' = 32.751, PctSE 9.46, bootstrap Q10-Q90 28.644-36.585
    lvp  <- log(71.874);  label("Peripheral1 volume of distribution V3 (L)")               # Suppl Table S3b, Theta 4 'V3 (L)' = 71.874, PctSE 12.69, bootstrap Q10-Q90 61.651-83.481
    lq2  <- log(31.088);  label("Intercompartmental clearance to peripheral2, Q4 (L/h)")   # Suppl Table S3b, Theta 5 'Q4 (L/h)' = 31.088, PctSE 7.23, bootstrap Q10-Q90 28.027-33.405
    lvp2 <- log(803.193); label("Peripheral2 volume of distribution V4 (L)")               # Suppl Table S3b, Theta 6 'V4 (L)' = 803.193, PctSE 6.17, bootstrap Q10-Q90 728.986-861.162

    # Absorption. Two separately estimated rate constants and two separately
    # estimated lag times, each selected by formulation.
    lka_fmi   <- log(0.420); label("Absorption rate constant, final market image formulation (1/h)")  # Suppl Table S3b, Theta 7 'KA.FMI (1/h)' = 0.420, PctSE 7.65, bootstrap Q10-Q90 0.375-0.457
    lka_csf   <- log(0.631); label("Absorption rate constant, clinical service formulation (1/h)")    # Suppl Table S3b, Theta 8 'KA.CSF (1/h)' = 0.631, PctSE 6.02, bootstrap Q10-Q90 0.584-0.682
    ltlag_fmi <- log(0.162); label("Absorption lag time, final market image formulation (h)")         # Suppl Table S3b, Theta 20 'LAG.FMI (h)' = 0.162, PctSE 5.81, bootstrap Q10-Q90 0.158-0.164
    ltlag_csf <- log(0.296); label("Absorption lag time, clinical service formulation (h)")           # Suppl Table S3b, Theta 21 'LAG.CSF (h)' = 0.296, PctSE 9.25, bootstrap Q10-Q90 0.285-0.353

    # Absolute oral bioavailability, shared by both oral formulations.
    lfdepot <- log(0.194); label("Absolute oral bioavailability (fraction)")  # Suppl Table S3b, Theta 9 'F1' = 0.194, PctSE 5.62, bootstrap Q10-Q90 0.181-0.209

    # Allometric exponents, held constant at the canonical values rather than
    # estimated. Supplementary Table S3b reports Thetas 10 and 11 as exactly
    # 0.750 and 1.000 with no standard error, no percent standard error and no
    # bootstrap interval, and the Results state that all clearances were assumed
    # proportional to weight^0.75 and all volumes to weight^1. The control stream
    # writes the same two numbers as literal constants on Q3, V3, Q4 and V4, so
    # these two parameters carry the exponents for every clearance and every
    # volume in the model.
    e_wt_cl <- fixed(0.750); label("Allometric exponent on weight relative to 70 kg for all clearances (unitless)")  # Suppl Table S3b, Theta 10 'CL.WT' = 0.750, SE and bootstrap reported as NA; same value hardcoded on Q3 and Q4 in Suppl Table S3a
    e_wt_vc <- fixed(1.000); label("Allometric exponent on weight relative to 70 kg for all volumes (unitless)")     # Suppl Table S3b, Theta 11 'V2.WT' = 1.000, SE and bootstrap reported as NA; same value hardcoded on V3 and V4 in Suppl Table S3a

    # Age effects, one per disposition parameter.
    e_age_cl  <- 0.137;  label("Power exponent on age relative to 61 years for clearance (unitless)")                          # Suppl Table S3b, Theta 12 'CL.AGE' = 0.137, PctSE 63.06, bootstrap Q10-Q90 0.018-0.243
    e_age_vc  <- -0.005; label("Power exponent on age relative to 61 years for central volume (unitless)")                     # Suppl Table S3b, Theta 13 'V2.AGE' = -0.005, PctSE 213.48, bootstrap Q10-Q90 -0.146 to 0.410
    e_age_q   <- 0.410;  label("Power exponent on age relative to 61 years for intercompartmental clearance Q3 (unitless)")    # Suppl Table S3b, Theta 22 'Q3.AGE' = 0.410, PctSE 59.43, bootstrap Q10-Q90 0.098-0.742
    e_age_vp  <- 0.713;  label("Power exponent on age relative to 61 years for peripheral1 volume V3 (unitless)")              # Suppl Table S3b, Theta 23 'V3.AGE' = 0.713, PctSE 43.37, bootstrap Q10-Q90 0.395-1.145
    e_age_q2  <- 0.212;  label("Power exponent on age relative to 61 years for intercompartmental clearance Q4 (unitless)")    # Suppl Table S3b, Theta 24 'Q4.AGE' = 0.212, PctSE 78.25, bootstrap Q10-Q90 0.008-0.345
    e_age_vp2 <- 0.530;  label("Power exponent on age relative to 61 years for peripheral2 volume V4 (unitless)")              # Suppl Table S3b, Theta 25 'V4.AGE' = 0.530, PctSE 28.50, bootstrap Q10-Q90 0.331-0.698

    # Race factors, on CL and V2 only.
    e_race_asian_cl <- 1.203; label("Multiplicative factor on clearance for Asian relative to Caucasian patients (unitless)")        # Suppl Table S3b, Theta 14 'CL.ASIAN' = 1.203, PctSE 8.44, bootstrap Q10-Q90 1.087-1.333
    e_race_asian_vc <- 2.060; label("Multiplicative factor on central volume for Asian relative to Caucasian patients (unitless)")   # Suppl Table S3b, Theta 15 'V2.ASIAN' = 2.060, PctSE 30.75, bootstrap Q10-Q90 1.253-2.817
    e_race_black_cl <- 0.941; label("Multiplicative factor on clearance for Black relative to Caucasian patients (unitless)")        # Suppl Table S3b, Theta 16 'CL.BLACK' = 0.941, PctSE 14.31, bootstrap Q10-Q90 0.802-1.125
    e_race_black_vc <- 1.817; label("Multiplicative factor on central volume for Black relative to Caucasian patients (unitless)")   # Suppl Table S3b, Theta 17 'V2.BLACK' = 1.817, PctSE 36.95, bootstrap Q10-Q90 1.080-2.928
    e_race_other_cl <- 0.665; label("Multiplicative factor on clearance for race category 'other' relative to Caucasian (unitless)")      # Suppl Table S3b, Theta 18 'CL.OTHER' = 0.665, PctSE 20.98, bootstrap Q10-Q90 0.513-0.877
    e_race_other_vc <- 0.835; label("Multiplicative factor on central volume for race category 'other' relative to Caucasian (unitless)") # Suppl Table S3b, Theta 19 'V2.OTHER' = 0.835, PctSE 30.76, bootstrap Q10-Q90 0.609-1.300

    # Interindividual variability. Exponential parameterization throughout
    # (Savelieva 2015 Methods), so the tabulated OMEGA entries are variances on
    # the log scale. This model adds multiplicative random effects to the four
    # peripheral parameters, arranged as three 2x2 blocks: CL with V2, Q3 with
    # V3, and Q4 with V4. The V2 variance of 1.668 is far larger than the 0.334
    # of the first final model, which fits the paper's note that this model did
    # not meet NONMEM's convergence criterion.
    etalcl + etalvc ~ c(0.407,
                        0.151, 1.668)  # Suppl Table S3b, Omega 1,1 'OM.CL' = 0.407; 2,1 'OM.CLV2' = 0.151; 2,2 'OM.V2' = 1.668
    etalq + etalvp ~ c(0.666,
                       0.505, 0.441)   # Suppl Table S3b, Omega 3,3 'OM.Q3' = 0.666; 4,3 'OM.Q3V3' = 0.505; 4,4 'OM.V3' = 0.441
    etalq2 + etalvp2 ~ c(0.497,
                         0.556, 0.700) # Suppl Table S3b, Omega 5,5 'OM.Q4' = 0.497; 6,5 'OM.Q4V4' = 0.556; 6,6 'OM.V4' = 0.700

    # Residual error. Suppl Table S3a: 'Y=F*(1+EPS(1))+EPS(2)', a combined
    # proportional plus additive model, with the SIGMA entries reported as
    # variances in Suppl Table S3b.
    propSd <- 0.42426; label("Proportional residual error (fraction)")   # Suppl Table S3b, Sigma 1 'VAR.PROP' = 0.180, PctSE 4.52; sqrt(0.180)
    addSd  <- 0.10488; label("Additive residual error (ng/mL)")          # Suppl Table S3b, Sigma 2 'VAR.ADD' = 0.011, PctSE 46.22; sqrt(0.011)
  })

  model({
    # Covariate multipliers, transcribed from Savelieva 2015 Supplementary
    # Table S3a:
    #   CL1 = THETA(1)*EXP(ETA(1))*(WT0/70)**THETA(10)
    #   CL  = CL1*(AGE0/61)**THETA(12)*THETA(14)**(AS)*THETA(16)**(BL)*THETA(18)**(OT)
    #   V21 = THETA(2)*EXP(ETA(2))*(WT0/70)**THETA(11)
    #   V2  = V21*(AGE0/61)**THETA(13)*THETA(15)**(AS)*THETA(17)**(BL)*THETA(19)**(OT)
    #   Q3  = THETA(3)*(WT0/70)**0.75*EXP(ETA(3))*(AGE0/61)**THETA(22)
    #   V3  = THETA(4)*(WT0/70)     *EXP(ETA(4))*(AGE0/61)**THETA(23)
    #   Q4  = THETA(5)*(WT0/70)**0.75*EXP(ETA(5))*(AGE0/61)**THETA(24)
    #   V4  = THETA(6)*(WT0/70)     *EXP(ETA(6))*(AGE0/61)**THETA(25)
    # The literal 0.75 and 1 on the peripheral parameters are the same two
    # allometric exponents carried by e_wt_cl and e_wt_vc.
    wt_cl   <- (WT / 70)^e_wt_cl
    wt_vc   <- (WT / 70)^e_wt_vc
    race_cl <- e_race_asian_cl^RACE_ASIAN * e_race_black_cl^RACE_BLACK * e_race_other_cl^RACE_OTHER
    race_vc <- e_race_asian_vc^RACE_ASIAN * e_race_black_vc^RACE_BLACK * e_race_other_vc^RACE_OTHER

    cl  <- exp(lcl  + etalcl)  * wt_cl * (AGE / 61)^e_age_cl * race_cl
    vc  <- exp(lvc  + etalvc)  * wt_vc * (AGE / 61)^e_age_vc * race_vc
    q   <- exp(lq   + etalq)   * wt_cl * (AGE / 61)^e_age_q
    vp  <- exp(lvp  + etalvp)  * wt_vc * (AGE / 61)^e_age_vp
    q2  <- exp(lq2  + etalq2)  * wt_cl * (AGE / 61)^e_age_q2
    vp2 <- exp(lvp2 + etalvp2) * wt_vc * (AGE / 61)^e_age_vp2

    # Formulation-dependent absorption, written as the source's linear switches
    # 'KA = THETA(7)*(1-FORM) + THETA(8)*FORM' and
    # 'ALAG1 = THETA(20)*(1-FORM) + THETA(21)*FORM'. The source additionally
    # forces KA = 0, ALAG1 = 0 and F1 = 1 for intravenous records; here that is
    # achieved instead by routing intravenous doses straight into central, which
    # never touches the depot, ka, the lag or f(depot).
    ka     <- exp(lka_fmi)   * (1 - FORM_PANO_CSF) + exp(lka_csf)   * FORM_PANO_CSF
    tlag   <- exp(ltlag_fmi) * (1 - FORM_PANO_CSF) + exp(ltlag_csf) * FORM_PANO_CSF
    fdepot <- exp(lfdepot)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    f(depot)    <- fdepot
    alag(depot) <- tlag

    # Doses are in mg and volumes in L, so central/vc is mg/L; the factor 1000
    # converts to the ng/mL scale the paper reports and on which the additive
    # residual error is expressed. This reproduces the source's scaling
    # 'S2=V2/1000' (Suppl Table S3a).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
