Muhamad_2025_cholecalciferol_pbpk <- function() {
  description <- paste(
    "PBPK (lumped whole-body, nlmixr2/rxode2). Serum 25-hydroxyvitamin D3",
    "in 77 healthy Cape Town schoolchildren aged 6-11 years given 250 ug",
    "(10,000 IU) oral vitamin D3 weekly for 3 years (Muhamad 2025, final",
    "Model 9d). Vitamin D3 moves from a gastrointestinal depot into the",
    "liver, where a fixed fraction (1/3) of hepatic elimination forms",
    "25(OH)D3; the remaining vitamin D3 distributes into venous blood and",
    "a lumped rest-of-body pool that is held in quasi-steady state with",
    "venous blood. 25(OH)D3 is carried as five compartments - venous",
    "blood, arterial blood, liver, fat mass and lean mass - and is",
    "eliminated from the liver by a sigmoidal, concentration-dependent",
    "clearance that rises with 25(OH)D3 and so produces the less than",
    "dose-proportional exposure that vitamin D is known for. A per-child",
    "constant input ENDOG, representing combined cutaneous synthesis and",
    "dietary intake, is back-calculated so the whole system starts exactly",
    "at that child's measured baseline 25(OH)D3. Compartment volumes scale",
    "linearly with body weight (30 kg reference) and the maximum clearance",
    "rate constant allometrically with an exponent of 0.75; the fat-mass",
    "and lean-mass split is driven by the BMI-for-age Z-score. Only the",
    "maximum clearance rate constant and the fat-mass partition",
    "coefficient were estimated - every other value was fixed from the",
    "authors' earlier healthy-adult model. Combined additive and",
    "proportional residual error; interindividual variability on the",
    "maximum clearance rate constant only."
  )
  reference <- paste(
    "Muhamad N, Walker N, Middelkoop K, Ganmaa D, Martineau AR, You T.",
    "Physiologically Based Pharmacokinetic Modelling of Serum",
    "25-Hydroxyvitamin D Concentrations in Schoolchildren Receiving Weekly",
    "Oral Vitamin D3 Supplementation. Nutrients. 2025;17(19):3028.",
    "doi:10.3390/nu17193028.",
    sep = " "
  )
  vignette <- "Muhamad_2025_cholecalciferol_pbpk"

  # Muhamad 2025 works in hours throughout (ka and CLmax in h^-1, blood
  # flows in L/h, the weight-interpolation denominator 24*365*3 h). Both
  # analytes are carried as nmol so that the 1/3 metabolic stoichiometry
  # and the nmol/L assay scale line up; an oral dose stated in ug must be
  # divided by the cholecalciferol molar mass before use (see the vignette).
  units <- list(
    time          = "h",
    dosing        = "nmol",
    concentration = "nmol/L"
  )

  # Issue #482. Vitamin D3 itself was never assayed in the ViDiKids Cape
  # Town children (Muhamad 2025 Section 4.3, fourth limitation), so its
  # three states are mechanistic and carry no measured matrix of their own;
  # the assayed quantity is serum 25(OH)D3, which the model reads off the
  # venous 25(OH)D3 compartment.
  compartmentData <- list(
    depot = list(
      analyte = "cholecalciferol", units = "nmol",
      specimen = "administration site", verified = TRUE
    ),
    venous = list(
      analyte = "cholecalciferol", units = "nmol",
      specimen = "whole blood", verified = TRUE
    ),
    liver = list(
      analyte = "cholecalciferol", units = "nmol",
      specimen = "tissue", verified = TRUE
    ),
    venous_25d3 = list(
      analyte = "25-hydroxyvitamin D3", units = "nmol",
      specimen = "serum", verified = TRUE
    ),
    arterial_25d3 = list(
      analyte = "25-hydroxyvitamin D3", units = "nmol",
      specimen = "whole blood", verified = TRUE
    ),
    liver_25d3 = list(
      analyte = "25-hydroxyvitamin D3", units = "nmol",
      specimen = "tissue", verified = TRUE
    ),
    adipose_25d3 = list(
      analyte = "25-hydroxyvitamin D3", units = "nmol",
      specimen = "tissue", verified = TRUE
    ),
    other_25d3 = list(
      analyte = "25-hydroxyvitamin D3", units = "nmol",
      specimen = "tissue", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Muhamad 2025 Section 2.3 and Table S5 (models 7-10)",
        "interpolate weight linearly between the baseline and the 3-year",
        "measurement, WTcur = WT + (WT3Y - WT) * time / (24*365*3), and use",
        "that current weight everywhere: all four lumped compartment volumes",
        "scale as WTcur/30 and CLmax as (WTcur/30)^0.75. The model file takes",
        "WTcur directly as the time-varying covariate column WT, so the",
        "caller supplies the growth trajectory in the event table (see the",
        "vignette); this is algebraically identical to the published",
        "interpolation and additionally lets the model be used with the",
        "sex-specific mean growth rates the paper applies to its simulated",
        "cohorts (boys 3.0 kg/year, girls 3.4 kg/year, Section 2.6; 3.16",
        "kg/year for the average child of Figure S16). Reference weight 30 kg",
        "is the lumping weight of Table S3, and equals the fitted cohort's",
        "mean baseline weight of 30.57 kg (Table S1)."
      ),
      source_name        = "WT"
    ),
    BMIZ = list(
      description        = "BMI-for-age Z-score",
      units              = "unitless (z-score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Computed per child in the WHO 2007 Z-score calculator (Muhamad 2025",
        "Section 2.7, reference [24]). Measured once, at baseline, and",
        "carried as time-constant. Drives the fraction of body weight that is",
        "fat mass through the Monasor-Ortola 2021 cubic (Muhamad 2025 Table",
        "S5 footnote, reference [20]), which in turn splits the non-",
        "eliminating body mass between the adipose_25d3 and other_25d3",
        "compartments. Source symbol ZBMI. Cohort distribution (Table S1):",
        "-3 to -2 1.30%, -2 to -1 9.09%, -1 to +1 63.64%, +1 to +2 18.18%,",
        "+2 to +3 7.79%."
      ),
      source_name        = "ZBMI"
    ),
    D25OH_BL = list(
      description        = "Baseline (pre-supplementation) serum 25-hydroxyvitamin D concentration",
      units              = "nmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per child. Plays the same dual role as the endogenous",
        "baseline in the authors' healthy-adult model (Muhamad 2025 Section",
        "3.1): it is the initial condition of every 25(OH)D3 compartment",
        "(each scaled by that compartment's volume and partition",
        "coefficient), and it fixes the per-child constant vitamin D3 input",
        "ENDOG, which is back-solved so that the sigmoidal clearance exactly",
        "balances production at baseline and the untreated system therefore",
        "holds flat. Measured by LC-MS/MS (Muhamad 2025 Section 2.1.1,",
        "reference [17]); fitted cohort mean 64.67 nmol/L, SD 14.85 (Table",
        "S1). Source symbol D25BASE."
      ),
      source_name        = "D25BASE"
    )
  )

  # Muhamad 2025 Section 2.6 compares the Cape Town and Ulaanbaatar cohorts
  # on sex, weight and compliance, and Section 2.5 samples sex when building
  # the European virtual population, but neither sex nor compliance enters
  # the final model as an estimated covariate effect: sex acts only through
  # the weight trajectory the caller supplies (see covariateData$WT), and
  # compliance is expressed by omitting dose records rather than by a
  # parameter. They are recorded here so the paper's covariate handling is
  # not lost, without being flagged as declared-but-unreferenced.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Section 2.6 fits baseline weight against age to infer sex-specific",
        "mean weight gain (boys 3.0 kg/year, girls 3.4 kg/year) used when",
        "simulating the Mongolian and European cohorts, and Table S5 model 5",
        "tested a sex-dependent weight trajectory - which fitted worse than",
        "model 4 (OBJFV 1965 vs 1825). In the final model (9d) sex enters",
        "only through whichever weight trajectory the caller supplies, never",
        "as a coefficient. Cohort: 30 of 77 male, 38.96% (Table S1)."
      )
    ),
    COMPLIANCE = list(
      description = "Received doses divided by number of weeks",
      units       = "fraction",
      type        = "continuous",
      notes       = paste(
        "Defined in Section 2.6 and Figure S13. Muhamad 2025 fits and",
        "simulates each child from that child's actual dosing record",
        "(Section 2.1.1: 'Dosing records for each participant were used for",
        "model fitting and simulation'), so compliance is represented by the",
        "dose records present in the event table and never as a model",
        "parameter. Section 3.5 reports that the 3-year change in serum",
        "25(OH)D did not correlate strongly with compliance."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 77,
    n_studies      = 1,
    age_range      = "6-11 years",
    age_median     = "8.94 years (mean; SD 1.29)",
    weight_median  = "30.57 kg (mean; SD 8.99)",
    height_median  = "149 cm (mean; SD 10.91)",
    bmi_median     = "17.75 kg/m^2 (mean; SD 3.46)",
    sex_female_pct = 61.04,
    disease_state  = "healthy schoolchildren (no vitamin D-related disease)",
    dose_range     = "250 ug (10,000 IU) oral vitamin D3 once weekly for 3 years",
    regions        = "Klipfontein district, Cape Town, South Africa",
    notes          = paste(
      "Baseline characteristics from Muhamad 2025 Table S1 (n = 77, the",
      "children with serum 25(OH)D at baseline and at 1, 2 and 3 years, who",
      "were used to fit the model). Baseline serum 25(OH)D 64.67 nmol/L, SD",
      "14.85. Drawn from the intervention arm of the ViDiKids Cape Town",
      "trial safety sub-study (NCT02880982); baseline samples were collected",
      "March-September 2017. The published model was additionally tested,",
      "without refitting, against 463 further Cape Town children from the",
      "same trial (Section 3.3) and against 1756 Mongolian schoolchildren",
      "given 350 ug weekly in the sister ViDiKids Ulaanbaatar trial",
      "(NCT02276755, Section 3.5), whose 25(OH)D rise the model",
      "overpredicted."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # ESTIMATED (Muhamad 2025 Table 2, "Parameter estimates for the
    # final model" = Table S5/S6 model 9d). Both are reported on the
    # natural-log scale, which is the scale used here verbatim.
    # ---------------------------------------------------------------
    lclmax <- -4.43
    label("Maximum 25(OH)D3 clearance rate constant at 30 kg (log h^-1)")  # Muhamad 2025 Table 2: est. in natural log -4.43, SE 0.0787, %RSE 1.78; linear 0.0119 (95% CI 0.0102, 0.0139) h^-1
    lkp_adipose_25d3 <- 1.54
    label("25(OH)D3 fat-mass:venous-blood partition coefficient Kp25fm (log, unitless)")  # Muhamad 2025 Table 2: est. in natural log 1.54, SE 0.136, %RSE 8.8; linear 4.66 (95% CI 3.6, 6.11)

    # ---------------------------------------------------------------
    # FIXED - vitamin D3 disposition. Muhamad 2025 Section 1 and 2.3:
    # "we changed the physiological parameters for children, fitted only
    # these two parameters, and fixed all other parameters from the
    # adult model [13]". Values are tabulated in Table S4 of this paper,
    # so they are sourced here and not from the upstream publication.
    # ---------------------------------------------------------------
    lka <- fixed(log(0.19))
    label("Vitamin D3 first-order absorption rate constant (log h^-1)")  # Muhamad 2025 Table S4: Ka = 0.19 h^-1
    lmclh <- fixed(log(0.222))
    label("Vitamin D3 hepatic clearance after repeated dosing (log L/h)")  # Muhamad 2025 Table S4: MCLH = 0.222 (see vignette Errata on the unit label)
    lkp_liver <- fixed(log(1))
    label("Vitamin D3 liver:venous-blood partition coefficient Kpl (log, unitless)")  # Muhamad 2025 Table S4: Kpl = 1
    lkp_other <- fixed(log(0.09))
    label("Vitamin D3 rest-of-body:venous-blood partition coefficient Kprb (log, unitless)")  # Muhamad 2025 Table S4: Kprb = 0.09

    # ---------------------------------------------------------------
    # FIXED - 25(OH)D3 disposition and the sigmoidal clearance shape.
    # Section 4.4: "we made a practical choice of fitting CLmax while
    # fixing C50 and gamma at the adult values", because separating all
    # three needs data at more than one dose level.
    # ---------------------------------------------------------------
    lkp_liver_25d3 <- fixed(log(1))
    label("25(OH)D3 liver:venous-blood partition coefficient Kp25l (log, unitless)")  # Muhamad 2025 Table S4: Kp25l = 1
    lkp_other_25d3 <- fixed(log(4))
    label("25(OH)D3 lean-mass:venous-blood partition coefficient Kp25lm (log, unitless)")  # Muhamad 2025 Section 3.1 and Table S6 model 9d: Kp25lm FIXED = 4, the optimum of the Figure S7 objective-function scan over 1, 2, 3, 4, 6, 8
    lc50 <- fixed(log(86.3))
    label("25(OH)D3 concentration at half-maximum clearance rate C50 (log nmol/L)")  # Muhamad 2025 Table S4: C50 = 86.3 nmol/L; Table S6 model 9d lists C50 as FIXED
    lhill <- fixed(log(5.64))
    label("Hill exponent of the sigmoidal 25(OH)D3 clearance (log, unitless)")  # Muhamad 2025 Table S4: gamma = 5.64; Table S6 model 9d lists gamma as FIXED
    fm_25d3 <- fixed(1 / 3)
    label("Fraction of hepatic vitamin D3 elimination forming 25(OH)D3 (unitless)")  # Muhamad 2025 Section 2.3 states Fm assumed to be 0.33; the Figure S2/S3 equations use exactly 1/3 (formation term x 1/3, ENDOG term x 3), so 1/3 is used here to keep the baseline balance exact

    # ---------------------------------------------------------------
    # FIXED - allometry. Section 2.3: "the maximum clearance rate for
    # serum 25(OH)D (CLMAX) was proportional to WT^0.75"; the exponent
    # carries no standard error anywhere in the paper.
    # ---------------------------------------------------------------
    e_wt_clmax <- fixed(0.75)
    label("Allometric exponent of body weight on CLmax (unitless)")  # Muhamad 2025 Section 2.3 and Table S5 (models 3-11): CLmax x (WTcur/30)^0.75

    # ---------------------------------------------------------------
    # IIV. Table 1 gives eta on CLmax only for models 9-10; Table 2
    # footnote back-solves the variance from the reported %CV.
    # ---------------------------------------------------------------
    etalclmax ~ 0.332  # Muhamad 2025 Table 2 footnote gives IIV CV% as sqrt(e^omega - 1) x 100% with omega the VARIANCE of the random effects, hence omega CLmax = 0.332; check: sqrt(exp(0.332) - 1) = 62.75% vs the reported 62.8 %CV

    # ---------------------------------------------------------------
    # Residual error. Section 2.3: "All models assumed combined additive
    # and proportional residual errors."
    # ---------------------------------------------------------------
    propSd <- 0.109
    label("Proportional residual error (fraction)")  # Muhamad 2025 Table 2: proportional error 0.109
    addSd <- 0.00249
    label("Additive residual error (nmol/L)")  # Muhamad 2025 Table 2: additive error 0.00249 nmol/L, called negligible in Section 3.1
  })

  model({
    # =================================================================
    # 1. PAEDIATRIC PHYSIOLOGY, Muhamad 2025 Table S3 ("Physiological
    #    parameters of the PBPK model for paediatrics after compartments
    #    are lumped together"), reported for the 30 kg reference child.
    #    Blood flows are not weight-scaled anywhere in the paper; only
    #    volumes and CLmax are (Table S5, models 2-10).
    # =================================================================
    q_co    <- 273.00   # L/h  cardiac output; Section 2.2 averages boys 4.8 and girls 4.3 L/min at age 7-9 to 4.55 L/min
    q_liver <- 61.97    # L/h  hepatic blood flow (22.7% of cardiac output, Table S2)
    q_other <- 211.03   # L/h  rest-of-body blood flow; Table S3 note "Qrb = Qco - Ql"

    # Fat-mass / lean-mass split of the rest-of-body flow, needed by the
    # models 9-10 structure (Figure S3) but not tabulated as such. Qfm is
    # the adipose row of Table S2 (5.2% of cardiac output = 14.20 L/h) and
    # Qlm takes the remainder, which keeps the paper's mass balance
    # Qco = Ql + Qfm + Qlm required by the arterial 25(OH)D3 equation.
    # This split is immaterial: sweeping Qfm from 2 to 200 L/h moves the
    # simulated 3-year serum 25(OH)D3 by less than 0.1 nmol/L (0.06%),
    # because both tissues equilibrate far faster than 25(OH)D3 is
    # cleared. The vignette reproduces that sweep.
    q_adipose <- 14.20        # L/h  Muhamad 2025 Table S2, Adipose row
    q_lean    <- q_other - q_adipose

    # Reference compartment volumes (L) at 30 kg, Table S3. Table S3 note:
    # "Assume apparent density = 1 kg/L. Vrb = 30L - Vl - Vven - Vart."
    v_venous_ref   <- 1.80
    v_liver_ref    <- 0.77
    v_arterial_ref <- 0.60
    v_other_ref    <- 26.83

    # =================================================================
    # 2. INDIVIDUAL PARAMETERS
    # =================================================================
    ka              <- exp(lka)
    mclh            <- exp(lmclh)
    kp_liver        <- exp(lkp_liver)
    kp_other        <- exp(lkp_other)
    kp_liver_25d3   <- exp(lkp_liver_25d3)
    kp_adipose_25d3 <- exp(lkp_adipose_25d3)
    kp_other_25d3   <- exp(lkp_other_25d3)
    c50             <- exp(lc50)
    hill            <- exp(lhill)

    # Table S5 (models 2-10): every lumped volume is directly proportional
    # to the current weight, referenced to the 30 kg lumping weight.
    v_venous   <- v_venous_ref   * WT / 30
    v_liver    <- v_liver_ref    * WT / 30
    v_arterial <- v_arterial_ref * WT / 30

    # Table S5 (models 3-10): CLmax = exp(TCLmax + etaCLmax) * (WTcur/30)^0.75.
    clmax <- exp(lclmax + etalclmax) * (WT / 30)^e_wt_clmax

    # =================================================================
    # 3. BODY COMPOSITION from the BMI-for-age Z-score
    #    Muhamad 2025 Table S5 footnote, after Monasor-Ortola 2021 [20]:
    #      fFM = max((28.61 + 7.82*Z - 0.91*Z^2 + 0.03*Z^3)/100, 0.05)  Z <= 6.19
    #      fFM = 0.49                                                    Z >  6.19
    #    The published Vfm / Vlm are linear interpolations between the
    #    baseline and 3-year fat and lean masses. Because ZBMI is measured
    #    once and carried as time-constant, that interpolation collapses
    #    exactly onto evaluating the same expressions at the current
    #    weight, which is the form used here.
    #
    #    The supplement prints LM3Y with the BASELINE weight,
    #      LM3Y = WT*(1 - fFM(3Y)) - (1.80 + 0.77 + 0.6)*WT/30,
    #    where FM3Y directly above it uses WT3Y. Taken literally lean mass
    #    would be frozen at its baseline value while fat mass grew, and
    #    the body-composition identity would break. Using the terminal
    #    weight in both makes Vfm + Vlm + Vven + Vl + Vart equal the
    #    current weight exactly at every weight and every ZBMI, so the
    #    published expression is a typographical slip. See the vignette.
    # =================================================================
    ffm <- max((28.61 + 7.82 * BMIZ - 0.91 * BMIZ^2 + 0.03 * BMIZ^3) / 100, 0.05)
    if (BMIZ > 6.19) {
      ffm <- 0.49
    }
    v_adipose <- WT * ffm
    v_lean    <- WT * (1 - ffm) -
      (v_venous_ref + v_liver_ref + v_arterial_ref) * WT / 30

    # =================================================================
    # 4. ENDOGENOUS VITAMIN D3 INPUT
    #    Muhamad 2025 Figures S2/S3, first displayed equation. ENDOG is
    #    the combined rate of cutaneous synthesis and dietary intake,
    #    back-solved per child so that the sigmoidal 25(OH)D3 clearance
    #    exactly balances 25(OH)D3 formation at the measured baseline.
    #    The published factor 3 is 1/Fm.
    # =================================================================
    endog <- clmax * D25OH_BL^hill / (c50^hill + D25OH_BL^hill) *
      D25OH_BL / fm_25d3

    # =================================================================
    # 5. VITAMIN D3
    #    Section 2.3 holds arterial blood and the rest of the body in
    #    quasi-steady state with venous blood, so both are algebraic
    #    functions of the venous amount (paper Equations 1 and 2). The
    #    forms below are the supplement's, which are the dimensionally
    #    consistent ones; the main-text renderings of Equations 1 and 2
    #    are garbled (see the vignette Errata).
    # =================================================================
    cven_vd  <- venous / v_venous
    aart_vd  <- cven_vd * q_co * v_arterial / (q_liver + q_other)
    cart_vd  <- aart_vd / v_arterial
    aother_vd <- aart_vd * kp_other * v_other_ref * (WT / 30) / v_arterial
    cliver_vd <- liver / (v_liver * kp_liver)

    d/dt(depot)  <- -ka * depot
    d/dt(venous) <- q_liver * cliver_vd +
      q_other * aother_vd / (v_other_ref * (WT / 30) * kp_other) -
      q_co * cven_vd
    d/dt(liver)  <- -q_liver * cliver_vd + q_liver * aart_vd / v_arterial -
      mclh * cliver_vd + endog + ka * depot

    # =================================================================
    # 6. 25(OH)D3 - Muhamad 2025 Figure S3 (models 9 and 10). Vitamin D3
    #    metabolised in the liver enters as fm_25d3 of hepatic vitamin D3
    #    elimination; 25(OH)D3 leaves only through the sigmoidal hepatic
    #    clearance, which increases with 25(OH)D3 and so makes exposure
    #    less than proportional to dose (Section 4.4).
    # =================================================================
    c25_liver    <- liver_25d3    / (v_liver   * kp_liver_25d3)
    c25_adipose  <- adipose_25d3  / (v_adipose * kp_adipose_25d3)
    c25_lean     <- other_25d3    / (v_lean    * kp_other_25d3)
    c25_venous   <- venous_25d3   / v_venous
    c25_arterial <- arterial_25d3 / v_arterial

    clr_25d3 <- clmax * c25_liver^hill / (c50^hill + c25_liver^hill) * c25_liver

    d/dt(venous_25d3) <- q_liver * c25_liver + q_adipose * c25_adipose +
      q_lean * c25_lean - q_co * c25_venous
    d/dt(liver_25d3) <- -q_liver * c25_liver + q_liver * c25_arterial +
      mclh * cliver_vd * fm_25d3 - clr_25d3
    d/dt(adipose_25d3) <- q_adipose * c25_arterial - q_adipose * c25_adipose
    d/dt(other_25d3)   <- q_lean * c25_arterial - q_lean * c25_lean
    d/dt(arterial_25d3) <- q_co * c25_venous - q_liver * c25_arterial -
      (q_adipose + q_lean) * c25_arterial

    # =================================================================
    # 7. INITIAL CONDITIONS - Muhamad 2025 Figure S3, "Initial
    #    conditions". Every 25(OH)D3 compartment starts at the child's
    #    measured baseline concentration times its own volume and
    #    partition coefficient, and the two vitamin D3 states start at
    #    the steady state that ENDOG sustains.
    # =================================================================
    depot(0)         <- 0
    venous(0)        <- endog * v_venous / mclh
    liver(0)         <- endog * v_liver / mclh * kp_liver
    venous_25d3(0)   <- D25OH_BL * v_venous
    liver_25d3(0)    <- D25OH_BL * v_liver * kp_liver_25d3
    adipose_25d3(0)  <- D25OH_BL * v_adipose * kp_adipose_25d3
    other_25d3(0)    <- D25OH_BL * v_lean * kp_other_25d3
    arterial_25d3(0) <- D25OH_BL * v_arterial

    # =================================================================
    # 8. OBSERVATION - serum 25(OH)D3, the only quantity the trial
    #    measured (Section 2.1.1).
    # =================================================================
    Cc <- c25_venous
    Cc ~ add(addSd) + prop(propSd)
  })
}
