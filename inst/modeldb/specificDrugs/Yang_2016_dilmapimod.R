Yang_2016_dilmapimod <- function() {
  description <- "Three-compartment IV population PK model for dilmapimod (SB-681323, a p38 MAPK inhibitor) coupled with an empirical indirect-response model for the inflammatory biomarker C-reactive protein (CRP) in severe-trauma adults at risk for acute respiratory distress syndrome (Yang 2016). BMI is a power covariate on CL and Q2. No statistically significant dilmapimod effect on CRP was retained in the final PD model, so the CRP component is an empirical post-injury production-decline / first-order-loss profile that is decoupled from dilmapimod exposure (Yang 2016 Results section 3.3.1)."
  reference   <- "Yang S, Pene Dumitrescu T. Population pharmacokinetics and pharmacodynamics modelling of dilmapimod in severe trauma subjects at risk for acute respiratory distress syndrome. Drugs R D. 2017;17(1):145-156. doi:10.1007/s40268-016-0161-9"
  vignette    <- "Yang_2016_dilmapimod"
  units       <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    BMI = list(
      description        = "Body mass index at baseline.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on CL and Q2 with reference value 27.4 kg/m^2 (population mean per Yang 2016 Results section 3.2). Applied as (BMI/27.4)^exponent. Time-fixed at baseline.",
      source_name        = "BMI"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 57L,
    n_studies      = 1L,
    age_range      = "median 39 years (adult trauma population)",
    age_median     = "39 years",
    weight_range   = "median 86 kg",
    weight_median  = "86 kg",
    sex_female_pct = NA_real_,
    disease_state  = "Severe major-trauma (non-head-injury) adults at risk for acute respiratory distress syndrome; median Injury Severity Score 25, median Glasgow Coma Score 15.",
    dose_range     = "Cohort 1: 3 mg IV over 4 h x 3 daily doses; cohort 2: 7.5 mg IV over 24 h x 3 daily doses; cohort 3: 7.5 mg IV over 4 h x 3 daily doses; cohort 4: 10 mg IV over 24 h x 3 daily doses (clinicaltrials.gov NCT00996840).",
    regions        = "United States (six sites)",
    bmi_median     = "27 kg/m^2 (population mean 27.4 used as the BMI covariate reference)",
    time_since_trauma_to_first_dose = "9-28 h (median 22 h); only 4 of 57 subjects were dosed within 12 h of injury, 41 between 12-24 h, and 28 (with overlap) between 24-36 h.",
    notes          = "Phase IIa randomised, double-blind, placebo-controlled, parallel-group study at six US sites (NCT00996840). 57 active-treatment subjects supplied 471 dilmapimod concentration records (40, or 8.5%, were excluded as outliers, missing-infusion, or implausible pre-dose). The combined CRP PK/PD analysis comprised 73 subjects (53 active + 20 placebo) and 651 CRP records. Baseline demographics from Yang 2016 section 3.1 and the cited primary clinical paper Christie et al. (reference [10] in Yang 2016)."
  )

  ini({
    # ===== Population PK (Yang 2016 Table 2; final 3-compartment IV model) =====
    # Each typical value is reported in the paper both as exp(theta_pop,k) and as
    # the underlying theta_pop,k on the natural-log scale.
    lcl    <- log(35.9)   ; label("Clearance CL (L/h)")                                                        # Yang 2016 Table 2: exp(theta_pop,1) = 35.9, theta_pop,1 = 3.58
    lvc    <- log(8.1)    ; label("Central volume of distribution Vc (V1, L)")                                 # Yang 2016 Table 2: exp(theta_pop,2) = 8.1,  theta_pop,2 = 2.09
    lq     <- log(28.2)   ; label("Inter-compartmental clearance to peripheral1 (Q1, L/h)")                    # Yang 2016 Table 2: exp(theta_pop,3) = 28.2, theta_pop,3 = 3.34
    lvp    <- log(35.9)   ; label("Peripheral volume of distribution peripheral1 (V2, L)")                     # Yang 2016 Table 2: exp(theta_pop,4) = 35.9, theta_pop,4 = 3.58
    lq2    <- log(5.7)    ; label("Inter-compartmental clearance to peripheral2 (Q2, L/h)")                    # Yang 2016 Table 2: exp(theta_pop,5) = 5.7,  theta_pop,5 = 1.74
    lvp2   <- log(115.6)  ; label("Peripheral volume of distribution peripheral2 (V3, L)")                     # Yang 2016 Table 2: exp(theta_pop,6) = 115.6, theta_pop,6 = 4.75

    # BMI covariate (power model centred at 27.4 kg/m^2 per Yang 2016 Results section 3.2)
    e_bmi_cl  <- 1.36     ; label("Power exponent of (BMI/27.4) on CL (unitless)")                              # Yang 2016 Table 2: BMI covariate on CL = 1.36 (95% CI 0.866-1.85)
    e_bmi_q2  <- 2.42     ; label("Power exponent of (BMI/27.4) on Q2 (unitless)")                              # Yang 2016 Table 2: BMI covariate on Q2 = 2.42 (95% CI 1.40-3.36)

    # Inter-individual variability (log-normal; final model retained IIV only on CL and Q2;
    # Yang 2016 Results section 3.2: data did not support adding IIV on other PK parameters).
    # NONMEM convention: omega reported as variance, eta ~ N(0, variance); approximate CV%
    # in Table 2 footnote c is computed as 100 * sqrt(variance).
    etalcl ~ 0.0991   # Yang 2016 Table 2: IIV CL  variance 0.0991, CV 31.5%
    etalq2 ~ 0.226    # Yang 2016 Table 2: IIV Q2  variance 0.226,  CV 47.5%

    # Residual error on dilmapimod plasma concentration. Yang 2016 Discussion section 3.2:
    # "The proportional error model was adequate". Sigma in Table 2 is the variance of the
    # proportional residual; SD = sqrt(variance) per the footnote-c %CV definition.
    propSd <- 0.379474   ; label("Proportional residual error on dilmapimod Cc (fraction)")                     # Yang 2016 Table 2: Sigma variance = 0.144 -> SD = sqrt(0.144) = 0.3795 (CV 37.9%)

    # ===== Base CRP indirect-response model (Yang 2016 Table 3; no drug effect retained) =====
    # dA/dt = K_in0 * exp(-K_decline * t) - K_out * A, with A(0) = 1.35 mg/L (CRP at injury),
    # where t is the time since injury (Yang 2016 Equations 3 and 4). The final drug-effect
    # search did not retain a statistically significant dilmapimod term (Yang 2016 Results
    # section 3.3.1), so the model below is identical for active-drug and placebo subjects.
    lkdecline <- log(0.008)   ; label("CRP production-rate decline constant K_decline (1/h)")                  # Yang 2016 Table 3: exp(theta_pop,1) = 0.008,  theta_pop,1 = -4.77
    lkin0     <- log(7.171)   ; label("CRP production rate at the moment of injury K_in0 (mg/L per h)")        # Yang 2016 Table 3: exp(theta_pop,2) = 7.171,  theta_pop,2 = 1.97. Units are stated as 1/h in the paper's Table 3 header, but the indirect-response mass balance dA/dt = K_in - K_out * A with A in mg/L requires K_in to have units of mg/L per h.
    lkout     <- log(0.026)   ; label("CRP first-order elimination rate K_out (1/h)")                          # Yang 2016 Table 3: exp(theta_pop,3) = 0.026,  theta_pop,3 = -3.64

    # IIV on CRP parameters (log-normal; reported variances in Table 3)
    etalkdecline ~ 0.377   # Yang 2016 Table 3: IIV K_decline variance 0.377, CV 61.4%
    etalkin0     ~ 0.230   # Yang 2016 Table 3: IIV K_in0     variance 0.230, CV 48.0%
    etalkout     ~ 0.373   # Yang 2016 Table 3: IIV K_out     variance 0.373, CV 61.1%

    # Residual error on CRP. The paper reports a single Sigma in Table 3 with %CV computed
    # as 100 * sqrt(variance) (same footnote-c convention as the PK model); proportional
    # error is the only interpretation consistent with that %CV definition.
    propSd_crpobs <- 0.464758  ; label("Proportional residual error on CRP observation (fraction)")            # Yang 2016 Table 3: Sigma variance = 0.216 -> SD = sqrt(0.216) = 0.4648 (CV 46.5%)
  })

  model({
    # Reference BMI for the centred-power covariate effect (population mean per Yang 2016
    # Results section 3.2; 27.4 kg/m^2 was used as the centring constant, not the rounded
    # 27 kg/m^2 quoted in section 3.1).
    bmi_ref <- 27.4

    # Individual PK parameters; BMI as a power covariate on CL and Q2 only (Yang 2016
    # Table 2: covariate effects on CL and Q2 only retained after backward elimination).
    cl   <- exp(lcl + etalcl) * (BMI / bmi_ref)^e_bmi_cl
    vc   <- exp(lvc)
    q    <- exp(lq)
    vp   <- exp(lvp)
    q2   <- exp(lq2 + etalq2) * (BMI / bmi_ref)^e_bmi_q2
    vp2  <- exp(lvp2)

    # Individual CRP parameters. No drug-effect link (Yang 2016 Results section 3.3.1: none
    # of the inhibitory / stimulatory drug-effect models on K_in0, K_decline, or K_out
    # produced an acceptable fit).
    kdecline <- exp(lkdecline + etalkdecline)
    kin0     <- exp(lkin0     + etalkin0)
    kout     <- exp(lkout     + etalkout)

    # Three-compartment IV PK with first-order elimination from the central compartment.
    # Dose lands in `central` via the user data set's cmt column; infusion duration and
    # rate are encoded on the dose record (4-h infusion for cohorts 1 and 3; 24-h infusion
    # for cohorts 2 and 4).
    d/dt(central)     <-  q  / vp  * peripheral1 + q2 / vp2 * peripheral2 -
                          (cl + q + q2) / vc * central
    d/dt(peripheral1) <-  q  / vc  * central     - q  / vp  * peripheral1
    d/dt(peripheral2) <-  q2 / vc  * central     - q2 / vp2 * peripheral2

    # CRP indirect response; simulation time `t` is interpreted as time since injury so
    # the production-rate decay K_in0 * exp(-K_decline * t) matches Yang 2016 Equation 4
    # without any offset. The baseline A(0) = 1.35 mg/L is the assumed healthy-subject CRP
    # level at the moment of injury (Yang 2016 section 2.5.2.1; reference [11] in the paper).
    crp(0)    <- 1.35
    d/dt(crp) <- kin0 * exp(-kdecline * t) - kout * crp

    # Observations.
    # Dilmapimod plasma concentration in ng/mL: dose in mg with vc in L gives mg/L; the
    # x1000 factor converts to ng/mL to match the units used in the paper's Methods
    # section 2.2 (LLOQ 0.1 ng/mL, ULOQ 100 ng/mL on the calibration curve).
    Cc      <- central / vc * 1000
    # CRP serum level (mg/L) as measured by the immunoturbidimetric assay; Yang 2016
    # Methods section 2.3 reports an assay linear range of 1.0-480 mg/L.
    crpobs  <- crp

    Cc     ~ prop(propSd)
    crpobs ~ prop(propSd_crpobs)
  })
}
