Zheng_2018_azithromycin <- function() {
  description <- "Pediatric population PK model for intravenous azithromycin in children with community-acquired pneumonia (Zheng 2018). Two-compartment model with linear elimination, allometric scaling on clearance and intercompartmental clearance (exponent 0.75 fixed) and on central and peripheral volumes (exponent 1.0 fixed) with reference body weight 21.5 kg, and a binary alanine aminotransferase covariate that reduces CL by 24 percent when ALT > 40 IU/L."
  reference <- "Zheng Y, Liu S-P, Xu B-P, Shi Z-R, Wang K, Yang J-B, Huang X, Tang B-H, Chen X-K, Shi H-Y, Zhou Y, Wu Y-E, Qi H, Jacqz-Aigrain E, Shen A-D, Zhao W. Population pharmacokinetics and dosing optimization of azithromycin in children with community-acquired pneumonia. Antimicrob Agents Chemother. 2018 Aug 27;62(9):e00686-18. doi:10.1128/AAC.00686-18"
  vignette <- "Zheng_2018_azithromycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at sampling time.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric exponents fixed a priori at 0.75 on CL and Q and at 1 on V1 and V2 (Zheng 2018 Methods 'Covariate analysis'). Reference weight 21.5 kg is the population median (Zheng 2018 Table 1 and Table 3 footnote).",
      source_name        = "WT"
    ),
    ALT = list(
      description        = "Serum alanine aminotransferase activity.",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Binarized inline as fliver <- (ALT > 40) per Zheng 2018 Table 3: F_liver = 0 if ALT <= 40, F_liver = 1 if ALT > 40. Multiplicative power-form effect on CL: 0.761^F_liver (a 24 percent reduction when ALT > 40 IU/L).",
      source_name        = "ALT"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 95,
    n_studies       = 1,
    age_range       = "2.1-11.7 years",
    age_median      = "5.9 years",
    weight_range    = "11.0-51.0 kg",
    weight_median   = "21.5 kg",
    sex_female_pct  = 44.2,
    disease_state   = "Hospitalized pediatric patients with suspected or confirmed community-acquired pneumonia (CAP).",
    dose_range      = "10 mg/kg intravenous infusion over 60 min once daily.",
    regions         = "China (multicenter: Beijing Children's Hospital, Children's Hospital of Hebei Province, Shandong Provincial Qianfoshan Hospital, Xintai People's Hospital).",
    n_observations  = "140 plasma azithromycin concentrations from 95 patients.",
    alt_range       = "3.3-374.4 IU/L (median 20.3; mean 26.8 +/- 40.2).",
    notes           = "Demographic counts and ranges reproduced from Zheng 2018 Table 1. Sex counts: 53 males / 42 females. An external validation cohort of 28 additional children (age 2.0-11.0 years, weight 10.0-50.0 kg) is reported separately in the source paper but was not used in model building."
  )

  ini({
    # Structural parameters at reference weight 21.5 kg and reference liver status ALT <= 40 IU/L.
    lcl <- log(27.8); label("Clearance, CL, at 21.5 kg with ALT <= 40 IU/L (L/h)")    # Zheng 2018 Table 3: theta1 = 27.8 L/h (RSE 7.1 percent)
    lvc <- log(39.5); label("Central volume of distribution, V1, at 21.5 kg (L)")      # Zheng 2018 Table 3: theta2 = 39.5 L (RSE 43.8 percent)
    lvp <- log(377);  label("Peripheral volume of distribution, V2, at 21.5 kg (L)")   # Zheng 2018 Table 3: theta3 = 377 L (RSE 13.4 percent)
    lq  <- log(55.7); label("Inter-compartmental clearance, Q, at 21.5 kg (L/h)")      # Zheng 2018 Table 3: theta4 = 55.7 L/h (RSE 22.4 percent)

    # Allometric exponents held fixed at the paper's a priori values (Zheng 2018 Methods 'Covariate analysis'); reported in Table 3 as the constants 0.75 and 1.
    allo_cl <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)")  # Zheng 2018 Methods 'Covariate analysis' and Table 3 equations
    allo_v  <- fixed(1.0);  label("Allometric exponent on V1 and V2 (unitless)") # Zheng 2018 Methods 'Covariate analysis' and Table 3 equations

    # Binarized ALT covariate effect on CL, expressed on the log scale: exp(e_alt_cl * fliver) reproduces the paper's power form theta5^F_liver.
    e_alt_cl <- log(0.761); label("Log-effect of ALT > 40 IU/L on CL (unitless)")  # Zheng 2018 Table 3: theta5 = 0.761 (RSE 14.7 percent; -24 percent on CL when ALT > 40)

    # IIV - Zheng 2018 Table 3 reports percent CV for an exponential IIV
    # model (Methods 'Model development'). Variances stored here use the
    # log-normal conversion omega^2 = log(CV^2 + 1).
    etalcl ~ 0.09805  # Zheng 2018 Table 3: 32.1 percent CV -> log(0.321^2 + 1) = 0.09805
    etalvc ~ 0.54299  # Zheng 2018 Table 3: 84.9 percent CV -> log(0.849^2 + 1) = 0.54299
    etalq  ~ 0.23624  # Zheng 2018 Table 3: 51.6 percent CV -> log(0.516^2 + 1) = 0.23624

    # Residual error - Zheng 2018 Methods reports an exponential residual
    # model fitted to the data (Y = IPRED * exp(EPS)); Table 3 reports
    # 5.7 percent CV with RSE 161.6 percent. Encoded as a log-normal residual
    # error with SD on the log scale.
    expSd <- 0.057; label("Log-normal residual error SD (unitless)")  # Zheng 2018 Table 3: 5.7 percent CV (exponential model; RSE 161.6 percent)
  })

  model({
    # Reference covariate values (Zheng 2018 Table 1 median weight, Table 3 ALT threshold).
    ref_wt  <- 21.5
    alt_thr <- 40

    # Binary liver-function indicator derived from continuous ALT - matches the paper's NONMEM coding F_liver = 0 if ALT <= 40, F_liver = 1 if ALT > 40.
    fliver <- (ALT > alt_thr)

    # Individual structural parameters with allometric scaling on body
    # weight and the multiplicative ALT effect on CL.
    cl <- exp(lcl + etalcl) * (WT / ref_wt)^allo_cl * exp(e_alt_cl * fliver)
    vc <- exp(lvc + etalvc) * (WT / ref_wt)^allo_v
    vp <- exp(lvp)          * (WT / ref_wt)^allo_v
    q  <- exp(lq  + etalq)  * (WT / ref_wt)^allo_cl

    # Micro-constants for the explicit two-compartment ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system - intravenous dosing routes directly into the central
    # compartment (no depot); the 60-min infusion duration is supplied by
    # the user via the rate column of the event table.
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
