Bijleveld_2017_gentamicin <- function() {
  description <- "Two-compartment population PK model of intravenous gentamicin in (pre)term neonates with suspected or proven Gram-negative sepsis (Bijleveld 2017), with fixed allometric body-weight scaling (exponents 0.75 on CL and Q, 1 on Vc and Vp) and an estimated postmenstrual-age power effect on CL."
  reference <- paste(
    "Bijleveld YA, van den Heuvel ME, Hodiamont CJ, Mathot RAA, de Haan TR.",
    "Population pharmacokinetics and dosing considerations for gentamicin in",
    "newborns with suspected or proven sepsis caused by gram-negative bacteria.",
    "Antimicrob Agents Chemother. 2017;61(1):e01304-16.",
    "doi:10.1128/AAC.01304-16.",
    sep = " "
  )
  vignette <- "Bijleveld_2017_gentamicin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at start of treatment",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Median 1.85 kg, range 0.76-4.67 kg in the source cohort (Table 1). Allometric scaling uses a 70 kg adult reference: fixed exponent 0.75 on CL and Q, fixed exponent 1 on Vc and Vp.",
      source_name        = "BW"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age plus postnatal age)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. The Bijleveld 2017 paper reports PMA in DAYS as the sum",
        "of gestational age (weeks * 7) and postnatal age (days); reference PMA",
        "in the paper's CL power equation is 225 days (median of the cohort).",
        "Canonical PAGE in nlmixr2lib carries months, so the model converts",
        "internally as PMA_days = PAGE * 30.4375 before evaluating",
        "(PMA_days / 225)^theta_clpma. Reference of 225 days corresponds to",
        "approximately 7.39 months."
      ),
      source_name        = "PMA"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 65L,
    n_studies       = 1L,
    n_observations  = 136L,
    ga_range        = "25-42 weeks (mean 32)",
    pna_range       = "0-31 days (median 1)",
    weight_range    = "0.76-4.67 kg at start of treatment (median 1.85)",
    birth_weight_range = "0.43-4.67 kg (mean 1.85)",
    sex_female_pct  = 43,
    disease_state   = paste(
      "(Pre)term neonates with suspected or proven Gram-negative sepsis",
      "treated in a neonatal intensive care unit. 75% born prematurely",
      "(GA < 37 weeks); 8% mortality during the observation window."
    ),
    dose_range      = paste(
      "4 mg/kg every 48 h (GA <= 32 wk, PNA < 4 wk), 4 mg/kg every 24 h",
      "(GA > 32 wk, PNA < 4 wk), or 5 mg/kg every 24 h (PNA > 4 wk),",
      "administered as a 30-min intravenous infusion."
    ),
    creatinine_summary = "Serum creatinine 26-112 umol/L (median 55.5).",
    regions         = "Single centre (Academic Medical Center NICU, Amsterdam, the Netherlands). Data collected October 2012 to January 2013.",
    notes           = paste(
      "Baseline demographics from Bijleveld 2017 Table 1. Three previously",
      "published gentamicin neonatal popPK models (Nielsen 2009, Fuchs 2014,",
      "De Cock 2014) were found to overpredict this cohort, motivating the new",
      "model. Final model retains only postmenstrual age as a covariate on CL;",
      "PNA, gender, and ibuprofen co-medication did not improve the fit. Fit in",
      "NONMEM 7.2 with FOCE-I."
    )
  )

  ini({
    # Structural population PK parameters (Bijleveld 2017 Table 3, Final
    # Model column). All four typical values are reported on the linear
    # scale per 70 kg reference; nlmixr2lib log-transforms positive-
    # constrained parameters.
    lcl <- log(1.00);  label("Clearance (CL, L/h per 70 kg)")                 # Bijleveld 2017 Table 3 Final model: 1.00 L/h/70 kg
    lvc <- log(37.4);  label("Central volume of distribution (Vc, L per 70 kg)") # Bijleveld 2017 Table 3 Final model: 37.4 L/70 kg
    lq  <- log(0.57);  label("Intercompartmental clearance (Q, L/h per 70 kg)")  # Bijleveld 2017 Table 3 Final model: 0.57 L/h/70 kg
    lvp <- log(17.1);  label("Peripheral volume of distribution (Vp, L per 70 kg)") # Bijleveld 2017 Table 3 Final model: 17.1 L/70 kg

    # Allometric body-weight exponents (Bijleveld 2017 Table 3 footnotes
    # b-e: TVCl uses (BW/70)^0.75, TVVc uses (BW/70)^1, TVQ uses (BW/70)^0.75,
    # TVVp uses (BW/70)^1). Reported without uncertainty; standard
    # West-style allometric exponents held constant.
    e_wt_cl_q   <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)") # Bijleveld 2017 Table 3 footnotes b, d
    e_wt_vc_vp  <- fixed(1);    label("Allometric exponent on Vc and Vp (unitless)") # Bijleveld 2017 Table 3 footnotes c, e

    # Postmenstrual-age power effect on CL (Bijleveld 2017 Table 3 Final
    # Model column: theta_ClPMA = 1.89, CV 26%). Enters as
    # (PMA_days / 225)^theta_ClPMA where 225 days is the cohort median.
    e_page_cl   <- 1.89;        label("Power exponent of (PMA/225 days) on CL (unitless)") # Bijleveld 2017 Table 3 Final model: theta_ClPMA = 1.89

    # Inter-individual variability (Bijleveld 2017 Table 3 Final Model
    # column). Reported as %CV; for log-normal eta the variance on the
    # log scale is omega^2 = log(1 + CV^2).
    #   IIV CL : 33.1% CV  -> log(1 + 0.331^2) = 0.10396
    #   IIV Vc : 24.4% CV  -> log(1 + 0.244^2) = 0.05784
    # CL and Vc are correlated (r = 0.08, Results "Pharmacokinetic model"
    # paragraph); the lower-triangular covariance is
    #   cov(CL, Vc) = r * omega_CL * omega_Vc
    #              = 0.08 * sqrt(0.10396) * sqrt(0.05784) = 0.006204.
    etalcl + etalvc ~ c(0.10396,
                        0.006204, 0.05784)  # Bijleveld 2017 Table 3: IIV CL 33.1%, IIV Vc 24.4%, corr 0.08
    # Q and Vp had no IIV in the final model (Bijleveld 2017 Results
    # "Pharmacokinetic model" paragraph: IIV could be estimated only for
    # CL and Vc).

    # Residual error. The paper describes residual variability on
    # log-transformed observations as an additive error (sigma = 0.43,
    # Bijleveld 2017 Table 3 Final Model column). The NONMEM LTBS pattern
    #   log(Cobs) = log(Cpred) + eps,  eps ~ N(0, expSd^2)
    # maps directly to nlmixr2's lnorm() residual with expSd = sigma.
    expSd <- 0.43;  label("Lognormal residual SD (additive on log-transformed concentration)") # Bijleveld 2017 Table 3 Final model: additive sigma on log-data = 0.43
  })

  model({
    # Convert canonical PAGE (months) to postmenstrual age in days for
    # the Bijleveld 2017 CL covariate equation, which uses a 225-day
    # cohort-median reference.
    pma_days <- PAGE * 30.4375

    # Individual PK parameters (Bijleveld 2017 Table 3 footnotes b-e).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q   * (pma_days / 225)^e_page_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp)          * (WT / 70)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV PK. Doses are administered as 30-min infusions
    # in the source paper; the library model does not hard-code the
    # infusion duration so users specify rate or dur per dose in their
    # event table.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Plasma concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
