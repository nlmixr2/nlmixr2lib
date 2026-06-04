Edstein_2001_tafenoquine <- function() {
  description <- "One-compartment population PK model for oral tafenoquine in 135 male Thai soldiers receiving 400 mg base for malaria prophylaxis (monthly n=104 or weekly n=31). The final model carries correlated IIV on apparent clearance and apparent volume of distribution (rho ~ 0.71) plus separate IIV on the first-order absorption rate constant; no covariates retained (centred age and weight on V/F and a prior-malaria indicator on CL/F were screened but not deemed to have sufficient clinical impact to alter the base model)."
  reference <- paste(
    "Edstein MD, Kocisko DA, Brewer TG, Walsh DS, Eamsila C, Charles BG (2001).",
    "Population pharmacokinetics of the new antimalarial agent tafenoquine in",
    "Thai soldiers.",
    "Br J Clin Pharmacol 52(6):663-670.",
    "doi:10.1046/j.0306-5251.2001.01482.x.",
    sep = " "
  )
  vignette <- "Edstein_2001_tafenoquine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  # Screened-but-not-retained covariates documented in covariatesDataExcluded
  # per references/parameter-names.md. WT and AGE (centred at the population
  # means 60.3 kg and 28.9 y) reduced the OFV when added to V/F (Table 2 model
  # numbers 4 and 5) and the malaria-status indicator MAL reduced the OFV when
  # added to CL/F (model number 6); the authors concluded "in view of the
  # relatively small changes in the pharmacokinetic parameter values, the base
  # model (model no. 1) was deemed to be adequate" (Results) and adopted the
  # no-covariate base model as final.
  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at the population mean 60.3 kg (Table 1 mean). Tested on V/F (Table 2 model number 4): reduction in OFV of 12 (P < 0.01) over base model. Not retained in the final model -- typical V/F at the weight extremes 45 kg and 90 kg differed by ~9 % (1633 L vs 1947 L) and was judged insufficient to alter dosing.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at the population mean 28.9 y (Table 1 mean). Tested on V/F (Table 2 model number 5): reduction in OFV of 8 (P < 0.01) over base model. Removed in backwards stepwise elimination (Table 2 model number 11) with OFV increase of only 2 vs the combined WT+AGE model.",
      source_name        = "AGE"
    ),
    MAL = list(
      description        = "Indicator of contracting malaria during the placebo phase of the trial (1 = subject contracted malaria and was re-treated; 0 = malaria-free subject on monthly prophylaxis throughout)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (malaria-free monthly cohort)",
      notes              = "Tested on CL/F as a piecewise effect 'CL/F = theta_1 * MAL + theta_2 * (1 - MAL)' (Table 2 models 6 and 8). Reduction in OFV of 58 when added (model 6), but the typical CL/F values for the two groups were similar (4.0 L/h malaria-affected vs 3.1 L/h not affected; ~28 % increase). Not retained because the authors adopted 'the parsimonious approach' (Discussion).",
      source_name        = "MAL"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 135L,
    n_studies       = 1L,
    age_range       = "21-46 years",
    age_median      = "28.9 years (mean; median not tabulated)",
    weight_range    = "45-90 kg",
    weight_median   = "60.3 kg (mean; median not tabulated)",
    sex_female_pct  = 0,
    race_ethnicity  = c(Asian = 100),
    disease_state   = "healthy adult male Thai soldiers, G6PD-normal, on tafenoquine malaria prophylaxis; 104 received monthly prophylaxis and 31 who contracted malaria during the placebo phase were re-treated and then maintained on weekly prophylaxis",
    dose_range      = "400 mg tafenoquine base orally: loading 400 mg daily for 3 days, then 400 mg monthly for 5 consecutive months (n = 104, monthly cohort); for the 31 subjects who contracted malaria during the placebo period, re-treatment with the same artesunate+doxycycline regime followed by 400 mg loading daily for 3 days and then 400 mg weekly (weekly cohort)",
    regions         = "Thailand (Thai-Cambodian border, Ubol Ratchatani province), Phase II clinical trial Apr-Oct 1998",
    notes           = "Demographics from Edstein 2001 Table 1 (mean and range). All subjects were deployed on security operations. Pre-treatment for parasite clearance: artesunate (300 mg day 1, 120 mg days 2-3) plus 200 mg daily doxycycline for 7 days before tafenoquine. Doses were taken with food (cake and biscuits) and 80-100 mL water. Plasma tafenoquine assayed by reversed-phase HPLC with fluorescence detection (LLOQ 10 ng/mL); interday and intraday CVs <= 8.4 %, mean recovery 81 %. Blood sampling was sparse and randomised in the field (mean 12.6 +/- 7.1 volunteers per collection)."
  )

  ini({
    # Structural parameters (Edstein 2001 Table 3, 'Parameter values for
    # population model number 1' -- the final model adopted from the forward /
    # backward covariate development summarised in Table 2). The reported CV %
    # next to each value is the RSE of the point estimate, not the IIV.
    lka  <- log(0.694) ; label("First-order absorption rate constant (Ka, 1/h)")     # Table 3 theta_3: Ka = 0.694 /h (RSE 12.8 %)
    lcl  <- log(3.20)  ; label("Apparent clearance CL/F (L/h)")                       # Table 3 theta_1: CL/F = 3.20 L/h (RSE 2.7 %)
    lvc  <- log(1820)  ; label("Apparent central volume of distribution V/F (L)")     # Table 3 theta_2: V/F = 1820 L (RSE 1.7 %)

    # Inter-individual variability. The paper writes the intersubject variance
    # model as exp(eta) on CL/F, V/F, and Ka and reports the random-effect
    # magnitudes as CV % (Table 3, 'Variance parameters' rows). Translate to
    # the internal log-scale variance via omega^2 = log(CV^2 + 1):
    #   CL/F : log(0.253^2 + 1) = 0.06204
    #   V/F  : log(0.148^2 + 1) = 0.02167
    #   Ka   : log(0.612^2 + 1) = 0.31791
    # A further run estimated the off-diagonal element via $OMEGA BLOCK(2)
    # (paper Results: OFV -4556 -> -4601); cov(CL/F, V/F) = 0.0265 reported
    # directly on the omega scale. Back-computed correlation from the values
    # above is 0.0265 / sqrt(0.06204 * 0.02167) = 0.72 (paper reports 0.71;
    # modest rounding from the CV % values quoted to three significant figures
    # in Table 3).
    etalcl + etalvc ~ c(0.06204,
                        0.0265, 0.02167)                                              # Table 3: 25.3 % CV CL/F (RSE 14.0); 14.8 % CV V/F (RSE 21.1); cov(CL/F,V/F) = 0.0265 (RSE 21.5) via $OMEGA BLOCK(2)
    etalka          ~ 0.31791                                                         # Table 3: 61.2 % CV Ka (RSE 167)

    # Residual variability. The paper Methods 'Population modelling' specifies
    # an exponential error model C_obs_ij = C_pred_ij * exp(eps_ij) with
    # variance sigma^2; Table 3 reports sigma = 17.9 % (CV) (RSE 8.60). In
    # nlmixr2, 'Cc ~ prop(propSd)' encodes additive-on-log-scale (== exponential
    # / log-normal proportional on linear scale), per the convention used by
    # Birgersson 2016 artemisinin and Birgersson 2019 artesunate. propSd
    # carries the proportional SD coefficient directly.
    propSd <- 0.179 ; label("Proportional residual error (fraction)")                 # Table 3 sigma: 17.9 % CV (RSE 8.60)
  })

  model({
    # Individual PK parameters. No covariates retained in the final model;
    # see Edstein 2001 Results: 'in view of the relatively small changes in
    # the pharmacokinetic parameter values, the base model (model no. 1) was
    # deemed to be adequate.'
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    # Elimination micro-constant from the single (central) disposition
    # compartment. Derived t_1/2 = log(2) * V/F / CL/F = 0.693 * 1820 / 3.20
    # = 394 h = 16.4 days (matches Table 3).
    kel <- cl / vc

    # ODE system: depot -> central -> elimination (one-compartment with
    # first-order absorption and elimination per Methods 'Population
    # modelling' and Results).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Plasma concentration. Dose is in mg, V/F is in L, so central / vc has
    # units of mg/L = 1000 ng/mL; multiply by 1000 to match the paper's
    # reported ng/mL (LLOQ = 10 ng/mL per Methods 'Drug analysis'; peak
    # ~260 ng/mL after the first 400 mg loading dose per Discussion).
    Cc <- central / vc * 1000

    # Residual error: NONMEM exponential / log-additive error maps to
    # proportional in nlmixr2's linear space.
    Cc ~ prop(propSd)
  })
}
