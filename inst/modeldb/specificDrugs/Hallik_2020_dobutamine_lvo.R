Hallik_2020_dobutamine_lvo <- function() {
  description <- "Simultaneous PKPD model for intravenous dobutamine and left ventricular cardiac output (LVO, echocardiography) in critically ill preterm and term neonates in the first 3 days of life (Hallik 2020, Table 4). PK is the paper's one-compartment linear model re-estimated jointly with LVO: clearance scales with birth weight^0.75 and a sigmoidal postmenstrual-age maturation function (PMA50 and Hill fixed from the PK-only fit), volume scales linearly with birth weight, both referenced to 1618 g, and one shared random effect enters volume multiplied by an estimated scale factor. LVO follows a sigmoidal Emax model in the plasma concentration, rising from a baseline of 131 mL/kg/min toward a maximum LEVEL of 157 mL/kg/min (the paper's Emax is the plateau value, not the increment) with EC50 117 ug/L and Hill 2.82, and no effect-compartment delay. Residual errors are proportional for both concentration and LVO. The PK-only model is Hallik_2020_dobutamine."
  reference <- paste(
    "Hallik M, Ilmoja M-L, Standing JF, Soeorg H, Jalas T, Raidmae M, Uibo K,",
    "Kobas K, Sonajalg M, Takkis K, Veigure R, Kipper K, Starkopf J, Metsvaht T.",
    "Population pharmacokinetics and pharmacodynamics of dobutamine in neonates",
    "on the first days of life.",
    "Br J Clin Pharmacol. 2020;86(2):318-328.",
    "doi:10.1111/bcp.14146.",
    sep = " "
  )
  vignette <- "Hallik_2020_dobutamine"

  # Dose in ug and time in hours, so an infusion of R ug/kg/min enters as a rate
  # of R * WT * 60 ug/h. Concentrations are ug/L; LVO is mL/kg/min.
  units <- list(time = "h", dosing = "ug", concentration = "ug/L", lvo = "mL/kg/min")

  compartmentData <- list(
    central = list(analyte = "dobutamine", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT_BIRTH = list(
      description = "Birth weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric size descriptor for CL (exponent 0.75, fixed) and V (exponent 1), normalised to the population median birth weight of 1618 g (Equations 1-2; Table 1). The paper enters grams; this file takes kilograms with a 1.618 kg reference, the same ratio. All infants were studied in the first 72 h of life, so birth weight stands in for current weight.",
      source_name = "BW"
    ),
    PAGE = list(
      description = "Postmenstrual age (gestational age at birth plus postnatal age)",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the sigmoidal maturation of CL (Equation 1). PMA50 = 37.4 weeks and Hill = 2.67 are fixed at the PK-only estimates of Table 2 (Methods 2.3). Declared in weeks, as the paper states it. Table 1 gestational age at birth median 30.4 weeks (22.7-41.0).",
      source_name = "PMA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 28,
    n_studies = 1,
    age_range = "Postnatal age at recruitment 2-28 h (median 6 h); gestational age at birth 22.7-41.0 weeks (median 30.4)",
    weight_range = "Birth weight 465-4380 g (median 1618 g)",
    sex_female_pct = 36,
    disease_state = "Critically ill preterm and term neonates in the first 72 h of life needing inotropic support on clinical grounds. Baseline LVO median 128 mL/kg/min (71-338), Table 1.",
    dose_range = "Continuous IV infusion started at 5 ug/kg/min and raised by 5 ug/kg/min roughly every 30 min to a maximum of 20 ug/kg/min",
    regions = "Estonia (Tallinn Children's Hospital and Tartu University Hospital)",
    notes = "Prospective 2-centre study (EU CTR 2015-004836-36); 28 of 31 recruited neonates analysed. LVO was measured by echocardiography before the infusion and about 20-30 min after each dose escalation (Methods 2.2). See Hallik_2020_dobutamine for the full demographic summary."
  )

  ini({
    # PK parameters re-estimated in the joint fit (Table 4, 'PKPD model for LVO
    # effect'); CL is the fully mature value for a 1618 g neonate.
    lcl <- log(40.7); label("Clearance for a 1618 g neonate at full maturation (L/h)")  # Table 4 LVO row 'CL (L h-1 1618-g-1)' = 40.7 (SE 3.03)
    lvc <- log(5.14); label("Volume of distribution for a 1618 g neonate (L)")  # Table 4 LVO row 'V (L 1618-g-1)' = 5.14 (SE 0.726)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of birth weight on CL (unitless)")  # Equation 1 literal exponent 0.75
    # Methods 2.3: PKPD models used 'the final linear PK structural model with
    # Hill coefficient and PMA50 fixed to values estimated from PK data (Table 2)'.
    ltmat50 <- fixed(log(37.4)); label("Postmenstrual age at 50% of mature CL (weeks)")  # Table 2 'PMA 50 (weeks)' = 37.4, fixed per Methods 2.3
    lhill_mat <- fixed(log(2.67)); label("Hill coefficient of the CL maturation function (unitless)")  # Table 2 'Hill' = 2.67, fixed per Methods 2.3
    vc_eta_scale <- 1.33; label("Scale factor applied to the shared CL random effect on V (unitless)")  # Table 4 LVO row 'Shared BSV scale factor' = 1.33 (SE 0.490)

    # PD parameters: sigmoidal Emax E = E0 + (Emax - E0) * C^g / (EC50^g + C^g)
    # (Equation 6). Emax is 'the estimated maximum HD parameter value', i.e.
    # the plateau LEVEL, hence lrmax_ rather than lemax.
    lrbase_lvo <- log(131); label("Baseline left ventricular output E0 (mL/kg/min)")  # Table 4 LVO row 'E0 (mL kg-1 min-1)' = 131 (SE 9.64)
    lrmax_lvo <- log(157); label("Maximum attainable left ventricular output Emax (mL/kg/min)")  # Table 4 LVO row 'Emax (mL kg-1 min-1)' = 157 (SE 21.8)
    lec50_lvo <- log(117); label("Plasma concentration at half-maximal LVO change (ug/L)")  # Table 4 LVO row 'EC50 (ug L-1)' = 117 (SE 37.0), no BSV
    lhill_lvo <- log(2.82); label("Hill coefficient of the LVO concentration-effect curve (unitless)")  # Table 4 LVO row 'gamma' = 2.82 (SE 1.25), no BSV

    # Table 4 footnote a: BSV CV = sqrt(omega2) x 100%, so omega2 = CV^2.
    etalcl ~ 0.0625  # Table 4 LVO CL BSV 25% (SE 17.5%), shrinkage 21%; omega2 = 0.25^2
    etalrbase_lvo ~ 0.1296  # Table 4 LVO E0 BSV 36% (SE 19.8%), shrinkage 3%; omega2 = 0.36^2
    etalrmax_lvo ~ 0.1936  # Table 4 LVO Emax BSV 44% (SE 39.7%), shrinkage 29%; omega2 = 0.44^2

    propSd <- 0.589; label("Proportional residual error, dobutamine concentration (fraction)")  # Table 4 LVO row 'Pharmacokinetic residual error (proportional)' = 0.589 (SE 0.053)
    propSd_lvo <- 0.167; label("Proportional residual error, LVO (fraction)")  # Table 4 LVO row 'Pharmacodynamic residual error (proportional)' = 0.167 (SE 0.015)
  })

  model({
    # Size and maturation (Equations 1-2), reference birth weight 1.618 kg.
    tmat50 <- exp(ltmat50)
    hill_mat <- exp(lhill_mat)
    fmat <- PAGE^hill_mat / (tmat50^hill_mat + PAGE^hill_mat)

    cl <- exp(lcl + etalcl) * (WT_BIRTH / 1.618)^e_wt_cl * fmat
    vc <- exp(lvc + vc_eta_scale * etalcl) * (WT_BIRTH / 1.618)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    Cc <- central / vc

    # Sigmoidal Emax driven by the plasma concentration (Equation 6); Table 3
    # selects this form for LVO, without an effect compartment.
    rbase_lvo <- exp(lrbase_lvo + etalrbase_lvo)
    rmax_lvo <- exp(lrmax_lvo + etalrmax_lvo)
    ec50_lvo <- exp(lec50_lvo)
    hill_lvo <- exp(lhill_lvo)
    lvo <- rbase_lvo + (rmax_lvo - rbase_lvo) * Cc^hill_lvo / (ec50_lvo^hill_lvo + Cc^hill_lvo)

    Cc ~ prop(propSd)
    lvo ~ prop(propSd_lvo)
  })
}
