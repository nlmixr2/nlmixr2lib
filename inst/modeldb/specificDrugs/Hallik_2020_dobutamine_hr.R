Hallik_2020_dobutamine_hr <- function() {
  description <- "Simultaneous PKPD model for intravenous dobutamine and heart rate (HR, continuously monitored and averaged over 1 min) in critically ill preterm and term neonates in the first 3 days of life (Hallik 2020, Table 5). PK is the paper's one-compartment linear model re-estimated jointly with HR: clearance scales with birth weight^0.75 and a sigmoidal postmenstrual-age maturation function (PMA50 and Hill fixed from the PK-only fit), volume scales linearly with birth weight, both referenced to 1618 g, and one shared random effect enters volume multiplied by an estimated scale factor. HR follows a sigmoidal Emax model in the EFFECT-COMPARTMENT concentration (keo 6.59 1/h, mean equilibration time about 9 min), rising from a baseline of 138 beats/min toward a maximum LEVEL of 172 beats/min (the paper's Emax is the plateau value, not the increment) with EC50 39.2 ug/L and Hill 3.36. Residual errors are proportional for both concentration and HR. The PK-only model is Hallik_2020_dobutamine."
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
  # of R * WT * 60 ug/h. Concentrations are ug/L; HR is beats/min.
  units <- list(time = "h", dosing = "ug", concentration = "ug/L", hr = "beats/min")

  # `effect` holds the effect-site CONCENTRATION (ug/L), per Equation 3, which
  # assumes 'the amount of drug distributing into the effect compartment does
  # not influence the overall PK'.
  compartmentData <- list(
    central = list(analyte = "dobutamine", units = "ug", specimen = "plasma", verified = TRUE),
    effect = list(analyte = "dobutamine", units = "ug/L", specimen = "not applicable", verified = TRUE)
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
    disease_state = "Critically ill preterm and term neonates in the first 72 h of life needing inotropic support on clinical grounds.",
    dose_range = "Continuous IV infusion started at 5 ug/kg/min and raised by 5 ug/kg/min roughly every 30 min to a maximum of 20 ug/kg/min",
    regions = "Estonia (Tallinn Children's Hospital and Tartu University Hospital)",
    notes = "Prospective 2-centre study (EU CTR 2015-004836-36); 28 of 31 recruited neonates analysed. HR was recorded every 2.5 s; the 15 min around each dose change (5 min before, 10 min after) was averaged over 1-min bins for the analysis (Methods 2.3). See Hallik_2020_dobutamine for the full demographic summary."
  )

  ini({
    # PK parameters re-estimated in the joint fit (Table 5, 'PKPD model for HR
    # effect'); CL is the fully mature value for a 1618 g neonate.
    lcl <- log(42.3); label("Clearance for a 1618 g neonate at full maturation (L/h)")  # Table 5 HR row 'CL (L h-1 1618-g-1)' = 42.3 (SE 3.22)
    lvc <- log(5.42); label("Volume of distribution for a 1618 g neonate (L)")  # Table 5 HR row 'V (L 1618-g-1)' = 5.42 (SE 0.766)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of birth weight on CL (unitless)")  # Equation 1 literal exponent 0.75
    # Methods 2.3: PKPD models used 'the final linear PK structural model with
    # Hill coefficient and PMA50 fixed to values estimated from PK data (Table 2)'.
    ltmat50 <- fixed(log(37.4)); label("Postmenstrual age at 50% of mature CL (weeks)")  # Table 2 'PMA 50 (weeks)' = 37.4, fixed per Methods 2.3
    lhill_mat <- fixed(log(2.67)); label("Hill coefficient of the CL maturation function (unitless)")  # Table 2 'Hill' = 2.67, fixed per Methods 2.3
    vc_eta_scale <- 1.72; label("Scale factor applied to the shared CL random effect on V (unitless)")  # Table 5 HR row 'Shared BSV scale factor' = 1.72 (SE 0.421)

    # Effect-compartment equilibration (Equation 3). Results 3.1: mean effect
    # time (keo^-1 x 60 min) of 9 min for HR; 60 / 6.59 = 9.1 min.
    lke0 <- log(6.59); label("Effect-compartment equilibration rate constant keo (1/h)")  # Table 5 HR row 'KE0 (h-1)' = 6.59 (SE 2.18), no BSV

    # PD parameters: sigmoidal Emax E = E0 + (Emax - E0) * Ce^g / (EC50^g + Ce^g)
    # (Equation 6). Emax is 'the estimated maximum HD parameter value', i.e.
    # the plateau LEVEL; Discussion 4.2: 'increase of mean HR from 138 to 172'.
    lrbase_hr <- log(138); label("Baseline heart rate E0 (beats/min)")  # Table 5 HR row 'E0 (min-1)' = 138 (SE 4.2)
    lrmax_hr <- log(172); label("Maximum attainable heart rate Emax (beats/min)")  # Table 5 HR row 'Emax (min-1)' = 172 (SE 2.5)
    lec50_hr <- log(39.2); label("Effect-site concentration at half-maximal HR change (ug/L)")  # Table 5 HR row 'EC50 (ug L-1)' = 39.2 (SE 5.56)
    lhill_hr <- log(3.36); label("Hill coefficient of the HR concentration-effect curve (unitless)")  # Table 5 HR row 'gamma' = 3.36 (SE 0.326), no BSV

    # Table 5 footnote a: BSV CV = sqrt(omega2) x 100%, so omega2 = CV^2.
    etalcl ~ 0.0729  # Table 5 HR CL BSV 27% (SE 18.0%), shrinkage 15%; omega2 = 0.27^2
    etalrbase_hr ~ 0.0225  # Table 5 HR E0 BSV 15% (SE 8.1%), shrinkage 4%; omega2 = 0.15^2
    etalec50_hr ~ 0.25  # Table 5 HR EC50 BSV 50% (SE 32.2%), shrinkage 21%; omega2 = 0.50^2
    etalrmax_hr ~ 0.0025  # Table 5 HR Emax BSV 5% (SE 3.1%), shrinkage 28%; omega2 = 0.05^2

    propSd <- 0.590; label("Proportional residual error, dobutamine concentration (fraction)")  # Table 5 HR row 'Pharmacokinetic residual error (proportional)' = 0.590 (SE 0.052)
    propSd_hr <- 0.051; label("Proportional residual error, HR (fraction)")  # Table 5 HR row 'Pharmacodynamic residual error (proportional)' = 0.051 (SE 0.001)
  })

  model({
    # Size and maturation (Equations 1-2), reference birth weight 1.618 kg.
    tmat50 <- exp(ltmat50)
    hill_mat <- exp(lhill_mat)
    fmat <- PAGE^hill_mat / (tmat50^hill_mat + PAGE^hill_mat)

    cl <- exp(lcl + etalcl) * (WT_BIRTH / 1.618)^e_wt_cl * fmat
    vc <- exp(lvc + vc_eta_scale * etalcl) * (WT_BIRTH / 1.618)
    ke0 <- exp(lke0)

    kel <- cl / vc

    d/dt(central) <- -kel * central
    Cc <- central / vc

    # Effect-site concentration (Equation 3): dCe/dt = keo * C - keo * Ce.
    d/dt(effect) <- ke0 * (Cc - effect)

    # Sigmoidal Emax driven by the effect-site concentration (Equation 6);
    # Table 3 selects 'Sigmoidal Emax' with KEO for HR.
    rbase_hr <- exp(lrbase_hr + etalrbase_hr)
    rmax_hr <- exp(lrmax_hr + etalrmax_hr)
    ec50_hr <- exp(lec50_hr + etalec50_hr)
    hill_hr <- exp(lhill_hr)
    hr <- rbase_hr + (rmax_hr - rbase_hr) * effect^hill_hr / (ec50_hr^hill_hr + effect^hill_hr)

    Cc ~ prop(propSd)
    hr ~ prop(propSd_hr)
  })
}
