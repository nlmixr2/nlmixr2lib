Hallik_2020_dobutamine_map <- function() {
  description <- "Simultaneous PKPD model for intravenous dobutamine and mean arterial blood pressure (MAP, continuously monitored and averaged over 1 min) in critically ill preterm and term neonates in the first 3 days of life (Hallik 2020, Table 5). PK is the paper's one-compartment linear model re-estimated jointly with MAP: clearance scales with birth weight^0.75 and a sigmoidal postmenstrual-age maturation function (PMA50 and Hill fixed from the PK-only fit), volume scales linearly with birth weight, both referenced to 1618 g, and one shared random effect enters volume multiplied by an estimated scale factor. MAP follows a steep sigmoidal Emax model in the plasma concentration (no effect-compartment delay), rising from a baseline of 39.7 mmHg toward a maximum LEVEL of 41.9 mmHg (the paper's Emax is the plateau value, not the increment) with EC50 25.4 ug/L and Hill 13.5. Residual errors are proportional for both concentration and MAP. The PK-only model is Hallik_2020_dobutamine."
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
  # of R * WT * 60 ug/h. Concentrations are ug/L; MAP is mmHg.
  units <- list(time = "h", dosing = "ug", concentration = "ug/L", map = "mmHg")

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
    disease_state = "Critically ill preterm and term neonates in the first 72 h of life needing inotropic support on clinical grounds.",
    dose_range = "Continuous IV infusion started at 5 ug/kg/min and raised by 5 ug/kg/min roughly every 30 min to a maximum of 20 ug/kg/min",
    regions = "Estonia (Tallinn Children's Hospital and Tartu University Hospital)",
    notes = "Prospective 2-centre study (EU CTR 2015-004836-36); 28 of 31 recruited neonates analysed. MAP was recorded every 2.5 s; the 15 min around each dose change (5 min before, 10 min after) was averaged over 1-min bins for the analysis (Methods 2.3). See Hallik_2020_dobutamine for the full demographic summary."
  )

  ini({
    # PK parameters re-estimated in the joint fit (Table 5, 'PKPD model for MAP
    # effect'); CL is the fully mature value for a 1618 g neonate.
    lcl <- log(37.2); label("Clearance for a 1618 g neonate at full maturation (L/h)")  # Table 5 MAP row 'CL (L h-1 1618-g-1)' = 37.2 (SE 2.88)
    lvc <- log(4.88); label("Volume of distribution for a 1618 g neonate (L)")  # Table 5 MAP row 'V (L 1618-g-1)' = 4.88 (SE 0.958)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of birth weight on CL (unitless)")  # Equation 1 literal exponent 0.75
    # Methods 2.3: PKPD models used 'the final linear PK structural model with
    # Hill coefficient and PMA50 fixed to values estimated from PK data (Table 2)'.
    ltmat50 <- fixed(log(37.4)); label("Postmenstrual age at 50% of mature CL (weeks)")  # Table 2 'PMA 50 (weeks)' = 37.4, fixed per Methods 2.3
    lhill_mat <- fixed(log(2.67)); label("Hill coefficient of the CL maturation function (unitless)")  # Table 2 'Hill' = 2.67, fixed per Methods 2.3
    vc_eta_scale <- 3.62; label("Scale factor applied to the shared CL random effect on V (unitless)")  # Table 5 MAP row 'Shared BSV scale factor' = 3.62 (SE 0.760)

    # PD parameters: sigmoidal Emax E = E0 + (Emax - E0) * C^g / (EC50^g + C^g)
    # (Equation 6). Emax is 'the estimated maximum HD parameter value', i.e.
    # the plateau LEVEL, hence lrmax_ rather than lemax.
    lrbase_map <- log(39.7); label("Baseline mean arterial pressure E0 (mmHg)")  # Table 5 MAP row 'E0 (mmHg)' = 39.7 (SE 1.75)
    lrmax_map <- log(41.9); label("Maximum attainable mean arterial pressure Emax (mmHg)")  # Table 5 MAP row 'Emax (mmHg)' = 41.9 (SE 2.19)
    lec50_map <- log(25.4); label("Plasma concentration at half-maximal MAP change (ug/L)")  # Table 5 MAP row 'EC50 (ug L-1)' = 25.4 (SE 2.00), no BSV
    lhill_map <- log(13.5); label("Hill coefficient of the MAP concentration-effect curve (unitless)")  # Table 5 MAP row 'gamma' = 13.5 (SE 2.94), no BSV

    # Table 5 footnote a: BSV CV = sqrt(omega2) x 100%, so omega2 = CV^2.
    etalcl ~ 0.0576  # Table 5 MAP CL BSV 24% (SE 16.6%), shrinkage 4%; omega2 = 0.24^2
    etalrbase_map ~ 0.0484  # Table 5 MAP E0 BSV 22% (SE 11.8%), shrinkage 2%; omega2 = 0.22^2
    etalrmax_map ~ 0.0676  # Table 5 MAP Emax BSV 26% (SE 13.9%), shrinkage 4%; omega2 = 0.26^2

    propSd <- 0.675; label("Proportional residual error, dobutamine concentration (fraction)")  # Table 5 MAP row 'Pharmacokinetic residual error (proportional)' = 0.675 (SE 0.062)
    propSd_map <- 0.065; label("Proportional residual error, MAP (fraction)")  # Table 5 MAP row 'Pharmacodynamic residual error (proportional)' = 0.065 (SE 0.001)
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
    # selects 'Sigmoidal Emax' without KEO for MAP.
    rbase_map <- exp(lrbase_map + etalrbase_map)
    rmax_map <- exp(lrmax_map + etalrmax_map)
    ec50_map <- exp(lec50_map)
    hill_map <- exp(lhill_map)
    map <- rbase_map + (rmax_map - rbase_map) * Cc^hill_map / (ec50_map^hill_map + Cc^hill_map)

    Cc ~ prop(propSd)
    map ~ prop(propSd_map)
  })
}
