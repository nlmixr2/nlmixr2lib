Hallik_2020_dobutamine_rvo <- function() {
  description <- "Simultaneous PKPD model for intravenous dobutamine and right ventricular cardiac output (RVO, echocardiography) in critically ill preterm and term neonates in the first 3 days of life (Hallik 2020, Table 4). PK is the paper's one-compartment linear model re-estimated jointly with RVO: clearance scales with birth weight^0.75 and a sigmoidal postmenstrual-age maturation function (PMA50 and Hill fixed from the PK-only fit), volume scales linearly with birth weight, both referenced to 1618 g, and one shared random effect enters volume multiplied by an estimated scale factor. RVO rises linearly with the plasma concentration from an estimated baseline of 151 mL/kg/min (slope 0.214 mL/kg/min per ug/L), with no effect-compartment delay. Residual errors are proportional for both concentration and RVO. The PK-only model is Hallik_2020_dobutamine."
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
  # of R * WT * 60 ug/h. Concentrations are ug/L; RVO is mL/kg/min.
  units <- list(time = "h", dosing = "ug", concentration = "ug/L", rvo = "mL/kg/min")

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
    disease_state = "Critically ill preterm and term neonates in the first 72 h of life needing inotropic support on clinical grounds. Baseline RVO median 136 mL/kg/min (75-306), Table 1.",
    dose_range = "Continuous IV infusion started at 5 ug/kg/min and raised by 5 ug/kg/min roughly every 30 min to a maximum of 20 ug/kg/min",
    regions = "Estonia (Tallinn Children's Hospital and Tartu University Hospital)",
    notes = "Prospective 2-centre study (EU CTR 2015-004836-36); 28 of 31 recruited neonates analysed. RVO was measured by echocardiography before the infusion and about 20-30 min after each dose escalation (Methods 2.2). See Hallik_2020_dobutamine for the full demographic summary."
  )

  ini({
    # PK parameters re-estimated in the joint fit (Table 4, 'PKPD model for RVO
    # effect'); CL is the fully mature value for a 1618 g neonate.
    lcl <- log(41.0); label("Clearance for a 1618 g neonate at full maturation (L/h)")  # Table 4 RVO row 'CL (L h-1 1618-g-1)' = 41.0 (SE 3.15)
    lvc <- log(5.31); label("Volume of distribution for a 1618 g neonate (L)")  # Table 4 RVO row 'V (L 1618-g-1)' = 5.31 (SE 0.753)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of birth weight on CL (unitless)")  # Equation 1 literal exponent 0.75
    # Methods 2.3: PKPD models used 'the final linear PK structural model with
    # Hill coefficient and PMA50 fixed to values estimated from PK data (Table 2)'.
    ltmat50 <- fixed(log(37.4)); label("Postmenstrual age at 50% of mature CL (weeks)")  # Table 2 'PMA 50 (weeks)' = 37.4, fixed per Methods 2.3
    lhill_mat <- fixed(log(2.67)); label("Hill coefficient of the CL maturation function (unitless)")  # Table 2 'Hill' = 2.67, fixed per Methods 2.3
    vc_eta_scale <- 1.50; label("Scale factor applied to the shared CL random effect on V (unitless)")  # Table 4 RVO row 'Shared BSV scale factor' = 1.50 (SE 0.479)

    # PD parameters: linear model E = E0 + SL * C (Equation 4).
    lrbase_rvo <- log(151); label("Baseline right ventricular output E0 (mL/kg/min)")  # Table 4 RVO row 'E0 (mL kg-1 min-1)' = 151 (SE 12.8)
    lslope_rvo <- log(0.214); label("Slope of RVO on plasma dobutamine concentration ((mL/kg/min) per ug/L)")  # Table 4 RVO row 'SL' = 0.214 (SE 0.067), no BSV

    # Table 4 footnote a: BSV CV = sqrt(omega2) x 100%, so omega2 = CV^2.
    etalcl ~ 0.0729  # Table 4 RVO CL BSV 27% (SE 18.7%), shrinkage 18%; omega2 = 0.27^2
    etalrbase_rvo ~ 0.1681  # Table 4 RVO E0 BSV 41% (SE 22.2%), shrinkage 3%; omega2 = 0.41^2

    propSd <- 0.583; label("Proportional residual error, dobutamine concentration (fraction)")  # Table 4 RVO row 'Pharmacokinetic residual error (proportional)' = 0.583 (SE 0.052)
    propSd_rvo <- 0.184; label("Proportional residual error, RVO (fraction)")  # Table 4 RVO row 'Pharmacodynamic residual error (proportional)' = 0.184 (SE 0.014)
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

    # Linear concentration-effect model driven by the plasma concentration
    # (Equation 4); Table 3 selects the linear form for RVO.
    rbase_rvo <- exp(lrbase_rvo + etalrbase_rvo)
    slope_rvo <- exp(lslope_rvo)
    rvo <- rbase_rvo + slope_rvo * Cc

    Cc ~ prop(propSd)
    rvo ~ prop(propSd_rvo)
  })
}
