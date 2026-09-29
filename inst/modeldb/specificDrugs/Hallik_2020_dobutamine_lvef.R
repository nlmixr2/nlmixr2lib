Hallik_2020_dobutamine_lvef <- function() {
  description <- "Simultaneous PKPD model for intravenous dobutamine and left ventricular ejection fraction (LVEF, echocardiography) in critically ill preterm and term neonates in the first 3 days of life (Hallik 2020, Table 4). PK is the paper's one-compartment linear model re-estimated jointly with LVEF: clearance scales with birth weight^0.75 and a sigmoidal postmenstrual-age maturation function (PMA50 and Hill fixed from the PK-only fit), volume scales linearly with birth weight, both referenced to 1618 g, and one shared random effect enters volume multiplied by an estimated scale factor. LVEF rises linearly with the plasma concentration from an estimated baseline of 63.5% (slope 0.0285 percentage points per ug/L), with no effect-compartment delay. Residual errors are proportional for both concentration and LVEF. The PK-only model is Hallik_2020_dobutamine."
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
  # of R * WT * 60 ug/h. Concentrations are ug/L; LVEF is in percent.
  units <- list(time = "h", dosing = "ug", concentration = "ug/L", lvef = "%")

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
    disease_state = "Critically ill preterm and term neonates in the first 72 h of life needing inotropic support on clinical grounds. Baseline LVEF median 64% (51-79), Table 1.",
    dose_range = "Continuous IV infusion started at 5 ug/kg/min and raised by 5 ug/kg/min roughly every 30 min to a maximum of 20 ug/kg/min",
    regions = "Estonia (Tallinn Children's Hospital and Tartu University Hospital)",
    notes = "Prospective 2-centre study (EU CTR 2015-004836-36); 28 of 31 recruited neonates analysed. LVEF was measured by echocardiography before the infusion and about 20-30 min after each dose escalation (Methods 2.2). See Hallik_2020_dobutamine for the full demographic summary."
  )

  ini({
    # PK parameters re-estimated in the joint fit (Table 4, 'PKPD model for
    # LVEF effect'); CL is the fully mature value for a 1618 g neonate.
    lcl <- log(41.2); label("Clearance for a 1618 g neonate at full maturation (L/h)")  # Table 4 LVEF row 'CL (L h-1 1618-g-1)' = 41.2 (SE 3.22)
    lvc <- log(5.26); label("Volume of distribution for a 1618 g neonate (L)")  # Table 4 LVEF row 'V (L 1618-g-1)' = 5.26 (SE 0.753)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of birth weight on CL (unitless)")  # Equation 1 literal exponent 0.75
    # Methods 2.3: PKPD models used 'the final linear PK structural model with
    # Hill coefficient and PMA50 fixed to values estimated from PK data (Table 2)'.
    ltmat50 <- fixed(log(37.4)); label("Postmenstrual age at 50% of mature CL (weeks)")  # Table 2 'PMA 50 (weeks)' = 37.4, fixed per Methods 2.3
    lhill_mat <- fixed(log(2.67)); label("Hill coefficient of the CL maturation function (unitless)")  # Table 2 'Hill' = 2.67, fixed per Methods 2.3
    vc_eta_scale <- 1.38; label("Scale factor applied to the shared CL random effect on V (unitless)")  # Table 4 LVEF row 'Shared BSV scale factor' = 1.38 (SE 0.434)

    # PD parameters: linear model E = E0 + SL * C (Equation 4).
    lrbase_lvef <- log(63.5); label("Baseline left ventricular ejection fraction E0 (%)")  # Table 4 LVEF row 'E0 (%)' = 63.5 (SE 1.46)
    lslope_lvef <- log(0.0285); label("Slope of LVEF on plasma dobutamine concentration (% per ug/L)")  # Table 4 LVEF row 'SL' = 0.0285 (SE 0.0145), no BSV

    # Table 4 footnote a: BSV CV = sqrt(omega2) x 100%, so omega2 = CV^2.
    etalcl ~ 0.0784  # Table 4 LVEF CL BSV 28% (SE 19.3%), shrinkage 17%; omega2 = 0.28^2
    etalrbase_lvef ~ 0.0081  # Table 4 LVEF E0 BSV 9% (SE 5.6%), shrinkage 10%; omega2 = 0.09^2

    propSd <- 0.580; label("Proportional residual error, dobutamine concentration (fraction)")  # Table 4 LVEF row 'Pharmacokinetic residual error (proportional)' = 0.580 (SE 0.051)
    propSd_lvef <- 0.098; label("Proportional residual error, LVEF (fraction)")  # Table 4 LVEF row 'Pharmacodynamic residual error (proportional)' = 0.098 (SE 0.007)
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
    # (Equation 4); Table 3 selects the linear form for LVEF.
    rbase_lvef <- exp(lrbase_lvef + etalrbase_lvef)
    slope_lvef <- exp(lslope_lvef)
    lvef <- rbase_lvef + slope_lvef * Cc

    Cc ~ prop(propSd)
    lvef ~ prop(propSd_lvef)
  })
}
