Hallik_2020_dobutamine_ftoe_cerebral <- function() {
  description <- "Simultaneous PKPD model for intravenous dobutamine and cerebral fractional tissue oxygen extraction (cFTOE = (SaO2 - rScO2) / SaO2, from continuous pulse oximetry and near-infrared spectroscopy, averaged over 1 min) in critically ill preterm and term neonates in the first 3 days of life (Hallik 2020, Table 5). PK is the paper's one-compartment linear model re-estimated jointly with cFTOE: clearance scales with birth weight^0.75 and a sigmoidal postmenstrual-age maturation function (PMA50 and Hill fixed from the PK-only fit), volume scales linearly with birth weight, both referenced to 1618 g, and one shared random effect enters volume multiplied by an estimated scale factor. cFTOE follows a sigmoidal Emax model in the plasma concentration (no effect-compartment delay), FALLING from a baseline of 0.227 toward a plateau LEVEL of 0.206 (the paper's Emax is the plateau value, not the increment, so the typical effect is negative) with EC50 52.9 ug/L and Hill 3.65. Residual errors are proportional for both concentration and cFTOE. The PK-only model is Hallik_2020_dobutamine."
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
  # of R * WT * 60 ug/h. Concentrations are ug/L; cFTOE is a unitless fraction.
  units <- list(time = "h", dosing = "ug", concentration = "ug/L", ftoe_cerebral = "fraction")

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
    notes = "Prospective 2-centre study (EU CTR 2015-004836-36); 28 of 31 recruited neonates analysed. SaO2 and cerebral regional oxygen saturation (rScO2) were recorded every 2.5 s; the 15 min around each dose change (5 min before, 10 min after) was averaged over 1-min bins and cFTOE computed as (SaO2 - rScO2) / SaO2 (Methods 2.3). See Hallik_2020_dobutamine for the full demographic summary."
  )

  ini({
    # PK parameters re-estimated in the joint fit (Table 5, 'PKPD model for
    # cFTOE effect'); CL is the fully mature value for a 1618 g neonate.
    lcl <- log(37.2); label("Clearance for a 1618 g neonate at full maturation (L/h)")  # Table 5 cFTOE row 'CL (L h-1 1618-g-1)' = 37.2 (SE 3.61)
    lvc <- log(4.80); label("Volume of distribution for a 1618 g neonate (L)")  # Table 5 cFTOE row 'V (L 1618-g-1)' = 4.80 (SE 0.993)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of birth weight on CL (unitless)")  # Equation 1 literal exponent 0.75
    # Methods 2.3: PKPD models used 'the final linear PK structural model with
    # Hill coefficient and PMA50 fixed to values estimated from PK data (Table 2)'.
    ltmat50 <- fixed(log(37.4)); label("Postmenstrual age at 50% of mature CL (weeks)")  # Table 2 'PMA 50 (weeks)' = 37.4, fixed per Methods 2.3
    lhill_mat <- fixed(log(2.67)); label("Hill coefficient of the CL maturation function (unitless)")  # Table 2 'Hill' = 2.67, fixed per Methods 2.3
    vc_eta_scale <- 2.42; label("Scale factor applied to the shared CL random effect on V (unitless)")  # Table 5 cFTOE row 'Shared BSV scale factor' = 2.42 (SE 0.333)

    # PD parameters: sigmoidal Emax E = E0 + (Emax - E0) * C^g / (EC50^g + C^g)
    # (Equation 6). Emax is 'the estimated maximum HD parameter value', i.e.
    # the plateau LEVEL, hence lrmax_ rather than lemax; here it lies below E0,
    # matching Discussion 4.2 'Decrease in cFTOE with dobutamine'.
    lrbase_ftoe_cerebral <- log(0.227); label("Baseline cerebral fractional tissue oxygen extraction E0 (fraction)")  # Table 5 cFTOE row 'E0' = 0.227 (SE 0.023)
    lrmax_ftoe_cerebral <- log(0.206); label("Plateau cerebral fractional tissue oxygen extraction Emax (fraction)")  # Table 5 cFTOE row 'Emax' = 0.206 (SE 0.027)
    lec50_ftoe_cerebral <- log(52.9); label("Plasma concentration at half-maximal cFTOE change (ug/L)")  # Table 5 cFTOE row 'EC50 (ug L-1)' = 52.9 (SE 7.26), no BSV
    lhill_ftoe_cerebral <- log(3.65); label("Hill coefficient of the cFTOE concentration-effect curve (unitless)")  # Table 5 cFTOE row 'gamma' = 3.65 (SE 0.573), no BSV

    # Table 5 footnote a: BSV CV = sqrt(omega2) x 100%, so omega2 = CV^2.
    etalcl ~ 0.1225  # Table 5 cFTOE CL BSV 35% (SE 21.5%), shrinkage 9%; omega2 = 0.35^2
    etalrbase_ftoe_cerebral ~ 0.25  # Table 5 cFTOE E0 BSV 50% (SE 26.8%), shrinkage 2%; omega2 = 0.50^2
    etalrmax_ftoe_cerebral ~ 0.36  # Table 5 cFTOE Emax BSV 60% (SE 36.2%), shrinkage 9%; omega2 = 0.60^2

    propSd <- 0.653; label("Proportional residual error, dobutamine concentration (fraction)")  # Table 5 cFTOE row 'Pharmacokinetic residual error (proportional)' = 0.653 (SE 0.061)
    propSd_ftoe_cerebral <- 0.181; label("Proportional residual error, cFTOE (fraction)")  # Table 5 cFTOE row 'Pharmacodynamic residual error (proportional)' = 0.181 (SE 0.004)
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
    # selects 'Sigmoidal Emax' without KEO for cFTOE (OFV -5567 vs -5152 with
    # KEO). Because E0 and Emax carry independent etas, an individual's
    # plateau may lie above or below their own baseline.
    rbase_ftoe_cerebral <- exp(lrbase_ftoe_cerebral + etalrbase_ftoe_cerebral)
    rmax_ftoe_cerebral <- exp(lrmax_ftoe_cerebral + etalrmax_ftoe_cerebral)
    ec50_ftoe_cerebral <- exp(lec50_ftoe_cerebral)
    hill_ftoe_cerebral <- exp(lhill_ftoe_cerebral)
    ftoe_cerebral <- rbase_ftoe_cerebral + (rmax_ftoe_cerebral - rbase_ftoe_cerebral) * Cc^hill_ftoe_cerebral / (ec50_ftoe_cerebral^hill_ftoe_cerebral + Cc^hill_ftoe_cerebral)

    Cc ~ prop(propSd)
    ftoe_cerebral ~ prop(propSd_ftoe_cerebral)
  })
}
