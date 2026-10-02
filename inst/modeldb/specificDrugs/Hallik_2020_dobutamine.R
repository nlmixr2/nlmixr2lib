Hallik_2020_dobutamine <- function() {
  description <- "One-compartment population PK model for intravenous dobutamine in critically ill preterm and term neonates in the first 3 days of life, given as a continuous infusion titrated from 5 to at most 20 ug/kg/min (Hallik 2020, final linear PK model of Table 2). Clearance is allometrically scaled to birth weight with a fixed exponent of 0.75 and multiplied by a sigmoidal postmenstrual-age maturation function (PMA50 = 37.4 weeks, Hill = 2.67); volume scales linearly with birth weight. Both are referenced to the 1618 g cohort median birth weight. Clearance and volume share a single random effect, which enters volume multiplied by an estimated scale factor of 1.34. Residual error is proportional (58.1%). The six simultaneous PKPD models the paper fitted on top of this PK structure (right and left ventricular output, ejection fraction, heart rate, mean arterial pressure, cerebral fractional tissue oxygen extraction) are the companion Hallik_2020_dobutamine_* models."
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
  # of R * WT * 60 ug/h. Concentrations are ug/L (the paper's unit; LLOQ
  # 0.97 ug/L, maximum measured 330 ug/L).
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  compartmentData <- list(
    central = list(analyte = "dobutamine", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT_BIRTH = list(
      description = "Birth weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric size descriptor for both CL (exponent 0.75, fixed) and V (exponent 1), normalised to the population median birth weight of 1618 g (Equations 1-2; Table 1 median 1618 g, range 465-4380 g). The paper enters birth weight in grams; this file takes kilograms and uses a 1.618 kg reference, which is the same ratio. Every infant was studied within the first 72 h of life (age at recruitment median 6 h), so birth weight stands in for current weight.",
      source_name = "BW"
    ),
    PAGE = list(
      description = "Postmenstrual age (gestational age at birth plus postnatal age)",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the sigmoidal maturation of CL, PAGE^Hill / (PMA50^Hill + PAGE^Hill) (Equation 1), with PMA50 = 37.4 weeks and Hill = 2.67 both estimated from the PK data (Table 2). Declared in WEEKS, as the paper states it, rather than the register's default of months, because PMA50 is only meaningful on the weeks scale. In this first-days-of-life cohort PAGE is essentially the gestational age at birth (Table 1 median 30.4 weeks, range 22.7-41.0 weeks).",
      source_name = "PMA"
    )
  )

  # Covariates screened in the PK model but not retained. Methods: 'In
  # covariate analysis parameterization of PK model with postnatal age,
  # antenatal glucocorticoid hormone administration, coadministration of
  # dopamine, blood haemoglobin and albumin concentration, patent ductus
  # arteriosus diameter, baseline LVEF and baseline RVO was tested.' Results:
  # allometry plus maturation 'lowered BSV by 62%, without further improvement
  # by other covariates'. Antenatal glucocorticoids, patent ductus arteriosus
  # diameter, baseline LVEF and baseline RVO were screened the same way; they
  # are not given register entries here because none has a canonical column.
  covariatesDataExcluded <- list(
    PNA = list(
      description = "Postnatal age",
      units = "months",
      type = "continuous",
      notes = "Screened on the PK parameters and not retained (Methods, covariate analysis; Results 3.1). Age at recruitment was median 6 h (range 2-28 h), Table 1."
    ),
    CONMED_DOPAMINE = list(
      description = "Concomitant dopamine indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened on the PK parameters and not retained. 5 of 28 neonates (18%) received dopamine (Table 1)."
    ),
    HGB = list(
      description = "Blood haemoglobin concentration",
      units = "g/L",
      type = "continuous",
      notes = "Screened on the PK parameters and not retained. Median 163 g/L (range 117-203) at recruitment (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/L",
      type = "continuous",
      notes = "Screened on the PK parameters and not retained. Median 27.7 g/L (range 21.6-41.3) at recruitment (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 28,
    n_studies = 1,
    n_observations = 119,
    age_range = "Postnatal age at recruitment 2-28 h (median 6 h); gestational age at birth 22.7-41.0 weeks (median 30.4)",
    weight_range = "Birth weight 465-4380 g (median 1618 g)",
    sex_female_pct = 36,
    disease_state = "Critically ill preterm and term neonates in the first 72 h of life needing inotropic support on clinical grounds (main diagnoses respiratory distress syndrome, early-onset sepsis, perinatal asphyxia, meconium aspiration, foeto-foetal transfusion). Congenital heart disease, hydrops and therapeutic hypothermia excluded.",
    dose_range = "Continuous IV infusion started at 5 ug/kg/min and raised by 5 ug/kg/min roughly every 30 min to a maximum of 20 ug/kg/min; maximal dose 10, 15 and 20 ug/kg/min in 1, 17 and 10 neonates; infusion duration median 3.5 days (1.4-17.7)",
    regions = "Estonia (Tallinn Children's Hospital and Tartu University Hospital)",
    notes = "Prospective 2-centre study, April 2016 to December 2017 (EU CTR 2015-004836-36). 31 recruited, 28 analysed. Gestational-age bands: <28 weeks 7 (25%), 28-32 weeks 9 (32%), 32-37 weeks 7 (25%), >37 weeks 5 (18%). 119 plasma samples; 9 below the LLOQ of 0.97 ug/L were set to 0.5 ug/L and 2 had no detectable drug (Methods 2.1; Results). Table 1."
  )

  ini({
    # Structural parameters (Table 2, 'The linear pharmacokinetic model').
    # CL and V are typical values for a neonate of the median birth weight
    # (1618 g); CL is the fully mature value, before the maturation factor.
    lcl <- log(41.2); label("Clearance for a 1618 g neonate at full maturation (L/h)")  # Table 2 row 'CL (L h-1 1618-g-1)' = 41.2 (SE 44.5)
    lvc <- log(5.29); label("Volume of distribution for a 1618 g neonate (L)")  # Table 2 row 'V (L 1618-g-1)' = 5.29 (SE 0.821); Discussion '5.29 L 1618-g-1 or 3.27 L kg-1'

    # Allometric exponent on CL. Discussion: 'allometric scaling to population
    # median BW with power coefficient of 0.75'; Equation 1 prints the 0.75 as a
    # literal and Table 2 gives it no estimate, so it is fixed. V scales
    # linearly (Equation 2 has no exponent).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of birth weight on CL (unitless)")  # Equation 1 literal exponent 0.75; Discussion 4.1

    # Maturation of CL with postmenstrual age (Equation 1). Methods: 'PMA50
    # and Hill's coefficient for the maturation function of dobutamine CL were
    # estimated from PK data.'
    ltmat50 <- log(37.4); label("Postmenstrual age at 50% of mature CL (weeks)")  # Table 2 row 'PMA 50 (weeks)' = 37.4 (SE 30.6)
    lhill_mat <- log(2.67); label("Hill coefficient of the CL maturation function (unitless)")  # Table 2 row 'Hill' = 2.67 (SE 1.90)

    # Shared random effect. Methods: 'Individual estimates for BSV for CL and V
    # were highly correlated (r = 1.0), so a shared BSV was used with an
    # estimated scale factor applied for V.' The eta sits on CL; V receives it
    # multiplied by the scale factor.
    vc_eta_scale <- 1.34; label("Scale factor applied to the shared CL random effect on V (unitless)")  # Table 2 row 'Shared BSV scale factor' = 1.34 (SE 0.373)

    # Table 2 footnote a: BSV 'presented as coefficient of variation,
    # calculated as: (square root of omega2) x 100%', so omega2 = CV^2.
    etalcl ~ 0.0841  # Table 2 row CL BSV 29% (SE 17.2%), shrinkage 17%; omega2 = 0.29^2 = 0.0841

    # Results 3.1: 'proportional error model for residual variability';
    # Results: residual variability 'remaining >50% in PK observations'.
    propSd <- 0.581; label("Proportional residual error (fraction)")  # Table 2 row 'Residual error (proportional)' = 0.581 (SE 0.229)
  })

  model({
    # Size and maturation (Equations 1-2), reference birth weight 1.618 kg.
    tmat50 <- exp(ltmat50)
    hill_mat <- exp(lhill_mat)
    fmat <- PAGE^hill_mat / (tmat50^hill_mat + PAGE^hill_mat)

    cl <- exp(lcl + etalcl) * (WT_BIRTH / 1.618)^e_wt_cl * fmat
    vc <- exp(lvc + vc_eta_scale * etalcl) * (WT_BIRTH / 1.618)

    kel <- cl / vc

    # One-compartment disposition with first-order elimination; dobutamine is
    # infused directly into the central compartment.
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
