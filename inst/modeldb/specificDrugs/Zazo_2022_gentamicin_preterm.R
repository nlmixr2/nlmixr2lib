Zazo_2022_gentamicin_preterm <- function() {
  description <- paste(
    "PBPK (minimal, PhysPK / EcosimPro). Gentamicin in preterm neonates",
    "(gestational age 32-37 weeks, postnatal age < 1 week) after a 30-min",
    "intravenous infusion. Serum (volume of distribution Vd) is connected",
    "to a well-stirred heart chamber perfused at cardiac output, a",
    "flow-limited eliminating kidney (elimination flow = TBW * GFR * 0.79,",
    "i.e. glomerular filtration net of 21% tubular reabsorption) and a",
    "flow-limited rest-of-tissues compartment. Built entirely from",
    "literature values with no fitting; between-subject variability is the",
    "authors' log-normal Monte Carlo over GFR, Vd, cardiac output and kidney",
    "blood flow. CAUTION: the printed model does not reproduce the paper's",
    "own concentrations. Clearance matches the paper's Table 4, but",
    "simulated serum concentrations at 2 h are 30-40% below the",
    "published predictions (Table 4, Figure 2), because the printed serum",
    "Vd caps the peak at dose/Vd. The printed serum equation (Eq. 5) is",
    "also not mass-conserving and is encoded here in corrected form, and",
    "the heart-chamber volume is not reported and is a maintainer",
    "assumption. See the vignette Errata before use. Companion model for",
    "term neonates: Zazo_2022_gentamicin_term."
  )
  reference <- paste(
    "Zazo H, Lagarejos E, Prado-Velasco M, Sanchez-Herrero S, Serna J,",
    "Rueda-Ferreiro A, Martin-Suarez A, Calvo MV, Perez-Blanco JS, Lanao JM.",
    "Physiologically-based pharmacokinetic modelling and dosing evaluation",
    "of gentamicin in neonates using PhysPK.",
    "Front Pharmacol. 2022;13:977372. doi:10.3389/fphar.2022.977372.",
    sep = " "
  )
  vignette <- "Zazo_2022_gentamicin"

  # Eq. 1 heart chamber: a well-stirred blood volume through which the
  # cardiac output recirculates ("the heart is not modelled as a tissue but
  # as the source of cardiac output", source paper section 2.1.2). It holds
  # blood, not myocardium, so the canonical `heart` organ compartment does
  # not apply.
  paper_specific_compartments <- c("heart_chamber")

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "gentamicin", units = "mg", specimen = "serum", verified = TRUE),
    heart_chamber = list(analyte = "gentamicin", units = "mg", specimen = "whole blood", verified = TRUE),
    kidney = list(analyte = "gentamicin", units = "mg", specimen = "tissue", verified = TRUE),
    other = list(analyte = "gentamicin", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight (TBW)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales serum Vd (L/kg), GFR, cardiac output and both tissue blood",
        "flows (all mL/min/kg in source paper Tables 1-2). The organ",
        "volumes of Table 1 are absolute and are not weight-scaled. The",
        "authors' Monte Carlo sampled TBW log-normally with mean 1.73 kg and",
        "SD 0.35 kg (Table 2); the observed preterm cohort weighed",
        "1.96 +/- 0.41 kg (Table 3)."
      ),
      source_name = "TBW"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 31L,
    n_studies = 1L,
    age_range = "postnatal age 2-7 days (mean 3.03 +/- 1.02)",
    ga_range = "gestational age 32-37 weeks (mean 33.7 +/- 1.56)",
    weight_range = "1.16-3.00 kg (mean 1.96 +/- 0.41)",
    sex_female_pct = NA_real_,
    disease_state = paste(
      "Preterm neonates in the first week of life on gentamicin for",
      "suspected Gram-negative infection, monitored by routine TDM."
    ),
    dose_range = paste(
      "6 mg/kg nominal as a 0.5-h IV infusion every 48 h (observed",
      "5.74 +/- 0.80 mg/kg, range 4.01-6.88); serum sampled at 2 and 24 h",
      "after the start of the infusion."
    ),
    regions = "Spain (University Hospital of Salamanca, 1999-2003)",
    notes = paste(
      "Demographics from source paper Table 3. The model was not fitted to",
      "these data: it was built from literature physiology (Tables 1-2) and",
      "validated against the 31 preterm patients (the paper's 47-neonate",
      "cohort, 91 samples in total, split into preterm GA <= 37 weeks and",
      "term GA > 37 weeks). Sex distribution not reported."
    )
  )

  ini({
    # Every value below is a literature physiological value implemented in
    # PhysPK without estimation (source paper section 2.1.3), so all are fixed.
    # Parameters the authors varied in their Monte Carlo are reported as
    # arithmetic mean +/- SD of a log-normal distribution (Table 2); each
    # typical value here is the corresponding log-normal MEDIAN,
    # mean / sqrt(1 + CV^2), and each omega^2 is log(1 + CV^2). The paper
    # confirms this reading for GFR: the Discussion quotes the implemented
    # values as 1.21 (preterm) and 1.59 (term) mL/min/kg, which are exactly the
    # log-normal medians of Table 2's 1.31 +/- 0.53 and 1.72 +/- 0.70.

    # Serum
    lvc <- fixed(log(0.497005)); label("Gentamicin serum volume of distribution Vd (L/kg)") # Table 2 preterm 'Gentamicin volume of distribution' 0.52 +/- 0.16 L/kg (mean +/- SD); median 0.52 / sqrt(1 + (0.16/0.52)^2)

    # Kidney (flow-limited, eliminating; Eqs. 2-3)
    lgfr <- fixed(log(1.214377)); label("Glomerular filtration rate GFR (mL/min/kg)") # Table 2 preterm 'Glomerular filtration rate' 1.31 +/- 0.53 mL/min/kg; median 1.31 / sqrt(1 + (0.53/1.31)^2) = 1.21 as quoted in the Discussion
    lq_kidney <- fixed(log(9.168429)); label("Kidney blood flow Qk (mL/min/kg)") # Table 1 'Blood flow rates' kidney 9.35 mL/min/kg; section 2.1.3 CV 20 percent; median 9.35 / sqrt(1 + 0.2^2)
    lv_kidney <- fixed(log(0.03)); label("Kidney volume Vk (mL)") # Table 1 'Organ volumes (ml)' kidney 0.03, encoded as printed (see vignette Errata on the units)
    lkp_kidney <- fixed(log(10)); label("Kidney:serum partition coefficient Pk (unitless)") # Table 1 'Partition coefficient' kidney 10.0

    # Rest of tissues (flow-limited; Eq. 4)
    lq_other <- fixed(log(217)); label("Rest-of-tissues blood flow Qr (mL/min/kg)") # Table 1 'Blood flow rates' rest of tissues 217 mL/min/kg; not varied in the Monte Carlo
    lv_other <- fixed(log(0.87)); label("Rest-of-tissues volume Vr (mL)") # Table 1 'Organ volumes (ml)' rest of tissues 0.87, encoded as printed (see vignette Errata on the units)
    lkp_other <- fixed(log(1)); label("Rest-of-tissues:serum partition coefficient Pr (unitless)") # Table 1 'Partition coefficient' rest of tissues 1.00

    # Heart chamber (Eq. 1)
    lqc <- fixed(log(165.007811)); label("Cardiac output CO used as heart flow Qh (mL/min/kg)") # Table 2 'Cardiac output' 172 +/- 50.6 mL/min/kg; median 172 / sqrt(1 + (50.6/172)^2)
    # Not reported anywhere in the source paper (a PhysPK built-in). Assumed by
    # the maintainers; the vignette quantifies the sensitivity (C2h moves by
    # under 5 percent across 0.5-50 mL).
    lv_heart_chamber <- fixed(log(5)); label("Heart chamber volume Vh (mL); maintainer assumption, not in the source") # NOT FROM THE PAPER: maintainer assumption, see vignette Errata

    # Monte Carlo variability (section 2.1.3 and Table 2), log-normal
    etalvc ~ fixed(0.0904571) # Table 2 preterm Vd SD 0.16 / mean 0.52 = CV 30.8 percent; log(1 + CV^2)
    etalgfr ~ fixed(0.1515918) # Table 2 preterm GFR SD 0.53 / mean 1.31 = CV 40.5 percent; log(1 + CV^2)
    etalq_kidney ~ fixed(0.0392207) # Section 2.1.3: Qk CV 20 percent; log(1 + 0.2^2)
    etalqc ~ fixed(0.0830033) # Table 2 CO SD 50.6 / mean 172 = CV 29.4 percent (section 2.1.3 text says 24 percent; see vignette Errata); log(1 + CV^2)

    # No residual error model: the paper's Monte Carlo varied physiology only
    propSd <- fixed(0); label("Proportional residual error (fraction); not part of the source model") # Not reported: the validation Monte Carlo had no residual error
  })
  model({
    # Individual physiology. Flows mL/min/kg -> L/h (x WT x 60 / 1000);
    # organ volumes mL -> L (/ 1000); serum Vd L/kg -> L (x WT).
    vc <- exp(lvc + etalvc) * WT
    gfr <- exp(lgfr + etalgfr)
    q_kidney <- exp(lq_kidney + etalq_kidney) * WT * 60 / 1000
    q_other <- exp(lq_other) * WT * 60 / 1000
    qc <- exp(lqc + etalqc) * WT * 60 / 1000
    v_kidney <- exp(lv_kidney) / 1000
    v_other <- exp(lv_other) / 1000
    v_heart_chamber <- exp(lv_heart_chamber) / 1000
    kp_kidney <- exp(lkp_kidney)
    kp_other <- exp(lkp_other)

    # Eq. 3: elimination flow, GFR x TBW net of 21% tubular reabsorption
    # (Contrepois 1985); L/h
    qe <- WT * gfr * 0.79 * 60 / 1000

    # Steady-state systemic clearance of a flow-limited eliminating kidney,
    # the quantity the paper compares with MAP-Bayesian CL in Table 4 (L/h)
    cl <- q_kidney * qe / (q_kidney + qe)

    Cc <- central / vc
    c_heart_chamber <- heart_chamber / v_heart_chamber
    c_kidney <- kidney / v_kidney
    c_other <- other / v_other

    # Eq. 1
    d/dt(heart_chamber) <- qc * Cc - qc * c_heart_chamber
    # Eq. 2
    d/dt(kidney) <- q_kidney * Cc - q_kidney * c_kidney / kp_kidney - qe * c_kidney / kp_kidney
    # Eq. 4
    d/dt(other) <- q_other * Cc - q_other * c_other / kp_other
    # Eq. 5, corrected: the printed equation omits the -(Qk + Qr) * Cb outflow
    # to the tissues that Eqs. 2 and 4 take up, so as printed it creates mass.
    # The 30-min infusion G * P(t) is supplied by the dosing record (dur = 0.5).
    d/dt(central) <- qc * c_heart_chamber + q_kidney * c_kidney / kp_kidney + q_other * c_other / kp_other - qc * Cc - (q_kidney + q_other) * Cc

    Cc ~ prop(propSd)
  })
}
