Rantanen_2026_methylprednisolone <- function() {
  description <- "One-compartment population PK model for intravenous methylprednisolone in neonates and infants undergoing cardiac surgery with cardiopulmonary bypass (Rantanen 2026). Fixed allometric body-weight scaling (exponent 0.75 on CL, 1 on V) centred on 5.7 kg, an empirical power effect of postmenstrual age on CL centred on 424 days, log-additive operation-group effects on V (tetralogy of Fallot, VSD/AVSD and Glenn cohorts versus the neonatal cohort), and a proportional residual error estimated separately for each operation group."
  reference <- "Rantanen M, Keski-Nisula J, Olkkola KT, Neuvonen PJ, Pesonen E, Suominen P, Saari TI. Population pharmacokinetics of methylprednisolone in neonates and paediatric patients undergoing cardiac surgery with cardiopulmonary bypass. Eur J Clin Pharmacol. 2026;82(9):235. doi:10.1007/s00228-026-04168-7"
  vignette <- "Rantanen_2026_methylprednisolone"
  units <- list(time = "min", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed at surgery. Fixed allometric scaling, (WT / 5.7)^0.75 on CL and (WT / 5.7)^1 on V, centred on the 5.7 kg study-population median (Methods 'Covariate analysis'; Table 2 'Fixed scaling parameters'). Cohort mean (SD) 5.8 (1.9) kg; neonatal cohort 3.6 (0.41) kg (Table 1).",
      source_name = "WT"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = "Canonical PAGE is carried in months; the source expresses PMA in DAYS and centres it on the 424-day median, so model() recovers days as PAGE * 30.4375 (the Bijleveld_2017_gentamicin convention) before forming (PMA_days / 424)^0.503 on CL. For the tetralogy of Fallot, VSD/AVSD and Glenn cohorts the authors assumed a gestational age of 280 days when calculating PMA (Table 1 footnote a), i.e. PMA_days = 280 + postnatal age in days. Observed range 270-854 days (8.9-28.1 months).",
      source_name = "PMA"
    ),
    SURG_TOF = list(
      description = "Tetralogy of Fallot total-correction cohort indicator (operation group OPER 1)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (neonatal cardiac-surgery cohort, OPER 0, when all three SURG_ indicators are 0)",
      notes = "Time-fixed. 1 = infant undergoing total correction of tetralogy of Fallot (Keski-Nisula 2016 study, n = 28). Multiplies V by exp(-0.404) = 0.668 and selects the OPER 1 proportional residual error. SURG_TOF, SURG_VSD_AVSD and SURG_GLENN are mutually exclusive. Each operation group is also one source study, so the authors caution that the group effects are empirical descriptors of the pooled data, combining developmental, anatomical and perioperative differences, not isolated effects of the procedure (Discussion).",
      source_name = "OPER = 1"
    ),
    SURG_VSD_AVSD = list(
      description = "Ventricular or atrioventricular septal defect repair cohort indicator (operation group OPER 2)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (neonatal cardiac-surgery cohort, OPER 0, when all three SURG_ indicators are 0)",
      notes = "Time-fixed. 1 = infant undergoing VSD or AVSD repair (Keski-Nisula 2015 study, n = 30). Multiplies V by exp(-0.475) = 0.622 and selects the OPER 2 proportional residual error. Half of this cohort received the dose into the CPB-circuit priming solution rather than intravenously; the source model treats both routes as an intravenous bolus (Methods 'Pharmacometric analysis'), and route was screened but not retained.",
      source_name = "OPER = 2"
    ),
    SURG_GLENN = list(
      description = "Bidirectional Glenn procedure cohort indicator (operation group OPER 3)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (neonatal cardiac-surgery cohort, OPER 0, when all three SURG_ indicators are 0)",
      notes = "Time-fixed. 1 = infant with tricuspid atresia undergoing the bidirectional Glenn procedure (Keski-Nisula 2020 study, n = 15). Multiplies V by exp(-0.884) = 0.413 and selects the OPER 3 proportional residual error.",
      source_name = "OPER = 3"
    )
  )

  covariatesDataExcluded <- list(
    NEONATE = list(
      description = "Neonatal status indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened; removing it from the model raised the OFV by only 3.18, so the simpler model without it was accepted (Results 'Population pharmacokinetic model'). Neonatal status is collinear with the reference operation group (OPER 0)."
    ),
    T_CPB = list(
      description = "Duration of cardiopulmonary bypass",
      units = "min",
      type = "continuous",
      notes = "Screened as a centred continuous covariate (Methods 'Covariate analysis'); not retained in the final model."
    ),
    TEMP_MIN_INTRAOP = list(
      description = "Lowest intraoperative body temperature",
      units = "degrees C",
      type = "continuous",
      notes = "Screened; did not improve the model and was not retained (Results 'Population pharmacokinetic model')."
    ),
    T_CIRC_ARREST = list(
      description = "Duration of circulatory arrest",
      units = "min",
      type = "continuous",
      notes = "Screened; did not improve the model and was not retained (Results 'Population pharmacokinetic model')."
    ),
    ROUTE_CPB = list(
      description = "Methylprednisolone given into the CPB-circuit priming solution rather than intravenously",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Methods 'Covariate analysis'); not retained. Both routes are modelled as an intravenous bolus."
    )
  )

  compartmentData <- list(
    central = list(analyte = "methylprednisolone", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 93L,
    n_studies = 4L,
    n_observations = 278L,
    age_range = "postnatal age 1-574 days (median 140 days)",
    age_median = "140 days postnatal",
    page_range = "postmenstrual age 270-854 days (median 424 days)",
    weight_range = "neonatal cohort mean 3.6 (SD 0.41) kg; tetralogy of Fallot cohort 7.0 (SD 1.9) kg",
    weight_median = "5.7 kg (mean 5.8, SD 1.9 kg)",
    sex_female_pct = 57,
    disease_state = "Congenital heart disease in neonates (n = 20) and infants undergoing open-heart surgery with cardiopulmonary bypass: total correction of tetralogy of Fallot (n = 28), ventricular or atrioventricular septal defect repair (n = 30), and the bidirectional Glenn procedure for tricuspid atresia (n = 15).",
    dose_range = "Single intravenous bolus of 30 mg/kg after anaesthesia induction (neonatal and Glenn cohorts); 5 or 30 mg/kg randomized (tetralogy of Fallot cohort); 30 mg/kg either intravenously or into the CPB-circuit priming solution, randomized (VSD/AVSD cohort). Mean (SD) dose 159 (58) mg.",
    regions = "Finland (single-country trials approved by the Ethics Committee of Helsinki University Central Hospital).",
    notes = "Pooled from four prospective trials (Keski-Nisula 2013, 2015, 2016, 2020; Table 1). Arterial samples at 30 min after CPB initiation, 5 min after protamine and 6 h after CPB cessation (median 3 per patient); sampling 23-758 min post-dose. Total plasma methylprednisolone by LC-MS/MS, LLOQ 2 ng/mL; one BLQ sample excluded. Sex 53 female / 40 male (Table 1)."
  )

  ini({
    # Structural parameters at the 5.7 kg / PMA 424 day reference, neonatal
    # operation group (OPER 0). The source fits time in minutes, so CL is in
    # L/min (Methods 'Pharmacometric analysis'); 0.0761 L/min = 4.57 L/h.
    lcl <- log(0.0761); label("Clearance at 5.7 kg and PMA 424 days (L/min)") # Table 2 THETA(1) = 0.0761 L/min (RSE 6.38%)
    lvc <- log(27.1); label("Volume of distribution at 5.7 kg, neonatal cohort (L)") # Table 2 THETA(2) = 27.1 L (RSE 7.90%)

    # Fixed allometric exponents
    e_wt_cl <- fixed(0.75); label("Allometric body-weight exponent on CL (unitless)") # Table 2 'Fixed scaling parameters', 0.750 (fixed)
    e_wt_vc <- fixed(1); label("Allometric body-weight exponent on V (unitless)") # Table 2 'Fixed scaling parameters', 1.00 (fixed)

    # Covariate effects
    e_page_cl <- 0.503; label("Power exponent of postmenstrual age on CL (unitless)") # Table 2 THETA(3) = 0.503 (RSE 45.0%)
    e_surg_tof_vc <- -0.404; label("Log-scale effect of the tetralogy of Fallot cohort on V (unitless)") # Table 2 THETA(4) = -0.404 (V ratio 0.668)
    e_surg_vsd_avsd_vc <- -0.475; label("Log-scale effect of the VSD/AVSD cohort on V (unitless)") # Table 2 THETA(5) = -0.475 (V ratio 0.622)
    e_surg_glenn_vc <- -0.884; label("Log-scale effect of the Glenn cohort on V (unitless)") # Table 2 THETA(6) = -0.884 (V ratio 0.413)

    # Interindividual variability (variances; Table 2 footnote 2 gives
    # CV% = 100 * sqrt(exp(omega^2) - 1))
    etalcl ~ 0.208 # Table 2 OMEGA(1,1) = 0.208 (CV 48.1%)
    etalvc ~ 0.0388 # Table 2 OMEGA(2,2) = 0.0388 (CV 19.9%)

    # Proportional residual error, one per operation group. Table 2 reports
    # the SIGMA variances; the SD is sqrt(variance) (Table 2 footnote 2).
    propSd_neonate <- 0.2917; label("Proportional residual error, neonatal cohort OPER 0 (fraction)") # Table 2 SIGMA(1,1) = 0.0851 -> sqrt = 0.2917 (29.2%)
    propSd_tof <- 0.4231; label("Proportional residual error, tetralogy of Fallot cohort OPER 1 (fraction)") # Table 2 SIGMA(2,2) = 0.179 -> sqrt = 0.4231 (42.3%)
    propSd_vsd_avsd <- 0.2069; label("Proportional residual error, VSD/AVSD cohort OPER 2 (fraction)") # Table 2 SIGMA(3,3) = 0.0428 -> sqrt = 0.2069 (20.7%)
    propSd_glenn <- 0.1404; label("Proportional residual error, Glenn cohort OPER 3 (fraction)") # Table 2 SIGMA(4,4) = 0.0197 -> sqrt = 0.1404 (14.0%)
  })

  model({
    # Postmenstrual age in days (canonical PAGE is in months)
    pma_days <- PAGE * 30.4375

    cl <- exp(lcl + etalcl) * (WT / 5.7)^e_wt_cl * (pma_days / 424)^e_page_cl
    vc <- exp(
      lvc +
        etalvc +
        e_surg_tof_vc * SURG_TOF +
        e_surg_vsd_avsd_vc * SURG_VSD_AVSD +
        e_surg_glenn_vc * SURG_GLENN
    ) *
      (WT / 5.7)^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # mg / L -> ng/mL
    Cc <- 1000 * central / vc

    # Operation-group-specific proportional residual error
    sd_prop <- propSd_neonate *
      (1 - SURG_TOF - SURG_VSD_AVSD - SURG_GLENN) +
      propSd_tof * SURG_TOF +
      propSd_vsd_avsd * SURG_VSD_AVSD +
      propSd_glenn * SURG_GLENN
    Cc ~ prop(sd_prop)
  })
}
