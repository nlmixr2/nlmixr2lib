Abdulla_2020_ciprofloxacin <- function() {
  description <- "Two-compartment population PK model for intravenous ciprofloxacin in adult ICU patients (Abdulla 2020). Linear elimination from the central compartment, IIV on CL and Vc, and a combined additive + proportional residual error. The model was fitted to protein-binding-corrected (unbound) plasma concentrations, so the observation Cc = central / vc is the UNBOUND concentration; total plasma ciprofloxacin is reconstructed as Ctot = Cc / fu with the paper's assumed 30% plasma protein binding (fu = 0.7). Serum creatinine, eGFR, albumin, BMI, weight, sex, renal replacement therapy and age were screened but none was retained."
  reference <- "Abdulla A, Rogouti O, Hunfeld NGM, Endeman H, Dijkstra A, van Gelder T, Muller AE, de Winter BCM, Koch BCP. Population pharmacokinetics and target attainment of ciprofloxacin in critically ill patients. Eur J Clin Pharmacol. 2020;76(7):957-967. doi:10.1007/s00228-020-02873-5"
  vignette <- "Abdulla_2020_ciprofloxacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. central holds UNBOUND drug amount (the fit was on
  # protein-binding-corrected concentrations; see the Cc comment in
  # model()).
  compartmentData <- list(
    central = list(analyte = "ciprofloxacin (unbound)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ciprofloxacin (unbound)", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened in the forward-inclusion / backward-elimination covariate search (Methods 'Covariate analysis'); not retained. Cohort median 90 umol/L (IQR 70-153), Table 1."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (MDRD, BSA-normalized)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Screened as eGFR (MDRD); not retained. Cohort median 58.5 (IQR 32-101), Table 1.",
      source_name = "eGFR"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened; not retained. Cohort median 25 g/L (IQR 22-29), Table 1."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened; not retained. Cohort median 26 (range 17.8-46.3), Table 1."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened; not retained. Cohort median 80 kg (IQR 64-90), Table 1."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained. 17 of 42 female, Table 1."
    ),
    RRT_CRRT_STATUS = list(
      description = "Renal replacement therapy (continuous venovenous haemofiltration)",
      units = "(binary)",
      type = "binary",
      notes = "Screened as renal replacement therapy; not retained. 10 of 42 on CVVH, Table 1.",
      source_name = "RRT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened; not retained. Cohort median 65.5 years (IQR 56-71), Table 1."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 42L,
    n_studies = 1L,
    n_observations = 204L,
    age_median = "65.5 years (IQR 56-71)",
    weight_median = "80 kg (IQR 64-90)",
    sex_female_pct = 40.5,
    disease_state = "Critically ill adult ICU patients receiving intravenous ciprofloxacin (primary diagnosis respiratory 47.6%, sepsis 19.0%, cardiovascular 11.9%, gastrointestinal 11.9%, neurological 4.8%, other 4.8%); median APACHE II 22, SOFA 13; 23.8% on continuous venovenous haemofiltration.",
    dose_range = "400 mg IV q24h (n = 3), q12h (n = 25) or q8h (n = 14), infused over 30-60 min; sampled on day 2 of therapy.",
    regions = "The Netherlands (Erasmus Medical Centre and Maasstad Hospital, Rotterdam; EXPAT study).",
    renal_function = "eGFR (MDRD) median 58.5 mL/min/1.73 m^2 (IQR 32-101); serum creatinine median 90 umol/L (IQR 70-153).",
    notes = "Prospective two-centre observational study, 2016. 204 plasma concentrations (4.9 per patient: pre-dose, 15-30 min after end of infusion, 1 and 3 h after infusion, and pre-next-dose). UPLC-MS/MS, linear 0.04-5.0 mg/L. NONMEM 7.2 FOCE-I. Demographics from Table 1."
  )

  ini({
    # Structural parameters: Table 2 'Final model' column (RSE% in brackets).
    lcl <- log(25.4); label("Clearance of unbound ciprofloxacin (L/h)")                              # Table 2, CL = 25.4 L/h (RSE 11%)
    lvc <- log(91.1); label("Central volume of distribution, unbound-referenced (L)")                # Table 2, Vc = 91.1 L (RSE 13%)
    lvp <- log(164);  label("Peripheral volume of distribution, unbound-referenced (L)")             # Table 2, Vp = 164 L (RSE 15%)
    lq  <- log(91.9); label("Intercompartmental clearance (L/h)")                                     # Table 2, Q = 91.9 L/h (RSE 10%)

    # Unbound fraction: not estimated. Methods 'Blood sampling and assays':
    # observed concentrations were corrected for protein binding using an
    # average plasma protein binding of 30% (fAUC = AUC * 0.7).
    fu <- fixed(0.7); label("Fraction of ciprofloxacin unbound in plasma (unitless)")                # Methods, average PPB 30% -> fu = 0.7

    # IIV reported as CV%; log-normal omega^2 = log(1 + CV^2). The CL-Vc
    # omega block was retained (Results 'Final model') but its covariance is
    # not reported, so the etas are carried as diagonal.
    etalcl ~ 0.378220  # Table 2, IIV CL = 67.8% -> log(1 + 0.678^2)
    etalvc ~ 0.231191  # Table 2, IIV Vc = 51.0% -> log(1 + 0.510^2)

    # Residual error on the unbound concentration (the fitted DV). Table 2
    # prints both rows under 'Residual variability (%)'; the additive row
    # 14.3 is read as SD = 0.143 mg/L (14.3 mg/L would exceed every observed
    # concentration; see vignette 'Assumptions and deviations').
    propSd <- 0.153; label("Proportional residual error on unbound Cc (fraction)")    # Table 2, Proportional = 15.3%
    addSd  <- 0.143; label("Additive residual error on unbound Cc (mg/L)")           # Table 2, Additional = 14.3 (read as 0.143 mg/L)
  })

  model({
    # Individual parameters (no covariates retained; Results 'Covariate analysis').
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # IV infusion into central (rate / dur on the dose record).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # The fitted concentrations were the protein-binding-corrected (unbound)
    # values: the model's AUC = dose / CL reproduces the paper's reported
    # fAUC0-24 (Table 1 and Fig. 4), not the total AUC. Total plasma
    # ciprofloxacin is reconstructed with the paper's fixed fu (same
    # Cc-unbound / Ctot-total convention as Beijer_2026_cloxacillin).
    Cc   <- central / vc
    Ctot <- Cc / fu

    Cc ~ add(addSd) + prop(propSd)
  })
}
