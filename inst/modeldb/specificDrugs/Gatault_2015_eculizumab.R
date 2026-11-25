Gatault_2015_eculizumab <- function() {
  description <- "One-compartment population PK model of eculizumab with parallel linear and Michaelis-Menten elimination and a body-weight power effect on clearance, in adults on maintenance therapy for atypical hemolytic-uremic syndrome or paroxysmal nocturnal hemoglobinuria, with the in vitro exponential relationship between free eculizumab and terminal complement activity (Gatault 2015)"
  reference <- "Gatault P, Brachet G, Ternant D, Degenne D, Recipon G, Barbet C, Gyan E, Gouilleux-Gruart V, Bordes C, Farrell A, Halimi JM, Watier H. Therapeutic drug monitoring of eculizumab: Rationale for an individualized dosing schedule. MAbs. 2015;7(6):1205-1211. doi:10.1080/19420862.2015.1086049"
  vignette <- "Gatault_2015_eculizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL centred on the median weight of the nine patients, 63 kg (Table 2 weights 50-92 kg; the Results also call 63 kg the 'medium-weight' patient). Methods print the covariate model as 'ln(theta_i) = ln(theta) = beta_WT . (WT/med(WT))', which is garbled; it is read as the Monolix log-transformed continuous covariate, ln(CL_i) = ln(CL) + beta_WT * ln(WT/63), i.e. CL = 0.18 * (WT/63)^2.5. Back-solving the typical CL at each weight from the authors' simulated Figure 4 gives a log-log slope of 2.44-2.56 across the three dosing intervals, confirming the power form (a linear-exponential exp(beta * (WT/63 - 1)) fits that figure about 90-fold worse in sum of squared log errors). Weight enters CL only; it was the only covariate tested. One weight per patient (Table 2).",
      source_name = "WT"
    )
  )

  compartmentData <- list(
    central = list(analyte = "eculizumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 9L,
    n_studies = 1L,
    age_range = "23.9-67.0 years",
    age_median = "40.9 years",
    weight_range = "50-92 kg",
    weight_median = "63 kg",
    sex_female_pct = 77.8,
    race_ethnicity = "All Caucasian",
    disease_state = "Seven adults with atypical hemolytic-uremic syndrome (aHUS; two of them kidney transplant recipients) and two with paroxysmal nocturnal hemoglobinuria (PNH), all on eculizumab maintenance therapy.",
    dose_range = "aHUS: 900 mg IV weekly for 4 weeks, 1200 mg at week 5, then 1200 mg every 2 weeks. PNH: 600 mg IV weekly for 4 weeks, 900 mg at week 5, then 900 mg every 2 weeks.",
    regions = "France (Tours University Hospital, single centre)",
    notes = "Retrospective therapeutic-drug-monitoring cohort (Gatault 2015 Table 1 demographics; Table 2 weights and trough summaries). 97 stored pre-infusion serum samples assayed for free eculizumab by an in-house C5-capture ELISA, plus a dense inter-infusion profile in Patient 6 and post-discontinuation samples in Patient 2. Fitted by SAEM in Monolix 4.3.2. The terminal-complement-activity relationship was fitted separately, by regression on seven spiked pooled healthy-donor sera (20-166 ug/mL)."
  )

  ini({
    # Structural PK parameters, final model (Gatault 2015 Results,
    # 'Eculizumab clearance'), typical values at the 63-kg median weight.
    lvc <- log(3.2); label("Volume of distribution V (L)") # Results: V = 3.2 L (RSE 10%)
    lcl <- log(0.18); label("Linear clearance CL at 63 kg (L/day)") # Results: CL = 0.18 L/day (RSE 7%)
    lvmax <- log(5.0); label("Maximum Michaelis-Menten elimination rate Vm (mg/day)") # Results: Vm = 5.0 (RSE 4%); printed as mg/L/day, encoded as an amount rate (Monolix library convention; see vignette)
    lkm <- log(1.1); label("Michaelis-Menten constant Km (ug/mL)") # Results: Km = 1.1 mg/L (RSE 33%)

    # Body-weight effect on CL (Results: beta_WT = 2.5, p < 0.0001)
    e_wt_cl <- 2.5; label("Power exponent of body weight on CL, reference 63 kg (unitless)") # Results: beta_WT = 2.5

    # Inter-individual variability, exponential model; Monolix reports omega
    # as a standard deviation (Results: 15.0% for V, 14.6% for CL).
    # IIV on Vm and Km was fixed to 0 (Results; Supplement 'Results for the
    # base model'), so those parameters carry no eta.
    etalvc ~ 0.0225 # Results: omega V = 0.150 (SD), variance 0.150^2
    etalcl ~ 0.021316 # Results: omega CL = 0.146 (SD), variance 0.146^2

    # Residual error: proportional (Supplement: 'The best residual model was proportional')
    propSd <- 0.25; label("Proportional residual error (fraction)") # Results: proportional SD = 0.25 (RSE 1.6%)

    # Terminal complement activity (C5b-9 / C9 neoantigen ELISA, arbitrary
    # units) versus free eculizumab concentration, fitted in vitro on spiked
    # sera: A = beta0 * exp(beta1 * C) with A = 20.6 * exp(-0.084 * C)
    # (Methods, 'Measurement of eculizumab trough levels'; Figure 5).
    lrbase_tca <- log(20.6); label("Terminal complement activity at zero free eculizumab, beta0 (arbitrary units)") # Methods: beta0 = 20.6
    lslope_tca <- log(0.084); label("Exponential decline rate of terminal complement activity with free eculizumab, -beta1 (mL/ug)") # Methods: beta1 = -0.084
  })

  model({
    WT_ref <- 63

    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl) * (WT / WT_ref)^e_wt_cl
    vmax <- exp(lvmax)
    km <- exp(lkm)
    kel <- cl / vc

    Cc <- central / vc

    d/dt(central) <- -kel * central - vmax * Cc / (km + Cc)

    # In vitro concentration-effect relationship applied to the predicted
    # free concentration, as in the authors' Figure 4 projections.
    rbase_tca <- exp(lrbase_tca)
    slope_tca <- exp(lslope_tca)
    tca <- rbase_tca * exp(-slope_tca * Cc)

    Cc ~ prop(propSd)
  })
}
