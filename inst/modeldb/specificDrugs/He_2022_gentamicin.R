He_2022_gentamicin <- function() {
  description <- "One-compartment IV population PK model for gentamicin in 14 critically ill adults on continuous renal replacement therapy (CVVH or CVVHDF), pooled from two published studies (He 2022). Total clearance is the sum of an estimated endogenous (body) clearance of 1.20 L/h, with log-normal IIV, and the individually calculated CRRT clearance supplied as the data column QEFF; volume of distribution 27.6 L with log-normal IIV; combined additive and proportional residual error. Body weight, age, sex and CRRT modality were screened but not retained."
  reference <- "He S, Cheng Z, Xie F. Population pharmacokinetics and dosing optimization of gentamicin in critically ill patients undergoing continuous renal replacement therapy. Drug Des Devel Ther. 2022;16:13-22. doi:10.2147/DDDT.S343385"
  vignette <- "He_2022_gentamicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "gentamicin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    QEFF = list(
      description = "Individually calculated CRRT (CVVH or CVVHDF) clearance of gentamicin",
      units = "L/h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "He 2022 Methods 'Calculation of CRRT Clearance and Exploratory Data Analysis', Equations 2-8:",
        "CL_CRRT is computed from the circuit settings (effluent flow QE = Qf + Qd in mL/h, times the",
        "sieving coefficient Sc for CVVH or the saturation coefficient Sd for CVVHDF, divided by 1000;",
        "Equations 3 and 6 add a Qb / (Qb + Qrep) pre-dilution factor). It is supplied to the model as a",
        "per-subject data column and ADDED to the estimated endogenous clearance (Equation 1; Methods,",
        "'Population PK Model Development': 'CLbody*EXP(ETA(1)) + CLCRRT'). Every Table 1 value is",
        "reproduced to its printed precision by the post-dilution form CRRT dose [mL/kg/h] * WT [kg] *",
        "Sc-or-Sd / 1000, e.g. patient 1: 45 * 80 * 0.75 / 1000 = 2.70 L/h. Table 1 median 2.6 L/h",
        "(range 1.6-3.5). The paper's Monte Carlo simulations fix Sd at the observed median 0.77 and set",
        "the CRRT dose to 30, 40 or 50 mL/kg/h, i.e. QEFF = 0.77 * dose * WT / 1000. Convert a",
        "weight-normalized prescription to L/h on ingestion (multiply by WT and by the coefficient,",
        "divide by 1000)."
      ),
      source_name = "CLCRRT"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "He 2022 Results/Abstract: body weight screened, 'did not influence the PK parameters of gentamicin based on our dataset'. Table 1 median 72 kg (range 48-102). Weight nevertheless matters OUTSIDE the model: gentamicin doses are mg/kg and the CRRT dose is mL/kg/h, so weight sets both the dose and QEFF in the paper's simulations."
    ),
    AGE = list(
      description = "Age",
      units = "year",
      type = "continuous",
      notes = "He 2022 Results: screened, not retained. Table 1 median 68.5 years (range 34-79)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = NA_character_,
      type = "categorical",
      notes = "He 2022 Results: gender screened, not retained. Table 1: 9 male, 5 female."
    ),
    RRT_CVVHDF_STATUS = list(
      description = "CVVHDF (versus CVVH) continuous-RRT modality indicator",
      units = NA_character_,
      type = "categorical",
      notes = "He 2022 Results: CRRT modality screened, not retained. Table 1: 7 patients on CVVH (all from Petejova et al; the reference level in this paper) and 7 on CVVHDF (all from D'Arcy et al). The CVVHDF saturation coefficient was not measured in D'Arcy et al and was imputed as 0.77, the median of the Petejova sieving coefficients (Table 1 footnote a). The modality enters the model only through how QEFF is computed."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 14L,
    n_studies = 2L,
    n_concentrations = 151L,
    age_median = "68.5 years (range 34-79)",
    weight_median = "72 kg (range 48-102)",
    sex_female_pct = 35.7,
    disease_state = "Critically ill adults with acute kidney injury and septic shock receiving continuous renal replacement therapy and gentamicin. Median APACHE II 28.5 (range 24-42).",
    renal_function = "All on CRRT: 7 CVVH, 7 CVVHDF; blood flow 200 mL/min in every patient; CRRT dose 36.6-67.7 mL/kg/h (median 45); CRRT clearance 1.6-3.5 L/h (median 2.6). Residual endogenous clearance was judged negligible (CL_body / CL_total < 0.05) in 4 of 14 patients.",
    dose_range = "Gentamicin 2.4-5 mg/kg q24h as a 30-minute intravenous infusion (Petejova et al 2.4-3.3 mg/kg; D'Arcy et al 5 mg/kg).",
    regions = "Pooled individual concentration data from two published studies: Petejova et al (n = 7, CVVH) and D'Arcy et al (n = 7, CVVHDF); countries not stated in He 2022.",
    notes = "Data from a PubMed literature search (to 14 May 2021) for gentamicin PK studies in CRRT with individual concentrations and full dosing / CRRT settings; 2 of 4 candidate studies qualified. NONMEM 7.3 FOCE-I via PsN 5.0.0 / Pirana 2.9.9. Model evaluated by 1000-sample bootstrap (Table 2) and pcVPC. Dosing simulations used 10,000 virtual patients drawn from MIMIC-III, weight restricted to 48-102 kg, Sd fixed at 0.77, 3-7 mg/kg q24h and CRRT dose 30-50 mL/kg/h."
  )

  ini({
    # Structural parameters: He 2022 Table 2, 'Final PK Model Estimate (RSE%)'
    # column. Equation 1: CL_total = CL_body + CL_CRRT, with the individual
    # clearance written as CLbody*EXP(ETA(1)) + CLCRRT (Methods, 'Population
    # PK Model Development'). CL_CRRT is calculated per subject (Equations 2-8)
    # and enters as the data column QEFF -- it is not estimated.
    # Table 2: CL_body = 1.20 L/h (RSE 20.3%; bootstrap median 1.29, 95% CI 0.89-1.93)
    lcl <- log(1.20); label("Endogenous (body) clearance (L/h)")
    # Table 2: V_d = 27.60 L (RSE 5.8%; bootstrap median 27.13, 95% CI 24.13-30.77)
    lvc <- log(27.60); label("Volume of distribution (L)")

    # Inter-individual variability, Table 2 reported as CV% with the footnote
    # CV(%) = sqrt(exp(omega^2) - 1) * 100, so omega^2 = log(CV^2 + 1).
    # Table 2: IIV CL_body 69.3 %CV (RSE 35.6%, shrinkage 17.8%)
    etalcl ~ 0.392211 # log(0.693^2 + 1)
    # Table 2: IIV V_d 22.1 %CV (RSE 37.0%, shrinkage 0.1%)
    etalvc ~ 0.047686 # log(0.221^2 + 1)

    # Combined additive + proportional residual error (Results, 'Population
    # PK': 'the additive and proportional errors were 0.156 mg/L and 8.0%').
    # Table 2: proportional error 8.0% (RSE 50.8%, shrinkage 7.5%)
    propSd <- 0.080; label("Proportional residual error (fraction)")
    # Table 2: additive error 0.156 mg/L (RSE 30.6%, shrinkage 7.5%)
    addSd <- 0.156; label("Additive residual error (mg/L)")
  })
  model({
    # Endogenous ("body") clearance: the estimated arm of Equation 1, with the
    # log-normal random effect on this arm only.
    cl_body <- exp(lcl + etalcl)

    # Total clearance, Equation 1. QEFF is the per-subject CRRT clearance in
    # L/h. The TOTAL must be the variable named `cl`: rxode2 recognises a
    # `cl` / `vc` pair and solves the one-compartment system analytically
    # from those two variables, ignoring the explicit d/dt(central) below, so
    # naming the endogenous arm `cl` would silently drop QEFF from the solve.
    cl <- cl_body + QEFF
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # One-compartment model with linear elimination (Results, 'Population
    # PK': a two-compartment model 'does not significantly improve the model
    # fit').
    d/dt(central) <- -kel * central

    # Dose in mg and vc in L give mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
