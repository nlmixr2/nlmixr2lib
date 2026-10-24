Crass_2025_pegcetacoplan_intravitreal <- function() {
  description <- "One-compartment population PK model for serum pegcetacoplan after intravitreal injection in adults with geographic atrophy or neovascular age-related macular degeneration (Crass 2025). The intravitreal dose enters a vitreous compartment and is absorbed into the systemic (serum) central compartment by a first-order rate ka that is slower than systemic elimination (flip-flop kinetics); systemic bioavailability is assumed to be 1, so CL and Vc are the apparent values CL/F and Vc/F. CL/F carries power effects of age (reference 80 years) and baseline serum complement C3 (reference 1.2 g/L) and a log-linear female-sex effect. Inter-individual variability is on CL/F (Manly / Box-Cox transformed eta) and ka; Vc/F has no random effect. Residual error is additive on the log scale (log-normal). Vitreous concentration is derived from the vitreous amount over an assumed 4 mL vitreous volume and is not observed."
  reference <- paste(
    "Crass RL, Prem K, Gaudreault F, Ribeiro R, Baumal CR, Smith B, Epling D,",
    "Chapel S. Population Pharmacokinetics of Pegcetacoplan in Patients with",
    "Geographic Atrophy or Neovascular Age-related Macular Degeneration.",
    "Ophthalmol Sci. 2025;5(2):100657. doi:10.1016/j.xops.2024.100657"
  )
  vignette <- "Crass_2025_pegcetacoplan_intravitreal"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL/F referenced to 80 years: (AGE / 80)^e_age_cl, written in the paper as theta6 * (log(Age) - log(80 years)) inside the log-scale covariate sum COVCL (Table 2 footnote equations). Analysis-population mean (SD) 79.6 (7.69) years (Table 1).",
      source_name = "Age"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Enters the log-scale covariate sum COVCL as theta8 * FemaleSex, so female CL/F = male CL/F * exp(theta8) = 0.751-fold (Table 2). The typical patient of Table 2 is male. 56.7% female in the analysis population (Table 1).",
      source_name = "FemaleSex"
    ),
    C3_BL = list(
      description = "Baseline serum complement component 3 (C3) concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL/F referenced to 1.2 g/L: (C3_BL / 1.2)^e_c3_bl_cl, written in the paper as theta7 * (log(Baseline C3) - log(1.2 g/L)) (Table 2 footnote equations). Analysis-population mean (SD) 1.20 (0.241) g/L (Table 1). Time-fixed per subject.",
      source_name = "Baseline C3"
    )
  )

  covariatesDataExcluded <- list(
    DIS_NAMD = list(
      description = "Disease type indicator, neovascular age-related macular degeneration (1) versus geographic atrophy (0)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL/F, VC/F and KA (Methods, Covariate Parameterization) but not retained in the final model."
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened on CL/F and VC/F but not retained in the final model."
    ),
    TBILI = list(
      description = "Baseline total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened on CL/F but not retained in the final model."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL/F but not retained in the final model."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened on CL/F but not retained in the final model."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (BSA-normalized)",
      units = "mL/min/1.73m^2",
      type = "continuous",
      notes = "Screened on CL/F as eGFR but not retained in the final model."
    ),
    CONMED_ANTIVEGF = list(
      description = "Concomitant intravitreal anti-VEGF medication (time-varying)",
      units = "(binary)",
      type = "binary",
      notes = "Tested only in a separate post hoc covariate model (anti-VEGF on CL/F, VC/F and KA) that is not the final model. The post hoc model's parameter estimates are not reported; only the 95% CI of the KA multiplier (0.980, 1.37) and the exposure-ratio forest plot (Figure 3: steady-state serum AUC ratio 0.831, Cmax 0.788, Cmin 0.903) are printed."
    ),
    ADA_POS = list(
      description = "Anti-pegcetacoplan peptide antibody or anti-PEG antibody status",
      units = "(binary)",
      type = "binary",
      notes = "Explored post hoc (Figure 3); not part of the final model and no parameter estimates are reported."
    )
  )

  compartmentData <- list(
    vitreous = list(analyte = "pegcetacoplan", units = "mg", specimen = "vitreous", verified = TRUE),
    central = list(analyte = "pegcetacoplan", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 261L,
    n_studies = 4L,
    n_observations = 2064L,
    age_mean = "79.6 years (SD 7.69)",
    sex_female_pct = 56.7,
    race_ethnicity = c(White = 97.3, Black = 0.8, Other_or_missing = 1.9),
    disease_state = "Geographic atrophy secondary to age-related macular degeneration (n = 231, 88.5%) or neovascular age-related macular degeneration (n = 30, 11.5%).",
    dose_range = "Single intravitreal injections of 4, 10 or 20 mg, and repeated intravitreal injections of 15 mg monthly or every other month.",
    renal_function = "Normal 19.5%, mild impairment 42.9%, moderate 31.0%, severe 0.8%, end stage 4.2% (Table 1).",
    c3_mean = "1.20 g/L (SD 0.241)",
    studies = "POT-CP043014 (n = 13, nAMD), APL2-203 (n = 17, nAMD), POT-CP121614 / FILLY (n = 164, GA), APL2-303 / DERBY (n = 67, GA).",
    notes = "Serum PK only; 1581 of 2064 samples (76.6%) were quantifiable and 483 (23.4%) were below the 0.10 ug/mL lower limit of quantification, handled with the M3 method. Baseline demographics from Table 1."
  )

  ini({
    # Structural parameters (Table 2). The 'Estimate' column is the log-scale
    # NONMEM THETA; the 'Transformed Estimate' column is its exponential. The
    # log of the three-significant-figure transformed value is used.
    lka <- log(0.0528); label("Vitreous-to-serum absorption rate constant ka (1/day)") # Table 2 KA: Estimate -2.94, Transformed 0.0528 1/day
    lcl <- log(0.325); label("Apparent systemic clearance from serum CL/F (L/day)") # Table 2 CL/F: Estimate -1.12, Transformed 0.325 L/day
    lvc <- log(1.83); label("Apparent central serum volume Vc/F (L)") # Table 2 VC/F: Estimate 0.603, Transformed 1.83 L

    # Covariate effects on CL/F (Table 2, theta6-theta8; all three log-scale
    # estimates carry a minus sign in the PDF glyph stream, confirmed by the
    # direction of every Figure 2 forest-plot ratio).
    e_age_cl <- -0.559; label("Power exponent of age on CL/F, reference 80 years (unitless)") # Table 2 'Age on CL/F (theta6)' = -0.559
    e_c3_bl_cl <- -0.350; label("Power exponent of baseline serum C3 on CL/F, reference 1.2 g/L (unitless)") # Table 2 'C3 on CL/F (theta7)' = -0.350
    e_sexf_cl <- -0.286; label("Log-scale effect of female sex on CL/F (unitless)") # Table 2 'Female sex on CL/F (theta8)' = -0.286; transformed 0.751

    # Manly (Box-Cox) shape of the CL/F random effect (Table 2 'Lambda' and
    # footnote equation TETACL = (exp(eta_CL/F * lambda) - 1) / lambda).
    boxcox_lcl <- 2.57; label("Manly / Box-Cox shape parameter of the CL/F random effect (unitless)") # Table 2 'Lambda' = 2.57

    # Inter-individual variability (Table 2). No eta on Vc/F (Table 2
    # footnote: VC/Fi = exp(VC/F)).
    etalcl ~ 0.0587 # Table 2 IIV omega^2 CL/F = 0.0587 (before the Manly transform)
    etalka ~ 0.122 # Table 2 IIV omega^2 KA = 0.122 (36.0% CV)

    # Residual error: log-transform-both-sides additive SD (Table 2
    # 'Residual variability, % Log additive' = 0.300, i.e. 30.0%).
    expSd <- 0.300; label("Additive residual SD on the log scale (unitless)") # Table 2 log-additive residual 0.300
  })

  model({
    # Log-scale covariate sum on CL/F (Table 2 footnote, COVCL)
    covcl <- e_age_cl * log(AGE / 80) + e_c3_bl_cl * log(C3_BL / 1.2) + e_sexf_cl * SEXF

    # Manly-transformed CL/F random effect (Table 2 footnote, TETACL)
    etatrcl <- (exp(etalcl * boxcox_lcl) - 1) / boxcox_lcl

    ka <- exp(lka + etalka)
    cl <- exp(lcl + covcl + etatrcl)
    vc <- exp(lvc)

    kel <- cl / vc

    # vitreous: intravitreal dose administration compartment; central: serum
    # (Figure 1). Systemic bioavailability is assumed to be 1 (Methods).
    d / dt(vitreous) <- -ka * vitreous
    d / dt(central) <- ka * vitreous - kel * central

    # Vitreous concentration (ug/mL) over the assumed 4 mL vitreous volume
    # (Methods, Model Development); mg/mL * 1000 = ug/mL. Not observed.
    Cvit <- vitreous / 4 * 1000

    # Serum concentration: mg / L = ug/mL
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
