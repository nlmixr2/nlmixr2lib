Kruizinga_2022_clonazepam <- function() {
  description <- paste(
    "Two-compartment oral population PK model for clonazepam in plasma and saliva in healthy",
    "adults aged 18-30 years given a single 0.5 or 1.0 mg oral solution (Kruizinga 2022).",
    "First-order absorption with a two-class mixture on the absorption rate constant: 75% of",
    "subjects absorb slowly (ka = 1.106 1/h with IIV) and 25% fast (ka fixed to 100 1/h), encoded",
    "by the binary covariate MIX_FAST_ABS. Clearance and inter-compartmental clearance scale",
    "allometrically with weight (exponent 0.75) and both volumes linearly, referenced to 70 kg.",
    "IIV is carried on the slow-class ka, Q and the relative bioavailability (F fixed to 1).",
    "The saliva concentration is the sum of an oral-contamination term and a plasma-driven term.",
    "The contamination term is a 1 mL saliva compartment that receives a fraction (0.033 per",
    "mille) of the dose and empties first-order (kel_saliva = 1.95 1/h). The plasma-driven term is",
    "Cc times a saliva:plasma ratio that rises with plasma concentration by a saturable function",
    "(maximum 0.195, half-saturation 2.581 ug/L). A dose must therefore be given TWICE in the event",
    "table: once to depot and once, with the same amount, to saliva; the model applies the",
    "deposited fraction to the saliva dose through f(saliva).",
    sep = " "
  )
  reference <- paste(
    "Kruizinga MD, Zuiker RGJA, Bergmann KR, Egas AC, Cohen AF, Santen GWE, van Esdonk MJ.",
    "Population pharmacokinetics of clonazepam in saliva and plasma: Steps towards noninvasive",
    "pharmacokinetic studies in vulnerable populations.",
    "Br J Clin Pharmacol. 2022;88(5):2236-2245. doi:10.1111/bcp.15152",
    sep = " "
  )
  vignette <- "Kruizinga_2022_clonazepam"

  # Concentrations are reported in ug/L throughout the paper (Table 2 RatioKM,
  # Figure 1). Dose is in mg and volumes in L, so mg/L is converted to ug/L by
  # a factor of 1000 in model().
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  compartmentData <- list(
    depot = list(
      analyte = "clonazepam",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "clonazepam",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "clonazepam",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    # The paper's 'saliva contamination compartment' (Results 3.1, Equation 1):
    # residue of the oral solution left in the mouth. It is dosed directly and
    # is not connected to central.
    saliva = list(
      analyte = "clonazepam",
      units = "mg",
      specimen = "saliva",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Fixed allometric scaling referenced to 70 kg: exponent 0.75 on CL and Q, 1 on Vc and Vp",
        "(Kruizinga 2022 Methods 2.4 and Table 2 footnote). Cohort mean (SD) 67.8 (8.3) kg (Table 1)."
      ),
      source_name = "WGT"
    ),
    MIX_FAST_ABS = list(
      description = paste(
        "Latent mixture-model class indicator for absorption speed: 1 = fast-absorption",
        "subpopulation (ka fixed to 100 1/h), 0 = slow-absorption subpopulation (ka = 1.106 1/h",
        "with IIV)."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Not a measured covariate: the class is the NONMEM mixture assignment. The estimated",
        "probability of the slow class is 0.75 (Kruizinga 2022 Table 2, 'Prob. slow group', RSE",
        "13.0%), so for population simulation draw MIX_FAST_ABS ~ Bernoulli(0.25) per subject; for",
        "typical-value simulation of the majority phenotype set MIX_FAST_ABS = 0. No covariates",
        "predicted class membership (Results 3.1)."
      ),
      source_name = "mixture subpopulation"
    )
  )

  covariatesDataExcluded <- list(
    ALB = list(
      description = "Serum albumin.",
      units = "g/L",
      type = "continuous",
      notes = paste(
        "Tested as a covariate on the saliva:plasma relationship as a surrogate for unbound plasma",
        "concentration; no correlation was identified (Kruizinga 2022 Discussion)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20,
    n_studies = 1,
    age_range = "18-30 years (mean 22.4, SD 2.8)",
    weight_range = "mean 67.8 kg (SD 8.3)",
    sex_female_pct = 55,
    race_ethnicity = c(Caucasian = 100),
    disease_state = "Healthy volunteers",
    dose_range = "Single oral dose of 0.5 mg (n = 10) or 1.0 mg (n = 10) clonazepam solution (Rivotril) in lemonade",
    regions = "Netherlands (Centre for Human Drug Research, Leiden)",
    notes = paste(
      "Kruizinga 2022 Methods 2.2 and Table 1: 9 of 20 male, height 175.1 (7.3) cm, BMI 22.2",
      "(2.4) kg/m^2. Paired plasma and saliva samples at 0.5, 1, 2, 4, 6, 8, 24 and 48 h post-dose;",
      "160 plasma and 154 analysable saliva samples, none below the LOQ. Saliva was collected",
      "with a SalivaBio Infant Swab after a mouth rinse 10 minutes before sampling."
    )
  )

  ini({
    # Plasma kinetics (Kruizinga 2022 Table 2)
    lka <- log(1.106); label("Absorption rate constant, slow-absorption class (1/h)") # Table 2 'theta ka - slow group' = 1.106 (RSE 12.8%)
    lka_fastabs <- fixed(log(100)); label("Absorption rate constant, fast-absorption class (1/h)") # Results 3.1: fast absorption population had a fixed ka of 100/h
    lcl <- log(2.98); label("Clearance at 70 kg (L/h)") # Table 2 'theta clearance' = 2.98 (RSE 5.5%)
    lvc <- log(109.5); label("Central volume of distribution at 70 kg (L)") # Table 2 'theta VD central' = 109.5 (RSE 10.4%)
    lq <- log(61.37); label("Intercompartmental clearance at 70 kg (L/h)") # Table 2 'theta Q' = 61.37 (RSE 18.3%)
    lvp <- log(130.6); label("Peripheral volume of distribution at 70 kg (L)") # Table 2 'theta VD peripheral' = 130.6 (RSE 9.9%)
    lfdepot <- fixed(log(1)); label("Relative bioavailability of the oral solution (fraction)") # Methods 2.4: no bioavailability could be estimated from oral-only data; only IIV on F plasma
    e_wt_cl <- fixed(0.75); label("Allometric exponent of weight on CL and Q (unitless)") # Methods 2.4 and Table 2 footnote: (WGT/70)^0.75
    e_wt_vc <- fixed(1); label("Allometric exponent of weight on Vc and Vp (unitless)") # Methods 2.4 and Table 2 footnote: (WGT/70)^1

    # Saliva kinetics (Kruizinga 2022 Table 2 and Equations 1-3)
    lfcontam_saliva <- log(0.033e-3); label("Fraction of the dose deposited in the saliva contamination compartment (fraction)") # Table 2 'theta F saliva' = 0.033 per mille (RSE 14.5%)
    lkel_saliva <- log(1.95); label("First-order elimination rate constant of saliva contamination (1/h)") # Table 2 'theta Kel saliva' = 1.95 (RSE 6.0%); Equation 1
    lvsaliva <- fixed(log(0.001)); label("Volume of the saliva contamination compartment (L)") # Results 3.1: volume fixed to 1 mL; Equation 3 divisor 0.001 L
    lfsaliva_max <- log(0.195); label("Maximum saliva:plasma concentration ratio (unitless)") # Table 2 'theta RatioMAX' = 0.195 (RSE 16.0%); Equation 2
    lkm_fsaliva <- log(2.581); label("Plasma concentration at half-maximal saliva:plasma ratio (ug/L)") # Table 2 'theta RatioKM' = 2.581 ug/L (RSE 28.5%); Equation 2

    # IIV (Kruizinga 2022 Table 2; variances, shrinkage in the paper's brackets)
    etalka ~ 0.16 # Table 2 'omega2 ka - slow group' = 0.16 (shrinkage 2.24%)
    etalq ~ 0.25 # Table 2 'omega2 Q' = 0.25 (shrinkage 0.33%)
    etalfdepot ~ 0.026 # Table 2 'omega2 F plasma' = 0.026 (shrinkage 3.98%)
    etalfcontam_saliva ~ 0.28 # Table 2 'omega2 F saliva' = 0.28 (shrinkage 5.30%)
    etalkel_saliva ~ 0.056 # Table 2 'omega2 Kel saliva' = 0.056 (shrinkage 17.2%)

    # Residual error (Kruizinga 2022 Table 2 reports variances; SD = sqrt(variance))
    propSd <- sqrt(0.0058); label("Proportional residual error, plasma (fraction)") # Table 2 'sigma2 proportional plasma' = 0.0058 -> SD 0.0762
    propSd_Csaliva <- sqrt(0.057); label("Proportional residual error, saliva (fraction)") # Table 2 'sigma2 proportional saliva' = 0.057 -> SD 0.239
  })
  model({
    # Absorption-rate mixture: IIV applies to the slow class only
    ka_slow <- exp(lka + etalka)
    ka_fastabs <- exp(lka_fastabs)
    ka <- ka_slow * (1 - MIX_FAST_ABS) + ka_fastabs * MIX_FAST_ABS

    cl <- exp(lcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    q <- exp(lq + etalq) * (WT / 70)^e_wt_cl
    vp <- exp(lvp) * (WT / 70)^e_wt_vc
    fdepot <- exp(lfdepot + etalfdepot)

    fcontam_saliva <- exp(lfcontam_saliva + etalfcontam_saliva)
    kel_saliva <- exp(lkel_saliva + etalkel_saliva)
    vsaliva <- exp(lvsaliva)
    fsaliva_max <- exp(lfsaliva_max)
    km_fsaliva <- exp(lkm_fsaliva)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    # Equation 1: dContamination/dt = -kel * Contamination
    d/dt(saliva) <- -kel_saliva * saliva

    f(depot) <- fdepot
    f(saliva) <- fcontam_saliva

    # mg/L -> ug/L
    Cc <- 1000 * central / vc
    # Equation 2: saliva:plasma ratio = RatioMAX * Cplasma / (Cplasma + RatioKM)
    fsaliva <- fsaliva_max * Cc / (Cc + km_fsaliva)
    # Equation 3: Csaliva = Contamination / 0.001 + Cplasma * ratio
    Csaliva <- 1000 * saliva / vsaliva + Cc * fsaliva

    Cc ~ prop(propSd)
    Csaliva ~ prop(propSd_Csaliva)
  })
}
