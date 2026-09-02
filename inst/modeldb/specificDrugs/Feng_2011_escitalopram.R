Feng_2011_escitalopram <- function() {
  description <- "Two-compartment population PK model with first-order absorption for escitalopram, used as the long half-life drug scenario in a trial-simulation study of Cpred/Cobs and Cipred/Cobs concentration ratios as adherence metrics (Feng 2011)"
  reference <- paste(
    "Feng Y, Gastonguay MR, Pollock BG, Frank E, Kepple GH, Bies RR.",
    "Performance of Cpred/Cobs concentration ratios as a metric reflecting",
    "adherence to antidepressant drug therapy.",
    "Neuropsychiatr Dis Treat. 2011;7:117-125. doi:10.2147/NDT.S15921.",
    "The parameter set in Feng 2011 Table 2 was adapted by those authors from",
    "the published escitalopram literature (their refs 15-16:",
    "Sogaard B, Mengel H, Rao N, Larsen F. J Clin Pharmacol. 2005;45:1400-1406;",
    "Gutierrez MM, Rosenberg J, Abramowitz W. Clin Ther. 2003;25:1200-1210);",
    "Feng 2011 is the transcription source for every value below.",
    sep = " "
  )
  vignette <- "Feng_2011_escitalopram"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Feng 2011: the Table 2 abbreviations
  # footnote defines V2 as the "volume of distribution of central compartment"
  # and V3 as the "volume of distribution of peripheral compartment", and the
  # Methods describe "a simulated (observed) plasma drug concentration (Cobs)"
  # after "oral administration".
  compartmentData <- list(
    depot       = list(analyte = "escitalopram", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "escitalopram", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "escitalopram", units = "mg", specimen = "plasma", verified = TRUE)
  )

  # Feng 2011 Table 2 reports structural and variance parameters only; the
  # simulation model carries no covariate effects. The MEMS study cohort that
  # supplied the dosing histories is described in the Results ("Subjects and
  # MEMS data") but no demographic covariate enters the pharmacokinetic model.
  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 65,
    n_studies      = 1,
    disease_state  = "Chronic psychiatric disorders; unipolar depression cohort of the NIMH-sponsored 'Depression: the search for treatment relevant phenotypes' trial (Frank 2008, Feng 2011 ref 12)",
    dose_range     = "10 mg oral once daily (Feng 2011 Table 3, 'Dose (mg): Long half-life drug (10)')",
    regions        = "United States (University of Pittsburgh) and Italy (Pisa); Pittsburgh-Pisa phenotypes study",
    notes          = paste(
      "Feng 2011 is a trial-simulation study. The 65 subjects supplied the",
      "'true' dosage histories via Medication Event Monitoring System (MEMS)",
      "caps -- 863 clinic-visit records over the first 6 months -- and set the",
      "simulation sample size (Table 3, 'Sample size for each simulation",
      "replicate: 65'; 100 replicates). They did NOT supply concentration data:",
      "the pharmacokinetic parameters in Table 2 were adapted from the",
      "published escitalopram literature and used as fixed simulation inputs,",
      "so every ini() value below is wrapped in fixed(). Age, weight, sex and",
      "race distributions for the MEMS cohort are not reported in Feng 2011.",
      sep = " "
    ),
    parameter_provenance = paste(
      "Structural and variance parameters are literature-adapted, not fitted in",
      "Feng 2011: 'An established escitalopram model was utilized to represent",
      "the long half-life drug. A two-compartment model with additive and",
      "proportional residual error models was adapted from the literature",
      "reports for escitalopram as a long half-life drug (Table 2).' No",
      "standard errors, RSEs or confidence intervals are reported for any",
      "value.",
      sep = " "
    ),
    sampling_design = paste(
      "One plasma sample per simulated clinic visit; 18 clinic-visit records per",
      "subject on average (range 2-43). Sampling time drawn uniformly between",
      "8 am and 6 pm (clinic hours); nominal dose time 9 pm with SD 1 h",
      "(Feng 2011 Table 3).",
      sep = " "
    )
  )

  ini({
    # Structural parameters. Every value is a FIXED simulation input in
    # Feng 2011 -- adapted from the prior escitalopram literature rather than
    # estimated here -- and is reported without any uncertainty, so all are
    # wrapped in fixed(). Parameterisation is NONMEM ADVAN4 TRANS4
    # (KA, CL, V2, Q, V3) per the Methods: "the NONMEM program
    # (two-compartment, ADVAN4 TRANS4)". Data are oral-only, so CL and the
    # volumes are apparent (CL/F, V/F) although Feng 2011 labels them CL, V2, V3.
    lka <- fixed(log(0.16));  label("Absorption rate constant (Ka, 1/h)")                          # Feng 2011 Table 2, row Ka (/h) = 0.16
    lcl <- fixed(log(24.5));  label("Apparent oral clearance (CL/F, L/h)")                         # Feng 2011 Table 2, row CL (L/h) = 24.5
    lvc <- fixed(log(417));   label("Apparent central volume of distribution (V2/F, L)")           # Feng 2011 Table 2, row V2 (L) = 417
    lq  <- fixed(log(35.7));  label("Apparent inter-compartmental clearance (Q/F, L/h)")           # Feng 2011 Table 2, row Q (L/h) = 35.7
    lvp <- fixed(log(541));   label("Apparent peripheral volume of distribution (V3/F, L)")        # Feng 2011 Table 2, row V3 (L) = 541

    # Inter-individual variability. Feng 2011 states the IIV terms are
    # "log-normally distributed, with a mean of zero and variance of omega^2"
    # and tabulates omega as a percentage, so the variance used here is
    # (omega%/100)^2 -- the paper's own arithmetic. Table 2 reports no IIV on
    # Ka, so none is encoded. See the vignette Errata for the alternative
    # exact-CV reading omega^2 = log(CV^2 + 1), which differs by <11% in
    # variance and is NOT what this paper's text defines.
    #
    # Source trace for the four omegas below. These comments sit ABOVE their
    # ini() lines rather than trailing them: rxode2 promotes a trailing comment
    # on an ini() entry that has no label() into a label(), and a quoted or
    # multi-token trailing comment then fails to re-parse.
    #   etalcl -- Feng 2011 Table 2, row omega CL% = 50 -> 0.50^2 = 0.2500
    #   etalvc -- Feng 2011 Table 2, row omega v2% = 35 -> 0.35^2 = 0.1225
    #   etalq  -- Feng 2011 Table 2, row omega Q%  = 30 -> 0.30^2 = 0.0900
    #   etalvp -- Feng 2011 Table 2, row omega v3% = 30 -> 0.30^2 = 0.0900
    etalcl ~ fixed(0.2500)
    etalvc ~ fixed(0.1225)
    etalq  ~ fixed(0.0900)
    etalvp ~ fixed(0.0900)

    # Residual unexplained variability. The Methods describe "additive and
    # proportional residual error models", but Table 2 tabulates only a single
    # residual term, "sigma 1% 30", which the Table 2 footnote defines as a
    # "coefficient of variation of residual error" -- i.e. the proportional
    # component. No value is given anywhere in the paper for the additive
    # component, so it is carried at fixed(0) to preserve the stated combined
    # structure without inventing a number. See the vignette Errata.
    propSd <- fixed(0.30);  label("Proportional residual error (fraction)")   # Feng 2011 Table 2, row sigma 1% = 30
    addSd  <- fixed(0);     label("Additive residual error (ug/mL)")          # Feng 2011 Methods states a combined additive + proportional error model; no additive value is reported in Table 2 or anywhere else in the paper.
  })

  model({
    # Individual parameters. Exponential (log-normal) IIV per the Methods:
    # "CL = TVCL x exp(eta CL)".
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # Micro-constants (ADVAN4 TRANS4 -> ADVAN4 TRANS3 equivalents).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Observation: dose in mg, volume in L -> mg/L == ug/mL. At the 10 mg
    # once-daily regimen of Feng 2011 Table 3 this is ~0.02-0.05 ug/mL
    # (20-50 ng/mL); the vignette reports ng/mL for readability.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
