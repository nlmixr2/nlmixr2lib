Richardson_2025_ribociclib <- function() {
  description <- paste(
    "Automated model-search (pyDarwin) two-compartment population PK model",
    "for oral ribociclib, fitted to a synthetic 96-subject data set that was",
    "simulated from a previously published ribociclib popPK model at 100,",
    "200, 400 and 600 mg (Richardson 2025). This is the top-ranked structure",
    "returned by the paper's automated Bayesian-optimisation model search,",
    "and the paper reports that it recovered the data-generating structure",
    "exactly. Absorption is zero-order directly into the central compartment",
    "over a modelled duration D with an absorption lag time; the search space",
    "carries a depot compartment but fixes its transfer rate constant KA to",
    "zero for the zero-order absorption option, so the depot is inactive and",
    "is omitted here. Disposition is two-compartment linear with apparent",
    "clearance and volumes (bioavailability fixed to 1, so all disposition",
    "parameters are apparent CL/F and V/F). Dose enters as a power-model",
    "covariate on CL, Vp and Q, normalised to a reference dose. IIV is",
    "diagonal log-normal on CL, Vc, Vp and Q, with the lag time and the",
    "zero-order duration carrying IIV fixed at 10% CV. Residual error is",
    "combined additive plus proportional. Note that this is a machine-search",
    "structure fitted to synthetic data, not the manually developed expert",
    "model; the expert ribociclib model is Lu et al. (2021).",
    sep = " "
  )
  reference <- paste(
    "Richardson S, Irurzun Arana I, Nowojewski A, Zhou D, Leander J,",
    "Tang W, Dearden R, Gibbs M.",
    "A machine learning approach to population pharmacokinetic modelling",
    "automation. Commun Med. 2025;5:325. doi:10.1038/s43856-025-01054-8.",
    "Parameter estimates from Supplementary Table 2; omega values and dose",
    "effects from the authors' published NONMEM output FinalResultFile.lst",
    "in the Code Availability repository",
    "https://github.com/samjrrr/autopk_synthetic_example.",
    "The synthetic data were generated from the expert ribociclib model of",
    "Lu Y, Yang S, Ho Y-Y, Ji Y. J Clin Pharmacol. 2021;61:1054-1068.",
    sep = " "
  )
  vignette <- "Richardson_2025_automated_popPK"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    central     = list(analyte = "ribociclib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ribociclib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    DOSE = list(
      description        = "Administered dose amount of the current dosing record",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters as a power model on CL, Vp and Q, normalised to a reference",
        "dose: P = TVP * (DOSE / 600)^e_dose_P. The authors' NONMEM data set",
        "carried DOSE and a separate DOSEREF column and the model term was",
        "MU(P) = LOG(TVP) + LOG(DOSE/DOSEREF) * THETA(cov_P) (template.txt and",
        "FinalControlFile.mod in the Code Availability repository). The numeric",
        "value of DOSEREF is not reported for any data set in the paper or its",
        "supplement; 600 mg (the highest of the four simulated dose levels and",
        "the therapeutic ribociclib dose) is used here as the reference so that",
        "the tabulated typical values apply at 600 mg. See the validation",
        "vignette Assumptions and deviations section - a different reference",
        "dose D rescales each affected typical value by (600 / D)^e_dose_P and",
        "leaves the shape of the dose dependence unchanged.",
        sep = " "
      ),
      source_name        = "DOSE"
    )
  )

  population <- list(
    species     = "human (synthetic data simulated from a published human popPK model)",
    n_subjects  = 96,
    n_studies   = 1,
    dose_range  = "100, 200, 400 and 600 mg oral",
    notes       = paste(
      "Synthetic data set, not a clinical cohort: 96 individuals simulated",
      "from the published ribociclib popPK model of Lu et al. (2021) with the",
      "covariate effects excluded, using sampling that mimicked the Phase 1",
      "study (pre-dose and 0.5, 1, 2, 4, 6, 8 and 24 h after dosing on day 1",
      "and day 21). Richardson 2025 Methods 'Data' and Table 1; 1,344",
      "observation records over 96 individuals per Supplementary Table 7 and",
      "the published NONMEM output. No demographic covariates were simulated,",
      "so no baseline demographics are reported.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- apparent (F fixed to 1) typical values at the
    # 600 mg reference dose. Supplementary Table 2 reports the
    # back-transformed values; the corresponding log-scale THETAs are in the
    # authors' published FinalResultFile.lst (THETA1, THETA2, THETA7, THETA8).
    lcl  <- log(23.40)  ; label("Apparent clearance CL/F at the reference dose (L/h)")                     # Suppl Table 2 THETA1 = 23.40 L/h
    lvc  <- log(244.38) ; label("Apparent central volume V2/F (L)")                                        # Suppl Table 2 THETA2 = 244.38 L
    lvp  <- log(714.62) ; label("Apparent peripheral volume V3/F at the reference dose (L)")               # Suppl Table 2 THETA7 = 714.62 L
    lq   <- log(68.67)  ; label("Apparent inter-compartmental clearance Q/F at the reference dose (L/h)")  # Suppl Table 2 THETA8 = 68.67 L/h
    ltlag <- log(0.36)  ; label("Absorption lag time ALAG (h)")                                           # Suppl Table 2 THETA11 = 0.36 h
    ld1  <- log(3.47)   ; label("Duration of zero-order input into the central compartment D (h)")         # Suppl Table 2 THETA12 = 3.47 h

    # Bioavailability was fixed to 1 in the model space (template.txt: F1 = 1),
    # so every disposition parameter above is apparent.
    lfcentral <- fixed(log(1)) ; label("Bioavailability F (unitless)")                        # template.txt F1 = 1   ; INPUT1

    # Dose effects: power model on (DOSE / 600). Taken from the published
    # NONMEM output, which carries one more significant figure than the
    # rounded values (-0.49, -0.61, -0.75) in Supplementary Table 2.
    e_dose_cl <- -0.493 ; label("Power exponent on (DOSE / 600) for CL (unitless)")  # FinalResultFile.lst THETA3 = -4.93E-01; Suppl Table 2 -0.49
    e_dose_vp <- -0.606 ; label("Power exponent on (DOSE / 600) for Vp (unitless)")  # FinalResultFile.lst THETA9 = -6.06E-01; Suppl Table 2 -0.61
    e_dose_q  <- -0.747 ; label("Power exponent on (DOSE / 600) for Q (unitless)")   # FinalResultFile.lst THETA10 = -7.47E-01; Suppl Table 2 -0.75

    # IIV -- diagonal log-normal variances, taken directly from the OMEGA
    # matrix of the published NONMEM output. Supplementary Table 2 reports
    # these only as CV%, computed as sqrt(exp(omega^2) - 1) * 100; the values
    # below round-trip to the tabulated 52%, 121%, 53%, 61%, 10% and 10%.
    etalcl   ~ 0.243        # FinalResultFile.lst OMEGA(1,1) = 2.43E-01 -> 52% CV (Suppl Table 2)
    etalvc   ~ 0.903        # FinalResultFile.lst OMEGA(2,2) = 9.03E-01 -> 121% CV (Suppl Table 2)
    etalvp   ~ 0.244        # FinalResultFile.lst OMEGA(4,4) = 2.44E-01 -> 53% CV (Suppl Table 2)
    etalq    ~ 0.315        # FinalResultFile.lst OMEGA(5,5) = 3.15E-01 -> 61% CV (Suppl Table 2)
    etaltlag ~ fixed(0.01)  # FinalResultFile.lst OMEGA(6,6) = 1.00E-02 FIX -> 10% CV (Suppl Table 2 10% (fix))
    etald1   ~ fixed(0.01)  # FinalResultFile.lst OMEGA(7,7) = 1.00E-02 FIX -> 10% CV (Suppl Table 2 10% (fix))

    # Residual error -- combined additive plus proportional on the linear
    # concentration scale. $SIGMA was 1 FIX and the error was carried in
    # THETA space: W = SQRT(add^2 + (prop * IPRED)^2); Y = IPRED + W * ERR(1).
    addSd  <- 3.42  ; label("Additive residual error (ng/mL)")            # Suppl Table 2 THETA5 = 3.42
    propSd <- 0.217 ; label("Proportional residual error (fraction)")     # FinalResultFile.lst THETA6 = 2.17E-01; Suppl Table 2 0.22
  })

  model({
    # Individual parameters. The dose effect is a power model on the dose
    # normalised to the 600 mg reference dose, equivalent to the authors'
    # MU(P) = LOG(TVP) + LOG(DOSE/DOSEREF) * THETA(cov_P).
    cl <- exp(lcl + etalcl) * (DOSE / 600)^e_dose_cl
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp) * (DOSE / 600)^e_dose_vp
    q  <- exp(lq  + etalq)  * (DOSE / 600)^e_dose_q

    tlag <- exp(ltlag + etaltlag)
    d1   <- exp(ld1   + etald1)

    # Micro-constants, matching the ADVAN5 parameterisation of the source
    # control stream (k20 = CL/V2, k23 = Q/V2, k32 = Q/V3).
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central)     <- k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Zero-order input straight into the central compartment after the lag
    # time. Dose records use rate = -2 so rxode2 takes the duration from
    # dur(central), matching RATE = -2 in the authors' NONMEM data set.
    f(central)   <- exp(lfcentral)
    dur(central) <- d1
    alag(central) <- tlag

    # central is in mg and vc in L, so central / vc is mg/L; the factor of
    # 1000 puts the observation on the ng/mL scale that the additive residual
    # error (3.42) is expressed in.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
