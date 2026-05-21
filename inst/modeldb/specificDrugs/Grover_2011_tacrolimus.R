Grover_2011_tacrolimus <- function() {
  description <- "Two-compartment population PK model for oral tacrolimus in adult Native American kidney transplant recipients (Grover 2011), with first-order absorption after a lag time, no covariate effects (the Native American cohort showed no association of age, sex, weight, BMI, or post-transplant duration with PK parameters), and a placeholder proportional residual error model (residual error was not reported in the short communication)."
  reference <- "Grover A, Frassetto LA, Benet LZ, Chakkera HA. Pharmacokinetic Differences Corroborate Observed Low Tacrolimus Dosage in Native American Renal Transplant Patients. Drug Metab Dispos. 2011 Nov;39(11):2017-2019. doi:10.1124/dmd.111.041350."
  vignette <- "Grover_2011_tacrolimus"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 1L,
    n_observations = NA_integer_,
    age_range      = "Adults; mean (SD) 52 (13) years",
    age_median     = "52 years (mean)",
    weight_range   = "Mean (SD) 83.6 (19.3) kg",
    weight_median  = "83.6 kg (mean)",
    sex_female_pct = 37,
    race_ethnicity = "100% Native American (Navajo 58%, Hopi 21%, Other tribal affiliation 21%); all subjects had both parents and both sets of grandparents of American Indian tribal-group descent (Table 1).",
    disease_state  = "Adult Native American renal transplant recipients on stable oral tacrolimus dosing for at least one month post-transplant. 67% had diabetes mellitus before transplant; 25% had a history of acute rejection; 13% had new-onset diabetes mellitus post-transplant; mean (SD) duration post-transplant 30 (23) months.",
    dose_range     = "Twice-daily oral tacrolimus capsules; mean (SD) twice-daily dose 1.27 (0.644) mg [0.016 (0.010) mg/kg]; mean (SD) total daily dose 2.54 (1.22) mg [0.033 (0.021) mg/kg/day, excluding one subject on therapy < 4 months with a higher target trough]; titrated to hospital target trough levels of 10-12 ng/mL within the first month post-transplant, 8-10 ng/mL between months 1 and 4, and 5-8 ng/mL after month 4; mean (SD) measured trough 6.53 (2.43) ng/mL.",
    regions        = "United States (Mayo Clinic Arizona, Phoenix, AZ).",
    bmi_mean_sd    = "29.9 (5) kg/m^2",
    tribal_affiliation = "Navajo 58%, Hopi 21%, Other 21%",
    sampling_design = "Single 12-hour steady-state pharmacokinetic profile per subject: blood drawn pre-dose and at 0.5, 1, 2, 4, 6, 8, and 12 hours after the morning capsule dose, following an overnight fast. EDTA whole-blood samples analysed by the Architect tacrolimus immunoassay (Bazin 2010) at Mayo Clinic Arizona.",
    co_medication  = "No concomitant medications, supplements, or foods known to interact with tacrolimus (antifungals, antiepileptics, macrolides, St. John's wort, grapefruit) at the time of the PK profile.",
    notes          = "Population estimates and IIV (%CV) come from Grover 2011 Table 2 'NONMEM Parameter Estimates' row; the model was fit by NONMEM 7.1 using an empirical Bayesian approach on the 24-subject single-dose steady-state PK profiles. The 'Mean estimate' and 'S.D.' rows in Table 2 are summaries of individual Bayesian posterior modes and are not what populates `ini()` -- the population estimate row is. IIV in Vss/F was not estimable from the 24-subject dataset (Table 2: 'n.e.' for Vss/F); all subjects therefore share the same population Vss/F = 462 L, and Vp/F = Vss/F - V/F = 462 - 73.3 = 388.7 L inherits no IIV in the source model."
  )

  ini({
    # Structural PK parameters from Grover 2011 Table 2 (NONMEM Parameter
    # Estimates row, p. 2018). All clearances and volumes are apparent values
    # (CL/F, V/F, Q/F, Vss/F) because the source data are oral. Bioavailability
    # is folded into the apparent values and is not separately parameterised.
    lka   <- log(1.38)  ; label("First-order absorption rate constant ka (1/h)")             # Table 2 NONMEM Parameter Estimates: ka = 1.38 1/h
    ltlag <- log(0.573) ; label("Absorption lag time tlag (h)")                                # Table 2 NONMEM Parameter Estimates: tlag = 0.573 h
    lcl   <- log(10.1)  ; label("Apparent oral clearance CL/F (L/h)")                          # Table 2 NONMEM Parameter Estimates: CL/F = 10.1 L/h
    lvc   <- log(73.3)  ; label("Apparent central volume V/F (L)")                             # Table 2 NONMEM Parameter Estimates: V/F = 73.3 L
    lq    <- log(27.1)  ; label("Apparent inter-compartmental clearance Q/F (L/h)")            # Table 2 NONMEM Parameter Estimates: Q/F = 27.1 L/h
    # Apparent peripheral volume Vp/F is fixed by the source's Vss/F
    # parameterisation: the paper reports Vss/F = 462 L and notes that all
    # subjects share this population value (IIV n.e.), so the peripheral
    # volume Vp/F = Vss/F - V/F = 462 - 73.3 = 388.7 L is held at this
    # derived constant with no IIV. This matches the secondary calculated
    # V2/F mean of 391 L in Table 2.
    lvp   <- fixed(log(388.7)); label("Apparent peripheral volume Vp/F = Vss/F - V/F = 462 - 73.3 (L; fixed)") # Table 2 NONMEM Parameter Estimates: Vss/F = 462 L, V/F = 73.3 L; Vp/F derived and fixed (Vss IIV not estimable)

    # Inter-individual variability. Grover 2011 Table 2 reports IIV as %CV;
    # convert to log-scale variance via omega^2 = log(CV^2 + 1):
    #   CL/F  CV 43.5% -> log(0.435^2 + 1) = 0.17326
    #   V/F   CV 36.9% -> log(0.369^2 + 1) = 0.12772
    #   Q/F   CV 50.4% -> log(0.504^2 + 1) = 0.22641
    #   ka    CV 46.0% -> log(0.460^2 + 1) = 0.19196
    #   tlag  CV 13.3% -> log(0.133^2 + 1) = 0.01754
    # Vss/F IIV was not estimable in NONMEM with the 24-subject dataset
    # (Table 2: 'n.e.'); no eta is carried on lvp.
    etalka   ~ 0.19196   # Table 2 NONMEM Parameter Estimates: IIV ka = 46.0% CV
    etaltlag ~ 0.01754   # Table 2 NONMEM Parameter Estimates: IIV tlag = 13.3% CV
    etalcl   ~ 0.17326   # Table 2 NONMEM Parameter Estimates: IIV CL/F = 43.5% CV
    etalvc   ~ 0.12772   # Table 2 NONMEM Parameter Estimates: IIV V/F = 36.9% CV
    etalq    ~ 0.22641   # Table 2 NONMEM Parameter Estimates: IIV Q/F = 50.4% CV

    # Residual error: NOT reported by Grover 2011. The Materials and Methods
    # section states "Pharmacokinetic parameters were estimated using NONMEM
    # (version 7.1) ... using an empirical Bayesian approach" but neither the
    # main text, Table 2, nor any supplement reports the residual-error model
    # form or magnitude. A small fixed proportional placeholder (1%) is
    # supplied so the observation has a valid residual model and typical-value
    # simulations are effectively deterministic. Users who need realistic
    # stochastic prediction intervals should substitute a tacrolimus-class
    # typical value (e.g. 15-30% CV proportional) and document the choice.
    # Precedent: Jelliffe_2014_digoxin.R and Taylor_2020_methotrexate.R apply
    # the same placeholder pattern.
    propSd <- fixed(0.01); label("Proportional residual error (placeholder; not reported in source)")  # not in Grover 2011 -- see vignette Assumptions and deviations
  })

  model({
    # Individual PK parameters. No covariate effects: Grover 2011 Results
    # state "Within all 24 patients, no other demographic characteristics
    # (age, gender, weight, BMI, or total days on therapy) were associated
    # with clearance or other pharmacokinetic parameters."
    ka   <- exp(lka   + etalka)
    tlag <- exp(ltlag + etaltlag)
    cl   <- exp(lcl   + etalcl)
    vc   <- exp(lvc   + etalvc)
    q    <- exp(lq    + etalq)
    vp   <- exp(lvp)

    # Two-compartment oral PK with first-order absorption + lag time. Dose
    # lands in `depot`; bioavailability is folded into the apparent CL/V/Q
    # estimates so f(depot) is not assigned (F = 1 structurally).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Tacrolimus is reported in whole blood in ng/mL (Architect immunoassay).
    # Dose units mg, volumes L -> central / vc returns mg/L = ug/mL; multiply
    # by 1000 to express the prediction in ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
