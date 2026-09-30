Bai_2021_mkAsodnLiposome_human_pbpk <- function() {
  description <- paste(
    "PBPK-derived reduced one-compartment intravenous model for MK-ASODN",
    "nanoliposomes (a nanoliposome-encapsulated antisense oligonucleotide",
    "against midkine mRNA, developed for hepatocellular carcinoma) in humans,",
    "predicted from monkey data ahead of first-in-human dosing. The source",
    "paper extrapolated its GastroPlus 8.0 monkey PBPK model to a human",
    "(Chinese male) GastroPlus PBPK model, substituting human plasma protein",
    "binding, predicting human Vss with the Rodgers-Single (Lukacova) method",
    "and scaling clearance from monkey by single-species allometry with a",
    "fixed exponent of 0.8. The 14-tissue whole-body structure, organ volumes,",
    "blood flows and Kp values are platform internals that the paper does not",
    "print and are not reproduced here. This file encodes the reduced",
    "disposition card the paper does print: CL = 4 L/h and Vss = 7.89 L. The",
    "reduction reproduces the paper's predicted AUC0-inf after 90 mg to 0.1%",
    "and the terminal half-life of the simulated curve in Figure 3 to within",
    "3%. It does not reproduce the tabulated predicted Cmax of 49.98 ug/mL,",
    "which is the platform's instantaneous-bolus spike into a small plasma",
    "volume and lasts under about 10 minutes; the reduction starts at dose /",
    "Vss (11.4 ug/mL), which is where the simulated curve settles after the",
    "spike. This is a prediction, not a fit to human data, and the paper",
    "reports no variability, so no etas are declared and both residual error",
    "terms are fixed(0).",
    sep = " "
  )
  reference <- paste(
    "Bai H, Cheng Y, Che J. Pharmacokinetics and Disposition of Heparin-Binding",
    "Growth Factor Midkine Antisense Oligonucleotide Nanoliposomes in",
    "Experimental Animal Species and Prediction of Human Pharmacokinetics Using",
    "a Physiologically Based Pharmacokinetic Model.",
    "Front Pharmacol. 2021;12:769538. doi:10.3389/fphar.2021.769538.",
    "PMCID PMC8595129.",
    "The human PBPK method, allometric equation and 90 mg dose derivation are",
    "Methods 2.7; human CL and Vss are Results 3.6; predicted human Cmax and",
    "AUC0-inf are Table 4; the simulated human plasma profile is Figure 3.",
    sep = " "
  )
  vignette <- "Bai_2021_mkAsodnLiposome_pbpk"
  # The paper reports concentrations in ug/mL and AUC in ug.h/mL; mg/L is the
  # numerically identical form.
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "MK-ASODN", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 0,
    age_range = NA_character_,
    weight_range = NA_character_,
    sex_female_pct = 0,
    race_ethnicity = "Chinese (GastroPlus virtual subject)",
    disease_state = paste(
      "Predicted, not observed: a single GastroPlus virtual Chinese male",
      "(Methods 2.7). No human data were available.",
      sep = " "
    ),
    dose_range = paste(
      "90 mg single intravenous injection, derived from the 46 mg/kg monkey",
      "dose by body-surface-area conversion with a safety factor of 10",
      "(Methods 2.7).",
      sep = " "
    ),
    regions = "China",
    notes = paste(
      "Prediction for first-in-human dose selection. Clearance scaled from a",
      "monkey clearance by CL_human = CL_monkey * (BW_human / BW_monkey)^0.8",
      "(Methods 2.7); neither body weight is printed. Human plasma protein",
      "binding was 96-99% (Table 3).",
      sep = " "
    )
  )

  ini({
    # Results 3.6: 'The human CL, extrapolated from monkeys, was 4 L h-1, and
    # the Vss was 7.89 L.' Dose / CL = 90 / 4 = 22.50 mg.h/L against the
    # Table 4 predicted AUC0-inf of 22.48 ug.h/mL.
    lcl <- fixed(log(4)) ; label("Clearance (L/h)")  # Results 3.6, human CL = 4 L/h
    lvc <- fixed(log(7.89)) ; label("Volume of distribution, Vss used as the single volume (L)")  # Results 3.6, human Vss = 7.89 L

    # Residual error. A prediction for one virtual subject with no reported
    # variability, so both terms are fixed at zero rather than invented.
    addSd <- fixed(0) ; label("Additive residual error (mg/L; not published by Bai 2021)")
    propSd <- fixed(0) ; label("Proportional residual error (fraction; not published by Bai 2021)")
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)

    kel <- cl / vc

    # Intravenous injection only (Methods 2.7); no absorption model.
    d/dt(central) <- -kel * central

    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
