Kaikousidis_2023_donepezil_oral_meandata_mouse <- function() {
  description <- "Preclinical (mouse). Three-compartment blood + brain PK model for donepezil after an oral solution (first-order absorption from a gut depot, blood-peripheral exchange, unidirectional blood-to-brain transfer and first-order brain elimination), fitted by naive averaging to the time-point MEAN serum and brain concentrations (Kaikousidis 2023 Table 2, 'Per Os - Mean Data')"
  reference <- "Kaikousidis C, Papakyriakopoulou P, Dokoumetzidis A, Valsami G. Donepezil Brain and Blood Pharmacokinetic Modeling after Nasal Film and Oral Solution Administration in Mice. Pharmaceutics. 2023;15(5):1409. doi:10.3390/pharmaceutics15051409"
  vignette <- "Kaikousidis_2023_donepezil_mouse"
  units <- list(time = "min", dosing = "mg", concentration = "ug/mL")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "donepezil", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "donepezil", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "donepezil", units = "mg", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "donepezil", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "mouse (C57BL/6J)",
    n_subjects = NA_integer_,
    n_studies = 1,
    age_range = "8 weeks",
    weight_range = "about 25 g (10 mg/kg for a 0.25 mg oral dose; 4 mg/kg for a 0.1 mg nasal film)",
    sex_female_pct = NA_real_,
    disease_state = "healthy",
    dose_range = "single oral dose of 0.2 mL of a 1.25 mg/mL donepezil solution by gastric gavage (0.25 mg, 10 mg/kg)",
    notes = "Destructive (one sample per animal) design: serum and perfused brain sampled at 5, 10, 15, 30, 45, 60, 120, 240 and 360 min after the oral dose, n >= 5 animals per time point (Section 2.1). This file is the fit to the per-time-point MEAN concentrations (Section 2.2 'mean data'); its raw-data sibling is Kaikousidis_2023_donepezil_oral_rawdata_mouse. Values below the LOQ (serum 0.043 ug/mL, brain 0.038 ug/mL) were set to LOQ/2. Sex and the exact number of animals are not reported. No between-animal variability is estimated; the model is a typical-value model with residual error only."
  )

  ini({
    # Kaikousidis 2023 Table 2, 'Per Os - Mean Data' block. The paper numbers
    # its states 1 = depot, 2 = central/blood, 3 = brain, 4 = peripheral
    # (Figure 2), so k12 is the absorption rate constant, k24/k42 the
    # blood-peripheral pair (canonical k12/k21) and k23 the blood-to-brain
    # transfer.
    lka <- log(0.059); label("First-order absorption rate constant, paper k12 (1/min)") # Table 2 Per Os Mean 'k12 (1/min)' 0.059
    lkel <- log(0.012); label("Elimination rate constant from blood, paper k20 (1/min)") # Table 2 Per Os Mean 'k20 (1/min)' 0.012
    lk_central_brain <- log(0.014); label("Blood-to-brain transfer rate constant, paper k23 (1/min)") # Table 2 Per Os Mean 'k23 (1/min)' 0.014
    lk12 <- log(0.034); label("Blood-to-peripheral rate constant, paper k24 (1/min)") # Table 2 Per Os Mean 'k24 (1/min)' 0.034
    lk21 <- log(0.011); label("Peripheral-to-blood rate constant, paper k42 (1/min)") # Table 2 Per Os Mean 'k42 (1/min)' 0.011
    lkeff_brain <- log(0.325); label("Lumped brain elimination rate constant, paper k30 (1/min)") # Table 2 Per Os Mean 'k30 (1/min)' 0.325
    lvc <- fixed(log(0.279)); label("Apparent blood volume of distribution Vblood (L)") # Table 1 'Vblood (L)' 0.279; Section 3.1 fixes it in the combined Per Os model
    lv_brain <- log(0.018); label("Apparent brain volume of distribution Vbrain (L)") # Table 2 Per Os Mean 'Vbrain (L)' 0.018

    # Residual error, Table 2 'Per Os - Mean Data' b1 (blood) and b2 (brain),
    # both proportional (Section 3.1). Read as NONMEM variances; see the
    # vignette for the evidence.
    propSd <- sqrt(0.012); label("Proportional residual error on serum concentration (fraction)") # Table 2 Per Os Mean 'b1' 0.012 (variance)
    propSd_Cbrain <- sqrt(0.062); label("Proportional residual error on brain concentration (fraction)") # Table 2 Per Os Mean 'b2' 0.062 (variance)
  })
  model({
    ka <- exp(lka)
    kel <- exp(lkel)
    k_central_brain <- exp(lk_central_brain)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    keff_brain <- exp(lkeff_brain)
    vc <- exp(lvc)
    v_brain <- exp(lv_brain)

    # Figure 2: blood-to-brain transfer is unidirectional and depletes blood;
    # keff_brain lumps every loss from brain (return to blood, lymphatic
    # clearance, metabolism), as in a parent-metabolite model (Section 3.1).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + k12 + k_central_brain) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(brain) <- k_central_brain * central - keff_brain * brain

    Cc <- central / vc
    Cbrain <- brain / v_brain
    Cc ~ prop(propSd)
    Cbrain ~ prop(propSd_Cbrain)
  })
}
