Morschhauser_2018_ibritumomabTiuxetan_fraction2 <- function() {
  description <- paste(
    "Five-compartment dosimetry PK model for 111In-ibritumomab tiuxetan (the",
    "imaging surrogate of 90Y-ibritumomab tiuxetan) in adults with follicular",
    "lymphoma after the second radioimmunotherapy fraction (F2) of the FIZZ",
    "trial, built from SPECT/CT organ activity and blood samples. Blood is the",
    "central compartment, exchanging first-order with the L2-L4 vertebral",
    "marrow and the liver; the spleen takes up drug from blood and drains to",
    "the liver; elimination is from blood only. Typical values are the means of",
    "the per-patient Kinetica fits (Table 3) and the k's act on amounts. The",
    "deep compartment (remaining organs and tumor) is carried with its",
    "unpublished exchange constants fixed to zero, so the model does not",
    "reproduce the paper's printed blood MRT, half-life or Figure 2 profile.",
    sep = " "
  )
  reference <- paste(
    "Morschhauser F, Dekyndt B, Baillet C, Barthelemy C, Malek E, Fulcrand J,",
    "Bigot P, Huglo D, Decaudin B, Simon N, Odou P. A new pharmacokinetic model",
    "for 90Y-ibritumomab tiuxetan based on 3-dimensional dosimetry.",
    "Sci Rep. 2018;8:14860. doi:10.1038/s41598-018-33160-0",
    sep = " "
  )
  vignette <- "Morschhauser_2018_ibritumomabTiuxetan"
  units <- list(time = "h", dosing = "ug", concentration = "ug/mL")

  compartmentData <- list(
    central = list(
      analyte = "111In-ibritumomab tiuxetan",
      units = "ug",
      specimen = "whole blood",
      verified = TRUE
    ),
    bonemarrow = list(
      analyte = "111In-ibritumomab tiuxetan",
      units = "ug",
      specimen = "tissue",
      verified = TRUE
    ),
    liver = list(
      analyte = "111In-ibritumomab tiuxetan",
      units = "ug",
      specimen = "tissue",
      verified = TRUE
    ),
    spleen = list(
      analyte = "111In-ibritumomab tiuxetan",
      units = "ug",
      specimen = "tissue",
      verified = TRUE
    ),
    other = list(
      analyte = "111In-ibritumomab tiuxetan",
      units = "ug",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 13L,
    n_studies = 1L,
    age_range = "28-80 years",
    age_mean = "59.6 years",
    weight_range = "49-83 kg",
    weight_mean = "72.9 kg",
    sex_female_pct = 76.9,
    disease_state = paste(
      "Previously untreated CD20-positive follicular lymphoma needing",
      "treatment (GELF/BNLI criteria); FLIPI low / intermediate / high risk",
      "4 / 3 / 6 patients"
    ),
    dose_range = paste(
      "138 MBq 111In-ibritumomab tiuxetan (specific activity 138 MBq/mg, i.e.",
      "1 mg) IV with the first of two pre-dosing rituximab infusions",
      "(250 mg/m2) of fraction 2, given 8-12 weeks after fraction 1; the therapeutic 90Y-ibritumomab tiuxetan",
      "fraction was 11.1 MBq/kg (maximum 888 MBq)"
    ),
    regions = "France (Lille university hospital; dosimetry centralised in Nantes)",
    notes = paste(
      "Dosimetric subset of the FIZZ phase II trial (NCT01493479); Table 1,",
      "F2 column. Six of the 19 F1 patients did not receive F2 (haematologic",
      "toxicity or HAMA). Mean tumor volume 1850 mm3 (range 0-8859). The first",
      "fraction is a separate fit in",
      "Morschhauser_2018_ibritumomabTiuxetan_fraction1."
    )
  )

  ini({
    # Blood (central). Table 3 reports CL in mL/min and the rate constants in
    # min^-1; each value is multiplied by 60 here so the model runs in hours.
    lcl <- log(0.74 * 60); label("Clearance from blood (mL/h)") # Table 3, Blood, F2 'Clearance' = 0.74 mL/min
    lvc <- log(4064.8); label("Blood volume of distribution (mL)") # Table 3, Blood, F2 'Vd' = 4064.8 mL

    # Fig. 3 compartments: 1 blood, 2 L2-L4 vertebrae, 3 liver, 4 spleen, 5 deep
    lkin_bonemarrow <- log(0.60e-4 * 60); label("Blood to L2-L4 vertebral marrow rate constant, k12 (1/h)") # Table 3, Vertebrae L2-L4, F2 k12 = 0.60E-4 min^-1
    lkout_bonemarrow <- log(0.96e-4 * 60); label("L2-L4 vertebral marrow to blood rate constant, k21 (1/h)") # Table 3, Vertebrae L2-L4, F2 k21 = 0.96E-4 min^-1
    lkin_liver <- log(2.05e-4 * 60); label("Blood to liver rate constant, k13 (1/h)") # Table 3, Liver, F2 k13 = 2.05E-4 min^-1
    lkout_liver <- log(1.77e-4 * 60); label("Liver to blood rate constant, k31 (1/h)") # Table 3, Liver, F2 k31 = 1.77E-4 min^-1
    lkin_spleen <- log(1.17e-4 * 60); label("Blood to spleen rate constant, k14 (1/h)") # Table 3, Spleen, F2 row printed 'k13' = 1.17E-4 min^-1; Fig. 3 arrow k14
    lk_spleen_liver <- log(1.64e-4 * 60); label("Spleen to liver rate constant, k43 (1/h)") # Table 3, Spleen, F2 row printed 'k31' = 1.64E-4 min^-1; Fig. 3 arrow k43

    # Deep compartment: drawn in Fig. 3 (k15, k51) but no value printed anywhere
    kin_other <- fixed(0); label("Blood to deep compartment rate constant, k15 (1/h)") # Fig. 3 only; unprinted, held at 0
    kout_other <- fixed(0); label("Deep compartment to blood rate constant, k51 (1/h)") # Fig. 3 only; unprinted, held at 0
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    kin_bonemarrow <- exp(lkin_bonemarrow)
    kout_bonemarrow <- exp(lkout_bonemarrow)
    kin_liver <- exp(lkin_liver)
    kout_liver <- exp(lkout_liver)
    kin_spleen <- exp(lkin_spleen)
    k_spleen_liver <- exp(lk_spleen_liver)

    kel <- cl / vc

    d/dt(central) <- -(kel + kin_bonemarrow + kin_liver + kin_spleen + kin_other) * central +
      kout_bonemarrow * bonemarrow + kout_liver * liver + kout_other * other
    d/dt(bonemarrow) <- kin_bonemarrow * central - kout_bonemarrow * bonemarrow
    d/dt(liver) <- kin_liver * central + k_spleen_liver * spleen - kout_liver * liver
    d/dt(spleen) <- kin_spleen * central - k_spleen_liver * spleen
    d/dt(other) <- kin_other * central - kout_other * other

    Cc <- central / vc
  })
}
