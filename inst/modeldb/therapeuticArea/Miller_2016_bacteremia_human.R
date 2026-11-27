Miller_2016_bacteremia_human <- function() {
  description <- paste(
    "Mechanistic five-compartment systems model of symptomatic",
    "A. baumannii bacteremia in an adult human with an optional",
    "extracorporeal magnetic bacterial-separation device (Miller 2016),",
    "extrapolated from the authors' rodent model by substituting human",
    "organ volumes, blood flows and device size. Bacterial concentrations",
    "(CFU/mL) in lungs, spleen, liver, blood and the separation device are",
    "linked by organ blood flow, organ volume and a bacterial blood:organ",
    "partition coefficient; each organ has a first-order net bacterial",
    "growth rate. Deterministic: no IIV, no residual error. Default",
    "parameters are the as-run immunocompromised + colistin scenario",
    "(paper Figure 6) with the device switched off (fsep = 0); set",
    "fsep = 1 for the 100% efficient device, and halve every kgrow_* for",
    "the multi-drug-resistant scenario of Figure 7.",
    sep = " "
  )
  reference <- paste(
    "Miller SE, Bell CS, McClain MS, Cover TL, Giorgio TD. (2016).",
    "Dynamic Computational Model of Symptomatic Bacteremia to Inform",
    "Bacterial Separation Treatment Requirements.",
    "PLoS One 11(9):e0163167. doi:10.1371/journal.pone.0163167.",
    "Structure and constants verified against the authors' MATLAB source",
    "code (S11 File).",
    sep = " "
  )
  vignette <- "Miller_2016_bacteremia"

  paper_specific_compartments <- c(
    "bact_lung",
    "bact_spleen",
    "bact_liver",
    "bact_device",
    "bact_blood"
  )

  units <- list(
    time = "h",
    dosing = "none (initial-condition model; the bacterial burden at presentation is set through bact0_* and colistin acts through the kgrow_* values)",
    concentration = "CFU/mL"
  )

  compartmentData <- list(
    bact_lung = list(analyte = "Acinetobacter baumannii", units = "CFU/mL", specimen = "tissue", verified = TRUE),
    bact_spleen = list(analyte = "Acinetobacter baumannii", units = "CFU/mL", specimen = "tissue", verified = TRUE),
    bact_liver = list(analyte = "Acinetobacter baumannii", units = "CFU/mL", specimen = "tissue", verified = TRUE),
    bact_device = list(
      analyte = "Acinetobacter baumannii captured in the extracorporeal separation device",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    bact_blood = list(analyte = "Acinetobacter baumannii", units = "CFU/mL", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    weight_range = "Adult with 5000 mL total blood volume (no body weight given)",
    disease_state = paste(
      "Symptomatic A. baumannii bacteremia (blood 10 CFU/mL, lungs 1e3",
      "CFU/mL at presentation) in an immunocompromised adult treated with",
      "colistin methanesulfonate, with or without extracorporeal bacterial",
      "separation. A hypothetical extrapolation; no human data were fitted."
    ),
    dose_range = "Colistin methanesulfonate (represented only through the net growth rates; no human dose given).",
    notes = paste(
      "Organ volumes and specific blood flows are literature human values",
      "(Table 8). The net growth rates and partition coefficients are the",
      "rodent values; the paper says the human growth rates were",
      "extrapolated from the rodent model but prints none, and the rodent",
      "Table 12 values reproduce the human AUCs of Table 6 to three",
      "significant figures."
    )
  )

  ini({
    # Net bacterial growth rates, 1/h: A. baumannii, immunocompromised +
    # colistin. No human-specific values are printed; the rodent values ran.
    kgrow_lung <- fixed(-0.24)
    label("Net bacterial growth rate in lungs (1/h)") # Table 12, Lung (GL) = -0.24 (rodent value used for human)
    kgrow_spleen <- fixed(-0.07)
    label("Net bacterial growth rate in spleen (1/h)") # Table 12, Spleen (GS) = -0.07 (rodent value used for human)
    kgrow_liver <- fixed(-0.18)
    label("Net bacterial growth rate in liver (1/h)") # Table 12, Liver (GH) = -0.18 (rodent value used for human)
    kgrow_blood <- fixed(-0.05)
    label("Net bacterial growth rate in blood (1/h)") # Table 12, Blood (GB) = -0.05 (rodent value used for human)

    # Bacterial organ:blood partition coefficients (unitless).
    lkp_lung <- fixed(log(3))
    label("Bacterial partition coefficient, lungs (unitless)") # Table 9, Lungs, Immunocompromised = 3
    lkp_spleen <- fixed(log(28))
    label("Bacterial partition coefficient, spleen (unitless)") # Table 9, Spleen, Immunocompromised = 28
    lkp_liver <- fixed(log(79))
    label("Bacterial partition coefficient, liver (unitless)") # Table 9, Liver, Immunocompromised = 79
    lkp_device <- fixed(log(1))
    label("Bacterial partition coefficient, separation device (unitless)") # Table 9, Device = 1 (both columns)

    # Organ volumes (g, taken as mL) and specific blood flows (mL/min/g).
    v_lung <- fixed(1315)
    label("Lung volume (mL)") # Table 8, Lungs V = 1315
    v_spleen <- fixed(127)
    label("Spleen volume (mL)") # Table 8, Spleen V = 127
    v_liver <- fixed(830)
    label("Liver volume (mL)") # Table 8, Liver V = 830
    v_device <- fixed(90)
    label("Separation-device volume (mL)") # Table 8, Device V = 90; S11 code V_E = 90
    qperf_lung <- fixed(1.14)
    label("Lung blood flow per gram (mL/min/g)") # Table 8, Lungs Q = 1.14
    qperf_spleen <- fixed(1.97)
    label("Spleen blood flow per gram (mL/min/g)") # Table 8, Spleen Q = 1.97
    qperf_liver <- fixed(0.97)
    label("Liver blood flow per gram (mL/min/g)") # Table 8, Liver Q = 0.97

    # Extracorporeal separation device.
    qbfr <- fixed(1000)
    label("Blood flow through the separation device (mL/min, as run)") # S11 code Q_E = 1000, multiplied by 60 like the organ flows; text: one-fifth of 5000 mL blood volume
    fsep <- fixed(0)
    label("Bacterial separation efficiency f(rf) (fraction; 0 = device off, 1 = 100%)") # eqs 4-5, Fig 6 (100%); default 0 = antibiotic alone

    # Bacterial concentrations at presentation (simulation start), CFU/mL.
    bact0_lung <- fixed(1000)
    label("Initial lung bacterial concentration (CFU/mL)") # Results: symptomatic bacteremia at 1e3 CFU/mL in the lungs
    bact0_spleen <- fixed(757)
    label("Initial spleen bacterial concentration (CFU/mL)") # back-solved from the eight AUCs of Table 6 (not printed)
    bact0_liver <- fixed(761)
    label("Initial liver bacterial concentration (CFU/mL)") # back-solved from the eight AUCs of Table 6 (not printed)
    bact0_blood <- fixed(10)
    label("Initial blood bacterial concentration (CFU/mL)") # Results: symptomatic bacteremia at 10 CFU/mL in the bloodstream
  })

  model({
    kp_lung <- exp(lkp_lung)
    kp_spleen <- exp(lkp_spleen)
    kp_liver <- exp(lkp_liver)
    kp_device <- exp(lkp_device)

    # Organ and device blood flows in mL/h; the source code multiplies each
    # Q/V term by 60 to convert the per-minute flows to the per-hour time base.
    q_lung <- qperf_lung * v_lung * 60
    q_spleen <- qperf_spleen * v_spleen * 60
    q_liver <- qperf_liver * v_liver * 60
    q_device <- qbfr * 60

    # Equations 1-5. Exchange terms are flow over the DONOR or RECEIVING
    # organ volume, as printed and as coded; blood volume does not appear.
    d/dt(bact_lung) <- kgrow_lung * bact_lung +
      q_lung / v_lung * bact_blood -
      q_lung / (v_lung * kp_lung) * bact_lung
    d/dt(bact_spleen) <- kgrow_spleen * bact_spleen +
      q_spleen / v_spleen * bact_blood -
      q_spleen / (v_spleen * kp_spleen) * bact_spleen
    d/dt(bact_liver) <- kgrow_liver * bact_liver +
      q_liver / v_liver * bact_blood +
      q_spleen / (v_spleen * kp_spleen) * bact_spleen -
      (q_liver + q_spleen) / (v_liver * kp_liver) * bact_liver
    d/dt(bact_device) <- q_device / v_device * bact_blood * fsep -
      q_device / (v_device * kp_device) * bact_device * (1 - fsep)
    d/dt(bact_blood) <- kgrow_blood * bact_blood +
      (q_liver + q_spleen) / (v_liver * kp_liver) * bact_liver +
      q_lung / (v_lung * kp_lung) * bact_lung +
      q_device / (v_device * kp_device) * bact_device * (1 - fsep) -
      (q_liver / v_liver + q_spleen / v_spleen + q_lung / v_lung +
        fsep * q_device / v_device) * bact_blood

    bact_lung(0) <- bact0_lung
    bact_spleen(0) <- bact0_spleen
    bact_liver(0) <- bact0_liver
    bact_device(0) <- 0
    bact_blood(0) <- bact0_blood
  })
}
