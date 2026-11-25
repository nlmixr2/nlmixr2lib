Miller_2016_bacteremia_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). Mechanistic five-compartment systems model of",
    "symptomatic Gram-negative bacteremia (Acinetobacter baumannii or",
    "Klebsiella pneumoniae) with an optional extracorporeal magnetic",
    "bacterial-separation device (Miller 2016). Bacterial concentrations",
    "(CFU/mL) in lungs, spleen, liver, blood and the separation device are",
    "linked by organ blood flow, organ volume and a bacterial blood:organ",
    "partition coefficient; each organ has a first-order net bacterial",
    "growth rate (proliferation minus immune clearance) that changes with",
    "immune status and colistin treatment. Deterministic: no IIV, no",
    "residual error. Default parameters are the as-run A. baumannii",
    "immunocompromised + colistin scenario (paper Figure 3) with the device",
    "switched off (fsep = 0); set fsep = 1 for the 100% efficient device",
    "(Figure 4). Other scenarios (immunonormal, untreated neutropenic,",
    "K. pneumoniae) are reproduced by overriding kgrow_*, lkp_* and bact0_*;",
    "see the vignette.",
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
    dosing = "none (initial-condition model; the bacterial inoculum is set through bact0_* and colistin acts through the kgrow_* values)",
    concentration = "CFU/mL"
  )

  compartmentData <- list(
    bact_lung = list(
      analyte = "Gram-negative bacteria (A. baumannii or K. pneumoniae)",
      units = "CFU/mL",
      specimen = "tissue",
      verified = TRUE
    ),
    bact_spleen = list(
      analyte = "Gram-negative bacteria (A. baumannii or K. pneumoniae)",
      units = "CFU/mL",
      specimen = "tissue",
      verified = TRUE
    ),
    bact_liver = list(
      analyte = "Gram-negative bacteria (A. baumannii or K. pneumoniae)",
      units = "CFU/mL",
      specimen = "tissue",
      verified = TRUE
    ),
    bact_device = list(
      analyte = "Gram-negative bacteria captured in the extracorporeal separation device",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    bact_blood = list(
      analyte = "Gram-negative bacteria (A. baumannii or K. pneumoniae)",
      units = "CFU/mL",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "mouse",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    weight_range = "31 g reference mouse (total blood volume 1.80 mL at 58.5 mL/kg)",
    disease_state = paste(
      "Intratracheal Gram-negative challenge (1e7 CFU/mL A. baumannii or",
      "K. pneumoniae) in immunonormal or neutropenic (immunocompromised)",
      "rodents, untreated or treated with colistin methanesulfonate 3 mg/kg."
    ),
    dose_range = "Colistin methanesulfonate 3 mg/kg (represented only through the net growth rates of Tables 12-13).",
    notes = paste(
      "Not fitted to individual data: organ volumes, blood flows and",
      "partition coefficients are literature physiological values (Tables",
      "7 and 9, from reference 18) and the net growth rates are slopes of",
      "published log10 CFU time courses (Tables 10-13; S5 Fig). Initial",
      "conditions are those of the simulation start, 3 h after",
      "inoculation."
    )
  )

  ini({
    # Net bacterial growth rates (proliferation minus immune clearance), 1/h.
    # Default: A. baumannii, immunocompromised + colistin (Table 12).
    kgrow_lung <- fixed(-0.24)
    label("Net bacterial growth rate in lungs (1/h)") # Table 12, Lung (GL) = -0.24
    kgrow_spleen <- fixed(-0.07)
    label("Net bacterial growth rate in spleen (1/h)") # Table 12, Spleen (GS) = -0.07
    kgrow_liver <- fixed(-0.18)
    label("Net bacterial growth rate in liver (1/h)") # Table 12, Liver (GH) = -0.18
    kgrow_blood <- fixed(-0.05)
    label("Net bacterial growth rate in blood (1/h)") # Table 12, Blood (GB) = -0.05

    # Bacterial organ:blood partition coefficients (unitless). The colistin
    # runs (Figures 3-5, Tables 1-5, S2-S4 Figs, S6-S10 Tables) were made
    # with the column Table 9 labels 'Non-Immunocompromised'; see vignette.
    lkp_lung <- fixed(log(93))
    label("Bacterial partition coefficient, lungs (unitless)") # Table 9, Lungs, Non-Immunocompromised = 93
    lkp_spleen <- fixed(log(59))
    label("Bacterial partition coefficient, spleen (unitless)") # Table 9, Spleen, Non-Immunocompromised = 59
    lkp_liver <- fixed(log(749))
    label("Bacterial partition coefficient, liver (unitless)") # Table 9, Liver, Non-Immunocompromised = 749
    lkp_device <- fixed(log(1))
    label("Bacterial partition coefficient, separation device (unitless)") # Table 9, Device = 1 (both columns)

    # Organ volumes (g, taken as mL) and specific blood flows (mL/min/g).
    v_lung <- fixed(2.10)
    label("Lung volume (mL)") # Table 7, Lungs V = 2.10
    v_spleen <- fixed(0.15)
    label("Spleen volume (mL)") # Table 7, Spleen V = 0.15
    v_liver <- fixed(1.99)
    label("Liver volume (mL)") # Table 7, Liver V = 1.99
    v_device <- fixed(0.0324)
    label("Separation-device channel volume (mL)") # S11 code V_E = 0.002*0.0006*0.027 m^3 = 0.0324 mL; Table 7 rounds to 0.032
    qperf_lung <- fixed(6.84)
    label("Lung blood flow per gram (mL/min/g)") # Table 7, Lungs Q = 6.84
    qperf_spleen <- fixed(1.20)
    label("Spleen blood flow per gram (mL/min/g)") # Table 7, Spleen Q = 1.20
    qperf_liver <- fixed(1.50)
    label("Liver blood flow per gram (mL/min/g)") # Table 7, Liver Q = 1.50

    # Extracorporeal separation device.
    qbfr <- fixed(0.36)
    label("Blood flow through the separation device (mL/min, as run)") # S11 code Q_E = .36, multiplied by 60 like the organ flows; text: one-fifth of 1.80 mL blood volume
    fsep <- fixed(0)
    label("Bacterial separation efficiency f(rf) (fraction; 0 = device off, 1 = 100%)") # eqs 4-5, Table 2 (100%); default 0 = antibiotic alone

    # Initial bacterial concentrations at the simulation start (3 h after
    # inoculation), CFU/mL. Default: A. baumannii, immunocompromised + colistin.
    bact0_lung <- fixed(1e7)
    label("Initial lung bacterial concentration (CFU/mL)") # Results: treatment starts when lung burden reaches 1e7 CFU/mL
    bact0_spleen <- fixed(7.6e6)
    label("Initial spleen bacterial concentration (CFU/mL)") # digitised from Fig 3 at t = 3 h
    bact0_liver <- fixed(1.1e7)
    label("Initial liver bacterial concentration (CFU/mL)") # digitised from Fig 3 at t = 3 h
    bact0_blood <- fixed(6.4e4)
    label("Initial blood bacterial concentration (CFU/mL)") # digitised from Fig 3 at t = 3 h
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
