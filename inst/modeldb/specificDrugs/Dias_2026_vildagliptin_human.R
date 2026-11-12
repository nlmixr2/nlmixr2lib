Dias_2026_vildagliptin_human <- function() {
  description <- paste(
    "Allometric human projection of the Dias 2026 rat vildagliptin model",
    "for patients with type 2 diabetes, with an Imax model of DPP-4",
    "inhibition by unbound plasma, skeletal-muscle and liver concentrations.",
    "Four-compartment disposition (central, peripheral, muscle and liver,",
    "the tissues linked to central by separate influx and efflux",
    "clearances) whose clearances and volumes are the diabetic-rat per-kg",
    "estimates scaled by Equation 2 (exponent 0.75 for clearances, 1 for",
    "volumes; rat reference 0.35 kg), with oral absorption through two",
    "sequential first-order steps and bioavailability 77.2% taken from",
    "Landersdorfer 2012. Typical-value simulation model: no IIV and no",
    "residual error.",
    sep = " "
  )
  reference <- paste(
    "Dias BB, Olivo LB, de Araujo BV. Preclinical Modeling and Simulation to",
    "Explore the Tissue/Plasma Exposure and Pharmacodynamic Effect of",
    "Vildagliptin in Diabetes Treatment. CPT Pharmacometrics Syst Pharmacol.",
    "2026;15(2):e70165. doi:10.1002/psp4.70165. PMCID: PMC12896384.",
    "Rat parameters: Table 1. Allometric scaling: Equation 2. Human",
    "parameters: Supporting Information Table S2. Absorption: Methods",
    "'Allometry' (from Landersdorfer 2012, doi:10.1111/j.1365-2125.2011.04108.x).",
    "Imax model: Equation 3 and Figure 5.",
    sep = " "
  )
  vignette <- "Dias_2026_vildagliptin"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "vildagliptin", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "vildagliptin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "vildagliptin (unbound)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vildagliptin (unbound)", units = "mg", specimen = "tissue", verified = TRUE),
    muscle = list(
      analyte = "vildagliptin (unbound, skeletal muscle interstitial fluid)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    liver = list(
      analyte = "vildagliptin (unbound, liver interstitial fluid)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight of the patient",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      source_name = "Weight_human",
      notes = paste(
        "Equation 2: Parameter_human = Parameter_rat x (Weight_human /",
        "Weight_rat)^b, b = 0.75 for clearances and 1 for volumes, applied",
        "to the absolute rat parameter (per-kg value x 0.35 kg). The paper",
        "does not print either weight; Supporting Information Table S2 is",
        "reproduced to four significant figures by Weight_human = 91 kg and",
        "Weight_rat = 0.35 kg (91 kg is the mean weight of the Landersdorfer",
        "2012 patients). Simulate the published projection at WT = 91."
      )
    )
  )

  population <- list(
    species = "human (allometric projection from rat)",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    weight_mean = "91 kg (simulated typical patient)",
    disease_state = "Type 2 diabetes mellitus (simulated; no human data were fit).",
    dose_range = "Simulated oral 25, 50 and 100 mg once (q24h) or twice (q12h) daily to steady state.",
    notes = paste(
      "No human data: every disposition parameter is the diabetic-rat",
      "estimate (Table 1) scaled allometrically, assuming diffusion-governed",
      "tissue distribution and the same diabetes effect on the",
      "intercompartmental clearances in rats and humans (Discussion).",
      "Absorption rate constants and bioavailability are the Landersdorfer",
      "2012 estimates in patients with type 2 diabetes."
    )
  )

  ini({
    # Diabetic-rat disposition parameters per kg (Table 1), fixed in the
    # projection and scaled to humans inside model() by Equation 2.
    lcl <- fixed(log(2.72)); label("Rat clearance CL per kg (L/h/kg)") # Table 1: CL = 2.72 L/h/kg
    lvc <- fixed(log(1.11)); label("Rat central volume V1 per kg (L/kg)") # Table 1: V1 = 1.11 L/kg
    lq <- fixed(log(2.61)); label("Rat central-peripheral clearance Q1, diabetic, per kg (L/h/kg)") # Table 1: Q1,diabetic = 2.61 L/h/kg
    lvp <- fixed(log(1.85)); label("Rat peripheral volume V2 per kg (L/kg)") # Table 1: V2 = 1.85 L/kg
    lclin_muscle <- fixed(log(2.09)); label("Rat influx clearance central to muscle per kg (L/h/kg)") # Table 1: Qin,muscle = 2.09 L/h/kg
    lclef_muscle <- fixed(log(14.3)); label("Rat efflux clearance muscle to central per kg (L/h/kg)") # Table 1: Qout,muscle = 14.3 L/h/kg
    lv_muscle <- fixed(log(2.55)); label("Rat muscle compartment volume V3 per kg (L/kg)") # Table 1: V3 = 2.55 L/kg
    lclin_liver <- fixed(log(2.5)); label("Rat influx clearance central to liver per kg (L/h/kg)") # Table 1: Qin,liver = 2.5 L/h/kg
    lclef_liver <- fixed(log(19.7)); label("Rat efflux clearance liver to central, diabetic, per kg (L/h/kg)") # Table 1: Qout,liver,diabetic = 19.7 L/h/kg
    lv_liver <- fixed(log(0.856)); label("Rat liver compartment volume V4 per kg (L/kg)") # Table 1: V4 = 0.856 L/kg

    # Allometric exponents (Methods 'Allometry', after Equation 2).
    e_wt_cl <- fixed(0.75); label("Allometric exponent on every clearance (unitless)") # Methods 'Allometry': b = 0.75 for clearance parameters
    e_wt_vc <- fixed(1); label("Allometric exponent on every volume (unitless)") # Methods 'Allometry': b = 1 for volumes

    # Oral absorption from Landersdorfer 2012 (Methods 'Allometry').
    lktr <- fixed(log(1.26)); label("First absorption step depot to transit, Ka1 (1/h)") # Methods 'Allometry': Ka1 = 1.26 h-1 (Landersdorfer 2012)
    lka <- fixed(log(1.05)); label("Second absorption step transit to central, Ka2 (1/h)") # Methods 'Allometry': Ka2 = 1.05 h-1 (Landersdorfer 2012)
    lfdepot <- fixed(log(0.772)); label("Absolute oral bioavailability Fabs (fraction)") # Methods 'Allometry': Fabs = 77.2% (Landersdorfer 2012)

    # DPP-4 inhibition, Equation 3.
    limax <- fixed(log(100)); label("Maximum DPP-4 inhibition Imax (%)") # Not printed; Figure 5 inhibition plateaus at ~100% for Cu ~ 200 ng/mL
    lic50 <- fixed(log(1.35)); label("Unbound concentration for 50% DPP-4 inhibition IC50 (ng/mL)") # Figure 5 caption: IC50 = 1.35 ng/mL; Methods: 4.5 nmol/L
  })

  model({
    # Equation 2 on the absolute rat parameter (per-kg value x wt_rat).
    # wt_rat = 0.35 kg reproduces Table S2; see covariateData$WT.
    wt_rat <- 0.35
    scale_cl <- wt_rat * (WT / wt_rat)^e_wt_cl
    scale_v <- wt_rat * (WT / wt_rat)^e_wt_vc

    cl <- exp(lcl) * scale_cl
    vc <- exp(lvc) * scale_v
    q <- exp(lq) * scale_cl
    vp <- exp(lvp) * scale_v
    clin_muscle <- exp(lclin_muscle) * scale_cl
    clef_muscle <- exp(lclef_muscle) * scale_cl
    v_muscle <- exp(lv_muscle) * scale_v
    clin_liver <- exp(lclin_liver) * scale_cl
    clef_liver <- exp(lclef_liver) * scale_cl
    v_liver <- exp(lv_liver) * scale_v

    ktr <- exp(lktr)
    ka <- exp(lka)
    imax <- exp(limax)
    ic50 <- exp(lic50)

    # Two sequential first-order absorption steps, then Equations 4-7.
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ka * transit1
    d/dt(central) <- ka * transit1 + q / vp * peripheral1 + clef_muscle / v_muscle * muscle + clef_liver / v_liver * liver - (cl / vc + q / vc + clin_muscle / vc + clin_liver / vc) * central
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1
    d/dt(muscle) <- clin_muscle / vc * central - clef_muscle / v_muscle * muscle
    d/dt(liver) <- clin_liver / vc * central - clef_liver / v_liver * liver

    f(depot) <- exp(lfdepot)

    # Unbound concentrations (mg/L x 1000 = ng/mL) and Equation 3.
    Cu <- central / vc * 1000
    Cmuscle <- muscle / v_muscle * 1000
    Cliver <- liver / v_liver * 1000
    inh_plasma <- imax * Cu / (Cu + ic50)
    inh_muscle <- imax * Cmuscle / (Cmuscle + ic50)
    inh_liver <- imax * Cliver / (Cliver + ic50)
  })
}
