Lu_2024_desflurane_pig_endtidal <- function() {
  description <- paste(
    "Preclinical (pig). Two-compartment model with zero-order input for the",
    "end-tidal desflurane concentration (percent vapour) in mechanically",
    "ventilated piglets breathing 6 percent inspiratory desflurane for 2 h,",
    "followed by a 2-h washout (Lu 2024). The compartments are those of a",
    "descriptive curve fit observed in exhaled gas, not anatomical spaces.",
    "Lu 2024 fitted each piglet individually in WinNonlin with constant-rate",
    "(zero-order) input over the 120-min uptake phase and reported the mean of",
    "the individual biexponential estimates (Table 3); the volumes and",
    "clearances carried here are reconstructed from those macro-constants",
    "with the clearance anchored to the printed mean, so no between-subject",
    "variability or residual-error model is available and the residual SD is",
    "fixed at zero. The end-tidal means are dominated by between-animal",
    "spread (beta half-life 61.90 +/- 105.93 min), so this typical-value",
    "curve over-predicts the observed mean end-tidal plateau by about half;",
    "prefer the blood models for simulation. The input is the nominal amount",
    "of the WinNonlin fit, in percent x litre: a 120-min zero-order infusion",
    "of 76.953 percent x litre into central corresponds to the 6 percent",
    "inspiratory exposure (the amount scales only the clearance and volumes,",
    "not the predicted concentrations). Lu 2024 fitted the arterial and",
    "pulmonary-artery curves as separate models; see",
    "modellib('Lu_2024_desflurane_pig_arterial') and",
    "modellib('Lu_2024_desflurane_pig_mixedvenous')."
  )
  reference <- paste(
    "Lu CC, Ho ST, Hu OYP, Hsiong CH, Cheng YC, Hsu CH, Lin TC (2024).",
    "Pharmacokinetics of desflurane uptake and disposition in piglets.",
    "Front Pharmacol 15:1339690. doi:10.3389/fphar.2024.1339690"
  )
  vignette <- "Lu_2024_desflurane"
  units <- list(time = "min", dosing = "%*L", concentration = "%")

  # The end-tidal fit's compartments are curve-fit spaces observed through
  # exhaled gas; the specimen vocabulary has no exhaled-gas member.
  compartmentData <- list(
    central = list(analyte = "desflurane", units = "%*L", specimen = "not applicable", verified = TRUE),
    peripheral1 = list(analyte = "desflurane", units = "%*L", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "pig (domestic piglets)",
    n_subjects = 7,
    n_studies = 1,
    age_range = "5 months",
    weight_mean = "28 +/- 5 kg",
    sex_female_pct = 0,
    disease_state = paste(
      "Healthy piglets under general anaesthesia (tiletamine/zolazepam and",
      "atropine induction; fentanyl, propofol and pancuronium maintenance),",
      "tracheostomised and mechanically ventilated with 40 percent oxygen,",
      "tidal volume 8 mL/kg, PEEP 5 cmH2O"
    ),
    dose_range = paste(
      "Desflurane 6 percent inspiratory via a 5 L/min fresh-gas flow for",
      "120 min (uptake), then 40 percent oxygen alone for 120 min (disposition)"
    ),
    regions = "Taiwan",
    notes = paste(
      "Lu 2024 Sections 2.2-2.4: seven male domestic piglets, fasted",
      "overnight. Inspiratory and end-tidal desflurane were recorded",
      "continuously with an infrared multigas analyser; end-tidal values below",
      "0.2 percent were discarded, so the end-tidal curve is truncated after",
      "about 60 min of washout (Section 4.5). No covariates were investigated."
    )
  )

  ini({
    # Two-compartment zero-order-input fit. Lu 2024 Table 3 (end-tidal
    # column) gives, as means of the seven individual fits, the biexponential
    # Y = A*exp(-alpha*t) + B*exp(-beta*t) with A = 731.26 %, B = 3.45 %,
    # alpha half-life 0.75 min and beta half-life 61.90 min, and the
    # clearance CL = 0.07 L/min. The micro-constants follow from the
    # macro-constants (alpha = ln(2)/0.75, beta = ln(2)/61.90):
    #   k21 = (A*beta + B*alpha)/(A + B) = 0.0154851 /min
    #   k10 = alpha*beta/k21             = 0.668323 /min
    #   k12 = alpha + beta - k21 - k10   = 0.251586 /min
    # and the clearance anchors the volumes: vc = CL/k10, q = k12*vc,
    # vp = q/k21. The implied nominal input is vc*(A + B) = 76.953 % x L.
    # These regenerate Table 3's central volume 0.10 L (0.105) but not its
    # Vss = 1.29 L (1.81) or MRT = 15.51 min (25.8): the end-tidal means are
    # pulled by one or more animals with a long beta phase.
    lcl <- fixed(log(0.07)); label("Clearance (L/min)") # Table 3, end-tidal 'Clearance, l.min-1' = 0.07
    lvc <- fixed(log(0.1047398)); label("Central volume of distribution (L)") # CL/k10 from Table 3 end-tidal A, B, alpha and beta half-lives; Table 3 prints 0.10
    lq <- fixed(log(0.02635106)); label("Intercompartmental clearance (L/min)") # k12*vc from Table 3 end-tidal macro-constants
    lvp <- fixed(log(1.70171)); label("Peripheral volume of distribution (L)") # q/k21 from Table 3 end-tidal macro-constants

    addSd <- fixed(0); label("Additive residual SD on end-tidal desflurane, not reported (%)") # Lu 2024 reports no residual-error model
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc

    Cc ~ add(addSd)
  })
}
