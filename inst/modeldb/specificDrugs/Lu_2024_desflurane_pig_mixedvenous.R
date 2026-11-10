Lu_2024_desflurane_pig_mixedvenous <- function() {
  description <- paste(
    "Preclinical (pig). Two-compartment model with zero-order input for the",
    "pulmonary-artery (mixed venous) blood desflurane concentration (percent",
    "vapour, i.e. partial pressure as volume percent) in mechanically",
    "ventilated piglets breathing 6 percent inspiratory desflurane for 2 h,",
    "followed by a 2-h washout (Lu 2024). Lu 2024 fitted each piglet",
    "individually in WinNonlin with constant-rate (zero-order) input over the",
    "120-min uptake phase and reported the mean of the individual",
    "biexponential estimates (Table 3); the volumes and clearances carried",
    "here are reconstructed from those macro-constants with the clearance",
    "anchored to the printed mean, so no between-subject variability or",
    "residual-error model is available and the residual SD is fixed at zero.",
    "The input is the nominal amount of the WinNonlin fit, in percent x litre:",
    "a 120-min zero-order infusion of 52.291 percent x litre into central",
    "reproduces the 6 percent inspiratory exposure (the amount scales only the",
    "clearance and volumes, not the predicted concentrations). Lu 2024 fitted",
    "the end-tidal and femoral-arterial curves as separate models; see",
    "modellib('Lu_2024_desflurane_pig_endtidal') and",
    "modellib('Lu_2024_desflurane_pig_arterial')."
  )
  reference <- paste(
    "Lu CC, Ho ST, Hu OYP, Hsiong CH, Cheng YC, Hsu CH, Lin TC (2024).",
    "Pharmacokinetics of desflurane uptake and disposition in piglets.",
    "Front Pharmacol 15:1339690. doi:10.3389/fphar.2024.1339690"
  )
  vignette <- "Lu_2024_desflurane"
  units <- list(time = "min", dosing = "%*L", concentration = "%")

  compartmentData <- list(
    central = list(analyte = "desflurane", units = "%*L", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "desflurane", units = "%*L", specimen = "whole blood", verified = TRUE)
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
      "Lu 2024 Sections 2.2-2.5: seven male domestic piglets, fasted",
      "overnight. Mixed venous blood was drawn through a pulmonary-artery",
      "catheter (with simultaneous femoral-arterial samples) at 1, 3, 5, 10, 20,",
      "30, 40, 50, 60, 80, 100 and 120 min of each phase and assayed by",
      "headspace gas chromatography; concentrations are reported as percent",
      "vapour (partial pressure). The blood-gas partition coefficient measured",
      "before dosing was 0.38 +/- 0.01. No covariates were investigated."
    )
  )

  ini({
    # Two-compartment zero-order-input fit. Lu 2024 Table 3 (pulmonary artery
    # blood column) gives, as means of the seven individual fits, the
    # biexponential Y = A*exp(-alpha*t) + B*exp(-beta*t) with A = 147.32 %,
    # B = 2.10 %, alpha half-life 1.86 min and beta half-life 61.29 min, and
    # the clearance CL = 0.09 L/min. The micro-constants follow from the
    # macro-constants (alpha = ln(2)/1.86, beta = ln(2)/61.29):
    #   k21 = (A*beta + B*alpha)/(A + B) = 0.0163879 /min
    #   k10 = alpha*beta/k21             = 0.257174 /min
    #   k12 = alpha + beta - k21 - k10   = 0.110408 /min
    # and the clearance anchors the volumes: vc = CL/k10, q = k12*vc,
    # vp = q/k21. The implied nominal input is vc*(A + B) = 52.291 % x L.
    # These regenerate Table 3's central volume 0.35 L (0.350), Vss = 2.75 L
    # (2.708) and MRT = 31.26 min (30.09).
    lcl <- fixed(log(0.09)); label("Clearance (L/min)") # Table 3, pulmonary artery 'Clearance, l.min-1' = 0.09
    lvc <- fixed(log(0.3499581)); label("Central volume of distribution (L)") # CL/k10 from Table 3 pulmonary artery A, B, alpha and beta half-lives; Table 3 prints 0.35
    lq <- fixed(log(0.03863803)); label("Intercompartmental clearance (L/min)") # k12*vc from Table 3 pulmonary artery macro-constants
    lvp <- fixed(log(2.357724)); label("Peripheral volume of distribution (L)") # q/k21 from Table 3 pulmonary artery macro-constants

    addSd <- fixed(0); label("Additive residual SD on mixed venous blood desflurane, not reported (%)") # Lu 2024 reports no residual-error model
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
