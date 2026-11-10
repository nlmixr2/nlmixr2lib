Yang_2025b_iruplinalkib_human <- function() {
  description <- paste(
    "Animal-scaled human projection (not fitted to human data). Two-compartment",
    "PK model with first-order absorption and first-order elimination for oral",
    "iruplinalkib (WX-0593, a selective ALK/ROS1 tyrosine kinase inhibitor),",
    "built before the first clinical data were available to predict the human",
    "effective dose. Mouse, rat, dog and monkey intravenous profiles were",
    "superimposed with Wajima's Css-MRT normalisation and rescaled to a",
    "60 kg human using predicted Vss (7.9 L/kg) and CL (311 mL/min) averaged",
    "over several interspecies methods; a two-compartment model was fitted to",
    "the rescaled profile and combined with the mean animal oral",
    "bioavailability (59.2%) and the monkey Ka (0.341 1/h) to give apparent",
    "oral parameters. Typical-value model with no variability."
  )

  reference <- paste(
    "Yang G, Wang Y, Zhao H, Zheng Q, Wang X, Jia L, Xin Q, Ma C, Zhang Y,",
    "Zheng S, Kang X.",
    "Translational Model-Informed Dose Selection for Iruplinalkib, a",
    "Selective Oral ALK/ROS1 Tyrosine Kinase Inhibitor.",
    "Clin Transl Sci. 2025;18(7):e70287.",
    "doi:10.1111/cts.70287"
  )

  vignette <- "Yang_2025b_iruplinalkib"

  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "iruplinalkib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "iruplinalkib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "iruplinalkib", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human (projected from mouse, rat, dog and monkey data)",
    n_subjects = NA_integer_,
    weight_range = "60 kg typical human assumed for every calculation (Methods S1)",
    disease_state = "Hypothetical typical adult; no human data were used",
    dose_range = paste(
      "Projected oral once-daily regimens of 1, 2.4 and 3 mg/kg (60, 144 and",
      "180 mg for a 60 kg subject; Figure 3C and Table S2)"
    ),
    notes = paste(
      "Human Vss was the mean of four predictions (Oie-Tozer 6.2, monkey",
      "proportionality 7.7, allometric scaling 8.4 and Wajima 9.2 L/kg) and",
      "human CL the mean of four (hepatic-blood-flow-corrected monkey 336,",
      "fu-corrected intercept method in dogs and monkeys 349 and 357, and",
      "rat-dog allometry 203 mL/min). For a 60 kg human the predicted Css was",
      "0.1266 ug/mL and MRT 25.4 h after a 1 mg/kg IV dose. Fitted in Phoenix",
      "NLME 8.3.4 (FOCE-ELS). The 59.2% bioavailability is the mean gavage",
      "bioavailability in mice, rats and monkeys."
    )
  )

  ini({
    # Table 2 'Predicted PK parameters of orally administered iruplinalkib in
    # humans based on animal study data'. The disposition parameters are
    # APPARENT oral values (CL/F, V1/F, Q/F, V2/F): the 59.2% bioavailability
    # is already divided into them, so no separate F is applied. Simulating
    # 180 mg QD with these values reproduces the paper's preclinical
    # projection in Table S2 (Ctrough,ss 126 and Cmax,ss 305 ng/mL); see the
    # vignette.
    lka <- fixed(log(0.341)); label("First-order absorption rate constant Ka (1/h)")    # Table 2: Ka = 0.341 1/h (FIX); Methods S1: monkey Ka (0.34 1/h) used as human Ka
    lcl <- log(33.7); label("Apparent oral clearance CL/F (L/h)")                       # Table 2: CL/F = 33.7 L/h (RSE 1.7%)
    lvc <- log(379); label("Apparent central volume of distribution V1/F (L)")          # Table 2: V1/F = 379 L (RSE 9.6%)
    lq <- log(134); label("Apparent intercompartmental clearance Q/F (L/h)")            # Table 2: Q/F = 134 L/h (RSE 29.7%)
    lvp <- log(239); label("Apparent peripheral volume of distribution V2/F (L)")       # Table 2: V2/F = 239 L (RSE 8.5%)
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L = ug/mL
    Cc <- central / vc
  })
}
