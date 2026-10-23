Porchet_1987_nicotine_human_pbpk <- function() {
  description <- paste(
    "PBPK (physiologically parameterised three-compartment mammillary model).",
    "Intravenous nicotine in healthy adult male cigarette smokers, scaled from",
    "the same authors' rabbit model: arterial blood is the central compartment,",
    "exchanging with a lumped peripheral (musculoskeletal and visceral) tissue",
    "compartment and a brain compartment. Each tissue volume is its physical",
    "volume times a tissue-to-blood partition coefficient, each transfer rate",
    "is the organ blood flow divided by the volume it leaves, and elimination",
    "is from arterial blood. Volumes and flows are literature values for a",
    "70 kg man, the brain partition coefficient is carried over from the",
    "rabbit, clearance is dose / AUC of forearm venous blood, and only the",
    "peripheral partition coefficient was fitted to the mean forearm venous",
    "concentrations of seven men. Forearm venous blood is the peripheral",
    "venous output Cother; Cbrain is the predicted 'brain' concentration that",
    "the paper uses as the effect-site driver of heart rate. No between-subject",
    "or residual variability was estimated."
  )
  reference <- paste(
    "Porchet HC, Benowitz NL, Sheiner LB, Copeland JR.",
    "Apparent tolerance to the acute effect of nicotine results in part from",
    "distribution kinetics. J Clin Invest. 1987;80(5):1466-1471.",
    "doi:10.1172/JCI113227.",
    "Brain partition coefficient from the rabbit fit; see",
    "modellib('Porchet_1987_nicotine_rabbit_pbpk')."
  )
  vignette <- "Porchet_1987_nicotine_pbpk"
  units <- list(time = "min", dosing = "ug", concentration = "ng/mL")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "nicotine", units = "ug", specimen = "whole blood", verified = TRUE),
    other = list(analyte = "nicotine", units = "ug", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "nicotine", units = "ug", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 7L,
    n_studies = 1L,
    age_range = "24-60 years",
    age_mean = "37 years",
    sex_female_pct = 0,
    disease_state = "healthy regular cigarette smokers (30-40 cigarettes/day), normal liver and kidney function",
    dose_range = paste(
      "IV nicotine base (as bitartrate) on two days: rapid loading 17.5 ug/kg/min",
      "for 1.5 min then 1.75 ug/kg/min to 30 min; slow loading 2.5 ug/kg/min for",
      "30 min; both followed by 0.35 ug/kg/min from 30 to 180 min"
    ),
    regions = "United States (San Francisco General Hospital)",
    notes = paste(
      "Methods, Subjects. Body weights are not reported; the physiological",
      "volumes and flows (Table I) are those of a 70 kg man. The model was",
      "fitted to the mean concentration-time curves, so it describes a typical",
      "subject."
    )
  )

  ini({
    # Table II, human column. Clearance is '*' (primary) but was fixed to dose
    # / AUC(0-inf) of forearm venous blood rather than fitted (Methods); the
    # Results give 1,323 and 1,260 ml/min for the two regimens.
    lcl <- log(1.300)
    label("Total (arterial-blood) clearance (L/min)") # Table II: total clearance 1,300 ml/min
    lkp_brain <- fixed(log(0.779))
    label("Brain-to-blood partition coefficient (unitless)") # Table II: PCbrain 0.779, assumed equal to the rabbit value (Methods)
    lkp_other <- log(1.44)
    label("Peripheral tissue-to-blood partition coefficient (unitless)") # Table II: PCperiph 1.44, the only fitted human parameter

    # Physical volumes and blood flows, fixed at literature values (Table I,
    # human column).
    lvc <- fixed(log(5.6))
    label("Total blood volume, the central (arterial) compartment volume (L)") # Table I: total blood volume 5,600 ml
    lv_other <- fixed(log(68.6))
    label("Peripheral tissue physical volume (L)") # Table I: peripheral tissue volume 68,600 ml (total body 70,000 minus brain 1,400)
    lv_brain <- fixed(log(1.4))
    label("Brain physical volume (L)") # Table I: brain volume 1,400 ml
    lq_other <- fixed(log(4.65))
    label("Peripheral tissue blood flow (L/min)") # Table I: peripheral blood flow 4,650 ml/min
    lq_brain <- fixed(log(0.75))
    label("Brain blood flow (L/min)") # Table I: brain blood flow 750 ml/min
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    kp_brain <- exp(lkp_brain)
    kp_other <- exp(lkp_other)
    v_other <- exp(lv_other)
    v_brain <- exp(lv_brain)
    q_other <- exp(lq_other)
    q_brain <- exp(lq_brain)

    # Apparent volumes: physical volume times partition coefficient (Methods,
    # Pharmacokinetic analysis).
    vapp_other <- v_other * kp_other
    vapp_brain <- v_brain * kp_brain

    # k_ij = blood flow from i to j / apparent volume of i; k10 = CL / V1.
    kel <- cl / vc
    k12 <- q_other / vc
    k21 <- q_other / vapp_other
    k13 <- q_brain / vc
    k31 <- q_brain / vapp_brain

    d/dt(central) <- -(kel + k12 + k13) * central + k21 * other + k31 * brain
    d/dt(other) <- k12 * central - k21 * other
    d/dt(brain) <- k13 * central - k31 * brain

    # ug / L = ng/mL. Under venous equilibration the blood leaving a tissue
    # carries amount / apparent volume, so Cother is the forearm venous blood
    # the human model was fitted to and Cbrain is the paper's predicted
    # 'brain' concentration (Figures 6 and 7).
    Cc <- central / vc
    Cother <- other / vapp_other
    Cbrain <- brain / vapp_brain
  })
}
