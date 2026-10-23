Porchet_1987_nicotine_rabbit_pbpk <- function() {
  description <- paste(
    "Preclinical (rabbit). PBPK (physiologically parameterised three-compartment",
    "mammillary model). Intravenous nicotine in anaesthetised New Zealand white",
    "rabbits: arterial blood is the central compartment, exchanging with a",
    "lumped peripheral (musculoskeletal and visceral) tissue compartment and a",
    "brain compartment. Each tissue volume is its physical volume times a",
    "tissue-to-blood partition coefficient, each transfer rate is the organ",
    "blood flow divided by the volume it leaves, and elimination is from",
    "arterial blood. Volumes and flows are literature values; total clearance",
    "and the two partition coefficients were fitted to the mean of six rabbits.",
    "Concentrations in the two tissues are reported as the concentration in the",
    "venous blood draining them (femoral and internal jugular vein). No",
    "between-subject or residual variability was estimated."
  )
  reference <- paste(
    "Porchet HC, Benowitz NL, Sheiner LB, Copeland JR.",
    "Apparent tolerance to the acute effect of nicotine results in part from",
    "distribution kinetics. J Clin Invest. 1987;80(5):1466-1471.",
    "doi:10.1172/JCI113227."
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
    species = "rabbit (New Zealand white)",
    n_subjects = 6L,
    n_studies = 1L,
    weight_range = "3.0-4.2 kg",
    weight_mean = "3.6 kg",
    sex_female_pct = 0,
    disease_state = "healthy, anaesthetised with ketamine and diazepam",
    dose_range = "25 ug/kg nicotine base (as bitartrate) intravenously over 1 min",
    notes = paste(
      "Methods, Subjects and Experimental protocol. Blood sampled simultaneously",
      "from the femoral artery, femoral vein and internal jugular vein at 1-60",
      "min. The model was fitted to the mean of the six animals, not to",
      "individual data, so it describes a typical 3.6 kg rabbit."
    )
  )

  ini({
    # Fitted parameters (Table II, rabbit column, '*' = primary parameter
    # obtained by fitting the data). Volumes in L and flows in L/min so that
    # ug / L = ng/mL.
    lcl <- log(0.055)
    label("Total (arterial-blood) clearance (L/min)") # Table II: total clearance 55 ml/min
    lkp_brain <- log(0.779)
    label("Brain-to-blood partition coefficient (unitless)") # Table II: PCbrain 0.779
    lkp_other <- log(0.729)
    label("Peripheral tissue-to-blood partition coefficient (unitless)") # Table II: PCperiph 0.729

    # Physical volumes and blood flows, fixed at literature values (Table I,
    # rabbit column).
    lvc <- fixed(log(0.225))
    label("Total blood volume, the central (arterial) compartment volume (L)") # Table I: total blood volume 225 ml
    lv_other <- fixed(log(3.587))
    label("Peripheral tissue physical volume (L)") # Table I: peripheral tissue volume 3,587 ml (total body 3,600 minus brain 13)
    lv_brain <- fixed(log(0.013))
    label("Brain physical volume (L)") # Table I: brain volume 13 ml
    lq_other <- fixed(log(0.760))
    label("Peripheral tissue blood flow (L/min)") # Table I: peripheral blood flow 760 ml/min
    lq_brain <- fixed(log(0.010))
    label("Brain blood flow (L/min)") # Table I: brain blood flow 10 ml/min
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
    # carries amount / apparent volume, so Cother is femoral venous blood and
    # Cbrain is internal jugular venous blood (the paper's 'brain' curve).
    Cc <- central / vc
    Cother <- other / vapp_other
    Cbrain <- brain / vapp_brain
  })
}
