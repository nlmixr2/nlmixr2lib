Xiao_2016_tazobactam_pneumonia <- function() {
  description <- paste(
    "Simulation model for tazobactam plasma and epithelial lining fluid",
    "(ELF) concentrations in adults with nosocomial pneumonia, used by Xiao",
    "2016 to justify the phase 3 dose. It combines the two-compartment",
    "plasma population PK of complicated intra-abdominal infection (cIAI)",
    "patients (Chandorkar 2015) with the plasma-to-ELF transfer kinetics",
    "fitted in healthy volunteers (Xiao_2016_tazobactam_healthy): a one-way",
    "central-to-ELF clearance Qce, a one-way ELF-to-central clearance Qec",
    "carrying between-subject variability, and an ELF volume fixed at 1 L.",
    "Creatinine clearance acts on CL through a power function; the cIAI",
    "shift on Vc is built into the typical value. Every value is carried",
    "from the two component models, so all are fixed; the model itself was",
    "never fitted. No residual error was specified for the simulations."
  )
  reference <- paste(
    "Xiao AJ, Miller BW, Huntington JA, Nicolau DP.",
    "Ceftolozane/tazobactam pharmacokinetic/pharmacodynamic-derived dose",
    "justification for phase 3 studies in patients with nosocomial",
    "pneumonia. J Clin Pharmacol. 2016;56(1):56-66. doi:10.1002/jcph.566.",
    "Values: Results, the equations following 'The final PK model used in",
    "the Monte Carlo simulation'. The printed Vc equation is garbled; the",
    "cIAI shift and the Vc variance are restored from the source cIAI model,",
    "Chandorkar G et al., J Clin Pharmacol. 2015;55(2):230-239",
    "(doi:10.1002/jcph.395), Table 3B (see the Vc notes in ini())."
  )
  vignette <- "Xiao_2016_ceftolozane_tazobactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description = "Baseline creatinine clearance",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Xiao 2016 simulation equation 'CL=18.0*(CrCL/115)^0.670*exp[...]',",
        "with CrCL in mL/min, the unit of the Chandorkar 2015 source model,",
        "which estimated CrCL by Cockcroft-Gault on actual body weight. The",
        "normal-renal-function scenario sampled CrCL uniformly on 90-150",
        "mL/min; moderate and severe impairment used 30-50 and 15-30 mL/min."
      ),
      source_name = "CrCL"
    )
  )

  compartmentData <- list(
    central = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tazobactam", units = "mg", specimen = "tissue", verified = TRUE),
    elf = list(analyte = "tazobactam", units = "mg", specimen = "epithelial lining fluid", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = "1000 simulated patients per scenario (no fitted cohort)",
    n_studies = 0L,
    disease_state = paste(
      "Projected adults with nosocomial (hospital-acquired or",
      "ventilator-associated) pneumonia. The plasma PK is that of cIAI",
      "patients from the phase 2 cIAI study pooled in Chandorkar 2015,",
      "chosen because cIAI patients had the larger clearance, volume and",
      "variability; the ELF kinetics are from 25 healthy volunteers."
    ),
    dose_range = paste(
      "Simulated ceftolozane/tazobactam 1.5 g (1000/500 mg) and 3 g",
      "(2000/1000 mg) as a 60-minute infusion every 8 hours for normal renal",
      "function; the tazobactam component is 500 or 1000 mg."
    ),
    renal_function = "Simulated CrCL uniform 90-150 mL/min (normal), 30-50 (moderate), 15-30 (severe)",
    notes = paste(
      "Not fitted to a cohort. The plasma component came from the pooled",
      "phase 1/2 analysis of Chandorkar 2015 and the ELF component from 25",
      "healthy volunteers (Xiao 2016 Table 1), assumed independent of the",
      "plasma PK. Simulations in SAS 9.3 at a 0.001-hour time step, 1000",
      "patients per case. PK/PD index: free-drug %fT>MEC with unbound",
      "fraction 0.70, targets 20-50%."
    )
  )

  ini({
    lcl <- fixed(log(18.0)); label("Tazobactam clearance at CrCL = 115 mL/min (L/h)") # Simulation equation: 'CL=18.0*(CrCL/115)^0.670*exp[N(0, 0.252)]'
    lvc <- fixed(log(14.2)); label("Tazobactam central volume before the cIAI shift (L)") # Simulation equation: 'Vc=14.2*exp[N(0, 0.387)]'
    lq <- fixed(log(3.13)); label("Tazobactam central-peripheral inter-compartmental clearance Q2 (L/h)") # Simulation equation: 'Q2=3.13'
    lvp <- fixed(log(4.29)); label("Tazobactam peripheral volume of distribution V2 (L)") # Simulation equation: 'V2=4.29'
    lq_c_elf <- fixed(log(0.792)); label("Tazobactam one-way central-to-ELF clearance Qce (L/h)") # Simulation equation: 'Qce=0.792'
    lq_elf_c <- fixed(log(1.72)); label("Tazobactam one-way ELF-to-central clearance Qec (L/h)") # Simulation equation: 'Qec=1.72*exp[N(0, 0.229)]'
    lv_elf <- fixed(log(1)); label("ELF compartment volume Ve (L)") # Simulation equation: 'Ve=1'

    e_crcl_cl <- fixed(0.670); label("Power exponent of CrCL on tazobactam CL, CRCL/115 (unitless)") # Simulation equation: '(CrCL/115)^0.670'
    e_ciai_vc <- fixed(0.387); label("cIAI-patient shift on log tazobactam Vc, applied to every simulated patient (log)") # Printed in the variance slot 'N(0, 0.387)'; exp(0.387) = 1.47-fold, the cIAI Vc factor of Chandorkar 2015 Table 3B; Figure 4B median peak needs it

    etalcl ~ fixed(0.252) # Simulation equation: 'N(0, 0.252)' on CL (= 0.502^2, Chandorkar 2015 BSV 50.2%)
    etalvc ~ fixed(0.275625) # Not legible in the printed Vc equation; 0.525^2 from Chandorkar 2015 Table 3B BSV 52.5%, the CV^2 convention of the other three printed variances
    etalq_elf_c ~ fixed(0.229) # Simulation equation: 'N(0, 0.229)' on Qec

    propSd <- fixed(0); label("Proportional residual error, plasma (fraction; none specified for the simulations)") # not reported for the simulation model
    propSd_Celf <- fixed(0); label("Proportional residual error, ELF (fraction; none specified for the simulations)") # not reported for the simulation model
  })

  model({
    cl <- exp(lcl + etalcl) * (CRCL / 115)^e_crcl_cl
    vc <- exp(lvc + etalvc) * exp(e_ciai_vc)
    q <- exp(lq)
    vp <- exp(lvp)
    q_c_elf <- exp(lq_c_elf)
    q_elf_c <- exp(lq_elf_c + etalq_elf_c)
    v_elf <- exp(lv_elf)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k_central_elf <- q_c_elf / vc
    k_elf_central <- q_elf_c / v_elf

    d/dt(central) <- -(kel + k12 + k_central_elf) * central + k21 * peripheral1 + k_elf_central * elf
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(elf) <- k_central_elf * central - k_elf_central * elf

    Cc <- central / vc
    Celf <- elf / v_elf

    Cc ~ prop(propSd)
    Celf ~ prop(propSd_Celf)
  })
}
