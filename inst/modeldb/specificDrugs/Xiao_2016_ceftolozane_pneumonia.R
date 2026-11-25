Xiao_2016_ceftolozane_pneumonia <- function() {
  description <- paste(
    "Simulation model for ceftolozane plasma and epithelial lining fluid",
    "(ELF) concentrations in adults with nosocomial pneumonia, used by Xiao",
    "2016 to justify the phase 3 dose. It combines the two-compartment",
    "plasma population PK of complicated intra-abdominal infection (cIAI)",
    "patients (Chandorkar 2015) with the plasma-to-ELF transfer kinetics",
    "fitted in healthy volunteers (Xiao_2016_ceftolozane_healthy): a",
    "one-way central-to-ELF clearance Qce carrying between-subject",
    "variability, a one-way ELF-to-central clearance Qec, and an ELF volume",
    "fixed at 1 L. Creatinine clearance acts on CL through a power function",
    "and body weight proportionally on central volume; the cIAI shifts on CL",
    "and Vc are built into the typical values. Every value is carried from",
    "the two component models, so all are fixed; the model itself was never",
    "fitted. No residual error was specified for the simulations."
  )
  reference <- paste(
    "Xiao AJ, Miller BW, Huntington JA, Nicolau DP.",
    "Ceftolozane/tazobactam pharmacokinetic/pharmacodynamic-derived dose",
    "justification for phase 3 studies in patients with nosocomial",
    "pneumonia. J Clin Pharmacol. 2016;56(1):56-66. doi:10.1002/jcph.566.",
    "All values: Results, the equations following 'The final PK model used",
    "in the Monte Carlo simulation'. The plasma part reproduces the cIAI",
    "model of Chandorkar G et al., J Clin Pharmacol. 2015;55(2):230-239",
    "(doi:10.1002/jcph.395), Table 3A."
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
        "Xiao 2016 simulation equation 'CL=5.11*(CrCL/109)^0.715*exp[...]',",
        "where 'CrCL represents renal clearance in millimeters per minute' --",
        "read as mL/min, the unit of the Chandorkar 2015 source model, which",
        "estimated CrCL by Cockcroft-Gault on actual body weight. The",
        "normal-renal-function scenario sampled CrCL uniformly on 90-150",
        "mL/min; moderate and severe impairment used 30-50 and 15-30 mL/min."
      ),
      source_name = "CrCL"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Simulation equation 'Vc=11.4*(WT/74)*exp[...]'. Body weight was",
        "sampled from a normal distribution with mean 74 kg and a printed",
        "spread of 0.205 kg ('mean (standard error)'). The source cIAI model",
        "(Chandorkar 2015) switches the weight effect off in cIAI patients;",
        "Xiao 2016 printed and simulated it with the weight term retained,",
        "which is what this model encodes."
      ),
      source_name = "WT"
    )
  )

  compartmentData <- list(
    central = list(analyte = "ceftolozane", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ceftolozane", units = "mg", specimen = "tissue", verified = TRUE),
    elf = list(analyte = "ceftolozane", units = "mg", specimen = "epithelial lining fluid", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = "1000 simulated patients per scenario (no fitted cohort)",
    n_studies = 0L,
    weight_median = "74 kg (simulated)",
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
      "function; 1.5 g and 750 mg for moderate and severe renal impairment."
    ),
    renal_function = "Simulated CrCL uniform 90-150 mL/min (normal), 30-50 (moderate), 15-30 (severe)",
    notes = paste(
      "Not fitted to a cohort. The plasma component came from 8 phase 1 and",
      "2 phase 2 studies (Chandorkar 2015) and the ELF component from 25",
      "healthy volunteers (Xiao 2016 Table 1), assumed independent of the",
      "plasma PK. Simulations in SAS 9.3 at a 0.001-hour time step, 1000",
      "patients per case. PK/PD target: free-drug %fT>MIC with unbound",
      "fraction 0.79; 24.8% (stasis), 32.2% (1-log kill), 40% and 50%."
    )
  )

  ini({
    lcl <- fixed(log(5.11)); label("Ceftolozane clearance at CrCL = 109 mL/min before the cIAI shift (L/h)") # Simulation equation: 'CL=5.11*(CrCL/109)^0.715*exp[0.195+N(0, 109)]'
    lvc <- fixed(log(11.4)); label("Ceftolozane central volume at 74 kg before the cIAI shift (L)") # Simulation equation: 'Vc=11.4*(WT/74)*exp[0.464+N(0, 0.158)]'
    lq <- fixed(log(1.19)); label("Ceftolozane central-peripheral inter-compartmental clearance Q2 (L/h)") # Simulation equation: 'Q2=1.19'
    lvp <- fixed(log(2.88)); label("Ceftolozane peripheral volume of distribution V2 (L)") # Simulation equation: 'V2=2.88'
    lq_c_elf <- fixed(log(0.850)); label("Ceftolozane one-way central-to-ELF clearance Qce (L/h)") # Simulation equation: 'Qce=0.850*exp[N(0, 0.174)]'
    lq_elf_c <- fixed(log(1.66)); label("Ceftolozane one-way ELF-to-central clearance Qec (L/h)") # Simulation equation: 'Qec=1.66'
    lv_elf <- fixed(log(1)); label("ELF compartment volume Ve (L)") # Simulation equation: 'Ve=1'

    e_crcl_cl <- fixed(0.715); label("Power exponent of CrCL on ceftolozane CL, CRCL/109 (unitless)") # Simulation equation: '(CrCL/109)^0.715'
    e_wt_vc <- fixed(1); label("Power exponent of body weight on ceftolozane Vc, WT/74 (unitless)") # Simulation equation: '(WT/74)', no exponent
    e_ciai_cl <- fixed(0.195); label("cIAI-patient shift on log ceftolozane CL, applied to every simulated patient (log)") # Simulation equation: 'exp[0.195+...]' on CL; exp(0.195) = 1.22-fold (Chandorkar 2015 Table 3A)
    e_ciai_vc <- fixed(0.464); label("cIAI-patient shift on log ceftolozane Vc, applied to every simulated patient (log)") # Simulation equation: 'exp[0.464+...]' on Vc; exp(0.464) = 1.59-fold (Chandorkar 2015 Table 3A)

    etalcl ~ fixed(0.109) # Simulation equation: 'N(0, 109)' on CL, read as 0.109 = 0.330^2 (Chandorkar 2015 BSV 33.0%)
    etalvc ~ fixed(0.158) # Simulation equation: 'N(0, 0.158)' on Vc
    etalq_c_elf ~ fixed(0.174) # Simulation equation: 'N(0, 0.174)' on Qce

    propSd <- fixed(0); label("Proportional residual error, plasma (fraction; none specified for the simulations)") # not reported for the simulation model
    propSd_Celf <- fixed(0); label("Proportional residual error, ELF (fraction; none specified for the simulations)") # not reported for the simulation model
  })

  model({
    cl <- exp(lcl + etalcl) * exp(e_ciai_cl) * (CRCL / 109)^e_crcl_cl
    vc <- exp(lvc + etalvc) * exp(e_ciai_vc) * (WT / 74)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp)
    q_c_elf <- exp(lq_c_elf + etalq_c_elf)
    q_elf_c <- exp(lq_elf_c)
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
