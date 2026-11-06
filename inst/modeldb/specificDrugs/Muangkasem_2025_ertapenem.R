Muangkasem_2025_ertapenem <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenously infused",
    "ertapenem in critically ill adults, combined with a saturable",
    "(Michaelis-Menten) albumin-binding model that converts the total",
    "plasma concentration to the unbound concentration (Muangkasem 2025).",
    "The structural PK parameters and their interindividual variability",
    "are carried from the Burkhardt 2007 model in critically ill patients",
    "with ventilator-associated pneumonia; clearance scales as a power",
    "function of Cockcroft-Gault creatinine clearance, (CRCL / 78)^0.535,",
    "with the exponent fitted by Muangkasem 2025 to the Mistry 2006",
    "renal-impairment data. The binding capacity Bmax = 0.316 * ALB^1.611",
    "and the dissociation constant KD = 3.301 * exp(0.082 * ALB) (both",
    "mg/L, serum albumin in g/L) were fitted to the Liebchen 2014",
    "per-patient bound-versus-unbound ertapenem data from six critically",
    "ill patients. Built for Monte Carlo dosing simulation; no residual",
    "error is reported.",
    sep = " "
  )
  reference <- paste(
    "Muangkasem A, Nosoongnoen W, Sririttha S, Aramruang T, Montakantikul P",
    "(2025). Dosage Regimen Optimization of Ertapenem",
    "Against ESBL-Producing Enterobacterales Infection in Critically Ill",
    "Patients Using Monte Carlo Simulation. Clinical and Translational",
    "Science 18(7):e70274. doi:10.1111/cts.70274.",
    "Structural PK model from Burkhardt O, Kumar V, Katterwe D, et al.",
    "(2007). J Antimicrob Chemother 59(2):277-284; binding data from",
    "Liebchen U, et al. (2014). J Antimicrob Chemother 69(11):3108-3111;",
    "renal-function data from Mistry GC, et al. (2006). J Clin Pharmacol",
    "46(10):1128-1138.",
    sep = " "
  )
  vignette <- "Muangkasem_2025_ertapenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "ertapenem (total, bound + unbound)", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ertapenem", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation (raw, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Raw mL/min, NOT BSA-normalized (Methods: 'CrCl is the creatinine",
        "clearance calculated by the Cockcroft-Gault equation'). Enters",
        "clearance as the power function (CRCL / 78)^0.535 (Equation 1 with",
        "the Results power coefficient 0.535). The 78 mL/min reference is",
        "printed in Equation 1 without saying which cohort statistic it is.",
        "The 0.535 exponent was fitted by least squares to clearance versus",
        "creatinine clearance digitised from Mistry 2006 (r^2 = 0.81,",
        "Figure S1). The source simulations span 15-200 mL/min."
      ),
      source_name = "CrCl"
    ),
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Not a covariate on any PK parameter. Sets the albumin-binding",
        "capacity Bmax = 0.316 * ALB^1.611 (power fit, r^2 = 0.99,",
        "Figure S3) and the dissociation constant KD = 3.301 *",
        "exp(0.082 * ALB) (exponential fit, r^2 = 1, Figure S4), both in",
        "mg/L, which convert total to unbound ertapenem. The binding data",
        "came from six critically ill patients with albumin 24-30 g/L",
        "(Discussion, limitations); the source simulations extrapolate to",
        "20-40 g/L."
      ),
      source_name = "serum albumin"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 3L,
    disease_state = paste(
      "Critically ill adults. The model is a composite assembled from three",
      "published data sets rather than fitted to one cohort: the Burkhardt",
      "2007 two-compartment model in critically ill patients with",
      "early-onset ventilator-associated pneumonia (structural parameters",
      "and IIV), the Mistry 2006 study in patients with varying degrees of",
      "renal insufficiency (clearance-creatinine clearance exponent), and",
      "the Liebchen 2014 unbound-fraction data in six ICU patients",
      "(albumin-binding model)."
    ),
    dose_range = paste(
      "Simulated: 0.5 g or 1 g every 12 or 24 h and 2 g every 24 h, each as",
      "a 0.5 h or 3 h IV infusion (10 regimens)."
    ),
    renal_function = "Simulated Cockcroft-Gault creatinine clearance 15-200 mL/min.",
    albumin_range = "Simulated serum albumin 20-40 g/L; binding data observed at 24-30 g/L.",
    regions = "Authors in Thailand; the source cohorts' locations are not stated in the paper.",
    notes = paste(
      "The paper does not tabulate the demographics of the source cohorts;",
      "it reports only the ranges of the simulated population (Methods,",
      "Dosage Regimen Simulations). n_subjects is left NA because the model",
      "pools parameters from separate studies."
    )
  )

  ini({
    # Structural PK. Results 'Clearance and Creatinine Clearance': 'The final
    # population pharmacokinetic parameter estimates and inter-individual
    # variability in parenthesis were as follows: 2.63 L/h (31.4%) for
    # clearance, 10.6 L (23.7%) for central volume of distribution, 4.24 L
    # (50.6%) for peripheral volume of distribution, and 12.32 L/h (0%) for
    # intercompartment clearance.' These are carried from Burkhardt 2007
    # (Methods: 'A previous ertapenem population pharmacokinetic model in
    # critically ill patients was developed by Burkhardt et al.') and were not
    # re-estimated, so they are fixed.
    lcl <- fixed(log(2.63)); label("Clearance at CRCL = 78 mL/min (L/h)") # Results; Equation 1
    lvc <- fixed(log(10.6)); label("Central volume of distribution (L)") # Results
    lvp <- fixed(log(4.24)); label("Peripheral volume of distribution (L)") # Results
    lq <- fixed(log(12.32)); label("Intercompartmental clearance (L/h)") # Results

    # Covariate effect of creatinine clearance on CL, estimated in this paper
    # by least-squares power regression on the Mistry 2006 data.
    e_crcl_cl <- 0.535; label("Power exponent of CRCL on CL (unitless)") # Results: 'A power coefficient of 0.535 was incorporated into Equation 1'

    # Albumin-binding model (Equation 2, Cbound = Bmax * Cfree / (KD + Cfree)).
    # Results 'Bmax and Serum Albumin Concentration': power law with
    # 'constant = 0.316, power coefficient = 1.611'; 'KD and Serum Albumin
    # Concentration': exponential with 'constant = 3.301, exponential
    # coefficient = 0.082'. ALB in g/L, Bmax and KD in mg/L (reproduces the
    # reported median Bmax 75.74 and KD 38.64 mg/L at ALB 30 g/L).
    lbmax <- log(0.316); label("Albumin-binding capacity Bmax at ALB = 1 g/L (mg/L)") # Results: Bmax power-law constant 0.316
    e_alb_bmax <- 1.611; label("Power exponent of ALB on Bmax (unitless)") # Results: Bmax power coefficient 1.611
    lkd <- log(3.301); label("Equilibrium dissociation constant KD at ALB = 0 g/L (mg/L)") # Results: KD exponential constant 3.301
    e_alb_kd <- 0.082; label("Exponential coefficient of ALB on KD (per g/L)") # Results: KD exponential coefficient 0.082

    # IIV: exponential (log-normal) random effects, omega^2 = log(1 + CV^2).
    # No eta on Q (reported IIV 0%).
    etalcl ~ fixed(0.0940330) # Results: CL IIV 31.4% -> log(1 + 0.314^2)
    etalvc ~ fixed(0.0546482) # Results: Vc IIV 23.7% -> log(1 + 0.237^2)
    etalvp ~ fixed(0.2279610) # Results: Vp IIV 50.6% -> log(1 + 0.506^2)

    # Residual error: not reported (simulation-only paper). Zero rather than
    # invented; see the vignette Assumptions and deviations.
    propSd <- fixed(0); label("Proportional residual error on total concentration (fraction; ZERO - not reported in source)") # not reported
  })
  model({
    cl <- exp(lcl + etalcl) * (CRCL / 78)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Dose is an IV infusion into central (rate or duration on the dose record).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Total plasma concentration (mg/L).
    Cc <- central / vc

    # Unbound concentration from the albumin-binding model. Solving
    # Cc = Cunbound + bmax * Cunbound / (kd + Cunbound) for Cunbound gives
    # Cunbound^2 + b * Cunbound - kd * Cc = 0 with b = kd + bmax - Cc; the
    # positive root is written in the rationalised form, which stays finite
    # and non-negative for every Cc >= 0.
    bmax <- exp(lbmax) * ALB^e_alb_bmax
    kd <- exp(lkd) * exp(e_alb_kd * ALB)
    bterm <- kd + bmax - Cc
    Cunbound <- 2 * kd * Cc / (bterm + sqrt(bterm * bterm + 4 * kd * Cc))

    Cc ~ prop(propSd)
  })
}
