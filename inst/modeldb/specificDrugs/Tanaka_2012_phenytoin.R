Tanaka_2012_phenytoin <- function() {
  description <- "Two-compartment population PK model for phenytoin after IV fosphenytoin sodium administration in Japanese healthy volunteers and adult / pediatric patients (Tanaka 2012). The fosphenytoin compartment converts first-order (K12) to the phenytoin central compartment; phenytoin is cleared from central and exchanges with a peripheral compartment via Q."
  reference <- "Tanaka J, Kasai H, Shimizu K, Shimasaki S, Kumagai Y. Population pharmacokinetics of phenytoin after intravenous administration of fosphenytoin sodium in pediatric patients, adult patients, and healthy volunteers. Eur J Clin Pharmacol. 2013;69(3):489-497. doi:10.1007/s00228-012-1373-8"
  vignette <- "Tanaka_2012_phenytoin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as a power covariate on CL, V2 (central), and V3 (peripheral). Reference 60 kg (Tanaka 2012 Results: 'The average BW of adult Japanese men (60 kg) was selected as the standard value'). The WT exponent on V2 was fixed at 1.0 (Results: 'The influence factor of V2 was fixed to 1 on the basis of statistical significance'). Paper uses the alias 'BW' interchangeably.",
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 71L,
    n_studies      = 3L,
    age_range      = "2-86 years",
    weight_range   = "7.8-74.4 kg",
    sex_female_pct = 38.0,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Pooled cohort of 24 healthy adult volunteers (Phase I crossover + Phase I dose-escalation), 14 adult patients (Phase III; status epilepticus, acute repetitive seizures, or post-neurosurgical seizure prophylaxis), and 33 pediatric patients (Phase III; same indications).",
    dose_range     = "IV fosphenytoin sodium 375-750 mg (healthy volunteers, fixed doses) or 15-22.5 mg/kg (patients) infused at 8.3-75 mg/min (healthy volunteers) or 1-3 mg/kg/min capped at 150 mg/min (patients). Phase I crossover arm also dosed phenytoin sodium 250 mg IV directly.",
    regions        = "Japan",
    notes          = "Tanaka 2012 Table 2 baseline demographics. Three pooled studies: two Phase I in healthy adult Japanese males (n = 24 total, all male) and one Phase III in neurosurgical and epileptic patients (adult n = 14, pediatric n = 33; combined adult+pediatric sex 20/27 M/F). Overall sex 44/27 M/F. 923 phenytoin plasma concentrations. Reference body weight 60 kg (average Japanese adult male)."
  )

  ini({
    # Structural parameters - Tanaka 2012 Table 3 (final population PK model);
    # typical values reported for a 60-kg adult Japanese male.
    lcl <- log(1.61);  label("Phenytoin clearance (L/h)")                          # Table 3 theta_1: CL = 1.61 L/h
    lvc <- log(20.8);  label("Phenytoin central volume of distribution V2 (L)")    # Table 3 theta_3: V2 = 20.8 L
    lq  <- log(53.0);  label("Phenytoin intercompartmental clearance Q (L/h)")     # Table 3 theta_4: Q  = 53.0 L/h
    lvp <- log(26.0);  label("Phenytoin peripheral volume of distribution V3 (L)") # Table 3 theta_5: V3 = 26.0 L
    lka <- log(5.02);  label("Fosphenytoin to phenytoin conversion rate K12 (1/h)") # Table 3 theta_7: K12 = 5.02 1/h

    # Allometric exponents on body weight (Tanaka 2012 Equation 1, p494).
    # Reference weight 60 kg (average Japanese adult male; Results p493-494).
    e_wt_cl <- 0.569;      label("Body-weight exponent on CL (unitless)")  # Table 3 theta_2: 0.569
    e_wt_vc <- fixed(1);   label("Body-weight exponent on V2, fixed (unitless)")  # Results p494: "The influence factor of V2 was fixed to 1"; no theta entry in Table 3
    e_wt_vp <- 0.584;      label("Body-weight exponent on V3 (unitless)")  # Table 3 theta_6: 0.584

    # Inter-subject variability (variances reported on omega^2 scale in Table 3).
    etalcl ~ 0.194    # Table 3 omega_{CL,CL}  = 0.194 (variance)
    etalvc ~ 0.161    # Table 3 omega_{V2,V2}  = 0.161
    etalq  ~ 0.271    # Table 3 omega_{Q,Q}    = 0.271
    etalvp ~ 0.0430   # Table 3 omega_{V3,V3}  = 0.0430
    etalka ~ 0.106    # Table 3 omega_{K12,K12} = 0.106

    # Combined exponential + additive residual error.
    # Tanaka 2012 p491: Y_ij = F_ij * exp(eps1_ij) + eps2_ij, sigma_1^2 and sigma_2^2
    # are reported on the variance scale (Table 3). Converting variance -> SD for
    # nlmixr2: propSd = sqrt(sigma_1^2); addSd = sqrt(sigma_2^2). The exponential
    # form linearizes to (1 + eps) for small eps, mapping to nlmixr2's prop() arm.
    propSd <- sqrt(0.00148);  label("Proportional residual SD (fraction)")  # Table 3 sigma_{1,1} = 0.00148 (variance)
    addSd  <- sqrt(0.317);    label("Additive residual SD (ug/mL)")         # Table 3 sigma_{2,2} = 0.317  (variance)
  })

  model({
    # Individual PK parameters with allometric weight scaling (reference 60 kg)
    cl  <- exp(lcl + etalcl) * (WT / 60)^e_wt_cl
    vc  <- exp(lvc + etalvc) * (WT / 60)^e_wt_vc
    q   <- exp(lq  + etalq)
    vp  <- exp(lvp + etalvp) * (WT / 60)^e_wt_vp
    ka  <- exp(lka + etalka)

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system (Tanaka 2012 Figure 1):
    #   depot       = fosphenytoin (IV dosing compartment for fosphenytoin sodium)
    #   central     = phenytoin central (V2); converted from fosphenytoin at rate K12 (= ka here)
    #   peripheral1 = phenytoin peripheral (V3)
    # IV fosphenytoin doses enter depot; if direct IV phenytoin is dosed (Phase I
    # crossover arm only), the dataset cmt column targets central.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Observation: total plasma phenytoin concentration (ug/mL = mg/L)
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
