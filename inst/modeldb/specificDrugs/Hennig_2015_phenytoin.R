Hennig_2015_phenytoin <- function() {
  description <- "One-compartment population PK model for phenytoin in critically ill children with a linear partition coefficient describing protein binding to albumin (Hennig 2015)."
  reference <- "Hennig S, Norris R, Tu Q, van Breda K, Riney K, Foster K, Lister B, Charles B. Population Pharmacokinetics of Phenytoin in Critically Ill Children. J Clin Pharmacol. 2015;55(3):355-364. doi:10.1002/jcph.417"
  vignette <- "Hennig_2015_phenytoin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Allometric scaling on CL (exponent 0.75) and on V2 / V3 (exponent 1) referenced to 70 kg per Hennig 2015 Methods (p357) and final model equations (p360).",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Reference 35 g/L. Linear effect on the unbound-bound partition coefficient PUB per Hennig 2015 Table 2 and the equation on p360: PUB = 8.23 * (1 + 0.00737 * (ALB - 35)).",
      source_name        = "ALB"
    )
  )

  population <- list(
    n_subjects     = 32L,
    n_studies      = 1L,
    age_range      = "0.08-17 years (1 month to 17 years)",
    age_median     = "mean 6.9 years (SD 5.9)",
    weight_range   = "4.0-80.0 kg",
    weight_median  = "mean 27.9 kg (SD 21.2)",
    sex_female_pct = 34.4,
    race_ethnicity = "Not reported",
    disease_state  = "Critically ill children admitted to a paediatric intensive care unit and treated with phenytoin for prevention or control of seizures (including status epilepticus). Underlying etiologies included acute brain injury, hypoxic ischemic encephalopathy, traumatic and infective brain injury, brain injury from motor vehicle accidents, brain tumour, haemorrhage (e.g., cerebral aneurysm), and longstanding structural abnormalities (e.g., cerebral palsy).",
    dose_range     = "IV 1.7-22.0 mg/kg per dose (mean 5.4); PO 0.7-10.0 mg/kg per dose (mean 3.3). 19.5 percent of doses IV, 80.5 percent PO (Dilantin 30 mg/5 mL pediatric oral suspension administered mainly via nasogastric or transpyloric tube; IV via 100 mg/2 mL ampoules). 8 children received both routes at different times.",
    regions        = "Single-centre paediatric intensive care unit at Mater Children's Hospital, Brisbane, Queensland, Australia (November 2006 to October 2009).",
    notes          = "Hennig 2015 Table 1 baseline demographics. 11 patients aged <2 years, 11 patients 2-11 years, 10 patients >11 years (boys/girls 21/11). Albumin range 11.0-48.0 g/L (mean 35.0, SD 9.3); 70 percent of patients had a change in serum albumin during the study with mean (range) absolute change of 26 percent (2-54). Total of 292 paired phenytoin concentrations: 146 protein-unbound + 146 protein-bound. Fosphenytoin not used in the study (not registered in Australia at the time)."
  )

  ini({
    # Structural PK parameters - Hennig 2015 Table 2 final-model column.
    # Reference body weight 70 kg; concentration in mg/L, time in hours.
    lcl  <- log(14.0); label("Clearance for unbound phenytoin at 70 kg (CL, L/h)")        # Table 2: CL = 14.0 L/h/70 kg
    lvc  <- log(447);  label("Volume of distribution for unbound phenytoin at 70 kg (V2, L)")  # Table 2: V2 = 447 L/70 kg
    lvp  <- fixed(log(2.8));        label("Volume of distribution for bound phenytoin at 70 kg, fixed to plasma albumin volume (V3, L)")  # Table 2: V3 = 2.8 L/70 kg, fixed
    lka  <- fixed(log(0.225));      label("First-order oral absorption rate constant, fixed (Ka, 1/h)")  # Methods p359 / Table 2: ka = 0.225 1/h, fixed after sensitivity analysis
    lpub <- log(8.23); label("Unbound-bound partition coefficient at reference albumin 35 g/L (PUB, dimensionless)")  # Table 2: PUB = 8.23 at ALB = 35 g/L
    lteq <- fixed(log(1.1e-4));     label("Half-life of unbound-bound equilibration, fixed to 0.4 s (Teq, h)")  # Methods p358: Teq fixed to 0.4 s = 1.1e-4 h
    logitfdepot <- log(0.63 / (1 - 0.63));  label("Logit-transformed oral bioavailability (logit-scale; F = 0.63)")  # Table 2: F1 = 63 percent (logit-transformed in NONMEM run186 to keep F in [0,1])

    # Covariate effects
    e_wt_cl   <- fixed(0.75); label("Allometric exponent on CL, fixed (unitless)")  # Methods p357 and p360 final equation: CL_i = CL * (WT/70)^0.75
    e_wt_vc   <- fixed(1);    label("Allometric exponent on V2, fixed (unitless)")  # Methods p360: V2_i = V2 * (WT/70)
    e_wt_vp   <- fixed(1);    label("Allometric exponent on V3, fixed (unitless)")  # Methods p360: V3_i = 2.8 * (WT/70)
    e_alb_pub <- 0.00737;     label("Linear slope of PUB on (ALB - 35), per g/L")   # Table 2 and equation p360: PUB = 8.23 * (1 + 0.00737 * (ALB - 35))

    # IIV. Hennig 2015 reports CV percent equivalent to 100 * sqrt(omega^2) for
    # log-normal etas; values below are omega^2 derived from Table 2 final-model
    # column (squared fractional CV).
    etalcl         ~ 0.2275  # Table 2 IIV CL = 47.7 percent CV
    etalvc         ~ 0.7174  # Table 2 IIV V2 = 84.7 percent CV
    etalka         ~ 3.27    # Table 2 IIV ka = 180.8 percent CV (typical ka fixed; IIV estimated)
    etalpub        ~ 0.0098  # Table 2 IIV PUB = 9.9 percent CV
    etalogitfdepot ~ 3.30    # NONMEM run186 OMEGA on the logit-F eta; paper Table 2 reports the back-transformed effective IIV F = 87.7 percent CV via the propagation IIV_F (CV) ~ sqrt(F * (1 - F)) * omega_F * 100

    # Residual error. The original NONMEM model used two independent proportional
    # EPS terms per output (analytical-assay error + a shared common-proportional
    # error EPS3). In nlmixr2 these collapse into a single proportional SD per
    # output: total prop SD = sqrt(assay_SD^2 + common_SD^2).
    # Cu (unbound): sqrt(0.136^2 + 0.174^2) = 0.221
    # Cb (bound):   sqrt(0.103^2 + 0.174^2) = 0.202
    # The within-pair correlation between Cu and Cb residuals induced by the
    # shared EPS3 in NONMEM is not preserved; marginal proportional SD per
    # output is.
    propSd_Cu <- 0.221; label("Proportional residual SD for unbound phenytoin (Cu, fraction)")  # Table 2: PHYu assay error 13.6 percent combined with common prop 17.4 percent
    propSd_Cb <- 0.202; label("Proportional residual SD for bound phenytoin (Cb, fraction)")    # Table 2: PHYb assay error 10.3 percent combined with common prop 17.4 percent
  })

  model({
    # Individual PK parameters with allometric weight scaling (reference 70 kg)
    cl  <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc  <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp  <- exp(lvp)          * (WT / 70)^e_wt_vp
    ka  <- exp(lka + etalka)
    teq <- exp(lteq)
    keq <- log(2) / teq
    q   <- vp * keq

    # Albumin effect on the partition coefficient (linear about reference 35 g/L)
    pub <- exp(lpub + etalpub) * (1 + e_alb_pub * (ALB - 35))

    # Logit-transformed bioavailability (kept in [0, 1] regardless of eta value).
    # Split the mu-reference onto its own line so nlmixr2 recognises etalogitfdepot
    # as mu-referenced.
    phi_fdepot <- logitfdepot + etalogitfdepot
    fdepot <- 1 / (1 + exp(-phi_fdepot))

    # Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system. depot = oral gut compartment; central = unbound phenytoin;
    # peripheral1 = compartment whose mass / V3 multiplied by PUB equals the
    # measured plasma bound phenytoin concentration (Hennig 2015 Figure 1).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Bioavailability applies to the depot (oral); IV doses are administered
    # directly to central via the dataset cmt column and bypass fdepot.
    f(depot) <- fdepot

    # Observations (Hennig 2015 Figure 1 caption):
    #   Cu = A_central / V2
    #   Cb = (A_peripheral1 / V3) * PUB
    Cu <- central / vc
    Cb <- (peripheral1 / vp) * pub

    Cu ~ prop(propSd_Cu)
    Cb ~ prop(propSd_Cb)
  })
}
