Caldes_2009_ganciclovir <- function() {
  description <- "Two-compartment population PK model for ganciclovir after IV ganciclovir and oral valganciclovir administration in solid organ transplant patients infected with cytomegalovirus, with first-order absorption, lag time, logit-transformed bioavailability, and creatinine-clearance scaling on CL (Caldes 2009)"
  reference   <- "Caldes A, Colom H, Armendariz Y, Garrido MJ, Troconiz IF, Gil-Vernet S, Lloberas N, Pou L, Peraire C, Grinyo JM. Population pharmacokinetics of ganciclovir after intravenous ganciclovir and oral valganciclovir administration in solid organ transplant patients infected with cytomegalovirus. Antimicrob Agents Chemother. 2009;53(11):4816-4824. doi:10.1128/AAC.00085-09"
  vignette    <- "Caldes_2009_ganciclovir"
  units       <- list(time = "hr", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance (Cockcroft-Gault formula, raw mL/min, NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-scaling effect on CL with reference 57 mL/min (mean population CLCR). Caldes 2009 reports CL = 7.49 * (CLCR / 57) L/h (Table 3); covariate enters as a linear ratio (power exponent fixed at 1), not as a power-law exponent. Cockcroft-Gault was computed with the manufacturer's renal dose-adjustment table (Table 1); per the Methods, CrCl was used in raw mL/min without BSA normalization. The canonical CRCL register's units are mL/min/1.73 m^2, but per Henin 2009 / Goel 2016 precedent the canonical name is reused for raw-CG-CrCl with the unit difference documented here.",
      source_name        = "CLCR"
    )
  )

  population <- list(
    n_subjects     = 20L,
    n_studies      = 1L,
    age_range      = "mean 55.7 (SD 11.8) years; per-transplant subgroup means 54.6 (kidney), 53.2 (liver), 54.4 (heart)",
    weight_range   = "mean 66.2 (SD 12.9) kg; per-transplant subgroup means 65.9 (kidney), 69.6 (liver), 65.8 (heart)",
    sex_female_pct = 50,
    race_ethnicity = "All 20 patients reported as Caucasian (Results, Patients section).",
    disease_state  = "Solid organ transplant recipients (10 kidney, 5 liver, 5 heart) with established CMV infection (positive CMV antigenemia pp65 >= 20 positive cells per 10^5 peripheral blood mononuclear cells); excluded if positive CMV tissue-invasive disease, ANC < 500/mm^3, platelets < 25,000/mm^3, hemoglobin < 80 g/L, or estimated GFR < 10 mL/min.",
    dose_range     = "IV ganciclovir 5 mg/kg BID for 5 days, then oral valganciclovir 900 mg BID for 16 days; doses adjusted for renal function per Table 1 (ganciclovir 5/2.5/2.5 mg/kg at q12/q24/q48 for CLCR >50 / 25-50 / <25 mL/min; valganciclovir 900/450/450/450 mg at q12/q12/q24/q48 for CLCR bands >=60 / 40-59 / 25-39 / 10-24 mL/min).",
    regions        = "Single center, Barcelona, Spain (Hospital Universitari de Bellvitge); enrollment March 2004-February 2006.",
    crcl_range     = "mean 57.0 (SD 25.3) mL/min (Cockcroft-Gault); per-transplant subgroup means 39.9 (kidney), 75.1 (liver), 73.2 (heart).",
    co_medications = "Concomitant immunosuppression: mycophenolate mofetil 16/20, cyclosporine 11/20, tacrolimus 8/20, sirolimus 1/20.",
    trial_id       = "ClinicalTrials.gov NCT00730769",
    n_observations = 382L,
    notes          = "21 patients enrolled; one liver-transplant patient excluded from the population PK analysis due to pancytopenia (treatment not completed). Final dataset: 382 ganciclovir serum concentrations (190 IV + 192 oral, 6 below LOQ of 5 ug/L). Sampling on day 5 (last IV dose) and day 15 (oral): 0, 0.5, 1, 1.5, 2, 3, 4, 8, 12 h post-dose, with extended sampling to 24 h for patients with CLCR < 10 mL/min (IV) or 16 mL/min (oral). Demographics from Table 2."
  )

  ini({
    # Structural parameters at the reference covariate set (CRCL = 57 mL/min).
    # All structural PK parameters in linear units (L, L/h, 1/h, h); concentrations
    # in ug/mL = mg/L given mg-doses and L-volumes.
    lcl    <- log(7.49);  label("Clearance at reference CRCL = 57 mL/min (CL, L/h)")    # Caldes 2009 Table 3 final CL
    lvc    <- log(31.90); label("Central volume of distribution (V1, L)")               # Caldes 2009 Table 3 final V1
    lq     <- log(10.20); label("Inter-compartmental clearance (CL_D, L/h)")            # Caldes 2009 Table 3 final CL_D
    lvp    <- log(32.0);  label("Peripheral volume of distribution (V2, L)")            # Caldes 2009 Table 3 final V2
    lka    <- log(0.895); label("First-order oral absorption rate constant (Ka, 1/h)")  # Caldes 2009 Table 3 final Ka
    ltlag  <- log(0.382); label("Absorption lag time (Tlag, h)")                        # Caldes 2009 Table 3 final Lag time

    # Bioavailability is logit-transformed in NONMEM to keep F in [0,1].
    # logit(0.825) = log(0.825 / (1 - 0.825)) = 1.5505978...
    logitfdepot <- logit(0.825); label("Logit of oral ganciclovir-equivalent bioavailability (F = 0.825)")  # Caldes 2009 Table 3 final F

    # Linear CRCL effect on CL: CL = 7.49 * (CRCL / 57). Stored as a power exponent
    # of 1 on the (CRCL / 57) ratio so the covariate-effect parameter is identifiable
    # by checkModelConventions() and the form mirrors other CL ~ renal-function models.
    e_crcl_cl <- 1.0; label("Power exponent of CRCL on CL (unitless; reference 57 mL/min, fixed at 1 per linear paper form)")  # Caldes 2009 Table 3 footnote c (linear scaling: CL = 7.49 * (CLCR/57))

    # Inter-individual variability. Caldes 2009 Table 3 reports the IIV variances
    # directly (omega^2 column, expressed as variance). NONMEM exponential-on-log-scale
    # parametrisation: CL_i = CL_pop * exp(eta_CL); F_i is on the logit scale with
    # logit(F_i) = logit(F_pop) + eta_F.
    etalcl         ~ 0.107  # Caldes 2009 Table 3 final IIV CL = 0.107 (32.71% CV)
    etalvc         ~ 0.227  # Caldes 2009 Table 3 final IIV V1 = 0.227 (53.00% CV)
    etalka         ~ 0.464  # Caldes 2009 Table 3 final IIV Ka = 0.464 (36.85% CV; SE on CV scale)
    etalogitfdepot ~ 0.049  # Caldes 2009 Table 3 final IIV F = 0.049 (logit-scale variance; back-transformed CV reported as 63.90%)

    # Residual error: combined additive + proportional (Caldes 2009 Methods, "RE was
    # described by a combined error model"). Table 3 reports sigma_1 as the additive
    # SD (ug/mL) and sigma_2^2 as the proportional variance.
    propSd <- sqrt(0.143); label("Proportional residual error SD (fraction; sqrt of variance 0.143 reported in Table 3)")  # Caldes 2009 Table 3 final sigma_2^2
    addSd  <- 0.465;        label("Additive residual error SD (ug/mL)")                                                     # Caldes 2009 Table 3 final sigma_1
  })

  model({
    # Renal-function effect on CL (linear scaling per Caldes 2009 footnote c).
    cl <- exp(lcl + etalcl) * (CRCL / 57)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q  <- exp(lq)
    ka <- exp(lka + etalka)

    tlag_depot <- exp(ltlag)
    fdepot     <- expit(logitfdepot + etalogitfdepot)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)    <- fdepot
    alag(depot) <- tlag_depot

    # Concentration: dose in mg, volume in L -> central / vc in mg/L = ug/mL.
    # IV doses target the central compartment; oral valganciclovir doses target
    # depot and must be supplied as ganciclovir-equivalent (multiply the
    # valganciclovir dose by 0.720, the ganciclovir/valganciclovir mass ratio,
    # per Caldes 2009 Methods "Pharmacokinetic data analysis").
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
