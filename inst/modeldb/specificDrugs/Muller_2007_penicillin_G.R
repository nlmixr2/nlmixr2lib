Muller_2007_penicillin_G <- function() {
  description <- "Two-compartment IV bolus population PK model for penicillin G (benzylpenicillin) in 20 very preterm neonates with gestational age less than 32 weeks studied on day 3 of life (Muller 2007). Clearance is linearly scaled to current body weight with reference 1.195 kg (cohort mean); central volume, peripheral volume, and intercompartmental clearance are not weight-scaled in the final model."
  reference <- paste(
    "Muller AE, DeJongh J, Bult Y, Goessens WHF, Mouton JW, Danhof M, van den Anker JN.",
    "Pharmacokinetics of penicillin G in infants with a gestational age of less than 32 weeks.",
    "Antimicrob Agents Chemother. 2007 Oct;51(10):3720-5.",
    "doi:10.1128/AAC.00318-07.",
    sep = " "
  )
  vignette <- "Muller_2007_penicillin_G"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Current body weight on the PK sampling day (day 3 of life)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linearly scaled on CL with reference 1.195 kg (cohort mean per Muller 2007 Table 1). The paper retained body weight on CL in the final model (P < 0.01, Results page 3724 and Fig. 3) but did not report the functional form or coefficient; per operator sidecar-001 response Q2=B, a linear-with-weight scaling (exponent 1) at the cohort-mean reference is imputed. CL at reference subject = 0.103 L/h.",
      source_name        = "WT"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 20L,
    n_studies       = 1L,
    age_range       = "Day 3 of life (postnatal age); gestational age range 26 3/7 to 32 0/7 weeks at birth",
    age_median      = "Gestational age at birth: median 29 5/7 weeks (SD 1 5/7 weeks)",
    weight_range    = "650 to 2,030 g birth weight",
    weight_median   = "mean 1,195 g (SD 387 g)",
    sex_female_pct  = 40,
    race_ethnicity  = "Not reported (single-centre Dutch cohort; Erasmus MC-Sophia, Rotterdam)",
    disease_state   = "Preterm neonates with suspected or documented septicemia or invasive infection (no positive blood cultures in this cohort; one positive superficial culture for Streptococcus agalactiae). Hemodynamically stable, normal liver function, no nephrotoxic drugs, no intracranial hemorrhage beyond grade II.",
    dose_range      = "Penicillin G 50,000 U/kg as IV bolus every 12 h; ~30 mg/kg q12h using the conventional conversion 1 IU = 0.6 mg",
    regions         = "Netherlands (Erasmus MC-Sophia, Sophia Children's Hospital, Rotterdam)",
    gestational_age_range = "26 3/7 to 32 0/7 weeks at birth (Muller 2007 Table 1)",
    samples_plasma  = "167 samples; arterial-line draws pre-dose and at 0.03, 0.5, 1, 2.5, 4, 8, and 12 h after a dose, plus a 24 h sample in subjects who skipped the next dose (n = 9 of 20)",
    notes           = "Hematocrit median 46% (range 33-63), platelets median 203 x10^3/mm3 (range 74-497), creatinine median 46 (range 10-82, units not reported), Apgar 1 min median 6 (range 1-10), Apgar 5 min median 8 (range 6-10). 9 of 20 ventilated. Half of the cohort were born to mothers with preeclampsia or HELLP syndrome. Per Methods 'Estimation of fT>MIC', protein binding was estimated at 40% +/- 2.5% (Ebert 1988 ref. 14) but the paper notes this is likely an overestimate in neonates; protein binding is NOT used inside the structural ODEs (the model is on total drug). PK estimates were obtained with NONMEM v.V ADVAN5 (general linear), FOCE+I; combined additive + proportional residual error."
  )

  ini({
    # Structural parameters (Muller 2007 Table 2 'Structural model
    # parameters' column; mean +/- standard error). Reference subject
    # is the cohort-mean 1.195 kg neonate (Muller 2007 Table 1).
    lcl <- log(0.103); label("Clearance at reference weight (L/h at 1.195 kg)")        # Muller 2007 Table 2: CL = 0.103 L/h (SE 0.0104)
    lvc <- log(0.359); label("Central volume of distribution (L)")                     # Muller 2007 Table 2: V1 = 0.359 L (SE 0.0558)
    lvp <- log(0.152); label("Peripheral volume of distribution (L)")                  # Muller 2007 Table 2: V2 = 0.152 L (SE 0.0312)
    lq  <- log(0.774); label("Intercompartmental clearance (L/h)")                     # Muller 2007 Table 2: Q  = 0.774 L/h (SE 0.277)

    # Allometric / size-scaling exponent on CL. Operator sidecar-001
    # response Q2=B selected a linear-with-weight scaling (exponent 1)
    # at the cohort-mean reference because Muller 2007 retains body
    # weight on CL in the final model (P < 0.01, Results page 3724;
    # Fig. 3) but does not report the functional form or coefficient.
    e_wt_cl <- fixed(1.0); label("Exponent of (WT / 1.195 kg) on CL (unitless; imputed linear scaling)")  # Muller 2007 Results p. 3724 + Fig. 3 (form not reported; linear exponent imputed per operator sidecar-001 Q2=B)

    # Inter-individual variability (Muller 2007 Table 2 'Variance model
    # parameters'; reported as omega^2 of the log-normal eta).
    # The IIV on the second volume term is encoded on V1 (central) per
    # operator sidecar-001 response Q1=A: Table 2 explicitly labels the
    # row 'Interindividual variability in V1' with subscript 1 and the
    # Table-2 footnote defines V1 as the central compartment, while the
    # paper body text says V2 / peripheral. See vignette Errata for
    # the table/text discrepancy.
    etalcl ~ 0.164  # Muller 2007 Table 2: omega^2(CL) = 0.164 (SE 0.0865); text page 3723 reports 34.5% (interpreted as eta-shrinkage)
    etalvc ~ 0.39   # Muller 2007 Table 2: omega^2(V1) = 0.39  (SE 0.126);  text page 3723 reports 17.1% (interpreted as eta-shrinkage); see Errata for V1-vs-V2 discrepancy

    # Residual error (Muller 2007 Table 2 'Variance model parameters';
    # reported as variance of the combined additive + proportional
    # error model. Converted to SDs:
    #   proportional sigma^2 = 0.104 -> SD = sqrt(0.104) = 0.3225 (fraction)
    #   additive     sigma^2 = 1.12  -> SD = sqrt(1.12)  = 1.058  (mg/L)
    propSd <- 0.3225; label("Proportional residual error (fraction)")                 # Muller 2007 Table 2: sigma^2_prop = 0.104 (SE 0.0316); SD = sqrt(0.104)
    addSd  <- 1.058;  label("Additive residual error (mg/L)")                          # Muller 2007 Table 2: sigma^2_add  = 1.12  (SE 0.891);  SD = sqrt(1.12)
  })

  model({
    # Individual PK parameters. Reference subject WT = 1.195 kg
    # (cohort mean). Linear scaling on CL only; V1, V2, and Q
    # are not weight-scaled in the final model.
    cl <- exp(lcl + etalcl) * (WT / 1.195)^e_wt_cl
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q  <- exp(lq)

    # Micro-constants for the two-compartment ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Observation. Dose in mg, volumes in L -> concentrations in mg/L
    # (= ug/mL). Penicillin G dose conversion: 1 IU = 0.6 mg; the
    # paper's 50,000 U/kg q12h corresponds to ~30 mg/kg q12h.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
