Cella_2012_midazolam_children_adolescents <- function() {
  description <- paste(
    "Two-compartment population PK model for midazolam in children and",
    "adolescents (Cella 2012 Model 2), IV bolus only, with per-kg linear",
    "scaling of the central volume and a linear normalisation of the",
    "peripheral volume by age (months) at a 74-month reference. Cohort of",
    "18 paediatric oncology patients (ages 3.2 to 16.2 years, body weights",
    "12.6 to 60.1 kg) dosed at 0.12 mg/kg IV before invasive procedures.",
    "Fitted with informative priors from De Wildt 2002 via the NONMEM PRIOR /",
    "Wishart subroutine."
  )
  reference <- paste(
    "Cella M, Knibbe C, de Wildt SN, Van Gerven J, Danhof M, Della Pasqua O",
    "(2012). Scaling of pharmacokinetics across paediatric populations: the",
    "lack of interpolative power of allometric models. Br J Clin Pharmacol",
    "74(3):525-535. doi:10.1111/j.1365-2125.2012.04206.x.",
    sep = " "
  )
  vignette <- "Cella_2012_midazolam_paediatric_scaling"
  units    <- list(time = "minute", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. Used for per-kg linear scaling of the",
        "central volume (Vc = 1.95 * WT). The per-kg interpretation follows",
        "operator sidecar response 001 Q2 = A; the published Table 2 unit",
        "label 'l' (absolute) is kinetically implausible for a typical",
        "29 kg child (k10 ~ 0.097/min, t1/2 ~ 7 min). See the vignette",
        "Errata for the reasoning."
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. The model uses age in months (AGE * 12)",
        "normalised to a 74-month reference for the linear scaling of the",
        "peripheral volume per Cella 2012 Table 2 ('Vp (l x months/74)' =",
        "7.14)."
      ),
      source_name        = "AGE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 18L,
    n_studies      = 1L,
    age_range      = "3.2 - 16.2 years",
    age_median     = "7.7 years (mean); ~6.17 years (74 months) used as Vp reference",
    weight_range   = "12.6 - 60.1 kg",
    weight_median  = "29.0 kg (mean)",
    sex_female_pct = NA_real_,
    disease_state  = paste(
      "Paediatric oncology patients undergoing invasive procedures. Sparse",
      "PK sampling (mean 4.6 samples per subject)."
    ),
    dose_range     = "0.12 mg/kg IV bolus (mean dose; per-subject doses vary slightly)",
    regions        = "Collaborative study between Purdue University (USA) and Sophia Children's Hospital (Netherlands)",
    notes          = paste(
      "Per Cella 2012 Table 1. Informative priors from De Wildt (CL = 5",
      "mL/kg/min, Vc = 0.38 L/kg, Vp = 1.7 L/kg) were used to stabilise the",
      "model under sparse sampling per Cella 2012 Methods; the paper notes",
      "the priors did not appreciably affect the final parameter estimates."
    )
  )

  ini({
    # Structural parameters from Table 2 ('Mean' column for Model 2
    # 'Children and adolescents'). Vc per-kg encoding (sidecar response
    # 001 Q2 = A) gives physically plausible kinetics (k10 ~ 0.0034/min,
    # t1/2 ~ 206 min for a typical 29 kg child); the literal absolute
    # reading gives t1/2 ~ 7 min, ~25x too fast. Vp uses the published
    # 'Vp (l x months/74)' linear normalisation, with the 74 months
    # representing the population-median age in months.
    lcl <- log(0.19);   label("Clearance (CL, L/min)")                                              # Cella 2012 Table 2 (CL = 0.19)
    lvc <- log(1.95);   label("Central volume of distribution per kg body weight (Vc/WT, L/kg)")    # Cella 2012 Table 2 (Vc = 1.95; per-kg interpretation per sidecar response 001 Q2)
    lq  <- log(0.105);  label("Inter-compartmental clearance (Q, L/min)")                           # Cella 2012 Table 2
    lvp <- log(7.14);   label("Peripheral volume of distribution at AGE = 74 months reference (Vp, L)") # Cella 2012 Table 2 (Vp = 7.14 with linear age scaling)

    # Inter-individual variability (Table 2 reports CV%; omega^2 = log(CV^2 + 1))
    #   CL : 32.7% CV -> log(0.327^2 + 1) = 0.10157
    #   Vc : 31.5% CV -> log(0.315^2 + 1) = 0.09456
    etalcl ~ 0.10157                                                                                # Cella 2012 Table 2 (IIV CL 32.7% CV)
    etalvc ~ 0.09456                                                                                # Cella 2012 Table 2 (IIV Vc 31.5% CV)

    # Residual error (39.0% CV proportional, linear-scale SD = 0.39)
    propSd <- 0.39;     label("Proportional residual error (SD, fraction)")                         # Cella 2012 Table 2 (residual 39.0% CV, proportional model)
  })
  model({
    # AGE-in-months term for the Vp linear-age scaling.
    age_months <- AGE * 12

    # Individual PK parameters.
    # CL: constant population-typical value (no WT or age scaling in Model 2).
    # Vc: per-kg linear scaling (Vc/WT is the estimated per-kg coefficient).
    # Q : constant population-typical value.
    # Vp: linear normalisation by age (months) at a 74-month reference.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * WT
    q  <- exp(lq)
    vp <- exp(lvp) * (age_months / 74)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # IV bolus only (no oral arm in Model 2). Doses go directly to central.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Plasma midazolam concentration. With dose in mg and Vc in L, Cc is in
    # mg/L (numerically equivalent to ug/mL and to ng/mL * 1e-3).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
