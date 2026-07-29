Cella_2012_midazolam_infants_adults <- function() {
  description <- paste(
    "Two-compartment population PK model for midazolam in infants, toddlers,",
    "and adults (Cella 2012 Model 1), with first-order absorption supporting",
    "intravenous and oral dosing, body-weight allometric scaling of clearance",
    "(exponent fixed to 0.75 at a 70 kg reference), per-kg linear scaling of",
    "the central volume, and a constant peripheral volume. Pooled cohort of",
    "23 infants and toddlers in a paediatric surgical ICU and 34 healthy",
    "adult volunteers."
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
        "Time-fixed at baseline. Enters the model with two distinct",
        "functional forms: allometric (exponent 0.75) on clearance with a",
        "70 kg reference, and linear per-kg scaling on the central volume",
        "(Vc = 0.312 * WT). The per-kg Vc interpretation follows operator",
        "sidecar response 001 Q2 = A; the published Table 2 unit label",
        "'l' is kinetically implausible if taken as an absolute typical",
        "value (k10 ~ 0.75/min for a 70 kg adult). See the vignette",
        "Errata for the full reasoning."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 57L,
    n_studies      = 4L,
    age_range      = "3.2 months - 29.7 years (infants/toddlers 3.2-24.7 months; adults 19.9-29.7 years)",
    age_median     = "Mixed: 10.9 months in the infant/toddler arm; 23.8 years in the adult arm",
    weight_range   = "5.1 - 91 kg (infants/toddlers 5.1-12 kg; adults 59-91 kg)",
    weight_median  = "Mixed: 9.2 kg infants/toddlers, 72.3 kg adults; 70 kg used as the allometric CL reference (see vignette Errata)",
    sex_female_pct = NA_real_,
    disease_state  = paste(
      "Healthy adult volunteers (n=34) from three Centre for Human Drug",
      "Research crossover trials (89110-pilot, 89110, 94113) plus infants",
      "and toddlers (n=23) admitted to a paediatric surgical ICU at",
      "Erasmus MC - Sophia's Children's Hospital after elective",
      "craniofacial surgery."
    ),
    dose_range     = paste(
      "Infants/toddlers: 0.1 mg/kg IV bolus followed by 0.05 mg/kg/h IV",
      "infusion. Adults: 0.1-0.15 mg/kg IV (bolus or 15-20 min infusion)",
      "or 5/7.5/10 mg fixed oral doses by body-weight band (<60 / 60-80 /",
      ">80 kg)."
    ),
    regions        = "Netherlands",
    notes          = paste(
      "Per Cella 2012 Table 1. Sex split not reported across all four",
      "studies; Study 94113 is explicitly 20 healthy males.",
      "Apparent oral clearance (CL/F): the paper does not separately",
      "estimate bioavailability; the IV and oral data share a single CL",
      "parameter and the depot bioavailability defaults to 1 (true F is",
      "absorbed into CL)."
    )
  )

  ini({
    # Structural parameters from Table 2 ('Mean' column for Model 1
    # 'Infants, toddlers and adults'). The published table units are
    # 'CL (l min-1 kg^0.75)', 'Vc (l)', 'Q (l min-1)', 'Vp (l)' and
    # 'Ka (h-1)'. Two encoding decisions deviate from the literal table:
    #   (a) the median weight for the CL allometric normalisation is not
    #       stated in the paper; sidecar response 001 Q1 instructed to
    #       validate against Table 3, and the validation showed Wmed =
    #       70 kg reproduces the published Model-1-extrapolated AUC0-180
    #       with a geometric-mean pred/obs ratio of 0.98 across 16 non-
    #       outlier rows (vs 0.74 for Wmed = 31 kg);
    #   (b) Vc per-kg encoding (sidecar response 001 Q2 = A): the
    #       literal absolute reading gives k10 ~ 0.75/min for 70 kg
    #       adults, ~100x too fast vs published midazolam PK. The per-kg
    #       reading (Vc = 0.312 * WT L) is consistent with the De Wildt
    #       2002 prior the authors cited for Model 2 (Vc = 0.38 L/kg)
    #       and recovers physically plausible kinetics for both infants
    #       and adults.
    lka <- log(8.21 / 60);   label("First-order absorption rate constant (Ka, 1/min)")             # Cella 2012 Table 2 (Ka = 8.21 h^-1, converted to /min)
    lcl <- log(0.234);       label("Clearance at WT = 70 kg reference (CL, L/min)")                # Cella 2012 Table 2 (CL = 0.234 with allometric WT effect)
    lvc <- log(0.312);       label("Central volume of distribution per kg body weight (Vc/WT, L/kg)") # Cella 2012 Table 2 (Vc = 0.312; per-kg interpretation per sidecar response 001 Q2)
    lq  <- log(1.34);        label("Inter-compartmental clearance (Q, L/min)")                     # Cella 2012 Table 2
    lvp <- log(16.5);        label("Peripheral volume of distribution (Vp, L)")                    # Cella 2012 Table 2

    allo_cl <- fixed(0.75);  label("Allometric exponent of WT on CL (unitless)")                   # Cella 2012 Methods + Mahmood 1996 ref [46] (fixed at classic 0.75)

    # Inter-individual variability (Table 2 reports CV%; omega^2 = log(CV^2 + 1))
    #   CL : 39.9% CV -> log(0.399^2 + 1) = 0.14776
    #   Vp : 58.5% CV -> log(0.585^2 + 1) = 0.29446
    etalcl ~ 0.14776                                                                                # Cella 2012 Table 2 (IIV CL 39.9% CV)
    etalvp ~ 0.29446                                                                                # Cella 2012 Table 2 (IIV Vp 58.5% CV)

    # Residual error (40.0% CV proportional, linear-scale SD = 0.40)
    propSd <- 0.40;          label("Proportional residual error (SD, fraction)")                   # Cella 2012 Table 2 (residual 40.0% CV, proportional model)
  })
  model({
    # Individual PK parameters.
    # CL: allometric WT-scaling, reference 70 kg.
    # Vc: per-kg linear scaling (Vc/WT is the estimated per-kg coefficient).
    # Vp: constant population-typical absolute volume (no covariate effect).
    # Q and Ka: constant population-typical values (no covariate effect).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^allo_cl
    vc <- exp(lvc) * WT
    q  <- exp(lq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Plasma midazolam concentration. With dose in mg and Vc in L, Cc is in
    # mg/L (numerically equivalent to ug/mL and to ng/mL * 1e-3).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
