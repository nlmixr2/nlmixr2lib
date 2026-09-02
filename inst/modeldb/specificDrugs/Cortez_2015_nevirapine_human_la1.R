Cortez_2015_nevirapine_human_la1 <- function() {
  description <- paste(
    "Human infant projection of the NVP LA-1 long-acting injectable",
    "nevirapine microparticle suspension for prevention of mother-to-child",
    "HIV transmission in breastfeeding neonates. Two-compartment model with",
    "first-order absorption in which the SYSTEMIC parameters (Vc/F and CL/F)",
    "are the published breastfeeding-infant oral nevirapine population",
    "estimates from two Ugandan studies, while the PRESYSTEMIC absorption",
    "rate constant and the peripheral distribution parameters are carried",
    "over from the authors' own Sprague-Dawley rat fit",
    "(Cortez_2015_nevirapine_rat_la1) and allometrically scaled from the",
    "0.33 kg rat reference weight with fixed exponents (-0.25 on ka, 1 on",
    "volumes, 0.75 on clearances). Deterministic naive-pooled simulation",
    "model: the publication reports neither between-subject variability nor",
    "a residual-error model for the human projection, so none is encoded.",
    "Reproduces Cortez 2015 Figure 4A (200 mg, 2.0 kg) and 4B (200 mg,",
    "3.9 kg).",
    sep = " ")
  reference <- paste(
    "Cortez JM Jr, Quintero R, Moss JA, Beliveau M, Smith TJ, Baum MM.",
    "Pharmacokinetics of injectable, long-acting nevirapine for HIV",
    "prophylaxis in breastfeeding infants.",
    "Antimicrob Agents Chemother. 2015;59(1):59-66.",
    "doi:10.1128/AAC.03906-14.",
    "Presystemic and peripheral parameters carried from the companion rat",
    "fit; see modellib('Cortez_2015_nevirapine_rat_la1').",
    sep = " ")
  vignette <- "Cortez_2015_nevirapine"
  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  compartmentData <- list(
    depot       = list(analyte = "nevirapine", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "nevirapine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "nevirapine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric size covariate on the RAT-DERIVED parameters only (ka,",
        "Vp, CLp), referenced to the 0.33 kg mean weight of the study-1",
        "rats. The human systemic parameters Vc/F and CL/F are the",
        "published infant population estimates and are NOT weight-scaled --",
        "this composition was recovered by reproducing Cortez 2015 Figure",
        "4A and 4B and is documented in the vignette Assumptions and",
        "deviations. Simulated infant weights were 2.0 and 3.9 kg, the",
        "range of the Ugandan newborn cohorts."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 2L,
    age_range      = "Neonates from birth (WHO prophylaxis window is birth to 4-6 weeks of age)",
    weight_range   = "2.0-3.9 kg",
    disease_state  = paste(
      "HIV-uninfected breastfeeding infants born to HIV-1-infected mothers,",
      "receiving nevirapine as prophylaxis against mother-to-child",
      "transmission."
    ),
    dose_range     = "Single 200 mg subcutaneous injection of NVP LA-1 (Cortez 2015 Figure 4A and 4B)",
    regions        = "Uganda (source of the systemic infant PK estimates)",
    notes          = paste(
      "Cortez 2015 Methods 'Human simulations': \"Infant parameters for the",
      "human simulations were compiled from recent clinical studies in",
      "Uganda examining the population PK of single-dose NVP in",
      "breastfeeding infants of HIV-positive mothers (20, 21). The",
      "newborns' weights ranged from 2.0 to 3.9 kg, and the following",
      "systemic parameters were determined from population PK estimates:",
      "V/bioavailability (F), 20.0 liters; CL/F, 0.21 liters h-1; and",
      "t1/2, 66 h. These systemic parameters were combined with presystemic",
      "rat parameters, essentially the absorption rate constant, ka, in",
      "order to simulate human exposure. Simulations were performed using a",
      "naive-pooled approach to estimate the mean plasma drug",
      "concentrations in humans.\" The bioavailable fraction was assumed to",
      "be 100% (Cortez 2015 Results 'Human simulation'). The prophylactic",
      "target window motivating the simulation is 0.2 to 3.0 ug/mL",
      "(Cortez 2015 Discussion 'Implications for HIV prevention in",
      "breastfeeding infants')."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # SYSTEMIC parameters -- human infant population estimates, taken
    # as-is from Cortez 2015 Methods 'Human simulations'. They are NOT
    # weight-scaled: reproducing Figure 4A and 4B requires Vc/F = 20 L
    # and CL/F = 5.04 L/day at both 2.0 and 3.9 kg (see the vignette
    # Assumptions and deviations for the digitised check).
    #
    # Internal consistency check on the transcription: the paper also
    # reports t1/2 = 66 h for these infants, and
    #   log(2) * 20 L / 0.21 L/h = 66.0 h.
    # ------------------------------------------------------------------
    lvc <- log(20)   ; label("Apparent central volume V/F (L), human infant")                  # Cortez 2015 Methods 'Human simulations': V/bioavailability (F), 20.0 liters
    lcl <- log(5.04) ; label("Apparent clearance CL/F (L/day), human infant")                  # Cortez 2015 Methods 'Human simulations': CL/F = 0.21 L/h; 0.21 * 24 = 5.04 L/day in this model's day time base

    # ------------------------------------------------------------------
    # PRESYSTEMIC and PERIPHERAL parameters -- carried from the study-1
    # rat fit (Cortez 2015 Table 2, study 1 column) and allometrically
    # scaled in model() from the 0.33 kg rat reference weight. ka is the
    # rate-limiting release of drug from the microparticles, which is a
    # formulation property and therefore the one quantity the rat study
    # is meant to supply (Cortez 2015 Discussion 'Implications for HIV
    # prevention in breastfeeding infants').
    # ------------------------------------------------------------------
    lka <- log(1.89)  ; label("Absorption rate constant ka (1/day) at the 0.33 kg rat reference weight")             # Cortez 2015 Table 2 study 1: Ka = 1.89 /day
    lvp <- log(14.3)  ; label("Apparent peripheral volume Vp/F (L) at the 0.33 kg rat reference weight")             # Cortez 2015 Table 2 study 1: Vp/F = 14.3 L
    lq  <- log(0.541) ; label("Apparent inter-compartmental clearance CLp/F (L/day) at the 0.33 kg rat reference")   # Cortez 2015 Table 2 study 1: CLp/F = 0.541 L/day

    lfdepot <- fixed(log(1)) ; label("Subcutaneous bioavailability F (fraction)")                             # Cortez 2015 Results 'Human simulation': "The bioavailable fraction was assumed to be 100% in these simulations."

    # ------------------------------------------------------------------
    # Allometric exponents -- Cortez 2015 Methods 'Structural PK model
    # buildup': "The nominal covariate effect was -0.25 on absorption
    # rate, 1 on volumes, and 0.75 on clearance (19)." Imposed, not
    # estimated, hence fixed().
    # ------------------------------------------------------------------
    e_wt_ka <- fixed(-0.25) ; label("Allometric exponent on ka (unitless)")    # Cortez 2015 Methods 'Structural PK model buildup'
    e_wt_vp <- fixed(1)     ; label("Allometric exponent on Vp (unitless)")    # Cortez 2015 Methods 'Structural PK model buildup'
    e_wt_q  <- fixed(0.75)  ; label("Allometric exponent on CLp (unitless)")   # Cortez 2015 Methods 'Structural PK model buildup'

    # ------------------------------------------------------------------
    # No IIV and no residual-error model: Cortez 2015 simulated mean
    # plasma concentrations by a naive-pooled approach and reports
    # neither an omega nor a sigma for the human projection.
    # ------------------------------------------------------------------
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters.
    #    Human systemic parameters used at their published values;
    #    rat-derived presystemic / peripheral parameters allometrically
    #    scaled from the 0.33 kg rat reference weight.
    # ------------------------------------------------------------------
    vc <- exp(lvc)
    cl <- exp(lcl)
    ka <- exp(lka) * (WT / 0.33)^e_wt_ka
    vp <- exp(lvp) * (WT / 0.33)^e_wt_vp
    q  <- exp(lq)  * (WT / 0.33)^e_wt_q

    # ------------------------------------------------------------------
    # 2. Micro-constants.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ------------------------------------------------------------------
    # 3. ODE system.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # ------------------------------------------------------------------
    # 4. Bioavailability (assumed 100%).
    # ------------------------------------------------------------------
    f(depot) <- exp(lfdepot)

    # ------------------------------------------------------------------
    # 5. Observation. Deterministic simulation; no residual error.
    # ------------------------------------------------------------------
    Cc <- central / vc
  })
}
