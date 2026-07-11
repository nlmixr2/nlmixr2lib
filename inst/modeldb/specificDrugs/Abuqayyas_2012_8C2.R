Abuqayyas_2012_8C2 <- function() {
  description <- "Preclinical (mouse). Two-compartment IV bolus population PK model of 8C2, a murine IgG1 anti-topotecan monoclonal antibody, in C57BL/6 wild-type, Fc-gamma-RI/RIII knockout, and Fc-gamma-RIIb knockout mice (Abuqayyas 2012). Data were pooled across all three strains and three IV bolus dose levels (0.04, 0.1, 0.4 mg/kg); strain was screened and found not statistically significant on any structural parameter, so no covariate effects are retained. Structural PK parameters were reported per kg body weight (L/kg and L/day/kg), which corresponds to linear body-weight scaling with a fixed exponent of 1."
  reference <- "Abuqayyas L, Balthasar JP. Application of knockout mouse models to investigate the influence of Fc gamma R on the tissue distribution and elimination of 8C2, a murine IgG1 monoclonal antibody. Int J Pharm. 2012 Nov 15;439(1-2):8-16. doi:10.1016/j.ijpharm.2012.09.042. PMID: 23018115."
  vignette <- "Abuqayyas_2012_8C2"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Body weight enters the model as linear scaling on volumes and clearances (fixed exponent of 1). Abuqayyas 2012 reported the structural parameters in per-kg form (Table 2: L/kg for volumes, L/day/kg for clearances); body weight itself was not evaluated as a discrete covariate (mice were 20-38 g).",
      source_name        = "WT"
    )
  )

  covariatesDataExcluded <- list(
    STRAIN = list(
      description = "Mouse strain: C57BL/6 wild-type, B6.129P2-Fcer1g<tm1Rav> (Fc-gamma-RI/RIII gamma-chain knockout), or B6.129S4-Fcgr2b<tm1TtK> (Fc-gamma-RIIb knockout).",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened via forward selection with backward elimination on all structural parameters (Section 2.6). Pair-wise MVOF testing of CLc between each knockout strain and WT (Table 3): WT vs Fc-gamma-RI/RIII delta-MVOF = 0.22 (p = 0.639), WT vs Fc-gamma-RIIb delta-MVOF = 0.02 (p = 0.888). Not retained in the final model. Reference category (had the effect been retained): C57BL/6 wild-type."
    )
  )

  population <- list(
    species        = "mouse (C57BL/6 wild-type; B6.129P2-Fcer1g<tm1Rav> Fc-gamma-RI/RIII knockout; B6.129S4-Fcgr2b<tm1TtK> Fc-gamma-RIIb knockout)",
    n_subjects     = "n = 3-4 per dose per strain (three strains x three dose levels)",
    n_studies      = 1,
    weight_range   = "20-38 g (plasma PK study); 20-25 g (tissue distribution study)",
    disease_state  = "Naive laboratory mice (no disease model); 8C2 has no known murine target so behaves as a non-target-binding tracer IgG1.",
    dose_range     = "0.04, 0.1, and 0.4 mg/kg IV bolus of 125I-labelled 8C2 (tracer ~10 uCi/mouse).",
    regions        = "United States (University at Buffalo; mice from Taconic Laboratories, Hudson NY).",
    notes          = "Plasma sampled at 1, 3, 8 h and 1, 2, 4, 7, 10 days from retro-orbital plexus or sub-mandibular vein; radioactivity counted by gamma counter (LKB Wallac 1272) and decay-corrected. Concentrations reported in nM assuming an IgG molecular weight (see vignette for the nM <-> mg/L conversion). Mice were maintained on autoclaved KI-water (0.2 g/L) starting 2 days pre-injection to block thyroidal uptake of free iodine."
  )

  ini({
    # Structural PK parameters -- Abuqayyas 2012 Table 2 (all per kg body weight)
    lvc <- log(0.057)   ; label("Central volume of distribution (Vc, L/kg)")            # Abuqayyas 2012 Table 2: Vc = 0.057 L/kg (%SEM 4.0)
    lvp <- log(0.0996)  ; label("Peripheral (tissue) volume of distribution (Vt, L/kg)") # Abuqayyas 2012 Table 2: Vt = 0.0996 L/kg (%SEM 6.4)
    lcl <- log(0.00543) ; label("Apparent clearance (CLc, L/day/kg)")                    # Abuqayyas 2012 Table 2: CLc = 0.00543 L/day/kg (%SEM 17.5)
    lq  <- log(0.0598)  ; label("Distribution clearance (CLd, L/day/kg)")                # Abuqayyas 2012 Table 2: CLd = 0.0598 L/day/kg (%SEM 9.7)

    # Inter-animal variability -- exponential model, variances directly from Table 2
    # (Abuqayyas 2012 reports CV% = sqrt(omega^2) as a small-omega approximation)
    # No IIV was reported on Vt.
    etalvc ~ 0.0633  # Abuqayyas 2012 Table 2: omega^2 Vc  = 0.0633 (25.2% CV, %SEM 16.1)
    etalq  ~ 0.218   # Abuqayyas 2012 Table 2: omega^2 CLd = 0.218  (46.7% CV, %SEM 21.2)
    etalcl ~ 0.387   # Abuqayyas 2012 Table 2: omega^2 CLc = 0.387  (62.6% CV, %SEM 40.3)

    # Residual variability -- proportional (constant CV) error model
    # Table 2 reports sigma^2_prop = 0.0181 (variance); nlmixr2 propSd is a standard deviation
    propSd <- sqrt(0.0181); label("Proportional residual error (fraction)")   # Abuqayyas 2012 Table 2: sigma^2 prop = 0.0181 (13.4% CV, %SEM 13.7)
  })

  model({
    # Individual PK parameters -- per-kg reference values scaled linearly by body weight (kg)
    vc <- exp(lvc + etalvc) * WT
    vp <- exp(lvp)          * WT
    cl <- exp(lcl + etalcl) * WT
    q  <- exp(lq  + etalq)  * WT

    # Micro-constants for two-compartment mammillary model (Fig. 1, Section 2.6)
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Concentration -- dose in mg, volume in L -> mg/L (equivalent to ug/mL)
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
