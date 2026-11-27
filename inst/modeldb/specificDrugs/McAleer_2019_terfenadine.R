McAleer_2019_terfenadine <- function() {
  description <- paste(
    "QSP (in vitro microphysiological organ-chip). Terfenadine",
    "pharmacokinetic-pharmacodynamic model for a pumpless recirculating",
    "heart:liver (and heart-only) human organ-chip, built from iPSC-derived",
    "cardiomyocytes and primary human hepatocytes (McAleer 2019). Five-state",
    "PK: terfenadine distributes between a recirculating medium compartment",
    "(where drug is dosed), a PDMS-adsorption sink, a cardiac compartment and",
    "a liver compartment, with a saturable efflux out of the cardiac",
    "compartment; terfenadine is metabolised to fexofenadine in the cardiac",
    "and liver compartments, which accumulates in a fifth fexofenadine state.",
    "An indirect-response PD state driven by the cardiac terfenadine",
    "concentration produces the field-potential-duration (FPD) prolongation",
    "response, FPD = 1 + effect (fold of vehicle control). The heart-only chip",
    "is the same model with the liver-transfer rate constants (kl1, kl2, kfl)",
    "set to zero. Deterministic typical-value model; the source reports no",
    "between-chip random effects but tabulates additive residual-error SDs,",
    "which are encoded.",
    sep = " "
  )
  reference <- paste(
    "McAleer CW, Pointon A, Long CJ, Brighton RL, Wilkin BD, Bridges LR,",
    "Narasimhan Sriram N, Fabre K, McDougall R, Muse VP, Mettetal JT,",
    "Srivastava A, Williams D, Schnepper MT, Roles JL, Shuler ML, Hickman JJ,",
    "Ewart L. On the potential of in vitro organ-chip models to define",
    "temporal pharmacokinetic-pharmacodynamic relationships. Sci Rep.",
    "2019;9(1):9619. doi:10.1038/s41598-019-45656-4. PMCID: PMC6610665.",
    "Model structure (ODEs) and final parameter values from the deposited",
    "Phoenix NLME code (Supplementary Information, Supplementary Text 1);",
    "equations 1-2 of the main text; parameter estimates from the",
    "Supplementary Information 'Parameter values' block (the self-consistent",
    "final estimates -- see the vignette Errata on why the code block is used",
    "rather than the Supplementary Table 1 Estimate column).",
    sep = " "
  )
  vignette <- "McAleer_2019_terfenadine"
  units <- list(
    time = "min",
    dosing = "umol/L (nominal terfenadine concentration added to the 2 mL recirculating medium)",
    concentration = "umol/L"
  )

  covariateData <- list()

  # Issue #482: what each ODE state holds. All terfenadine states are
  # concentrations (umol/L) in medium-equivalent units (the published Phoenix
  # deriv transcript carries no compartment volumes); `fexo` is the
  # fexofenadine metabolite; `effect` is the dimensionless FPD indirect
  # response (fold of vehicle control). Verified against the deposited Phoenix
  # code and Fig 4a of McAleer 2019.
  compartmentData <- list(
    media = list(analyte = "terfenadine", units = "umol/L", specimen = "culture medium", verified = TRUE),
    pdms = list(analyte = "terfenadine", units = "umol/L", specimen = "not applicable", verified = TRUE),
    cardiac = list(analyte = "terfenadine", units = "umol/L", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "terfenadine", units = "umol/L", specimen = "tissue", verified = TRUE),
    fexo = list(analyte = "fexofenadine", units = "umol/L", specimen = "culture medium", verified = TRUE),
    effect = list(
      analyte = "field-potential-duration response",
      units = "fold of vehicle control (dimensionless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  # Paper-mechanistic in-vitro organ-chip states. `liver` and `effect` are
  # canonical (compartment-names.md); the recirculating medium, the PDMS
  # adsorption sink, the cardiomyocyte compartment and the fexofenadine
  # metabolite pool have no canonical PK role, so they are declared here.
  # Every state holds a concentration (umol/L); the published Phoenix deriv
  # transcript carries no explicit compartment volumes, so the states are in
  # medium-equivalent concentration units.
  paper_specific_compartments <- c("media", "pdms", "cardiac", "fexo")

  population <- list(
    species = paste(
      "in vitro (human iPSC-derived cardiomyocytes + primary human",
      "hepatocytes; pumpless recirculating heart:liver / heart-only",
      "microphysiological organ-chip)"
    ),
    n_subjects = NA,
    disease_state = "in vitro cardiac electrophysiology / hepatic metabolism assay",
    dose_range = "5 umol/L nominal terfenadine added to 2 mL recirculating medium",
    notes = paste(
      "Chip replicates: n = 2 for the PK (concentration) measurements and",
      "n >= 3 for the FPD measurements, at 0, 0.5, 1, 4, 12, 16 and 24 hours.",
      "Total system volume 2 mL; drug dosed at the medium reservoir closest",
      "to the liver chamber. Parameters were fit simultaneously across the",
      "heart:liver and heart-only chips in Phoenix NLME 7.0."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # All eleven structural rate / saturation parameters are the FINAL
    # estimates from the deposited Phoenix code 'Parameter values' block
    # (Supplementary Text 1). They are preferred over the Supplementary
    # Table 1 'Estimate' column because several Table 1 estimates (KP1,
    # KL1, KC2, IC50) fall OUTSIDE their own tabulated 95% CIs, whereas
    # every code-block value falls inside the Table 1 CI -- i.e. the Table 1
    # Estimate column is partially misaligned and the code block holds the
    # true fitted values. See the vignette Errata. Rate constants are 1/min;
    # the model time unit is minutes. No between-chip IIV was reported, so
    # these are deterministic typical values (no eta partners).
    # ------------------------------------------------------------------
    kp1 <- 0.0035097   ; label("Medium -> PDMS adsorption rate constant (1/min)")               # Supp Text 1 KP1
    kp2 <- 0.000526081 ; label("PDMS -> medium desorption rate constant (1/min)")               # Supp Text 1 KP2
    kc1 <- 0.00185601  ; label("Medium -> cardiac uptake rate constant (1/min)")                # Supp Text 1 KC1
    kc2 <- 0.0780807   ; label("Cardiac -> medium saturable efflux maximal rate constant (1/min)") # Supp Text 1 KC2
    kl1 <- 0.0256725   ; label("Medium -> liver uptake rate constant (1/min)")                  # Supp Text 1 KL1
    kl2 <- 0.0429636   ; label("Liver -> medium efflux rate constant (1/min)")                  # Supp Text 1 KL2
    kfl <- 0.0167355   ; label("Liver terfenadine -> fexofenadine metabolism rate constant (1/min)") # Supp Text 1 Kfl
    kfc <- 0.000912019 ; label("Cardiac terfenadine -> fexofenadine metabolism rate constant (1/min)") # Supp Text 1 Kfc
    ic50 <- 0.0699277  ; label("Cardiac terfenadine concentration at half-maximal efflux saturation (umol/L)") # Supp Text 1 IC50
    kon <- 0.147464    ; label("FPD effect production rate constant, driven by cardiac terfenadine (1/min)") # Supp Text 1 kon
    koff <- 0.00413283 ; label("FPD effect loss rate constant (1/min)")                         # Supp Text 1 koff

    # ------------------------------------------------------------------
    # Additive residual-error SDs from the deposited Phoenix 'error()'
    # statements for the HEART:LIVER chip (Supplementary Text 1). The
    # heart-only chip used slightly different SDs (cardiac 0.01338, medium
    # 1.07091, fexo 0.0323472, FPD 0.38903); the heart:liver values are
    # carried here as the default full model. Units match each output.
    # ------------------------------------------------------------------
    addSd_CterfMedia   <- 1.17403    ; label("Additive residual SD, medium terfenadine (umol/L)")   # Supp Text 1 CEps_media
    addSd_CterfCardiac <- 0.00642974 ; label("Additive residual SD, cardiac terfenadine (umol/L)")  # Supp Text 1 CEps_cardiac
    addSd_CterfLiver   <- 0.0101382  ; label("Additive residual SD, liver terfenadine (umol/L)")    # Supp Text 1 CEps_liver
    addSd_Cfexo        <- 1.84938    ; label("Additive residual SD, fexofenadine (umol/L)")         # Supp Text 1 CEps_fexo
    addSd_FPD          <- 0.995963   ; label("Additive residual SD, FPD (fold of control)")         # Supp Text 1 FPDEps_CL
  })

  model({
    # Saturable fraction of the cardiac efflux (main-text Eqs 1-2:
    # Kmediumtoheart = KC1; Khearttomedium = KC2 * 1/(1 + Cardiac/IC50)).
    cardiac_sat <- 1 / (1 + cardiac / ic50)

    # --- Terfenadine PK (deposited Phoenix deriv block, heart:liver) ---
    # Medium is the dosing compartment. The heart-only chip is recovered by
    # setting kl1 = kl2 = kfl = 0 (the liver chamber is eliminated).
    d/dt(media) <- -(kp1 + kc1 + kl1) * media +
      kp2 * pdms + kc2 * cardiac_sat * cardiac + kl2 * liver
    d/dt(pdms) <- kp1 * media - kp2 * pdms
    d/dt(cardiac) <- kc1 * media - kc2 * cardiac_sat * cardiac - kfc * cardiac
    d/dt(liver) <- kl1 * media - kl2 * liver - kfl * liver

    # --- Fexofenadine metabolite accumulation ---
    d/dt(fexo) <- kfl * liver + kfc * cardiac

    # --- FPD indirect-response PD, driven by cardiac terfenadine ---
    d/dt(effect) <- kon * cardiac - koff * effect

    # --- Observations ---
    CterfMedia <- media
    CterfCardiac <- cardiac
    CterfLiver <- liver
    Cfexo <- fexo
    FPD <- 1 + effect

    CterfMedia ~ add(addSd_CterfMedia)
    CterfCardiac ~ add(addSd_CterfCardiac)
    CterfLiver ~ add(addSd_CterfLiver)
    Cfexo ~ add(addSd_Cfexo)
    FPD ~ add(addSd_FPD)
  })
}
