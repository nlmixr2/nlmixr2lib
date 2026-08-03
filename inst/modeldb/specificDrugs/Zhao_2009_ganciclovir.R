Zhao_2009_ganciclovir <- function() {
  description <- paste(
    "Two-compartment population PK model for ganciclovir after oral valganciclovir",
    "in pediatric renal transplant patients (Zhao 2009), parameterized in apparent",
    "oral clearances and volumes, with first-order absorption after a lag time and",
    "a two-component additive apparent clearance built from a steep",
    "creatinine-clearance power term plus a weight-proportional term.",
    "Parameters transcribed from the Yang 2023 ganciclovir / valganciclovir",
    "population-PK model repository (Table 3), not from the primary publication;",
    "re-verify against Zhao 2009 when the primary is obtained -- the transcribed",
    "creatinine-clearance exponent of 2.93 is unusually steep and is flagged in the",
    "validation vignette.",
    sep = " "
  )
  reference <- paste(
    "Zhao W, Baudouin V, Zhang D, Deschenes G, Le Guellec C, Jacqz-Aigrain E.",
    "Population pharmacokinetics of ganciclovir following administration of",
    "valganciclovir in paediatric renal transplant patients.",
    "Clin Pharmacokinet. 2009;48(5):321-328. doi:10.2165/00003088-200948050-00004.",
    "Parameters transcribed from Yang W, Mak W, Gwee A, Gu M, Wu Y, Shi Y, He Q,",
    "Xiang X, Han B, Zhu X. Establishment and Evaluation of a Parametric Population",
    "Pharmacokinetic Model Repository for Ganciclovir and Valganciclovir.",
    "Pharmaceutics. 2023;15(7):1801. doi:10.3390/pharmaceutics15071801 (Table 3).",
    sep = " "
  )
  vignette <- "Yang_2023_ganciclovir_model_repository"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives the second, additive component of apparent oral clearance:",
        "3.62 * (BW/28) L/h, referenced to the cohort median 28 kg. Weight enters",
        "linearly (exponent 1) and additively rather than as an allometric",
        "multiplier, and does not affect the volumes. Cohort median 28 kg,",
        "range 12-76 kg; mean 34 (SD 19) kg (Yang 2023 Table 2).",
        sep = " "
      ),
      source_name        = "BW"
    ),
    CRCL = list(
      description        = "Creatinine clearance",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Yang 2023 Table 3 footnote defines CLcr as creatinine clearance in mL/min",
        "(raw, NOT BSA-normalized), distinct from the BSA-normalized `CrCL` used by",
        "Franck 2021 in the same repository. Drives the first, additive component of",
        "apparent oral clearance: 8.04 * (CLcr/89)^2.93 L/h, referenced to 89 mL/min.",
        "The exponent 2.93 as printed in Yang 2023 Table 3 is far steeper than any",
        "other renal-function effect in the repository (the others range from 0.29",
        "to 1.08) and produces a roughly 3.7-fold swing in CL/F across a two-fold",
        "CLcr range; it is encoded verbatim here per the trust-the-printed-equation",
        "rule, and flagged in the validation vignette for verification against the",
        "primary publication.",
        sep = " "
      ),
      source_name        = "CLcr"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 22L,
    n_studies      = 1L,
    n_observations = 164L,
    age_median     = "9 years (range 3-17); mean 10 (SD 5) years",
    weight_median  = "28 kg (range 12-76); mean 34 (SD 19) kg",
    sex_female_pct = 50,
    race_ethnicity = "Not reported.",
    disease_state  = "Pediatric renal transplant patients.",
    co_medication  = paste(
      "Prednisone and mycophenolate mofetil were tested as covariates and not",
      "retained. Yang 2023 Section 4.2 attributes this model's higher simulated",
      "concentrations relative to the other pediatric valganciclovir models to the",
      "cohort's mycophenolate mofetil immunosuppression, which is known to reduce",
      "the renal clearance of ganciclovir.",
      sep = " "
    ),
    dose_range     = paste(
      "Prophylactic therapy: oral valganciclovir 900 mg/24 h. Pre-emptive therapy:",
      "IV ganciclovir 5 mg/kg/12 h for 15 days followed by oral valganciclovir",
      "10 mg/kg/12 h for 3 months.",
      sep = " "
    ),
    regions        = "France (prospective study).",
    bioassay       = "HPLC, LLOQ 0.25 ug/mL.",
    notes          = paste(
      "Demographics and dosing from Yang 2023 Table 2. Intensive sampling at",
      "0 (pre-dose), 1, 2, 4, 6, 8 and 12 h post dose. Covariates tested: weight,",
      "age, height, CLcr, AST, ALT, serum protein, prednisone and mycophenolate",
      "mofetil; CLcr and weight were retained (Yang 2023 Table 4). This study did",
      "NOT convert valganciclovir doses to ganciclovir equivalents; oral doses are",
      "supplied as valganciclovir milligrams and the molecular-weight ratio plus",
      "the true bioavailability are absorbed into the apparent parameters.",
      "Simulation targets in the primary study were an AUC0-24h of 45 mg*h/L and",
      "trough concentrations of 0.5 and 1 mg/L.",
      sep = " "
    )
  )

  ini({
    # Structural PK -- Yang 2023 Table 3, Zhao et al. (2009) row. Apparent oral
    # clearance is the SUM of two components (Yang 2023 Table 3 prints
    # CL/F = 8.04 * (CLcr/89)^2.93 + 3.62 * (BW/28)), encoded with the registered
    # multi-component clearance names: `lcl_renal` for the creatinine-clearance-driven
    # arm and `lcl_nonren` for the weight-proportional arm. Reference subject:
    # CLcr = 89 mL/min, WT = 28 kg -> CL/F = 8.04 + 3.62 = 11.66 L/h.
    lcl_renal  <- log(8.04); label("Renal-function arm of apparent oral clearance at CLcr = 89 mL/min (L/h)") # Yang 2023 Table 3 (Zhao 2009): 8.04 * (CLcr/89)^2.93
    lcl_nonren <- log(3.62); label("Weight-proportional arm of apparent oral clearance at WT = 28 kg (L/h)")  # Yang 2023 Table 3 (Zhao 2009): + 3.62 * (BW/28)
    lvc        <- log(5.2) ; label("Apparent central volume Vc/F (L)")                                        # Yang 2023 Table 3 (Zhao 2009): Vc/F = 5.2
    lq         <- log(3.97); label("Apparent inter-compartmental clearance Q/F (L/h)")                        # Yang 2023 Table 3 (Zhao 2009): Q/F = 3.97
    lvp        <- log(30.7); label("Apparent peripheral volume Vp/F (L)")                                     # Yang 2023 Table 3 (Zhao 2009): Vp/F = 30.7
    lka        <- log(0.369); label("First-order oral absorption rate constant (ka, 1/h)")                    # Yang 2023 Table 3 (Zhao 2009): Ka = 0.369
    ltlag      <- log(0.743); label("Absorption lag time (Tlag, h)")                                          # Yang 2023 Table 3 (Zhao 2009): Tlag = 0.743

    # Covariate effects. The CLcr exponent 2.93 is transcribed verbatim from
    # Yang 2023 Table 3 (see covariateData[[CRCL]]$notes for the flag). The weight
    # exponent is 1 (linear) as printed and is encoded as fixed.
    e_crcl_cl <- 2.93    ; label("Power exponent of CLcr on the renal arm of CL/F (unitless; reference 89 mL/min)")   # Yang 2023 Table 3 (Zhao 2009): (CLcr/89)^2.93
    e_wt_cl   <- fixed(1); label("Power exponent of WT on the weight arm of CL/F (unitless; linear, reference 28 kg)") # Yang 2023 Table 3 (Zhao 2009): (BW/28)

    # Between-subject variability. Yang 2023 Methods: %CV = sqrt(omega^2) * 100%,
    # so variance = (BSV% / 100)^2. The source reports a single BSV for CL/F, which
    # is applied multiplicatively to the summed two-component clearance. Q/F, Vp/F
    # and Tlag carry no BSV in the source table.
    etalcl ~ 0.05678689  # Yang 2023 Table 3 (Zhao 2009): BSV CL/F = 23.83% -> 0.2383^2
    etalvc ~ 0.33895684  # Yang 2023 Table 3 (Zhao 2009): BSV Vc/F = 58.22% -> 0.5822^2
    etalka ~ 0.10400625  # Yang 2023 Table 3 (Zhao 2009): BSV ka   = 32.25% -> 0.3225^2

    # Residual unexplained variability. Yang 2023 Table 3 reports 20.93%
    # "exponential error"; an exponential residual on the natural scale is
    # equivalent to an additive residual on the log scale, which nlmixr2 expresses
    # as a proportional error in linear space (matching the encoding used for the
    # other exponential-residual models in this repository).
    propSd <- 0.2093; label("Proportional residual error (fraction; exponential residual in the source)")  # Yang 2023 Table 3 (Zhao 2009): 20.93% exponential error
  })

  model({
    # Two-component additive apparent clearance; the single reported CL/F eta is
    # applied to the sum.
    cl_renal  <- exp(lcl_renal)  * (CRCL / 89)^e_crcl_cl
    cl_nonren <- exp(lcl_nonren) * (WT / 28)^e_wt_cl
    cl <- (cl_renal + cl_nonren) * exp(etalcl)

    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Bioavailability is absorbed into the apparent parameters (CL/F, Vc/F, Q/F,
    # Vp/F), so no explicit f(depot) is applied. Oral valganciclovir doses enter
    # `depot` as valganciclovir milligrams.
    alag(depot) <- tlag

    # Dose in mg, apparent volume in L -> concentration in mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
