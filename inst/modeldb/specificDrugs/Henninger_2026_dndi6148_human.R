Henninger_2026_dndi6148_human <- function() {
  description <- paste(
    "Human translation of the Henninger 2026 murine DNDI-6148 target-site PK/PD model,",
    "used to predict an efficacious oral dose for cutaneous leishmaniasis. Plasma PK is",
    "one-compartment with first-order absorption and LINEAR elimination, its three typical",
    "values obtained by single-species allometric scaling from the mouse fit; the murine",
    "saturable clearance and dose-dependent bioavailability were judged mouse-specific and",
    "are deliberately absent. The murine skin penetration coefficient and the murine",
    "sigmoidal Emax parasite-killing component are carried over unchanged, so free skin",
    "concentration drives parasite clearance from a baseline burden set equal to the mouse.",
    "Every parameter is fixed: none was estimated from human data. See",
    "modellib('Henninger_2026_dndi6148_mouse') for the murine model these values come from."
  )
  reference <- paste(
    "Henninger RH, Schouten WM, Arana B, Gillon JY, Mowbray CE, Kratz JM,",
    "Van Bocxlaer K, Dorlo TPC. (2026). Translational Pharmacokinetic-Pharmacodynamic",
    "Modeling and Efficacious Human Dose Prediction of DNDI-6148 for the Treatment of",
    "Cutaneous Leishmaniasis. Clin Transl Sci 19(4):e70535.",
    "doi:10.1111/cts.70535."
  )
  vignette <- "Henninger_2026_dndi6148"

  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list()

  # Body weight is NOT a covariate here. The allometric exponents 0.75 / 1 /
  # -0.25 were used once, to bridge the 22 g mouse to a 70 kg human; the human
  # simulations were then run for 70 kg adults only (Methods 2.2.6), with dose
  # levels expressed in mg/kg. Carrying WT scaling inside the human model would
  # be an extrapolation the authors did not make.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Used only for the mouse-to-human allometric bridge (Methods 2.2.5, fixed",
        "exponents 0.75 for CL/F, 1.00 for Vd/F and -0.25 for ka) and to convert the",
        "mg/kg dose levels to milligrams for a 70 kg adult. Not retained as a",
        "within-human covariate: all human simulations used a single 70 kg subject",
        "weight, so no human weight effect is identified by this paper."
      )
    )
  )

  compartmentData <- list(
    depot     = list(analyte = "DNDI-6148", units = "mg",     specimen = "administration site", verified = TRUE),
    central   = list(analyte = "DNDI-6148", units = "mg",     specimen = "plasma",              verified = TRUE),
    skin      = list(analyte = "DNDI-6148", units = "ug/mL",  specimen = "tissue",              verified = TRUE),
    parasites = list(analyte = "Leishmania major bioluminescence", units = "photons/s", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species       = "human",
    n_subjects    = NA_integer_,
    n_studies     = 0L,
    weight_median = "70 kg (the single simulated adult body weight)",
    disease_state = "Cutaneous leishmaniasis (Leishmania major), simulated",
    dose_range    = paste(
      "Simulated only: 3-20 mg/kg once daily for 7-14 days in 1.0 mg/kg increments,",
      "and 1.5-7.0 mg/kg twice daily for 7-10 days in 0.5 mg/kg increments",
      "(Methods 2.2.6). 10,000 subjects per dose scenario."
    ),
    notes         = paste(
      "This is a SIMULATION model, not a fit: no human PK or PD data were analysed.",
      "The three typical PK values are single-species allometric projections from the",
      "murine estimates and every parameter is therefore fixed. The projections were",
      "checked by the authors against an independent Phase 1 study of DNDI-6148",
      "(10-380 mg single doses, n = 8 healthy males): predicted CL/F 3.44 L/h vs an",
      "observed 2.26-3.28 L/h, and predicted Vd/F 79.5 L vs an observed 53.2-104 L,",
      "all inside the 0.5-2-fold range (Discussion). Only three variance components",
      "were carried into the simulation and all were ASSUMED rather than estimated in",
      "humans; residual variability was excluded entirely (Methods 2.2.6)."
    )
  )

  ini({
    # =================================================================
    # Plasma PK. Single-species allometric projection from the murine fit,
    # Methods 2.2.5 and Results 3.2.6: fixed exponents 0.75 (CL/F), 1.00
    # (Vd/F) and -0.25 (ka) applied to the 22 g -> 70 kg size ratio of
    # 3181.8. Because the mouse and human unbound fractions are the same
    # 6.6% (Methods 2.2.3), the stated fu correction leaves the values
    # unchanged. Each is reproduced exactly by the printed murine estimate:
    #   CL = (65 / 8000) * 3181.8^0.75  = 3.44 L/h
    #   Vd =  0.025      * 3181.8^1     = 79.5 L
    #   ka =  2.7        * 3181.8^-0.25 = 0.360 /h
    # Saturable clearance is dropped, so CL for scaling was taken as
    # Vmax/F / Km (Discussion). Dose-dependent bioavailability is dropped
    # too; both were judged mouse-specific (Methods 2.2.5).
    # =================================================================
    lka <- fixed(log(0.360)) ; label("Oral absorption rate constant (ka, 1/h)")                        # Results 3.2.6: 0.360
    lcl <- fixed(log(3.44))  ; label("Apparent clearance for a 70 kg adult (CL/F, L/h)")               # Results 3.2.6: 3.44
    lvc <- fixed(log(79.5))  ; label("Apparent central volume for a 70 kg adult (Vd/F, L)")            # Results 3.2.6: 79.5

    # =================================================================
    # Skin distribution -- carried over unchanged from the murine model
    # ("Predicted human PK parameters were combined with tissue distribution
    # characteristics and PD components from the murine PK/PD model").
    # =================================================================
    lke0     <- fixed(log(20))   ; label("Plasma-to-tissue equilibration rate constant (k_plasma-tissue, 1/h)")  # Table 1: 20, Fixed; transferred from the murine model
    lkp_skin <- fixed(log(0.56)) ; label("Skin-to-plasma penetration coefficient (R_skin-plasma, unitless)")     # Table 1: 0.56 [0.49-0.68]; transferred from the murine model

    fu <- fixed(0.07) ; label("Fraction unbound in plasma and in skin (unitless; assumed equal in both matrices)")  # Methods 2.2.6: "An fu of 0.07 in plasma was applied, and a similar fu in tissue was assumed"

    # =================================================================
    # Parasite killing -- murine PD component transferred unchanged.
    # =================================================================
    lemax <- fixed(log(0.049)) ; label("Maximal rate of parasite elimination (Emax, 1/h)")                     # Table 2: 0.049 [0.042-0.058]; transferred from the murine model
    lec50 <- fixed(log(165))   ; label("Free skin concentration achieving half of Emax (fEC50, ug/L)")         # Table 2: 165 [125-236]; transferred from the murine model
    lhill <- fixed(log(1.6))   ; label("Sigmoidicity factor on free skin concentration (gamma, unitless)")     # Table 2: 1.6 [1.1-2.8]; transferred from the murine model

    lrbase_parasites <- fixed(log(10^8.01)) ; label("Baseline skin parasite bioluminescence (photons/s)")      # Methods 2.2.6: "Initial parasite loads in humans were set to the same as in mice"; Table 2 value 8.01 log10 photons/s

    # =================================================================
    # Between-subject variability. All three were ASSUMED for the human
    # simulation rather than estimated in humans (Methods 2.2.6). CV% is
    # converted with OMEGA = log(CV^2 + 1), the convention stated in the
    # Table 1 / Table 2 footnotes.
    # =================================================================
    etalcl              ~ fixed(0.086178)  # Methods 2.2.6, an inflated 30% CV on CL; log(0.30^2 + 1)
    etalkp_skin         ~ fixed(0.128309)  # Methods 2.2.6, 37% CV on R_skin-plasma; log(0.37^2 + 1)
    etalrbase_parasites ~ fixed(0.307485)  # Methods 2.2.6, 60% CV on parasite load at baseline; log(0.60^2 + 1)

    # Methods 2.2.6: "Residual variability was excluded." No residual-error
    # parameter is therefore defined -- the model is a typical-value plus BSV
    # simulation tool, and adding an invented residual SD would misstate it.
  })

  model({
    # ---- 1. Individual parameters ---------------------------------------
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)

    ke0     <- exp(lke0)
    kp_skin <- exp(lkp_skin + etalkp_skin)

    emax <- exp(lemax)
    hill <- exp(lhill)

    # fEC50 is published in ug/L (Table 2) and kept at that literal value in
    # ini(); units$concentration for this model is ug/mL, so convert here
    # rather than storing a pre-divided number that no longer matches the paper.
    ec50 <- exp(lec50) / 1000

    rbase_parasites <- exp(lrbase_parasites + etalrbase_parasites)

    # ---- 2. Plasma PK ----------------------------------------------------
    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc <- central / vc

    # ---- 3. Skin distribution (murine Equation 3) ------------------------
    d/dt(skin) <- ke0 * (kp_skin * Cc - skin)
    Cskin      <- skin

    # ---- 4. Parasite killing (murine Equations 4-5) ----------------------
    fCskin <- Cskin * fu
    kkill  <- emax * fCskin^hill / (ec50^hill + fCskin^hill)

    d/dt(parasites) <- -kkill * parasites
    parasites(0)    <-  rbase_parasites

    log10BLI <- log10(parasites)

    # The paper's efficacy endpoint: percent reduction in parasite load from
    # baseline, against which the 95% and 99% PTA targets are evaluated
    # (Figures 3-4, Table S1).
    parasiteReduction <- 100 * (1 - parasites / rbase_parasites)
  })
}
