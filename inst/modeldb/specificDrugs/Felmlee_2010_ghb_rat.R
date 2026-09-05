Felmlee_2010_ghb_rat <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Mechanistic hybrid-physiological toxicokinetic model for",
    "gamma-hydroxybutyric acid (GHB) in male Sprague-Dawley rats after",
    "intravenous bolus doses of 200 to 1,000 mg/kg, fitting plasma",
    "concentrations and cumulative urinary excretion simultaneously.",
    "Disposition is a plasma compartment exchanging with two tissue",
    "compartments (fast and slow) by distributional clearances, plus a",
    "kidney compartment perfused at renal blood flow. Capacity-limited",
    "metabolic elimination (Michaelis-Menten) acts on the first tissue",
    "compartment. Drug is filtered from the kidney compartment at the",
    "glomerular filtration rate into a proximal-tubule ultrafiltrate",
    "compartment, from which it is actively reabsorbed back into the",
    "kidney by a saturable (Michaelis-Menten) monocarboxylate-transporter",
    "process; the surviving filtrate flows at the urine flow rate through",
    "a distal-tubule ultrafiltrate compartment into the cumulative urine",
    "compartment. Plasma volume, kidney volume, both ultrafiltrate",
    "volumes, renal blood flow and urine flow were fixed to rat",
    "physiological values. The saturable reabsorption term carries an",
    "optional monocarboxylate-transporter inhibition extension (source",
    "paper Eqs. 11-13) that reproduces the paper's competitive,",
    "noncompetitive and uncompetitive inhibition simulations; with the",
    "two inhibition covariates at their default of zero the term reduces",
    "exactly to the fitted model (source paper Eq. 4). Fit by NONMEM VI",
    "ADVAN9 with FOCE."
  )
  reference <- paste(
    "Felmlee MA, Wang Q, Cui D, Roiko SA, Morris ME.",
    "Mechanistic toxicokinetic model for gamma-hydroxybutyric acid:",
    "inhibition of active renal reabsorption as a potential therapeutic",
    "strategy.",
    "The AAPS Journal. 2010;12(3):407-416.",
    "doi:10.1208/s12248-010-9197-x.",
    sep = " "
  )
  vignette <- "Felmlee_2010_ghb_rat"

  # Aulf1 / Aulf2 of source paper Eqs. 4-8: the ultrafiltrate contained in
  # the proximal and the distal renal tubule respectively. Active
  # (transporter-mediated) reabsorption occurs only from the proximal one;
  # the distal one exists to reproduce the delay between glomerular
  # filtration and the appearance of drug in voided urine (source paper
  # Discussion, "The two ultrafiltrate compartments in the model represent
  # the proximal and distal tubules of the kidney").
  paper_specific_compartments <- c("ulf1", "ulf2")

  units <- list(time = "min", dosing = "mg", concentration = "mg/mL")

  compartmentData <- list(
    central     = list(analyte = "GHB", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "GHB", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "GHB", units = "mg", specimen = "tissue", verified = TRUE),
    kidney      = list(analyte = "GHB", units = "mg", specimen = "tissue", verified = TRUE),
    ulf1        = list(analyte = "GHB", units = "mg", specimen = "urine", verified = TRUE),
    ulf2        = list(analyte = "GHB", units = "mg", specimen = "urine", verified = TRUE),
    urine       = list(analyte = "GHB", units = "mg", specimen = "urine", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only Table I parameter reported per kilogram is the glomerular",
        "filtration rate (10 mL/min/kg), so body weight enters the model",
        "only through gfr = exp(lgfr) * WT. Every other volume and flow in",
        "Table I is an absolute per-rat value and is NOT weight-scaled. The",
        "study rats weighed 280-320 g (source paper Methods, Animals and",
        "Surgery); 0.3 kg is the value at which the model reproduces the",
        "published Table II simulations."
      ),
      source_name        = "WT"
    ),
    INH_MCT_KM_RATIO = list(
      description        = "MCT-inhibitor concentration ratio [I]/Ki acting on the reabsorption Michaelis-Menten constant",
      units              = "unitless",
      type               = "continuous",
      reference_category = "0 = no inhibitor present",
      notes              = paste(
        "Source paper Eqs. 11-13 modify the Eq. 4 reabsorption term for a",
        "steady-state transport inhibitor with R = [I]/Ki. Set together with",
        "INH_MCT_CONC_RATIO to select the mechanism:",
        "competitive (Eq. 11) INH_MCT_KM_RATIO = R, INH_MCT_CONC_RATIO = 0;",
        "noncompetitive (Eq. 12) both = R;",
        "uncompetitive (Eq. 13) INH_MCT_KM_RATIO = 0, INH_MCT_CONC_RATIO = R.",
        "Both zero recovers the fitted model, Eq. 4."
      ),
      source_name        = "R"
    ),
    INH_MCT_CONC_RATIO = list(
      description        = "MCT-inhibitor concentration ratio [I]/Ki acting on the proximal-tubule substrate concentration",
      units              = "unitless",
      type               = "continuous",
      reference_category = "0 = no inhibitor present",
      notes              = paste(
        "Companion to INH_MCT_KM_RATIO; see that entry for the",
        "mechanism-selection table (source paper Eqs. 11-13)."
      ),
      source_name        = "R"
    )
  )

  population <- list(
    species        = "rat (male Sprague-Dawley)",
    n_subjects     = 34L,
    n_studies      = 1L,
    age_range      = "adult (exact age not reported)",
    weight_range   = "280-320 g",
    sex_female_pct = 0,
    disease_state  = paste(
      "Healthy male Sprague-Dawley rats (Harlan, Indianapolis, IN) with",
      "cannulas implanted in the right jugular vein, housed individually",
      "for 2-3 days after surgery and then in metabolic cages for the",
      "duration of the study."
    ),
    dose_range     = paste(
      "GHB 200, 400, 600 or 1,000 mg/kg as a single intravenous bolus into",
      "the jugular vein cannula, N = 7-10 rats per dose group. Plasma",
      "sampled at 0, 5, 10, 20, 30, 60, 90, 120, 180, 240, 300 and 360 min;",
      "urine collected over 0-60, 60-120, 120-240 and 240-360 min."
    ),
    regions        = "USA (University at Buffalo)",
    notes          = paste(
      "n_subjects is the toxicokinetic dataset used to fit the model (four",
      "dose groups of 7-10 rats each; the source paper reports the range",
      "rather than the exact per-group counts, so 34 is the midpoint of the",
      "28-40 implied range). A separate GHB / L-lactate interaction study",
      "(GHB 600 mg/kg alone, N = 3, versus GHB 600 mg/kg plus L-lactate",
      "330 mg/kg IV bolus and 121 mg/kg/h IV infusion, N = 5) was used to",
      "confirm the inhibition simulations and did not contribute to the",
      "model fit. Animal protocols approved by the University at Buffalo",
      "Institutional Animal Care and Use Committee. GHB in plasma and urine",
      "measured by LC/MS/MS. Model fit in NONMEM VI (level 1.1) with FOCE",
      "and ADVAN9, plasma and urine data from all doses fitted",
      "simultaneously. See source paper Methods and Table I."
    )
  )

  ini({
    # ================================================================
    # Volumes (source paper Table I). Table I footnote a marks the
    # parameters that were fixed to rat physiological values: Vplasma,
    # Vkidney, Vulf1, Vulf2, QR and UF.
    # ================================================================
    lvc       <- fixed(log(10.5)) ; label("Plasma volume Vplasma (mL)")                          # Table I: Vplasma = 10.5 mL, footnote a (fixed to physiological value)
    lvp       <- log(75.9)        ; label("Volume of the first (fast) tissue compartment Vtissue1 (mL)")   # Table I: Vtissue1 = 75.9 mL
    lvp2      <- log(26.8)        ; label("Volume of the second (slow) tissue compartment Vtissue2 (mL)")  # Table I: Vtissue2 = 26.8 mL
    lv_kidney <- fixed(log(4.0))  ; label("Kidney volume Vkidney (mL)")                          # Table I: Vkidney = 4.0 mL, footnote a
    lv_ulf1   <- fixed(log(3.0))  ; label("Volume of the proximal-tubule ultrafiltrate compartment Vulf1 (mL)")  # Table I: Vulf1 = 3.0 mL, footnote a
    lv_ulf2   <- fixed(log(1.0))  ; label("Volume of the distal-tubule ultrafiltrate compartment Vulf2 (mL)")    # Table I: Vulf2 = 1.0 mL, footnote a

    # ================================================================
    # Flows and distributional clearances (source paper Table I).
    # ================================================================
    lq        <- log(26.9)        ; label("Distributional clearance to the first tissue compartment CLD (mL/min)")   # Table I: CLD = 26.9 mL/min
    lq2       <- log(3.07)        ; label("Distributional clearance to the second tissue compartment CLD2 (mL/min)") # Table I: CLD2 = 3.07 mL/min
    lq_kidney <- fixed(log(12.5)) ; label("Renal blood flow QR (mL/min)")                        # Table I: QR = 12.5 mL/min, footnote a
    # GFR is encoded as fixed although Table I does NOT mark it with footnote a
    # (see vignette Errata). Two Table I features say it is a looked-up
    # physiological constant rather than an estimate: it is the only row whose
    # units are per kilogram, and it is the only non-error row printed as a bare
    # integer while every estimated parameter carries three significant figures.
    # The paper also states it fixed parameters to physiological values citing
    # refs 27-29, one of which is Davies and Morris 1993, the standard
    # laboratory-animal physiology compendium.
    lgfr      <- fixed(log(10))   ; label("Glomerular filtration rate GFR (mL/min/kg)")          # Table I: GFR = 10; Table I prints the unit as mg/min/kg, a typo for mL/min/kg (see vignette Errata)
    luf       <- fixed(log(0.1))  ; label("Urine flow UF (mL/min)")                              # Table I: UF = 0.1 mL/min, footnote a

    # ================================================================
    # Capacity-limited metabolic elimination from the first tissue
    # compartment (source paper Eq. 2 and Table I).
    # ================================================================
    lvmax     <- log(0.581)       ; label("Maximum metabolic elimination rate Vmax,m (mg/min)")  # Table I: Vmax,m = 0.581 mg/min
    lkm       <- log(0.054)       ; label("Metabolic Michaelis-Menten constant Km,m (mg/mL)")    # Table I: Km,m = 0.054 mg/mL

    # ================================================================
    # Saturable active renal reabsorption from the proximal-tubule
    # ultrafiltrate compartment (source paper Eq. 4 and Table I).
    # ================================================================
    lvmax_reab <- log(2.34)       ; label("Maximum renal reabsorption rate Vmax,R (mg/min)")            # Table I: Vmax,R = 2.34 mg/min
    lkm_reab   <- log(0.46)       ; label("Renal reabsorption Michaelis-Menten constant Km,R (mg/mL)")  # Table I: Km,R = 0.46 mg/mL

    # ================================================================
    # Between-subject variability. Source paper: "an exponential model
    # was used to describe between-subject variability, Pi = theta *
    # exp(eta_i)", with Table I reporting CV% for between-subject
    # variability on exactly two parameters, Vtissue1 (16%) and UF
    # (114%). omega^2 = log(CV^2 + 1) for the exponential model, the
    # same convention used for the companion rat GHB model
    # Fung_2008_butanediol_rat.
    # ================================================================
    etalvp ~ 0.025278             # Table I: Vtissue1 CV 16 percent -> log(1 + 0.16^2) = 0.025278
    etaluf ~ 0.832909             # Table I: UF CV 114 percent -> log(1 + 1.14^2) = 0.832909

    # ================================================================
    # Residual unexplained variability. Source paper Eqs. 9-10 use a
    # proportional log error model, Y = Log(prediction) + epsilon, i.e.
    # additive on the natural-log scale, which is the lognormal
    # residual-error model in nlmixr2.
    # ================================================================
    expSd       <- 0.025          ; label("Plasma proportional (log-scale additive) residual error (fraction)")            # Table I: epsilon_plasma = 2.5 percent
    expSd_urine <- 0.46           ; label("Urinary-excretion proportional (log-scale additive) residual error (fraction)") # Table I: epsilon_urine = 46 percent
  })

  model({
    # ================================================================
    # 1. Individual parameters
    # ================================================================
    vc       <- exp(lvc)
    vp       <- exp(lvp + etalvp)
    vp2      <- exp(lvp2)
    v_kidney <- exp(lv_kidney)
    v_ulf1   <- exp(lv_ulf1)
    v_ulf2   <- exp(lv_ulf2)

    q        <- exp(lq)
    q2       <- exp(lq2)
    q_kidney <- exp(lq_kidney)
    uf       <- exp(luf + etaluf)

    # Table I reports GFR per kilogram of body weight; every other Table I
    # flow and volume is an absolute per-rat value.
    gfr      <- exp(lgfr) * WT

    vmax      <- exp(lvmax)
    km        <- exp(lkm)
    vmax_reab <- exp(lvmax_reab)
    km_reab   <- exp(lkm_reab)

    # ================================================================
    # 2. Compartment concentrations
    # ================================================================
    Cc    <- central / vc
    Cp1   <- peripheral1 / vp
    Cp2   <- peripheral2 / vp2
    Ck    <- kidney / v_kidney
    Culf1 <- ulf1 / v_ulf1
    Culf2 <- ulf2 / v_ulf2

    # ================================================================
    # 3. Saturable active renal reabsorption (source paper Eq. 4),
    # generalised to the transport-inhibition forms of Eqs. 11-13.
    # "Reabsorption" has units of clearance (mL/min); the reabsorbed
    # mass flux is Reabsorption * Culf1, which is the Michaelis-Menten
    # rate Vmax,R * Culf1 / (Km,R + Culf1).
    #
    #   INH_MCT_KM_RATIO = INH_MCT_CONC_RATIO = 0  -> Eq. 4  (fitted model)
    #   INH_MCT_KM_RATIO = R, INH_MCT_CONC_RATIO = 0 -> Eq. 11 (competitive)
    #   INH_MCT_KM_RATIO = INH_MCT_CONC_RATIO = R  -> Eq. 12 (noncompetitive)
    #   INH_MCT_KM_RATIO = 0, INH_MCT_CONC_RATIO = R -> Eq. 13 (uncompetitive)
    # ================================================================
    reabsorption <- vmax_reab /
      (km_reab * (1 + INH_MCT_KM_RATIO) + Culf1 * (1 + INH_MCT_CONC_RATIO))

    # ================================================================
    # 4. ODE system (source paper Eqs. 1-3 and 5-8). All compartments
    # start empty; the IV bolus enters `central`.
    # ================================================================
    d/dt(central)     <- -(q_kidney + q + q2) * Cc + q_kidney * Ck + q * Cp1 + q2 * Cp2  # Eq. 1
    d/dt(peripheral1) <- q * Cc - (q + vmax / (km + Cp1)) * Cp1                          # Eq. 2
    d/dt(peripheral2) <- q2 * Cc - q2 * Cp2                                              # Eq. 3
    d/dt(kidney)      <- q_kidney * Cc - (q_kidney + gfr) * Ck + reabsorption * Culf1     # Eq. 5
    d/dt(ulf1)        <- gfr * Ck - uf * Culf1 - reabsorption * Culf1                     # Eq. 6
    d/dt(ulf2)        <- uf * Culf1 - uf * Culf2                                          # Eq. 7
    d/dt(urine)       <- uf * Culf2                                                       # Eq. 8

    # ================================================================
    # 5. Observations. Two simultaneously fitted outputs: the plasma
    # concentration (source paper Eq. 9) and the cumulative amount
    # excreted in urine (source paper Eq. 10). Both use a proportional
    # log error model, which is nlmixr2's lognormal residual error.
    # ================================================================
    Cc ~ lnorm(expSd)
    urine ~ lnorm(expSd_urine)
  })
}
