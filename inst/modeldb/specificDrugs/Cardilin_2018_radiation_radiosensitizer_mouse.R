Cardilin_2018_radiation_radiosensitizer_mouse <- function() {
  # Preclinical (mouse FaDu xenograft) tumor-growth-inhibition (TGI) model for
  # combination therapy with ionizing radiation and a radiosensitizing compound
  # (anonymized "RS1" in the source). Untreated growth is logistic in a
  # proliferating compartment with a three-stage natural-death transit chain.
  # Radiation acts through the linear-quadratic (LQ) surviving-fraction theory as
  # an (approximately) instantaneous mass transfer of proliferating cells into an
  # irradiated chain that divides at most once more before dying. The
  # radiosensitizer (one-compartment oral PK) stimulates both the natural-death
  # process and the radiation-induced kill (synergy).
  description <- paste(
    "Preclinical (mouse, FaDu head-and-neck xenograft).",
    "Tumor growth inhibition model for combination therapy with ionizing",
    "radiation and a radiosensitizer (linear-quadratic radiation kill with a",
    "damage-compartment transit chain, driven by a one-compartment",
    "radiosensitizer PK)."
  )
  reference <- paste(
    "Cardilin T, Almquist J, Jirstrand M, Zimmermann A, El Bawab S,",
    "Gabrielsson J. Model-Based Evaluation of Radiation and Radiosensitizing",
    "Agents in Oncology. CPT Pharmacometrics Syst Pharmacol. 2018;7(1):51-58.",
    "doi:10.1002/psp4.12268. PMID 29218836; PMCID PMC5784742."
  )
  vignette <- "Cardilin_2018_radiation_radiosensitizer_mouse"
  # Tumor volume is in mm^3; radiosensitizer concentration (Cc) is in ug/mL;
  # radiation dose is in Gy; RS1 dose is in mg/kg.
  units <- list(time = "day", dosing = "mg/kg", concentration = "ug/mL")

  # No patient covariates: this is a fixed-design xenograft experiment. The
  # radiation dose per fraction is encoded as the parameter `radDose` (Gy) so it
  # can be overridden per scenario / subject without a covariate column.
  covariateData <- list()

  population <- list(
    species       = "mouse (female CD1 nu/nu or NMRI nu/nu; FaDu head-and-neck squamous-cell-carcinoma xenograft)",
    n_subjects    = 80,
    n_studies     = 1,
    disease_state = "FaDu xenograft (human head and neck squamous cell carcinoma)",
    dose_range    = "RS1 10-200 mg/kg PO once daily on days 3-7; radiation 2 Gy/fraction on days 3-7",
    regions       = "Preclinical (Merck, Darmstadt, Germany)",
    notes         = paste(
      "Primary tumor dataset: 80 mice in 8 groups of 10 -- (A) vehicle,",
      "(B) radiation 2 Gy, (C-E) RS1 10/50/200 mg/kg, (F-H) combination of",
      "2 Gy + RS1 10/50/200 mg/kg, dosed once daily on days 3-7",
      "(Supplementary Material B). RS1 was given orally 10 min before",
      "irradiation. A separate 10-mouse dataset (2 Gy + RS1 25 mg/kg, 5 days/week",
      "for 6 weeks) was used for prediction (Figure 5). Radiosensitizer PK was",
      "estimated from 147 tumor-bearing female mice across 5 studies",
      "(Supplementary Table SI). Between-subject variability was estimated for",
      "initial tumor volume, tumor capacity, and the two PK parameters."
    )
  )

  ini({
    # --- Tumor model (Table 1 of the reference) ---
    lkg    <- log(0.50)   ; label("Tumor natural growth rate kg (1/day)")              # Table 1
    lkk    <- log(0.28)   ; label("Tumor natural kill rate kk (1/day)")                # Table 1
    lcap   <- log(2200)   ; label("Tumor capacity K (mm^3)")                           # Table 1
    lv0    <- log(40.0)   ; label("Initial volume of the proliferating compartment V0 (mm^3)")  # Table 1

    # Linear-quadratic radiation parameters. alpha/beta was FIXED to 10
    # (Methods; Table 1 footnote), so beta is derived as alpha/abratio.
    lalpha  <- log(0.08)  ; label("LQ linear radiation parameter alpha (1/Gy)")        # Table 1
    abratio <- fixed(10)  ; label("LQ alpha/beta ratio (Gy) - fixed")                  # Methods (fixed at 10)

    # Radiosensitizer pharmacodynamic coefficients. The Table 1 letters and the
    # Mathematica supplement code-variable names for these two coefficients are
    # swapped; the placement here follows the paper's quantitative text and is
    # confirmed by the reported combination kill fractions (2 Gy + 10/50/200
    # mg/kg => 24/46/85% of proliferating cells killed per fraction):
    #   * natural-death stimulation = 0.09 mL/ug ("9% increase per ug/mL"),
    #   * radiation synergy          = 0.45 mL/ug ("45% increase per ug/mL").
    laDeath <- log(0.09)  ; label("RS1 stimulation of natural cell death (mL/ug)")     # Table 1 (paper symbol a)
    lbRad   <- log(0.45)  ; label("RS1 stimulation of radiation-induced kill (mL/ug)") # Table 1 (paper symbol b)

    # --- Radiosensitizer (RS1) one-compartment oral PK (Supplementary Table SI) ---
    # ke reported as 0.21 1/h; converted to 1/day (x 24) to match the day-based
    # tumor model. V reported as 0.0095 (apparent V/F); the printed unit label
    # ("mL/kg") is a mis-transcription -- on a ng/mL concentration basis the
    # value is 0.009507 (mg/kg)/(ng/mL) = 9.507 L/kg. This reproduces the paper's
    # reported average plasma concentrations (25 mg/kg -> 0.52 ug/mL,
    # 48 mg/kg -> 0.9 ug/mL, 400 mg/kg -> 8.2 ug/mL) and the combination kill
    # fractions above.
    lke    <- log(0.21011834 * 24) ; label("RS1 elimination rate ke (1/day)")          # Suppl. Table SI (0.21/h x 24)
    lvf    <- log(9.506615)        ; label("RS1 apparent volume of distribution V/F (L/kg)")  # Suppl. Table SI

    # Radiation dose per fraction (Gy). An experimental input, not a fitted
    # parameter: fixed to the primary-experiment value of 2 Gy and overridable
    # per scenario / subject (e.g. the TSE-curve doses 0.67 / 2.17 Gy).
    radDose <- fixed(2) ; label("Radiation dose per fraction (Gy) - experimental input")  # Methods (2 Gy)

    # --- Between-subject variability ---
    # Reported as BSV% = sqrt(omega^2) x 100 (Table 1 / Table SI footnotes), so
    # the log-scale variance is (BSV%/100)^2.
    etalv0  ~ 0.0441   # V0  BSV 21% (Table 1):    0.21^2
    etalcap ~ 0.1024   # K   BSV 32% (Table 1):    0.32^2
    etalvf  ~ 0.7744   # V/F BSV 88% (Table SI):   0.88^2
    etalke  ~ 0.0225   # ke  BSV 15% (Table SI):   0.15^2

    # --- Residual error ---
    propSd <- 0.23 ; label("Proportional residual error on tumor volume (fraction)")  # Table 1 (sigma 23%)
  })

  model({
    # 1. Individual parameters
    kg     <- exp(lkg)
    kk     <- exp(lkk)
    cap    <- exp(lcap + etalcap)
    v0     <- exp(lv0 + etalv0)
    alpha  <- exp(lalpha)
    beta   <- alpha / abratio
    aDeath <- exp(laDeath)
    bRad   <- exp(lbRad)
    ke     <- exp(lke + etalke)
    vf     <- exp(lvf + etalvf)

    # 2. Numerical device for the (Dirac-delta) instantaneous radiation kill.
    # Radiation is delivered as a unit bolus (amt = 1) into `radDepot`, which
    # decays at the fast rate krad so that integral(krad * radDepot dt) = 1 per
    # fraction. The kill hazard below therefore integrates to the LQ lethal-lesion
    # number over each fraction, multiplying the proliferating pool by the
    # surviving fraction exp(-(1 + bRad*Cc)(alpha*D + beta*D^2)). krad is a
    # solver convenience (large => near-instantaneous), NOT a fitted parameter.
    krad <- 500

    # 3. Radiosensitizer plasma concentration (one-compartment, oral bolus input)
    Cc <- central / vf                       # ug/mL

    # 4. Tumor dynamics
    Vtot <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3 + irrad1 + irrad2
    logi <- 1 - Vtot / cap                    # logistic capacity factor

    lethal  <- (1 + bRad * Cc) * (alpha * radDose + beta * radDose * radDose)
    killHaz <- lethal * krad * radDepot       # near-instantaneous LQ kill at each fraction

    d/dt(central)  <- -ke * central           # RS1 PK (mg/kg)
    d/dt(radDepot) <- -krad * radDepot        # radiation timing trigger (unit area / fraction)

    d/dt(cycling_cells)  <- kg * cycling_cells * logi - (1 + aDeath * Cc) * kk * cycling_cells - killHaz * cycling_cells
    d/dt(damaged_cells1) <- (1 + aDeath * Cc) * kk * cycling_cells - kk * damaged_cells1
    d/dt(damaged_cells2) <- kk * damaged_cells1 - kk * damaged_cells2
    d/dt(damaged_cells3) <- kk * damaged_cells2 - kk * damaged_cells3
    d/dt(irrad1)  <- killHaz * cycling_cells - kk * irrad1 - kg * irrad1 * logi
    d/dt(irrad2)  <- 2 * kg * irrad1 * logi - kk * irrad2   # factor 2: one more division before death

    # 5. Initial conditions (Eq. 7; Mathematica supplement). The damage chain
    # starts at the pseudo-steady-state ratios (kk/kg)^n relative to the
    # proliferating compartment so that an untreated tumor grows exponentially at
    # rate kg - kk.
    cycling_cells(0)  <- v0
    damaged_cells1(0) <- v0 * (kk / kg)
    damaged_cells2(0) <- v0 * (kk / kg)^2
    damaged_cells3(0) <- v0 * (kk / kg)^3

    # 6. Observation: total tumor volume (sum of all six compartments)
    tumor_vol <- Vtot
    tumor_vol ~ prop(propSd)
  })
}
attr(Cardilin_2018_radiation_radiosensitizer_mouse, "message") <-
  "Radiation is given as a unit bolus (amt=1) into the radDepot compartment at each irradiation time; the per-fraction radiation dose (Gy) is the parameter radDose (default 2). RS1 is dosed (mg/kg) into the central compartment. Observation tumor_vol is total tumor volume (mm^3)."
Cardilin_2018_radiation_radiosensitizer_mouse
