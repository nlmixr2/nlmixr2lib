Cardilin_2019_radiation_radiosensitizer_mouse <- function() {
  # Preclinical (mouse FaDu xenograft) long-term tumor growth and kill model for
  # combination therapy with fractionated ionizing radiation and an anonymized
  # small-molecule radiosensitizer that interferes with DNA-damage repair.
  # Untreated growth is exponential in a proliferating compartment (V1) that
  # feeds a three-stage natural-death transit chain (V2-V4). Radiation has two
  # effects: (1) an instantaneous linear-quadratic (LQ) transfer of proliferating
  # cells into a lethally-irradiated pool that divides at most once more (U1 ->
  # U2) before dying, and (2) a permanent inhibition of the growth rate that
  # depends on the accumulated effective radiation dose IR_Tot. The
  # radiosensitizer's plasma concentration at the time of irradiation enhances
  # both effects.
  description <- paste(
    "Preclinical (mouse, FaDu head-and-neck xenograft).",
    "Long-term tumor growth, regrowth and eradication model for fractionated",
    "radiation combined with a radiosensitizer: exponential growth, a",
    "three-stage natural-death transit chain, linear-quadratic radiation kill",
    "into a one-more-division irradiated chain, and permanent growth",
    "inhibition by the accumulated radiation dose, both radiation effects",
    "enhanced by a one-compartment radiosensitizer PK."
  )
  reference <- paste(
    "Cardilin T, Almquist J, Jirstrand M, Zimmermann A, Lignet F, El Bawab S,",
    "Gabrielsson J. Modeling long-term tumor growth and kill after",
    "combinations of radiation and radiosensitizing agents. Cancer Chemother",
    "Pharmacol. 2019;83(6):1159-1173. doi:10.1007/s00280-019-03829-y.",
    "PMID 30976845; PMCID PMC6499765."
  )
  vignette <- "Cardilin_2019_radiation_radiosensitizer_mouse"
  # Tumor volume is in mm^3; radiosensitizer concentration (Cc) in ug/mL;
  # radiation dose in Gy; radiosensitizer dose in mg/kg.
  units <- list(time = "day", dosing = "mg/kg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "radiosensitizer", units = "mg/kg", specimen = "plasma", verified = TRUE),
    rad_depot = list(
      analyte = "radiation timing trigger",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    cycling_cells = list(analyte = "proliferating tumor cells", units = "mm^3", specimen = "tumor", verified = TRUE),
    irrad1 = list(
      analyte = "lethally irradiated cells (one division left)",
      units = "mm^3",
      specimen = "tumor",
      verified = TRUE
    ),
    irrad2 = list(analyte = "lethally irradiated daughter cells", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells1 = list(analyte = "dying tumor cells", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells2 = list(analyte = "dying tumor cells", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells3 = list(analyte = "dying tumor cells", units = "mm^3", specimen = "tumor", verified = TRUE),
    irtot = list(
      analyte = "accumulated effective radiation dose",
      units = "Gy",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  # No animal covariates: a fixed-design xenograft experiment. The radiation
  # dose per fraction is the parameter `radDose` (Gy), overridable per scenario.
  covariateData <- list()

  population <- list(
    species = "mouse (female CD1 nu/nu or NMRI nu/nu; FaDu head-and-neck squamous-cell-carcinoma xenograft)",
    n_subjects = 40,
    n_studies = 1,
    age_range = "5-6 weeks at purchase",
    sex_female_pct = 100,
    disease_state = "FaDu xenograft (human head and neck squamous cell carcinoma), 50-110 mm^3 at treatment start",
    dose_range = paste(
      "Radiation 2 Gy per fraction, alone or 30 min after an oral radiosensitizer",
      "dose of 25 or 100 mg/kg, 5 days a week (Mon-Fri) for 6 weeks (60 Gy total)"
    ),
    regions = "Preclinical (Merck KGaA, Darmstadt, Germany)",
    notes = paste(
      "Tumor dataset: 40 mice, 10 per arm -- vehicle, radiation 2 Gy, and",
      "radiation 2 Gy + radiosensitizer 25 or 100 mg/kg; calipers twice a",
      "week for up to 12 weeks after treatment stopped (Methods, Experimental",
      "data). Radiosensitizer PK from 8 further mice (4 at 25 and 4 at 100",
      "mg/kg PO; samples 1, 2 and 6 h post-dose). Mixed-effects FOCE fit in a",
      "Mathematica framework (Fraunhofer-Chalmers Centre)."
    )
  )

  ini({
    # --- System (tumor) parameters: Table 1 ---
    lkg <- log(0.40); label("Tumor growth rate kg (1/day)") # Table 1: kg = 0.40 day^-1 (RSE 4%)
    lkk <- log(0.26); label("Natural kill (transit) rate kk (1/day)") # Table 1: kk = 0.26 day^-1 (RSE 5%)
    lv0 <- log(27.0); label("Initial proliferating volume V1(0) = V0 (mm^3)") # Table 1: V0 = 27.0 mm^3 (RSE 7%)

    # --- Radiation parameters: Table 1 (printed in kGy^-1; converted to Gy^-1) ---
    lgammaRad <- log(0.0040); label("Long-term radiation growth-inhibition coefficient gamma (1/Gy)") # Table 1: gamma = 4.0 kGy^-1 (RSE 14%) = 0.0040 Gy^-1
    lalpha <- log(0.054); label("LQ linear radiation parameter alpha (1/Gy)") # Table 1: alpha = 54.0 kGy^-1 (RSE 18%) = 0.054 Gy^-1
    abratio <- fixed(10); label("LQ alpha/beta ratio (Gy)") # Methods, Computational methods: alpha/beta fixed to 10

    # --- Radiosensitizer potency parameters: Table 1 ---
    # The Table 1 letters are swapped relative to Eqs. 3 and 7. The value whose
    # effect sets the long-term (IR_Tot) enhancement is 0.42 mL/ug: it
    # reproduces the TSE curve of Fig. 6 (107.7 Gy at C = 0, ~58 Gy at C = 2,
    # ~15 Gy at C = 15 ug/mL) and the end-of-treatment net growth rates of
    # Fig. 5 (curve C ~ -0.006, curve D ~ -0.125 day^-1). The other value,
    # 0.15 mL/ug, is therefore the short-term (LQ kill) enhancement; it
    # reproduces the reported 25% kill per 2 Gy fraction at 100 mg/kg.
    lbRad <- log(0.15); label("Radiosensitizer enhancement of the LQ radiation kill (mL/ug)") # Table 1 row 'b' = 0.15 mL/ug (RSE 35%); Eq. 3 symbol a
    lbRadTot <- log(0.42); label("Radiosensitizer enhancement of the accumulated radiation dose (mL/ug)") # Table 1 row 'a' = 0.42 mL/ug (RSE 31%); Eq. 7 symbol b; Figs. 5-6

    # --- Radiosensitizer one-compartment PK (Eq. 1; summary values only) ---
    # The paper gives no PK parameter table. ke is from the reported 3 h
    # half-life; V/F is from the simulated peaks of Fig. 3 (2.09 and 8.39 ug/mL
    # after 25 and 100 mg/kg, digitised by the maintainers; the text rounds them
    # to 2 and 8 ug/mL). Held fixed in the tumor model (sequential PK-PD).
    lke <- fixed(log(log(2) / 3 * 24)); label("Radiosensitizer elimination rate ke (1/day)") # Results, Exposure: half-life 3 h -> ln(2)/3 h^-1 x 24
    lvf <- fixed(log(11.95)); label("Radiosensitizer apparent volume V/F (L/kg)") # Fig. 3 peaks 2.09 / 8.39 ug/mL after 25 / 100 mg/kg (digitised); text: 2 and 8 ug/mL

    # Radiation dose per fraction (Gy): an experimental input, not a fitted
    # parameter. Override per scenario (e.g. for TSE-curve exploration).
    radDose <- fixed(2); label("Radiation dose per fraction (Gy) - experimental input") # Methods: 2 Gy per dose

    # --- Between-subject variability (Table 1, CV%) ---
    # Variance = (CV/100)^2, the sqrt(omega^2) x 100 convention this group's
    # Mathematica mixed-effects framework uses (Cardilin 2018 Table 1 footnote).
    etalv0 ~ 0.25 # Table 1: V0 BSV 50% -> 0.50^2
    etalgammaRad ~ 0.2116 # Table 1: gamma BSV 46% -> 0.46^2
    etalalpha ~ 0.1764 # Table 1: alpha BSV 42% -> 0.42^2

    # --- Residual error (Table 1) ---
    propSd <- 0.24; label("Proportional residual error on tumor volume (fraction)") # Table 1: sigma_a = 24.0% (proportional)
    addSd <- 6.9; label("Additive residual error on tumor volume (mm^3)") # Table 1: sigma_b = 6.9 mm^3 (additive)
  })

  model({
    # 1. Individual parameters
    kg <- exp(lkg)
    kk <- exp(lkk)
    v0 <- exp(lv0 + etalv0)
    gammaRad <- exp(lgammaRad + etalgammaRad)
    alpha <- exp(lalpha + etalalpha)
    beta <- alpha / abratio
    bRad <- exp(lbRad)
    bRadTot <- exp(lbRadTot)
    ke <- exp(lke)
    vf <- exp(lvf)

    # 2. Instantaneous radiation effects at each fraction (Eqs. 2, 6, 7).
    # Each fraction is a unit bolus (amt = 1) into rad_depot, which decays at
    # krad so that integral(krad * rad_depot dt) = 1 per fraction. The kill
    # hazard then integrates to (1 + bRad*C)(alpha*D + beta*D^2), multiplying
    # V1 by 1 - F(D, C) (Eq. 3) and moving F(D, C) * V1 into U1; the IR_Tot
    # input integrates to (1 + bRadTot*C) * D (Eq. 7). krad is a solver device
    # (large => near-instantaneous), not a fitted parameter.
    krad <- 5000

    # 3. Radiosensitizer plasma concentration (Eq. 1)
    Cc <- central / vf

    lethal <- (1 + bRad * Cc) * (alpha * radDose + beta * radDose * radDose)
    killHaz <- lethal * krad * rad_depot

    # 4. Long-term growth inhibition by the accumulated dose (Eq. 4)
    inhib <- exp(-gammaRad * irtot)

    d/dt(central) <- -ke * central
    d/dt(rad_depot) <- -krad * rad_depot

    # Eq. 2 (V1), Eq. 6 (U1, U2), Eq. 5 (V2-V4), Eq. 7 (IR_Tot)
    d/dt(cycling_cells) <- kg * inhib * cycling_cells - kk * cycling_cells - killHaz * cycling_cells
    d/dt(irrad1) <- killHaz * cycling_cells - kg * irrad1 - kk * irrad1
    d/dt(irrad2) <- 2 * kg * irrad1 - kk * irrad2
    d/dt(damaged_cells1) <- kk * cycling_cells + kk * irrad1 + kk * irrad2 - kk * damaged_cells1
    d/dt(damaged_cells2) <- kk * damaged_cells1 - kk * damaged_cells2
    d/dt(damaged_cells3) <- kk * damaged_cells2 - kk * damaged_cells3
    d/dt(irtot) <- (1 + bRadTot * Cc) * radDose * krad * rad_depot

    # 5. Initial conditions (Eq. 8): V_i(0) = V0 (kk/kg)^(i-1), U_j(0) = 0,
    # IR_Tot(0) = 0 -- the pseudo-steady split that makes untreated growth
    # strictly exponential at kg - kk.
    cycling_cells(0) <- v0
    damaged_cells1(0) <- v0 * (kk / kg)
    damaged_cells2(0) <- v0 * (kk / kg)^2
    damaged_cells3(0) <- v0 * (kk / kg)^3

    # Net growth rate of the proliferating pool (Eq. 9; plotted in Fig. 5)
    knet <- kg * inhib - kk

    # 6. Observation: total tumor volume
    tumor_vol <- cycling_cells + irrad1 + irrad2 + damaged_cells1 + damaged_cells2 + damaged_cells3
    tumor_vol ~ add(addSd) + prop(propSd)
  })
}
attr(Cardilin_2019_radiation_radiosensitizer_mouse, "message") <-
  "Radiation is given as a unit bolus (amt=1) into the rad_depot compartment at each irradiation time; the per-fraction radiation dose (Gy) is the parameter radDose (default 2). The radiosensitizer is dosed (mg/kg) into the central compartment. Observation tumor_vol is total tumor volume (mm^3); knet is the net growth rate of Eq. 9 (1/day)."
Cardilin_2019_radiation_radiosensitizer_mouse
