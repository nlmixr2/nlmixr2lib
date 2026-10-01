Kosinsky_2022_camostat_pkpd <- function() {
  description <- paste(
    "Semi-mechanistic (QSP-style) PK/PD model for oral camostat mesylate",
    "and its active metabolite FOY-251 against SARS-CoV-2 (Kosinsky 2022).",
    "Camostat is hydrolysed too rapidly to quantify, so the modelled PK",
    "analyte is FOY-251: a one-compartment disposition (first-order",
    "absorption, linear elimination) fit to digitised human IV-infusion",
    "data (Midgley 1994) and oral data (FOIPAN package insert); the",
    "camostat dose is treated as equimolar to FOY-251. FOY-251 plasma",
    "concentration drives reversible covalent inhibition of the host",
    "serine protease TMPRSS2: a two-state target-turnover model in which",
    "active enzyme (target, relative activity 1 at baseline) is",
    "covalently bound by FOY-251 (rate kcat, half-saturation ki) into an",
    "inactive complex that slowly recovers (kdis, 14 h complex half-life),",
    "with enzyme synthesis/degradation (kdeg, assumed 12 h enzyme",
    "half-life) maintaining the baseline. Remaining TMPRSS2 activity is",
    "linked to SARS-2-S-driven viral entry rate by an empirical",
    "Hill-Langmuir relationship (ki50 = 0.047 activity at half-maximal",
    "entry, Hill 0.59). The PD parameters are carried from the in-vitro",
    "fit of Kosinsky_2022_camostat_invitro; only the PK parameters are",
    "estimated here. No between-subject variability was reported (the",
    "model was fit to digitised literature profiles and used for",
    "deterministic dose simulations)."
  )
  reference <- paste(
    "Kosinsky Y, Peskov K, Stanski DR, Wetmore D, Vinetz J.",
    "Semi-Mechanistic Pharmacokinetic-Pharmacodynamic Model of Camostat",
    "Mesylate-Predicted Efficacy against SARS-CoV-2 in COVID-19.",
    "Microbiol Spectr. 2022;10(2):e02167-21.",
    "doi:10.1128/spectrum.02167-21.",
    "PK parameters in Table 1; PD parameters in Table 2; the in-vivo",
    "PK/PD ODEs are Equations 1 and 5 of Materials and Methods.",
    sep = " "
  )
  vignette <- "Kosinsky_2022_camostat"
  units <- list(
    time = "h",
    dosing = "nmol",
    concentration = "nM"
  )

  # The modelled dose is FOY-251 in nmol (camostat mesylate is considered
  # equimolar to FOY-251, Methods). Camostat mesylate (salt) MW = 494.52
  # g/mol, so a clinical mg dose converts as nmol = mg / 494.52 * 1e6
  # (e.g. 200 mg = 404431 nmol). The MW is a standard chemical constant
  # (PubChem CID 5284360), documented in the vignette; the paper reports
  # only FOY-251 MW = 313 g/mol (Methods).
  compartmentData <- list(
    depot = list(
      analyte = "FOY-251 (camostat-equimolar)",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "FOY-251 (active camostat metabolite)",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    target = list(
      analyte = "TMPRSS2 (active host serine protease)",
      units = "fraction of baseline activity",
      specimen = "not applicable",
      verified = TRUE
    ),
    complex = list(
      analyte = "FOY-251-TMPRSS2 covalent complex (inactive enzyme)",
      units = "fraction of baseline activity",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    disease_state = "healthy adults (PK); model simulations target COVID-19",
    dose_range = "camostat mesylate 200-600 mg orally q6-q8h (simulated)",
    notes = paste(
      "The PK model was fit to digitised FOY-251 plasma profiles: a 12-h",
      "IV infusion of 40 mg camostat (Midgley et al. 1994, Xenobiotica",
      "24:79-92) and a single 200 mg oral dose (FOIPAN package insert);",
      "see Results 'Camostat mesylate pharmacokinetic model in humans' and",
      "Figure S1. The PD parameters were estimated in vitro (Hoffmann et",
      "al. 2020 data) and reused here. No subject-level cohort or",
      "between-subject variability was reported; the model is used for",
      "deterministic dose-regimen simulations (Table 3, Figure 3)."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # FOY-251 one-compartment PK (Table 1; all estimated).
    # ODEs are Equation 1 of Methods: dAd/dt = -ka*Ad;
    # dAc/dt = ka*Ad*F - kel*Ac; Cc (nM) = Ac/Vd.
    # ---------------------------------------------------------------------
    lka <- log(0.67)
    label("FOY-251 first-order absorption rate ka (1/h)") # Table 1 ka = 0.67 (RSE 5.78%)
    lfdepot <- log(0.051)
    label("Camostat oral bioavailability F (fraction)") # Table 1 Fbio = 0.051 (RSE 4.66%)
    lvc <- log(22.36)
    label("FOY-251 volume of distribution Vd (L)") # Table 1 Vd = 22.36 L (RSE 7.07%)
    lkel <- log(1.22)
    label("FOY-251 linear elimination rate kel (1/h)") # Table 1 kel = 1.22 (RSE 4.70%)

    # ---------------------------------------------------------------------
    # TMPRSS2 covalent-inhibition turnover PD, carried (fixed) from the
    # in-vitro fit (Kosinsky_2022_camostat_invitro; Table 2) plus the
    # assumed enzyme turnover. In-vivo ODEs are Equation 5 of Methods.
    # ---------------------------------------------------------------------
    lki <- fixed(log(45638.51))
    label("FOY-251 TMPRSS2 half-saturation / inhibition constant ki (nM)") # Table 2 Ki = 45638.51 nM (= 45.6 uM; RSE 7.46%)
    kcat <- fixed(400)
    label("FOY-251 covalent-binding catalytic rate kcat (1/h)") # Table 2 kcat = 400 (fixed)
    kdis <- fixed(0.049)
    label("TMPRSS2 activity recovery rate kdis from covalent complex (1/h)") # Table 2 kdis = 0.049 (fixed; 14 h complex half-life)
    kdeg <- fixed(0.0575)
    label("TMPRSS2 degradation rate kdeg (1/h)") # Methods kdeg = log(2)/12 h = 0.0575 (assumed 12 h enzyme half-life)
    lki50 <- fixed(log(0.047))
    label("TMPRSS2 activity at half-maximal viral entry ki50 (fraction)") # Table 2 Ksp = 0.047 (RSE 42.4%)
    hill <- fixed(0.59)
    label("Hill coefficient of viral entry vs TMPRSS2 activity (unitless)") # Table 2 h = 0.59 (RSE 16.6%)

    # ---------------------------------------------------------------------
    # Residual error (PK only; Table 1). The in-vivo TMPRSS2 and viral
    # entry outputs are deterministic predictions and carry no residual
    # error in the paper.
    # ---------------------------------------------------------------------
    propSd <- 0.14
    label("Proportional residual error on FOY-251 concentration (fraction)") # Table 1 proportional residual error = 0.14 (RSE 17.3%)
  })

  model({
    # 1-2. Individual PK parameters (no IIV reported).
    ka <- exp(lka)
    f <- exp(lfdepot)
    vc <- exp(lvc)
    kel <- exp(lkel)

    # PD constants.
    ki <- exp(lki)
    ksp <- exp(lki50)

    # 4. PK ODE system (Equation 1). FOY-251 amount in nmol.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- f

    # FOY-251 plasma concentration in nM (Methods: Cc = Ac/Vd).
    Cc <- central / vc

    # TMPRSS2 covalent inhibition with turnover (Equation 5). target and
    # complex are relative to baseline enzyme activity (target(0) = 1).
    bind <- kcat * target * Cc / (Cc + ki)
    d/dt(target) <- kdeg * (1 - target) - bind + kdis * complex
    d/dt(complex) <- -kdeg * complex + bind - kdis * complex
    target(0) <- 1
    complex(0) <- 0

    # Predicted TMPRSS2 activity (% of baseline) and SARS-2-S viral entry
    # rate (% of baseline) by the empirical Hill-Langmuir link (Equation 4).
    tmprss2 <- 100 * target
    viralentry <- 100 * target^hill / (target^hill + ksp^hill)

    Cc ~ prop(propSd)
  })
}
