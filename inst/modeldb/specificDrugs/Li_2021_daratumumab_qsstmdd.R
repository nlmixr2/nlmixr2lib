Li_2021_daratumumab_qsstmdd <- function() {
  description <- "Two-compartment semi-mechanistic target-mediated drug-disposition (TMDD) population PK model for IV daratumumab (anti-CD38 IgG1) in adults with multiple myeloma, with parallel non-specific linear clearance and CD38-mediated saturable clearance under the quasi-steady-state (QSS) approximation of Gibiansky 2008. The TMDD/QSS form supersedes an earlier empirical Michaelis-Menten parameterisation with time-dependent Vmax: receptor (CD38) turnover and complex internalisation reproduce mechanistically the observed Vmax time-dependency. PAGE 29 (2021) abstract II-52 by Li, Perez Ruixo, Zhou, Perez Ruixo, and Dosne (Janssen R and D, Beerse). Distinct from Xu 2020 daratumumab, which uses the empirical 2-cmt parallel-linear / time-dependent Vmax form."
  reference <- "Li X, Perez Ruixo C, Zhou H, Perez Ruixo JJ, Dosne AG. Population-based Target-Mediated Drug Disposition (TMDD) Pharmacokinetics Model of Daratumumab in Patients With Multiple Myeloma Following Intravenously Daratumumab Monotherapy. PAGE 29 (2021) Abstr 9701. www.page-meeting.org/?abstract=9701"
  vignette <- "Li_2021_daratumumab_qsstmdd"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  # Covariates that the PAGE 29 abstract names as statistically significant
  # but for which neither the functional form nor the coefficient value is
  # reported in the source. They are documented here (not in covariateData)
  # so checkModelConventions() does not flag them as declared-but-unused;
  # promote each entry into covariateData when the full parameter table is
  # eventually published.
  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reported by Li 2021 as a statistically significant covariate on both linear CL and on V1; the abstract does not report the functional form or coefficient values. Effect omitted from model() pending publication of the full parameter table. See vignette Errata.",
      source_name        = "body weight"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reported by Li 2021 as a statistically significant covariate on linear CL; functional form and coefficient value are NOT in the abstract. Effect omitted from model() pending publication of the full parameter table. See vignette Errata.",
      source_name        = "albumin concentration"
    ),
    MM_NIGG = list(
      description        = "Multiple-myeloma immunoglobulin type indicator: 1 = non-IgG MM, 0 = IgG MM",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (IgG MM)",
      notes              = "Reported by Li 2021 as a statistically significant covariate on linear CL (described as 'type of myeloma [IgG versus non-IgG]'). Coefficient value and reference orientation are NOT in the abstract. Effect omitted from model() pending publication of the full parameter table. See vignette Errata. The canonical MM_NIGG column orientation (1 = non-IgG, 0 = IgG) is preserved.",
      source_name        = "type of myeloma (IgG vs non-IgG)"
    ),
    SEXF = list(
      description        = "Biological sex indicator: 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Reported by Li 2021 as a statistically significant covariate on V1. Coefficient value is NOT in the abstract. Effect omitted from model() pending publication of the full parameter table. See vignette Errata.",
      source_name        = "sex"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 223L,
    n_observations = 2572L,
    n_studies      = 2L,
    age_range      = "not reported in abstract",
    weight_range   = "not reported in abstract",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported in abstract",
    disease_state  = "Multiple myeloma. Pooled GEN501 + MMY2002 daratumumab monotherapy populations: GEN501 was a Phase 1/2 dose-escalation (0.005-24 mg/kg) + single-arm (8 or 16 mg/kg) study; MMY2002 (SIRIUS) was a Phase 2 study of heavily pretreated MM. Patients in both studies received daratumumab monotherapy with no other background therapies.",
    dose_range     = "0.005-24 mg/kg IV. Maintenance schedules included 8 mg/kg or 16 mg/kg IV once weekly (QW) for 8-16 weeks, every 2 weeks (Q2W) for 16 weeks, then every 4 weeks (Q4W) thereafter (GEN501 single-arm and MMY2002).",
    regions        = "Multinational (GEN501 + MMY2002 / SIRIUS clinical trial sites).",
    notes          = "Pooled GEN501 + MMY2002 daratumumab-monotherapy populations. 2,572 PK samples from 223 subjects (PAGE 29 abstract Results section). Demographic detail (age, weight distribution, sex split, race) is not reported in the abstract; the abstract names body weight, albumin, myeloma type (IgG vs non-IgG), and sex as statistically significant covariates without reporting coefficient values. The TMDD/QSS structural form replaces the previous empirical 2-cmt parallel-linear + Michaelis-Menten model with time-dependent Vmax (Xu 2017 / Xu 2020 lineage; see modellib('Xu_2020_daratumumab') for the empirical parameterisation)."
  )

  ini({
    # Structural PK parameters (typical-value estimates from Li 2021 PAGE 29
    # abstract Results section). Reference covariate values are not stated in
    # the abstract; the typical-value estimates below are therefore not
    # anchored to a specific reference WT / ALB / sex / myeloma-type. The
    # abstract reports CL = 0.0051 L/h, V1 = 3.76 L, Q = 0.0329 L/h, V2 = 3.31 L
    # ("L/h" in the printed abstract for V2 is a typo: V2 is a VOLUME, L).
    lcl <- log(0.0051); label("Non-specific linear (FcRn-mediated) clearance CL (L/hour)")           # Li 2021 PAGE 29 Results: CL 0.0051 L/h
    lvc <- log(3.76);   label("Central volume of distribution V1 (L)")                               # Li 2021 PAGE 29 Results: V1 3.76 L
    lq  <- log(0.0329); label("Inter-compartmental clearance Q (L/hour)")                            # Li 2021 PAGE 29 Results: Q 0.0329 L/h
    lvp <- log(3.31);   label("Peripheral volume of distribution V2 (L)")                            # Li 2021 PAGE 29 Results: V2 3.31 L (printed as 'L/h' in abstract; V2 is a volume)

    # TMDD / QSS parameters (Gibiansky 2008 parameterisation).
    # Ksyn is the zero-order target-synthesis rate constant (mg/L/h);
    # Kdeg is the first-order degradation rate of free target;
    # Kint is the first-order internalisation rate of the drug-target complex;
    # Kss is the QSS dissociation constant (mg/L), Kss = (koff + kint) / kon.
    # The abstract prints Kss in 1/h, which is dimensionally inconsistent with
    # the QSS expression and is treated here as a transcription error;
    # mg/L is the canonical concentration unit for Kss in the Gibiansky form.
    lksyn <- log(0.178);  label("Target (CD38) zero-order synthesis rate constant Ksyn (mg/L/hour)")  # Li 2021 PAGE 29 Results: Ksyn 0.178 mg/L/h
    lkdeg <- log(0.0082); label("Free-target first-order degradation rate constant Kdeg (1/hour)")    # Li 2021 PAGE 29 Results: Kdeg 0.0082 1/h
    lkint <- log(0.092);  label("Drug-target complex internalisation rate constant Kint (1/hour)")    # Li 2021 PAGE 29 Results: Kint 0.092 1/h
    lkss  <- log(1.79);   label("QSS steady-state dissociation constant Kss (mg/L; abstract prints 1/h, treated as typo)") # Li 2021 PAGE 29 Results: Kss 1.79 (units printed as 1/h; encoded as mg/L per Gibiansky 2008 QSS form)

    # Residual error. Li 2021 PAGE 29 reports RUV = 22.7% (proportional model
    # assumed, the standard form for daratumumab population PK).
    propSd <- 0.227; label("Proportional residual error (fraction)") # Li 2021 PAGE 29 Results: RUV 22.7%

    # Inter-individual variability (IIV) is NOT reported in the PAGE 29
    # abstract; no eta declarations are included here. The IIV magnitudes
    # are flagged in the vignette Errata so simulations with the loaded
    # model are typical-value trajectories only. Do NOT substitute
    # class-typical IIV placeholders.
  })
  model({
    # Individual structural parameters. No IIV or covariate effects are applied
    # because (a) the abstract does not report IIV magnitudes, and (b) the
    # abstract names body weight, albumin, myeloma type, and sex as significant
    # covariates without giving coefficient values. Both gaps are listed in the
    # vignette Errata.
    cl   <- exp(lcl)
    vc   <- exp(lvc)
    q    <- exp(lq)
    vp   <- exp(lvp)
    ksyn <- exp(lksyn)
    kdeg <- exp(lkdeg)
    kint <- exp(lkint)
    kss  <- exp(lkss)

    # Total target baseline (mg/L) at steady state with no drug.
    t0 <- ksyn / kdeg

    # Initial condition: free target at the synthesis/degradation equilibrium.
    total_target(0) <- t0

    # QSS algebra (Gibiansky et al. 2008 Eq. 7).
    # central is the total drug AMOUNT in the central compartment (mg, free +
    # complex); total_target is the total target CONCENTRATION (mg/L, free +
    # complex). Peripheral distribution acts on free drug only.
    ctot     <- central / vc
    ttot     <- total_target
    qss_disc <- (ctot - ttot)^2 + 2 * (ctot + ttot) * kss + kss^2
    complex  <- ((ctot + ttot + kss) - sqrt(qss_disc)) / 2
    cfree    <- ctot - complex
    tfree    <- ttot - complex

    # Two-compartment IV-input PK. Linear (non-specific, FcRn-mediated) CL
    # operates on the free drug concentration; the drug-target complex is
    # eliminated via internalisation at rate kint. Peripheral distribution
    # uses free drug only.
    d/dt(central)      <- -cl * cfree - q * cfree + q * (peripheral1 / vp) - kint * complex * vc
    d/dt(peripheral1)  <-               q * cfree - q * (peripheral1 / vp)
    d/dt(total_target) <-  ksyn - kdeg * tfree - kint * complex

    # Observation: total drug concentration in the central compartment.
    Cc <- ctot
    Cc ~ prop(propSd)
  })
}
