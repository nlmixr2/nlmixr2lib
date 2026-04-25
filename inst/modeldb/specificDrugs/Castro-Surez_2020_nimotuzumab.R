`Castro-Surez_2020_nimotuzumab` <- function() {
  description <- "Semi-mechanistic two-compartment QSS TMDD population PK model for nimotuzumab (anti-EGFR humanized IgG1) in adults with autosomal dominant polycystic kidney disease (Castro-Suarez 2020); EGFR binding represented in both central (Rtot) and peripheral (Rtotp) compartments under quasi-steady-state, plus a turnover mediator that stimulates non-specific clearance via a sigmoid Emax of free central nimotuzumab."
  reference <- "de Castro-Suarez N, Trame MN, Mangas-Sanjuan V, Garcia-Cremades M, Boix-Montanes A, Fernandez-Teruel C, Munoz-Camara A, Martin-Suarez A, Rebollo-Fernandez G, Lleonart-Vidal R. Semi-Mechanistic Pharmacokinetic Model to Guide the Dose Selection of Nimotuzumab in Patients with Autosomal Dominant Polycystic Kidney Disease. Pharmaceutics. 2020;12(12):1147. doi:10.3390/pharmaceutics12121147"
  vignette <- "Castro-Surez_2020_nimotuzumab"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    n_subjects     = 20L,
    n_observations = 422L,
    n_studies      = 1L,
    age_range      = "21-61 years (median 42; mean 39, SD 11)",
    age_median     = "42 years",
    weight_range   = "approx 38-94 kg (median 65.7; mean 66.98, SD 14.69)",
    weight_median  = "65.7 kg",
    height_median  = "163.5 cm (mean 163.60, SD 8.99)",
    bsa_median     = "1.70 m^2 (mean 1.72, SD 0.21)",
    serum_creatinine_median = "0.72 mg/dL (mean 0.77, SD 0.14)",
    crcl_median    = "105.7 mL/min/1.73 m^2 (mean 103.43, SD 22.63)",
    tkv_mean       = "Men 822.18 mL (SD 486.22; median 678.85). Women 924.14 mL (SD 404.27; median 846.55).",
    tcv_mean       = "339.93 mL (SD 201.19; median 310.3)",
    sex_female_pct = 70,
    race_ethnicity = c(Caucasian = 75, AfroAmerican = 5, Other = 20),
    disease_state  = "Adult patients with autosomal dominant polycystic kidney disease (ADPKD); inclusion required GFR >= 50 mL/min/1.73 m^2 and urinary protein excretion < 1 g/24 h.",
    dose_range     = "Single 30-min IV infusion at one of four fixed doses: 50, 100, 200, or 400 mg (n=5 per cohort; total n=20).",
    regions        = "Single-center Cuban phase I clinical trial (Cuban National Regulatory Agency CECMED 442/05.014.08-B); demographics from Castro-Suarez 2020 Table 1.",
    notes          = "All 422 serum nimotuzumab concentrations were quantifiable (none below LOQ). Concentrations measured by ELISA capturing nimotuzumab via recombinant EGFR extracellular domain. Covariate search (body weight, height, age, BSA, CrCL, serum creatinine, total kidney volume, total cyst volume, sex, race) found no covariate statistically significant on any PK parameter; the model is therefore reported without covariate effects."
  )

  ini({
    # Structural PK parameters - Castro-Suarez 2020 Table 2 (final model). Reference body weight 65 kg.
    lcl    <- log(9.64e-3); label("Linear (non-target-mediated) clearance CL (L/h) at typical mediator A3 = 1") # Castro-Suarez 2020 Table 2: CL = 9.64e-3 L/h
    lvc    <- log(2.63);    label("Apparent central volume of distribution V1 (L)")                               # Castro-Suarez 2020 Table 2: V1 = 2.63 L
    lq     <- log(2.88e-2); label("Inter-compartmental clearance Q (L/h)")                                       # Castro-Suarez 2020 Table 2: Q = 2.88e-2 L/h
    lvp    <- log(9.92e-3); label("Apparent peripheral volume of distribution V2 (L)")                            # Castro-Suarez 2020 Table 2: V2 = 9.92e-3 L

    # QSS TMDD parameters - shared Kss in both compartments, separate Rtot / Rtotp.
    lkss   <- log(15.5);    label("Quasi-steady-state binding constant of nimotuzumab to EGFR, Kss (mg/L)")        # Castro-Suarez 2020 Table 2: Kss = 15.5 mg/L
    lkint  <- log(4.94e-3); label("First-order internalization rate of nimotuzumab-EGFR complex, kint (1/h)")     # Castro-Suarez 2020 Table 2: kint = 4.94e-3 1/h
    lrtot  <- log(1.05e-2); label("Apparent EGFR concentration in central compartment, Rtot (mg/L)")              # Castro-Suarez 2020 Table 2: Rtot = 1.05e-2 mg/L
    lrtotp <- log(956);     label("Apparent EGFR concentration in peripheral compartment, Rtotp (mg/L)")          # Castro-Suarez 2020 Table 2: Rtotp = 956 mg/L

    # Mediator turnover parameters - sigmoid Emax stimulation of non-specific CL by free central nimotuzumab.
    # Hill coefficient gamma is fixed at 1 (it appears in equation (3) of the paper but is absent from
    # Table 2; operator confirmed gamma = 1, reducing the sigmoid to a hyperbolic Emax). With A3(0) = 1,
    # ksyn = kout * A3(0) = kout so kin and kout share the same value.
    lkout  <- log(1.33e-2); label("First-order degradation rate of mediator, kout (1/h); kin = kout x A3(0) = kout") # Castro-Suarez 2020 Table 2: kout = 1.33e-2 1/h
    lsmax  <- log(3.18);    label("Maximal stimulation of non-specific CL by mediator, Smax (unitless)")          # Castro-Suarez 2020 Table 2: Smax = 3.18
    ls50   <- log(8.57);    label("Free nimotuzumab concentration achieving half of Smax, S50 (mg/L)")             # Castro-Suarez 2020 Table 2: S50 = 8.57 mg/L

    # Inter-individual variability - statistically significant only on Rtotp and Kout. Exponential
    # (log-normal) IIV. omega^2 = log(CV^2 + 1):
    #   135% CV -> omega^2 = log(1.35^2 + 1) = 1.0379
    #   197% CV -> omega^2 = log(1.97^2 + 1) = 1.5853
    etalrtotp ~ 1.0379  # Castro-Suarez 2020 Table 2: Rtotp IIV 135% CV (eta-shrinkage 14%)
    etalkout  ~ 1.5853  # Castro-Suarez 2020 Table 2: Kout  IIV 197% CV (eta-shrinkage 21%)

    # Residual error - paper reports an additive RUV on the log-transformed observations of 48% (4),
    # which maps to proportional error in nlmixr2's linear-space concentration parameterization.
    propSd <- 0.48; label("Proportional residual error (fraction); reported by paper as additive on log scale = 48%") # Castro-Suarez 2020 Table 2: 48%
  })

  model({
    # Individual PK parameters. Only Rtotp and Kout carry IIV per the published final model.
    cl    <- exp(lcl)
    vc    <- exp(lvc)
    q     <- exp(lq)
    vp    <- exp(lvp)
    kss   <- exp(lkss)
    kint  <- exp(lkint)
    rtot  <- exp(lrtot)
    rtotp <- exp(lrtotp + etalrtotp)
    kout  <- exp(lkout  + etalkout)
    smax  <- exp(lsmax)
    s50   <- exp(ls50)

    # kin from steady-state condition on the mediator (Castro-Suarez 2020 §2.3, equation (3) initial
    # conditions): A3(0) = 1 and dA3/dt = 0 at baseline imply kin = kout * A3(0) = kout.
    kin <- kout

    # Initial mediator (effect) state (dimensionless; baseline = 1 per paper). The canonical
    # nlmixr2lib compartment name "effect" is used for the indirect-response turnover species the
    # paper calls A3 / "mediator"; its level multiplies the linear CL term in central.
    # central and peripheral1 default to 0.
    effect(0) <- 1

    # ---- QSS TMDD: solve free nimotuzumab concentration in each compartment ----
    # Total drug concentration in the central / peripheral compartments (mg/L = ug/mL):
    #   ctot = Atotal / V
    # Free concentration cfree solves: ctot = cfree + R * cfree / (Kss + cfree), giving the closed form
    #   cfree = 0.5 * ((ctot - R - Kss) + sqrt((ctot - R - Kss)^2 + 4 * Kss * ctot))
    # which is the standard QSS-TMDD algebraic free-drug solution (Gibiansky et al. 2008).
    ctot1 <- central / vc
    disc1 <- ctot1 - rtot - kss
    c1f   <- 0.5 * (disc1 + sqrt(disc1 * disc1 + 4 * kss * ctot1))

    ctot2 <- peripheral1 / vp
    disc2 <- ctot2 - rtotp - kss
    c2f   <- 0.5 * (disc2 + sqrt(disc2 * disc2 + 4 * kss * ctot2))

    # Internalized (bound) drug amount in central compartment, used by the kint elimination term.
    bound1 <- vc * rtot * c1f / (kss + c1f)

    # ---- ODEs (Castro-Suarez 2020 equations (1)-(3)) ----
    # central tracks the TOTAL nimotuzumab amount in central compartment (mg); the dose enters here.
    # Linear elimination is on free drug, modulated by the mediator state A3 (=effect compartment).
    # Distribution between central <-> peripheral is on FREE drug only.
    # kint elimination is on the bound complex amount in central; no kint elimination from peripheral
    # (paper §2.3: "No elimination process of the complex was assumed from the peripheral compartment").
    d/dt(central)     <- -cl * effect * c1f - q * c1f + q * c2f - kint * bound1
    d/dt(peripheral1) <-  q * c1f - q * c2f

    # Mediator (effect): zero-order synthesis kin stimulated by a sigmoid Emax of free central drug
    # (gamma fixed at 1 -> hyperbolic Emax), and first-order degradation by kout.
    d/dt(effect)      <-  kin * (1 + smax * c1f / (s50 + c1f)) - kout * effect

    # Observation: total nimotuzumab concentration in serum (free + EGFR-bound), measured by the
    # receptor-binding ELISA. Total = ctot1 by construction of the QSS solution.
    Cc <- ctot1
    Cc ~ prop(propSd)
  })
}
