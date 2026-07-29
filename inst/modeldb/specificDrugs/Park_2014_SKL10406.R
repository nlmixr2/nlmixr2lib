Park_2014_SKL10406 <- function() {
  description <- "Two-compartment first-order oral absorption population PK with effect-compartment Emax PK-PD model for striatal serotonin transporter (SERT) occupancy by SKL10406 (a triple monoamine reuptake inhibitor candidate) in healthy adult volunteers (Park 2014; EME variant, Table 3)"
  reference <- "Park JS, Lee J, Meyer J, Ilankumaran P, Han S, Yim DS. Serotonin transporter occupancy of SKL10406 in humans: comparison of pharmacokinetic-pharmacodynamic modeling methods for estimation of occupancy parameters. Transl Clin Pharmacol. 2014;22(2):83-91. doi:10.12793/tcp.2014.22.2.83"
  vignette <- "Park_2014_SKL10406"
  units <- list(time = "hr", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 11L,
    n_studies      = 1L,
    age_range      = "18-50 years (inclusion criterion); 42.5 +/- 4.51 in 50 mg bid cohorts, 44.8 +/- 11.82 in 75 mg bid cohorts",
    age_median     = "43.5 years (pooled mean across 11 completers)",
    weight_range   = "Body weight < 125.0 kg (inclusion criterion); BMI 19.0-30.0 kg/m^2 (inclusion), observed BMI 25.7 +/- 2.90 kg/m^2",
    weight_median  = NULL,
    sex_female_pct = 9.1,
    race_ethnicity = c(White = 90.9, Black_or_African_American = 9.1),
    disease_state  = "Healthy adult volunteers; key exclusions included use of medications affecting SERT or DAT binding within 1 week prior to first dosing, and history of serious medical or psychiatric illness.",
    dose_range     = "100 mg/day SKL10406 (50 mg PO every 12 h) for 6 days in Cohorts 1 and 2; 100 mg/day for 4 days followed by 150 mg/day (75 mg PO every 12 h) for 6 days in Cohorts 3 and 4. PD model was fit to the SERT cohorts only (Cohorts 1 and 3; n = 6 subjects).",
    regions        = "Canada (single-centre PET study at the Centre for Addiction and Mental Health, Toronto; PK / PD sampling at Kendle Early Stage, Toronto).",
    n_observations = "Full PK sampling for all 11 completers at 0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 12, and 16 h post-dose on the steady-state dosing day. Two post-dose SERT-occupancy PET scans per subject (approximately 4 h and 16 h post-dose) in Cohorts 1 and 3 only.",
    notes          = "Fifteen subjects enrolled; 11 completed both PK and PD assessments (six SERT cohorts 1 and 3, five DAT cohorts 2 and 4). One subject in the SERT analysis refused the 16 h Day 11 PET scan. One subject was suspected to be a poor metabolizer based on elevated plasma SKL10406 levels and was retained in the analysis. Baseline demographics per Park 2014 Table 1."
  )

  ini({
    # Structural PK parameters - Park 2014 Table 3, ME row (one PK model shared across all
    # three PD methods). Two-compartment first-order oral absorption with absorption lag.
    lcl   <- log(49.6);  label("Apparent clearance CL/F (L/h)")                              # Park 2014 Table 3 ME: CL = 49.6 L/h (RSE 19.0%; bootstrap 43.7, 95% CI 32.9-68.9)
    lvc   <- log(176);   label("Apparent central volume of distribution V2/F (L)")          # Park 2014 Table 3 ME: V2 = 176 L (RSE 14.4%; bootstrap 169, 95% CI 15.9-221)
    lka   <- log(5.05);  label("First-order absorption rate constant ka (1/h)")             # Park 2014 Table 3 ME: ka = 5.05 1/h (RSE 57.6%; bootstrap 18.6, 95% CI 0.51-32.4); no IIV
    lvp   <- log(63.7);  label("Apparent peripheral volume of distribution V3/F (L)")       # Park 2014 Table 3 ME: V3 = 63.7 L (RSE 33.2%; bootstrap 51.2, 95% CI 38.8-133); no IIV
    lq    <- log(21.1);  label("Apparent inter-compartmental clearance Q/F (L/h)")          # Park 2014 Table 3 ME: Q = 21.1 L/h (RSE 18.3%; bootstrap 19.0, 95% CI 15.0-32.8); no IIV
    ltlag <- log(0.336); label("Absorption lag time ALAG (h)")                              # Park 2014 Table 3 ME: ALAG = 0.336 h (RSE 27.0%; bootstrap 0.429, 95% CI 0.18-0.48); no IIV. Inclusion of ALAG decreased OFV by 7.95.

    # PD parameters - Park 2014 Table 3, EME row (effect-compartment Emax model on
    # striatal SERT occupancy). Emax is reported in percent occupancy (0-100); EC50 is in
    # ng/mL on the effect compartment concentration scale; ke0 is the first-order
    # equilibration rate (equilibration t1/2 = log(2) / 0.288 = 2.41 h).
    lemax <- log(68.6);  label("Maximum striatal SERT occupancy Emax (%)")                  # Park 2014 Table 3 EME: Emax = 68.6% (RSE 4.62%; bootstrap 69.9, 95% CI 50.3-104.8). IIV could not be estimated; the estimation step did not converge with eta on Emax (paper Discussion).
    lec50 <- log(40.2);  label("Effect-compartment concentration at 50% Emax, EC50 (ng/mL)") # Park 2014 Table 3 EME: EC50 = 40.2 ng/mL (RSE 34.3%; bootstrap 44.2, 95% CI 12.9-140)
    lke0  <- log(0.288); label("Effect-compartment equilibration rate constant ke0 (1/h)")  # Park 2014 Table 3 EME: ke0 = 0.288 1/h (RSE 9.65%; bootstrap 0.286, 95% CI 0.202-0.409). Equilibrium half-life = 2.41 h (paper Results).

    # IIV - log-normal CV% reported in Park 2014 Table 3.
    # omega^2 = log(CV^2 + 1).
    # log(0.622^2 + 1) = 0.3270; log(0.416^2 + 1) = 0.1598; log(0.683^2 + 1) = 0.3826.
    etalcl   ~ 0.3270  # Park 2014 Table 3 ME: IIV CL/F = 62.2% CV
    etalvc   ~ 0.1598  # Park 2014 Table 3 ME: IIV V2/F = 41.6% CV
    etalec50 ~ 0.3826  # Park 2014 Table 3 EME: IIV EC50 = 68.3% CV; effect-compartment scale.

    # Residual error.
    # Park 2014 Methods states a "combined" (additive + proportional) residual error model
    # for plasma SKL10406, and a "proportional" residual error model for the EME SERT
    # occupancy fit ("successful convergence was not achieved with either the additive or
    # combined error models for both DME and EME"). The numeric magnitudes of these residual
    # errors are NOT tabulated in Park 2014 Table 3 nor in the paper text. The values below
    # are pragmatic library defaults wrapped in fixed() so that downstream users can simulate
    # plausible prediction intervals; they are NOT the paper's estimated residual error
    # magnitudes. See the vignette Errata for details.
    propSd    <- fixed(0.20); label("Proportional residual error on plasma SKL10406 (fraction; library default; not reported in Park 2014)") # Library default; see vignette Errata
    addSd     <- fixed(1.0);  label("Additive residual error on plasma SKL10406 (ng/mL; library default; not reported in Park 2014)")        # Library default; see vignette Errata
    propSd_TO <- fixed(0.15); label("Proportional residual error on SERT occupancy TO (fraction; library default; not reported in Park 2014)") # Library default; see vignette Errata
  })
  model({
    # Individual PK parameters.
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc)
    ka   <- exp(lka)
    vp   <- exp(lvp)
    q    <- exp(lq)
    tlag <- exp(ltlag)

    # Individual PD parameters. Emax has no IIV (paper Discussion: estimation did not
    # converge with eta on Emax in the EME variant). EC50 carries log-normal IIV. ke0 has
    # no IIV.
    emax <- exp(lemax)
    ec50 <- exp(lec50 + etalec50)
    ke0  <- exp(lke0)

    # Micro-rate constants for the explicit two-compartment ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Plasma concentration: dose in mg, volume in L -> central / vc in mg/L. Multiply by
    # 1000 to express Cc in ng/mL to match Park 2014 Table 3 units (EC50 in ng/mL).
    Cc <- central / vc * 1000

    # ODE system: depot -> central -> peripheral1; effect compartment driven by plasma
    # concentration Cc, with first-order equilibration rate ke0 (Park 2014 Eq. 3).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(effect)      <-  ke0 * (Cc - effect)

    # Absorption lag time on the depot compartment.
    alag(depot) <- tlag

    # Striatal SERT occupancy (% target sites blocked). Emax is in percent units; effect
    # holds the effect-compartment SKL10406 concentration in ng/mL.
    TO <- emax * effect / (ec50 + effect)

    Cc ~ add(addSd) + prop(propSd)
    TO ~ prop(propSd_TO)
  })
}
