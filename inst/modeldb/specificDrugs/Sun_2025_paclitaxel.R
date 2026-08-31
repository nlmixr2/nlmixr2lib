Sun_2025_paclitaxel <- function() {
  description <- "Two-compartment population PK model of intravenous paclitaxel coupled to a pharmacodynamic model of chemotherapy-induced peripheral neuropathy (CIPN), in which central-compartment concentrations above an estimated threshold concentration drive a pure-integrator effect compartment that stimulates the input rate of an EORTC QLQ-CIPN20 sensory-subscale (CIPN8) turnover pool, in adult women with breast cancer receiving weekly 1-h paclitaxel infusions (Sun 2025)"
  reference <- "Sun Y, Pai MP, Henry NL, Hertz DL. Pharmacokinetic-Pharmacodynamic Model of Paclitaxel-Induced Peripheral Neuropathy. Clin Transl Sci. 2025;18(11):e70404. doi:10.1111/cts.70404"
  vignette <- "Sun_2025_paclitaxel"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # The CIPN8 turnover state is the paper's patient-reported efficacy/toxicity
  # endpoint: the 8-item sensory subscale of the EORTC QLQ-CIPN20, summed and
  # then offset by -8 so that 0 = no symptoms and 24 = maximal symptoms
  # (Methods 2.1). It is the direct analogue of the registered `das28` and
  # `cdai` disease-activity-score output compartments (precedent:
  # Zhang_2023_brazikumab_crp.R) but is not itself a registered canonical, so
  # it is declared here rather than silently introduced.
  paper_specific_compartments <- c("cipn8")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: analyte and specimen checked against
  # Figure 2 (model schematic) and Methods 2.1/2.4 of the source paper.
  compartmentData <- list(
    central     = list(analyte = "paclitaxel", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "paclitaxel", units = "mg", specimen = "plasma", verified = TRUE),
    # The effect compartment is a latent kinetic state, not a biological
    # matrix. Table 2's footnote labels its value in ng/mL ("The effect
    # compartment concentration across individual is 1.56 +/- 1.18 ng/mL"),
    # which is the paper's own label; see the dimensional note on `lke0` in
    # ini() and the vignette Errata.
    effect      = list(analyte = "paclitaxel", units = "ng/mL", specimen = "not applicable", verified = TRUE),
    cipn8       = list(analyte = "EORTC QLQ-CIPN20 sensory subscale score (CIPN8)", units = "score", specimen = "not applicable", verified = TRUE)
  )

  # No covariate effects were retained in the final model: Table 2 reports only
  # structural fixed effects, IIV standard deviations and residual-error terms,
  # and neither Methods 2.2 nor 2.4 describes a covariate search. Body surface
  # area determines the administered dose (Methods 2.1) but enters the model as
  # the `amt` of the infusion record, not as a covariate on any parameter.
  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 60,
    n_studies      = 1,
    age_range      = "28-71 years",
    age_median     = "52.30 years (mean)",
    sex_female_pct = 100,
    race_ethnicity = c(White = 93.3, Black = 3.3, Asian = 3.3),
    disease_state  = "stage I-III or oligometastatic breast cancer",
    dose_range     = "80 mg/m2 IV weekly for 12 weeks (1-h infusion; 90-min for the first dose)",
    regions        = "United States (University of Michigan Rogel Cancer Center)",
    notes          = paste(
      "UMCC2014.002 observational study (NCT02338115); 65 enrolled, 60 analysed after 5",
      "exclusions for withdrawal or protocol violation (Results 3.1). Baseline demographics,",
      "dosing and exposure summarised in Table 1: mean Cmax 2364.16 ng/mL (range 907-4340),",
      "mean 18-26 h concentration 22.02 ng/mL (range 10-51.80), mean cumulative dose 1617.20 mg",
      "(range 456-2412), mean 11 of 12 planned doses received (range 3-12) with 68% completing",
      "all 12. Mean baseline CIPN8 0.3 (range 0-3) with 83% scoring 0; mean change in CIPN8",
      "6.03 (range 0-21). Two plasma samples per patient from the first infusion only (final",
      "10 min of the infusion, and 16-26 h after infusion start), so the PK model is informed by",
      "sparse peak-and-trough data; 57 patients had an end-of-infusion sample and 59 had a",
      "post-infusion sample. External evaluation used published CIPN incidences from the CALGB",
      "C9840 trial (no PK data). Generalisability is limited to weekly dosing in females who",
      "were 93% white with breast cancer (Discussion)."
    )
  )

  ini({
    # ---- PK: two-compartment model, IV infusion (Table 2, "PK model" panel) ----
    # Reported as untransformed population values; log-transformed here because
    # Table 2's footnote states the omega terms model IIV as log-normal.
    lcl <- log(26.04); label("Clearance (L/h)")                                                   # Table 2 Cl_pop = 26.04 L/hour (SE 1.17, RSE 4.51%)
    lvc <- log(34.08); label("Central volume of distribution (L)")                                # Table 2 V_pop = 34.08 L (SE 2.44, RSE 7.15%)
    lq  <- log(4.2);   label("Intercompartmental clearance (L/h)")                                 # Table 2 Q_pop = 4.2 L/hour (SE 0.17, RSE 4.01%)
    lvp <- log(67.67); label("Peripheral volume of distribution (L)")                              # Table 2 V2_pop = 67.67 L (SE 6.2, RSE 9.17%)

    # ---- PD: threshold + effect compartment + turnover (Table 2, "PD model" panel) ----
    # `cthres` is the paper's threshold concentration `x` (Equation 1). It is a
    # NEW canonical registered with this extraction: an estimated plasma
    # concentration that hard-gates a downstream PD effect. Distinct from the
    # occupancy-scoped `thres` (a percent) and from the measured `mic`.
    lcthres <- log(1735.75); label("Threshold paclitaxel concentration in the central compartment above which CIPN accrues (ng/mL)")  # Table 2 x_pop = 1735.75 ng/mL (SE 1.84, RSE 0.106%); restated in Results 3.2 and Discussion
    # Applied verbatim in the paper's printed units. Equation 3 is
    # d ce/dt = cact * ke0 with cact in ng/mL and ke0 in mL/ng, which makes the
    # right-hand side dimensionless rather than a rate; the paper's unit labels
    # for ke0 and for the effect compartment are not mutually consistent. The
    # numeric value and the equation are reproduced exactly as printed, and the
    # resulting effect-compartment magnitude is verified against Table 2's
    # footnote value in the vignette. See vignette Errata.
    lke0    <- log(0.00022); label("Effect-compartment input rate constant (mL/ng)")                # Table 2 ke0_pop = 0.00022 mL/ng (SE 0.000063, RSE 28.4%)
    lkout   <- log(0.0027);  label("First-order output rate constant of the CIPN8 turnover pool (1/h)")  # Table 2 kout_pop = 0.0027 (SE 0.0016, RSE 58.9%); units not printed, inferred as 1/h from the L/hour PK panel that shares the time axis
    lslope  <- log(21.65);   label("Linear coefficient of the effect compartment on the CIPN8 input rate (per effect-compartment unit)")  # Table 2 k_pop = 21.65 (SE 9.52, RSE 44.0%); the `k` of Equation 5
    # R0_pop is reported with SE and RSE both NA, i.e. not estimated: it was
    # held at the cohort mean baseline CIPN8 score of 0.3 (Table 1).
    lrbase  <- fixed(log(0.3)); label("Baseline CIPN8 score (0-24 scale)")                          # Table 2 R0_pop = 0.3 with SE = NA and RSE = NA -> fixed; equals the Table 1 mean baseline CIPN8 of 0.3

    # ---- IIV: log-normal (Table 2 footnote: "The omega values are standard
    # deviation", and IIV is modelled "in log normal distribution"). nlmixr2
    # eta declarations take VARIANCES, so each entry is the reported omega
    # squared. The PD omegas are large and carry high shrinkage (reported in
    # Table 2), which the vignette Errata discusses.
    etalcl    ~ 0.0625   # Table 2 Omega_Cl = 0.25 (SD, RSE 15.6%) -> variance 0.25^2
    etalvc    ~ 0.2601   # Table 2 Omega_V1 = 0.51 (SD, RSE 9.54%) -> variance 0.51^2
    etalq     ~ 0.04     # Table 2 Omega_Q = 0.2 (SD, RSE 20.0%) -> variance 0.2^2
    etalvp    ~ 0.0529   # Table 2 Omega_V2 = 0.23 (SD, RSE 36.9%) -> variance 0.23^2
    etalkout  ~ 6.0025   # Table 2 Omega_kout = 2.45 (SD, RSE 21.5%, shrinkage -75%) -> variance 2.45^2
    etalslope ~ 3.6481   # Table 2 Omega_k = 1.91 (SD, RSE 14.7%, shrinkage 52%) -> variance 1.91^2
    etalke0   ~ 0.7225   # Table 2 Omega_ke0 = 0.85 (SD, RSE 40.2%, shrinkage 90%) -> variance 0.85^2

    # ---- Residual error ----
    # PK: constant (additive) error model. Table 2 footnote:
    # "Concentration = Cc + a x e", with e ~ N(0, 1).
    addSd <- 9.11; label("Additive residual error on paclitaxel concentration (ng/mL)")             # Table 2 PK error model a = 9.11 ng/mL (SE 1.27, RSE 13.9%); Results 3.2 "constant error model"
    # PD: Monolix combined1 with a power term. Table 2 footnote:
    # "CIPN8 = R + (a1 + b1 x R^0.74) x e", i.e. SD = a1 + b1 * R^0.74, a LINEAR
    # sum of an additive term and a power term -> combined1() in nlmixr2.
    addSd_cipn8  <- 1.62;  label("Additive component of the CIPN8 residual error (score)")          # Table 2 PD error model a1 = 1.62 (SE 0.059, RSE 3.66%)
    propSd_cipn8 <- 0.056; label("Power-term coefficient of the CIPN8 residual error")              # Table 2 PD error model b1 = 0.056 (SE 0.018, RSE 32.1%)
    # The exponent 0.74 appears only inside the printed error expression and is
    # absent from Table 2's estimated-parameter rows (which list a1 and b1
    # only, each with an SE and RSE), so it was held fixed rather than
    # estimated.
    powExp_cipn8 <- fixed(0.74); label("Power exponent on the CIPN8 score in the residual error")   # Table 2 footnote "CIPN8 = R + (a1 + b1 x R^0.74) x e"; not reported as an estimated parameter
  })

  model({
    # ---- Individual parameters ----
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    ke0    <- exp(lke0 + etalke0)
    kout   <- exp(lkout + etalkout)
    slope  <- exp(lslope + etalslope)
    cthres <- exp(lcthres)
    rbase  <- exp(lrbase)

    # ---- Micro-constants ----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- PK: two-compartment model; states are amounts in mg (Figure 2) ----
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Central concentration in ng/mL. central is in mg and vc in L, so
    # central/vc is mg/L = ug/mL; the factor of 1000 converts to the ng/mL used
    # throughout Table 1, Table 2 and Figure 1.
    Cc <- central / vc * 1000

    # ---- Threshold gate (Equation 1) ----
    # f(x): cact = 0 for c < x, and cact = c - x for c >= x, where x is the
    # threshold concentration in the central compartment. A hard switch, not a
    # saturable function.
    cact <- max(Cc - cthres, 0)

    # ---- Effect compartment (Equations 2 and 3) ----
    # A pure integrator with no loss term, exactly as printed: the effect
    # compartment only accumulates. The Discussion confirms this is intended
    # ("using the effect compartment alone is not optimal since the CIPN will
    # always increase"), which is why the turnover pool downstream supplies the
    # resolution behaviour.
    effect(0)    <- 0
    d/dt(effect) <- cact * ke0

    # ---- CIPN8 turnover (Equations 4 and 5) ----
    kin <- rbase * kout
    d/dt(cipn8) <- kin * (1 + slope * effect) - kout * cipn8
    cipn8(0) <- rbase

    # ---- Observations ----
    Cc ~ add(addSd)
    cipn8 ~ add(addSd_cipn8) + pow(propSd_cipn8, powExp_cipn8) + combined1()
  })
}
