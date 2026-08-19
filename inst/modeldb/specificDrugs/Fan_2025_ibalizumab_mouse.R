Fan_2025_ibalizumab_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). Mechanism-based quasi-equilibrium (QE)",
    "target-mediated drug disposition (TMDD) PK model coupled to an HIV-1",
    "viral-dynamic PD model for ibalizumab (anti-CD4 humanized IgG4, about",
    "150 kDa) run as the positive control arm of the Fan 2025 anti-CD4",
    "nanobody study in HIV-1 CH058-infected humanized NDG-HuPBL mice.",
    "Structure is identical to the companion Fan_2025_nb457trimer_mouse.R",
    "except that ibalizumab disposition is TWO-compartment: 'The PK data of",
    "Ibalizumab were fitted to a two-compartment disposition model. The other",
    "model structure and parameters are the same as Nb457-NbHSA-Nb457'",
    "(Table 1 footnote). Parallel intraperitoneal and subcutaneous first-order",
    "absorption depots feed a central state carrying TOTAL drug (free plus",
    "CD4-bound); the free concentration Cc is recovered at every integration",
    "step from the QE quadratic (eq 5), is cleared linearly (CL), exchanges",
    "with peripheral1 through Q, and is removed as drug-receptor complex at",
    "kint. Total CD4 receptor turns over with kdeg fixed equal to kint, so",
    "total receptor stays at RTOT0 and ksyn = kint * RTOT0. The PD layer is",
    "the same two-state viral-dynamic model (infected CD4+ T-cell pool",
    "expanding at kPV, lost at kCV, drug acting as an Imax inhibition of the",
    "kPV term, feeding a free virus pool through beta toward",
    "lambda * infected). Ibalizumab's fitted Imax = 0.42 is below 1, so the",
    "drug slows but never reverses infected-cell expansion - the mechanistic",
    "contrast with the nanobody, whose Imax = 1.963. The ibalizumab target",
    "parameters are very poorly identified (RTOT0 %RSE 249.6, K_M 308.5,",
    "kint 1121); see the validation vignette Errata."
  )
  reference <- paste(
    "Fan X, Cao K, Wu X, Yan X. Pharmacokinetic-pharmacodynamic modeling of a",
    "highly potent and broadly neutralizing anti-CD4 trimeric nanobody to",
    "inhibit HIV-1 infection. Microbiol Spectr. 2025;13(11):e00805-25.",
    "doi:10.1128/spectrum.00805-25. PMCID: PMC12584636.",
    "Structural equations from eqs 1-7 (pp. 3-4) plus the Table 1 footnote",
    "declaring a two-compartment ibalizumab disposition; PK parameter values",
    "from Table 1 (Ibalizumab rows); PD parameter values from Table 2",
    "(Ibalizumab rows)."
  )
  vignette <- "Fan_2025_antiCD4_nanobody_HIV"
  units <- list(time = "h", dosing = "ug", concentration = "ug/mL")

  compartmentData <- list(
    depot_ip     = list(analyte = "ibalizumab", units = "ug", specimen = "administration site", verified = TRUE),
    depot_sc     = list(analyte = "ibalizumab", units = "ug", specimen = "administration site", verified = TRUE),
    central      = list(analyte = "ibalizumab (total: free plus CD4-bound)", units = "ug", specimen = "serum", verified = TRUE),
    peripheral1  = list(analyte = "ibalizumab", units = "ug", specimen = "tissue", verified = FALSE),
    total_target = list(analyte = "CD4 receptor (total, free plus drug-bound; carried as a drug-equivalent concentration)", units = "ug/mL", specimen = "blood cell", verified = FALSE),
    infected     = list(analyte = "HIV-1-infected CD4+ T cells", units = "cells/mL", specimen = "whole blood", verified = FALSE),
    virus        = list(analyte = "HIV-1 RNA", units = "copies/mL", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "mouse (NDG-HuPBL humanized, female)",
    n_subjects     = 8L,
    n_studies      = 1L,
    age_range      = "Not reported.",
    weight_range   = "Not reported. A reference body weight of 0.02 kg is implied by the paper's own allometric projection (V3 = 0.54 mL reported as 26.9 mL/kg in humans with a volume exponent of 1).",
    sex_female_pct = 100,
    race_ethnicity = "Not applicable (murine study).",
    disease_state  = "HIV-1 CH058 infection established by challenge with 10 ng p24 on day 0.",
    dose_range     = "400 ug per mouse, identical to the nanobody arms. PK arms: a single 400 ug intraperitoneal (IP) OR subcutaneous (SC) dose on day 1 post-infection, blood sampled at 0, 1, 4, 8, 12, 24, 72 and 120 h. PD arm: 400 ug IP on day 1 followed by 400 ug SC on days 3, 5 and 7, viral load sampled at weeks 1, 2, 3 and 4 post-infection.",
    regions        = "Nanjing, People's Republic of China.",
    notes          = "Ibalizumab was administered as the positive control on the same regimen as Nb457-NbHSA-Nb457 (Materials and Methods, Data; Table S1); n = 4 female mice per arm. Fitted by a naive pooled-data approach, yet Table 1 still reports estimated IIV variances for CL and V3, which this model file carries; see the vignette Errata."
  )

  ini({
    # ---- Absorption (Table 1, Ibalizumab rows) --------------------------
    lka_ip <- log(0.26)
    label("First-order absorption rate constant after intraperitoneal injection (1/h)")  # Table 1, Ibalizumab row "k a (IP) (1/h)" = 0.26 (%RSE 46.85)
    lka_sc <- log(0.22)
    label("First-order absorption rate constant after subcutaneous injection (1/h)")     # Table 1, Ibalizumab row "k a (SC) (1/h)" = 0.22 (%RSE 58.57)

    # ---- Linear two-compartment disposition (Table 1, Ibalizumab rows) --
    lcl <- log(0.024)
    label("Linear clearance of free drug from the central compartment (mL/h)")   # Table 1, Ibalizumab row "CL (mL/h)" = 0.024 (%RSE 20.29)
    lvc <- log(0.54)
    label("Volume of distribution of the central compartment (mL)")              # Table 1, Ibalizumab row "V 3 (mL)" = 0.54 (%RSE 35.27)
    lq  <- log(1.86)
    label("Inter-compartmental (tissue distribution) clearance (mL/h)")          # Table 1, Ibalizumab row "Q (mL/h)" = 1.86 (%RSE 35.27)
    lvp <- log(0.12)
    label("Volume of distribution of the peripheral compartment V4 (mL)")        # Table 1, Ibalizumab row "V 4 (mL)" = 0.12 (%RSE 14.25)

    # ---- Target binding and turnover (Table 1, Ibalizumab rows) ---------
    lkm <- log(12.73)
    label("Quasi-equilibrium drug-receptor dissociation constant K_M (ug/mL)")  # Table 1, Ibalizumab row "K M (ug/mL)" = 12.73 (%RSE 308.5)

    lrbase_total_target <- log(11.73)
    label("Baseline total CD4 receptor RTOT0, in drug-equivalent concentration (ug/mL)")  # Table 1, Ibalizumab row "RTOT 0 (ug/mL)" = 11.73 (%RSE 249.6)

    lkint <- log(0.0036)
    label("Drug-receptor complex internalization rate constant kint, also used for free-receptor degradation kdeg (1/h)")  # Table 1, Ibalizumab row "k INT (1/h)" = 0.0036 (%RSE 1121)

    # ---- Viral dynamics (Table 2, Ibalizumab rows) ----------------------
    lkpv <- log(0.027)
    label("Viral infection / infected-cell expansion rate constant kPV (1/h)")  # Table 2, Ibalizumab row "k PV (1/h)" = 0.027 (%RSE 12.68)
    lkcv <- log(0.0072)
    label("Death rate constant of infected cells kCV (1/h)")                    # Table 2, Ibalizumab row "k CV (1/h)" = 0.0072 (%RSE 38.42)
    lbeta_virus <- log(2.09e-7)
    label("First-order rate constant for HIV-1 release from infected cells beta (1/h)")  # Table 2, Ibalizumab row "beta (1/h)" = 2.09e-7 (%RSE 58)
    llambda_virus <- log(2.001)
    label("Scaling ratio lambda relating the infected-cell pool to viral load (unitless)")  # Table 2, Ibalizumab row "lambda" = 2.001 (%RSE 38.67)
    lrbase_infected <- log(1225)
    label("Estimated initial condition of the infected CD4+ T-cell pool BASE_T (cells/mL)")  # Table 2, Ibalizumab row "BASE T" = 1,225 (%RSE 36.38)

    # ---- Drug effect on viral dynamics (Table 2, Ibalizumab rows) -------
    limax <- log(0.42)
    label("Maximum inhibitory effect of ibalizumab on HIV-1 infection Imax (unitless)")  # Table 2, Ibalizumab row "I max" = 0.42 (%RSE 51.79)

    # Fixed, not estimated - Table 2 footnote c "The IC 50 was fixed to the
    # observed value". The Ibalizumab row is labelled SC50 (the symbol used in
    # eq 6) whereas the Nb457 row is labelled IC50; they are the same
    # quantity.
    lic50 <- fixed(log(0.1449))
    label("Free concentration producing half the maximum inhibition, SC50 in eq 6, set to the observed value (ug/mL)")  # Table 2, Ibalizumab row "SC 50 (ug/mL)" = 0.1449, FIXED (footnote c)

    # ---- Between-subject variability (Table 1, Ibalizumab rows) ---------
    # Table 1 footnote a establishes that the tabulated omega and sigma are
    # NONMEM VARIANCES.
    etalcl ~ 0.29   # Table 1, Ibalizumab row omega CL = 0.29 (variance; %RSE 22.95)
    etalvc ~ 0.68   # Table 1, Ibalizumab row omega V3 = 0.68 (variance; %RSE 40.53)

    # ---- Residual error -------------------------------------------------
    # Table 1 reports only a proportional sigma for ibalizumab (no additive
    # counterpart to the nanobody's sigma 2). The row's Description cell reads
    # "Proportional error of Nb 457 -Nb HSA -Nb 457", which is a copy-paste
    # slip in the published table - the row sits in the Ibalizumab block and
    # carries the ibalizumab estimate. See the vignette Errata.
    propSd <- sqrt(0.053)
    label("Proportional residual error SD on serum concentration (fraction)")  # Table 1, Ibalizumab row "sigma 1" = 0.053 (variance; %RSE 8.69) -> SD = sqrt(0.053)
    propSd_virus <- sqrt(0.94)
    label("Proportional residual error SD on viral load (fraction)")           # Table 2, Ibalizumab row "sigma 3" = 0.94 (variance; %RSE 23.09) -> SD = sqrt(0.94)
  })

  model({
    # ---- Individual parameters ------------------------------------------
    ka_ip <- exp(lka_ip)
    ka_sc <- exp(lka_sc)
    cl    <- exp(lcl + etalcl)
    vc    <- exp(lvc + etalvc)
    q     <- exp(lq)
    vp    <- exp(lvp)
    km    <- exp(lkm)
    kint  <- exp(lkint)

    # Paper assumption (text after eq 5): kint == kdeg.
    kdeg  <- kint

    rbase_total_target <- exp(lrbase_total_target)
    ksyn  <- kdeg * rbase_total_target

    kpv            <- exp(lkpv)
    kcv            <- exp(lkcv)
    beta_virus     <- exp(lbeta_virus)
    lambda_virus   <- exp(llambda_virus)
    rbase_infected <- exp(lrbase_infected)
    imax           <- exp(limax)
    ic50           <- exp(lic50)

    # ---- Quasi-equilibrium free-drug concentration (eq 5) ---------------
    ctot   <- central / vc
    qeDisc <- ctot - total_target - km
    Cc     <- 0.5 * (qeDisc + sqrt(qeDisc * qeDisc + 4 * km * ctot))

    # Peripheral distribution is driven by the FREE central concentration,
    # the only species available to leave the vascular space; the paper adds
    # the second compartment via the Table 1 footnote without printing the
    # modified eq 3. Documented in the vignette Assumptions section.
    cp <- peripheral1 / vp

    inhVirus <- imax * Cc / (ic50 + Cc)

    # ---- ODE system (eqs 1-4, 6-7; eq 3 extended per Table 1 footnote) --
    d/dt(depot_ip) <- -ka_ip * depot_ip                                            # eq 1
    d/dt(depot_sc) <- -ka_sc * depot_sc                                            # eq 2
    d/dt(central)  <- ka_ip * depot_ip + ka_sc * depot_sc -
      cl * Cc - (central - Cc * vc) * kint - q * Cc + q * cp                       # eq 3 + two-compartment exchange
    d/dt(peripheral1)  <- q * Cc - q * cp
    d/dt(total_target) <- ksyn - (kint - kdeg) * (ctot - Cc) - kdeg * total_target # eq 4
    total_target(0)    <- rbase_total_target                                       # eq 4, A4(0) = RTOT0

    d/dt(infected) <- kpv * infected * (1 - inhVirus) - kcv * infected             # eq 6
    infected(0)    <- rbase_infected                                               # Table 2 BASE_T
    d/dt(virus)    <- beta_virus * (lambda_virus * infected - virus)               # eq 7

    # ---- Observations ---------------------------------------------------
    Cc    ~ prop(propSd)
    virus ~ prop(propSd_virus)
  })
}
