Fan_2025_ibalizumab_human <- function() {
  description <- paste(
    "Human allometric projection, PK-validated. Mouse-to-human scaled",
    "quasi-equilibrium (QE) target-mediated drug disposition (TMDD) PK model",
    "with two-compartment disposition, coupled to the HIV-1 viral-dynamic PD",
    "model, for ibalizumab (anti-CD4 humanized IgG4) as scaled in Fan 2025,",
    "Microbiol Spectr. The structure is identical to the fitted murine model",
    "Fan_2025_ibalizumab_mouse.R; only clearances and volumes are scaled",
    "(exponent 1 on volumes, 0.85 on clearances). The paper prints the two",
    "scaled central parameters directly - a human central volume of",
    "26.9 mL/kg and a human linear clearance of 0.35 mL/h/kg - anchored here",
    "at a 70 kg reference adult (vc = 1.883 L, cl = 0.0245 L/h). The",
    "peripheral parameters Q and V4 are NOT printed for the human; they are",
    "derived here by applying the same published exponents to the murine",
    "values (Q = 1.86 mL/h * (70/0.02)^0.85 = 1.915 L/h;",
    "V4 = 0.12 mL * (70/0.02)^1 = 0.420 L). All CD4 target parameters (RTOT0,",
    "K_M, kint = kdeg) and all viral-dynamic parameters are held at their",
    "murine values per the paper's explicit assumption. This is the leg of",
    "the extrapolation the authors validated: predicted concentrations after",
    "the published 10 mg/kg and 15 mg/kg IV regimens fell almost entirely",
    "within a two-fold error margin of observed clinical data (Fig. S7,",
    "Fig. S8). The PD extrapolation is unvalidated for both drugs. The",
    "underlying murine target parameters are very poorly identified (RTOT0",
    "%RSE 249.6, K_M 308.5, kint 1121), so the target-mediated limb of this",
    "projection should be treated as weakly determined."
  )
  reference <- paste(
    "Fan X, Cao K, Wu X, Yan X. Pharmacokinetic-pharmacodynamic modeling of a",
    "highly potent and broadly neutralizing anti-CD4 trimeric nanobody to",
    "inhibit HIV-1 infection. Microbiol Spectr. 2025;13(11):e00805-25.",
    "doi:10.1128/spectrum.00805-25. PMCID: PMC12584636.",
    "Structural equations from eqs 1-7 (pp. 3-4) plus the Table 1 footnote",
    "declaring a two-compartment ibalizumab disposition; human CL and V3 from",
    "the Results section 'Extrapolation of the PK/PD Model to Humans'",
    "(26.9 mL/kg and 0.35 mL/h/kg); allometric exponents from Materials and",
    "Methods; Q and V4 derived by the same exponents from the Table 1",
    "Ibalizumab murine values; all remaining parameters from Table 1 and",
    "Table 2 (Ibalizumab rows) unchanged. Clinical validation regimens are",
    "summarised in supplementary Table S2."
  )
  vignette <- "Fan_2025_antiCD4_nanobody_HIV"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling covariate, reference 70 kg. Exponent 1 on both volumes and 0.85 on both clearances, per Materials and Methods 'Extrapolation of the mouse model to humans'. The 70 kg reference weight is not printed; it is recovered from the paper's own arithmetic (a 0.02 kg mouse V3 of 0.54 mL scaled with exponent 1 gives 27.0 mL/kg against the reported 26.9 mL/kg, and CL of 0.024 mL/h scaled with exponent 0.85 gives 0.353 mL/h/kg against the reported 0.35 mL/h/kg).",
      source_name        = "Not a fitted covariate; the paper applies the scaling to a per-kilogram basis, and the validation regimens are themselves mg/kg."
    )
  )

  compartmentData <- list(
    depot_ip     = list(analyte = "ibalizumab", units = "mg", specimen = "administration site", verified = TRUE),
    depot_sc     = list(analyte = "ibalizumab", units = "mg", specimen = "administration site", verified = TRUE),
    central      = list(analyte = "ibalizumab (total: free plus CD4-bound)", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1  = list(analyte = "ibalizumab", units = "mg", specimen = "tissue", verified = FALSE),
    total_target = list(analyte = "CD4 receptor (total, free plus drug-bound; carried as a drug-equivalent concentration)", units = "ug/mL", specimen = "blood cell", verified = FALSE),
    infected     = list(analyte = "HIV-1-infected CD4+ T cells", units = "cells/mL", specimen = "whole blood", verified = FALSE),
    virus        = list(analyte = "HIV-1 RNA", units = "copies/mL", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human (allometric projection from mouse; PK validated against published clinical data)",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    age_range      = "Not reported in Fan 2025; the validation PK were taken from a previously published clinical study (reference 14, summarised in Table S2).",
    weight_range   = "Not reported; simulations are anchored at a 70 kg reference adult and the clinical regimens are mg/kg.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported.",
    disease_state  = "Multidrug-resistant HIV-1 infection (the approved ibalizumab indication); the validation data are trough and profile concentrations from the cited clinical study.",
    dose_range     = "Validation regimens (Table S2, IV infusion): 10 mg/kg weekly for the first 9 doses then 10 mg/kg every 2 weeks; and 15 mg/kg alternating weekly with placebo for the first 9 doses (up to week 8) then 15 mg/kg every 2 weeks. Trough samples were collected at all visits from day 1, weekly through week 4, then every 2 weeks through week 48.",
    regions        = "Not reported.",
    notes          = "No human data were fitted; the human parameters are a forward allometric projection from the murine fit, checked against the published clinical concentrations. 'Nearly all model-estimated concentrations varied within a twofold error margin compared to the observed data' (Results). The PD layer was not validated in humans because human viral-dynamic data for these scenarios were unavailable."
  )

  ini({
    # ---- Absorption - carried over unchanged from the murine fit --------
    # Retained so the packaged structure matches the murine model; the human
    # validation and simulations use intravenous dosing into central.
    lka_ip <- log(0.26)
    label("First-order absorption rate constant after intraperitoneal injection (1/h)")  # Table 1, Ibalizumab row "k a (IP) (1/h)" = 0.26; not scaled
    lka_sc <- log(0.22)
    label("First-order absorption rate constant after subcutaneous injection (1/h)")     # Table 1, Ibalizumab row "k a (SC) (1/h)" = 0.22; not scaled

    # ---- Allometrically scaled two-compartment disposition --------------
    # 0.35 mL/h/kg * 70 kg = 24.5 mL/h = 0.0245 L/h
    lcl <- log(0.0245)
    label("Linear clearance of free drug for a 70 kg adult (L/h)")  # Results, human extrapolation: "a projected human linear clearance of 0.35 mL/h/kg"
    # 26.9 mL/kg * 70 kg = 1883 mL = 1.883 L
    lvc <- log(1.883)
    label("Central volume of distribution for a 70 kg adult (L)")   # Results, human extrapolation: "predicting a human central compartment volume of distribution of 26.9 mL/kg"
    # Derived, not printed: 1.86 mL/h * (70/0.02)^0.85 = 1914.96 mL/h
    lq  <- log(1.915)
    label("Inter-compartmental clearance for a 70 kg adult (L/h)")   # Derived: Table 1 Ibalizumab Q = 1.86 mL/h scaled with the published CL exponent 0.85
    # Derived, not printed: 0.12 mL * (70/0.02)^1 = 420 mL
    lvp <- log(0.420)
    label("Peripheral volume of distribution V4 for a 70 kg adult (L)")  # Derived: Table 1 Ibalizumab V4 = 0.12 mL scaled with the published volume exponent 1

    e_wt_cl <- fixed(0.85)
    label("Allometric exponent on linear clearance for log(WT/70) (unitless)")           # Methods: "the exponent for linear clearance (CL) was set to 0.85"
    e_wt_vc <- fixed(1)
    label("Allometric exponent on central volume for log(WT/70) (unitless)")             # Methods: "The allometric exponent for the volume of distribution (V3) was set to 1"
    e_wt_q  <- fixed(0.85)
    label("Allometric exponent on inter-compartmental clearance for log(WT/70) (unitless)")  # Assumed equal to the published CL exponent; Q is not named in the Methods
    e_wt_vp <- fixed(1)
    label("Allometric exponent on peripheral volume for log(WT/70) (unitless)")          # Assumed equal to the published volume exponent; V4 is not named in the Methods

    # ---- Target binding and turnover - assumed equal to the mouse -------
    lkm <- log(12.73)
    label("Quasi-equilibrium drug-receptor dissociation constant K_M (ug/mL)")  # Table 1, Ibalizumab row "K M (ug/mL)" = 12.73; assumed unchanged in humans
    lrbase_total_target <- log(11.73)
    label("Baseline total CD4 receptor RTOT0, in drug-equivalent concentration (ug/mL)")  # Table 1, Ibalizumab row "RTOT 0 (ug/mL)" = 11.73; assumed unchanged in humans
    lkint <- log(0.0036)
    label("Drug-receptor complex internalization rate constant kint, also used for kdeg (1/h)")  # Table 1, Ibalizumab row "k INT (1/h)" = 0.0036; assumed unchanged in humans

    # ---- Viral dynamics - assumed equal to the humanized-mouse fit ------
    lkpv <- log(0.027)
    label("Viral infection / infected-cell expansion rate constant kPV (1/h)")  # Table 2, Ibalizumab row "k PV (1/h)" = 0.027
    lkcv <- log(0.0072)
    label("Death rate constant of infected cells kCV (1/h)")                    # Table 2, Ibalizumab row "k CV (1/h)" = 0.0072
    lbeta_virus <- log(2.09e-7)
    label("First-order rate constant for HIV-1 release from infected cells beta (1/h)")  # Table 2, Ibalizumab row "beta (1/h)" = 2.09e-7
    llambda_virus <- log(2.001)
    label("Scaling ratio lambda relating the infected-cell pool to viral load (unitless)")  # Table 2, Ibalizumab row "lambda" = 2.001
    lrbase_infected <- log(1225)
    label("Initial condition of the infected CD4+ T-cell pool BASE_T (cells/mL)")  # Table 2, Ibalizumab row "BASE T" = 1,225

    limax <- log(0.42)
    label("Maximum inhibitory effect on HIV-1 infection Imax (unitless)")  # Table 2, Ibalizumab row "I max" = 0.42
    lic50 <- fixed(log(0.1449))
    label("Free concentration producing half the maximum inhibition, SC50 in eq 6 (ug/mL)")  # Table 2, Ibalizumab row "SC 50 (ug/mL)" = 0.1449, FIXED (footnote c)

    # ---- Between-subject variability, carried over from the murine fit --
    etalcl ~ 0.29   # Table 1, Ibalizumab row omega CL = 0.29 (variance)
    etalvc ~ 0.68   # Table 1, Ibalizumab row omega V3 = 0.68 (variance)

    # ---- Residual error, carried over from the murine fit ---------------
    propSd <- sqrt(0.053)
    label("Proportional residual error SD on serum concentration (fraction)")  # Table 1, Ibalizumab row "sigma 1" = 0.053 (variance) -> SD = sqrt(0.053)
    propSd_virus <- sqrt(0.94)
    label("Proportional residual error SD on viral load (fraction)")           # Table 2, Ibalizumab row "sigma 3" = 0.94 (variance) -> SD = sqrt(0.94)
  })

  model({
    # ---- Individual parameters ------------------------------------------
    ka_ip <- exp(lka_ip)
    ka_sc <- exp(lka_sc)
    cl    <- exp(lcl + e_wt_cl * log(WT / 70) + etalcl)
    vc    <- exp(lvc + e_wt_vc * log(WT / 70) + etalvc)
    q     <- exp(lq  + e_wt_q  * log(WT / 70))
    vp    <- exp(lvp + e_wt_vp * log(WT / 70))
    km    <- exp(lkm)
    kint  <- exp(lkint)
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
    cp     <- peripheral1 / vp

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
