Fan_2025_nb457trimer_human <- function() {
  description <- paste(
    "Human allometric projection (no human data). Mouse-to-human scaled",
    "quasi-equilibrium (QE) target-mediated drug disposition (TMDD) PK model",
    "coupled to the HIV-1 viral-dynamic PD model for the anti-CD4 trimeric",
    "nanobody Nb457-NbHSA-Nb457 (Fan 2025, Microbiol Spectr). The structure is",
    "identical to the fitted murine model Fan_2025_nb457trimer_mouse.R; only",
    "clearance and central volume are scaled. The paper scaled V3 with an",
    "allometric exponent of 1 and CL with an exponent of 0.85 from the single",
    "mouse species, giving a human central volume of 36.26 mL/kg and a human",
    "linear clearance of 10.21 mL/h/kg (Results, Extrapolation of the PK/PD",
    "Model to Humans); this file anchors those at a 70 kg reference adult",
    "(vc = 2.538 L, cl = 0.7147 L/h) and keeps the exponents as fixed WT",
    "covariate effects so other body weights can be simulated. Every remaining",
    "parameter - absorption rate constants, the CD4 target parameters (RTOT0,",
    "K_M, kint = kdeg) and all viral-dynamic parameters - is held at its mouse",
    "value, which is the paper's explicit assumption ('all other human",
    "CD4-targeted parameters were assumed to be equivalent to those observed",
    "in mice ... all the human viral dynamic-related parameters were assumed",
    "same as the mouse parameters because the experiments were performed in",
    "humanized mouse models'). The companion ibalizumab projection",
    "(Fan_2025_ibalizumab_human.R) is the only leg of the extrapolation the",
    "authors could check against clinical data; the nanobody projection is",
    "unvalidated, and the PD extrapolation could not be validated for either",
    "drug. Intended use is simulation of human SC dosing regimens such as the",
    "published 20 mg/kg once every 2 days."
  )
  reference <- paste(
    "Fan X, Cao K, Wu X, Yan X. Pharmacokinetic-pharmacodynamic modeling of a",
    "highly potent and broadly neutralizing anti-CD4 trimeric nanobody to",
    "inhibit HIV-1 infection. Microbiol Spectr. 2025;13(11):e00805-25.",
    "doi:10.1128/spectrum.00805-25. PMCID: PMC12584636.",
    "Structural equations from eqs 1-7 (pp. 3-4); human CL and V3 from the",
    "Results section 'Extrapolation of the PK/PD Model to Humans' (36.26 mL/kg",
    "and 10.21 mL/h/kg); allometric exponents from Materials and Methods",
    "'Extrapolation of the mouse model to humans'; all remaining parameters",
    "from Table 1 and Table 2 (Nb457 rows) unchanged."
  )
  vignette <- "Fan_2025_antiCD4_nanobody_HIV"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling covariate, reference 70 kg. The exponent on the central volume is fixed at 1 and the exponent on linear clearance at 0.85, both taken from Materials and Methods 'Extrapolation of the mouse model to humans'. The 70 kg reference weight is not printed in the paper; it is recovered from the paper's own arithmetic (a 0.02 kg mouse V3 of 0.73 mL scaled with exponent 1 gives 36.5 mL/kg against the reported 36.26 mL/kg, and CL of 0.69 mL/h scaled with exponent 0.85 gives 10.15 mL/h/kg against the reported 10.21 mL/h/kg).",
      source_name        = "Not a fitted covariate; the paper applies the scaling to a per-kilogram basis."
    )
  )

  compartmentData <- list(
    depot_ip     = list(analyte = "Nb457-NbHSA-Nb457", units = "mg", specimen = "administration site", verified = TRUE),
    depot_sc     = list(analyte = "Nb457-NbHSA-Nb457", units = "mg", specimen = "administration site", verified = TRUE),
    central      = list(analyte = "Nb457-NbHSA-Nb457 (total: free plus CD4-bound)", units = "mg", specimen = "serum", verified = TRUE),
    total_target = list(analyte = "CD4 receptor (total, free plus drug-bound; carried as a drug-equivalent concentration)", units = "ug/mL", specimen = "blood cell", verified = FALSE),
    infected     = list(analyte = "HIV-1-infected CD4+ T cells", units = "cells/mL", specimen = "whole blood", verified = FALSE),
    virus        = list(analyte = "HIV-1 RNA", units = "copies/mL", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human (allometric projection from mouse; no human PK or PD data)",
    n_subjects     = NA_integer_,
    n_studies      = 0L,
    age_range      = "Not applicable; no human subjects contributed to this model.",
    weight_range   = "Not applicable; simulations are anchored at a 70 kg reference adult.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not applicable.",
    disease_state  = "Intended target population is HIV-1 infection; no human data were fitted.",
    dose_range     = "Simulated regimens only: 20 mg/kg SC once every 2 days (Fig. 4 and Fig. 5), plus 20 mg/kg once daily and 40 mg/kg once every 2 days in the Fig. 5 dose-intensity comparison, each for one month.",
    regions        = "Not applicable.",
    notes          = "This is a forward projection, not a fit. Only CL and the volumes were scaled; all CD4-target and viral-dynamic parameters were carried over from the murine fit unchanged. The authors validated the scaling procedure on ibalizumab, for which clinical PK were available (Fig. S7 and S8), and state that the PD extrapolation could not be validated for either drug because human viral-dynamic data under these scenarios do not exist. The Discussion also flags that interspecies differences in albumin binding by the NbHSA arm were not modelled, and that CD4 downregulation in chronic HIV-1 infection would reduce target-mediated clearance in patients."
  )

  ini({
    # ---- Absorption - carried over unchanged from the murine fit --------
    lka_ip <- log(0.46)
    label("First-order absorption rate constant after intraperitoneal injection (1/h)")  # Table 1, Nb457 row "k a (IP) (1/h)" = 0.46; not scaled
    lka_sc <- log(0.37)
    label("First-order absorption rate constant after subcutaneous injection (1/h)")     # Table 1, Nb457 row "k a (SC) (1/h)" = 0.37; not scaled

    # ---- Allometrically scaled disposition ------------------------------
    # 10.21 mL/h/kg * 70 kg = 714.7 mL/h = 0.7147 L/h
    lcl <- log(0.7147)
    label("Linear clearance of free drug for a 70 kg adult (L/h)")  # Results, human extrapolation: "human linear clearance of 10.21 mL/h/kg"
    # 36.26 mL/kg * 70 kg = 2538.2 mL = 2.5382 L
    lvc <- log(2.5382)
    label("Central volume of distribution for a 70 kg adult (L)")   # Results, human extrapolation: "scaled human central compartment volume of distribution was 36.26 mL/kg"

    e_wt_cl <- fixed(0.85)
    label("Allometric exponent on linear clearance for log(WT/70) (unitless)")  # Methods, human extrapolation: "the exponent for linear clearance (CL) was set to 0.85"
    e_wt_vc <- fixed(1)
    label("Allometric exponent on central volume for log(WT/70) (unitless)")    # Methods, human extrapolation: "The allometric exponent for the volume of distribution (V3) was set to 1"

    # ---- Target binding and turnover - assumed equal to the mouse -------
    lkm <- log(0.14)
    label("Quasi-equilibrium drug-receptor dissociation constant K_M (ug/mL)")  # Table 1, Nb457 row "K M (ug/mL)" = 0.14; assumed unchanged in humans
    lrbase_total_target <- log(174.7)
    label("Baseline total CD4 receptor RTOT0, in drug-equivalent concentration (ug/mL)")  # Table 1, Nb457 row "RTOT 0 (ug/mL)" = 174.7; assumed unchanged in humans
    lkint <- log(0.0027)
    label("Drug-receptor complex internalization rate constant kint, also used for kdeg (1/h)")  # Table 1, Nb457 row "k int (1/h)" = 0.0027; assumed unchanged in humans

    # ---- Viral dynamics - assumed equal to the humanized-mouse fit ------
    lkpv <- log(0.041)
    label("Viral infection / infected-cell expansion rate constant kPV (1/h)")  # Table 2, Nb457 row "k PV (1/h)" = 0.041
    lkcv <- log(0.025)
    label("Death rate constant of infected cells kCV (1/h)")                    # Table 2, Nb457 row "k CV (1/h)" = 0.025
    lbeta_virus <- log(8.12e-7)
    label("First-order rate constant for HIV-1 release from infected cells beta (1/h)")  # Table 2, Nb457 row "beta (1/h)" = 8.12e-7
    llambda_virus <- log(0.81)
    label("Scaling ratio lambda relating the infected-cell pool to viral load (unitless)")  # Table 2, Nb457 row "lambda" = 0.81
    # Table 2's units cell reads cells/mL, kept verbatim here. Note that
    # this value in cells/mL is 0.8 cells/uL, far below any viable CD4+
    # count, while in cells/uL it is unremarkable; the units label is
    # very likely off by 1000. Because eqs 6-7 are linear in this state,
    # `virus` is exactly proportional to BASE_T, and the vignette Errata
    # shows that a x1000 rescale closes most of the gap to the published
    # viral-load figures. Not applied: the error is inferred, not stated.
    lrbase_infected <- log(814.2)
    label("Initial condition of the infected CD4+ T-cell pool BASE_T (cells/mL)")  # Table 2, Nb457 row "BASE T (cells/mL)" = 814.2

    limax <- log(1.963)
    label("Maximum inhibitory effect on HIV-1 infection Imax (unitless)")  # Table 2, Nb457 row "I max" = 1.963
    lic50 <- fixed(log(0.0341))
    label("Free concentration producing half the maximum inhibition, SC50 in eq 6 (ug/mL)")  # Table 2, Nb457 row "IC 50 (ug/mL)" = 0.0341, FIXED (footnote c)

    # ---- Between-subject variability, carried over from the murine fit --
    # Table 1 footnote a establishes that the tabulated omega values are
    # NONMEM VARIANCES. Carrying murine IIV into a human projection is an
    # assumption of this packaged file, not a published human estimate; use
    # rxode2::zeroRe() for the typical-value projection the paper plotted.
    etalcl ~ 0.87    # Table 1, Nb457 row omega CL = 0.87 (variance)
    etalvc ~ 0.082   # Table 1, Nb457 row omega V3 = 0.082 (variance)

    # ---- Residual error, carried over from the murine fit ---------------
    propSd <- sqrt(0.82)
    label("Proportional residual error SD on serum concentration (fraction)")  # Table 1, Nb457 row "sigma 1" = 0.82 (variance) -> SD = sqrt(0.82)
    addSd <- sqrt(0.04)
    label("Additive residual error SD on serum concentration (ug/mL)")         # Table 1, Nb457 row "sigma 2" = 0.04 (variance) -> SD = sqrt(0.04)
    propSd_virus <- sqrt(0.21)
    label("Proportional residual error SD on viral load (fraction)")           # Table 2, Nb457 row "sigma 3" = 0.21 (variance) -> SD = sqrt(0.21)
  })

  model({
    # ---- Individual parameters ------------------------------------------
    ka_ip <- exp(lka_ip)
    ka_sc <- exp(lka_sc)
    cl    <- exp(lcl + e_wt_cl * log(WT / 70) + etalcl)
    vc    <- exp(lvc + e_wt_vc * log(WT / 70) + etalvc)
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
    # central is in mg and vc in L, so ctot and Cc are in mg/L, numerically
    # identical to the ug/mL used by the target parameters.
    ctot   <- central / vc
    qeDisc <- ctot - total_target - km
    Cc     <- 0.5 * (qeDisc + sqrt(qeDisc * qeDisc + 4 * km * ctot))

    inhVirus <- imax * Cc / (ic50 + Cc)

    # ---- ODE system (eqs 1-4, 6-7) --------------------------------------
    d/dt(depot_ip) <- -ka_ip * depot_ip                                            # eq 1
    d/dt(depot_sc) <- -ka_sc * depot_sc                                            # eq 2
    d/dt(central)  <- ka_ip * depot_ip + ka_sc * depot_sc -
      cl * Cc - (central - Cc * vc) * kint                                         # eq 3
    d/dt(total_target) <- ksyn - (kint - kdeg) * (ctot - Cc) - kdeg * total_target # eq 4
    total_target(0)    <- rbase_total_target                                       # eq 4, A4(0) = RTOT0

    d/dt(infected) <- kpv * infected * (1 - inhVirus) - kcv * infected             # eq 6
    infected(0)    <- rbase_infected                                               # Table 2 BASE_T
    d/dt(virus)    <- beta_virus * (lambda_virus * infected - virus)               # eq 7

    # ---- Observations ---------------------------------------------------
    Cc    ~ prop(propSd) + add(addSd)
    virus ~ prop(propSd_virus)
  })
}
