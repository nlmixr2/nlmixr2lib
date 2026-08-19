Fan_2025_nb457trimer_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). Mechanism-based quasi-equilibrium (QE)",
    "target-mediated drug disposition (TMDD) PK model coupled to an HIV-1",
    "viral-dynamic PD model for the anti-CD4 trimeric nanobody",
    "Nb457-NbHSA-Nb457 (about 40 kDa) in HIV-1 CH058-infected humanized",
    "NDG-HuPBL mice (Fan 2025, Microbiol Spectr). Disposition is",
    "one-compartment with parallel intraperitoneal and subcutaneous",
    "first-order absorption depots (separate ka per route, no bioavailability",
    "term); the central state carries TOTAL drug (free plus CD4-bound) and the",
    "free concentration Cc is recovered at every integration step from the QE",
    "quadratic (eq 5). Free drug is cleared linearly (CL) and the",
    "drug-receptor complex is internalized at kint. Total CD4 receptor",
    "(total_target, expressed in drug-equivalent ug/mL) turns over with",
    "zero-order synthesis ksyn and first-order degradation kdeg; the authors",
    "fixed kdeg = kint to avoid over-parameterization, which makes total",
    "receptor constant at its baseline RTOT0 and ksyn = kint * RTOT0. The PD",
    "layer is a two-state viral-dynamic model: an HIV-1-infected CD4+ T-cell",
    "pool (infected) that expands at kPV and is lost at kCV, with the drug",
    "acting as an Imax inhibition of the kPV expansion term, feeding a free",
    "virus pool (virus) through a slow first-order release constant beta that",
    "drives virus toward lambda * infected. Because the fitted Imax = 1.963",
    "exceeds 1, sustained free concentrations above IC50 / (Imax - 1) =",
    "0.0354 ug/mL drive net loss of the infected pool rather than merely",
    "slowing its growth. PK and PD were fitted sequentially by a naive pooled",
    "approach (n = 4 mice per arm); see the validation vignette Errata for the",
    "reproduction gap between these published parameters and the published",
    "PD figures."
  )
  reference <- paste(
    "Fan X, Cao K, Wu X, Yan X. Pharmacokinetic-pharmacodynamic modeling of a",
    "highly potent and broadly neutralizing anti-CD4 trimeric nanobody to",
    "inhibit HIV-1 infection. Microbiol Spectr. 2025;13(11):e00805-25.",
    "doi:10.1128/spectrum.00805-25. PMCID: PMC12584636.",
    "Structural equations from eqs 1-7 (pp. 3-4); PK parameter values from",
    "Table 1 (Nb457 rows); PD parameter values from Table 2 (Nb457 rows).",
    "The 19 September 2025 version was corrected on 4 November 2025 for an",
    "error in the Acknowledgments only; no model value is affected."
  )
  vignette <- "Fan_2025_antiCD4_nanobody_HIV"
  units <- list(time = "h", dosing = "ug", concentration = "ug/mL")

  compartmentData <- list(
    depot_ip     = list(analyte = "Nb457-NbHSA-Nb457", units = "ug", specimen = "administration site", verified = TRUE),
    depot_sc     = list(analyte = "Nb457-NbHSA-Nb457", units = "ug", specimen = "administration site", verified = TRUE),
    central      = list(analyte = "Nb457-NbHSA-Nb457 (total: free plus CD4-bound)", units = "ug", specimen = "serum", verified = TRUE),
    total_target = list(analyte = "CD4 receptor (total, free plus drug-bound; carried as a drug-equivalent concentration)", units = "ug/mL", specimen = "blood cell", verified = FALSE),
    infected     = list(analyte = "HIV-1-infected CD4+ T cells", units = "cells/mL", specimen = "whole blood", verified = FALSE),
    virus        = list(analyte = "HIV-1 RNA", units = "copies/mL", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "mouse (NDG-HuPBL humanized, female)",
    n_subjects     = 8L,
    n_studies      = 1L,
    age_range      = "Not reported.",
    weight_range   = "Not reported. A reference body weight of 0.02 kg is implied by the paper's own allometric projection (see notes) but is never printed.",
    sex_female_pct = 100,
    race_ethnicity = "Not applicable (murine study).",
    disease_state  = "HIV-1 CH058 infection established by challenge with 10 ng p24 on day 0.",
    dose_range     = "400 ug per mouse. PK arms: a single 400 ug intraperitoneal (IP) OR subcutaneous (SC) dose on day 1 post-infection, blood sampled at 0, 1, 4, 8, 12, 24, 72 and 120 h. PD arm: 400 ug IP on day 1 followed by 400 ug SC on days 3, 5 and 7, viral load sampled at weeks 1, 2, 3 and 4 post-infection.",
    regions        = "Nanjing, People's Republic of China.",
    notes          = "n = 4 female mice per arm for the PK study (one IP arm, one SC arm) and n = 4 for the PD arm (Materials and Methods, Data; Table S1). Because of the limited blood volume in mice, PK and PD were studied in separate animals, and the model was fitted by a naive pooled-data approach in which all individuals are treated as one subject (Materials and Methods, PD model). Despite that statement, Table 1 still reports estimated IIV variances for CL and V3, which this model file carries; see the vignette Errata. Ibalizumab was run as a positive control on the identical regimen and is packaged separately as Fan_2025_ibalizumab_mouse.R. Body weight is never printed; 0.02 kg is recoverable from the paper's own allometric projections (V3 = 0.73 mL scaled to 36.26 mL/kg, and 0.54 mL scaled to 26.9 mL/kg for ibalizumab, both with an exponent of 1)."
  )

  ini({
    # ---- Absorption (Table 1, Nb457 rows) -------------------------------
    lka_ip <- log(0.46)
    label("First-order absorption rate constant after intraperitoneal injection (1/h)")  # Table 1, Nb457 row "k a (IP) (1/h)" = 0.46 (%RSE 2.02)
    lka_sc <- log(0.37)
    label("First-order absorption rate constant after subcutaneous injection (1/h)")     # Table 1, Nb457 row "k a (SC) (1/h)" = 0.37 (%RSE 7.03)

    # ---- Linear disposition (Table 1, Nb457 rows) -----------------------
    lcl <- log(0.69)
    label("Linear clearance of free drug from the central compartment (mL/h)")  # Table 1, Nb457 row "CL (mL/h)" = 0.69 (%RSE 33.21)
    lvc <- log(0.73)
    label("Volume of distribution of the central compartment (mL)")             # Table 1, Nb457 row "V 3 (mL)" = 0.73 (%RSE 12.49)

    # ---- Target binding and turnover (Table 1, Nb457 rows) --------------
    # K_M is the paper's symbol for the drug-receptor complex equilibrium
    # dissociation constant under the quasi-equilibrium assumption (eq 5 and
    # the text following eq 5); Table 1 labels the same row "Michaelis
    # constant".
    lkm <- log(0.14)
    label("Quasi-equilibrium drug-receptor dissociation constant K_M (ug/mL)")  # Table 1, Nb457 row "K M (ug/mL)" = 0.14 (%RSE 19.25)

    lrbase_total_target <- log(174.7)
    label("Baseline total CD4 receptor RTOT0, in drug-equivalent concentration (ug/mL)")  # Table 1, Nb457 row "RTOT 0 (ug/mL)" = 174.7 (%RSE 12.09)

    # kdeg is set equal to kint in model() - an explicit paper assumption
    # ("To simplify the model and reduce the risk of overparameterization,
    # kint was assumed to be equal to kdeg", text after eq 5), so only one
    # rate constant is estimated.
    lkint <- log(0.0027)
    label("Drug-receptor complex internalization rate constant kint, also used for free-receptor degradation kdeg (1/h)")  # Table 1, Nb457 row "k int (1/h)" = 0.0027 (%RSE 42.4)

    # ---- Viral dynamics (Table 2, Nb457 rows) ---------------------------
    lkpv <- log(0.041)
    label("Viral infection / infected-cell expansion rate constant kPV (1/h)")  # Table 2, Nb457 row "k PV (1/h)" = 0.041 (%RSE 14.1)
    lkcv <- log(0.025)
    label("Death rate constant of infected cells kCV (1/h)")                    # Table 2, Nb457 row "k CV (1/h)" = 0.025 (%RSE 22.72)
    lbeta_virus <- log(8.12e-7)
    label("First-order rate constant for HIV-1 release from infected cells beta (1/h)")  # Table 2, Nb457 row "beta (1/h)" = 8.12e-7 (%RSE 21.1)
    llambda_virus <- log(0.81)
    label("Scaling ratio lambda relating the infected-cell pool to viral load (unitless)")  # Table 2, Nb457 row "lambda" = 0.81 (%RSE 36.1)
    # Table 2's units cell reads cells/mL, kept verbatim here. Note that
    # this value in cells/mL is 0.8 cells/uL, far below any viable CD4+
    # count, while in cells/uL it is unremarkable; the units label is
    # very likely off by 1000. Because eqs 6-7 are linear in this state,
    # `virus` is exactly proportional to BASE_T, and the vignette Errata
    # shows that a x1000 rescale closes most of the gap to the published
    # viral-load figures. Not applied: the error is inferred, not stated.
    lrbase_infected <- log(814.2)
    label("Estimated initial condition of the infected CD4+ T-cell pool BASE_T (cells/mL)")  # Table 2, Nb457 row "BASE T (cells/mL)" = 814.2 (%RSE 44.2)

    # ---- Drug effect on viral dynamics (Table 2, Nb457 rows) ------------
    # Imax > 1 is what the paper estimated; it lets the bracket
    # (1 - Imax * Cc / (ic50 + Cc)) go negative, so that a sustained free
    # concentration above ic50 / (Imax - 1) produces net loss of the
    # infected pool. Kept exactly as published.
    limax <- log(1.963)
    label("Maximum inhibitory effect of Nb457-NbHSA-Nb457 on HIV-1 infection Imax (unitless)")  # Table 2, Nb457 row "I max" = 1.963 (%RSE 16.24)

    # Fixed, not estimated: "The IC50 values for each drug were fixed based on
    # experimentally observed data because it would be inaccurate to estimate
    # IC50 with only one dose level" (Results); Table 2 footnote c "The IC 50
    # was fixed to the observed value". 0.0341 ug/mL is the measured IC50
    # against the challenge strain HIV-1 CH058.
    lic50 <- fixed(log(0.0341))
    label("Free concentration producing half the maximum inhibition, SC50 in eq 6, set to the observed HIV-1 CH058 IC50 (ug/mL)")  # Table 2, Nb457 row "IC 50 (ug/mL)" = 0.0341, FIXED (footnote c)

    # ---- Between-subject variability (Table 1, Nb457 rows) --------------
    # Table 1 footnote a - "RSE for omega and sigma are reported on the
    # approximate S.D. scale (standard error/variance estimate)/2" - states
    # that the tabulated omega and sigma are NONMEM VARIANCES, so they are
    # used here as variances directly.
    etalcl ~ 0.87    # Table 1, Nb457 row omega CL = 0.87 (variance; %RSE 27.73)
    etalvc ~ 0.082   # Table 1, Nb457 row omega V3 = 0.082 (variance; %RSE 42.75)

    # ---- Residual error -------------------------------------------------
    # Combined proportional plus additive on the serum concentration
    # (Table 1 reports both a proportional and an additive sigma for
    # Nb457-NbHSA-Nb457). Table values are variances (footnote a), so the
    # nlmixr2 SD-scale parameters are their square roots.
    propSd <- sqrt(0.82)
    label("Proportional residual error SD on serum concentration (fraction)")  # Table 1, Nb457 row "sigma 1" = 0.82 (variance; %RSE 11.75) -> SD = sqrt(0.82)
    addSd <- sqrt(0.04)
    label("Additive residual error SD on serum concentration (ug/mL)")         # Table 1, Nb457 row "sigma 2" = 0.04 (variance; %RSE 10.2) -> SD = sqrt(0.04)
    propSd_virus <- sqrt(0.21)
    label("Proportional residual error SD on viral load (fraction)")           # Table 2, Nb457 row "sigma 3" = 0.21 (variance; %RSE 14.5) -> SD = sqrt(0.21)
  })

  model({
    # ---- Individual parameters ------------------------------------------
    ka_ip <- exp(lka_ip)
    ka_sc <- exp(lka_sc)
    cl    <- exp(lcl + etalcl)
    vc    <- exp(lvc + etalvc)
    km    <- exp(lkm)
    kint  <- exp(lkint)

    # Paper assumption (text after eq 5): kint == kdeg.
    kdeg  <- kint

    rbase_total_target <- exp(lrbase_total_target)

    # Zero-order receptor synthesis implied by the eq 4 initial condition
    # A4(0) = RTOT0 held at steady state: ksyn = kdeg * RTOT0. Not tabulated
    # separately in Table 1; it is fully determined by kdeg and RTOT0.
    ksyn  <- kdeg * rbase_total_target

    kpv           <- exp(lkpv)
    kcv           <- exp(lkcv)
    beta_virus    <- exp(lbeta_virus)
    lambda_virus  <- exp(llambda_virus)
    rbase_infected <- exp(lrbase_infected)
    imax          <- exp(limax)
    ic50          <- exp(lic50)

    # ---- Quasi-equilibrium free-drug concentration (eq 5) ---------------
    # ctot is total drug concentration in the central compartment; Cc is the
    # positive root of the QE binding quadratic and is the quantity fitted
    # against the observed serum ELISA data.
    ctot   <- central / vc
    qeDisc <- ctot - total_target - km
    Cc     <- 0.5 * (qeDisc + sqrt(qeDisc * qeDisc + 4 * km * ctot))

    # ---- Drug effect on infected-cell expansion (eq 6 bracket) ----------
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
    # virus(0) = 0 by default ("Initial conditions were estimated for T cells
    # and zero for HIV-1", text after eq 7).

    # ---- Observations ---------------------------------------------------
    Cc    ~ prop(propSd) + add(addSd)
    virus ~ prop(propSd_virus)
  })
}
