Liang_2024_rituximab <- function() {
  description <- "Two-compartment quasi-steady-state target-mediated drug disposition (TMDD) population PK/PD model for rituximab in adults with primary membranous nephropathy, with CD20 as the target; outputs are rituximab concentration and total CD20 (the molar surrogate for CD20+ B cell count) (Liang 2024)"
  reference <- "Liang H, Deng Z, Niu S, Kong W, Liu Y, Wang S, Li H, Wang Y, Zheng D, Liu D. Dosing optimization of rituximab for primary membranous nephropathy by population pharmacokinetic and pharmacodynamic study. Front Pharmacol. 2024 Mar 26;15:1197651. doi:10.3389/fphar.2024.1197651"
  vignette <- "Liang_2024_rituximab"

  # The authors converted every observation to molar units before fitting: RTX
  # concentrations to molar units, and CD20+ B cell counts to a CD20 molar
  # concentration by multiplying 94000 CD20 molecules per cell and dividing by
  # Avogadro's constant (Liang 2024 Section 2.3). Doses must therefore be
  # supplied in umol; see the vignette for the mg -> umol conversion.
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  compartmentData <- list(
    # RTX was assayed by ELISA on blood samples and the Discussion refers to
    # "serum RTX concentration"; CD20+ B cells were counted in peripheral blood
    # (Liang 2024 Sections 2.2, 4).
    central      = list(analyte = "rituximab (total: free + CD20-bound)", units = "umol",   specimen = "serum",       verified = TRUE),
    peripheral1  = list(analyte = "rituximab (free)",                     units = "umol",   specimen = "tissue",      verified = TRUE),
    total_target = list(analyte = "CD20 (total: free + rituximab-bound)", units = "umol/L", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list()

  # More than 30 baseline physiological parameters were screened (Liang 2024
  # Section 3.4); age, sex and eGFR are named explicitly in Section 2.4, and
  # body weight / sex / body surface area are discussed as covariates reported
  # for rituximab in other diseases. None was retained: "no covariates were
  # identified in our study. Thus, the final model was replaced by the base
  # structural model" (Liang 2024 Section 3.4). They are recorded here as
  # documentation of the covariate screen, not as model inputs.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      notes = "Screened in the SCM covariate analysis (Liang 2024 Section 2.4) but not retained in the final model (Section 3.4). Population range 19-76 years (Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "binary",
      type = "categorical",
      notes = "Screened in the SCM covariate analysis (Liang 2024 Section 2.4) but not retained (Section 3.4); the Discussion notes sex has been identified as a clearance covariate for rituximab in other diseases but had little influence in PMN. 10/41 (24.4%) female (Table 1)."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (BSA-normalized renal function)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Screened in the SCM covariate analysis (Liang 2024 Section 2.4) but not retained (Section 3.4). Baseline 92.0 +/- 24.6 mL/min/1.73 m^2, range 32-151 (Table 1)."
    ),
    WT = list(
      description = "Body weight at baseline",
      units = "kg",
      type = "continuous",
      notes = "Discussed as a clearance covariate reported for rituximab in other diseases but with 'little influence on RTX PK profile' in PMN (Liang 2024 Discussion); not retained (Section 3.4). Population 75.4 +/- 11.3 kg, range 50-101 (Table 1)."
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Discussed as a central-volume covariate reported for rituximab in other diseases but with 'little influence on RTX PK profile' in PMN (Liang 2024 Discussion); not retained (Section 3.4). BSA is not tabulated in the paper."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 41L,
    n_studies      = 1L,
    age_range      = "19-76 years (mean 52.8 +/- 14.9)",
    weight_range   = "50-101 kg (mean 75.4 +/- 11.3)",
    sex_female_pct = 24.4,
    disease_state  = "Adults with primary membranous nephropathy (PMN) treated off-label with rituximab; baseline urine protein 8.0 +/- 3.6 g/day, serum albumin 25.8 +/- 5.7 g/L, anti-PLA2R antibody titer 260.3 +/- 453.2 U/mL (range 5.4-2695), 17/41 (41.5%) with high titer (>150 U/mL).",
    renal_function = "eGFR 92.0 +/- 24.6 mL/min/1.73 m^2 (range 32-151); eGFR remained relatively stable through 12 months of treatment.",
    dose_range     = "Most patients received a monthly mini-dose of 100 mg IV rituximab; 10/41 received 200-500 mg for some months. Cumulative dose 578 +/- 350 mg at month 6, 923 +/- 544 mg at month 12 and 1060 +/- 705 mg at last follow-up (Table 1).",
    regions        = "Single centre: Department of Nephrology, Peking University Third Hospital, Beijing, China; retrospective cohort treated March 2019 to December 2021. Registered as ChiCTR2200057381.",
    notes          = paste(
      "Retrospective study; 171 rituximab concentrations, 220 CD20+ B cell counts and 276 anti-PLA2R titers",
      "from 41 patients (Liang 2024 Section 3.2). Median follow-up 15.9 months (range 6-44).",
      "Rituximab concentrations were measured by ELISA with a 3 ng/mL limit of quantification; observations",
      "below LOQ were excluded. Four structural PK models were compared and the QSS TMDD model had by far the",
      "lowest OFV (-908.3 versus 1293.4 for the linear two-compartment model; Supplementary Table 2), so it was",
      "selected. No covariates were retained, so the reported final model is the base structural model."
    )
  )

  # ---------------------------------------------------------------------------
  # Units of kss and ksyn in Table 2 (the one non-obvious reading in this model)
  #
  # Table 2 prints the unit of kss as "umol^-1" and of ksyn as "umol^-1 hr^-1".
  # Both are dimensionally impossible: a quasi-steady-state constant is a
  # concentration and a zero-order synthesis rate is concentration per time. The
  # receptor scale is pinned by the paper's own numbers -- ksyn / kdeg =
  # 5.06e-8 / 1.63e-3 = 3.104e-5, and 3.104e-5 umol/L back-converts through the
  # paper's own 94000 CD20 molecules per cell and Avogadro's constant to 199
  # CD20+ B cells/uL, which is exactly the baseline plateau of Figure 4. So the
  # slip is a dropped "L": the intended units are umol/L and umol/L/h.
  #
  # For the NUMBER of kss, 6.21 umol/L is ruled out by the paper's own
  # simulations: at that value the model depletes CD20 only to 33-125 cells/uL
  # and never reaches the <5/uL depletion criterion, whereas Section 3.7 reports
  # depletion within 24 h under all three regimens and Figure 4 shows a nadir of
  # about 2 cells/uL. Reading the printed value as 6.21 nmol/L = 6.21e-3 umol/L
  # reproduces the paper's own reported output on three independent anchors:
  #   * nadir 1.5-1.9 cells/uL and <5/uL by 24 h under every regimen (Figure 4,
  #     Section 3.7);
  #   * recovery to 90% of baseline 12.8 / 13.1 months after the last standard
  #     dose versus 9.6 months after the last mini-dose, against the reported
  #     "12 months vs 10 months" (Section 3.7);
  #   * the five monthly pre-dose troughs of Figure 3C (0.33, 0.48, 0.55, 0.58,
  #     0.59 ug/mL simulated versus about 0.40, 0.52, 0.57, 0.59, 0.60
  #     digitised).
  # 6.21 nmol/L is also the only reading consistent with rituximab's reported
  # low-nanomolar CD20 affinity. The literal reciprocal reading (Kss = 1/6.21 =
  # 0.161 umol/L) fails: it leaves the 100 mg mini-dose arm at a nadir of
  # 14-17 cells/uL, far above the depletion criterion.
  #
  # This reading is the best supported of the three but not perfect: it lets CD20+
  # B cells drift back to ~8-11 cells/uL at each monthly pre-dose trough, so
  # simulated depletion durations run 30-40% short of Section 3.7 (though the
  # ordering across regimens is reproduced). kss has NOT been adjusted to close
  # that gap. kss is also by far the least well determined parameter in the model
  # (60.1% RSE, 95% CI 3.09-15.5). See the vignette Errata for the full account.
  # ---------------------------------------------------------------------------

  ini({
    # Rituximab disposition (Liang 2024 Table 2).
    lcl   <- log(0.0482);         label("Nonspecific clearance (CL, L/h)")                                        # Liang 2024 Table 2: CL = 0.0482 L/h (95% CI 0.0336, 0.0681)
    lvc   <- log(2.48);           label("Central volume of distribution (V, L)")                                   # Liang 2024 Table 2: V = 2.48 L (1.77, 2.96)
    lq    <- log(0.0073);         label("Intercompartmental clearance (Q, L/h)")                                   # Liang 2024 Table 2: Q = 0.0073 L/h (0.000679, 0.0372)
    lvp   <- log(4.68);           label("Peripheral volume of distribution (V2, L)")                               # Liang 2024 Table 2: V2 = 4.68 L (1.86, 21.9)

    # Target binding and turnover (Liang 2024 Table 2, Section 2.3).
    # The paper's ktmd is the rituximab-CD20 complex elimination rate constant.
    # Section 2.3 states explicitly that "the internalization rate constant kint
    # was replaced by the RTX-CD20 complex target-mediated elimination rate
    # constant (ktmd)" because CD20 is non-internalizing and the complex is
    # cleared by ADCC / CDC, so ktmd occupies the kint slot of the Gibiansky
    # 2008 QSS equations and is named with the library canonical kint here.
    lkint <- fixed(log(0.217));   label("Rituximab-CD20 complex elimination rate constant (paper's ktmd, 1/h)")    # Liang 2024 Table 2: ktmd = 0.217 hr-1, held constant at the maximal mAb elimination rate of Glassman and Balthasar 2017 (Section 2.3)
    lkss  <- log(6.21e-3);        label("Quasi-steady-state constant (kss, umol/L)")                               # Liang 2024 Table 2: kss = 6.21 (3.09, 15.5); see the units note below -- 6.21 nmol/L = 6.21e-3 umol/L
    lksyn <- fixed(log(5.06e-8)); label("CD20 synthesis rate (ksyn, umol/L/h)")                                    # Liang 2024 Table 2: ksyn = 5.06e-8, held constant
    lkdeg <- fixed(log(1.63e-3)); label("CD20 degradation rate constant (kdeg, 1/h)")                              # Liang 2024 Table 2: kdeg = 1.63e-3 hr-1, held constant at the ~3.9%/day B cell disappearance rate of Macallan 2005 (Section 2.3)

    # Between-subject variability, exponential (Liang 2024 Eq. 1); the Omega
    # column of Table 2 reports variances. Every other parameter is "0 (Fixed)".
    etalcl  ~ 0.49                                                                                                # Liang 2024 Table 2: IIV on CL = 0.49 (0.307, 0.738), 16.1% RSE, 5.3% shrinkage
    etalkss ~ 2.207                                                                                               # Liang 2024 Table 2: IIV on kss = 2.207 (1.498, 2.805), 24.7% RSE, 29.1% shrinkage

    # Residual error was proportional for both outputs (Liang 2024 Eq. 2 and
    # Supplementary Table 2), but no sigma estimate is reported anywhere in the
    # paper or the Supplementary Material, so both are encoded as zero rather
    # than invented. See the vignette Errata.
    propSd              <- fixed(0); label("Proportional residual standard deviation for rituximab concentration (fraction)")
    propSd_total_target <- fixed(0); label("Proportional residual standard deviation for total CD20 (fraction)")
  })

  model({
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc)
    q    <- exp(lq)
    vp   <- exp(lvp)
    kint <- exp(lkint)
    kss  <- exp(lkss + etalkss)
    ksyn <- exp(lksyn)
    kdeg <- exp(lkdeg)

    # CD20 starts at its drug-free steady state, ksyn / kdeg = 3.104e-5 umol/L,
    # which back-converts to 199 CD20+ B cells/uL - the baseline plateau of
    # Liang 2024 Figure 4.
    total_target(0) <- ksyn / kdeg

    # Quasi-steady-state TMDD (Gibiansky 2008), as described in Liang 2024
    # Section 2.3: free rituximab is eliminated by a first-order rate constant
    # (kel = CL/V), distributes to and from the peripheral compartment, and binds
    # CD20 to form a complex that is eliminated at kint (the paper's ktmd).
    # `central` carries the TOTAL (free + complex) drug amount in umol and
    # `total_target` the total (free + bound) CD20 concentration in umol/L; free
    # drug cfree is recovered from the QSS binding quadratic
    # ctot = cfree + total_target * cfree / (kss + cfree). Only free drug
    # distributes to the peripheral compartment.
    ctot  <- central / vc
    disc  <- ctot - total_target - kss
    cfree <- 0.5 * (disc + sqrt(disc * disc + 4 * kss * ctot))
    cplx  <- total_target * cfree / (kss + cfree)

    d/dt(central) <-
      -(cl / vc) * cfree * vc - (q / vc) * cfree * vc +
        (q / vp) * peripheral1 - kint * cplx * vc
    d/dt(peripheral1) <-
       (q / vc) * cfree * vc - (q / vp) * peripheral1
    d/dt(total_target) <-
       ksyn - kdeg * (total_target - cplx) - kint * cplx

    # Rituximab was assayed by ELISA (Liang 2024 Section 2.2), i.e. total drug.
    # Free and total differ by the complex concentration, which is at most
    # ksyn / kdeg = 3.1e-5 umol/L, about 0.01% of the peak concentration, so the
    # distinction is numerically immaterial here.
    Cc <- central / vc

    Cc              ~ prop(propSd)
    total_target    ~ prop(propSd_total_target)
  })
}
