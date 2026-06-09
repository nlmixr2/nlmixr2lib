Kamal_2015_oseltamivir <- function() {
  description <- paste(
    "Mechanistic drug-disease (viral-dynamics) model of influenza-virus",
    "progression and oseltamivir antiviral effect in adults with",
    "experimental and naturally-acquired influenza A (H1N1) virus",
    "infection (Kamal 2015). Builds on the Baccam et al. (2006)",
    "target-cell-limited viral-dynamics framework: uninfected target",
    "respiratory epithelial cells (target_cells) are infected by free",
    "virus (virus) at second-order rate beta_inf; infected cells",
    "(infected_cells) produce virus at rate p_prod per cell per day and",
    "die at rate delta_clr; free virus is cleared at rate c_clr.",
    "Oseltamivir inhibits viral production through an inhibitory Hill",
    "function acting on log10(p) (Equation 4 of Kamal 2015), parameterised",
    "so Emax is the maximum log10-fold reduction of p and ED50 is the",
    "dose producing a 2-fold (50%) reduction of p on the linear scale.",
    "Dose enters via the per-record DOSE covariate (mg per administered",
    "oseltamivir dose; 0 during placebo or outside the treatment window);",
    "no oseltamivir pharmacokinetics are modelled. Initial conditions",
    "are fixed per Baccam et al. (2006): target_cells(0) = 4e8 epithelial",
    "cells (from a 160 cm^2 upper-respiratory-tract surface area and",
    "2e-11 to 4e-11 m^2 per epithelial cell), infected_cells(0) = 0, and",
    "virus(0) = 10^0.25 TCID50/mL (the viral-titer lower limit of",
    "quantification, used as the inoculation viral titer). The viral",
    "load viralLoad (TCID50/mL of nasal wash, canonical PD-output name)",
    "is the single observed output with proportional residual error,",
    "equivalent to the paper's log10-transformed additive-error model.",
    "The three viral-dynamics compartments are declared paper-specific",
    "(see paper_specific_compartments)."
  )
  reference <- paste(
    "Kamal MA, Gieschke R, Lemenuel-Diot A, Beauchemin CAA, Smith PF,",
    "Rayner CR. (2015).",
    "A drug-disease model describing the effect of oseltamivir",
    "neuraminidase inhibition on influenza virus progression.",
    "Antimicrob Agents Chemother 59(9):5388-5395.",
    "doi:10.1128/AAC.00069-15. PMID 26100711; PMCID PMC4538529.",
    sep = " "
  )
  vignette <- "Kamal_2015_oseltamivir"

  paper_specific_compartments <- c("target_cells", "infected_cells", "virus")

  units <- list(
    time          = "day",
    dosing        = "mg per administered oseltamivir dose (per-record DOSE covariate; 0 during placebo and outside the 5-day b.i.d. treatment window)",
    concentration = "TCID50/mL of nasal wash (viralLoad)"
  )

  covariateData <- list(
    DOSE = list(
      description        = "Per-record administered oseltamivir dose level (mg) driving the inhibitory Hill function on viral production rate p. Set to 0 during placebo arms or outside the 5-day b.i.d. treatment window. The source paper does not include an oseltamivir PK ODE; dose enters the PD model directly (Equation 4).",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate matching the canonical DOSE entry's use case (b): time-varying current administered dose feeding a derived exposure term without an explicit PK compartment. For simulation of the paper's Figure 4 dose-response curves, set DOSE = 75 (or 20, 100, 150, 200, etc.) during the 5-day treatment window and 0 elsewhere. In study PV15616 treatment started 28 h (= 1.167 day) after intranasal inoculation, giving a treatment window of t in [1.167, 6.167] day in absolute time after infection.",
      source_name        = "Dose"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 208L,
    n_studies      = 4L,
    age_range      = "Adults (specific age range not tabulated in the source paper).",
    age_median     = "Not reported in the source paper.",
    weight_range   = "Not reported in the source paper.",
    weight_median  = "Not reported in the source paper.",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Experimental human inoculation with influenza A virus (H1N1: A/Hong Kong/123/77 in study Baccam; A/Texas/36/91 in studies PV15616 and PV15615) or naturally-acquired influenza-virus infection (study WV15670). All subjects gave informed consent under each study's institutional-review-board approval.",
    dose_range     = "Placebo (152 subjects across all four studies) and oral oseltamivir at 20, 100, or 200 mg b.i.d. or 200 mg q.d. for 5 days (56 subjects in PV15616). Simulation exercises in the paper also explored the 75 and 150 mg b.i.d. clinical doses.",
    regions        = "Not reported in the source paper.",
    notes          = "Pooled across four studies per Table 1 of Kamal 2015. Placebo data: 573 positive viral-titer time points across all four studies; oseltamivir treatment data: 298 positive viral-titer time points from PV15616 only (PV15616 was the sole dose-ranging study with the wide dose range and dense viral-titer sampling required for PD-parameter estimation). Viral titer was sampled in nasal washings as 50% tissue culture infective dose per mL (TCID50/mL) on MDCK cells and was assumed proportional to free-virus concentration at the site of infection."
  )

  ini({
    # Influenza viral-dynamics structural parameters (Kamal 2015 Table 2).
    # beta_inf (target-cell infection rate) and p_prod (viral production
    # rate per infected cell) were modelled as 10^theta in NONMEM
    # (theta a fixed-effects parameter) for numerical stability because
    # their linear-scale values are on the order of 10^-4; the Table 2
    # entries are the back-transformed linear-scale estimates.
    lbeta_inf  <- log(7.41e-4); label("Target-cell infection rate beta ((TCID50/mL)^-1 * day^-1)")        # Table 2: beta = 7.41E-4, %SEM 10
    lp_prod    <- log(2.0e-4);  label("Viral production rate per infected cell p (TCID50/mL * day^-1)")   # Table 2: p = 2.0E-4, %SEM 9
    lc_clr     <- log(3.33);    label("Free-virus clearance rate c (1/day)")                              # Table 2: c = 3.33, %SEM 22
    ldelta_clr <- log(2.49);    label("Infected-cell clearance rate delta (1/day)")                       # Table 2: delta = 2.49, %SEM 28

    # Oseltamivir-effect Hill function on log10(p) (Kamal 2015 Equation
    # 4 and footnote). Emax is the maximum log10-fold inhibition of p;
    # ED50 is the dose producing 50% inhibition of viral production per
    # infected cell, i.e. the dose at which the linear-scale p drops to
    # half its untreated value (log10(p) drops by log10(2)). ED50 is the
    # clinically meaningful parameter reported in Table 2; the directly
    # fit NONMEM parameter ED50* is derived inside model() from
    # ED50 = ED50* / (Emax / log10(2) - 1).
    lemax <- log(2.35); label("Maximum log10-fold inhibition of viral production rate p by oseltamivir (Emax, log10 units)")  # Table 2: Emax = 2.35, %SEM 25
    led50 <- log(3.2);  label("Oseltamivir dose producing 50% (2-fold) reduction of viral production rate p (ED50, mg)")       # Table 2: ED50 = 3.2 mg, %SEM 69

    # Inter-individual variability. Table 2 reports IIV on p_prod
    # (65% CV) and Emax (82% CV); no IIV on beta_inf, c_clr, delta_clr,
    # or ED50. The exponential variability model Pj = TVP * exp(eta_j)
    # places eta on the natural log of the linear-scale parameter so
    # omega^2 = log(CV^2 + 1).
    etalp_prod ~ log(0.65^2 + 1)  # Table 2: IIV p   = 65% CV  -> omega^2 = log(1.4225) = 0.3525
    etalemax   ~ log(0.82^2 + 1)  # Table 2: IIV Emax = 82% CV -> omega^2 = log(1.6724) = 0.5142

    # Residual error. The paper fit viral titer on the log10 scale with
    # additive error, which corresponds to a proportional error model
    # on untransformed viralLoad (Kamal 2015 Materials and Methods,
    # last paragraph of the "Influenza model and oseltamivir
    # pharmacodynamics" section). Reported sigma_error = 14% CV.
    propSd <- 0.14; label("Proportional residual error on viral load viralLoad (fraction)")  # Table 2: sigma_error = 14% CV
  })

  model({
    # Individual-level structural parameters. Only p_prod and emax
    # carry IIV per Table 2; the other four structural parameters are
    # typical-value-only.
    beta_inf  <- exp(lbeta_inf)
    p_prod    <- exp(lp_prod + etalp_prod)
    c_clr     <- exp(lc_clr)
    delta_clr <- exp(ldelta_clr)
    emax      <- exp(lemax + etalemax)
    ed50      <- exp(led50)

    # Convert the reported ED50 (50% inhibition in linear space) to the
    # directly fit NONMEM parameter ED50* (Kamal 2015 Equation 4
    # footnote): ED50 = ED50* / (Emax / log10(2) - 1)  =>
    # ED50* = ED50 * (Emax / log10(2) - 1). For the reported point
    # estimates Emax = 2.35 and ED50 = 3.2 mg this gives ED50* approx
    # 21.78 mg.
    ed50star <- ed50 * (emax / log10(2) - 1)

    # Inhibitory Hill function on log10(p) (Kamal 2015 Equation 4):
    # p_eff = 10^[log10(p_prod) - Emax * DOSE / (DOSE + ED50*)].
    # When DOSE = 0 (placebo or pre/post-treatment), inh = 0 and
    # p_eff = p_prod; when DOSE >> ED50*, inh -> Emax so
    # p_eff -> p_prod / 10^Emax (maximum suppression of virus
    # production per infected cell). The small 1e-12 in the denominator
    # guards against a 0/0 indeterminate when DOSE = 0 and ED50* = 0.
    inh   <- emax * DOSE / (DOSE + ed50star + 1e-12)
    p_eff <- p_prod * 10^(-inh)

    # Viral-dynamics ODE system (Kamal 2015 Equations 1-3; Baccam et al.
    # 2006 framework). target_cells (T), infected_cells (I), and virus
    # (V) are declared as paper_specific_compartments above.
    d/dt(target_cells)   <- -beta_inf * target_cells * virus
    d/dt(infected_cells) <-  beta_inf * target_cells * virus - delta_clr * infected_cells
    d/dt(virus)          <-  p_eff * infected_cells - c_clr * virus

    # Initial conditions per Kamal 2015 Methods (Influenza model and
    # oseltamivir pharmacodynamics section), inherited from Baccam et
    # al. (2006). T0 = 4 * 10^8 epithelial cells (160 cm^2 nasal surface
    # area at 2e-11 to 4e-11 m^2 per epithelial cell); I0 = 0 (no
    # infected cells at inoculation); V0 = 10^0.25 TCID50/mL (the
    # viral-titer lower limit of quantification, used as the
    # inoculation viral titer).
    target_cells(0)   <- 4e8
    infected_cells(0) <- 0
    virus(0)          <- 10^0.25

    # Observation: free viral load (TCID50/mL of nasal wash; canonical
    # PD-output name viralLoad). The source paper fit log10(V) with
    # additive error which corresponds to a proportional error model
    # on the untransformed viralLoad in nlmixr2.
    viralLoad <- virus
    viralLoad ~ prop(propSd)
  })
}
