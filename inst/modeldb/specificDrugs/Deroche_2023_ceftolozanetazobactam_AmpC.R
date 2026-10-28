Deroche_2023_ceftolozanetazobactam_AmpC <- function() {
  description <- "In vitro (Pseudomonas aeruginosa PAO1 carrying the single AmpC G183D cephalosporinase mutation). Semi-mechanistic sequential time-kill PK/PD model with adaptive resistance (Mouton-type) describing C/T killing and regrowth. Logistic bacterial growth gated by a steep lag-time switch is opposed by a drug effect that sums two Emax terms, one for the non-adapted fraction (ARoff, potency EC50,off) and one for the adapted fraction (ARon, potency EC50,on), with ARon = 1 - exp(-kon*t). One of 11 strain-by-drug files from this paper (this file: PAO1-AmpCG183D versus C/T). Emax and kon are shared across the PAO1-derived C/T strain group; all other parameters are strain-specific."
  reference <- "Deroche L, Plouzeau C, Bemer P, Dalichampt M, Michaud A, Chavanet P, Burucoa C, Rammaert B, Gregoire N, Marchand S, Cremniter J. Characterization of Pseudomonas aeruginosa resistance to ceftolozane-tazobactam due to ampC and/or ampD mutations observed during treatment using semi-mechanistic PKPD modeling. Antimicrob Agents Chemother. 2023 Oct 18;67(10):e00480-23. doi:10.1128/aac.00480-23. PMID:37695298; PMCID:PMC10583683. Structural model: effect Eq 1 (main text), adaptation Eqs 2-3 (main text), full differential equations Suppl Text S1 Eqs 2-6, schematic Fig 5. Parameter estimates: Table S2."
  vignette <- "Deroche_2023_ceftolozane_imipenem_resistance"
  units <- list(time = "h", dosing = "mg/L", concentration = "mg/L")

  # Drug concentration is a user-dosed state held constant within each
  # static-concentration time-kill; `cdrug` is not a canonical compartment name.
  paper_specific_compartments <- c("cdrug")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # matrix. Verified against the paper (Fig 5 schematic, Suppl Text S1).
  compartmentData <- list(
    bacteria = list(
      analyte = "Pseudomonas aeruginosa PAO1-AmpCG183D viable count",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    aroff = list(
      analyte = "non-adapted (drug-susceptible) bacterial fraction",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    aron = list(
      analyte = "adapted (adaptive-resistance) bacterial fraction",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    cdrug = list(
      analyte = "ceftolozane/tazobactam (C/T)",
      units = "mg/L",
      specimen = "administration site",
      verified = TRUE
    )
  )

  # OCC selects the time-kill-specific growth lag (1 = first TKC, 2 = second
  # TKC); OCC is the canonical occasion column.
  covariateData <- list(
    OCC = list(
      description = "Time-kill-curve index selecting the growth lag time: 1 = first TKC, 2 = second TKC.",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      source_name = "TKC index (Fig 1/3 sub-rows TKC1/TKC2)"
    )
  )

  population <- list(
    species = "in vitro (Pseudomonas aeruginosa PAO1-AmpCG183D)",
    n_subjects = 1L,
    n_studies = 1L,
    disease_state = "Catheter-related P. aeruginosa bacteremia; in-vitro sequential time-kill characterization of PAO1-derived C/T susceptibility",
    model_system = "Sequential 24-well static-concentration time-kill curves in cation-adjusted Mueller-Hinton broth II at 35 C; two successive 30-h TKCs, the second seeded at 24 h from regrowth wells of the first, sampled at 0, 3, 6, 24, 30 h with a 200 CFU/mL limit of quantification",
    notes = paste(
      "Semi-mechanistic PKPD model with adaptation (NONMEM 7.4, Laplacian, Beal M3 for BLQ) fit simultaneously to log10 CFU/mL counts across the PAO1-derived C/T strain group.",
      "Emax = 1.94 /h and kon = 0.0897 /h were estimated as one value per strain group and are reproduced here; all other parameters are specific to this strain.",
      "koff was fixed to 0: once adapted, bacteria do not revert within the experiment, so ARon = 1 - exp(-kon*t).",
      "The inoculum (1e6 CFU/mL) is a dose to `bacteria` (evid = 1 at t = 0; evid = 5 reinoculation at the start of the second TKC), so tad(bacteria) gives the per-TKC growth-lag clock.",
      "The second TKC used antibiotic up to 64x MIC; adaptive resistance develops continuously from t = 0 (ARon/ARoff are not reset between TKCs) while the growth lag resets per TKC (selected by OCC).",
      "Companion files (same paper): the other 10 Deroche_2023 strain-by-drug models."
    )
  )

  ini({
    # --- Bacterial growth (Suppl Text S1 Eqs 2, 6; Table S2) -----------------------
    knet        <- 0.823; label("Net bacterial growth rate constant k_net (1/h)")                 # Table S2, PAO1-AmpCG183D column: k_net = 0.823 /h
    log10bmax   <- 9.97; label("Maximal bacterial density B_max (log10 CFU/mL)")                 # Table S2, PAO1-AmpCG183D column: B_max = 9.97 log10 CFU/mL
    tlag1       <- fixed(0); label("Growth lag time during the first TKC (h)")           # Table S2, PAO1-AmpCG183D column: Tlag1 = 0 (fixed)
    tlag2       <- 3.09; label("Growth lag time during the second TKC (h)")             # Table S2, PAO1-AmpCG183D column: Tlag2 = 3.09 h

    # --- Drug effect with adaptive resistance (Eq 1; Suppl Eqs 3-5; Table S2) ------
    emax        <- 1.94; label("Maximum effect rate constant E_max (1/h); shared across the PAO1-derived C/T strain group")  # Table S2 (one value per PAO1 C/T strain group): E_max = 1.94 /h
    ec50off     <- 0.0606; label("C/T concn for 50% of E_max in non-adapted bacteria EC50,off (mg/L)")  # Table S2, PAO1-AmpCG183D column: EC50,off = 0.0606 mg/L
    ec50on      <- 19.4; label("C/T concn for 50% of E_max in adapted bacteria EC50,on (mg/L)")       # Table S2, PAO1-AmpCG183D column: EC50,on = 19.4 mg/L
    kon         <- 0.0897; label("Adaptive-resistance development rate constant k_on (1/h); shared across the PAO1-derived C/T strain group")  # Table S2 (one value per PAO1 C/T strain group): k_on = 0.0897 /h
    koff        <- fixed(0); label("Reversal rate constant of adaptive resistance k_off (1/h)")  # Table S2, PAO1-AmpCG183D column: k_off = 0 (fixed)

    # --- Residual error (additive on log10 CFU/mL; sigma2 is the variance) ---
    addSd       <- 1.0954; label("Additive residual SD on log10 CFU/mL scale")                      # Table S2, PAO1-AmpCG183D column: sigma^2 = 1.2 (variance); SD = sqrt(sigma^2) = 1.0954
  })

  model({
    # Growth lag switch (Suppl Eq 6: Delay = Time^20 / (Tlag^20 + Time^20)),
    # written in the algebraically identical reciprocal form so a fixed
    # Tlag = 0 gives Delay = 1 at t = 0 without a 0/0 singularity.
    # tad(bacteria) is the per-TKC clock (it restarts at the second-TKC
    # reinoculation event).
    tlag  <- tlag1 * (OCC == 1) + tlag2 * (OCC == 2)
    delay <- 1 / (1 + (tlag / (tad(bacteria) + 1e-10))^20)

    bmax <- 10^log10bmax

    # Effect = non-adapted Emax term + adapted Emax term (Eq 1 / Suppl Eq 5).
    effect <- emax * cdrug / (ec50off + cdrug) * aroff +
      emax * cdrug / (ec50on + cdrug) * aron

    d/dt(bacteria) <- knet * (1 - bacteria / bmax) * delay * bacteria -
      effect * bacteria
    d/dt(aroff)    <- koff * aron - kon * aroff
    d/dt(aron)     <- kon * aroff - koff * aron
    d/dt(cdrug)    <- 0

    # The inoculum is supplied as a dose to `bacteria` (evid = 1 at the start
    # of the first TKC; evid = 5 reinoculation at the start of the second),
    # not as an initial condition, so that tad(bacteria) provides the per-TKC
    # growth-lag clock. Methods: suspension adjusted to 1e6 CFU/mL at t = 0.
    bacteria(0) <- 0
    aroff(0)    <- 1
    aron(0)     <- 0

    log_cfu <- log10(bacteria)
    log_cfu ~ add(addSd)
  })
}
