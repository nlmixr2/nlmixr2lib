Dodds_2005_rfxiii_cyno <- function() {
  description <- paste(
    "Preclinical (cynomolgus monkey). Three-state mechanistic population PK model",
    "for recombinant Factor XIII A2 dimer (rA2) administered IV bolus, with",
    "endogenous constant-influx production of A2 dimer and B monomer, mass-action",
    "association A2 + 2 B -> A2B2 heterotetramer following the 1A2 + 2B -> 1A2B2",
    "stoichiometry, first-order elimination of each species, and three ELISA",
    "assay outputs (total A2 = A2 + A2B2; A2B2 tetramer; free B) with",
    "proportional + proportional + additive residual error. Parameters are",
    "weight-normalised throughout (mass / kg). Estimated in NONMEM V (Dodds 2005).",
    sep = " "
  )
  reference <- paste(
    "Dodds MG, Visich JE, Vicini P.",
    "Population Pharmacokinetics of Recombinant Factor XIII in Cynomolgus Monkeys.",
    "AAPS J. 2005 Oct 27;7(3):Article 70 (E693-E702).",
    "doi:10.1208/aapsj070370.",
    sep = " "
  )
  vignette <- "Dodds_2005_rfxiii_cyno"
  paper_specific_compartments <- c("A2", "A2B2", "B")

  units <- list(
    time          = "hour",
    dosing        = "mg/kg (IV bolus of rA2 dimer; doses are weight-normalised throughout)",
    concentration = "mg/L (each assay output: ELISA reading in mass/volume on plasma)"
  )

  covariateData <- list(
    # No covariates entered the final model. Body weight is baked into the
    # mg/kg parameterisation (doses, masses, influx rates, and apparent
    # volumes are all per-kg); the elimination rate constants carry the
    # only between-subject variability (Dodds 2005, Population/Statistical
    # Model section).
  )

  population <- list(
    species        = "cynomolgus monkey (Macaca fascicularis, Chinese origin)",
    n_subjects     = 30L,
    n_studies      = 2L,
    weight_typical = "~3.5 kg (Dodds 2005 Results: 'normal cynomolgus monkeys ~3.5 kg = 270 mL blood volume')",
    dose_range     = "Single-dose study: 0.5, 1.0, 5.0 mg/kg IV bolus rA2 (n = 4 per dose). Multiple-dose study: 0.3, 3.0, 6.0 mg/kg IV bolus rA2 once daily for 14 days (n = 6 per dose).",
    n_observations = 1674L,
    sampling       = "Single-dose: 0.25, 1, 2, 4, 8, 24, 72, 120, 168, 240, 336, 504, 672 h post-dose. Multiple-dose: pre-dose and 15 min post-dose on Days 2, 7, 10; 0, 0.25, 6, 24 h on Day 1; pre-dose and 0.25, 6, 24, 48 h post-dose on Day 14.",
    assays         = "Three validated ELISAs measured from the same plasma space: total A2 (anti-A2 capture and detect; sees A2 + A2B2; sensitivity 0.79 mg/L), A2B2 tetramer (anti-B capture, anti-A detect; sensitivity 1.5 mg/L), and free B (anti-mouse capture of monoclonal anti-free-B + biotinylated anti-B detect; sensitivity 0.313 mg/L).",
    disease_state  = "Healthy preclinical animals (no FXIII deficiency).",
    regions        = "USDA-regulated laboratory (location not stated in the paper).",
    notes          = "Subject counts derived from Dodds 2005 Study Design section. Weight range and per-subject ages / sex are not reported in the paper; the 3.5 kg figure is the paper's own typical-monkey reference used to compare model volumes against blood volume."
  )

  ini({
    # Structural parameters -- Table 1 of Dodds 2005 ("PK parameter (theta)" block).
    # Each value is the published point estimate; %RSEs from the same table.
    linf_dimer    <- log(0.000622); label("Endogenous influx of A2 dimer (mg/kg/hr)")    # Table 1: InfA2 = 0.000622 mg/kg/hr (%RSE 17.0%)
    linf_monomer  <- log(0.0121);   label("Endogenous influx of free B monomer (mg/kg/hr)")  # Table 1: InfB = 0.0121 mg/kg/hr (%RSE 8.76%)
    lvt_dimer     <- log(0.0407);   label("Apparent volume for total A2 assay (L/kg)")   # Table 1: VtA2 = 0.0407 L/kg (%RSE 3.93%)
    lv_tetramer   <- log(0.00934);  label("Apparent volume for A2B2 tetramer assay (L/kg)")  # Table 1: VA2B2 = 0.00934 L/kg (%RSE 3.62%)
    lv_monomer    <- log(0.0598);   label("Apparent volume for free B assay (L/kg)")     # Table 1: VB = 0.0598 L/kg (%RSE 12.3%)
    lke_dimer     <- log(0.208);    label("First-order elimination of A2 dimer (1/hr)")   # Table 1: keA2 = 0.208 1/hr (%RSE 9.09%); half-life 3.33 h
    lke_tetramer  <- log(0.0102);   label("First-order elimination of A2B2 tetramer (1/hr)") # Table 1: keA2B2 = 0.0102 1/hr (%RSE 10.5%); half-life 2.83 d
    lke_monomer   <- log(0.176);    label("First-order elimination of free B monomer (1/hr)") # Table 1: keB = 0.176 1/hr (%RSE 15.9%); half-life 3.94 h
    lk_assoc      <- log(6.59);     label("Bimolecular association rate constant (kg / mg / hr)") # Table 1: Ka = 6.59 mg^-1 kg /hr (%RSE 19.7%)

    # Between-subject variability -- Table 1 of Dodds 2005 ("Population
    # variability (Omega)" block). Paper reports variances on the
    # exponential-IIV scale (NONMEM Omega), interpretable as sqrt(variance) =
    # approximate CV; the reported BSV columns 25.1% / 16.4% / 38.9% are
    # sqrt(0.0629) / sqrt(0.0268) / sqrt(0.151). Only the three elimination
    # rate constants carry IIV; the other parameters were treated as
    # population-level fixed effects (Dodds 2005, Population/Statistical
    # Model section).
    etalke_dimer    ~ 0.0629   # Table 1: Var[keA2]    = 0.0629 (BSV 25.1%)
    etalke_tetramer ~ 0.0268   # Table 1: Var[keA2B2]  = 0.0268 (BSV 16.4%)
    etalke_monomer  ~ 0.151    # Table 1: Var[keB]     = 0.151  (BSV 38.9%)

    # Residual error -- Table 1 of Dodds 2005 ("Assay variability (Sigma)"
    # block). The paper reports variances; nlmixr2 takes residual SDs, so
    # values are sqrt of the reported variances. Total A2 and A2B2 assays
    # are proportional (unitless variance in Table 1); free B is additive
    # with variance in mg^2/L^2.
    propSd_total_a2_assay <- sqrt(0.0730); label("Proportional residual SD on total A2 assay (fraction)")  # Table 1: Var = 0.0730; RUV 27.0%
    propSd_a2b2_assay     <- sqrt(0.0751); label("Proportional residual SD on A2B2 assay (fraction)")      # Table 1: Var = 0.0751; RUV 27.4%
    addSd_b_assay         <- sqrt(0.0857); label("Additive residual SD on free B assay (mg/L)")            # Table 1: Var = 0.0857 mg^2/L^2; RUV 0.293 mg/L
  })

  model({
    # Individual (per-subject) elimination rates with log-normal IIV.
    ke_dimer    <- exp(lke_dimer    + etalke_dimer)
    ke_tetramer <- exp(lke_tetramer + etalke_tetramer)
    ke_monomer  <- exp(lke_monomer  + etalke_monomer)

    # Population-level structural parameters (no IIV in the source paper).
    inf_dimer   <- exp(linf_dimer)
    inf_monomer <- exp(linf_monomer)
    vt_dimer    <- exp(lvt_dimer)
    v_tetramer  <- exp(lv_tetramer)
    v_monomer   <- exp(lv_monomer)
    k_assoc     <- exp(lk_assoc)

    # Steady-state initial conditions from Dodds 2005 Equations 4-7
    # (analytical solution at t=0 of the system below with the
    # bolus dose D set to zero). psi has units of 1/hr^2; A2(0), A2B2(0),
    # and B(0) all in mg/kg.
    psi <- sqrt(4*inf_monomer*k_assoc*ke_dimer*ke_monomer +
                (k_assoc*(2*inf_dimer - inf_monomer) + ke_dimer*ke_monomer)^2)
    A2(0)   <- (2*inf_dimer*ke_monomer) /
               (k_assoc*(inf_monomer - 2*inf_dimer) + ke_dimer*ke_monomer + psi)
    A2B2(0) <- (3 / (4*k_assoc*ke_tetramer)) *
               (k_assoc*(inf_monomer + 2*inf_dimer) + ke_dimer*ke_monomer - psi)
    B(0)    <- (1 / (2*k_assoc*ke_monomer)) *
               (k_assoc*(inf_monomer - 2*inf_dimer) - ke_dimer*ke_monomer + psi)

    # ODE system from Dodds 2005 Equations 1-3 (the D*delta(t) bolus
    # term is supplied by the event table as an evid=1 amt dosed into
    # the A2 compartment, so it is not written explicitly here).
    # Coefficients 2 (on B loss) and 3 (on A2B2 gain) follow from the
    # 1A2 + 2B -> 1A2B2 stoichiometry; the "3" is the sum of the
    # A2 (1 unit mass) and B (2 unit masses) contributions per
    # reaction.
    d/dt(A2)   <- inf_dimer   - ke_dimer*A2     - k_assoc*A2*B
    d/dt(A2B2) <-              -ke_tetramer*A2B2 + 3*k_assoc*A2*B
    d/dt(B)    <- inf_monomer - ke_monomer*B   - 2*k_assoc*A2*B

    # Measurement equations -- Dodds 2005 Equations 8-10. Each assay
    # divides the relevant amount (mg/kg) by an apparent volume (L/kg),
    # giving concentration (mg/L). The total A2 assay sees both A2 and
    # A2B2; A2B2 and free B assays are species-specific.
    total_a2_assay <- (A2 + A2B2) / vt_dimer
    a2b2_assay     <- A2B2 / v_tetramer
    b_assay        <- B / v_monomer

    total_a2_assay ~ prop(propSd_total_a2_assay)
    a2b2_assay     ~ prop(propSd_a2b2_assay)
    b_assay        ~ add(addSd_b_assay)
  })
}
