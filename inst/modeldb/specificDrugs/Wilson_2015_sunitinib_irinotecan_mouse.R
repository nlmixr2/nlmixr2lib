Wilson_2015_sunitinib_irinotecan_mouse <- function() {
  description <- paste(
    "Preclinical (mouse with HT-29 colorectal-cancer xenograft).",
    "Mechanistic tumor-growth PD model for the antiangiogenic agent sunitinib",
    "(reduces vascular carrying capacity) combined with the cytotoxic agent",
    "irinotecan (three-stage transit-cell-death chain following Simeoni et al.",
    "2004) and an empirical interaction term (Wilson 2015 Equation 4) in which",
    "the irinotecan transit-death rate kC depends on the cumulative",
    "pre-irinotecan sunitinib exposure. Drug input is K-PD (no pharmacokinetic",
    "data; each oral sunitinib or 5-min IV irinotecan dose enters its",
    "drug-amount compartment with normalized magnitude 1)."
  )
  reference <- paste(
    "Wilson S, Tod M, Ouerdani A, Emde A, Yarden Y, Adda Berkane A,",
    "Kassour S, Wei MX, Freyer G, You B, Grenier E, Ribba B.",
    "Modeling and predicting optimal treatment scheduling between",
    "the antiangiogenic drug sunitinib and irinotecan in preclinical",
    "settings. CPT Pharmacometrics Syst Pharmacol. 2015;4(12):720-727.",
    "doi:10.1002/psp4.12045.",
    sep = " "
  )
  vignette <- "Wilson_2015_sunitinib_irinotecan_mouse"
  paper_specific_compartments <- c("intIrinotecan", "cumSunitinibFrozen", "sunitinib", "irinotecan")

  units <- list(
    time          = "day",
    dosing        = "unitless (K-PD; each oral sunitinib or IV irinotecan dose enters its drug-amount compartment with magnitude 1)",
    concentration = "mm (geometric mean of three orthogonal tumor diameters (l*w*h)^(1/3); not a drug concentration)"
  )

  covariateData <- list()

  population <- list(
    species        = "mouse (athymic nude male with subcutaneous HT-29 colorectal-adenocarcinoma xenograft)",
    n_subjects     = 105L,
    n_studies      = 2L,
    age_range      = "5-6 weeks at randomisation",
    weight_range   = "approximately 20 g each",
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state  = "subcutaneous HT-29 human colorectal-adenocarcinoma xenograft (3.0e6 cells in 200 uL inoculated into the flank)",
    dose_range     = "sunitinib 40 mg/kg oral gavage once daily for 12 consecutive days; irinotecan 90 mg/kg single 5-min IV infusion; in the model the magnitudes are normalized to 1 per dose (K-PD)",
    regions        = "preclinical (in-vivo xenograft) at CellVax Laboratory facility, Maisons Alfort, France",
    notes          = "1,371 longitudinal tumor-diameter observations across 105 mice in two model-building experiments (sunitinib monotherapy and combined sunitinib+irinotecan; see Supplemental Tables S1-S3). Wilson 2015 Table 1 (interaction-model parameters used here) are typical-value estimates from experiment #2 fit to the median tumor size per group via nonlinear least-squares; no IIV or residual-error structure is reported for this fit. Tumor diameter measured every 2-3 days by handheld caliper as the geometric mean of length, width, and height (l*w*h)^(1/3) in mm."
  )

  ini({
    # Estimated structural parameters (Wilson 2015 Table 1, interaction model).
    # Log-transformed in ini() so re-fits stay on the positive real line; the
    # paper reports the linear-scale value in Table 1.
    lk_growth <- log(1.34);    label("Tumor growth rate lambda (1/day)")                                       # Wilson 2015 Table 1: lambda = 1.34 (RSE 10%)
    lb_cap    <- log(0.0027);  label("Carrying-capacity rate constant b (1/(mm*day))")                         # Wilson 2015 Table 1: b = 0.0027 mm^-1 day^-1 (RSE 0.04%)
    lbeta_s   <- log(0.0317);  label("Sunitinib effect coefficient on carrying-capacity loss (conc.unit^-1)")  # Wilson 2015 Table 1: beta_S = 0.0317 conc.unit^-1 (RSE 0.31%)
    lbeta_c   <- log(0.3847);  label("Irinotecan effect coefficient on cycling-cell killing (conc.unit^-1)")   # Wilson 2015 Table 1: beta_C = 0.3847 conc.unit^-1 (RSE 5%)
    lks       <- log(0.155);   label("Sunitinib-irinotecan interaction coefficient k_S (1/day)")               # Wilson 2015 Table 1: k_S = 0.155 (RSE 6%); replaces the noninteraction kC

    # Fixed structural parameters (Wilson 2015 Table 1 marks D0/K0/p_S/p_C as
    # "fixed"; alpha = 0.1 is fixed per page 723 prose, "the growth law is
    # 'nearly' Gompertzian"). Linear-scale fixed() because they are not
    # re-estimated in any companion paper that this file inherits from.
    d0             <- fixed(0.29);   label("Initial mean tumor diameter D(0) (mm)")                          # Wilson 2015 Table 1: D(t=0) = 0.29 mm (fixed)
    k0             <- fixed(7.43);   label("Initial carrying capacity K(0) (mm)")                            # Wilson 2015 Table 1: K(t=0) = 7.43 mm (fixed)
    ps_elim        <- fixed(2.12);   label("Sunitinib K-PD elimination rate p_S (1/day)")                    # Wilson 2015 Table 1: p_S = 2.12 day^-1 (fixed)
    pc_elim        <- fixed(0.0850); label("Irinotecan K-PD elimination rate p_C (1/day)")                   # Wilson 2015 Table 1: p_C = 0.0850 day^-1 (fixed)
    alpha_logistic <- fixed(0.1);    label("Generalized-logistic exponent alpha (unitless)")                 # Wilson 2015 page 723 prose: alpha = 0.1 (fixed, "nearly Gompertzian")

    # Residual error -- placeholder. Wilson 2015 fit the interaction model
    # (Table 1) to the per-group median tumor diameter by nonlinear least
    # squares (Methods, "Parameter estimation"), and reports no residual-error
    # structure for that fit. Placeholder values chosen on the same order as
    # Simeoni 2004 (oncology_xenograft_simeoni_2004.R) so the model is usable
    # for stochastic simulation; a downstream user re-fitting the model should
    # re-estimate these.
    propSd_tumor_size <- 0.10;  label("Proportional residual error on tumor diameter (fraction; placeholder)")  # Not reported in Wilson 2015 Table 1; placeholder
    addSd_tumor_size  <- 0.50;  label("Additive residual error on tumor diameter (mm; placeholder)")            # Not reported in Wilson 2015 Table 1; placeholder
  })

  model({
    # --- Individual (typical-value) structural parameters -------------------
    # No IIV: Wilson 2015 Table 1 reports point estimates from nonlinear
    # least squares on median tumor data with no random effects.
    k_growth <- exp(lk_growth)
    b_cap    <- exp(lb_cap)
    beta_s   <- exp(lbeta_s)
    beta_c   <- exp(lbeta_c)
    k_s      <- exp(lks)

    # --- Interaction term (Wilson 2015 Equation 4) --------------------------
    # kC = k_S * exp(integral_0^T_C S(t) dt), where T_C is the time of the
    # irinotecan dose. The integral is evaluated *at* T_C, then kC is held
    # constant thereafter (Wilson 2015 page 724: "kC is proportional to the
    # cumulated exposure to sunitinib prior to the dose of irinotecan").
    #
    # Implementation: intIrinotecan accumulates the irinotecan K-PD amount
    # over time. Before the first irinotecan dose intIrinotecan == 0, so
    # gate_cumS == 1 and cumSunitinibFrozen integrates sunitinib normally.
    # The first irinotecan dose adds 1 to the irinotecan compartment, which
    # then begins contributing to intIrinotecan; intIrinotecan crosses 0
    # immediately, gate_cumS flips to 0, and cumSunitinibFrozen freezes at
    # its T_C value for the remainder of the simulation. The gate is
    # monotonic (intIrinotecan never decreases), so subsequent irinotecan
    # doses cannot reopen it.
    gate_cumS <- 1.0
    if (intIrinotecan > 0.0) {
      gate_cumS <- 0.0
    }
    d/dt(intIrinotecan)      <- irinotecan
    d/dt(cumSunitinibFrozen) <- sunitinib * gate_cumS

    # Effective transit-death rate constant for the irinotecan damage chain.
    # When no sunitinib has been given before T_C (irinotecan monotherapy),
    # cumSunitinibFrozen = 0 and k_c_eff = k_s; when sunitinib is given
    # before T_C the cumulative exposure scales k_c_eff exponentially.
    k_c_eff <- k_s * exp(cumSunitinibFrozen)

    # --- ODE system (Wilson 2015 Equation 3, with interaction kC from Eq. 4) ---
    # State variables and their paper names:
    #   cycling_cells     = D1 (proliferating tumor diameter compartment, mm)
    #   damaged_cells1    = D2 (first transit-death compartment)
    #   damaged_cells2    = D3
    #   damaged_cells3    = D4
    #   carrying_capacity = K  (vasculature-limited maximum tumor diameter, mm)
    #   sunitinib        = S  (K-PD drug-amount compartment, unitless)
    #   irinotecan       = C  (K-PD drug-amount compartment, unitless)
    #
    # Generalized-logistic tumor growth with varying carrying capacity follows
    # Hahnfeldt et al. 1999 (Wilson 2015 reference 14); with alpha = 0.1 the
    # growth law is "nearly Gompertzian" per Wilson 2015 page 723. Irinotecan
    # injures cycling cells at rate beta_C * p_C * C; injured cells pass
    # through the three damaged-cell compartments at rate k_c_eff before
    # leaving the system (Wilson 2015 Figure 3).
    d/dt(cycling_cells)     <- k_growth * cycling_cells *
                                (1 - (cycling_cells / carrying_capacity)^alpha_logistic) -
                              beta_c * pc_elim * irinotecan * cycling_cells
    d/dt(damaged_cells1)    <- beta_c * pc_elim * irinotecan * cycling_cells -
                              k_c_eff * damaged_cells1
    d/dt(damaged_cells2)    <- k_c_eff * damaged_cells1 - k_c_eff * damaged_cells2
    d/dt(damaged_cells3)    <- k_c_eff * damaged_cells2 - k_c_eff * damaged_cells3
    d/dt(carrying_capacity) <- b_cap * cycling_cells * cycling_cells -
                              beta_s * ps_elim * sunitinib * carrying_capacity
    d/dt(sunitinib)        <- -ps_elim * sunitinib
    d/dt(irinotecan)       <- -pc_elim * irinotecan

    # --- Initial conditions -------------------------------------------------
    cycling_cells(0)        <- d0
    carrying_capacity(0)    <- k0
    damaged_cells1(0)       <- 0.0
    damaged_cells2(0)       <- 0.0
    damaged_cells3(0)       <- 0.0
    sunitinib(0)           <- 0.0
    irinotecan(0)          <- 0.0
    intIrinotecan(0)       <- 0.0
    cumSunitinibFrozen(0)  <- 0.0

    # --- Observation and residual error -------------------------------------
    # Wilson 2015 fits the sum D = D1 + D2 + D3 + D4 to the observed mean
    # tumor diameter (Wilson 2015 Equation 3, final line). Combined
    # additive + proportional residual error is a placeholder (paper does
    # not specify a residual-error structure for the Table 1 fit; see
    # ini() comments).
    tumor_size <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3
    tumor_size ~ add(addSd_tumor_size) + prop(propSd_tumor_size)
  })
}
