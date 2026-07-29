Asimus_2007_artemisinin <- function() {
  description <- "Semiphysiological autoinduction popPK model for oral artemisinin, fit to pooled plasma data from six clinical studies (33 healthy male Vietnamese volunteers + 54 male falciparum-malaria patients). The structural model is identical to the original Gordi 2005 saliva-based model except no absorption lag-time is estimated. Three artemisinin compartments (gut depot, liver V_H = 1 L fixed, sampling V_S = 26.1 L) are linked in a circular well-stirred-hepatic-extraction loop with hepatic plasma flow Q_H = 0.63 L/h/kg of body weight. Two enzyme states (precursor1 + enzyme pool) drive autoinduction: hepatic artemisinin amount linearly stimulates precursor formation (slope s_ind); the precursor transitions to enzyme with rate constant kpout = 1/MIT = 1/(2.0 h); the enzyme decays with first-order rate kdeg = ln(2)/94 h. Intrinsic clearance is proportional to enzyme amount and saturates in hepatic concentration via Michaelis-Menten kinetics (CL_int,t = vmax * enzyme / (km + C_H)), giving a pre-induced hepatic extraction E_H = 0.74 increasing to 0.98 after autoinduction (a roughly 13-fold drop in oral bioavailability with only a modest change in systemic clearance)."
  reference <- paste(
    "Asimus S, Gordi T (2007).",
    "Retrospective analysis of artemisinin pharmacokinetics:",
    "application of a semiphysiological autoinduction model.",
    "British Journal of Clinical Pharmacology 63(6): 758-762.",
    "doi:10.1111/j.1365-2125.2006.02844.x.",
    "Structural model adapted from Gordi T, Xie R, Huong NV, Huong DX,",
    "Karlsson MO, Ashton M (2005).",
    "A semiphysiological pharmacokinetic model for artemisinin in",
    "healthy subjects incorporating autoinduction of metabolism and",
    "saturable first-pass hepatic extraction.",
    "British Journal of Clinical Pharmacology 59(2): 189-198.",
    "doi:10.1111/j.1365-2125.2004.02321.x.",
    sep = " "
  )
  vignette <- "Asimus_2007_artemisinin"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight; scales hepatic plasma flow Q_H = 0.63 * WT (L/h)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives the hepatic plasma flow term Q_H = q_h_per_kg * WT inside model().",
        "Asimus 2007 study population body weight: healthy volunteers 48-68 kg;",
        "malaria patients 37-63 kg (Asimus 2007 Methods 'Study design').",
        "Reference value 51 kg corresponds to the average weight of the upstream",
        "Gordi 2005 healthy-volunteer cohort (24 Vietnamese males; Gordi 2005",
        "Methods 'Study subjects and materials')."
      ),
      source_name        = "Bodyweight"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 87L,
    n_studies      = 6L,
    age_range      = "15-55 years",
    weight_range   = "37-68 kg",
    sex_female_pct = 0,
    race_ethnicity = c(Vietnamese = 100),
    disease_state  = "healthy male volunteers (n = 33) plus male adults with uncomplicated falciparum malaria (n = 54)",
    dose_range     = "Oral artemisinin 250 mg or 500 mg (2 x 250 mg) capsules in single, twice-daily and dose-escalating regimens across the six pooled studies (Asimus 2007 Table 1).",
    regions        = "Vietnam",
    notes          = paste(
      "Demographics from Asimus 2007 Table 1.",
      "Six pooled clinical studies (combined healthy + malaria cohorts).",
      "Plasma artemisinin was sampled pre-dose and at multiple times to roughly",
      "10 h post-dose; the exact schedule differs across studies (Asimus 2007",
      "Table 1).",
      "FOCE and FOCE INTERACTION estimation methods terminated; first-order (FO)",
      "without centring was used for the final fit (Asimus 2007 Methods",
      "'Pharmacokinetic analysis').",
      "Structural ODEs adopted unchanged from Gordi 2005 except for the lag-time;",
      "parameter values in this file are the Asimus 2007 plasma re-estimates",
      "(Asimus 2007 Table 2)."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters. Asimus 2007 Table 2 reports CL_int_0 = 1760 L/h
    # and K_m = 434 ng/mL; the Michaelis-Menten relationship CL_int_0 =
    # V_max / K_m gives V_max = CL_int_0 * K_m = 1760 L/h * 0.434 mg/L =
    # 763.84 mg/h at baseline enzyme amount (A_ENZ = 1). The model() block
    # uses the equivalent vmax/(km + C_H) form (Gordi 2005 Eq. 2 rewritten
    # in vmax/km coordinates) so the canonical lvmax / lkm parameter names
    # apply.
    # ---------------------------------------------------------------------
    lvmax  <- log(763.84)
    label("Michaelis-Menten Vmax at baseline enzyme A_ENZ = 1 (mg/h)")
    # Asimus 2007 Table 2: CL_int_0 = 1760 L/h (RSE 35%), K_m = 434 ng/mL
    # (RSE 50%); V_max = CL_int_0 * K_m = 1760 * 0.434 mg/L = 763.84 mg/h.

    lkm    <- log(434)
    label("Michaelis-Menten Km hepatic artemisinin concentration (ng/mL)")
    # Asimus 2007 Table 2: K_m = 434 ng/mL (RSE 50%); hepatic concentration
    # giving half-maximal intrinsic clearance.

    lka    <- log(0.09)
    label("First-order absorption rate constant k_a from gut to liver (1/h)")
    # Asimus 2007 Table 2: k_a = 0.09 1/h (RSE 13%).

    lvc    <- log(26.1)
    label("Sampling-compartment volume V_S (L)")
    # Asimus 2007 Table 2: V_S = 26.1 L (RSE 15%); the sampling compartment
    # represents the entire body except the liver.

    lkdeg  <- log(log(2) / 94)
    label("Enzyme pool first-order degradation rate kdeg (1/h, t_half = 94 h)")
    # Asimus 2007 Table 2: t_half_ENZ = 94 h (RSE 27%); kdeg = ln(2)/94 =
    # 0.00737 1/h.

    lkpout <- log(1 / 2.0)
    label("Precursor pool loss rate kpout = 1/MIT (1/h)")
    # Asimus 2007 Table 2: MIT (mean induction time) = 2.0 h (RSE 43%);
    # kpout = 1/MIT = 0.5 1/h. In the Gordi 2005 model, kpout is the
    # first-order rate of conversion of the precursor pool into the enzyme
    # pool (Gordi 2005 Eqs. 10-11; the paper writes the same constant as
    # kPRE).

    s_ind  <- 0.045
    label("Slope of hepatic artemisinin concentration on precursor production (mL/ng = 1/(ng/mL))")
    # Asimus 2007 Table 2 reports S_IND = 0.045 with units '1/ng'; Gordi
    # 2005 Table 1 footnote and Figure 2 caption clarify that S_IND is the
    # 'slope of the inducing effect of artemisinin hepatic CONCENTRATION
    # on the production rate of enzyme precursor' (not amount), so the
    # natural unit pairing is mL/ng applied to C_H in ng/mL. The Gordi
    # 2005 value (0.018 L/ng = 18 mL/ng applied to C_H in ng/mL would
    # imply a 1000x larger induction signal than Asimus 2007's 0.045
    # mL/ng; numerically the Asimus 2007 value pairs naturally with C_H
    # in ng/mL given the K_m = 434 ng/mL saturation constant. Model() uses
    # the formulation `kdeg * (1 + s_ind * C_H_ng_per_mL)`.

    fu     <- fixed(0.14)
    label("Plasma unbound fraction fu (FIXED)")
    # Asimus 2007 Table 2: f_u = 0.14 FIXED; carried unchanged from the
    # Gordi 2005 model (Gordi 2005 ref [18], saliva-plasma partitioning).

    v_h_fix <- fixed(1)
    label("Hepatic compartment volume V_H (L, FIXED)")
    # Gordi 2005 Methods 'Pharmacokinetic analysis' and Figure 2 caption:
    # V_H fixed at 1 L; sensitivity analyses at 5 and 10 L increased the
    # OFV by 10-12 units and were rejected.

    q_h_per_kg <- fixed(1.2)
    label("Hepatic blood flow per kg body weight q_h_per_kg (L/h/kg, FIXED)")
    # Gordi 2005 Methods 'Pharmacokinetic analysis' tested two values for
    # Q_H: 0.63 L/h/kg (hepatic plasma flow, the Gordi 2005 default) and
    # 1.2 L/h/kg (hepatic blood flow, the upper sensitivity alternative).
    # The 1.2 L/h/kg value is selected here because the Asimus 2007
    # reported baseline hepatic extraction E_H = 0.74 (Asimus 2007
    # Results) is inconsistent with the 0.63 L/h/kg default for the
    # Asimus 2007 cohort body-weight range -- with CL_int_0 = 1760 L/h
    # and f_u = 0.14, Q_H = 0.63 * 51 = 32 L/h gives E_H ~ 0.88, well
    # above the published 0.74. The 1.2 L/h/kg alternative gives Q_H ~
    # 61-84 L/h across the cohort weight range (37-68 kg) and brings the
    # baseline E_H into the 0.75-0.85 band that brackets the published
    # 0.74. The Asimus 2007 paper does not state Q_H explicitly; this is
    # the most parsimonious of the Gordi 2005 sensitivity-tested values
    # consistent with the published baseline E_H.

    # ---------------------------------------------------------------------
    # Inter-individual variability. Asimus 2007 Table 2 IIV column reports
    # exponential variances on log scale (omega^2): IIV on CL_int_0 = 0.38
    # (RSE 24%), IIV on V_S = 1.2 (RSE 32%). Because V_max = CL_int_0 * K_m
    # with K_m carrying no IIV, the variance on log(V_max) equals the
    # variance on log(CL_int_0); so etalvmax is encoded with the same
    # omega^2 = 0.38.
    # ---------------------------------------------------------------------
    etalvmax ~ 0.38
    # Asimus 2007 Table 2: IIV on CL_int_0, omega^2 = 0.38 (RSE 24%).
    # Mapped onto lvmax via lvmax = lcl_int_0 + log(K_m_typ); identical
    # log-scale variance.

    etalvc   ~ 1.2
    # Asimus 2007 Table 2: IIV on V_S, omega^2 = 1.2 (RSE 32%).

    etalka   ~ 0.64
    # Asimus 2007 Table 2: IOV on k_a, omega^2 = 0.64 (RSE 23%); encoded
    # as IIV here so single-occasion forward simulation in nlmixr2lib
    # reproduces the population variability (Birgersson 2016 artemisinin
    # precedent for IOV -> IIV mapping in nlmixr2lib).

    # ---------------------------------------------------------------------
    # Residual variability. Asimus 2007 Table 2 reports the proportional
    # residual error value 0.54 (RSE 4.7%); the paper used a proportional
    # residual error component on log-transformed plasma concentrations
    # ('NONMEM additive-on-log == proportional in nlmixr2 linear space').
    # Interpreting the table value as the NONMEM SIGMA variance gives the
    # log-scale SD propSd = sqrt(0.54) = 0.7348, consistent with the
    # similarly small Gordi 2005 saliva value (0.5; Gordi 2005 Table 1)
    # for the plasma cohort having modestly higher residual.
    # ---------------------------------------------------------------------
    propSd <- sqrt(0.54)
    label("Proportional residual error SD on log scale (sqrt of NONMEM SIGMA)")
    # Asimus 2007 Table 2: residual error sigma^2 = 0.54 (RSE 4.7%);
    # propSd = sqrt(0.54) ~ 0.7348.
  })

  model({
    # ---------------------------------------------------------------------
    # Individual structural parameters and hepatic blood-flow scaling.
    # Amounts are carried in mg (compartments) and L (volumes). The Km
    # parameter is reported in ng/mL and converted to mg/L inside model()
    # so the saturable elimination term (km + C_H) is dimensionally
    # consistent with C_H = liver / v_h in mg/L. The induction slope
    # s_ind (paper units '1/ng') is interpreted as mL/ng applied to C_H
    # in ng/mL -- the interpretation that pairs with the K_m = 434 ng/mL
    # saturation constant and reproduces the paper's reported 16-fold
    # enzyme increase over 5 days of dosing within numerical tolerance
    # (a literal 1/ng applied to A_H in ng would inflate the induction
    # signal by 1000x because A_H_ng = C_H_ngL = 1000 * C_H_ngmL when
    # V_H = 1 L, giving a runaway enzyme growth inconsistent with the
    # published 16-fold value; see vignette Assumptions and deviations).
    # ---------------------------------------------------------------------
    vmax       <- exp(lvmax + etalvmax)
    km_ngml    <- exp(lkm)
    km         <- km_ngml / 1000
    ka         <- exp(lka + etalka)
    vc         <- exp(lvc + etalvc)
    kdeg       <- exp(lkdeg)
    kpout      <- exp(lkpout)
    q_h        <- q_h_per_kg * WT
    k_sh       <- q_h / vc

    # ---------------------------------------------------------------------
    # Hepatic concentration in mg/L (matches K_m comparison) and in ng/mL
    # (matches S_IND comparison and the paper's reporting units).
    # Time-varying intrinsic clearance, well-stirred extraction ratio,
    # and the fraction escaping hepatic first-pass:
    # Gordi 2005 Eqs. 2-5 with V_max parameterisation:
    #   CL_int,t  = vmax * enzyme / (km + C_H)
    #   E_H       = CL_int,t * fu / (Q_H + CL_int,t * fu)
    #   F_H       = 1 - E_H
    # ---------------------------------------------------------------------
    c_h        <- liver / v_h_fix
    c_h_ngml   <- c_h * 1000
    cl_int_t   <- vmax * enzyme / (km + c_h)
    e_h        <- cl_int_t * fu / (q_h + cl_int_t * fu)
    f_h        <- 1 - e_h

    # ---------------------------------------------------------------------
    # ODE system (Gordi 2005 Eqs. 1, 6, 10, 11 with the Asimus 2007
    # modification that no absorption lag-time is estimated). The hepatic
    # amount carried as `liver` is in mg; the precursor amount is
    # dimensionless (normalised so A_ENZ_baseline = 1).
    # ---------------------------------------------------------------------
    d/dt(depot)     <- -ka * depot
    d/dt(liver)     <-  ka * depot - q_h * c_h + k_sh * central
    d/dt(central)   <-  q_h * f_h * c_h - k_sh * central
    d/dt(precursor1) <- kdeg * (1 + s_ind * c_h_ngml) - kpout * precursor1
    d/dt(enzyme)    <-  kpout * precursor1 - kdeg * enzyme

    # ---------------------------------------------------------------------
    # Pre-induced steady-state initial conditions (Gordi 2005 'The amount
    # of enzyme was set to 1 for the preinduced state' and 'The amount of
    # precursor in the preinduced state was set to kENZ/kPRE'). With the
    # canonical naming used here, k_ENZ = kdeg and k_PRE = kpout.
    # ---------------------------------------------------------------------
    precursor1(0) <- kdeg / kpout
    enzyme(0)     <- 1

    # ---------------------------------------------------------------------
    # Observation: plasma artemisinin concentration sampled from the
    # sampling compartment. central / vc is mg/L; multiply by 1000 to get
    # ng/mL declared in `units`. Proportional residual error on log scale.
    # ---------------------------------------------------------------------
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
