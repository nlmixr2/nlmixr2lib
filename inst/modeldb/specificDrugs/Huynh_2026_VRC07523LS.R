Huynh_2026_VRC07523LS <- function() {
  description <- "Two-compartment population PK model with zero-order subcutaneous absorption, allometric weight scaling, and binary effects of age (adult vs infant) and repeat dosing for the broadly neutralizing HIV-1 monoclonal antibody VRC07-523LS in healthy adults and HIV-exposed infants (Huynh 2026)."
  reference <- "Huynh D, Nikanjam M, Cunningham CK, McFarland EJ, Muresan P, Perlowski C, Yin DE, Moye J, Spiegel H, Gama L, Gaudinski M, Capparelli EV. Model-based assessment of VRC07-523LS dosing in infants through population pharmacokinetic-pharmacodynamic modelling in adults and infants. J Antimicrob Chemother. 2026; doi:10.1093/jac/dkaf449"
  vignette <- "Huynh_2026_VRC07523LS"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Updated at each PK visit in the infant cohort to track growth (Huynh 2026 Methods, Population pharmacokinetic analysis: 'PK parameters were scaled allometrically based on observed potential covariates ... Weight was updated as a PK study visit for both the modelling and simulation'). Used for allometric scaling on CL and Q with exponent 0.85 and on Vc and Vp with exponent 1.0; reference 70 kg.",
      source_name        = "WTKG"
    ),
    CHILD = list(
      description        = "Age-group indicator: 1 = infant, 0 = adult",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult)",
      notes              = "Time-fixed per subject. Encodes the IMPAACT P1112 infant cohort (median post-natal age 2.5 days) as CHILD = 1 and the VRC 605 adult cohort (median 29.4 years) as CHILD = 0. The published structural CL (47.76 mL/d/70 kg) and D1 (36 h) are reported with infants as the reference category, so the model file applies the effect via (1 + e_child_cl * (1 - CHILD)) (and analogously on D1) with a positive coefficient to preserve the paper's verbatim values; adults carry CL/F = 158.4 x (WTKG/70)^0.85 x 1.69 and D1 = 36 x 2.79 hours per Huynh 2026 Table 2 footnote.",
      source_name        = "AGEGRP (1 = infant, 0 = adult)"
    ),
    CYCLE = list(
      description        = "Dose-number indicator: 1 = first dose, 2+ = repeat dose",
      units              = "(count)",
      type               = "count",
      reference_category = "1 (first dose)",
      notes              = "Time-varying. Increment at the start of each new dose; supply CYCLE >= 1 at every observation row. Source paper applies a multiplicative factor of 1.49 on Vss (= Vc + Vp) for the second and subsequent doses (Huynh 2026 Table 2; 'V_ss increased by 49% in comparison to the first dose'). Mechanism is hypothesized to involve TMDD or FcRn-binding changes but was not pursued in the paper.",
      source_name        = "DOSEN (1 = first, >=2 = repeat)"
    )
  )

  population <- list(
    n_subjects     = 46L,
    n_studies      = 2L,
    age_range      = "Infants (IMPAACT P1112): 37-42 weeks gestational age, 0-5 days post-natal age at first dose. Adults (VRC 605): 22-48 years.",
    age_median     = "Infants 2.5 days post-natal age; adults 29.4 years",
    weight_range   = "Infants: 2.2-4.3 kg. Adults: 45-97 kg.",
    weight_median  = "Infants 2.8 kg; adults 71.1 kg",
    sex_female_pct = 50.0,
    disease_state  = "Healthy adults (VRC 605) and HIV-exposed (uninfected) infants (IMPAACT P1112) enrolled for HIV-1 prophylaxis evaluation",
    dose_range     = "Infants: 80 mg s.c. single dose (n=12) or 80 mg s.c. at birth + 100 mg s.c. at week 12 (n=9). Adults: 1, 5, 20, or 40 mg/kg i.v. single (n=12); 5 mg/kg s.c. single (n=3); 20 mg/kg i.v. q12w x3 (n=5); 5 mg/kg s.c. q12w x3 (n=5).",
    regions        = "United States (multi-site IMPAACT P1112 and VRC 605 trials)",
    n_infants      = 21L,
    n_adults       = 25L,
    sampling       = "638 VRC07-523LS serum concentrations across 84 weeks (211 infant samples; 480 adult samples after exclusion of one infant for tolerability). Quantification by anti-idiotype antibody capture ELISA.",
    notes          = "Per Huynh 2026 Table 1 (Baseline demographics and study designs) and Methods, Patient population. Sex split combines 12 male / 9 female infants and 11 male / 14 female adults."
  )

  ini({
    # Structural PK parameters - Huynh 2026 Table 2 final estimates.
    # Reference subject: 70 kg, infant (CHILD = 1), first dose (CYCLE = 1).
    # Vc and Vp are absolute volumes (L) for a 70 kg subject; CL and Q are
    # absolute clearances (L/day) for a 70 kg subject. Source values are
    # converted to L/day where the paper reports mL/day or L/h.
    lvc     <- log(1.48);          label("Central volume of distribution Vc for a 70 kg subject (L)")                                # Huynh 2026 Table 2: theta1 (Vc; L/70) = 1.48
    lvp     <- log(2.28);          label("Peripheral volume of distribution Vp for a 70 kg subject (L)")                             # Huynh 2026 Table 2: theta2 (Vp; L/70) = 2.28
    lcl     <- log(47.76 / 1000);  label("Clearance CL for a 70 kg infant (L/day) - paper reports 47.76 mL/d/70 kg")                 # Huynh 2026 Table 2: theta3 (CL; mL/d/70 kg) = 47.76
    lq      <- log(0.0243 * 24);   label("Inter-compartmental clearance Q for a 70 kg subject (L/day) - paper reports 0.0243 L/h/70 kg")  # Huynh 2026 Table 2: theta4 (Q; L/h/70 kg) = 0.0243
    lfdepot <- log(0.30);          label("SC bioavailability F (fraction)")                                                          # Huynh 2026 Table 2: theta6 (F) = 0.30
    ldur    <- log(36 / 24);       label("Zero-order SC absorption duration D1 for an infant (day) - paper reports 36 hours")        # Huynh 2026 Table 2: theta7 (D1; hours) = 36

    # Allometric exponents - shared between CL and Q, and between Vc and Vp.
    e_wt_cl_q  <- 0.85; label("Shared allometric exponent of WT on CL and Q (unitless)")                                              # Huynh 2026 Methods: "CL and Q by weight (WT^0.85)"
    e_wt_vc_vp <- 1.0;  label("Shared allometric exponent of WT on Vc and Vp (unitless)")                                             # Huynh 2026 Methods: "V_ss by weight (WT^1.0)"

    # Binary covariate effects. Paper's reference category for CL and D1 is the
    # infant cohort (CHILD = 1); to preserve the published structural values,
    # the effects are applied with (1 - CHILD) so the contrast multiplies the
    # adult cohort (CHILD = 0) by 1 + coefficient. The repeat-dose effect on
    # the volumes is applied with (CYCLE > 1).
    e_child_cl      <- 0.69; label("Adult-vs-infant fractional effect on CL (adult CL is 1.69x infant)")                              # Huynh 2026 Table 2: theta9 (adult factor CL) = 1.69
    e_child_dur     <- 1.79; label("Adult-vs-infant fractional effect on D1 (adult D1 is 2.79x infant)")                              # Huynh 2026 Table 2: theta10 (adult factor D1) = 2.79
    e_repdose_vc_vp <- 0.49; label("Repeat-dose fractional effect on Vc and Vp (Vss is 1.49x for second and later doses)")            # Huynh 2026 Table 2: theta8 (repeat dose factor Vc + Vp) = 1.49

    # Inter-individual variability (Huynh 2026 Table 2 BSV section).
    # Paper reports IIV as percent CV; convert to internal log-scale variance
    # via omega^2 = log(CV^2 + 1).
    # IIV on Vc + Vp = 35.7% -> omega^2 = log(0.357^2 + 1) = 0.11991. The paper
    # reports a single IIV applied jointly to Vc and Vp (i.e., a single eta in
    # the source NONMEM control stream assigned to both V2 and V3); a single
    # shared eta is added to both lvc and lvp inside model() so the volumes
    # vary proportionally and preserve the Vc:Vp ratio.
    # IIV on CL = 29.3% -> omega^2 = log(0.293^2 + 1) = 0.08234.
    # IIV on D1 = 10.0% -> omega^2 = log(0.100^2 + 1) = 0.00995.
    etalvc  ~ 0.11991  # Huynh 2026 Table 2: IIV on Vc + Vp = 35.7% (CV%); shared between vc and vp
    etalcl  ~ 0.08234  # Huynh 2026 Table 2: IIV on CL = 29.3% (CV%)
    etaldur ~ 0.00995  # Huynh 2026 Table 2: IIV on D1 = 10.0% (CV%)

    # Combined additive + proportional residual error (Huynh 2026 Methods,
    # Population pharmacokinetic analysis; Table 2 epsilon block).
    propSd <- 0.226; label("Proportional residual error (fraction)")                                                                  # Huynh 2026 Table 2: proportional epsilon = 22.6%
    addSd  <- 0.31;  label("Additive residual error (ug/mL)")                                                                         # Huynh 2026 Table 2: additive epsilon = 0.31 ug/mL
  })

  model({
    # Categorical-covariate multipliers.
    age_cl      <- 1 + e_child_cl      * (1 - CHILD)
    age_dur     <- 1 + e_child_dur     * (1 - CHILD)
    repdose_vss <- 1 + e_repdose_vc_vp * (CYCLE > 1)

    # Individual PK parameters with allometric WT scaling and categorical effects.
    cl     <- exp(lcl + etalcl)  * (WT / 70)^e_wt_cl_q  * age_cl
    vc     <- exp(lvc + etalvc)  * (WT / 70)^e_wt_vc_vp * repdose_vss
    vp     <- exp(lvp + etalvc)  * (WT / 70)^e_wt_vc_vp * repdose_vss
    q      <- exp(lq)            * (WT / 70)^e_wt_cl_q
    fdepot <- exp(lfdepot)
    dur1   <- exp(ldur + etaldur) * age_dur

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Short-term zero-order SC absorption. Each dose to 'depot' releases drug
    # into 'central' at a constant bioavailable rate F * Dose / D1 over the
    # absorption window of length D1. Outside the absorption window, or before
    # any depot dose has been given, kzero is held at 0 so the depot machinery
    # is inert. IV doses are given to 'central' and bypass this mechanism.
    # The default-and-conditional-upgrade pattern (rather than compute-then-
    # zero) is required because podo(depot) returns NA before the first depot
    # dose, and an NA seed propagates through the ODE solver.
    kzero <- 0.0
    if (tad(depot) <= dur1) kzero <- fdepot * podo(depot) / dur1

    d/dt(depot)       <- -kzero
    d/dt(central)     <-  kzero - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-          k12 * central - k21 * peripheral1

    # Concentration: dose in mg, vc in L gives central / vc in mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
