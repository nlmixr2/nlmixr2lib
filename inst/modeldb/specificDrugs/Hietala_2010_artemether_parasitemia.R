# Combined artemether (ARM) + dihydroartemisinin (DHA) population PK and
# semimechanistic parasite-dynamics PD model in Tanzanian children with
# uncomplicated falciparum malaria (Hietala 2010, Antimicrob Agents Chemother
# 54:4780-4788; doi:10.1128/AAC.00252-10). The PK structure is the same Table 1
# joint parent-metabolite model as modellib('Hietala_2010_artemether'); the PD
# parameters are from Table 3 (symptomatic-patient column). Lumefantrine PK
# is intentionally absent from this file because Hietala 2010 found "the
# introduction of an effect of LUM did not significantly improve the model fit"
# and the LUM PD-effect parameters could not be estimated; the companion
# Hietala_2010_lumefantrine PK model handles LUM in isolation.

Hietala_2010_artemether_parasitemia <- function() {
  description <- paste(
    "Joint artemether (ARM) + dihydroartemisinin (DHA) PK model coupled to a",
    "semimechanistic Plasmodium falciparum parasite life-cycle PD model in",
    "Tanzanian children (ages 1-10 years, weights 8-30 kg) with",
    "uncomplicated falciparum malaria (Hietala 2010). The PK is the same",
    "Table 1 two-compartment ARM (with time-dependent CL/F_ARM via the OCC",
    "dose-occasion covariate) and one-compartment DHA structure as",
    "modellib('Hietala_2010_artemether'). The PD (Table 3) is a five-stage",
    "parasite life-cycle model: parasites mature through tinyrings (PTR),",
    "smallrings (PSR), largerings (PLR), and mature trophozoites /",
    "schizonts (PMT), and parasites killed or injured by drug accumulate",
    "in a spleen compartment (Pspleen) before clearance at a fixed",
    "elimination rate k_spleen = 0.26 / h (Gordi et al. 2002, ref 9).",
    "Replication is encoded as a multiplication factor REPL_p applied to",
    "the PMT -> PTR transit (estimated to 4 in symptomatic patients;",
    "fixed to 1 in asymptomatic children, not encoded in this file).",
    "Drug killing is modelled on all visible developmental stages as",
    "k_ARM = S * log[ARM] and k_DHA = S * log[DHA] with shared slope",
    "S_ARMDHA = 0.073 (Table 3). Visible parasitemia is the sum of the",
    "ring-stage compartments plus the spleen pool. Lumefantrine effect",
    "was tested but not retained in the source paper and is intentionally",
    "absent from this model.",
    sep = " "
  )
  reference <- paste(
    "Friberg Hietala S, Martensson A, Ngasala B, Dahlstrom S, Lindegardh N,",
    "Annerberg A, Premji Z, Farnert A, Gil P, Bjorkman A, Ashton M (2010).",
    "Population pharmacokinetics and pharmacodynamics of artemether and",
    "lumefantrine during combination treatment in children with uncomplicated",
    "falciparum malaria in Tanzania. Antimicrob Agents Chemother",
    "54(11):4780-4788. doi:10.1128/AAC.00252-10. k_spleen fixed value from",
    "Gordi T, Xie R, Huong NV, Huong DX, Karlsson MO, Ashton M (2005). A",
    "semiphysiological pharmacokinetic model for artemisinin in healthy",
    "subjects incorporating autoinduction of metabolism and saturable",
    "first-pass hepatic extraction. Br J Clin Pharmacol 59:189-198.",
    "doi:10.1111/j.1365-2125.2004.02321.x (paper ref 9).",
    sep = " "
  )
  vignette <- "Hietala_2010_artemether_lumefantrine_malaria"
  units <- list(time = "hour", dosing = "mg", concentration = "nmol/L")

  paper_specific_compartments <- c(
    "parasite_tinyrings",
    "parasite_smallrings",
    "parasite_largerings",
    "parasite_matureschizonts",
    "parasite_spleen"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at admission",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear (per-kg) weight normalisation of CL / Q / V parameters",
        "for both ARM and DHA. Hietala 2010 reports all PK parameters",
        "on a per-kg basis (Table 1). Population mean WT 14 kg (range",
        "8-30 kg, Results)."
      ),
      source_name        = "WT"
    ),
    OCC = list(
      description        = "Integer-valued dose-occasion (dose number) indicator",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Values 1..6 identify the dose-occasion of the six-dose Coartem",
        "regimen (doses at 0, 8, 24, 36, 48, and 60 hours). Used to",
        "encode the linear time-dependent increase in apparent oral",
        "artemether clearance: CL/F_ARM = theta1 * (1 + theta2 * (OCC",
        "- 1)), reaching ~10 L/h/kg at OCC = 6, a ~3.4-fold rise over",
        "the regimen (Discussion: enzyme induction). Time-varying within",
        "subject; the per-dose-record OCC propagates through to the",
        "immediately-following sampling occasion."
      ),
      source_name        = "OCC"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 50L,
    n_studies       = 1L,
    n_observations  = 356L,
    age_range       = "1-10 years (mean 4; Results)",
    weight_range    = "8-30 kg (mean 14; Results)",
    sex_female_pct  = 62,
    disease_state   = paste(
      "Acute uncomplicated Plasmodium falciparum malaria with asexual",
      "parasite density 2,000-200,000 / microL at admission and either",
      "axillary temperature >= 37.5 degC or history of fever within 24 h",
      "(Methods)."
    ),
    dose_range      = paste(
      "Coartem (Novartis): 20 mg artemether + 120 mg lumefantrine per",
      "tablet. Weight-based dosing: 5-14 kg -> 1 tablet/dose,",
      "15-24 kg -> 2 tablets/dose, 25-34 kg -> 3 tablets/dose. Six",
      "doses (oral) at 0, 8, 24, 36, 48, and 60 hours."
    ),
    regions            = "Tanzania (Fukayosi Primary Health Care Centre, Bagamoyo District)",
    trial_registration = "ClinicalTrials.gov NCT00336375",
    notes              = paste(
      "PD-model dataset: 356 peripheral parasite counts from the 50",
      "symptomatic patients plus 104 counts from 11 asymptomatic children",
      "from the same coastal region (refs 8, 8a). This file encodes the",
      "symptomatic-patient parameterisation only (REPL_p = 4, A_p = 0,",
      "P_init = 1); the asymptomatic alternative (REPL_a = 1 fixed,",
      "A_a = 1.01 estimated, sine oscillation active) is described in",
      "the vignette Errata for reference but not in model(). Companion",
      "PK-only models for ARM/DHA and LUM:",
      "modellib('Hietala_2010_artemether'),",
      "modellib('Hietala_2010_lumefantrine')."
    )
  )

  ini({
    # === PK parameters (Hietala 2010 Table 1) -- identical to Hietala_2010_artemether ===
    lka      <- fixed(log(1))    ; label("Absorption rate constant ka (1/h); fixed")                  # Hietala 2010 Table 1: ka = 1 (fixed)
    lcl      <- log(2.6)         ; label("Apparent oral artemether clearance at OCC = 1, theta1 (L/h/kg)") # Hietala 2010 Table 1: theta1 = 2.6 (95% CI 1.5-2.6)
    e_occ_cl <- 0.57             ; label("Fractional increase in CL/F_ARM per dose occasion, theta2") # Hietala 2010 Table 1: theta2 = 0.57 (95% CI 0.39-0.75)
    lvc      <- log(5.2)         ; label("Apparent artemether central volume Vc/F_ARM (L/kg)")        # Hietala 2010 Table 1: Vc/F_ARM = 5.2 (95% CI 3.5-7.1)
    lq       <- log(1.4)         ; label("Apparent artemether intercompartmental clearance Q/F_ARM (L/h/kg)") # Hietala 2010 Table 1: Q_ARM = 1.4 (95% CI 1.1-1.8); see Hietala_2010_artemether comment on the Table caption unit typo
    lvp      <- log(41.4)        ; label("Apparent artemether peripheral volume Vp/F_ARM (L/kg)")     # Hietala 2010 Table 1: Vp/F_ARM = 41.4 (95% CI 29.0-58.1)
    lcl_dihydroart  <- log(6.8)         ; label("Apparent DHA metabolite clearance CL/F_DHA (L/h/kg); F_DHA fixed to 1") # Hietala 2010 Table 1: CL/F_DHA = 6.8 (95% CI 5.8-8.0)
    lvc_dihydroart  <- log(3.7)         ; label("Apparent DHA metabolite volume V/F_DHA (L/kg); F_DHA fixed to 1")      # Hietala 2010 Table 1: V/F_DHA = 3.7 (95% CI 2.3-8.7)

    # === PD parameters (Hietala 2010 Table 3, symptomatic-patient column) ===
    # P_init: per-individual initial parasitemia introduced into the tiny-ring
    # (PTR) compartment at simulation start. The point value is fixed at 1
    # parasite/microL with log-normal IIV CV = 119.2% (omega^2 = log(1.192^2
    # + 1) = 1.005). The infection is assumed to have started 4 cycles (4 *
    # MTT = 194 h) before the first sample in symptomatic patients (Results,
    # 'Pharmacodynamic model'), so a typical simulation starts at -194 h
    # with parasite_tinyrings(0) = p_init and runs forward; over 4 cycles
    # the parasite population grows REPL_p ^ 4 ~ 200-fold (Discussion: REPL_p
    # = 4 corresponds to a ~12-fold per-cycle multiplication factor) to reach
    # clinically realistic admission parasitemia of order 10^3-10^5 / microL.
    p_init   <- fixed(1)         ; label("Initial parasitemia introduced into PTR at simulation start (parasites/microL); fixed") # Hietala 2010 Table 3: P_init = 1 (fixed)

    # Mean visible parasite transit time (VPT, h) and mean intraerythrocytic
    # cycle time (MTT, h). Inside model() the visible-stage rate constant is
    # k_VPT = 3 / VPT (three visible compartments PTR -> PSR -> PLR each
    # passing on at the same rate) and the sequestered-to-PTR rate constant
    # is k_IPT = 1 / (MTT - VPT) (Methods / Results, 'Pharmacodynamic model').
    lvpt     <- log(15.5)        ; label("Mean visible parasite transit time VPT (h)")                # Hietala 2010 Table 3: VPT = 15.5 (95% CI 9.7-21.5)
    lmtt     <- log(48.5)        ; label("Mean intraerythrocytic cycle time MTT (h)")                 # Hietala 2010 Table 3: MTT = 48.5 (95% CI 48.0-49.1)

    # Replication factor at merogony for symptomatic patients. REPL_p = 4 is
    # estimated (asymptomatic REPL_a = 1 is fixed for the steady-state
    # asymptomatic cohort, not encoded here). REPL appears multiplicatively
    # in the PMT -> PTR feedback term in the dPTR/dt equation.
    repl     <- 4                ; label("Replication factor at merogony, REPL_p (unitless)")        # Hietala 2010 Table 3: REPL_p = 4 (95% CI 3.6-4.4)

    # Amplitude of the sine oscillation modulating PMT -> PTR transitions.
    # A_p = 0 (fixed) for symptomatic patients per Table 3 ('the sine function
    # produced a significant drop in the OFV in asymptomatic individuals but
    # was not supported by the sparse data in symptomatic patients'). With
    # A_p = 0 the sine modulation collapses to a unity factor.
    amp      <- fixed(0)         ; label("Amplitude of sine oscillation on PMT -> PTR transit, A_p (unitless); fixed") # Hietala 2010 Table 3: A_p = 0 (fixed)

    # Spleen elimination rate. Fixed to Gordi 2002 / Gordi 2005 (paper ref 9).
    lkspleen <- fixed(log(0.26)) ; label("Spleen elimination rate k_spleen (1/h); fixed")             # Hietala 2010 Table 3: k_spleen = 0.26 (fixed; ref 9 = Gordi 2005)

    # Slope of the shared ARM / DHA log-concentration killing effect. The
    # paper writes k_ARM = S * log[ARM] and k_DHA = S * log[DHA] applied to
    # all visible stages (PTR, PSR, PLR, and the new tiny-ring birth rate
    # from PMT); the parallel concentration profiles of ARM and DHA preclude
    # estimating separate effect parameters so S is shared (Methods).
    s_armdha <- 0.073            ; label("Slope of shared ARM / DHA log-concentration kill rate (1/h per log(nM))") # Hietala 2010 Table 3: S_ARMDHA = 0.073 (95% CI 0.049-0.423)

    # === PK IIV (Hietala 2010 Table 1) ===
    etalcl     ~ 0.155649    # Hietala 2010 Table 1: IIV on CL/F_ARM = 41% CV; omega^2 = log(0.41^2 + 1)
    etalcl_dihydroart ~ 0.198101    # Hietala 2010 Table 1: IIV on CL/F_DHA = 47% CV; omega^2 = log(0.47^2 + 1)

    # === PD IIV (Hietala 2010 Table 3) ===
    # P_init carries the only PD-side IIV; CV = 119.2% -> omega^2 = log(1.192^2 + 1).
    etap_init ~ 1.005        # Hietala 2010 Table 3: IIV on P_init = 119.2% CV (95% CI 98-161)

    # === Residual error (Hietala 2010 Tables 1 and 3) ===
    # ARM / DHA plasma concentrations: combined additive + proportional, with
    # the additive components fixed (Table 1).
    propSd     <- 0.61                                                                                 ; label("Proportional residual SD for artemether plasma concentration (fraction)") # Hietala 2010 Table 1: sigma_prop_ARM = 61%
    addSd      <- fixed(2)                                                                             ; label("Additive residual SD for artemether plasma concentration (nM); fixed")    # Hietala 2010 Table 1: sigma_add_ARM = 2 nM (fixed)
    propSd_dihydroart <- 0.82                                                                                 ; label("Proportional residual SD for DHA plasma concentration (fraction)")       # Hietala 2010 Table 1: sigma_prop_DHA = 82%
    addSd_dihydroart  <- fixed(3)                                                                             ; label("Additive residual SD for DHA plasma concentration (nM); fixed")          # Hietala 2010 Table 1: sigma_add_DHA = 3 nM (fixed)

    # Parasitemia residual error. Table 3 reports sigma = 138% as a
    # proportional residual on the log-transformed parasite-density
    # observation (the visual predictive check in Fig. 7 plots log-
    # transformed parasitemia). In nlmixr2's linear-concentration space the
    # equivalent error model is proportional with propSd = 1.38.
    propSd_visibleParasitemia <- 1.38                                                                  ; label("Proportional residual SD for visible parasitemia (fraction; log-scale sigma)") # Hietala 2010 Table 3: sigma_PD = 138% (95% CI 125-258)
  })

  model({
    # === ARM + DHA PK structure (identical to Hietala_2010_artemether) ===
    mw_arm <- 298.4   # g/mol
    mw_dihydroart <- 284.3   # g/mol

    # Time-dependent CL/F_ARM via the OCC dose-occasion covariate
    # (CL/F_ARM = theta1 * (1 + theta2 * (OCC - 1)) * exp(eta)). All
    # PK parameters are per-kg in the source; multiply by WT.
    ka     <-  exp(lka)
    cl     <- (exp(lcl     + etalcl)     * (1 + e_occ_cl * (OCC - 1))) * WT
    vc     <-  exp(lvc)                                                * WT
    q      <-  exp(lq)                                                 * WT
    vp     <-  exp(lvp)                                                * WT
    cl_dihydroart <-  exp(lcl_dihydroart + etalcl_dihydroart)                               * WT
    vc_dihydroart <-  exp(lvc_dihydroart)                                            * WT

    kel     <- cl     / vc
    k12     <- q      / vc
    k21     <- q      / vp
    kel_dihydroart <- cl_dihydroart / vc_dihydroart

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                              k12 * central - k21 * peripheral1
    d/dt(central_dihydroart) <-  kel * central * (mw_dihydroart / mw_arm) - kel_dihydroart * central_dihydroart

    # Plasma concentrations in nM
    Cc     <- 1e6 * central     / vc     / mw_arm
    Cc_dihydroart <- 1e6 * central_dihydroart / vc_dihydroart / mw_dihydroart

    # === PD: parasite life-cycle structure ===
    # Per-individual P_init (with log-normal IIV) and back-transformed PD
    # rate constants. k_VPT = 3 / VPT because PTR -> PSR -> PLR -> (exit to
    # PMT) is three visible-stage transitions each at the same rate;
    # k_IPT = 1 / (MTT - VPT) is the single mean rate from sequestered PMT
    # back to PTR (Methods / Results).
    p_init_i <- p_init * exp(etap_init)
    vpt      <- exp(lvpt)
    mtt      <- exp(lmtt)
    kspleen  <- exp(lkspleen)

    k_vpt <- 3 / vpt
    k_ipt <- 1 / (mtt - vpt)

    # Drug-kill rates from the shared log-concentration slope. The published
    # equation k_ARM = S * log[ARM] is undefined at ARM = 0; the dataset has
    # no parasitemia observations before t = 0 (first dose) where Cc would
    # otherwise be zero, but the simulation runs from t = -4 * MTT to allow
    # the parasite population to build up to admission levels before drug
    # exposure begins. During that pre-drug window the ARM / DHA central
    # compartments are still empty and Cc = 0. To keep the right-hand side
    # finite, the killing rate is gated by an indicator that the relevant
    # plasma concentration is above 1 nM (the natural-log argument crosses
    # unity at Cc = 1 nM, yielding zero killing at that crossover and a
    # smooth onset thereafter); below 1 nM the contribution is forced to 0.
    # This matches the standard NONMEM IF (CARM .GT. 1) form used to make
    # the published equation simulator-safe; the choice is documented in
    # the vignette's Assumptions and deviations section.
    k_arm <- (Cc     > 1) * s_armdha * log(Cc     + (Cc     <= 1))
    k_dihydroart <- (Cc_dihydroart > 1) * s_armdha * log(Cc_dihydroart + (Cc_dihydroart <= 1))

    # Sine modulation on the PMT -> PTR feedback. With amp fixed to 0 for
    # symptomatic patients the factor collapses to unity; the algebraic
    # form is preserved here so the file can be edited to A_p > 0 if needed.
    pi_const <- 3.141592653589793
    sin_mod  <- 1 + amp * sin(2 * pi_const * t / mtt)

    # Visible-stage decay rate, common to PTR, PSR, PLR (Methods: drug-
    # induced killing acts on all visible developmental stages).
    decay_visible <- k_vpt + k_arm + k_dihydroart

    # Five-compartment parasite life cycle (Equations 1-5 of the paper).
    # PTR is fed by sequestered PMT cells released at rate k_IPT and
    # amplified by the merogony multiplication factor REPL.
    d/dt(parasite_tinyrings)       <-  k_ipt * parasite_matureschizonts * repl * sin_mod - decay_visible * parasite_tinyrings
    d/dt(parasite_smallrings)      <-  k_vpt * parasite_tinyrings                        - decay_visible * parasite_smallrings
    d/dt(parasite_largerings)      <-  k_vpt * parasite_smallrings                       - decay_visible * parasite_largerings
    d/dt(parasite_matureschizonts) <-  k_vpt * parasite_largerings                       - parasite_matureschizonts * (k_ipt * sin_mod + k_arm + k_dihydroart)
    d/dt(parasite_spleen)          <- (k_arm + k_dihydroart) * (parasite_tinyrings + parasite_smallrings + parasite_largerings + parasite_matureschizonts) - kspleen * parasite_spleen

    # Initial condition: per-individual P_init parasites/microL placed in the
    # tiny-ring compartment at simulation start (Results: 'PD model ... is
    # initiated by the introduction of a certain parasitemia (P_init) into
    # the compartment denoted by tiny rings'). All other parasite
    # compartments start empty. Note: the published OCR of equation 1 carries
    # P_init as a continuous additive term, which is dimensionally and
    # biologically inconsistent (P_init has units of parasites/microL, not
    # parasites/microL/h); the text describes it as an introduced initial
    # parasitemia, which is the interpretation encoded here.
    parasite_tinyrings(0) <- p_init_i

    # Visible parasitemia observation: ring stages plus drug-injured pool in
    # spleen (Methods: 'Only the parasites in the ring stages and those
    # injured by drugs were assumed to be visible through microscopy').
    visibleParasitemia <- parasite_tinyrings + parasite_smallrings + parasite_largerings + parasite_spleen

    # Residual error
    Cc                 ~ add(addSd)     + prop(propSd)
    Cc_dihydroart             ~ add(addSd_dihydroart) + prop(propSd_dihydroart)
    visibleParasitemia ~ prop(propSd_visibleParasitemia)
  })
}
