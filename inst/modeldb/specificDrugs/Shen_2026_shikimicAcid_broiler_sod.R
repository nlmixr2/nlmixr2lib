Shen_2026_shikimicAcid_broiler_sod <- function() {
  description <- paste(
    "Preclinical (broiler chicken). Effect-compartment-link indirect-response",
    "pharmacokinetic-pharmacodynamic model for the effect of shikimic acid (SA)",
    "on serum superoxide dismutase activity (SOD) in 5-week-old yellow-feathered broilers after a",
    "single 50 mg/kg oral (gavage) dose (Shen 2026). A one-compartment",
    "first-order-absorption plasma model drives an effect compartment with",
    "equilibration rate constant ke0, and the effect-site concentration then",
    "acts linearly on an indirect-response turnover model in which SA",
    "stimulates the zero-order production kin. The biomarker baseline is the",
    "drug-free turnover steady state E0 = kin / kout. All states are expressed",
    "per kg body weight (volume in mL/kg, clearance in mL/h/kg, amounts in",
    "ug/kg), so plasma and effect-site concentrations are in ug/mL - the units",
    "in which the linear stimulation coefficient was estimated. IMPORTANT: Shen",
    "2026 reports the pharmacodynamic parameters in full (Table 5) but reports",
    "NO parameter of the one-compartment plasma model it says it fitted, so the",
    "three pharmacokinetic values here are reconstructed from the paper's own",
    "non-compartmental results and do NOT reproduce the published Tmax; see the",
    "ini() comments and the vignette Errata before using the plasma layer",
    "quantitatively. Shen 2026 fitted each bird individually in Phoenix",
    "WinNonlin and reported only the mean and SD of the individual estimates,",
    "so no between-subject variability is available; every parameter is fixed",
    "at the published mean and both residual SDs are fixed at zero. Companion",
    "models for the other two biomarkers of the same experiment are",
    "modellib('Shen_2026_shikimicAcid_broiler_taoc') and",
    "modellib('Shen_2026_shikimicAcid_broiler_mda'); for the same drug in",
    "growing pigs see modellib('Mo_2024_shikimicAcid_pig_oral')."
  )
  reference <- paste(
    "Shen Y, Mo K, Zhao H, Zhang Y, Huang X (2026).",
    "Pharmacokinetic-pharmacodynamic integration of shikimic acid in broilers:",
    "An indirect effect modeling approach. Poultry Science 105:106959.",
    "doi:10.1016/j.psj.2026.106959"
  )
  vignette <- "Shen_2026_shikimicAcid"
  units <- list(time = "h", dosing = "ug/kg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(
      analyte = "shikimic acid", units = "ug", specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "shikimic acid", units = "ug", specimen = "plasma",
      verified = TRUE
    ),
    effect = list(
      analyte = "shikimic acid", units = "ug/mL", specimen = "not applicable",
      verified = TRUE
    ),
    sod = list(
      analyte = "superoxide dismutase (SOD) activity", units = "U/mL", specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species       = "chicken (yellow-feathered broiler)",
    n_subjects    = 8,
    n_studies     = 1,
    age_median    = "5 weeks",
    weight_range  = "mean initial body weight 1.98 +/- 0.11 kg (all 16 birds)",
    sex           = "equal numbers of males and females",
    disease_state = "healthy broilers; no oxidative-stress challenge was applied",
    dose_range    = "shikimic acid 50 mg/kg as a single oral (gavage) dose",
    regions       = "China",
    notes         = paste(
      "Shen 2026 'Animals' and 'Dosing and sample collection': sixteen healthy",
      "5-week-old yellow-feathered broilers, acclimated for 1 week, housed four",
      "per cage, fasted 12 h before dosing, randomised by SAS to an intravenous",
      "or an oral group of n = 8 each. This model describes the ORAL group",
      "(n = 8) only: Shen 2026 states that antioxidant biomarkers were",
      "'evaluated only following oral administration to better reflect",
      "practical exposure conditions', so no pharmacodynamic data exist for the",
      "intravenous arm. Plasma was sampled at 15, 30 and 45 min and 1, 2, 3, 4,",
      "6, 8, 12 and 24 h; serum biomarkers at 0, 1, 2, 4, 6, 8 and 12 h",
      "(Figure 3). By 24 h plasma concentrations were below the limit of",
      "quantification in all birds. The 16-bird intravenous plus oral cohort",
      "supports the non-compartmental results of Tables 3 and 4 from which this",
      "model's plasma layer is reconstructed. Both sexes were studied but no",
      "sex covariate was fitted, and no covariate model of any kind was",
      "reported. No race/ethnicity data apply."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # PLASMA LAYER -- RECONSTRUCTED, NOT PUBLISHED. READ THIS BEFORE USE.
    #
    # Shen 2026 states in "Statistical analysis and PK-PD modeling" that "the
    # pharmacokinetic component was described by a one-compartment model with
    # extravascular input and clearance parameterization", but reports NO
    # parameter of that fit anywhere: Table 5 is explicitly and exclusively the
    # PD parameter table, and Tables 3 and 4 are Phoenix WinNonlin
    # non-compartmental analysis, not the compartmental fit. There is no ka, no
    # V/F, no CL/F and no lag time in the paper or in its supplement (which
    # holds only the analytical stability Tables S1/S2). The three values below
    # are therefore RECONSTRUCTED by us from the paper's own published
    # non-compartmental results, under an operator decision recorded in the
    # vignette Errata. They are not the authors' estimates.
    #
    # Two of the three are exact identities on published numbers:
    #   CL/F = Dose / AUC_po = 50,000 / 132.04 = 378.6731 mL/h/kg   (Table 4)
    #   V/F  = Vz / F        = 350 / 0.4037    = 866.9804 mL/kg     (Tables 3, 4)
    # These give kel = CL/F / (V/F) = 0.43677 /h, which agrees independently
    # with the intravenous CL/Vz = 0.15/0.35 = 0.42857 /h of Table 3 and with
    # the terminal slope of the digitised Figure 2B oral curve (0.43-0.44 /h).
    #
    # ka is then the single remaining degree of freedom and is set so the model
    # reproduces the published Cmax of 32.68 ug/mL (Table 4).
    #
    # WHY Tmax IS NOT REPRODUCED. The published oral triplet
    # (Cmax 32.68, Tmax 3.38 h, AUC 132.04) is INFEASIBLE for any lag-free
    # one-compartment first-order-absorption model, so no choice of parameters
    # can match all three. For that structure Tmax = ln(ka/kel)/(ka - kel) is
    # bounded above by 1/kel, so honouring Tmax = 3.38 h forces
    # kel <= 0.2959 /h; combined with CL/F fixed by the published AUC this
    # forces V/F >= 1279.9 mL/kg, and Cmax is then bounded by (Dose/V)/e =
    # 14.37 ug/mL -- a factor of 2.3 BELOW the observed 32.68. A lag time does
    # not rescue it either: Figure 2B shows a non-zero mean concentration
    # already at 0.25 h (1.68 ug/mL), so no absorption delay is admissible.
    # The mean Figure 2B curve peaks at 4 h with a broad 2-4 h plateau because
    # individual birds peak at 2, 3 or 4 h (Table 4) and several (inset birds 5
    # and 7) show frankly double-peaked absorption, consistent with the crop
    # reservoir the authors invoke in the Discussion. Averaging those profiles
    # produces a shape no single Bateman function has.
    #
    # We therefore honour Cmax and AUC -- the two exposure quantities the
    # downstream indirect-response model actually depends on -- and accept a
    # too-early Tmax (1.30 h modelled vs 3.38 h published). The resulting model
    # matches four published quantities exactly or near-exactly: Cmax 32.68,
    # AUC 132.04, V/F = Vz/F, and a terminal half-life of 1.587 h against the
    # 1.61 h intravenous and 1.80 h oral values of Tables 3 and 4. The vignette
    # shows that the full PK-PD system built on this layer reproduces the
    # observed biomarker time course of Figure 3.
    lka <- fixed(log(1.237772)); label("First-order absorption rate constant ka (1/h)")  # reconstructed: value giving the Table 4 Cmax of 32.68 ug/mL at the CL/F and V/F below; NOT published
    lcl <- fixed(log(378.6731)); label("Apparent clearance CL/F (mL/h/kg)")  # Dose / AUC_po = 50000 ug/kg / 132.04 h*ug/mL (Table 4)
    lvc <- fixed(log(866.9804)); label("Apparent central volume of distribution V/F (mL/kg)")  # Vz / F = 350 mL/kg (Table 3) / 0.4037 (Table 4)
    # ---------------------------------------------------------------------
    # PHARMACODYNAMICS -- published in full, Table 5 ("PD parameters of
    # indirect effect models for T-AOC, SOD and MDA in broilers after oral
    # administration of SA"), as the mean +/- SD of the eight individual fits.
    # The SOD row is carried here; its model was selected on AIC 33.50 +/- 10.11.
    #
    # Shen 2026 defines E0 = kin/kout, and that ratio, 8.98 / 0.31 = 28.97 U/mL, reproduces
    # the 0 h bar of Figure 3B (digitised at 30.13 U/mL). Because Table 5 states no
    # units of its own, that agreement confirms both the transcription and the
    # biomarker units simultaneously - kin carries the biomarker's units per
    # hour, so it can only be right if the units are right.
    #
    # The effect-site link is dCe/dt = ke0 * (Cp - Ce) and the linear
    # stimulation coefficient S enters as the factor (1 + S * Ce), so S carries
    # units of mL/ug (the reciprocal of the plasma concentration units).
    lke0 <- fixed(log(0.89)); label("Effect-compartment equilibration rate constant ke0 (1/h)")  # Table 5, SOD row: ke0 = 0.89 +/- 0.74
    lslope <- fixed(log(0.04)); label("Linear stimulation coefficient S of SA on SOD production (mL/ug)")  # Table 5, SOD row: S = 0.04 +/- 0.02
    lkin <- fixed(log(8.98)); label("SOD zero-order production rate kin (U/mL/h)")  # Table 5, SOD row: kin = 8.98 +/- 4.55
    lkout <- fixed(log(0.31)); label("SOD first-order loss rate constant kout (1/h)")  # Table 5, SOD row: kout = 0.31 +/- 0.16

    # Shen 2026 specifies the error STRUCTURE -- "a multiplicative error model
    # was used for plasma concentration data, whereas an additive error model
    # was used for pharmacodynamic response data" -- but the only magnitudes it
    # gives (0.1 for PK and 1.0 for PD) are explicitly the INITIAL standard
    # deviations handed to the estimator, not final estimates. The structures
    # are therefore encoded faithfully and both magnitudes are fixed at zero.
    propSd <- fixed(0); label("Proportional residual SD on shikimic acid plasma concentration (fraction; final value not reported)")  # Shen 2026 reports only the initial SD of 0.1
    addSd <- fixed(0); label("Additive residual SD on serum SOD (U/mL; final value not reported)")  # Shen 2026 reports only the initial SD of 1.0
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    ke0 <- exp(lke0)
    slope <- exp(lslope)
    kin <- exp(lkin)
    kout <- exp(lkout)

    kel <- cl / vc

    # One-compartment plasma model with first-order absorption. Amounts are
    # per kg body weight (ug/kg) and vc is in mL/kg, so Cc is in ug/mL.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc

    # Effect-compartment link, Shen 2026 "Statistical analysis and PK-PD
    # modeling": dCe/dt = ke0 * (Cp - Ce). The state IS the effect-site
    # concentration in ug/mL (no volume term is involved). It starts at zero
    # because dosing starts from a drug-free state.
    d/dt(effect) <- ke0 * (Cc - effect)

    # Indirect response, Shen 2026 "Statistical analysis and PK-PD modeling":
    # for T-AOC and SOD, SA stimulates the zero-order production,
    #   dE/dt = kin * (1 + S * Ce) - kout * E.
    # The baseline is the drug-free steady state E0 = kin/kout, stated
    # explicitly by Shen 2026.
    d/dt(sod) <- kin * (1 + slope * effect) - kout * sod
    sod(0) <- kin / kout

    Cc ~ prop(propSd)
    sod ~ add(addSd)
  })
}
