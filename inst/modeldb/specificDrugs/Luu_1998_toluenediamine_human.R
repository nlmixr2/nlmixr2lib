Luu_1998_toluenediamine_human <- function() {
  description <- paste(
    "PBPK (whole-body, human). 2,4-toluenediamine (2,4-TDA), a rodent",
    "carcinogen leached from the polyester-urethane foam cover of",
    "Meme/Replicon silicone breast implants (Luu, Hutter and Bushar",
    "1998, Environ Health Perspect). Human scale-up of the rat model:",
    "five flow-limited tissue compartments (slowly-perfused,",
    "richly-perfused, liver, kidney, GI tract) for the parent 2,4-TDA",
    "and for its liver metabolite 4-AAT (4-N-acetyl-2-amino toluene),",
    "coupled through a blood mix point. Parent 2,4-TDA is extensively",
    "plasma-protein bound (9.7% free); liver metabolism is",
    "Michaelis-Menten and a fixed selectivity (0.32) of metabolised",
    "2,4-TDA becomes 4-AAT. The default scenario is a 4.872 g foam",
    "breast implant releasing 2,4-TDA by zero-order degradation in a",
    "58 kg patient; an intravenous forcing function and first-order oral",
    "release are also available. Human physiology from Table 1 and the",
    "metabolism Vmax scaled from the rat by (58/0.25)^0.7 = 45.27. Used",
    "to predict the lifetime breast-implant exposure and cancer risk.",
    "Deterministic, no random effects. See the Luu_1998_toluenediamine_rat",
    "companion for the calibrated rat model."
  )
  reference <- paste(
    "Luu HMD, Hutter JC, Bushar HF. A physiologically based",
    "pharmacokinetic model for 2,4-toluenediamine leached from",
    "polyurethane foam-covered breast implants. Environ Health Perspect.",
    "1998;106(7):393-400. doi:10.1289/ehp.98106393. Physiologic volumes",
    "and flows from Table 1 (human column); metabolism, excretion,",
    "protein-binding and distribution coefficients from Table 3 (human",
    "column); coupled material-balance ODEs from the Appendix",
    "(equations 1-15)."
  )
  vignette <- "Luu_1998_toluenediamine"
  units <- list(time = "h", dosing = "mol", concentration = "ng/mL")

  # Material-balance states: Cs..Cg are 2,4-TDA tissue concentrations
  # (mol/mL); Ms..Mg are 4-AAT tissue concentrations (mol/mL); the four
  # A* states are cumulative excreted amounts (mol). None map to a
  # canonical compartment role. The blood mix point (Cv, Cai, Mai) is
  # algebraic, not a state.
  paper_specific_compartments <- c(
    "Cs",
    "Cr",
    "Cl",
    "Ck",
    "Cg",
    "Ms",
    "Mr",
    "Ml",
    "Mk",
    "Mg",
    "AurineTDA",
    "AfecesTDA",
    "AurineAAT",
    "AfecesAAT"
  )

  compartmentData <- list(
    Cs = list(analyte = "2,4-TDA", units = "mol/mL", specimen = "tissue", verified = TRUE),
    Cr = list(analyte = "2,4-TDA", units = "mol/mL", specimen = "tissue", verified = TRUE),
    Cl = list(analyte = "2,4-TDA", units = "mol/mL", specimen = "tissue", verified = TRUE),
    Ck = list(analyte = "2,4-TDA", units = "mol/mL", specimen = "tissue", verified = TRUE),
    Cg = list(analyte = "2,4-TDA", units = "mol/mL", specimen = "tissue", verified = TRUE),
    Ms = list(analyte = "4-AAT", units = "mol/mL", specimen = "tissue", verified = TRUE),
    Mr = list(analyte = "4-AAT", units = "mol/mL", specimen = "tissue", verified = TRUE),
    Ml = list(analyte = "4-AAT", units = "mol/mL", specimen = "tissue", verified = TRUE),
    Mk = list(analyte = "4-AAT", units = "mol/mL", specimen = "tissue", verified = TRUE),
    Mg = list(analyte = "4-AAT", units = "mol/mL", specimen = "tissue", verified = TRUE),
    AurineTDA = list(analyte = "2,4-TDA", units = "mol", specimen = "urine", verified = TRUE),
    AfecesTDA = list(analyte = "2,4-TDA", units = "mol", specimen = "faeces", verified = TRUE),
    AurineAAT = list(analyte = "4-AAT", units = "mol", specimen = "urine", verified = TRUE),
    AfecesAAT = list(analyte = "4-AAT", units = "mol", specimen = "faeces", verified = TRUE)
  )

  # Body weight is not a model covariate: the tissue volumes and blood
  # flows are the fixed physiologic values of Table 1, and the dose is
  # supplied directly as a molar amount. WT is recorded only to document
  # the 0.25 kg rat the Table 3 scale-up factor references.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "58 kg patient (the Table 2 human breast-implant case). Not used",
        "in model(): physiology is the fixed Table 1 (human) values and",
        "the implant source is k0 * mi. WT is recorded only because the",
        "human metabolism Vmax was scaled from the rat by",
        "(58 kg / 0.25 kg)^0.7 = 45.27."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1L,
    n_studies = 1L,
    age_range = "adult",
    weight_range = "58 kg (Table 2 breast-implant case)",
    sex_female_pct = 100,
    disease_state = paste(
      "A 58 kg patient carrying a 4.87 g polyester-urethane foam Meme",
      "breast implant for 10 years. The human model is an allometric",
      "scale-up of the calibrated rat model (no human 2,4-TDA kinetic",
      "data were fitted); it predicts lifetime low-dose exposure and",
      "the associated cancer risk."
    ),
    dose_range = "4.872 g foam breast implant (zero-order 2,4-TDA release)",
    regions = NA_character_,
    notes = paste(
      "Toxicokinetic PBPK for 2,4-TDA, a medical-device-leachate",
      "carcinogen, not a therapeutic drug. Human physiologic parameters",
      "from Table 1 (human); distribution/binding coefficients shared",
      "with the rat; metabolism Vmax/Vmmax scaled by (58/0.25)^0.7 =",
      "45.27 (Table 3 human column). All parameters are FIXED; the model",
      "is deterministic."
    )
  )

  ini({
    # -------------------------------------------------------------------
    # Calibrated 2,4-TDA / 4-AAT parameters, female F-344 rat. Source:
    # Luu, Hutter and Bushar 1998 Table 3 (rat column). All FIXED;
    # deterministic model solved via rxSolve. Concentrations are in
    # mol/mL and amounts in mol throughout, matching the paper's units.
    # -------------------------------------------------------------------

    # Michaelis-Menten metabolism (Table 3). Vmax/Vmmax are whole-liver
    # maximal rates (mol/hr); km/kmm are the half-saturation free-liver
    # concentrations (mol/mL). At the low doses modelled the liver free
    # concentration is << km, so metabolism is effectively first order.
    # Human Vmax = rat 0.12 x (58/0.25)^0.7 = 0.12 x 45.27 = 5.43 mol/hr.
    # Table 3 prints 0.54, but that contradicts the paper's own stated
    # scaling rule and the metabolite Vmmax (rat 0.01 -> human 0.46 IS the
    # 45.27-fold scale-up); 0.54 is a transcription typo in the scanned
    # table. The scaling-consistent 5.43 is used; see vignette Errata.
    lvmax <- fixed(log(5.43)); label("2,4-TDA hepatic Vmax (mol/hr)") # rat 0.12 x 45.27 (Table 3 human '0.54' is a typo)
    lkm <- fixed(log(7e-7)); label("2,4-TDA hepatic km (mol/mL)") # Table 3 'km' = 7e-7 mol/mL
    lvmmax <- fixed(log(0.46)); label("4-AAT hepatic Vmax (mol/hr)") # Table 3 'Vmmax' human = 0.46 mol/hr
    lkmm <- fixed(log(3e-7)); label("4-AAT hepatic km (mol/mL)") # Table 3 'kmm' = 3e-7 mol/mL

    # Free fraction and metabolite selectivity (Methods / Table 3 footnote)
    fu <- fixed(0.097); label("Free fraction of 2,4-TDA in blood (unitless)") # Methods: 9.7% of 2,4-TDA free (90.3% protein bound)
    sel <- fixed(0.32); label("Fraction of metabolised 2,4-TDA forming 4-AAT (unitless)") # Appendix 'Sel' = 0.32

    # 2,4-TDA tissue:blood distribution coefficients (Table 3). Dj =
    # tissue concentration / emergent venous blood concentration.
    lDs <- fixed(log(1.0)); label("2,4-TDA slowly-perfused:blood distribution coefficient") # Table 3 = 1.0
    lDr <- fixed(log(1.0)); label("2,4-TDA richly-perfused:blood distribution coefficient") # Table 3 = 1.0
    lDl <- fixed(log(5.0)); label("2,4-TDA liver:blood distribution coefficient") # Table 3 = 5.0
    lDk <- fixed(log(3.0)); label("2,4-TDA kidney:blood distribution coefficient") # Table 3 = 3.0
    lDg <- fixed(log(5.0)); label("2,4-TDA GI tract:blood distribution coefficient") # Table 3 = 5.0

    # 4-AAT tissue:blood distribution coefficients (Table 3).
    lDms <- fixed(log(1.0)); label("4-AAT slowly-perfused:blood distribution coefficient") # Table 3 = 1.0
    lDmr <- fixed(log(0.4)); label("4-AAT richly-perfused:blood distribution coefficient") # Table 3 = 0.4
    lDml <- fixed(log(1.0)); label("4-AAT liver:blood distribution coefficient") # Table 3 = 1.0
    lDmk <- fixed(log(10.0)); label("4-AAT kidney:blood distribution coefficient") # Table 3 = 10.0
    lDmg <- fixed(log(10.0)); label("4-AAT GI tract:blood distribution coefficient") # Table 3 = 10.0

    # Excretion flows and urine/feces distribution coefficients (Table 3).
    lQurine <- fixed(log(90)); label("Urine flow (mL/hr)") # Table 3 'Urine flow' human = 90 mL/hr
    lQfeces <- fixed(log(18)); label("Feces flow (mL/hr)") # Table 3 'Feces flow' human = 18 mL/hr
    lDurine <- fixed(log(0.5)); label("2,4-TDA urine:kidney distribution coefficient") # Table 3 '2,4-TDA Durine' = 0.5
    lDmurine <- fixed(log(10.0)); label("4-AAT urine:kidney distribution coefficient") # Table 3 '4-AAT Durine' = 10.0
    # Feces distribution coefficients are NOT tabulated in the paper;
    # feces excretion is < 0.1% of dose and the validated outputs are
    # insensitive to it (Sensitivity analysis, Results). Set to a neutral
    # 1.0 as a documented nominal; see vignette Errata.
    lDfeces <- fixed(log(1.0)); label("2,4-TDA feces:GI distribution coefficient (nominal, untabulated)") # not reported; nominal 1.0
    lDmfeces <- fixed(log(1.0)); label("4-AAT feces:GI distribution coefficient (nominal, untabulated)") # not reported; nominal 1.0

    # Dosing inputs; overridden per scenario in the vignette via rxSolve.
    # Default scenario = 4.872 g foam breast implant (reproduces Table 2
    # human implant column).
    iv_dose <- fixed(0); label("Intravenous 2,4-TDA dose (mol)") # human was not dosed iv; left at 0
    tinj <- fixed(0.01); label("Duration of the iv forcing function (hr)") # near-bolus; individual Qinj/Cinj/tinj not reported, see Errata
    oral_dose <- fixed(0); label("Oral 2,4-TDA dose released into the GI tract (mol)") # not used for the human implant case
    lka <- fixed(log(0.9)); label("First-order oral release rate constant (1/hr)") # Results: 0.015 /min = 0.9 /hr
    k0 <- fixed(3.001e-11); label("Zero-order implant degradation rate (mol/hr/g foam)") # 88 ng/g/day / 122.17 g/mol / 24 = 3.001e-11 mol/hr/g
    mi <- fixed(4.872); label("Implant foam mass (g)") # Table 2 human case = 4.872 g Meme foam implant

    # No residual-error model is reported for this deterministic PBPK.
    # FIXED placeholder, carries no information from the paper; see Errata.
    propSd <- fixed(0.1); label("Proportional residual error placeholder (fraction)") # not reported in Luu 1998; placeholder only
  })

  model({
    # ===================================================================
    # 0. Back-transforms (calibrated Table 3 parameters)
    # ===================================================================
    vmax <- exp(lvmax)
    km <- exp(lkm)
    vmmax <- exp(lvmmax)
    kmm <- exp(lkmm)
    Ds <- exp(lDs)
    Dr <- exp(lDr)
    Dl <- exp(lDl)
    Dk <- exp(lDk)
    Dg <- exp(lDg)
    Dms <- exp(lDms)
    Dmr <- exp(lDmr)
    Dml <- exp(lDml)
    Dmk <- exp(lDmk)
    Dmg <- exp(lDmg)
    Qurine <- exp(lQurine)
    Qfeces <- exp(lQfeces)
    Durine <- exp(lDurine)
    Dmurine <- exp(lDmurine)
    Dfeces <- exp(lDfeces)
    Dmfeces <- exp(lDmfeces)
    ka <- exp(lka)

    # Molar masses for the ng/mL observations.
    MW_TDA <- 122.17 # g/mol, 2,4-toluenediamine C7H10N2
    MW_AAT <- 164.20 # g/mol, 4-N-acetyl-2-amino toluene C9H12N2O

    # Bound:free ratio from the free fraction: (1 + K) = 1 / fu.
    K <- (1 - fu) / fu

    # ===================================================================
    # 1. Fixed physiology, human (Table 1). Flows in mL/hr, volumes in
    #    mL. Qs+Qr+Ql+Qk+Qg = Qi exactly (total cardiac output), which
    #    closes the blood mix balance.
    # ===================================================================
    Qs <- 26520 # Table 1 slowly-perfused flow 26,520 mL/hr
    Qr <- 123200 # Table 1 richly-perfused flow 123,200 mL/hr
    Ql <- 71760 # Table 1 liver flow 71,760 mL/hr
    Qk <- 46800 # Table 1 kidney flow 46,800 mL/hr
    Qg <- 43680 # Table 1 GI flow 43,680 mL/hr
    Qi <- 312000 # Table 1 total cardiac output 312,000 mL/hr

    Vs <- 13562 # Table 1 slowly-perfused volume 13,562 mL
    Vr <- 40672 # Table 1 richly-perfused volume 40,672 mL
    Vl <- 1419 # Table 1 liver volume 1,419 mL
    Vk <- 243 # Table 1 kidney volume 243 mL
    Vg <- 948 # Table 1 GI volume 948 mL

    # ===================================================================
    # 2. Blood mix point (algebraic). Cv is the flow-weighted free venous
    #    blood concentration (Appendix eq 14); Cai is the arterial input
    #    after the implant/iv source is added at the mix point (eq 1); Mai
    #    is the 4-AAT mix-point concentration, which is also the arterial
    #    4-AAT input (eq 15; 4-AAT is not protein bound).
    # ===================================================================
    Lambda <- (t < tinj) * (iv_dose / tinj) # iv forcing (Appendix eqs 2-3)
    Cv <- (Qs * Cs / Ds + Qr * Cr / Dr + Ql * Cl / Dl +
      Qk * Ck / Dk + Qg * Cg / Dg) / Qi # eq 14
    Cai <- (Qi * Cv * (1 + K) + k0 * mi + Lambda) / (Qi * (1 + K)) # eq 1
    Mai <- (Qs * Ms / Dms + Qr * Mr / Dmr + Ql * Ml / Dml +
      Qk * Mk / Dmk + Qg * Mg / Dmg) / Qi # eq 15

    # Michaelis-Menten rates (concentration rates, mol/mL/hr). The free
    # liver concentration is Cl/Dl (2,4-TDA) and Ml/Dml (4-AAT).
    metTDA <- (vmax / Vl) * (Cl / Dl) / (km + Cl / Dl) # 2,4-TDA metabolised in liver
    metAAT <- (vmmax / Vl) * (Ml / Dml) / (kmm + Ml / Dml) # 4-AAT metabolised in liver

    # ===================================================================
    # 3. Material-balance ODEs for 2,4-TDA (Appendix eqs 4-8). The
    #    (1 + K) factor scales the perfusion terms because blood carries
    #    total (free + bound) 2,4-TDA while Cj are free-basis.
    # ===================================================================
    d / dt(Cs) <- (Qs / Vs) * Cai * (1 + K) - (Qs * Cs / (Vs * Ds)) * (1 + K)
    d / dt(Cr) <- (Qr / Vr) * Cai * (1 + K) - (Qr * Cr / (Vr * Dr)) * (1 + K)
    d / dt(Cl) <- (Ql / Vl) * Cai * (1 + K) -
      (Ql * Cl / (Vl * Dl)) * (1 + K) - metTDA
    d / dt(Ck) <- (Qk / Vk) * Cai * (1 + K) -
      (Qk * Ck / (Vk * Dk)) * (1 + K) - (Qurine * Durine * Ck) / Vk
    d / dt(Cg) <- (Qg / Vg) * Cai * (1 + K) -
      (Qg * Cg / (Vg * Dg)) * (1 + K) +
      (ka / Vg) * oral_dose * exp(-ka * t) - (Qfeces * Dfeces * Cg) / Vg

    # ===================================================================
    # 4. Material-balance ODEs for 4-AAT (Appendix eqs 9-13). 4-AAT is
    #    produced in the liver (sel * metTDA), is not protein bound, and
    #    distributes with its own coefficients Dm*.
    # ===================================================================
    d / dt(Ms) <- (Qs / Vs) * Mai - (Qs * Ms / (Vs * Dms))
    d / dt(Mr) <- (Qr / Vr) * Mai - (Qr * Mr / (Vr * Dmr))
    d / dt(Ml) <- (Ql / Vl) * Mai - (Ql * Ml / (Vl * Dml)) +
      sel * metTDA - metAAT
    d / dt(Mk) <- (Qk / Vk) * Mai - (Qk * Mk / (Vk * Dmk)) -
      (Qurine * Dmurine * Mk) / Vk
    d / dt(Mg) <- (Qg / Vg) * Mai - (Qg * Mg / (Vg * Dmg)) -
      (Qfeces * Dmfeces * Mg) / Vg

    # ===================================================================
    # 5. Cumulative excreted amounts (mol): the excretion terms above
    #    times the compartment volume give amount rates.
    # ===================================================================
    d / dt(AurineTDA) <- Qurine * Durine * Ck
    d / dt(AfecesTDA) <- Qfeces * Dfeces * Cg
    d / dt(AurineAAT) <- Qurine * Dmurine * Mk
    d / dt(AfecesAAT) <- Qfeces * Dmfeces * Mg

    # ===================================================================
    # 6. Observations. Reported plasma concentrations are FREE 2,4-TDA
    #    and FREE 4-AAT; the free venous blood concentrations are Cv and
    #    Mai. Convert mol/mL to ng/mL with the molar masses.
    # ===================================================================
    Cc <- Cv * MW_TDA * 1e9 # free 2,4-TDA plasma (ng/mL)
    Caat <- Mai * MW_AAT * 1e9 # free 4-AAT plasma (ng/mL)

    Cc ~ prop(propSd)
  })
}
