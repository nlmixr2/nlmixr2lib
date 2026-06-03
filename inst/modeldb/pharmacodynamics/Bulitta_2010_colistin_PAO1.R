Bulitta_2010_colistin_PAO1 <- function() {
  description <- "In vitro (Pseudomonas aeruginosa PAO1; NONMEM VI primary fit). Mechanism-based population pharmacodynamic model for colistin static time-kill experiments across initial inocula of 10^4-10^9 CFU/mL. The bacterial system has three pre-existing subpopulations (susceptible / intermediate / least-susceptible) with different second-order killing rate constants by colistin, plus a lag compartment that initially holds the susceptible cells (no growth or natural death, but subject to colistin killing and signal-molecule inhibition). A first-order rate constant klag transfers cells from the lag compartment into the replicating susceptible compartment. Each replicating subpopulation grows via a saturable Michaelis-Menten-style function parameterised by POPmax (would-be plateau without signal-molecule inhibition) and a per-subpopulation low-density growth half-life t1/2(kg,low CFU); cells die first-order at kd = 0.3 /h. A receptor-occupancy submodel encodes competitive displacement of Mg2+/Ca2+ from outer-membrane LPS sites by colistin (Eq. 1) followed by a steep Hill function (Hill = 10 fixed) that maps the fraction of receptors not occupied by cations into the effective colistin concentration at the target site (Eq. 2). Signal molecules tracking CFUALL with first-order kinetics (kdeg) inhibit both bacterial replication (ImaxRep) and colistin killing (ImaxKill) via a Hill-1 function with the same IC50 (Eqs. 4, 5). Drug concentration Ccolistin (mg/L) and cation concentration Ccations (umol/L) are external time-varying covariates supplied from the data; the model has no human PK. Random effects (eta) are NOT included: the paper reports CV 24% IIV on the three growth half-lives only (jointly via difference-in-VGmax), reflecting between-experiment replicate variability; the packaged model is intended for typical-value simulation."
  reference <- paste(
    "Bulitta JB, Yang JC, Yohonn L, Ly NS, Brown SV, D'Hondt RE, Jusko WJ, Forrest A, Tsuji BT. (2010).",
    "Attenuation of colistin bactericidal activity by high inoculum of Pseudomonas aeruginosa",
    "characterized by a new mechanism-based population pharmacodynamic model.",
    "Antimicrobial Agents and Chemotherapy 54(5):2051-2062.",
    "doi:10.1128/AAC.00881-09.",
    sep = " "
  )
  vignette <- "Bulitta_2010_colistin"
  units <- list(time = "hour", dosing = "mg/L (colistin in broth)", concentration = "log10 CFU/mL (observation); mg/L (drug covariate); umol/L (cation covariate)")

  # Ccolistin and Ccations are experimentally-controlled in-vitro inputs
  # (static colistin concentration in broth + sum of Mg2+ and Ca2+ molar
  # concentration), not epidemiological subject covariates. Declaring them
  # in `depends` marks them as model inputs so the canonical-covariate
  # register check does not apply (they are documented in covariateData
  # below for provenance).
  depends <- c("Ccolistin", "Ccations")

  # Paper-mechanistic compartment names; declared so checkModelConventions()
  # accepts the non-canonical roles.
  paper_specific_compartments <- c("bact_slag", "bact_s", "bact_i", "bact_r", "signal")

  covariateData <- list(
    Ccolistin = list(
      description        = "Colistin concentration in growth medium",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Static concentration in supplemented LB broth held constant over the 24 h experiment (the in vitro PD design supplies colistin once at t = 0 and replication of the experiments confirmed negligible degradation; the paper's Discussion notes 'drug degradation was not the primary reason for the inoculum effect of colistin'). PAO1 was studied at 11 colistin concentrations up to 256 mg/L (64x the LB-broth MIC of 4 mg/L). In-vitro experimental input -- not in inst/references/covariate-columns.md (the canonical register is for human pop-PK covariates and does not apply to this in-vitro PD model).",
      source_name        = "Ccolistin (Bulitta 2010 Methods, paragraph on Time-kill experiments)"
    ),
    Ccations = list(
      description        = "Sum of Mg2+ and Ca2+ molar concentration in the growth medium",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Sum of the supplemented Mg2+ and Ca2+ concentrations in LB broth. Per Bulitta 2010 Table 1 footnote (f) the supplemented concentrations were 0.514 mmol/L Mg2+ and 0.624 mmol/L Ca2+ (sum 1.138 mmol/L = 1138 umol/L), and the sum was used in the receptor-occupancy model. Set to 1138 by default to recover the published fits; set to a smaller value (e.g. the EDTA-derived free-ion estimate of 0.06 mmol/L = 60 umol/L, paper Results paragraph 5) to simulate cation-chelated conditions. In-vitro experimental input -- not in inst/references/covariate-columns.md.",
      source_name        = "Ccations (Bulitta 2010 Methods receptor-occupancy paragraph + Table 1 footnote f)"
    )
  )

  population <- list(
    species             = "in vitro (Pseudomonas aeruginosa, PAO1 reference strain)",
    n_subjects          = NA_integer_,
    n_studies           = 1L,
    organism            = "Pseudomonas aeruginosa PAO1 (genetically characterised reference strain from the R.E.W. Hancock Laboratory, University of British Columbia; LB-broth MIC 4.0 mg/L, MHB MIC 2.0 mg/L; LB agar MIC 2 mg/L)",
    system              = "Static time-kill experiments at 37 C in supplemented cation-adjusted LB broth (25 mg/L Ca2+, 12.5 mg/L Mg2+); 20 mL cultures in constant shaking water bath",
    medium              = "Supplemented LB broth (cation-adjusted)",
    duration            = "24 h with sampling at 0, 0.25, 0.5, 1, 2, 3, 4, 8, 12, 16, 24 h",
    inoculum_range      = "10^4, 10^6, 10^8, and 10^9 CFU/mL",
    mic_values          = c(colistin_LB_broth = "4.0 mg/L", colistin_MHB = "2.0 mg/L", colistin_LB_agar = "2 mg/L"),
    regimens            = "11 colistin concentrations up to 256 mg/L (64x MIC in LB broth) at each initial inoculum; antibiotic-free growth controls; replicate viable counts on drug-free and drug-containing LB agar plates (0, 2, 4, 8, 16 mg/L)",
    notes               = "In-vitro pharmacodynamic study; no human or animal subjects. The model was fit simultaneously in NONMEM VI (results shown in the paper main text and figures) and S-ADAPT (results in Table 1 for cross-method comparison) using the first-order conditional estimation method with interaction in NONMEM and importance-sampling Monte-Carlo expectation-maximization (pmethod = 4 or 8) in S-ADAPT. PAO1 estimates packaged here are the NONMEM VI fit (Table 1, NONMEM column). The paper reports CV 24% inter-experiment variability on the three growth half-lives parameterised as a difference in VGmax (Table 1 footnote b); the packaged model omits these etas and provides typical-value simulation -- see vignette Assumptions and deviations. The same structural model was also fit to two clinical strains URMC1 and URMC2 with strain-specific S-ADAPT parameter estimates; those fits are packaged separately as Bulitta_2010_colistin_URMC1 and Bulitta_2010_colistin_URMC2."
  )

  ini({
    # =============================================================
    # Bacterial growth (low-density half-lives in MIN per Table 1)
    # =============================================================
    # Reparameterised growth function (Harigaya et al. 2009; Meagher
    # et al. 2004 -- ref 34 and 45 in the paper): VGmax = max growth
    # velocity (CFU/(mL*h)); CFUm = bacterial density at 50% of max
    # rate (CFU/mL). The paper re-expresses these as POPmax (the
    # would-be plateau without signal-molecule inhibition) and the
    # per-subpopulation low-density growth half-life t1/2(kg,low CFU).
    # See model() block for the back-conversion VGmax_X = kg_low_X *
    # CFUm and CFUm = POPmax * kd / (kg_low_S - kd).
    t12_kg_low_s <- 22.2
    label("Susceptible-population low-density growth half-life (min; t1/2(kg,low CFU)_S)")  # Bulitta 2010 Table 1, P. aeruginosa PAO1 / NONMEM
    t12_kg_low_i <- 28.7
    label("Intermediate-population low-density growth half-life (min; t1/2(kg,low CFU)_I)")  # Bulitta 2010 Table 1
    t12_kg_low_r <- 62.4
    label("Least-susceptible-population low-density growth half-life (min; t1/2(kg,low CFU)_R)")  # Bulitta 2010 Table 1

    # Natural death rate constant -- FIXED at the literature value
    # from Meagher et al. (ref 45) for a constant-drug-concentration
    # static time-kill setup; this anchors the equilibrium that
    # defines POPmax via VGmax_S = POPmax * kd.
    kd_nat <- fixed(0.3)
    label("First-order natural death rate constant (1/h; kd; FIXED to Meagher 2004 ref 45)")  # Bulitta 2010 Table 1 footnote (a)

    # Growth-lag half-life: cells start in bact_slag (no growth, no
    # natural death) and transit to bact_s at rate klag = ln(2)/t12.
    t12_klag <- 1.31
    label("Growth-lag half-life (h; t1/2(klag); klag = ln(2)/t12_klag)")  # Bulitta 2010 Table 1, NONMEM column

    # Maximum population size (would-be plateau without signal-mol
    # inhibition); footnote (c): the observed plateau is approximately
    # POPmax * (1 - ImaxRep).
    log10_popmax <- 9.59
    label("log10 maximum population size (log10 CFU/mL; POPmax without signal inhibition)")  # Bulitta 2010 Table 1, NONMEM column

    # =============================================================
    # Initial inoculum and subpopulation fractions
    # =============================================================
    # log10 CFUo for the published 10^6 inoculum (PAO1 NONMEM); the
    # user can override at simulation time via the dataset to reproduce
    # the 10^4, 10^8, and 10^9 inocula reported in Table 1.
    log10_cfu0 <- 6.17
    label("log10 initial inoculum (log10 CFUo; default = 10^6-class inoculum value, PAO1 NONMEM)")  # Bulitta 2010 Table 1, Log10 CFUo,6 / NONMEM

    # Fractions of the inoculum that are intermediate / least-susceptible.
    # Footnote-cited interpretation: I and R fractions held constant
    # across CFUo levels (see paper Results paragraph 7 + Discussion).
    log10_fr_i <- -3.43
    label("log10 intermediate-population fraction of initial inoculum (unitless; Log10 FrI)")  # Bulitta 2010 Table 1
    log10_fr_r <- -7.19
    label("log10 least-susceptible-population fraction of initial inoculum (unitless; Log10 FrR)")  # Bulitta 2010 Table 1

    # =============================================================
    # Signal molecules (the 'inoculum effect' driver)
    # =============================================================
    # FrSig: signal-molecule concentration at t = 0 divided by the
    # bacterial concentration at t = 0 (footnote d); the ODE for
    # d/dt(signal) is parameterised so that at steady state the signal
    # molecule numerical value equals CFUALL.
    log10_fr_sig <- -1.89
    label("log10 initial signal-molecule concentration relative to inoculum (ml/CFU; Log10 FrSig)")  # Bulitta 2010 Table 1, NONMEM column
    # IC50: signal-molecule concentration at which Imax (kill and rep)
    # is halved.
    log10_ic50_sig <- 6.17
    label("log10 signal-molecule concentration for 50% maximal inhibition (Log10 IC50)")  # Bulitta 2010 Table 1, NONMEM column
    t12_kdeg <- 0.970
    label("Signal-molecule degradation half-life (h; t1/2(kdeg))")  # Bulitta 2010 Table 1, NONMEM column
    imax_rep <- 0.422
    label("Maximal fractional inhibition of bacterial replication by signal molecules (ImaxRep, unitless)")  # Bulitta 2010 Table 1, NONMEM column
    imax_kill <- 0.992
    label("Maximal fractional inhibition of bacterial killing by signal molecules (ImaxKill, unitless)")  # Bulitta 2010 Table 1, NONMEM column

    # =============================================================
    # Receptor occupancy and bacterial killing
    # =============================================================
    # EC50_rec: fraction of LPS receptor sites NOT occupied by Mg2+
    # or Ca2+ at which the effective colistin concentration reaches
    # 50% of the broth concentration.
    ec50_rec <- 0.537
    label("Fraction of receptors not occupied by Mg2+/Ca2+ giving 50% effective colistin (EC50, unitless)")  # Bulitta 2010 Table 1, NONMEM column
    # Hill coefficient: initially estimated 10-20 then fixed to 10
    # for model stability (footnote e).
    hill_rec <- fixed(10)
    label("Receptor-occupancy Hill coefficient (unitless; gamma; FIXED at 10 per Table 1 footnote e)")  # Bulitta 2010 Table 1 footnote (e)
    # Receptor dissociation constants (footnote g) -- FIXED per
    # Schindler & Osborn 1979 (S. Typhimurium G30A) and Methods
    # Receptor occupancy paragraph. Units umol/L.
    kdiss_cation <- fixed(200)
    label("Receptor dissociation constant for Mg2+/Ca2+ (umol/L; KdCations; FIXED at 200)")  # Bulitta 2010 Methods, Table 1 footnote (g)
    kdiss_colistin <- fixed(0.3)
    label("Receptor dissociation constant for colistin (umol/L; KdColistin; FIXED at 0.3)")  # Bulitta 2010 Methods, Table 1 footnote (g)
    # Colistin molecular mass for mg/L <-> umol/L unit conversion.
    mw_colistin <- fixed(1.163)
    label("Mean colistin A+B molar mass (mg/umol = g/mmol; equivalent to 1163 g/mol; FIXED)")  # Bulitta 2010 Methods, paragraph after Eq. 1

    # Second-order killing rate constants per subpopulation.
    lk2s <- log(5.72)
    label("Susceptible-population second-order killing rate constant (L/(mg*h); k2S)")  # Bulitta 2010 Table 1, NONMEM column
    lk2i <- log(0.369)
    label("Intermediate-population second-order killing rate constant (L/(mg*h); k2I)")  # Bulitta 2010 Table 1, NONMEM column
    lk2r <- log(0.00210)
    label("Least-susceptible-population second-order killing rate constant (L/(mg*h); k2R)")  # Bulitta 2010 Table 1, NONMEM column

    # =============================================================
    # Residual error
    # =============================================================
    # Additive on log10 CFU/mL per Bulitta 2010 Table 1 footnote (h).
    # The observation Cc is the log10 of total viable count plus a 1-
    # CFU/mL floor (matches the Wicha 2017 / Landersdorfer 2018
    # in-vitro PD convention); the additive residual SD applies on
    # the log10 scale directly.
    addSd <- 0.474
    label("Additive residual SD on log10 total viable count (log10 CFU/mL; SD_CFU)")  # Bulitta 2010 Table 1 footnote (h), NONMEM
  })

  model({
    # ---- exponentiate log-transformed structural parameters ----
    k2s <- exp(lk2s)
    k2i <- exp(lk2i)
    k2r <- exp(lk2r)

    # ---- rate-constant derivations from half-lives ----
    # Growth half-lives reported in MIN; convert to /h by * 60 in the
    # numerator (kg = ln(2) / (t12_min/60) = ln(2)*60/t12_min).
    kg_low_s <- log(2) * 60 / t12_kg_low_s
    kg_low_i <- log(2) * 60 / t12_kg_low_i
    kg_low_r <- log(2) * 60 / t12_kg_low_r
    klag <- log(2) / t12_klag
    kdeg <- log(2) / t12_kdeg

    # ---- linear-scale population parameters ----
    popmax  <- 10 ^ log10_popmax
    cfu0    <- 10 ^ log10_cfu0
    fr_i    <- 10 ^ log10_fr_i
    fr_r    <- 10 ^ log10_fr_r
    fr_sig  <- 10 ^ log10_fr_sig
    ic50_sig <- 10 ^ log10_ic50_sig

    # ---- growth-function back-derivation (CFUm and VGmax_X) ----
    # The paper parameterises growth as VGmax_X / (CFUm + CFUALL); we
    # reconstruct CFUm and VGmax_X from POPmax and t1/2(kg,low CFU)_X.
    #
    # At plateau without signal-molecule inhibition, growth rate
    # equals natural-death rate for the susceptible population, so
    # POPmax + CFUm = VGmax_S / kd. Combining with VGmax_S = kg_low_S
    # * CFUm gives CFUm = POPmax * kd / (kg_low_S - kd) and
    # VGmax_X = kg_low_X * CFUm.
    cfum    <- popmax * kd_nat / (kg_low_s - kd_nat)
    vgmax_s <- kg_low_s * cfum
    vgmax_i <- kg_low_i * cfum
    vgmax_r <- kg_low_r * cfum

    # ---- total viable bacteria (Eq. 3) ----
    cfu_all <- bact_slag + bact_s + bact_i + bact_r

    # ---- receptor occupancy (Eq. 1) ----
    # Ccolistin is in mg/L; converted to umol/L via division by
    # mw_colistin (mg/umol). Kd ratio (cation/colistin) is unitless.
    ccolistin_um <- Ccolistin / mw_colistin
    frcations    <- Ccations / (kdiss_cation + Ccations + (kdiss_cation / kdiss_colistin) * ccolistin_um)

    # ---- effective colistin concentration at the target site (Eq. 2) ----
    # Hill function of the fraction of receptors NOT occupied by
    # cations, with Hill exponent 10 (fixed).
    not_occ      <- 1 - frcations
    not_occ_h    <- not_occ ^ hill_rec
    ec50_h       <- ec50_rec ^ hill_rec
    ccolistin_eff <- (not_occ_h / (ec50_h + not_occ_h)) * Ccolistin

    # ---- signal-molecule inhibition (Eqs. 4, 5) ----
    inh_factor <- signal / (ic50_sig + signal)
    inh_kill   <- 1 - imax_kill * inh_factor
    inh_rep    <- 1 - imax_rep * inh_factor

    # ---- per-subpopulation growth and killing rates ----
    growth_factor <- 1 / (cfum + cfu_all)
    growth_s <- inh_rep * vgmax_s * growth_factor
    growth_i <- inh_rep * vgmax_i * growth_factor
    growth_r <- inh_rep * vgmax_r * growth_factor

    kill_s <- inh_kill * k2s * ccolistin_eff
    kill_i <- inh_kill * k2i * ccolistin_eff
    kill_r <- inh_kill * k2r * ccolistin_eff

    # ---- ODE system (Eqs. 6-10) ----
    # Eq. 6: susceptible cells in lag compartment -- no growth, no
    # natural death, subject only to klag transfer and colistin
    # killing under signal-molecule inhibition.
    d/dt(bact_slag) <- (-klag - kill_s) * bact_slag

    # Eq. 7: replicating susceptible cells -- growth, natural death,
    # killing; gain from lag compartment.
    d/dt(bact_s)    <- (growth_s - kd_nat - kill_s) * bact_s + klag * bact_slag

    # Eq. 8: intermediate-susceptibility cells.
    d/dt(bact_i)    <- (growth_i - kd_nat - kill_i) * bact_i

    # Eq. 9: least-susceptible cells.
    d/dt(bact_r)    <- (growth_r - kd_nat - kill_r) * bact_r

    # Eq. 10: signal molecules track CFUALL with first-order kinetics
    # (numerical equivalence at steady state).
    d/dt(signal)    <- kdeg * (cfu_all - signal)

    # ---- initial conditions ----
    # Total inoculum cfu0 is partitioned: intermediate and least-
    # susceptible cells start at their fractions of the inoculum;
    # remaining cells start in the lag compartment.
    cfu_i0     <- fr_i * cfu0
    cfu_r0     <- fr_r * cfu0
    cfu_slag0  <- cfu0 - cfu_i0 - cfu_r0

    bact_slag(0) <- cfu_slag0
    bact_s(0)    <- 0
    bact_i(0)    <- cfu_i0
    bact_r(0)    <- cfu_r0
    signal(0)    <- fr_sig * cfu0

    # ---- observation (Eq. 11) ----
    # log10 of total viable count; the 1 CFU/mL floor is the
    # nlmixr2lib in-vitro PD convention (matches Wicha 2017 /
    # Landersdorfer 2018) for keeping log10 finite when bacteria are
    # driven below 1 CFU/mL. Additive residual error on the log10
    # scale per Bulitta 2010 Table 1 footnote (h).
    cfu_obs <- cfu_all + 1
    Cc      <- log10(cfu_obs)
    Cc      ~ add(addSd)
  })
}
