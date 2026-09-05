Krishna_2011_anacetrapib_ldlc <- function() {
  description <- paste0(
    "Proportional Emax exposure-response model for low-density lipoprotein ",
    "cholesterol (LDL-C) lowering by the cholesteryl ester transfer protein ",
    "(CETP) inhibitor anacetrapib, alone or with atorvastatin 20 mg once ",
    "daily, in healthy volunteers and patients with dyslipidemia (Krishna ",
    "2011, Table III population-predicted column and Eq. 9). The model is ",
    "algebraic and has no ODE states: LDL-C at steady state is a ",
    "population-specific baseline (separate for volunteers and patients) ",
    "multiplied by a bracket carrying an additive atorvastatin effect, an ",
    "Emax term in the population-predicted anacetrapib trough concentration ",
    "supplied as CSS_ANACETRAPIB, and an interaction parameter gamma that ",
    "scales how much atorvastatin modifies the anacetrapib effect. The ",
    "fitted gamma of 0.99 makes the bracket collapse to the exactly ",
    "multiplicative form, which is the paper's central pharmacological ",
    "conclusion: the two agents act independently. Drive it with the ",
    "companion popPK model Krishna_2011_anacetrapib.R; the HDL-C arm of the ",
    "same paper is Krishna_2011_anacetrapib_hdlc.R."
  )
  reference <- paste0(
    "Krishna R, Bergman AJ, Green M, Dockendorf MF, Wagner JA, Dykstra K. ",
    "Model-based development of anacetrapib, a novel cholesteryl ester ",
    "transfer protein inhibitor. AAPS J. 2011;13(2):179-190. ",
    "doi:10.1208/s12248-011-9254-0. Structural model: Eq. 9 and Table III, ",
    "'Estimate (SE) based on Population-predicted PK' column. PK layer: ",
    "modellib('Krishna_2011_anacetrapib')."
  )
  vignette <- "Krishna_2011_anacetrapib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL", ldlc = "mg/dL")

  covariateData <- list(
    CSS_ANACETRAPIB = list(
      description        = "Population-predicted anacetrapib steady-state trough (24 h post-dose) plasma concentration",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters Eq. 9 as Emax * CSS_ANACETRAPIB / (EC50 + ",
        "CSS_ANACETRAPIB). Same driver and same provenance rule as the ",
        "companion HDL-C model: a POPULATION-predicted trough, obtained by ",
        "solving Krishna_2011_anacetrapib.R under rxode2::zeroRe() at the ",
        "regimen, formulation and prandial state of interest and reading the ",
        "pre-dose concentration at steady state. Table III gives both a ",
        "population-predicted and an individual-predicted fit; the estimates ",
        "'do not differ meaningfully' and this file ships the ",
        "population-predicted column, the paper's selection. Do NOT ",
        "propagate PK inter-individual variability into this column. Set to ",
        "0 for placebo. Note that the LDL-C EC50 of 237 ng/mL is 76% higher ",
        "than the HDL-C EC50 of 135 ng/mL, so LDL-C is the less sensitive of ",
        "the two endpoints."
      ),
      source_name        = "C_trough (C24h)"
    ),
    DIS_HYPERLIP = list(
      description        = "Dyslipidemia diagnosis indicator, 1 = patient with dyslipidemia, 0 = normal healthy volunteer",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal healthy volunteer, NHV)",
      notes              = paste0(
        "Selects between the two baseline LDL-C estimates of Eq. 9, ",
        "rbase_nhv = 107 mg/dL and rbase_pts = 140 mg/dL. Krishna 2011 ",
        "Results: baseline LDL-C is 'about 33 mg/dL (27 to 40 mg/dL) higher ",
        "in patients than healthy subjects', and the Discussion puts this at ",
        "'approximately 30% higher'. The paper's Eq. 9 writes the two ",
        "baselines as a sum of complementary indicator terms ",
        "(BL_NHV * I_NHV + BL_Pts * I_Pts), which is reproduced here with ",
        "I_Pts = DIS_HYPERLIP and I_NHV = 1 - DIS_HYPERLIP. The indicator ",
        "acts on the baseline ONLY: no difference in drug treatment effect ",
        "between patients and volunteers was found for the proportional ",
        "model, which is precisely why the authors preferred it over the ",
        "additive Emax model that required population-specific Emax values. ",
        "The analysis set was 474 patients and 72 volunteers."
      ),
      source_name        = "I_Pts / I_NHV"
    ),
    CONMED_ATORVASTATIN = list(
      description        = "Concomitant atorvastatin indicator, 1 = receiving atorvastatin 20 mg once daily, 0 = no atorvastatin",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant atorvastatin)",
      notes              = paste0(
        "Acts twice in Eq. 9: as the direct proportional effect ",
        "e_conmed_atorvastatin_ldlc = -0.442 (a 44.2% LDL-C reduction from ",
        "baseline, 95% CI 42.5% to 46%, which the paper notes is close to a ",
        "published 42.7% for atorvastatin 20 mg), and inside the interaction ",
        "bracket as e_conmed_atorvastatin_emax (gamma) scaling that same ",
        "effect where it modifies the anacetrapib term. ONLY the 20 mg dose ",
        "was studied, in the phase IIb trial: the paper is explicit that it ",
        "'is unclear whether dose/response observed between anacetrapib and ",
        "20 mg atorvastatin will differ meaningfully relative to different ",
        "statins and different doses of statins', and its own higher-dose ",
        "extrapolations borrowed the Mandema 2005 statin dose-response ",
        "model rather than this indicator. Treat the indicator as specific ",
        "to atorvastatin 20 mg once daily and do not reinterpret it as a ",
        "generic statin flag; the register carries CONMED_STATIN and the ",
        "CONMED_STATIN_LI / _MI / _HI intensity strata for that purpose."
      ),
      source_name        = "I_Atorva"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 544L,
    n_observations   = 3078L,
    n_studies        = 5L,
    disease_state    = "Combined phase I and phase IIb dataset of 474 patients with dyslipidemia and 72 normal healthy volunteers (546 subjects); the LDL-C model converged on 544 subjects and 3,078 observations (Table III).",
    ldlc_baseline    = "Typical baseline LDL-C 107 mg/dL (SE 3) in normal healthy volunteers and 140 mg/dL (SE 1) in patients with dyslipidemia, with a shared 25 mg/dL between-subject SD (18% coefficient of variation on the patient baseline).",
    dose_range       = "Anacetrapib doses across the pooled phase I and phase IIb studies; the model was exercised over simulated arms of placebo and 10, 50, 100, 150, 250 and 300 mg once daily, with and without atorvastatin 20 mg once daily. 100 mg once daily as the hot-melt-extruded tablet was selected for phase III.",
    co_medication    = "Atorvastatin 20 mg once daily was the only statin regimen studied, in the phase IIb trial.",
    notes            = paste0(
      "Per-study designs are in Table IB of the Electronic Supplementary ",
      "Material, which is not on disk (see the vignette Errata). Krishna ",
      "2011 publishes no age, weight, sex or race distribution for the ",
      "exposure-response analysis set in the main text, and no demographic ",
      "covariate is in the model. The paper also fitted an additive Emax ",
      "model that scored better on -2 log-likelihood but required separate ",
      "population-specific Emax estimates; only the proportional model, ",
      "shipped here, was carried into the simulations."
    )
  )

  ini({
    # =========================================================================
    # Krishna 2011 Table III, "Estimate (SE) based on Population-predicted PK"
    # column, and Eq. 9:
    #
    #   LDL_ij = (BL_NHV * I_NHV_i + BL_Pts * I_Pts_i + eta_BL_i)
    #            * (1 + theta_A * I_Atorva_i
    #                 + Emax * C_i / (EC50 + C_i)
    #                   * (1 + gamma * theta_A * I_Atorva_i))
    #            + eps_ij
    #
    # THE INTERACTION BRACKET IS THE ONE PIECE OF PRINTED NOTATION THAT HAD
    # TO BE RECOVERED. The PDF's symbol font renders the second bracket as
    # "(1 + gA I_Atorva)", which reads literally as (1 + gamma * I_Atorva).
    # It is actually (1 + gamma * theta_A * I_Atorva) -- the font dropped the
    # theta glyph exactly as it did in the Table III row label, which prints
    # "Effect of atorvastatin, 20 mg/d ... A" for theta_A. Two independent
    # confirmations:
    #   (a) The literal reading predicts a 76.6% LDL-C decrease for the
    #       paper's own worked example (50 mg fasted plus atorvastatin) where
    #       the paper publishes 52.4%; the theta_A reading gives 53.5%.
    #   (b) With gamma = 1 the bracket becomes
    #       1 + theta_A + E * (1 + theta_A) = (1 + theta_A) * (1 + E),
    #       i.e. exactly multiplicative independence -- which is what the
    #       paper says gamma = 1 means, in four separate places. The literal
    #       reading has no such collapse.
    # See the vignette source-trace and Assumptions sections.
    #
    # Random effects are ADDITIVE on the natural scale, as printed and as
    # S-PLUS NLME parameterises them, so the omega and sigma rows are read as
    # STANDARD DEVIATIONS in mg/dL (25 mg/dL on a 140 mg/dL baseline is an
    # 18% CV; reading it as a variance would give 5 mg/dL, a 3.6% CV, far too
    # small for population LDL-C). This is the opposite convention from the
    # NONMEM-fitted popPK file, whose Table I omegas are provably variances.
    #
    # Emax is NEGATIVE, so no parameter that can take either sign is
    # log-transformed here; all thetas are on the natural scale, which also
    # keeps every eta mu-referenced as theta + eta.
    # =========================================================================
    rbase_nhv <- 107;  label("Typical baseline LDL-C in normal healthy volunteers (mg/dL)")     # Table III, population-predicted: BL_NHV = 107 mg/dL (SE 3); individual-predicted column gives 103
    rbase_pts <- 140;  label("Typical baseline LDL-C in patients with dyslipidemia (mg/dL)")    # Table III, population-predicted: BL_Pts = 140 mg/dL (SE 1); individual-predicted column gives 141
    emax      <- -0.80; label("Maximum proportional decrease in LDL-C by anacetrapib monotherapy (unitless)")  # Table III, population-predicted: Emax = -0.80 (SE 0.04); individual-predicted column gives -0.78
    lec50     <- log(237); label("Anacetrapib trough concentration giving half-maximal LDL-C effect (ng/mL)")   # Table III, population-predicted: EC50 = 237 ng/mL (SE 25); individual-predicted column gives 240

    e_conmed_atorvastatin_ldlc <- -0.442; label("Proportional effect of atorvastatin 20 mg/d on LDL-C (unitless)")            # Table III, population-predicted: theta_A = -0.442 (SE 0.009); Results: 44.2% reduction (95% CI 42.5% to 46%); individual-predicted column gives -0.445
    e_conmed_atorvastatin_emax <- 0.99;   label("Anacetrapib-atorvastatin interaction factor gamma on the anacetrapib effect (unitless)")  # Table III, population-predicted: gamma = 0.99 (SE 0.06); Results: 95% CI 0.88 to 1.1; individual-predicted column gives 0.95. gamma = 1 is exact pharmacologic independence, gamma > 1 synergy, gamma < 1 sub-additivity

    # ---- IIV, additive on the natural scale and shared across both baselines (Table III) ----
    etarbase ~ 625   # Table III, population-predicted: w_BL = 25 mg/dL (SE 0.9), read as an SD -> variance 25^2 = 625; individual-predicted column gives 24. Table III reports no IIV on Emax, unlike the HDL-C model

    # ---- Residual error (Table III) ----
    addSd <- 16; label("Additive residual error on LDL-C (mg/dL)")   # Table III, population-predicted: sigma = 16 mg/dL (SE 0.2); Eq. 9 adds eps_ij outside the proportional bracket, so the residual is purely additive
  })

  model({
    # Eq. 9. The two baselines are complementary indicator terms sharing one
    # additive random effect, exactly as printed.
    rbase_i <- rbase_nhv * (1 - DIS_HYPERLIP) + rbase_pts * DIS_HYPERLIP + etarbase

    ec50 <- exp(lec50)
    drug <- emax * CSS_ANACETRAPIB / (ec50 + CSS_ANACETRAPIB)

    # The interaction bracket: with e_conmed_atorvastatin_emax (gamma) = 1
    # this whole factor collapses to
    # (1 + e_conmed_atorvastatin_ldlc * CONMED_ATORVASTATIN) * (1 + drug).
    ldl <- rbase_i *
      (1 + e_conmed_atorvastatin_ldlc * CONMED_ATORVASTATIN +
         drug * (1 + e_conmed_atorvastatin_emax * e_conmed_atorvastatin_ldlc * CONMED_ATORVASTATIN))

    ldl ~ add(addSd)
  })
}
