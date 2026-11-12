Ooi_2023_bevurogant_pasi <- function() {
  description <- paste(
    "Bounded-integer PK/PD model of the Psoriasis Area and Severity Index",
    "(PASI) in patients with moderate-to-severe plaque psoriasis receiving",
    "oral bevurogant (BI 730357, a RORgammat antagonist) or placebo (Ooi",
    "2023). A probability-scale disease state follows a type I indirect",
    "response whose production is inhibited by an Imax function of the",
    "individual average steady-state bevurogant concentration (CAV); Imax",
    "depends on prior biologic and prior IL-17 inhibitor use. The logit of",
    "the state, minus an individual exponential placebo effect, is the latent",
    "value on the probit scale of the 73 PASI categories (0-72); the observed",
    "integer PASI is floor(73 * pnorm(latent + g * eps)) with an individual",
    "scaling g. The baseline is on the logit scale with a Box-Cox transformed",
    "random effect and a baseline sPGA shift."
  )
  reference <- paste(
    "Ooi QX, Wenning L, Plock N, Hoch M, Kerbusch T, Jacobs A, Zheng W, Ji Q,",
    "Gong J, Shah N, Joshi A, Jeon S, Sachs M, Kemmler M, Seeburger J,",
    "Hilger E, Pfeifer E. Bounded integer model-based analysis of psoriasis",
    "area and severity index in patients with moderate-to-severe plaque",
    "psoriasis receiving BI 730357. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(6):758-769. doi:10.1002/psp4.12948.",
    "Model equations from the Appendix S1 NONMEM control stream (run29).",
    sep = " "
  )
  vignette <- "Ooi_2023_bevurogant"
  units <- list(
    time = "h",
    dosing = "n/a (no dose events; exposure enters as the covariate CAV)",
    concentration = "pasiprobit (latent PASI value on the probit scale); pasi on the 0-73 scale; CAV and IC50 in nmol/L"
  )

  paper_specific_etas <- c("etaRUV", "etapmax")
  depends <- c("CAV")

  covariateData <- list(
    CAV = list(
      description = "Individual average steady-state bevurogant plasma concentration",
      units = "nmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the Imax function (IC50 490 nmol/L). In the analysis it is the individual model-predicted Cav,ss from the population PK model (Ooi_2023_bevurogant): F x dose / (CL x tau), with F and CL the individual absolute bioavailability and clearance, converted to nmol/L at 1877.5 nmol/mg. Time-fixed per subject; 0 for placebo.",
      source_name = "CAV"
    ),
    PRIOR_BIO = list(
      description = "Prior biologic use for psoriasis",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (biologic-naive)",
      notes = "Selects the Imax category together with PRIOR_IL17: Imax 0.411 when PRIOR_BIO = 0, 0.284 when PRIOR_BIO = 1 and PRIOR_IL17 = 0, 0.186 when PRIOR_BIO = 1 and PRIOR_IL17 = 1 (control stream 'IF(BIO.EQ.1.AND.IL17.EQ.0)' / 'IF(BIO.EQ.1.AND.IL17.EQ.1)').",
      source_name = "BIO"
    ),
    PRIOR_IL17 = list(
      description = "Prior IL-17 inhibitor use",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no prior IL-17 inhibitor)",
      notes = "A subset of PRIOR_BIO: PRIOR_IL17 = 1 requires PRIOR_BIO = 1, and the model applies the IL-17 Imax only when both are 1.",
      source_name = "IL17"
    ),
    SPGA_BL = list(
      description = "Baseline static Physician's Global Assessment score",
      units = "(score)",
      type = "count",
      reference_category = "3 (moderate; the most common category)",
      notes = "Additive shift of 0.154 on the logit baseline when SPGA_BL >= 4 (control stream 'IF(BSPGA.EQ.4) BASEBSPGA = THETA(11)'). Only scores 3 and 4 occur in the analysis (moderate-to-severe inclusion); other values are outside the data.",
      source_name = "BSPGA"
    )
  )

  compartmentData <- list(
    latent_pasi = list(
      analyte = "disease state on the probability scale; its logit is the latent PASI value",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 274L,
    n_studies = 1L,
    n_observations = 1554L,
    age_range = "18-75 years",
    sex_female_pct = 30.7,
    disease_state = "Moderate-to-severe plaque psoriasis (baseline PASI >= 12)",
    dose_range = "Placebo; 25, 50, 100 or 200 mg once daily fasted for 12 weeks (part 1); 400 mg once daily or 200 mg twice daily fed for 12 weeks (part 2)",
    regions = "Not reported",
    notes = "Phase II study 1407-0030 (Table 1); PASI assessed on days 1, 8, 15, 29, 57 and 85 within the 12-week analysis window. 69.3% men (Results). Prior biologic and IL-17 inhibitor use and baseline sPGA proportions are not reported."
  )

  ini({
    # Baseline on the logit scale, with a Box-Cox transformed random effect
    logitrbase <- -0.780; label("Logit of the baseline disease state, sPGA 3 (logit)") # Table 3 'Baseline' = -0.780
    e_spga_rbase <- 0.154; label("Additive shift of baseline sPGA 4 on the logit baseline (logit)") # Table 3 'Baseline sPGA effect on baseline' = 0.154
    boxcox_rbase <- 3.03; label("Box-Cox shape parameter of the baseline random effect (unitless)") # Table 3 'Box-Cox coefficient for baseline' = 3.03

    # Placebo effect: typical value 0, individual values only
    pmax <- fixed(0); label("Typical maximum placebo response on the latent scale (probit units)") # Table 3 'Maximum placebo response' = 0

    # Drug effect: Imax per prior-biologic category (probability scale -> logit)
    logitimax_bionaive <- logit(0.411); label("Logit of Imax, no prior biologic use (logit)") # Table 3 'I max for no prior biologic use' = 0.411
    logitimax_nonil17 <- logit(0.284); label("Logit of Imax, prior non-IL-17-inhibitor biologic use (logit)") # Table 3 'I max for prior non-IL-17 inhibitor use' = 0.284
    logitimax_il17 <- logit(0.186); label("Logit of Imax, prior IL-17 inhibitor use (logit)") # Table 3 'I max for prior IL-17 inhibitor use' = 0.186
    lic50 <- log(490); label("CAV at half-maximal inhibition (nmol/L)") # Table 3 'IC50 (nmol/L)' = 490
    lhill <- fixed(log(1)); label("Sigmoidicity of the Imax function (unitless)") # Table 3 'gamma' = 1.00
    lthalfrec <- log(33.8); label("Half-life of the disease-state turnover (day)") # Table 3 'Half-life to reach I max (days)' = 33.8

    # Scaling of the bounded-integer model (SD of the latent value)
    addSd_pasiprobit <- 0.0981; label("Scaling parameter g, SD on the probit scale of the PASI categories (probit units)") # Table 3 'Scaling parameter' = 0.0981

    # IIV (Table 3 SDs; control stream $OMEGA variances). The scaling and Imax
    # random effects are correlated (r = 0.709).
    etaRUV + etalogitimax ~ c(0.395917, 0.81062, 3.30427) # Table 3 IIV scaling 0.629 (SD), Imax 1.82 (SD), correlation 0.709; $OMEGA BLOCK(2) 0.395917 / 0.81062 3.30427
    etalogitrbase ~ 0.0308488 # Table 3 IIV baseline 0.176 (SD); $OMEGA OM3_BASE 0.0308488, Box-Cox transformed in model()
    etapmax ~ 0.0556852 # Table 3 IIV maximum placebo response 0.236 (SD); $OMEGA OM4_PLMAX 0.0556852
  })

  model({
    # Baseline: Box-Cox transformed eta, sPGA shift, logit scale
    eta_bc <- (exp(etalogitrbase * boxcox_rbase) - 1) / boxcox_rbase
    spga_severe <- 0
    if (SPGA_BL >= 4) spga_severe <- 1
    rbase <- expit(logitrbase + e_spga_rbase * spga_severe + eta_bc)

    # Imax by prior-biologic category, logit-scale IIV
    logitimax_typ <- logitimax_bionaive * (1 - PRIOR_BIO) +
      logitimax_nonil17 * PRIOR_BIO * (1 - PRIOR_IL17) +
      logitimax_il17 * PRIOR_BIO * PRIOR_IL17
    imax <- expit(logitimax_typ + etalogitimax)
    ic50 <- exp(lic50)
    hill <- exp(lhill)
    eff <- imax * CAV^hill / (ic50^hill + CAV^hill)

    # Indirect response on the probability scale (days -> h)
    thalfrec <- exp(lthalfrec) * 24
    kout <- log(2) / thalfrec
    kin <- rbase * kout
    latent_pasi(0) <- rbase
    d/dt(latent_pasi) <- kin * (1 - eff) - kout * latent_pasi

    # Individual placebo effect, rising with the same half-life
    pmax_i <- pmax + etapmax
    pl <- pmax_i * (1 - exp(-kout * t))

    # Latent value on the probit scale of the 73 categories: category j
    # (PASI = j) spans qnorm(j / 73) to qnorm((j + 1) / 73)
    pasiprobit <- logit(latent_pasi) - pl
    gi <- addSd_pasiprobit * exp(etaRUV)
    # Median PASI on the continuous 0-73 scale; the integer score is floor()
    pasi <- probitInv(pasiprobit, 0, 73)
    pasiprobit ~ add(gi)
  })
}
