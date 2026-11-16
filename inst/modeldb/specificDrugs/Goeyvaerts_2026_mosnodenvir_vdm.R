Goeyvaerts_2026_mosnodenvir_vdm <- function() {
  description <- paste(
    "QSP (within-host viral dynamics). Immune-response viral dynamic model",
    "(after Clapham 2014) of serum DENV-3 RNA after subcutaneous challenge",
    "with rDEN3delta30 in a controlled human infection model, with",
    "pre-exposure mosnodenvir prophylaxis (Goeyvaerts 2026, final model",
    "fitted to Cohorts 1 and 2 of the Phase 2a study). Four ODE states:",
    "uninfected target cells 'target' (produced at kT*T0, die at kT, infected",
    "at beta*V*T), infected cells 'infected' (die at kI and are removed by",
    "the immune response at alpha*A), free virus 'virus' (produced at",
    "(1 - eps)*pI*I, cleared at kV) and immune response 'immune' (grows at",
    "phi*I*A). Mosnodenvir inhibits virion production with a sigmoid Imax",
    "function eps = Imax*CAV^h/(IC50^h + CAV^h) of the individual",
    "steady-state average plasma concentration CAV (ng/mL), which is",
    "constant over the whole infection course; CAV comes from the",
    "Goeyvaerts_2026_mosnodenvir PopPK model as F*dose/(CL/F*tau).",
    "Time zero is inoculation. All parameters except IC50, h and their IIV",
    "were fixed to the placebo-model estimates. Derived outputs: basic",
    "(R0) and effective (Re) reproductive numbers and the reproduction",
    "minimum inhibitory concentration (RMIC). The interim Cohort-1-only",
    "model used for Cohort 2 dose selection is",
    "Goeyvaerts_2026_mosnodenvir_vdm_c1."
  )
  reference <- paste(
    "Goeyvaerts N, Woot de Trixhe X, Neyens M, Kakuda TN, Van Wesenbeeck L,",
    "Buelens A, Rasschaert F, Durbin A, Ackaert O. Viral Dynamic",
    "Model-Informed Dose Selection of Mosnodenvir Prophylaxis against DENV-3",
    "in a Controlled Human Infection Model. Clin Pharmacol Ther.",
    "2026;120(5):1379-1388. doi:10.1002/cpt.70439"
  )
  vignette <- "Goeyvaerts_2026_mosnodenvir"

  paper_specific_compartments <- c("immune")

  units <- list(
    time = "day",
    dosing = "no dosing (mosnodenvir enters as the CAV covariate)",
    concentration = "copies/mL for virus; cells/mL for target and infected cells; immune cells (arbitrary units) for the immune response; CAV in ng/mL"
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average mosnodenvir plasma concentration over the maintenance dosing interval at steady state (Cavg,SS)",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Constant within subject; 0 for placebo. The paper took individual Cavg,SS predictions from the Goeyvaerts_2026_mosnodenvir PopPK model (Results, Viral dynamic modeling); for a regimen with maintenance dose D every tau hours, Cavg,SS = 1000 * F * D / (CL/F * tau) ng/mL with the subject's F and CL/F from that model.",
      source_name = "Cavg,SS"
    )
  )

  compartmentData <- list(
    target = list(
      analyte = "uninfected DENV target cells (activated-macrophage-like)",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    infected = list(
      analyte = "DENV-infected cells",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    virus = list(
      analyte = "DENV-3 RNA (free virus)",
      units = "copies/mL",
      specimen = "serum",
      verified = TRUE
    ),
    immune = list(
      analyte = "host immune response (latent immune cells)",
      units = "immune cells (arbitrary units)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 54L,
    n_studies = 1L,
    n_observations = 1047L,
    age_range = "18-53 years (Cohort 1 median 32, Cohort 2 median 29; Table 1)",
    weight_range = "48.2-107 kg (Table 1)",
    sex_female_pct = 59.3,
    race_ethnicity = c(White = 50.0, Black = 35.2, Asian = 5.6, AmericanIndian = 1.9, Other = 7.4),
    disease_state = "Healthy dengue-seronegative adults challenged subcutaneously with 3 log10 PFU of the under-attenuated DENV-3 strain rDEN3delta30 (controlled human infection model).",
    dose_range = paste(
      "Cohort 1 (fasted capsule): 600/200, 200/50 or 40/10 mg QD",
      "loading (Day -5 to -1)/maintenance (Day 1-21); Cohort 2 (fed tablet):",
      "800 mg BID loading (Day -2 to -1) then 250 mg QD, or 450 or 250 mg BID",
      "loading then 1200 or 500 mg weekly (Days 1, 8, 15); placebo."
    ),
    regions = "United States (Phase 2a 64281802DNG2002, NCT05048875)",
    notes = paste(
      "Cohort 1 n = 29 (placebo 7, high 10, medium 6, low 6), Cohort 2 n = 25",
      "(placebo 6, daily 6, weekly high 7, weekly low 6) (Table 1). 1,047",
      "serum samples: 69 quantifiable, 71 below LLOQ (1,000 copies/mL, imputed",
      "at 500 copies/mL), 907 undetectable (below LOD 36 copies/mL, M3",
      "censoring). Fitted to natural-log DENV-3 RNA with NONMEM (Laplace)."
    )
  )

  ini({
    lrbase_target <- fixed(log(1e7))
    label("Baseline target cell density T0 (cells/mL)") # Table 2: T0 = Fixed 10^7 cells/mL (from Clapham 2014)
    lkdeg_target <- fixed(log(0.14))
    label("Target cell death rate kT (1/day)") # Table 2: kT = Fixed 0.14 1/day
    lkdeg_infected <- fixed(log(0.14))
    label("Infected cell death rate kI (1/day)") # Table 2: kI = Fixed 0.14 1/day
    lprod_virus <- fixed(log(1e2))
    label("Virion production rate pI (copies/mL per infected cell per day)") # Table 2: pI = Fixed 10^2
    lalpha_immune <- fixed(log(1e-3))
    label("Infected cell removal rate by the immune response alpha (1/(immune cell*day))") # Table 2: alpha = Fixed 10^-3
    lrbase_virus <- fixed(log(10^1.28))
    label("Baseline viral load V0 at inoculation (copies/mL)") # Table 2, C1 and C2 column: V0 = Fixed 10^1.28 copies/mL (placebo-model estimate)
    lbeta_infect <- fixed(log(10^-8.35))
    label("Target cell infection rate beta (1/(virion*day))") # Table 2, C1 and C2 column: beta = Fixed 10^-8.35
    lkdeg_virus <- fixed(log(1.46))
    label("Virion clearance rate kV (1/day)") # Table 2, C1 and C2 column: kV = Fixed 1.46 1/day
    lphi_immune <- fixed(log(10^-2.17))
    label("Immune cell proliferation rate phi (1/(infected cell*day))") # Table 2, C1 and C2 column: phi = Fixed 10^-2.17
    lrbase_immune <- fixed(log(1))
    label("Baseline immune response A0 (immune cells)") # Table 2: A0 = Fixed 1

    lic50 <- log(200)
    label("Cavg,SS giving half-maximal inhibition of virion production IC50 (ng/mL)") # Table 2, C1 and C2 column: IC50 = 200 ng/mL
    lhill <- log(1.55)
    label("Hill coefficient h (unitless)") # Table 2, C1 and C2 column: h = 1.55
    limax <- fixed(log(1))
    label("Maximal inhibition of virion production Imax (fraction)") # Table 2: Imax = Fixed 1

    # IIV was additive on the log10 scale; Table 2 footnote gives
    # CV% = sqrt(omega^2)/|theta| * 100 with theta the log10 estimate. The
    # natural-log variance is (log(10) * CV * |theta|)^2.
    etalrbase_virus ~ fixed(1.3484777) # Table 2, C1 and C2: IIV of log10 V0 = 39.4 %CV of 1.28 -> SD 0.50432 log10 units
    etalphi_immune ~ fixed(0.3715956) # Table 2, C1 and C2: IIV of log10 phi = 12.2 %CV of 2.17 -> SD 0.26474 log10 units
    etalic50 ~ 2.3446155 # Table 2, C1 and C2: IIV of log10 IC50 = 28.9 %CV of log10(200) = 2.30103 -> SD 0.66500 log10 units

    expSd <- 1.31
    label("Additive residual error SD on natural-log DENV-3 RNA (log copies/mL)") # Table 2, C1 and C2 column: additive residual error (SD) = 1.31; Methods: natural-log DENV-3 RNA
  })

  model({
    t0_target <- exp(lrbase_target)
    kdeg_target <- exp(lkdeg_target)
    kdeg_infected <- exp(lkdeg_infected)
    prod_virus <- exp(lprod_virus)
    alpha_immune <- exp(lalpha_immune)
    v0 <- exp(lrbase_virus + etalrbase_virus)
    beta_infect <- exp(lbeta_infect)
    kdeg_virus <- exp(lkdeg_virus)
    phi_immune <- exp(lphi_immune + etalphi_immune)
    a0 <- exp(lrbase_immune)
    ic50 <- exp(lic50 + etalic50)
    hill <- exp(lhill)
    imax <- exp(limax)

    # Zero-order target-cell production holds T at T0 before infection
    # (disease-free equilibrium, as in the R0 expression of the Supplement).
    kin_target <- kdeg_target * t0_target

    # Eq. 3: sigmoid Imax inhibition of virion production by Cavg,SS.
    eps_p <- imax * CAV^hill / (ic50^hill + CAV^hill)

    # Eq. 2 (immune response model) with pI -> (1 - eps_p) * pI.
    d/dt(target) <- kin_target - beta_infect * virus * target - kdeg_target * target
    d/dt(infected) <- beta_infect * virus * target - kdeg_infected * infected -
      alpha_immune * immune * infected
    d/dt(virus) <- (1 - eps_p) * prod_virus * infected - kdeg_virus * virus
    d/dt(immune) <- phi_immune * infected * immune
    target(0) <- t0_target
    infected(0) <- 0
    virus(0) <- v0
    immune(0) <- a0

    # Supplement, Reproductive numbers and RMIC (Jacqmin 2010).
    R0 <- prod_virus * beta_infect * t0_target /
      (kdeg_virus * (kdeg_infected + alpha_immune * a0))
    Re <- (1 - eps_p) * R0
    rmic <- (R0 - 1)^(1 / hill) * ic50

    viral_load <- virus
    viral_load ~ lnorm(expSd)
  })
}
