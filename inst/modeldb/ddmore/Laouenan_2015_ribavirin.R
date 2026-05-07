Laouenan_2015_ribavirin <- function() {
  description <- "Hemoglobin turnover (indirect-response) model describing ribavirin-induced anemia in HCV genotype-1 cirrhotic patients on telaprevir- or boceprevir-based triple therapy (Laouenan 2015). Hemoglobin (g/dL) follows a kin/kout indirect-response ODE in which ribavirin inhibits hemoglobin synthesis with an Imax = 1, EC50 form. The ribavirin concentration time-course is reconstructed analytically from per-subject empirical-Bayes regressors (CSS_RBV, K_RBV) supplied as data columns from a separately fitted Laouenan 2015 upstream ribavirin popPK fit; this PD model does not instantiate the PK ODE itself. Distributed in the DDMORE Foundation Model Repository as DDMODEL00000285; the linked publication fits the same equations to 15 ANRS-CO20-CUPIC patients (9 telaprevir, 6 boceprevir)."
  reference <- paste(
    "Laouenan C, Guedj J, Peytavin G, Nguyen TT, Lapalus M, Khelifa-Mouri F,",
    "Boyer N, Zoulim F, Serfaty L, Bronowicki JP, Martinot-Peignoux M,",
    "Lada O, Asselah T, Dorival C, Hezode C, Carrat F, Nicot F, Marcellin P,",
    "Mentre F. (2015).",
    "A Model-Based Illustrative Exploratory Approach to Optimize the Dosing",
    "of Peg-IFN/RBV in Cirrhotic Hepatitis C Patients Treated With Triple Therapy.",
    "CPT Pharmacometrics Syst Pharmacol 4(1):e00008.",
    "doi:10.1002/psp4.8.",
    "DDMORE Foundation Model Repository: DDMODEL00000285.",
    sep = " "
  )
  vignette <- "Laouenan_2015_ribavirin"
  units <- list(time = "day", dosing = "mg", concentration = "g/dL")
  ddmore_id    <- "DDMODEL00000285"
  replicate_of <- NULL

  covariateData <- list(
    CSS_RBV = list(
      description        = "Individual posthoc ribavirin steady-state trough plasma concentration from an upstream popPK fit (regressor input)",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject (time-fixed) empirical-Bayes (modal) estimate from the Laouenan 2015 upstream ribavirin popPK fit. Used together with K_RBV inside model() to reconstruct the individual ribavirin concentration time-course analytically: riba(t) = CSS_RBV * (1 - exp(-K_RBV * t)). The DDMORE bundle's Simulated_Laouenant_2015_CPTPSP_hb_RBV.txt carries this column as `css_mode`; rename `css_mode` -> CSS_RBV before passing the dataset to rxSolve. For typical-trajectory simulations a single CSS_RBV (e.g., the Table-1 cohort median ~3000 ng/mL) suffices.",
      source_name        = "css_mode"
    ),
    K_RBV = list(
      description        = "Individual posthoc ribavirin approach-to-steady-state rate constant from the same upstream popPK fit (regressor input)",
      units              = "1/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject (time-fixed) empirical-Bayes (modal) estimate from the Laouenan 2015 upstream ribavirin popPK fit. Apparent first-order rate constant of the lumped exponential trough-concentration model (NOT a structural elimination rate). Used together with CSS_RBV inside model() in the analytical expression riba(t) = CSS_RBV * (1 - exp(-K_RBV * t)). The DDMORE bundle's Simulated_Laouenant_2015_CPTPSP_hb_RBV.txt carries this column as `k_mode`; rename `k_mode` -> K_RBV before passing the dataset to rxSolve. Bundle range across 15 subjects: 0.013-0.47 day^-1 (approach-to-Css half-lives 1.5-55 days).",
      source_name        = "k_mode"
    )
  )

  population <- list(
    n_subjects     = 15L,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    disease_state  = "Adult HCV genotype 1 patients with cirrhosis (Metavir F4) and prior treatment failure to peg-interferon-alpha plus ribavirin, enrolled in the French ANRS-CO20-CUPIC compassionate-use cohort. Two treatment arms: peg-IFN-alpha2a + ribavirin + telaprevir (TVR; n = 9) and peg-IFN-alpha2a + ribavirin + boceprevir (BOC; n = 6). Baseline hemoglobin median 15.1 g/dL (range 10.8-16.0).",
    dose_range     = "Ribavirin 1000-1200 mg/day (weight-adjusted) plus telaprevir 750 mg q8h or boceprevir 800 mg q8h plus peg-IFN-alpha2a 180 ug/week SC. Doses are not consumed by this PD model; the ribavirin exposure is supplied as the per-subject CSS_RBV / K_RBV regressors from the upstream popPK fit.",
    regions        = "France (multicentre)",
    notes          = "Population n = 15 and treatment-arm split (9 TVR + 6 BOC) confirmed via the Laouenan 2015 publication (PMID 26225222, doi:10.1002/psp4.8) abstract retrieved from PubMed E-utilities. The DDMORE bundle's `Simulated_Laouenant_2015_CPTPSP_hb_RBV.txt` simulated dataset reproduces 15 subject identifiers from the same study. Demographic detail (age range, weight range, sex distribution, race) is reported in the publication's Table 1 but the publication PDF / PMC full text was not accessible from the worktree environment, so the corresponding fields are recorded as NA. The bundle's CAT column (BOC / TVR) records the protease-inhibitor arm but does not enter the PD model equations; the treatment effect on PK is absorbed into the per-subject CSS_RBV / K_RBV upstream-PK regressors."
  )

  ini({
    # Final population estimates from the DDMORE bundle's Output_real_Laouenant_2015_CPTPSP_hb_RBV
    # listing (Monolix mlxtran fit; "Estimation of the population parameters" block,
    # generated 2013-11-06). The publication-reported value of IC50RBV is 7,090 ng/mL
    # (Laouenan 2015 Results); the bundle's Output_real listing reports 8,280 ng/mL. The
    # difference is documented in the validation vignette's Errata and follows the
    # SKILL.md DDMORE-source guidance to use the .lst final estimates verbatim. The
    # publication does not report typical kout,Hb or hb0 values for direct comparison.
    lhb0  <- log(14.3)    ; label("Typical baseline hemoglobin Hb0 (g/dL)")                   # Output_real_Laouenant_2015_CPTPSP_hb_RBV: hb0 = 14.3
    lkout <- log(0.124)   ; label("Hemoglobin elimination (turnover) rate constant kout,Hb (1/day)")  # Output_real_Laouenant_2015_CPTPSP_hb_RBV: Kout = 0.124
    lec50 <- log(8.28e3)  ; label("Ribavirin half-maximal inhibitory concentration on Hb synthesis EC50 (ng/mL)")  # Output_real_Laouenant_2015_CPTPSP_hb_RBV: EC50 = 8.28e+003

    # Inter-individual variability -- Output_real_Laouenant_2015_CPTPSP_hb_RBV
    # "Estimation of the population parameters" block reports omegas as Monolix random-
    # effect SDs on the log scale (default log-normal parameterization for positive-
    # constrained parameters). The internal nlmixr2 ini() variance on the eta is omega^2.
    # The listing shows independent omegas (no BLOCK / correlated random-effect block);
    # the additional correlation matrices in the listing are correlations among the
    # parameter ESTIMATES from the linearized FIM, not structural random-effect
    # correlations.
    etalhb0  ~ 0.00728   # omega_hb0  = 0.0853 -> var = 0.0853^2 = 0.00728
    etalkout ~ 0.147     # omega_Kout = 0.383  -> var = 0.383^2  = 0.147
    etalec50 ~ 0.0906    # omega_EC50 = 0.301  -> var = 0.301^2  = 0.0906

    # Residual error - Output_real_Laouenant_2015_CPTPSP_hb_RBV reports a single
    # additive residual error parameter `a` = 0.737 g/dL on the hemoglobin observation.
    # Monolix's default residual model with only `a` reported is plain additive:
    # y = f + a * epsilon, epsilon ~ N(0, 1). nlmixr2 equivalent: hb ~ add(addSd_hb).
    addSd_hb <- 0.737 ; label("Additive residual error on hemoglobin (g/dL)")  # Output_real_Laouenant_2015_CPTPSP_hb_RBV: a = 0.737
  })

  model({
    # Individual parameters - log-normal etas around the typical values, matching
    # Monolix's default log-normal distribution for the structural parameters Hb0,
    # kout,Hb, and IC50RBV.
    hb0  <- exp(lhb0  + etalhb0)
    kout <- exp(lkout + etalkout)
    ec50 <- exp(lec50 + etalec50)

    # Predicted ribavirin plasma concentration as an analytical exponential approach
    # to steady state, using the per-subject empirical-Bayes regressors CSS_RBV
    # (steady-state trough) and K_RBV (apparent approach-to-Css rate constant)
    # supplied as data columns from the upstream Laouenan 2015 RBV popPK fit. At
    # t = 0 the formula evaluates to 0 (since 1 - exp(0) = 0), which makes hb start
    # at the per-subject baseline hb0 (no drug effect). Equivalent to the mlxtran
    # `if t<=0 then riba=0 else riba = css_mode * (1-exp(-k_mode*t))` block; the
    # conditional is unnecessary because the formula naturally gives riba(0) = 0
    # and simulation always starts at t >= 0.
    riba <- CSS_RBV * (1 - exp(-K_RBV * t))

    # Indirect-response turnover with Imax-style inhibition of hemoglobin synthesis
    # (Imax fixed at 1 per Laouenan 2015 Methods). The synthesis (kin) rate is set
    # to hb0 * kout so that hb = hb0 at steady state with no drug; the elimination
    # rate is kout * hb. The full ODE matches the publication's
    # `dHb/dt = kin,Hb * (1 - Imax * CRBV/(IC50RBV + CRBV)) - kout,Hb * Hb` and the
    # mlxtran `ddt_hb = hb0*Kout*(1-(riba/(riba+EC50))) - Kout*hb` verbatim.
    d/dt(hb) <- hb0 * kout * (1 - riba / (riba + ec50)) - kout * hb
    hb(0)    <- hb0

    # Observation - additive residual error on hemoglobin. The mlxtran model
    # exports `output = hb`; the listing's residual model is plain additive.
    hb ~ add(addSd_hb)
  })
}
