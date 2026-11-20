Zitzmann_2026_bemnifosbuvir_ruzasvir <- function() {
  description <- paste(
    "Multiscale HCV viral-kinetic model of once-daily bemnifosbuvir 550 mg",
    "plus ruzasvir 180 mg in adults with chronic hepatitis C (Zitzmann 2026,",
    "final model B24). The age-structured multiscale model (target cells,",
    "infected cells, intracellular HCV RNA, free virus) is encoded through",
    "its published analytic solution for the plasma viral load, a",
    "tri-exponential decline after a pharmacological delay tau, together",
    "with the analytic solution for serum ALT released by dying infected",
    "cells. The combination acts with constant effectiveness (no PK): it",
    "blocks intracellular HCV RNA production (eps_alpha) and virion",
    "assembly/secretion (eps_s), and multiplies intracellular RNA",
    "degradation by kappa; replication-complex decay gamma is fixed to 0.",
    "Compensated cirrhosis (F4 vs F0-F3) shifts eps_alpha and eps_s on the",
    "logit scale. Time is days since treatment start; the solution describes",
    "the on-treatment decline only. The paper does not report the",
    "random-effect variances or residual-error SDs, so every eta and both",
    "residual SDs are fixed to 0 and simulations are typical-value.",
    sep = " "
  )
  reference <- paste(
    "Zitzmann C, Ribeiro RM, Marc A, Zhou XJ, Perelson AS. Multiscale",
    "Modeling of an 8-Week Combination Regimen of Bemnifosbuvir and Ruzasvir",
    "in Patients with Chronic Hepatitis C Virus Infection. Clin Pharmacol",
    "Ther. 2026 (epub 2026-09-03). doi:10.1002/cpt.70454.",
    "Model equations from Supporting Information S1_Text (eqs 1-4); final",
    "estimates from Table 3 and Table 2 (model B24) of the main text and",
    "Table S5 of S1_Text.",
    sep = " "
  )
  vignette <- "Zitzmann_2026_bemnifosbuvir_ruzasvir"

  units <- list(
    time = "day",
    dosing = "n/a (no dosing events; drug effect is constant from time 0, the start of treatment)",
    concentration = "log10(IU/mL) for log10_viral_load; log(U/L) for lnalt"
  )

  covariateData <- list(
    DIS_CIRRHOSIS = list(
      description = "Compensated cirrhosis indicator: 1 = fibrosis stage F4 (compensated cirrhosis), 0 = fibrosis stages F0-F3 (no cirrhosis).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (F0-F3, non-cirrhotic; 181 of 215 patients, 84.2 percent)",
      notes = "Time-fixed per subject. Shifts the logit of eps_alpha and eps_s (Monolix categorical covariate on a logit-normal parameter). The paper reports the F0-F3 and F4 group values directly (Table 3, model B24); the shifts are the logit differences of those values. All F4 patients in the trial had compensated (Child-Pugh A) cirrhosis.",
      source_name = "Fibrosis stage group (F0-F3 vs F4)"
    )
  )

  covariatesDataExcluded <- list(
    HCV_GT2 = list(
      description = "HCV genotype-2 indicator (1 = genotype 2).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (genotype 1 was the reference in the G-model screen)",
      notes = "Genotype (GT1-GT4) was screened as a categorical covariate on eps_alpha, eps_s, kappa and gamma (S1_Text Table S2, G models) and on delta (Table S7, K24) but did not improve BICc and is not in the final model B24."
    ),
    HCV_GT3 = list(
      description = "HCV genotype-3 indicator (1 = genotype 3).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (genotype 1 was the reference in the G-model screen)",
      notes = "Screened as above (Tables S2 and S7); not retained. The individual (empirical Bayes) parameters still differ between GT1 and GT3 (Figure 3a), which the vignette reproduces from the S2_Text table."
    ),
    HCV_GT4 = list(
      description = "HCV genotype-4 indicator (1 = genotype 4).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (genotype 1 was the reference in the G-model screen)",
      notes = "Screened as above (Tables S2 and S7); not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 215L,
    n_studies = 1L,
    age_range = "20-85 years",
    age_mean = "51 years",
    sex_female_pct = 53.5,
    race_ethnicity = "White 72.6%, Asian 20.5%, Black/African American 2.8%, other 4.2%",
    disease_state = "DAA-naive adults with chronic HCV infection, genotype 1 (70.2%), 2 (1.9%), 3 (27.0%) or 4 (0.9%); fibrosis F0 22.8%, F1 27.9%, F2 23.7%, F3 9.8%, compensated cirrhosis F4 15.8%.",
    dose_range = "Bemnifosbuvir 550 mg plus ruzasvir 180 mg orally once daily for 8 weeks.",
    regions = NA_character_,
    notes = "Per-protocol, treatment-adherent subset (215 of 275 enrolled, adherence from pill counts and plasma drug levels) of the phase 2 single-arm trial NCT05904470 (Table 1). HCV RNA LLOQ 15 IU/mL. SVR12 97.7% (210/215)."
  )

  ini({
    # Viral-kinetic parameters (Table 3; days throughout)
    lrbase <- log(10^6.48)
    label("Baseline (pre-treatment) HCV RNA V0 (IU/mL)") # Table 3: V0 = 6.48 log10 IU/mL (RSE 0.79%); the S2_Text table reports it as Log10(V0)
    ldelta <- log(0.269)
    label("Loss rate of infected cells delta (1/day)") # Table 2 model B24: delta = 0.269 [0.245-0.296] /day; Table 3 rounds to 0.27 (RSE 4.90%)
    lc <- log(19.62)
    label("Plasma virion clearance rate c (1/day)") # Table 3: c = 19.62 /day (RSE 4.90%)
    ltau <- log(1.8 / 24)
    label("Pharmacological delay tau before the drug effect starts (day)") # Table 3: tau_v = 1.8 h (RSE 4.82%) = 0.075 day
    lalpha_r <- fixed(log(40))
    label("Maximum intracellular HCV RNA production rate alpha_R (RNA copies/day)") # Table 3: alpha_R = 40 /day, fixed (Methods: fixed from prior studies)
    lrho <- fixed(log(8.18))
    label("Virion assembly and secretion rate rho (1/day)") # Table 3: rho = 8.18 /day, fixed
    lmu <- fixed(log(1))
    label("Baseline intracellular HCV RNA degradation rate mu (1/day)") # Table 3: mu = 1 /day, fixed

    # Modes of action (Table 3 / Table 2 model B24)
    logiteps_alpha <- logit(0.9859)
    label("Drug effectiveness in blocking HCV RNA production eps_alpha, F0-F3 (fraction)") # Table 3: eps_alpha F0-F3 = 0.9859 (RSE 0.21%)
    e_dis_cirrhosis_eps_alpha <- logit(0.9935) - logit(0.9859)
    label("Logit-scale shift of eps_alpha for compensated cirrhosis F4 (unitless)") # Table 3: eps_alpha F4 = 0.9935 (RSE 0.22%), P = 0.00913
    logiteps_s <- logit(0.9978)
    label("Drug effectiveness in blocking virion assembly/secretion eps_s, F0-F3 (fraction)") # Table 3: eps_s F0-F3 = 0.9978 (RSE 0.04%)
    e_dis_cirrhosis_eps_s <- logit(0.9904) - logit(0.9978)
    label("Logit-scale shift of eps_s for compensated cirrhosis F4 (unitless)") # Table 3: eps_s F4 = 0.9904 (RSE 0.37%), P = 2.4e-5
    lkappa <- log(1.705)
    label("Drug-induced enhancement factor kappa of intracellular HCV RNA degradation (unitless)") # Table S5 model B24: kappa = 1.705 [1.516, 1.963]; Table 3 rounds to 1.70 (RSE 6.70%)
    gamma_rc <- fixed(0)
    label("Drug-induced replication-complex decay rate gamma (1/day)") # Table 3: gamma = 0 (model B24 fixes gamma to 0)

    # ALT sub-model (Table 3)
    lkin_alt <- log(5.23)
    label("Background ALT production rate s (U/L/day)") # Table 3: s = 5.23 (RSE 3.46%); printed unit 'IU/L', a rate per day by eq (3)
    lkout_alt <- log(0.30)
    label("ALT clearance rate c_A (1/day)") # Table 3: c_A = 0.30 /day (RSE 4.29%)
    lrbase_alt <- log(31.57)
    label("Baseline (pre-treatment) ALT A0 (U/L)") # Table 3: A0 = 31.57 U/L (RSE 5.68%)

    # Between-subject variability: the S2_Text individual estimates vary on
    # all ten estimated parameters, but no variance is reported.
    etalrbase ~ fixed(0)
    etaldelta ~ fixed(0)
    etalc ~ fixed(0)
    etaltau ~ fixed(0)
    etalogiteps_alpha ~ fixed(0)
    etalogiteps_s ~ fixed(0)
    etalkappa ~ fixed(0)
    etalkin_alt ~ fixed(0)
    etalkout_alt ~ fixed(0)
    etalrbase_alt ~ fixed(0)

    # Residual error: additive on log10 HCV RNA and on log ALT (Methods: logs
    # of the data were fitted); the SDs are not reported.
    addSd_log10_viral_load <- fixed(0)
    label("Additive residual SD on log10 HCV RNA (log10 IU/mL)") # not reported in the paper
    addSd_lnalt <- fixed(0)
    label("Additive residual SD on log ALT (log U/L)") # not reported in the paper
  })

  model({
    v0 <- exp(lrbase + etalrbase)
    delta <- exp(ldelta + etaldelta)
    c <- exp(lc + etalc)
    tau <- exp(ltau + etaltau)
    alpha_r <- exp(lalpha_r)
    rho <- exp(lrho)
    mu <- exp(lmu)
    eps_alpha <- expit(logiteps_alpha + e_dis_cirrhosis_eps_alpha * DIS_CIRRHOSIS + etalogiteps_alpha)
    eps_s <- expit(logiteps_s + e_dis_cirrhosis_eps_s * DIS_CIRRHOSIS + etalogiteps_s)
    kappa <- exp(lkappa + etalkappa)
    gamma <- gamma_rc

    kin_alt <- exp(lkin_alt + etalkin_alt)
    kout_alt <- exp(lkout_alt + etalkout_alt)
    a0 <- exp(lrbase_alt + etalrbase_alt)

    # S1_Text eq (2) and the definitions under it: N is the pre-treatment
    # steady-state intracellular-RNA scaling, A_R the residual RNA production
    # rate and B the on-treatment intracellular RNA loss rate.
    n_ss <- rho * (alpha_r + delta) / (delta * (rho + mu + delta))
    a_r <- alpha_r * (1 - eps_alpha)
    b_rna <- rho * (1 - eps_s) + kappa * mu
    pref <- (1 - eps_s) * c * rho / n_ss
    term3 <- a_r / ((b_rna - gamma) * delta * (c - delta - gamma))
    term2 <- (n_ss / rho - a_r / ((b_rna - gamma) * delta)) / (c - b_rna - delta)
    c3_vl <- pref * term3
    c2_vl <- pref * term2
    c1_vl <- 1 - c2_vl - c3_vl

    # Viral load is V0 until the pharmacological delay has elapsed; the three
    # coefficients sum to 1, so the solution is continuous at t = tau.
    t_eff <- (t - tau) * (t > tau)
    vl <- v0 * (c1_vl * exp(-c * t_eff) + c2_vl * exp(-(b_rna + delta) * t_eff) +
      c3_vl * exp(-(gamma + delta) * t_eff))

    # S1_Text eq (4): ALT relaxes from A0 to the uninfected level s / c_A as
    # infected cells are lost. The equation carries no delay term.
    alt_n <- kin_alt / kout_alt
    alt <- alt_n * (1 - (a0 - alt_n) / (alt_n * (kout_alt - delta)) *
      (delta * exp(-kout_alt * t) - kout_alt * exp(-delta * t)))

    log10_viral_load <- log10(vl)
    lnalt <- log(alt)

    log10_viral_load ~ add(addSd_log10_viral_load)
    lnalt ~ add(addSd_lnalt)
  })
}
