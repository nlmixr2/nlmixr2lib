Phee_2019_colistin_fusidicAcid <- function() {
  description <- "In vitro (Acinetobacter baumannii; ATCC 19606 and five multidrug-resistant clinical isolates). PK/PD model of static time-kill experiments with colistin and fusidic acid, alone and in combination: logistic bacterial growth whose net rate switches from -knet to +knet through a steep (power-20) lag function, sigmoidal Emax kill by each drug with additive combined effect, an empirical time- and concentration-dependent loss of drug effect (a Gompertz-type term, 1 - beta * (1 - exp(-t * sum(tau * C)))) representing resistance development, and categorical interaction terms through which fusidic acid lowers the colistin EC50 and colistin raises the fusidic acid Emax. Strain AB205 (colistin-resistant) and AB12 differ from the reference strains through categorical covariates. Drug concentrations are static covariates (no dosing events); the initial bacterial count is the observed inoculum."
  reference <- "Phee LM, Kloprogge F, Morris R, Barrett J, Wareham DW, Standing JF. Pharmacokinetic-pharmacodynamic modelling to investigate in vitro synergy between colistin and fusidic acid against MDR Acinetobacter baumannii. J Antimicrob Chemother. 2019;74(4):961. doi:10.1093/jac/dky524. Structural equations: Methods Eqs 3-8. Parameter estimates: Table 2."
  vignette <- "Phee_2019_colistin_fusidicAcid"
  units <- list(time = "h", dosing = "mg/L", concentration = "mg/L")

  compartmentData <- list(
    cfu = list(
      analyte = "Acinetobacter baumannii (total viable bacteria)",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CONC_COLISTIN_MGL = list(
      description = "Static (time-invariant) colistin concentration in the time-kill tube. Drives the sigmoidal Emax kill by colistin and its contribution to the time-dependent loss of effect; any colistin above zero also switches on the colistin-on-fusidic-acid Emax interaction.",
      units = "mg/L",
      type = "continuous",
      reference_category = "n/a -- 0 mg/L is the drug-free growth control; the time-kill grid spanned 0.25-8192 mg/L (Methods, Time-kill assays)",
      notes = "Colistin sulphate (Alfa Aesar). Held constant over the 24-h experiment: colistin degradation in broth (0.0175 1/h) did not improve the description of the CFU data and was not included in the final model (Results, Drug degradation experiments). Concentrations are in mg/L; Table 2 prints the EC50 units as ng/mL, which the vignette shows is inconsistent with the rest of the model.",
      source_name = "colistin concentration"
    ),
    CONC_FUS_MGL = list(
      description = "Static (time-invariant) fusidic acid concentration in the time-kill tube. Drives the sigmoidal Emax kill by fusidic acid and its contribution to the time-dependent loss of effect; any fusidic acid above zero also switches on the fusidic-acid-on-colistin EC50 interaction.",
      units = "mg/L",
      type = "continuous",
      reference_category = "n/a -- 0 mg/L is the drug-free growth control; the time-kill grid spanned 1-8192 mg/L (Methods, Time-kill assays)",
      notes = "Fusidic acid (Sigma-Aldrich). Constant in broth over 24 h by LC-MS/MS (Results, Drug degradation experiments).",
      source_name = "fusidic acid concentration"
    ),
    STRAIN_AB205 = list(
      description = "1 = the experiment used A. baumannii clinical isolate AB205 (colistin-resistant; colistin MIC > 256 mg/L, fusidic acid MIC 32 mg/L); 0 = any other strain.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 -- strains AB14 (OXA-23 UK clone I, the strain on which the model was developed), ATCC 19606, AB16 and AB315, which share the reference parameter values; AB12 carries its own indicator (STRAIN_AB12)",
      notes = "AB205 raises the colistin EC50 by (1 + 65.2), raises tau for fusidic acid by (1 + 2.12), lowers beta for fusidic acid monotherapy by (1 - 0.124), and in the presence of fusidic acid takes its own fusidic-acid-on-colistin EC50 interaction (-0.986) in place of the reference-strain value (-0.826). Mutually exclusive with STRAIN_AB12.",
      source_name = "strain"
    ),
    STRAIN_AB12 = list(
      description = "1 = the experiment used A. baumannii clinical isolate AB12 (SE clone; colistin MIC 1 mg/L, fusidic acid MIC 64 mg/L); 0 = any other strain.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 -- strains AB14, ATCC 19606, AB16 and AB315 (reference parameter values); AB205 carries its own indicator (STRAIN_AB205)",
      notes = "AB12 lowers beta for the colistin + fusidic acid combination by (1 - 0.0523) (Table 2, 'proportional decrease in theta_beta_comb with AB12'). Mutually exclusive with STRAIN_AB205.",
      source_name = "strain"
    ),
    INOC_LOG10CFU = list(
      description = "Observed starting inoculum of the time-kill experiment, used as the initial condition of the bacterial count.",
      units = "log10 CFU/mL",
      type = "continuous",
      reference_category = "n/a -- standard inoculum about 6 log10 CFU/mL; high-inoculum experiments with AB14 used 7-8 log10 CFU/mL (Methods, Time-kill assays)",
      notes = "Table 2 carries no inoculum parameter: the model starts each experiment from its observed baseline (Methods, Eq. 1 and the sentence that follows it). Eight high-inoculum experiments showed no inoculum effect on model fit (Results; Figures S5-S8).",
      source_name = "baseline cfu"
    )
  )

  population <- list(
    species = "in vitro (Acinetobacter baumannii; ATCC 19606 and clinical isolates AB12, AB14, AB16, AB205, AB315)",
    n_subjects = 312L,
    n_studies = 1L,
    disease_state = "Six A. baumannii strains: the antibiotic-susceptible type strain ATCC 19606 and five MDR clinical isolates (Public Health England) representing OXA-23 UK clone I (AB14, AB315), OXA-23 UK clone II (AB16), SE clone (AB12) and a colistin-resistant strain (AB205). Broth microdilution colistin MICs 0.25 (ATCC 19606), 0.5 (AB14, AB16), 1 (AB12, AB315) and > 256 mg/L (AB205); fusidic acid MICs 32 (AB205), 64 (AB12), 128 (AB14, ATCC 19606), 256 (AB315) and 512 mg/L (AB16).",
    model_system = "Static time-kill assays in 10 mL Iso-Sensitest broth at 37 degrees C, starting inoculum about 10^6 CFU/mL (10^7-10^8 CFU/mL in a separate AB14 inoculum experiment), sampled over 24 h; colony counts by plating 100 uL aliquots, lower limit 10 CFU/mL. 312 individual time-kill curves covering 114 unique pathogen-antimicrobial combinations (27 colistin only, 33 fusidic acid only, 54 colistin plus fusidic acid).",
    dose_range = "Static concentrations: colistin 0.25-8192 mg/L, fusidic acid 1-8192 mg/L, alone and in combination",
    n_observations = "312 time-kill curves (Results, Time-kill experiments)",
    notes = paste(
      "Estimated in NONMEM 7.3 with ADVAN9 and the Laplacian method; counts below 10 CFU/mL handled with the M3 method; CFU modelled as log10.",
      "The model was built stepwise with parameters fixed between steps: drug-free growth, then colistin and fusidic acid alone on AB14, then the combination on AB14, then the other strains added sequentially as categorical covariates.",
      "There is no between-experiment variability; the only random effect is an additive residual error on log10 CFU/mL with a separate magnitude for growth-control, colistin-only, fusidic-acid-only and combination experiments.",
      "Standard errors from a 1000-sample non-parametric bootstrap."
    )
  )

  ini({
    # --- Drug-free growth (Table 2; Eqs 3-4) ---
    kgrow <- 1.87; label("Net bacterial growth rate constant, theta_net (1/h)") # Table 2: theta_net 1.87 1/h (RSE 4.5%)
    log10bmax <- 9.41; label("Maximum carrying capacity, theta_cfumax (log10 CFU/mL)") # Table 2: theta_cfumax 9.41 (RSE 0.600%), printed in CFU/mL but a log10 value (Figure S1 plateaus near 10^9.4)
    tlag50 <- 0.352; label("Time of half-maximal lag, theta_lag (h)") # Table 2: theta_lag 0.352 h (RSE 14.4%)

    # --- Colistin effect (Table 2; Eq 5) ---
    # Table 2 prints the EC50 units as ng/mL; the model only reproduces the
    # Figure 2 kill-and-regrowth and the AB205 colistin MIC when the EC50 is
    # read in mg/L (see the vignette, 'EC50 units').
    ec50_col <- 9.86; label("Colistin concentration at half-maximal kill rate, theta_EC50col (mg/L)") # Table 2: theta_EC50col 9.86 (RSE 23.5%), printed as ng/mL
    emax_col <- 39.5; label("Colistin maximum kill rate, theta_Emaxcol (1/h)") # Table 2: theta_Emaxcol 39.5 1/h (RSE 10.6%)
    hill_col <- 0.855; label("Colistin sigmoidicity, theta_gamma_col (unitless)") # Table 2: theta_gamma_col 0.855 (RSE 4.60%)
    tau_col <- 0.177; label("Colistin rate of loss of effect per unit concentration, theta_tau_col (L/mg/h)") # Table 2: theta_tau_col 0.177 (RSE 38.6%)
    beta_col <- 0.895; label("Maximum fractional loss of colistin effect, theta_beta_col (fraction)") # Table 2: theta_beta_col 0.895 (RSE 1.10%)
    e_ab205_ec50_col <- 65.2; label("Proportional increase in the colistin EC50 for strain AB205 (fraction)") # Table 2: 'proportional increase in theta_EC50col with AB205' 65.2 (RSE 41.6%)

    # --- Fusidic acid effect (Table 2; Eq 5) ---
    ec50_fus <- 310; label("Fusidic acid concentration at half-maximal kill rate, theta_EC50fus (mg/L)") # Table 2: theta_EC50fus 310 (RSE 43.9%), printed as ng/mL
    emax_fus <- 23.2; label("Fusidic acid maximum kill rate, theta_Emaxfus (1/h)") # Table 2: theta_Emaxfus 23.2 1/h (RSE 6.30%)
    hill_fus <- 0.776; label("Fusidic acid sigmoidicity, theta_gamma_fus (unitless)") # Table 2: theta_gamma_fus 0.776 (RSE 8.30%)
    tau_fus <- 0.0102; label("Fusidic acid rate of loss of effect per unit concentration, theta_tau_fus (L/mg/h)") # Table 2: theta_tau_fus 0.0102 (RSE 37.2%)
    beta_fus <- fixed(0.895); label("Maximum fractional loss of fusidic acid effect, theta_beta_fus (fraction)") # Table 2: theta_beta_fus 0.895 (fixed); Methods: fixed to the colistin estimate
    e_ab205_tau_fus <- 2.12; label("Proportional increase in tau for fusidic acid for strain AB205 (fraction)") # Table 2: 'proportional increase in theta_tau_fus with AB205' 2.12 (RSE 48.6%)
    e_ab205_beta_fus <- -0.124; label("Proportional change in beta for fusidic acid for strain AB205 (fraction)") # Table 2: 'proportional decrease in theta_beta_fus with AB205' -0.124 (RSE 22.3%)

    # --- Colistin + fusidic acid interaction (Table 2; Eqs 7-8) ---
    e_fus_ec50_col <- -0.826; label("Proportional change in the colistin EC50 when fusidic acid is present, reference strains (fraction)") # Table 2: 'proportional decrease in theta_EC50col with fusidic acid' -0.826 (RSE 1.20%)
    e_fus_ec50_col_ab205 <- -0.986; label("Proportional change in the colistin EC50 when fusidic acid is present, strain AB205 (fraction)") # Table 2: 'proportional decrease in theta_EC50col with fusidic acid for AB205' -0.986 (RSE 0.200%)
    e_col_emax_fus <- 0.579; label("Proportional increase in the fusidic acid Emax when colistin is present (fraction)") # Table 2: 'proportional increase in theta_Emaxfus with colistin' 0.579 (RSE 52.2%)
    tau_comb <- 1.44; label("Rate of loss of combined effect per unit concentration product, theta_tau_comb (L^2/mg^2/h)") # Table 2: theta_tau_comb 1.44 (RSE 15.9%)
    beta_comb <- fixed(0.895); label("Maximum fractional loss of combined effect, theta_beta_comb (fraction)") # Table 2: theta_beta_comb 0.895 (fixed); Methods: fixed to the colistin estimate
    e_ab12_beta_comb <- -0.0523; label("Proportional change in beta for the combination for strain AB12 (fraction)") # Table 2: 'proportional decrease in theta_beta_comb with AB12' -0.0523 (RSE 22.6%)

    # --- Residual error on log10 CFU/mL (Table 2), one value per experiment type ---
    # Table 2 lists these under 'Fixed effect', i.e. as THETAs, read here as
    # standard deviations; the growth-control VPC band supports the SD reading
    # (see the vignette).
    addSd_growth <- 0.129; label("Additive residual SD, growth-control experiments (log10 CFU/mL)") # Table 2: 'Residual variability growth model' 0.129 (RSE 15.9%)
    addSd_col <- 2.8; label("Additive residual SD, colistin-only experiments (log10 CFU/mL)") # Table 2: 'Residual variability col' 2.8 (RSE 47.5%)
    addSd_fus <- 0.698; label("Additive residual SD, fusidic-acid-only experiments (log10 CFU/mL)") # Table 2: 'Residual variability fus' 0.698 (RSE 9.30%)
    addSd_comb <- 5.14; label("Additive residual SD, combination experiments (log10 CFU/mL)") # Table 2: 'Residual variability comb' 5.14 (RSE 8.8%)
  })

  model({
    bmax <- 10^log10bmax

    # Presence flags for the categorical interaction terms (Eq 7) and for the
    # choice of beta and residual SD.
    on_col <- 0
    if (CONC_COLISTIN_MGL > 0) on_col <- 1
    on_fus <- 0
    if (CONC_FUS_MGL > 0) on_fus <- 1
    on_comb <- on_col * on_fus

    # Growth rate with lag (Eq 4): -theta_net at t = 0, rising steeply through
    # zero at t = theta_lag to +theta_net; the exponent 20 is printed in Eq 4.
    kg <- -kgrow + 2 * kgrow * t^20 / (tlag50^20 + t^20)

    # Covariate-adjusted parameters, each as theta * (1 + theta_cov) (Eq 7).
    # For AB205 the fusidic-acid-on-colistin EC50 term (-0.986) replaces the
    # reference-strain term (-0.826); see the vignette for this reading.
    e_fus_ec50 <- e_fus_ec50_col * (1 - STRAIN_AB205) + e_fus_ec50_col_ab205 * STRAIN_AB205
    ec50_col_i <- ec50_col * (1 + e_ab205_ec50_col * STRAIN_AB205) * (1 + e_fus_ec50 * on_fus)
    emax_fus_i <- emax_fus * (1 + e_col_emax_fus * on_col)
    tau_fus_i <- tau_fus * (1 + e_ab205_tau_fus * STRAIN_AB205)
    beta_fus_i <- beta_fus * (1 + e_ab205_beta_fus * STRAIN_AB205)
    beta_comb_i <- beta_comb * (1 + e_ab12_beta_comb * STRAIN_AB12)

    # beta: colistin value alone, fusidic acid value alone, combination value
    # when both are present (Methods, text after Eqs 6 and 8).
    beta_i <- beta_col * on_col * (1 - on_fus) + beta_fus_i * on_fus * (1 - on_col) + beta_comb_i * on_comb

    # Sigmoidal Emax kill per drug (Eq 5), additive in combination.
    eff_col <- emax_col * CONC_COLISTIN_MGL^hill_col / (ec50_col_i^hill_col + CONC_COLISTIN_MGL^hill_col)
    eff_fus <- emax_fus_i * CONC_FUS_MGL^hill_fus / (ec50_fus^hill_fus + CONC_FUS_MGL^hill_fus)

    # Time- and concentration-dependent loss of effect (Eq 8; reduces to the
    # Eq 5 factor for either drug alone).
    tde <- 1 - beta_i * (1 - exp(-t * (tau_col * CONC_COLISTIN_MGL + tau_fus_i * CONC_FUS_MGL + tau_comb * CONC_COLISTIN_MGL * CONC_FUS_MGL)))
    kill <- (eff_col + eff_fus) * tde

    # Logistic growth and kill on the linear CFU scale (Eq 6).
    d/dt(cfu) <- (kg - kill) * cfu * (1 - cfu / bmax)
    cfu(0) <- 10^INOC_LOG10CFU

    # Residual SD by experiment type (Table 2).
    sd_i <- addSd_growth * (1 - on_col) * (1 - on_fus) + addSd_col * on_col * (1 - on_fus) + addSd_fus * on_fus * (1 - on_col) + addSd_comb * on_comb

    # log10 CFU/mL; the floor only guards log10() against a numerically
    # non-positive state after very fast kill.
    log_cfu <- log10(max(cfu, 1e-30))
    log_cfu ~ add(sd_i)
  })
}
