Zhao_2025_paracetamol <- function() {
  description <- "Parent-and-metabolites population PK model for oral paracetamol and its glucuronide, sulphate and combined oxidative (cysteine + mercapturate) metabolites in children and adults with spinal muscular atrophy (SMA) and healthy controls (Zhao 2025). One compartment per compound with first-order absorption through a depot, a fixed absorption lag time, and body-weight allometric scaling (exponent fixed at 0.75 on every clearance and 1 on every volume). Paracetamol leaves the central compartment by four parallel first-order routes: glucuronide formation, sulphate formation, oxidative-metabolite formation, and a leftover clearance covering unchanged drug plus any unaccounted route. Each metabolite has its own one-compartment plasma pool with an apparent volume fixed at 18% of the paracetamol volume and its own first-order elimination clearance. SMA disease status raises the paracetamol volume of distribution (and, through the fixed 0.18 ratio, every metabolite volume) by 58%; plasma myoglobin scales the paracetamol leftover clearance with a negative power exponent, and plasma total bilirubin scales sulphate-formation clearance positively and oxidative-metabolite elimination clearance negatively."
  reference <- paste(
    "Zhao Q, Naume MM, de Winter BCM, Krag T, Haslund-Krog SS, Revsbech KL,",
    "Vissing J, Holst H, Moller MH, Hornsyld TM, Duno M, Hoei-Hansen CE,",
    "Born AP, Jensen PB, Orngreen MC (2025).",
    "Paracetamol and its metabolites in children and adults with spinal",
    "muscular atrophy - a population pharmacokinetic model.",
    "Br J Clin Pharmacol 91(7):2045-2056.",
    "doi:10.1002/bcp.70028.",
    sep = " "
  )
  vignette <- "Zhao_2025_paracetamol"

  # Molar parameterisation, as in vanRongen_2016_acetaminophen.R. Paracetamol
  # transfers to each conjugate 1:1 on a molar basis, so a single set of
  # formation / elimination clearances describes the parent and all three
  # metabolite pools only when the amounts are moles. The scale is not stated
  # in Zhao 2025, but it is determined by the paper's own simulation output:
  # supplementary Tables S2 / S3 report steady-state AUC(0-24 h) for all four
  # analytes, and the AUC ratio of a metabolite to the parent is
  # volume-independent (it equals CL_formation / CL_elimination in whatever
  # amount unit the compartments carry). Reproducing the model with molar
  # amounts and converting each output to that species' own mass units gives
  # metabolite:paracetamol AUC(0-24 h) ratios of 2.12 / 0.615 / 0.077, against
  # the published 2.07 / 0.609 / 0.076 -- agreement to within 2.5%. A
  # mass-equivalent (paracetamol-mass) parameterisation gives 0.98 / 0.40 /
  # 0.043 on the same window, i.e. 34-53% low. Molecular weights for the mass
  # conversion (paracetamol 151.16, glucuronide 327.29, sulphate 231.23,
  # cysteine conjugate 270.30 g/mol) are tabulated in vanRongen 2016 Methods
  # and are applied in the validation vignette, not in model(); see the
  # vignette's Errata section.
  units <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot          = list(analyte = "paracetamol", units = "umol", specimen = "administration site", verified = TRUE),
    central        = list(analyte = "paracetamol", units = "umol", specimen = "plasma", verified = TRUE),
    central_gluc   = list(analyte = "paracetamol-glucuronide", units = "umol", specimen = "plasma", verified = TRUE),
    central_sulf   = list(analyte = "paracetamol-sulphate", units = "umol", specimen = "plasma", verified = TRUE),
    central_cysmer = list(analyte = "paracetamol oxidative metabolites (cysteine + mercapturate)", units = "umol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric size descriptor with reference weight 70 kg. Exponent fixed at 0.75 on every clearance and at 1 on every volume (Zhao 2025 Results section 3.2 and Table 2 caption: 'the exponent of body weight allometric scaling was fixed at 0.75 for clearance parameters and 1 for volume of distribution'). Estimating the exponent did not improve the fit, and fat-free mass performed worse than total body weight (section 3.2).",
      source_name        = "body weight"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-control indicator (1 = healthy control, 0 = patient with spinal muscular atrophy).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient with spinal muscular atrophy)",
      notes              = "Zhao 2025's source column is the reverse-coded 'disease' indicator (1 = SMA), entered on the paracetamol volume of distribution as the proportional form of supplementary Eq. 1 (Pi = theta_pop * theta_cov * exp(eta), theta_cov = 1 for the reference group). Table 2 gives theta_cov = 1.58 for SMA (bootstrap median 1.62, 95% CI 1.20-2.30) with the healthy controls as the paper's reference, and Results section 3.3 reports dOFV = 11.5. This file re-expresses the indicator on the canonical DIS_HEALTHY orientation (reference category 0 = patient), following Chen_2023_nemonoxacin.R and Cleary_2023_risdiplam.R -- the latter is the same SMA-vs-healthy contrast in a risdiplam popPK model. The structural typical is therefore shifted to the SMA state, lvc = log(63.5 * 1.58), and the effect parameter 1 / 1.58 restores the printed healthy-control typical of 63.5 L/70 kg at DIS_HEALTHY = 1. Because every metabolite volume is fixed at 0.18 * V_paracetamol, the same factor propagates to the glucuronide, sulphate and oxidative-metabolite volumes -- exactly as the abstract states ('SMA disease resulted in a 58% increase in the volume of distribution for paracetamol and its metabolites') and as Table 3's post-hoc medians confirm (V_glu HC 13.21 = 0.18 * 73.40; SMA 17.08 = 0.18 * 94.90).",
      source_name        = "disease"
    ),
    MYO = list(
      description        = "Plasma myoglobin concentration.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-model covariate on the paracetamol leftover clearance, normalised to the pooled population median of 25 ng/mL (Zhao 2025 Table 2 row label 'Myoglobin (mean 25 ng/mL)'; supplementary Eq. 2 gives the power form Pi = theta_pop * (COV_i / COV_median)^theta_cov * exp(eta)). Exponent -1.10 (bootstrap median -1.14, 95% CI -3.07 to -0.49); Results section 3.3 records dOFV = 18.5 and states the correlation is negative. Baseline medians: healthy controls 34 ng/mL, SMA 17 ng/mL (Table 1). Myoglobin is a skeletal-muscle protein used clinically as a muscle-damage marker; the authors interpret the effect as an indirect glutathione-availability / muscle-mass signal (Discussion).",
      source_name        = "myoglobin"
    ),
    TBILI = list(
      description        = "Plasma total bilirubin concentration.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-model covariate (supplementary Eq. 2) normalised to the pooled population median of 6 umol/L (Zhao 2025 Table 2 row label 'Bilirubin (mean 6 umol/L)'). Two separate effects: +0.18 on the sulphate-formation clearance (bootstrap median 0.16, 95% CI -0.02 to 0.39; positive correlation) and -0.177 on the oxidative-metabolite elimination clearance (bootstrap median -0.176, 95% CI -0.45 to -0.03; negative correlation). Both entered with dOFV = 7.3 (Results section 3.3). Baseline medians: healthy controls 7 umol/L, SMA 4 umol/L (Table 1); values were within the normal range in both groups.",
      source_name        = "bilirubin"
    )
  )

  # Screened in the covariate analysis (supplementary Methods, Covariate
  # Analysis) but not retained in the final model. Documented here so the
  # provenance of the covariate screen is preserved without declaring
  # covariates that model() never references. Urea and glomerular filtration
  # rate were also screened; they have no unambiguous canonical column in
  # inst/references/covariate-columns.md and are recorded in population$notes
  # instead.
  covariatesDataExcluded <- list(
    SEXF  = list(description = "Sex (female indicator)", units = "(binary)", type = "binary", notes = "Screened as a categorical covariate; not significant."),
    AGE   = list(description = "Age", units = "years", type = "continuous", notes = "Screened; not retained. The authors note myoglobin is itself age- and disease-dependent, so age may be a confounder of the retained myoglobin effect (Discussion)."),
    BMI   = list(description = "Body mass index", units = "kg/m^2", type = "continuous", notes = "Screened as an alternative size descriptor; 'BMI and FFM failed to become significant covariates' (Results section 3.3)."),
    LBM   = list(description = "Fat-free mass", units = "kg", type = "continuous", notes = "Tested both as an allometric size descriptor and as a covariate; 'using fat free mass (FFM) as an allometric scaling factor resulted in a less optimal model performance compared to using body weight' (Results section 3.2)."),
    ALT   = list(description = "Alanine aminotransferase", units = "U/L", type = "continuous", notes = "Screened as a hepatic-function marker; not retained."),
    AST   = list(description = "Aspartate aminotransferase", units = "U/L", type = "continuous", notes = "Screened as a hepatic-function marker; not retained."),
    ALP   = list(description = "Alkaline phosphatase", units = "U/L", type = "continuous", notes = "Screened; not retained."),
    LDH   = list(description = "Lactate dehydrogenase", units = "U/L", type = "continuous", notes = "Screened; not retained."),
    CPK   = list(description = "Creatine kinase", units = "U/L", type = "continuous", notes = "Screened as a muscle-damage marker alongside myoglobin; not retained."),
    CREAT = list(description = "Serum creatinine", units = "umol/L", type = "continuous", notes = "Screened; not retained. Median 79 umol/L in healthy controls vs 9 umol/L in SMA patients (Table 1)."),
    POT   = list(description = "Potassium", units = "mmol/L", type = "continuous", notes = "Screened; not retained."),
    SOD   = list(description = "Sodium", units = "mmol/L", type = "continuous", notes = "Screened; not retained.")
  )

  population <- list(
    species        = "human",
    n_subjects     = 23L,
    n_studies      = 1L,
    age_range      = "SMA patients 6-37 years (median 17); healthy controls 20-36 years (median 25)",
    age_median     = "17 years (SMA) / 25 years (healthy controls)",
    weight_range   = "SMA patients 22-57 kg (median 30.5); healthy controls 51-103 kg (median 78.0)",
    weight_median  = "30.5 kg (SMA) / 78.0 kg (healthy controls)",
    sex_female_pct = 43.5,
    race_ethnicity = NA_character_,
    disease_state  = "Six adults with spinal muscular atrophy, six children with spinal muscular atrophy, and 11 healthy controls (Zhao 2025 Methods section 2.1 and Results section 3.1). SMA is an inherited neuromuscular disorder caused by SMN1 mutations, giving progressive motor-neuron degeneration, muscle weakness and markedly reduced skeletal muscle mass. Baseline myoglobin and total bilirubin were both significantly lower in SMA patients than in healthy controls (Table 1).",
    dose_range     = "Oral paracetamol 15 mg/kg every 6 h for 3 days, capped at a maximum single dose of 1 g (Methods section 2.1). Blood samples hourly for 6-8 h on Days 1 and 3 after a pre-treatment baseline sample; 294 plasma samples per analyte entered the model.",
    regions        = "Denmark (Copenhagen University Hospital Rigshospitalet and Bispebjerg Hospital)",
    notes          = "Sex split 7 male / 5 female among SMA patients and 6 male / 5 female among healthy controls, i.e. 10 of 23 (43.5%) female. Baseline (pre-dose) samples were below the lower limit of quantification in most subjects and were fixed to their observed values during base-model building, so the model carries no baseline variability (Results section 3.2). Three cysteine (1.1%) and four mercapturic-acid (1.5%) on-treatment samples below the LLOQ were censored. Screened-but-not-retained covariates additionally included urea and glomerular filtration rate (supplementary Methods, Covariate Analysis); both are omitted from covariatesDataExcluded because neither has an unambiguous canonical column in the register (GFR was reported only as the interval '87 - >90 mL/min' for both groups). Missing covariate values were imputed by last observation carried forward. EudraCT 2018-002295-40; ethics approval H-18032928."
  )

  ini({
    # Final-model parameter estimates from Zhao 2025 Table 2, "Final model"
    # column (page 2049). The bootstrap median and 95% CI columns of the same
    # table are quoted alongside each value. NONMEM 7.4 FOCE-I, ADVAN5
    # (supplementary Methods, Base Model Development). All clearances are
    # reported per 70 kg and all volumes at 70 kg, so the allometric terms in
    # model() are written as (WT / 70)^exponent.

    # Absorption ----------------------------------------------------------
    ltlag <- fixed(log(0.153)) ; label("Absorption lag time (h)")                 # Table 2 t_lag = 0.153 h, flagged "(FIX)" in the row label and in the table footnote ("lag time was fixed to 0.153 in the model building"); bootstrap 95% CI 0.153-0.153 confirms it was not estimated
    lka   <- log(2.84)         ; label("First-order absorption rate constant (1/h)") # Table 2 k_a final 2.84 1/h (bootstrap median 2.44, 95% CI 1.43-5.64)

    # Paracetamol disposition ---------------------------------------------
    lvc        <- log(63.5 * 1.58) ; label("Apparent paracetamol volume of distribution V_pcm/F at 70 kg, SMA patient (L/70 kg)") # Table 2 V_pcm/F final 63.5 L/70 kg (bootstrap median 61.67, 95% CI 43.44-81.07), which is the healthy-control typical; multiplied by the Table 2 disease factor 1.58 to shift the structural typical onto the canonical DIS_HEALTHY = 0 (patient) reference state
    lcl        <- log(8.16) ; label("Apparent paracetamol leftover clearance CL_p/F at 70 kg, myoglobin 25 ng/mL (L/h/70 kg)") # Table 2 CL_p/F final 8.16 L/h/70 kg (bootstrap median 6.90, 95% CI 2.15-11.57); "leftover" = unchanged drug plus any route not captured by the three named pathways (Discussion)
    lcl_gluc   <- log(6.37) ; label("Apparent glucuronide formation clearance CL_pg/F at 70 kg (L/h/70 kg)")                   # Table 2 CL_pg/F final 6.37 L/h/70 kg (bootstrap median 6.61, 95% CI 4.86-8.46)
    lcl_sulf   <- log(8.43) ; label("Apparent sulphate formation clearance CL_ps/F at 70 kg, bilirubin 6 umol/L (L/h/70 kg)")  # Table 2 CL_ps/F final 8.43 L/h/70 kg (bootstrap median 9.00, 95% CI 7.25-11.80)
    lcl_cysmer <- log(0.20) ; label("Apparent oxidative-metabolite formation clearance CL_pox/F at 70 kg (L/h/70 kg)")         # Table 2 CL_pox/F final 0.20 L/h/70 kg (bootstrap median 0.21, 95% CI 0.170-0.26)

    # Metabolite elimination ----------------------------------------------
    lcle_gluc   <- log(5.69) ; label("Apparent glucuronide elimination clearance CL_gluc/F at 70 kg (L/h/70 kg)")                       # Table 2 CL_gluc/F final 5.69 L/h/70 kg (bootstrap median 5.73, 95% CI 4.78-7.17)
    lcle_sulf   <- log(20.4) ; label("Apparent sulphate elimination clearance CL_sulf/F at 70 kg (L/h/70 kg)")                          # Table 2 CL_sulf/F final 20.4 L/h/70 kg (bootstrap median 20.95, 95% CI 18.61-27.00)
    lcle_cysmer <- log(3.72) ; label("Apparent oxidative-metabolite elimination clearance CL_ox/F at 70 kg, bilirubin 6 umol/L (L/h/70 kg)") # Table 2 CL_ox/F final 3.72 L/h/70 kg (bootstrap median 3.65, 95% CI 2.94-4.80)

    # Allometric exponents, both fixed by the authors -----------------------
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of (WT / 70) on every clearance (unitless)") # Table 2 caption and Results section 3.2: "The allometric theory-based exponent (EXP) was fixed at 0.75 for clearance (CL) parameters and 1 for V d"
    e_wt_vc <- fixed(1)    ; label("Allometric exponent of (WT / 70) on every volume (unitless)")    # Table 2 caption and Results section 3.2, same sentence

    # Covariate effects ----------------------------------------------------
    e_dis_healthy_vc   <- 1 / 1.58 ; label("Multiplicative factor on V_pcm/F for healthy controls (DIS_HEALTHY = 1; unitless)") # Table 2 "Covariate effect on V_pcm/F -- Disease" final 1.58 (bootstrap median 1.62, 95% CI 1.20-2.30) with the healthy controls as the paper's reference; abstract quotes the same effect as "58% (95% CI 20%-130%) increase". Written as 1 / 1.58 = 0.633 so that the canonical DIS_HEALTHY = 1 state restores the printed healthy-control typical V_pcm/F of 63.5 L/70 kg
    e_myo_cl           <- -1.10  ; label("Power exponent of (MYO / 25) on CL_p/F (unitless)")                             # Table 2 "Covariate effect on CL_p/F -- Myoglobin (mean 25 ng/mL)" final 1.10, bootstrap median 1.14, bootstrap 95% CI printed as "3.07-0.49". The minus signs are lost in the typeset table (the CI is printed in descending order, which is only possible for a negative interval), and Results section 3.3 states the correlation is negative; the sign is confirmed numerically by supplementary Table S3, where the paracetamol half-life rises monotonically from 2.20 h at MYO = 15 ng/mL to 3.69 h at MYO = 70 ng/mL
    e_tbili_cl_sulf    <- 0.18   ; label("Power exponent of (TBILI / 6) on CL_ps/F (unitless)")                           # Table 2 "Covariate effect on CL_ps/F -- Bilirubin (mean 6 umol/L)" final 0.18 (bootstrap median 0.16, 95% CI -0.02 to 0.39); positive correlation per Results section 3.3
    e_tbili_cle_cysmer <- -0.177 ; label("Power exponent of (TBILI / 6) on CL_ox/F (unitless)")                           # Table 2 "Covariate effect on CL_ox/F -- Bilirubin (mean 6 umol/L)" final 0.177, bootstrap median 0.176, bootstrap 95% CI printed as "0.45-0.03" (descending, hence negative). Results section 3.3 states the correlation is negative and the Discussion adds "a decrease in bilirubin could ... increase CL_ox/F"; confirmed numerically by supplementary Table S3, where the oxidative-metabolite half-life rises from 3.18 h at 3 umol/L to 4.19 h at 24 umol/L

    # Between-subject variability. Zhao 2025 Table 2 reports a "BSV (%)"
    # block; the supplementary Methods state BSV was added "using an
    # exponential model", i.e. P_i = P_pop * exp(eta_i). The percentages are
    # the log-scale standard deviation omega * 100, NOT an apparent CV: the
    # Discussion quotes the CL_p/F BSV as 105%, and
    # sqrt(exp(0.862^2) - 1) = 1.050 reproduces that exactly from the Table 2
    # value of 86.2%. nlmixr2 eta values are variances, hence (BSV/100)^2.
    etalka         ~ 2.1609  # Table 2 K_a BSV 147%   -> omega = 1.47,  variance 1.47^2
    etalvc         ~ 0.0801  # Table 2 V_pcm/F BSV 28.3% -> omega = 0.283, variance 0.283^2
    etalcl         ~ 0.7430  # Table 2 CL_p/F BSV 86.2%  -> omega = 0.862, variance 0.862^2 (the 105% quoted in the Discussion is the corresponding CV)
    etalcl_gluc    ~ 0.0906  # Table 2 CL_pg/F BSV 30.1% -> omega = 0.301, variance 0.301^2
    etalcl_sulf    ~ 0.0650  # Table 2 CL_ps/F BSV 25.5% -> omega = 0.255, variance 0.255^2
    etalcle_cysmer ~ 0.0697  # Table 2 CL_ox/F BSV 26.4% -> omega = 0.264, variance 0.264^2

    # Residual error. Table 2's "Residual variability -- Proportional" block
    # gives one proportional term per analyte on the standard-deviation
    # (fraction) scale, matching the standard-deviation scale used for the BSV
    # block in the same table; Results section 3.4 quotes the same numbers as
    # "proportional error ... decreased from 0.283 to 0.280 ...". The
    # supplementary goodness-of-fit panel for paracetamol (observed vs
    # individual-predicted, alongside IWRES vs individual-predicted) shows
    # relative residuals whose spread is consistent with a ~0.28 SD and not
    # with a 0.53 SD (which is what reading 0.28 as a variance would imply).
    propSd        <- 0.28 ; label("Paracetamol proportional residual SD (fraction)")            # Table 2 final proportional residual, paracetamol = 0.28 (bootstrap 95% CI 0.23-0.33)
    propSd_gluc   <- 0.18 ; label("Glucuronide proportional residual SD (fraction)")            # Table 2 final proportional residual, paracetamol-glucuronide = 0.18 (bootstrap 95% CI 0.16-0.20)
    propSd_sulf   <- 0.22 ; label("Sulphate proportional residual SD (fraction)")               # Table 2 final proportional residual, paracetamol-sulfate = 0.22 (bootstrap 95% CI 0.18-0.24)
    propSd_cysmer <- 0.17 ; label("Oxidative-metabolite proportional residual SD (fraction)")   # Table 2 final proportional residual, paracetamol-oxidative pathway metabolites = 0.17 (bootstrap 95% CI 0.13-0.20)
  })

  model({
    # Absorption. The lag time is a fixed typical value with no IIV; the
    # absorption rate constant carries the largest BSV in the model (147%).
    tlag <- exp(ltlag)
    ka   <- exp(lka + etalka)

    # Paracetamol volume. Allometric exponent 1 on (WT / 70), then the
    # categorical disease effect written in the proportional form of
    # supplementary Eq. 1. lvc carries the SMA typical value, so the
    # multiplier is 1 for DIS_HEALTHY = 0 (the canonical patient reference)
    # and 1 / 1.58 for DIS_HEALTHY = 1, recovering the paper's printed
    # healthy-control V_pcm/F of 63.5 L/70 kg.
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc *
      (1 + (e_dis_healthy_vc - 1) * DIS_HEALTHY)

    # Metabolite volumes. Zhao 2025 Methods section 2.2: "The volume of
    # distribution (V d ) for all paracetamol metabolites cannot be
    # identified, so we fixed them to 18% of the central V d of paracetamol
    # in plasma based on literature." Table 2 confirms the arithmetic:
    # 0.18 * 63.5 = 11.43 L/70 kg, the reported V_glu/F = V_sulf/F = V_ox/F.
    # The SMA factor therefore propagates into every metabolite volume.
    vmet <- 0.18 * vc

    # Clearances. Every clearance carries the fixed 0.75 allometric exponent
    # on (WT / 70). The three retained continuous covariate effects use the
    # median-normalised power model of supplementary Eq. 2.
    cl         <- exp(lcl         + etalcl)         * (WT / 70)^e_wt_cl * (MYO / 25)^e_myo_cl
    cl_gluc    <- exp(lcl_gluc    + etalcl_gluc)    * (WT / 70)^e_wt_cl
    cl_sulf    <- exp(lcl_sulf    + etalcl_sulf)    * (WT / 70)^e_wt_cl * (TBILI / 6)^e_tbili_cl_sulf
    cl_cysmer  <- exp(lcl_cysmer)                   * (WT / 70)^e_wt_cl
    cle_gluc   <- exp(lcle_gluc)                    * (WT / 70)^e_wt_cl
    cle_sulf   <- exp(lcle_sulf)                    * (WT / 70)^e_wt_cl
    cle_cysmer <- exp(lcle_cysmer + etalcle_cysmer) * (WT / 70)^e_wt_cl * (TBILI / 6)^e_tbili_cle_cysmer

    # Micro-constants. Formation rates are partial paracetamol clearances
    # divided by the paracetamol volume; metabolite elimination rates are the
    # metabolite clearance divided by the (shared) metabolite volume.
    k_form_gluc   <- cl_gluc   / vc
    k_form_sulf   <- cl_sulf   / vc
    k_form_cysmer <- cl_cysmer / vc
    k_left        <- cl        / vc
    kel_gluc      <- cle_gluc   / vmet
    kel_sulf      <- cle_sulf   / vmet
    kel_cysmer    <- cle_cysmer / vmet

    # ODEs reproducing Figure 1 of Zhao 2025 (schematic of the structural
    # model): oral depot -> V_pcm with ka and a lag time; V_pcm drains to
    # V_gluc, V_sulf and V_ox through the three formation clearances and to
    # urine/bile through the leftover clearance CL_p; each metabolite pool
    # drains to urine/bile through its own elimination clearance.
    d/dt(depot)          <- -ka * depot
    alag(depot)          <-  tlag
    d/dt(central)        <-  ka * depot -
      (k_form_gluc + k_form_sulf + k_form_cysmer + k_left) * central
    d/dt(central_gluc)   <-  k_form_gluc   * central - kel_gluc   * central_gluc
    d/dt(central_sulf)   <-  k_form_sulf   * central - kel_sulf   * central_sulf
    d/dt(central_cysmer) <-  k_form_cysmer * central - kel_cysmer * central_cysmer

    # Plasma concentrations. Amounts are micromoles and volumes are litres,
    # so each ratio is umol/L. Users dosing in mg should convert with
    # dose_umol = dose_mg / 0.15116 (paracetamol MW 151.16 g/mol); to compare
    # a metabolite against a published mass concentration, multiply the molar
    # concentration by that metabolite's own molecular weight.
    Cc        <- central        / vc
    Cc_gluc   <- central_gluc   / vmet
    Cc_sulf   <- central_sulf   / vmet
    Cc_cysmer <- central_cysmer / vmet

    Cc        ~ prop(propSd)
    Cc_gluc   ~ prop(propSd_gluc)
    Cc_sulf   ~ prop(propSd_sulf)
    Cc_cysmer ~ prop(propSd_cysmer)
  })
}
