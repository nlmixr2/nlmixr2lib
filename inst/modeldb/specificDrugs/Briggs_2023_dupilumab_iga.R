Briggs_2023_dupilumab_iga <- function() {
  description <- "Exposure-response model for the Investigator's Global Assessment (IGA, ordinal 0-4) under dupilumab in children (6-11 years), adolescents and adults with moderate-to-severe atopic dermatitis (Briggs 2023). PD-only: the time-varying serum functional dupilumab concentration is an input column (CP_DUPILUMAB_MGL, mg/L), as in the paper's sequential PK-then-PD fit to individual predictions from the age-group dupilumab popPK models (see modellib('Kovalenko_2020_dupilumab_covariate'), modellib('Kovalenko_2021_dupilumab_adolescent_covariate') and modellib('Kovalenko_2021_dupilumab_children_covariate')). Cumulative-probit latent-variable model: probit P(IGA <= m) = Q_m + Pmax * t / (ET50 + t) + DSLP * R, where Q_3 is the baseline cutpoint (additive IIV, shifted for moderate baseline disease), the lower cutpoints are spaced by estimated gaps, and R is an indirect-response drug signal that falls from 0 towards -C / (IC50 + C). Covariates: age, topical corticosteroids, prior systemic immunotherapy and moderate disease on Pmax; body weight, Asian race, eosinophil count and TARC on the drug slope. Outputs the per-category probabilities p0-p4 and the IGA 0/1 responder probability prob_iga01. Companion EASI model from the same paper: modellib('Briggs_2023_dupilumab_easi')."

  reference <- paste(
    "Briggs E, Kamal MA, Kosloski MP, Linsmeier I, Jusko N, Dolphin N,",
    "Chittenden J, Simpson EL, Paller AS, Siegfried EC, Shumel B, Levit NA,",
    "Bansal A, Davis JD, Chapel S, Smith DE, Huniti N. Integrated",
    "Exposure-Response of Dupilumab in Children, Adolescents, and Adults With",
    "Atopic Dermatitis Using Categorical and Continuous Efficacy Assessments:",
    "A Population Analysis. Pharm Res. 2023;40(11):2653-2666.",
    "doi:10.1007/s11095-023-03616-8.",
    sep = " "
  )
  vignette <- "Briggs_2023_dupilumab"

  units <- list(
    time = "day",
    dosing = "(none; dupilumab exposure enters as the CP_DUPILUMAB_MGL concentration column, not as dosing events)",
    concentration = "(ordinal IGA score 0-4; the observation iga is a category and the outputs p0-p4 and prob_iga01 are probabilities; the drug input CP_DUPILUMAB_MGL is in mg/L)"
  )

  compartmentData <- list(
    latent_iga = list(
      analyte = "latent indirect-response variable driving the IGA drug effect",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CP_DUPILUMAB_MGL = list(
      description = "Serum functional dupilumab concentration",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying input supplied on every record: 0 on placebo. In the source analysis it was the individual prediction of the age-group dupilumab popPK model (Briggs 2023 Online Resource 'Population Pharmacokinetic Model' and Table S1, which reprints Kovalenko 2020 / 2021), evaluated with each patient's actual dosing history. The paper's covariate and exposure-response simulations hold it constant (68 mg/L in Figure 3; a range of trough values in Figure 5).",
      source_name = "C (dupilumab concentration)"
    ),
    IGA_MOD = list(
      description = "Moderate baseline disease indicator (1 = baseline IGA 3, moderate; 0 = baseline IGA 4, severe)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (severe, baseline IGA 4)",
      notes = "Additive +4.39 shift of the probit-scale baseline cutpoint for IGA <= 3 and linear effect (1 - 0.727 * IGA_MOD) on Pmax (Table III). Only baseline IGA 3 or 4 is allowed: every patient in the analysis had moderate or severe disease.",
      source_name = "Disease severity (IGA)"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (AGE / 30)^0.0550 on Pmax, centred on the 30-year Figure 3 reference patient (form and centring not printed in the paper; chosen to reproduce Figure 3).",
      source_name = "Age"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (WT / 70)^-0.174 on the drug slope DSLP, centred on the 70 kg Figure 3 reference patient (form and centring not printed in the paper; chosen to reproduce Figure 3).",
      source_name = "Baseline body weight"
    ),
    EOS = list(
      description = "Baseline blood eosinophil count",
      units = "cells/uL",
      type = "continuous",
      reference_category = NULL,
      notes = "The paper reports eosinophils in 10^9/L; multiply by 1000 to get cells/uL. Table III 'Log baseline EOS on drug effect' = -0.0776, encoded as (EOS / 481)^-0.0776 on DSLP; 481 cells/uL (0.481 x 10^9/L) is the median observed value used for the Figure 3 reference patient.",
      source_name = "EOS"
    ),
    TARC = list(
      description = "Baseline serum thymus and activation-regulated chemokine (TARC / CCL17)",
      units = "pg/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Table III 'Log baseline TARC on drug effect' = -0.0938, encoded as (TARC / 2033.6)^-0.0938 on DSLP; 2033.6 pg/mL is the median observed value used for the Figure 3 reference patient.",
      source_name = "TARC"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = any other race)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian; the Figure 3 reference patient is White)",
      notes = "Linear effect (1 - 0.222 * RACE_ASIAN) on DSLP. Black and Other race were screened (Table S2) but not retained.",
      source_name = "Asian race"
    ),
    CONMED_STEROID_TOPICAL = list(
      description = "Concomitant topical corticosteroid (TCS) therapy indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (dupilumab or placebo without TCS)",
      notes = "Subject-level grain, set by the study design: all patients in the children study (NCT03345914) and the adult CHRONOS study (NCT02260986) received TCS, no patient in the other four studies did (Briggs 2023 Table I and Table S3). Linear effect (1 + 0.229 * TCS) on Pmax.",
      source_name = "TCS co-administration"
    ),
    PRIOR_SYSTEMIC_IMMUNOTHERAPY = list(
      description = "Prior exposure to systemic immunotherapy for atopic dermatitis (1 = yes)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no prior systemic immunotherapy)",
      notes = "Linear effect (1 - 0.102 * PRIOR) on Pmax. Briggs 2023 Table S3 tabulates the column as 'Prior exposure to systemic immunotherapies' (31.8% yes overall).",
      source_name = "Prior exposure to systemic immunotherapies"
    )
  )

  covariatesDataExcluded <- list(
    RACE_BLACK = list(
      description = "Black race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened with Other race on the placebo and drug-effect parameters (Briggs 2023 Table S2) and removed by backward elimination (alpha = 0.001)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "The Methods list gender among the baseline demographics considered, but Table S2 does not list it among the covariates evaluated on IGA and it is not in the final model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2968L,
    n_studies = 6L,
    age_range = "6-75 years (children 6-11 n = 359, adolescents 12-17 n = 243, adults 18-75 n = 2366); mean 32.2 (SD 16.4) years",
    weight_range = "mean 69.8 (SD 23.2) kg overall; adults 76.1, adolescents 65.0, children 31.5 kg",
    sex_female_pct = 41.9,
    race_ethnicity = "White 67.2%, Asian 21.1%, Black 8.1%, Other / not reported 3.7%; 9.0% Japanese (all adults)",
    disease_state = "Moderate-to-severe atopic dermatitis (baseline IGA 3 in 45.8%, IGA 4 in 54.2%; only 1 child had moderate disease)",
    dose_range = "Dupilumab SC with a Day-1 loading dose of twice the maintenance dose: adults 100-300 mg q4w to 300 mg qw; adolescents 200 or 300 mg q2w or 300 mg q4w; children 100 or 200 mg q2w or 300 mg q4w with TCS; or placebo (with or without TCS)",
    regions = "Multinational Phase IIb / III programme (NCT01859988, NCT02277769, NCT02277743, NCT02260986, NCT03054428, NCT03345914)",
    co_medication = "Topical corticosteroids in 35.6% (all children and the adult CHRONOS study)",
    notes = "Demographics from Briggs 2023 Table S3; 29,420 IGA observations. 66.2% received dupilumab, 33.8% placebo. Estimated with the Laplacian method in NONMEM 7.3."
  )

  ini({
    # Briggs 2023 Table III (final IGA model) and Online Resource equations:
    #   probit P(IGA <= m) = Q_m + Pmax * t / (ET50 + t) + DSLP * R
    #   Q_3 = BASE; Q_m = Q_(m+1) - exp(beta_(m+1)) for m = 0, 1, 2
    #   dR/dt = kout * (1 - C / (IC50 + C)) - kout * (R + 1), R(0) = 0
    # Table III prints the cutpoint gaps exp(beta) directly; they are carried
    # here on the log scale, so exp(lgap_*) is the printed value.
    rbase <- -3.35
    label("Probit-scale cutpoint for IGA <= 3 at baseline in severe disease (BASE; unitless)") # Table III 'Baseline for IGA <= 3' = -3.35 (RSE 3.6%)
    lgap_iga2 <- log(2.64)
    label("Log gap between the IGA <= 3 and IGA <= 2 probit cutpoints (unitless)") # Table III 'Baseline adjustment for IGA <= 2' = 2.64 (95% CI 2.59-2.69)
    lgap_iga1 <- log(1.60)
    label("Log gap between the IGA <= 2 and IGA <= 1 probit cutpoints (unitless)") # Table III 'Baseline adjustment for IGA <= 1' = 1.60 (95% CI 1.56-1.63)
    lgap_iga0 <- log(1.75)
    label("Log gap between the IGA <= 1 and IGA <= 0 probit cutpoints (unitless)") # Table III 'Baseline adjustment for IGA <= 0' = 1.75 (95% CI 1.70-1.79)
    lpmax <- log(4.83)
    label("Maximum placebo effect on the probit scale (unitless)") # Table III 'Pmax' = 4.83 (RSE 2.5%)
    let50 <- log(6.37)
    label("Time to half the maximum placebo effect (day)") # Table III 'ET50 (days)' = 6.37 (RSE 3.0%)
    slope <- -2.20
    label("Drug effect slope DSLP on the indirect-response signal R (probit units)") # Table III 'Drug effect (DSLP)' = -2.20 (RSE 3.1%)
    lic50 <- log(27.1)
    label("Dupilumab concentration giving half the maximum drug signal (mg/L)") # Table III 'IC50 (mg/L)' = 27.1 (95% CI 21.3-34.4)
    lthalf_onset <- log(19.7)
    label("Half-life of drug effect onset THFD; kout = log(2) / THFD (day)") # Table III 'Half-life for drug effect onset (day)' = 19.7 (95% CI 18.2-21.4)

    e_iga_mod_rbase <- 4.39
    label("Additive shift of the IGA <= 3 probit cutpoint for moderate baseline disease (unitless)") # Table III 'Moderate IGA baseline additive shift' = 4.39 (RSE 3.0%)
    e_age_pmax <- 0.0550
    label("Power exponent of AGE/30 on Pmax (unitless)") # Table III 'Age on Pmax' = 0.0550 (RSE 21.3%)
    e_conmed_steroid_topical_pmax <- 0.229
    label("Fractional change in Pmax with topical corticosteroids (unitless)") # Table III 'TCS on Pmax' = 0.229 (RSE 8.1%)
    e_prior_systemic_immunotherapy_pmax <- -0.102
    label("Fractional change in Pmax with prior systemic immunotherapy (unitless)") # Table III 'Prior immunotherapy on Pmax' = -0.102 (RSE 13.1%)
    e_iga_mod_pmax <- -0.727
    label("Fractional change in Pmax with moderate baseline disease (unitless)") # Table III 'Moderate disease severity on Pmax' = -0.727 (RSE 1.7%)
    e_wt_slope <- -0.174
    label("Power exponent of WT/70 on DSLP (unitless)") # Table III 'Baseline body weight on drug effect' = -0.174 (RSE 26.1%)
    e_race_asian_slope <- -0.222
    label("Fractional change in DSLP for Asian race (unitless)") # Table III 'Asian race on drug effect' = -0.222 (RSE 14.2%)
    e_eos_slope <- -0.0776
    label("Power exponent of EOS/481 on DSLP (unitless)") # Table III 'Log baseline EOS on drug effect' = -0.0776 (RSE 19.3%)
    e_tarc_slope <- -0.0938
    label("Power exponent of TARC/2033.6 on DSLP (unitless)") # Table III 'Log baseline TARC on drug effect' = -0.0938 (RSE 13.4%)

    # IIV: additive on the probit-scale baseline cutpoint, shared by every
    # cutpoint because the lower ones are built from it.
    etarbase ~ 1.49 # Table III 'IGA baseline-additive (omega2)' = 1.49 (RSE 3.4%; shrinkage 8.62%)
  })

  model({
    # Probit-scale cutpoints (Online Resource, IGA model equations)
    rbase_i <- rbase + e_iga_mod_rbase * IGA_MOD + etarbase
    q3 <- rbase_i
    q2 <- q3 - exp(lgap_iga2)
    q1 <- q2 - exp(lgap_iga1)
    q0 <- q1 - exp(lgap_iga0)

    # Placebo (non-drug) time course
    pmax_i <- exp(lpmax) *
      (AGE / 30)^e_age_pmax *
      (1 + e_conmed_steroid_topical_pmax * CONMED_STEROID_TOPICAL) *
      (1 + e_prior_systemic_immunotherapy_pmax * PRIOR_SYSTEMIC_IMMUNOTHERAPY) *
      (1 + e_iga_mod_pmax * IGA_MOD)
    et50 <- exp(let50)
    placebo <- pmax_i * time / (et50 + time)

    # Indirect-response drug signal. The state is the paper's R + 1, so it
    # starts at 1 and the paper's R is latent_iga - 1.
    slope_i <- slope *
      (WT / 70)^e_wt_slope *
      (1 + e_race_asian_slope * RACE_ASIAN) *
      (EOS / 481)^e_eos_slope *
      (TARC / 2033.6)^e_tarc_slope
    ic50 <- exp(lic50)
    thalf_onset <- exp(lthalf_onset)
    kout <- log(2) / thalf_onset
    latent_iga(0) <- 1
    d / dt(latent_iga) <- kout * (1 - CP_DUPILUMAB_MGL / (ic50 + CP_DUPILUMAB_MGL)) - kout * latent_iga
    r_drug <- latent_iga - 1

    # Cumulative and per-category probabilities
    lv <- placebo + slope_i * r_drug
    pcum0 <- pnorm(q0 + lv)
    pcum1 <- pnorm(q1 + lv)
    pcum2 <- pnorm(q2 + lv)
    pcum3 <- pnorm(q3 + lv)
    p0 <- pcum0
    p1 <- pcum1 - pcum0
    p2 <- pcum2 - pcum1
    p3 <- pcum3 - pcum2
    p4 <- 1 - pcum3
    prob_iga01 <- pcum1

    # rxode2 ordinal endpoint: categories 1-5 carry probabilities p0 to p3
    # and 1 - p0 - p1 - p2 - p3, i.e. simulated category k is IGA k - 1.
    iga ~ c(p0, p1, p2, p3)
  })
}
