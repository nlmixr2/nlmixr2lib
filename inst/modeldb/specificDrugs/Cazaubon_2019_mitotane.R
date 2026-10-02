Cazaubon_2019_mitotane <- function() {
  description <- paste(
    "One-compartment population pharmacokinetic model for oral mitotane",
    "(o,p'-DDD) in adults with adrenocortical carcinoma (Cazaubon 2019;",
    "38 patients, 503 therapeutic-drug-monitoring plasma concentrations,",
    "Monolix 2019R1 SAEM). First-order absorption (ka fixed at 24 /day)",
    "with bioavailability fixed at 35 %, and linear elimination. Clearance",
    "carries power-form effects of serum triglycerides and HDL cholesterol",
    "(both lower clearance) and a two-class latent-covariate mixture: an",
    "'ultrafast metabolizer' subpopulation (11.5 % of subjects) with a",
    "3.06-fold higher clearance. Log-normal IIV on V and CL and a combined",
    "additive + proportional residual error. Time is in days."
  )
  reference <- paste(
    "Cazaubon Y, Talineau Y, Feliu C, Konecki C, Russello J, Mathieu O,",
    "Djerada Z. Population Pharmacokinetics Modelling and Simulation of",
    "Mitotane in Patients with Adrenocortical Carcinoma: An Individualized",
    "Dose Regimen to Target All Patients at Three Months?",
    "Pharmaceutics. 2019;11(11):566. doi:10.3390/pharmaceutics11110566."
  )
  vignette <- "Cazaubon_2019_mitotane"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "mitotane", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mitotane", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    TRIG = list(
      description = "Serum triglyceride concentration, per-subject median over the observation window (time-fixed).",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL: (TRIG / 1.56)^e_trig_cl with e_trig_cl = -0.526",
        "(Cazaubon 2019 Table 2 and footnote equation). The reference 1.56 g/L",
        "is the cohort median (Table 1, range 0.48-5.24 g/L). The paper used",
        "the median of each subject's values ('As variation of covariates were",
        "not significate, we took the median of each for each individual',",
        "Methods 2.5). Note the g/L unit: 1 g/L = 100 mg/dL ~= 1.13 mmol/L."
      ),
      source_name = "Tg"
    ),
    HDLC = list(
      description = "Serum high-density lipoprotein cholesterol concentration, per-subject median over the observation window (time-fixed).",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL: (HDLC / 0.65)^e_hdlc_cl with e_hdlc_cl = -0.344",
        "(Cazaubon 2019 Table 2 and footnote equation). The reference 0.65 g/L",
        "is the cohort median (Table 1, range 0.26-2.1 g/L). Per-subject median",
        "per Methods 2.5. Note the g/L unit: 1 g/L = 100 mg/dL ~= 2.59 mmol/L."
      ),
      source_name = "HDL"
    ),
    MIX_FAST_ELIM = list(
      description = paste(
        "Per-subject latent mixture-model class indicator: 1 = subject in the",
        "'ultrafast metabolizer' subpopulation (Monolix latent covariate",
        "modality lcat = 2, CL 3.06-fold higher); 0 = subject in the majority",
        "subpopulation (lcat = 1). Time-fixed."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (lcat = 1 majority subpopulation, 88.5 % of subjects)",
      notes = paste(
        "Not a measured clinical covariate -- it is the Monolix latent",
        "covariate (Cazaubon 2019 Methods 2.5 equation 3:",
        "log(Cl) = log(Clpop) + beta_Cl,lcat2 * 1[lcat = 2] + eta_Cl).",
        "Population probabilities: plcat_1 = 0.885, plcat_2 = 0.115 (Table 2).",
        "For population simulation draw MIX_FAST_ELIM ~ Bernoulli(0.115). The",
        "paper's dose-regimen simulations (Figure 4a-d, 4g-h) were run for the",
        "majority subpopulation (MIX_FAST_ELIM = 0), and Figure 4f for",
        "'plcat2 = 100%' (MIX_FAST_ELIM = 1); see the vignette. The paper",
        "speculates the class relates to the CYP2B6 G516T polymorphism but",
        "did not test it."
      ),
      source_name = "lcat (lcat2 modality)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (Methods 2.5), not retained in the final model."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Methods 2.5), not retained in the final model."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Methods 2.5), not retained in the final model."
    ),
    IBW = list(
      description = "Ideal body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Methods 2.5), not retained in the final model."
    ),
    LBW = list(
      description = "Lean body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Methods 2.5), not retained in the final model."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened (Methods 2.5), not retained in the final model."
    ),
    LDLC = list(
      description = "LDL cholesterol",
      units = "g/L",
      type = "continuous",
      notes = "Screened (Methods 2.5); negative trend on CL not significant (P = 0.056, Discussion); not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 38L,
    n_studies = 2L,
    n_observations = 503L,
    age_range = "14-76 years; median 51 (Table 1)",
    weight_range = "39-139 kg; median 71.7 (Table 1)",
    sex_female_pct = 100 * 11 / 38,
    disease_state = "Adrenocortical carcinoma treated with oral mitotane (therapeutic initiation), routine TDM",
    dose_range = "1-7.25 g/day orally in 2, 3 or 4 administrations (median 2.9 g/day); eight patients received 7.5-12 g/day during part of treatment (Table 1)",
    regions = "France (Reims University Hospital n = 25; Montpellier University Hospital n = 13)",
    notes = paste(
      "Retrospective pooled TDM data collected 2008-2016; at least four",
      "samples per patient (median 9, range 4-46). Covariate medians",
      "(Table 1): BMI 25.8 kg/m^2, CrCL 103 mL/min, HDL 0.65 g/L, LDL",
      "1.61 g/L, TG 1.56 g/L. Mitotane quantified by HPLC-UV."
    )
  )

  ini({
    # Structural parameters: Cazaubon 2019 Table 2 (final model, Monolix 2019R1).
    lka <- fixed(log(24)); label("Absorption rate constant (1/day)") # Table 2 'Ka (day-1) 24 FIX'; Methods 2.4
    lfdepot <- fixed(log(0.35)); label("Oral bioavailability (fraction)") # Table 2 'F (%) 35 FIX'; Methods 2.4
    lvc <- log(8900); label("Central volume of distribution (L)") # Table 2 'V (L) 8900 (18.2)'
    lcl <- log(70); label("Clearance at median TRIG and HDLC, majority subpopulation (L/day)") # Table 2 'Cl (L day-1) 70 (6.64)'

    # Covariate effects on CL: Table 2 footnote
    # Cli = Clpop * (Tgi/1.56)^beta_Tg * (HDLi/0.65)^beta_HDL * exp(beta_lcat2)
    e_trig_cl <- -0.526; label("Power exponent of TRIG/1.56 on CL (unitless)") # Table 2 'beta Tg -0.526 (25.5)'
    e_hdlc_cl <- -0.344; label("Power exponent of HDLC/0.65 on CL (unitless)") # Table 2 'beta HDL -0.344 (45.9)'
    e_mix_fast_elim_cl <- 1.12; label("Log-scale shift in CL for the lcat2 ultrafast subpopulation (unitless)") # Table 2 'beta lcat2 1.12 (20.1)'

    # IIV: Monolix reports omega as the SD of the log-normal random effect;
    # Table 2 prints them as percentages (omega x 100). Variance = omega^2.
    etalvc ~ 0.817216 # Table 2 'omega V (%) 90.4 (17.5)' -> 0.904^2
    etalcl ~ 0.085849 # Table 2 'omega Cl (%) 29.3 (16.7)' -> 0.293^2

    # Residual error: Table 2 'a (constant) 1.06 (13.5)' and 'b (proportional) 0.17 (8.51)'
    addSd <- 1.06; label("Additive residual error (mg/L)") # Table 2 'a (constant) 1.06'
    propSd <- 0.17; label("Proportional residual error (fraction)") # Table 2 'b (proportional) 0.17'
  })
  model({
    # Individual parameters (Methods 2.5 equations 1 and 3; Table 2 footnote)
    ka <- exp(lka)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + e_mix_fast_elim_cl * MIX_FAST_ELIM + etalcl) *
      (TRIG / 1.56)^e_trig_cl * (HDLC / 0.65)^e_hdlc_cl

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- exp(lfdepot)

    Cc <- central / vc
    # Monolix 2019R1 default 'combined1' error model: SD = a + b * f
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
