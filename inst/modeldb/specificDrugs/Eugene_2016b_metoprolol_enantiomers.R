Eugene_2016b_metoprolol_enantiomers <- function() {
  description <- "Enantiomer-resolved one-compartment PK model for oral metoprolol in healthy young adults, tracking R-metoprolol and S-metoprolol in parallel with first-order absorption and an absorption lag time and with every structural parameter stratified by sex. Typical values only (no inter-individual variability was estimated). The kinetics are flip-flop: absorption is rate-limiting, so the terminal slope is set by Ka rather than by CL/V. These are the parameters that generate every published Cmax, Tmax, AUC and half-life in Eugene 2016 (Eugene 2016)."
  reference <- paste(
    "Eugene AR. (2016).",
    "Metoprolol Dose Equivalence in Adult Men and Women Based on Gender",
    "Differences: Pharmacokinetic Modeling and Simulations.",
    "Med Sci (Basel) 4(4):18.",
    "doi:10.3390/medsci4040018.",
    sep = " "
  )
  vignette <- "Eugene_2016b_metoprolol_gender_dose_equivalence"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. No interconversion between the enantiomers is modelled;
  # each is dosed into, and eliminated from, its own pair of compartments.
  compartmentData <- list(
    depot_r_enant = list(analyte = "R-metoprolol", units = "mg", specimen = "administration site", verified = TRUE),
    central_r_enant = list(analyte = "R-metoprolol", units = "mg", specimen = "plasma", verified = TRUE),
    depot_s_enant = list(analyte = "S-metoprolol", units = "mg", specimen = "administration site", verified = TRUE),
    central_s_enant = list(analyte = "S-metoprolol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    SEXF = list(
      description = "Biological sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male). Eugene 2016 Table 1 reports a Female and a Male column for each enantiomer rather than a reference plus an effect; the male column is taken as the structural reference here and each female value is reproduced exactly by a log-additive offset.",
      notes = "Eugene 2016 Table 1 ('One-compartment pharmacokinetic parameters for R- and S-metoprolol for young men and women') prints four independent parameter sets -- R/Female, R/Male, S/Female, S/Male -- with no common reference value and no published covariate coefficient. To express this under the canonical SEXF (1 = female, 0 = male), each enantiomer's male column is the reference and each e_sexf_* term is computed as log(female value / male value), reproducing the published female columns to full printed precision. The sex effects are estimated separately per enantiomer and are NOT shared between them.",
      source_name = "gender"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    age_range = "20-36 years",
    age_median = "not reported",
    weight_range = "men 83.9 +/- 10.7 kg, women 62.0 +/- 7.3 kg (mean +/- SD, Luzier 1999 as quoted in Eugene 2016 Discussion)",
    weight_median = "not reported (means reported instead)",
    sex_female_pct = 50,
    race_ethnicity = "not reported",
    disease_state = "Healthy volunteers",
    dose_range = "Nine 100 mg oral doses of racemic metoprolol every 12 h (Eugene 2016 Methods 2.1, citing the Luzier 1999 protocol).",
    regions = "United States (Luzier 1999 cohort; Eugene 2016 analysis performed at Mayo Clinic, Rochester MN)",
    notes = "n = 20 healthy volunteers (10 men, 10 women) from Luzier 1999. Eugene 2016 did not have access to the individual-level data: the mean R-metoprolol and S-metoprolol concentration-time curves were DIGITIZED from the Luzier 1999 figures, giving four distinct profiles (enantiomer x sex), and each was fit in MONOLIX 4.3.3 (SAEM-MCMC). Because the fits are to mean curves rather than individual subjects, no inter-individual variability and no residual-error magnitude are reported for this model -- see the vignette Errata. Eugene 2016 undertook the fit because the original Luzier 1999 publication did not report Ka or Tlag (Methods 2.2)."
  )

  ini({
    # ---- R-metoprolol, Eugene 2016 Table 1, Male column (structural reference).
    ltlag_r_enant <- log(0.59); label("R-metoprolol absorption lag time in men (h)") # Eugene 2016 Table 1, R-Metoprolol Male: Tlag = 0.59 h
    lka_r_enant <- log(0.234); label("R-metoprolol first-order absorption rate constant in men (1/h)") # Eugene 2016 Table 1, R-Metoprolol Male: Ka = 0.234 1/h
    lvc_r_enant <- log(63.9); label("R-metoprolol apparent central volume Vc/F in men (L)") # Eugene 2016 Table 1, R-Metoprolol Male: V = 63.9 L
    lcl_r_enant <- log(316); label("R-metoprolol apparent clearance CL/F in men (L/h)") # Eugene 2016 Table 1, R-Metoprolol Male: CL = 316 L/h

    # ---- R-metoprolol female-sex offsets, reproducing Table 1's Female column.
    e_sexf_ltlag_r_enant <- log(0.39 / 0.59); label("Log-additive female-sex effect on R-metoprolol lag time") # Eugene 2016 Table 1, R-Metoprolol Female: Tlag = 0.39 h
    e_sexf_lka_r_enant <- log(0.165 / 0.234); label("Log-additive female-sex effect on R-metoprolol Ka") # Eugene 2016 Table 1, R-Metoprolol Female: Ka = 0.165 1/h
    e_sexf_lvc_r_enant <- log(38.1 / 63.9); label("Log-additive female-sex effect on R-metoprolol Vc/F") # Eugene 2016 Table 1, R-Metoprolol Female: V = 38.1 L
    e_sexf_lcl_r_enant <- log(120 / 316); label("Log-additive female-sex effect on R-metoprolol CL/F") # Eugene 2016 Table 1, R-Metoprolol Female: CL = 120 L/h

    # ---- S-metoprolol, Eugene 2016 Table 1, Male column (structural reference).
    ltlag_s_enant <- log(0.67); label("S-metoprolol absorption lag time in men (h)") # Eugene 2016 Table 1, S-Metoprolol Male: Tlag = 0.67 h
    lka_s_enant <- log(0.241); label("S-metoprolol first-order absorption rate constant in men (1/h)") # Eugene 2016 Table 1, S-Metoprolol Male: Ka = 0.241 1/h
    lvc_s_enant <- log(55.3); label("S-metoprolol apparent central volume Vc/F in men (L)") # Eugene 2016 Table 1, S-Metoprolol Male: V = 55.3 L
    lcl_s_enant <- log(253); label("S-metoprolol apparent clearance CL/F in men (L/h)") # Eugene 2016 Table 1, S-Metoprolol Male: CL = 253 L/h

    # ---- S-metoprolol female-sex offsets, reproducing Table 1's Female column.
    e_sexf_ltlag_s_enant <- log(0.38 / 0.67); label("Log-additive female-sex effect on S-metoprolol lag time") # Eugene 2016 Table 1, S-Metoprolol Female: Tlag = 0.38 h
    e_sexf_lka_s_enant <- log(0.161 / 0.241); label("Log-additive female-sex effect on S-metoprolol Ka") # Eugene 2016 Table 1, S-Metoprolol Female: Ka = 0.161 1/h
    e_sexf_lvc_s_enant <- log(34.9 / 55.3); label("Log-additive female-sex effect on S-metoprolol Vc/F") # Eugene 2016 Table 1, S-Metoprolol Female: V = 34.9 L
    e_sexf_lcl_s_enant <- log(101 / 253); label("Log-additive female-sex effect on S-metoprolol CL/F") # Eugene 2016 Table 1, S-Metoprolol Female: CL = 101 L/h

    # Residual error is fixed to zero because Eugene 2016 Table 1 reports no
    # residual-error magnitude for these fits: they are fits to DIGITIZED MEAN
    # concentration-time curves, so there is no within-subject residual to
    # report. Encoded as fixed(0) rather than invented, per the no-fabrication
    # rule; see the vignette Errata.
    propSd_r_enant <- fixed(0); label("R-metoprolol proportional residual SD (fraction; 0 -- not reported in the source)")
    propSd_s_enant <- fixed(0); label("S-metoprolol proportional residual SD (fraction; 0 -- not reported in the source)")
  })

  model({
    # ---- R-metoprolol individual parameters (typical values; no IIV reported).
    tlag_r_enant <- exp(ltlag_r_enant + e_sexf_ltlag_r_enant * SEXF)
    ka_r_enant <- exp(lka_r_enant + e_sexf_lka_r_enant * SEXF)
    vc_r_enant <- exp(lvc_r_enant + e_sexf_lvc_r_enant * SEXF)
    cl_r_enant <- exp(lcl_r_enant + e_sexf_lcl_r_enant * SEXF)
    kel_r_enant <- cl_r_enant / vc_r_enant

    # ---- S-metoprolol individual parameters (typical values; no IIV reported).
    tlag_s_enant <- exp(ltlag_s_enant + e_sexf_ltlag_s_enant * SEXF)
    ka_s_enant <- exp(lka_s_enant + e_sexf_lka_s_enant * SEXF)
    vc_s_enant <- exp(lvc_s_enant + e_sexf_lvc_s_enant * SEXF)
    cl_s_enant <- exp(lcl_s_enant + e_sexf_lcl_s_enant * SEXF)
    kel_s_enant <- cl_s_enant / vc_s_enant

    d/dt(depot_r_enant) <- -ka_r_enant * depot_r_enant
    d/dt(central_r_enant) <- ka_r_enant * depot_r_enant - kel_r_enant * central_r_enant

    d/dt(depot_s_enant) <- -ka_s_enant * depot_s_enant
    d/dt(central_s_enant) <- ka_s_enant * depot_s_enant - kel_s_enant * central_s_enant

    lag(depot_r_enant) <- tlag_r_enant
    lag(depot_s_enant) <- tlag_s_enant

    # Each enantiomer's CL/F and V/F were fit against the FULL administered
    # metoprolol dose, not against half of it: Eugene 2016 Results 3.2 reports
    # AUC = 394 ng*h/mL for men from a 100 mg dose, and 100 mg / 253 L/h =
    # 395 ng*h/mL reproduces it. Dose each depot with the full dose amount.
    # Dose is in mg and vc is in L, so central / vc is mg/L = 1000 ng/mL.
    Cc_r_enant <- 1000 * central_r_enant / vc_r_enant
    Cc_s_enant <- 1000 * central_s_enant / vc_s_enant

    Cc_r_enant ~ prop(propSd_r_enant)
    Cc_s_enant ~ prop(propSd_s_enant)
  })
}
