Ahmad_2021_homovanillicAcid <- function() {
  description <- "Coupled turnover model for the endogenous renal OAT1/OAT3 biomarker homovanillic acid (HVA) in healthy adults (Ahmad 2021), fit simultaneously to HVA plasma concentrations and urinary amounts with and without the OAT1/3 inhibitor probenecid. HVA is produced at a zero-order synthesis rate ksyn and eliminated by renal clearance CLr (about 94% of total clearance) and nonrenal clearance CLnr into a one-compartment plasma pool; the renally cleared amount accumulates in a urine compartment. Probenecid competitively inhibits CLr through its total-plasma OAT1/3 inhibition constant Ki, driven by probenecid's own one-compartment first-order-absorption PK carried in the same file; synthesis and nonrenal clearance are unaffected. With no probenecid dosed the model sits at its steady-state baseline of about 9.8 ng/mL (typical value). Probenecid doses are given in umol (500 mg = 1752 umol). Companion to Ahmad_2021_pyridoxicAcid, which shares the probenecid PK."
  reference <- paste(
    "Ahmad A, Ogungbenro K, Kunze A, Jacobs F, Snoeys J, Rostami-Hodjegan A,",
    "Galetin A. Population pharmacokinetic modeling and simulation to support",
    "qualification of pyridoxic acid as endogenous biomarker of OAT1/3 renal",
    "transporters. CPT Pharmacometrics Syst Pharmacol. 2021;10(5):467-477.",
    "doi:10.1002/psp4.12610.",
    "Structural equations from Eqs 2-3 of the main text (stated to apply",
    "equally to HVA); all parameter values from Table 1; clinical-study",
    "design and probenecid model structure from the Supplementary Material.",
    sep = " "
  )
  vignette <- "Ahmad_2021_oatBiomarkers"
  units <- list(time = "h", dosing = "umol", concentration = "ug/L")

  # What each ODE state holds, in what amount units, in what biological
  # matrix. Checked against Ahmad 2021 Figure 1 and Eqs 2-3.
  compartmentData <- list(
    central = list(analyte = "homovanillic acid", units = "ug", specimen = "plasma", verified = TRUE),
    urine = list(analyte = "homovanillic acid", units = "ug", specimen = "urine", verified = TRUE),
    depot_prob = list(analyte = "probenecid", units = "umol", specimen = "administration site", verified = TRUE),
    central_prob = list(analyte = "probenecid", units = "umol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    OCC = list(
      description = "Study occasion (crossover phase): 1 = control phase (no probenecid), 2 = probenecid interaction phase.",
      units = "(count)",
      type = "categorical",
      reference_category = "n/a -- decomposed into occasion indicators that select the per-occasion eta on CLr",
      notes = "Carries only the interoccasion variability on the HVA renal clearance CLr (Ahmad 2021 Table 1, HVA IOV column). Two occasions, one per crossover phase, separated by a 21-day washout. Set OCC = 1 throughout for single-occasion simulation; an OCC value other than 1 or 2 switches the IOV term off (typical-value occasion).",
      source_name = "phase (control / probenecid)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 6L,
    n_studies = 1L,
    age_range = "48-58 years",
    weight_range = "(not reported; body mass index 22-28 kg/m^2)",
    sex_female_pct = 100,
    race_ethnicity = "Caucasian (6 of 6)",
    disease_state = "Healthy volunteers. HVA was monitored as an endogenous biomarker of renal OAT1/OAT3 activity during a transporter drug-drug-interaction study.",
    dose_range = "Endogenous biomarker (no exogenous HVA dose). Perpetrator: probenecid 500 mg orally at 18:00 and 23:00 on the day before the interaction phase, then 500 mg about every 6 h (07:00, 13:00, 18:00, 23:00) for 7 days.",
    regions = "Europe (EudraCT 2016-003923-49)",
    renal_function = "CKD-EPI glomerular filtration rate 3.04-4.59 mL/min/kg (mean 3.92 mL/min/kg)",
    notes = "Two-phase crossover (Willemin et al. 2021) with a 21-day washout: phase I a single dose of a Janssen victim compound alone (used as the HVA baseline), phase II the same compound co-administered with multiple-dose probenecid. 132 HVA plasma samples (66 per phase) and 108 HVA urine samples (54 per phase), plus 36 probenecid plasma trough samples. The probenecid model was fit first; its individual empirical Bayes estimates were then fixed while HVA plasma and urine data from both phases were fit simultaneously in NONMEM (FOCE-I)."
  )

  ini({
    # ---------------------------------------------------------------
    # HVA biomarker -- Ahmad 2021 Table 1, HVA block. Turnover model
    # of Eq. 2 (plasma) and Eq. 3 (urine) applied to HVA ('A similar
    # procedure was done in the case of HVA').
    # ---------------------------------------------------------------
    lksyn <- log(212)
    label("HVA zero-order synthesis rate ksyn (ug/h)")
    # Table 1, HVA row 'k syn, ug/h' = 212 (SE 16%).

    lvc <- log(105)
    label("HVA volume of distribution V (L)")
    # Table 1, HVA row 'V, L' = 105 (SE 35%).

    lcl_renal <- log(20.4)
    label("HVA renal clearance CLr (L/h)")
    # Table 1, HVA row 'CL r, L/h' = 20.4 (SE 14%).

    lcl_nonren <- log(1.29)
    label("HVA nonrenal clearance CLnr (L/h)")
    # Table 1, HVA row 'CL nr, L/h' = 1.29 (SE 69%). CLr / (CLr + CLnr) =
    # 20.4 / 21.69 = 94%, the renal fraction stated in the Results.

    lki_oat13 <- log(137)
    label("Probenecid total-plasma OAT1/3 inhibition constant Ki, HVA as probe (umol/L)")
    # Table 1, Probenecid row 'K i, uM (HVA data)' = 137 (SE 24%). The
    # Table 1 note states this is the TOTAL Ki; the unbound value 8.5 uM
    # (fu = 0.062) is reported only for comparison with in vitro data and
    # is not used by the model.

    # ---------------------------------------------------------------
    # Probenecid perpetrator PK -- Table 1, Probenecid block (shared
    # with the PDA model). One compartment with first-order absorption
    # and linear elimination (Figure 1; Supplementary Material).
    # ---------------------------------------------------------------
    lka_prob <- log(0.74)
    label("Probenecid absorption rate constant ka (1/h)")
    # Table 1, Probenecid row 'k a /h' = 0.74 (SE 45%).

    lvc_prob <- log(15)
    label("Probenecid apparent volume of distribution V (L)")
    # Table 1, Probenecid row 'V, L' = 15 (SE 22%).

    lcl_prob <- log(0.82)
    label("Probenecid apparent clearance CL (L/h)")
    # Table 1, Probenecid row 'CL, L/h' = 0.82 (SE 9%).

    # ---------------------------------------------------------------
    # Between-subject and between-occasion variability. Table 1
    # reports IIV / IOV as percentages of an exponential (log-normal)
    # random effect; read as CV% and converted with
    # omega^2 = log(1 + CV^2).
    # ---------------------------------------------------------------
    etalksyn ~ 0.0703190 # Table 1, HVA IIV ksyn = 27% (SE 19%); log(1 + 0.27^2)
    etalcl_renal ~ 0.0515460 # Table 1, HVA IIV CLr = 23% (SE 26%); log(1 + 0.23^2)
    etalvc_prob ~ 0.1484200 # Table 1, Probenecid IIV V = 40% (SE 29%); log(1 + 0.40^2)
    etalcl_prob ~ 0.0515460 # Table 1, Probenecid IIV CL = 23% (SE 23%); log(1 + 0.23^2)
    etaiov_cl_renal_1 ~ 0.0024969 # Table 1, HVA IOV CLr = 5% (SE 48%); log(1 + 0.05^2)
    etaiov_cl_renal_2 ~ fixed(0.0024969) # Table 1, HVA IOV CLr = 5%; shared variance, occasion 2

    # ---------------------------------------------------------------
    # Residual error -- Table 1. Combined proportional + additive on
    # every output (Methods, 'Structural and statistic models').
    # ---------------------------------------------------------------
    propSd <- 0.173
    label("Proportional residual error, HVA plasma (fraction)")
    # Table 1, HVA row 'sigma prop (%) - plasma' = 17.3 (SE 8%).

    addSd <- fixed(0.001)
    label("Additive residual error, HVA plasma (ug/L = ng/mL)")
    # Table 1, HVA row 'sigma add, ng/ml - plasma' = 0.001 Fixed.

    propSd_Uhva <- 0.269
    label("Proportional residual error, HVA urine amount (fraction)")
    # Table 1, HVA row 'sigma prop (%) - urine' = 26.9 (SE 14%).

    addSd_Uhva <- 373
    label("Additive residual error, HVA urine amount (ug)")
    # Table 1, HVA row 'sigma add, ug - urine' = 373 (SE 20%).

    propSd_prob <- 0.188
    label("Proportional residual error, probenecid plasma (fraction)")
    # Table 1, Probenecid row 'sigma prop, (%)' = 18.8 (SE 30%).

    addSd_prob <- fixed(0.001)
    label("Additive residual error, probenecid plasma (umol/L)")
    # Table 1, Probenecid row 'sigma add, uM' = 0.001 (Fixed).
  })

  model({
    # --- Occasion indicators (IOV on CLr) ---------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl_renal <- oc1 * etaiov_cl_renal_1 + oc2 * etaiov_cl_renal_2

    # --- Individual parameters -------------------------------------
    ksyn <- exp(lksyn + etalksyn)
    vc <- exp(lvc)
    cl_renal <- exp(lcl_renal + etalcl_renal + iov_cl_renal)
    cl_nonren <- exp(lcl_nonren)
    ki_oat13 <- exp(lki_oat13)

    ka_prob <- exp(lka_prob)
    vc_prob <- exp(lvc_prob + etalvc_prob)
    cl_prob <- exp(lcl_prob + etalcl_prob)

    # --- Probenecid PK (umol, umol/L) -------------------------------
    Cc_prob <- central_prob / vc_prob

    # --- HVA turnover ----------------------------------------------
    # Ahmad 2021 Eqs 2-3: probenecid competitively inhibits CLr only;
    # with no probenecid on board cl_renal_eff = cl_renal.
    cl_renal_eff <- cl_renal / (1 + Cc_prob / ki_oat13)

    Cc <- central / vc

    # ODE declaration order sets compartment numbering: biomarker
    # states first (central = 1, urine = 2), then the perpetrator
    # (depot_prob = 3, central_prob = 4).
    d/dt(central) <- ksyn - (cl_renal_eff + cl_nonren) * Cc
    d/dt(urine) <- cl_renal_eff * Cc
    d/dt(depot_prob) <- -ka_prob * depot_prob
    d/dt(central_prob) <- ka_prob * depot_prob - cl_prob / vc_prob * central_prob

    # Inhibitor-free steady-state baseline, Eq. 2 with dC/dt = 0:
    # Css = ksyn / (CLr + CLnr); typical 212 / 21.69 = 9.77 ug/L.
    central(0) <- ksyn / (cl_renal + cl_nonren) * vc

    # --- Outputs ----------------------------------------------------
    # Uhva is the cumulative amount excreted since time 0; an observed
    # collection-interval amount is the difference of Uhva between the
    # interval end points.
    Uhva <- urine

    Cc ~ add(addSd) + prop(propSd)
    Uhva ~ add(addSd_Uhva) + prop(propSd_Uhva)
    Cc_prob ~ add(addSd_prob) + prop(propSd_prob)
  })
}
