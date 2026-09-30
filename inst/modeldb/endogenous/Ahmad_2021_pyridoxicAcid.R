Ahmad_2021_pyridoxicAcid <- function() {
  description <- "Coupled turnover model for the endogenous renal OAT1/OAT3 biomarker pyridoxic acid (PDA) in healthy adults (Ahmad 2021), fit simultaneously to PDA plasma concentrations and urinary amounts with and without the OAT1/3 inhibitor probenecid. PDA is produced at a zero-order synthesis rate ksyn and eliminated by renal clearance CLr (about 82% of total clearance) and nonrenal clearance CLnr into a one-compartment plasma pool; the renally cleared amount accumulates in a urine compartment. Probenecid competitively inhibits CLr through its total-plasma OAT1/3 inhibition constant Ki, driven by probenecid's own one-compartment first-order-absorption PK carried in the same file; synthesis and nonrenal clearance are unaffected. With no probenecid dosed the model sits at its steady-state baseline of about 3.2 ng/mL (typical value). Probenecid doses are given in umol (500 mg = 1752 umol)."
  reference <- paste(
    "Ahmad A, Ogungbenro K, Kunze A, Jacobs F, Snoeys J, Rostami-Hodjegan A,",
    "Galetin A. Population pharmacokinetic modeling and simulation to support",
    "qualification of pyridoxic acid as endogenous biomarker of OAT1/3 renal",
    "transporters. CPT Pharmacometrics Syst Pharmacol. 2021;10(5):467-477.",
    "doi:10.1002/psp4.12610.",
    "Structural equations from Eqs 2-3 of the main text; all parameter values",
    "from Table 1; clinical-study design and probenecid model structure from",
    "the Supplementary Material.",
    sep = " "
  )
  vignette <- "Ahmad_2021_oatBiomarkers"
  units <- list(time = "h", dosing = "umol", concentration = "ug/L")

  # What each ODE state holds, in what amount units, in what biological
  # matrix. Checked against Ahmad 2021 Figure 1 and Eqs 2-3.
  compartmentData <- list(
    central = list(analyte = "pyridoxic acid", units = "ug", specimen = "plasma", verified = TRUE),
    urine = list(analyte = "pyridoxic acid", units = "ug", specimen = "urine", verified = TRUE),
    depot_prob = list(analyte = "probenecid", units = "umol", specimen = "administration site", verified = TRUE),
    central_prob = list(analyte = "probenecid", units = "umol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    OCC = list(
      description = "Study occasion (crossover phase): 1 = control phase (no probenecid), 2 = probenecid interaction phase.",
      units = "(count)",
      type = "categorical",
      reference_category = "n/a -- decomposed into occasion indicators that select the per-occasion eta on ksyn",
      notes = "Carries only the interoccasion variability on the PDA synthesis rate ksyn (Ahmad 2021 Table 1, IOV column; 'IOV between the two interaction phases was estimated only for ksyn'). Two occasions, one per crossover phase, separated by a 21-day washout. Set OCC = 1 throughout for single-occasion simulation; an OCC value other than 1 or 2 switches the IOV term off (typical-value occasion).",
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
    disease_state = "Healthy volunteers. PDA was monitored as an endogenous biomarker of renal OAT1/OAT3 activity during a transporter drug-drug-interaction study.",
    dose_range = "Endogenous biomarker (no exogenous PDA dose). Perpetrator: probenecid 500 mg orally at 18:00 and 23:00 on the day before the interaction phase, then 500 mg about every 6 h (07:00, 13:00, 18:00, 23:00) for 7 days.",
    regions = "Europe (EudraCT 2016-003923-49)",
    renal_function = "CKD-EPI glomerular filtration rate 3.04-4.59 mL/min/kg (mean 3.92 mL/min/kg)",
    notes = "Two-phase crossover (Willemin et al. 2021) with a 21-day washout: phase I a single dose of a Janssen victim compound alone (used as the PDA baseline), phase II the same compound co-administered with multiple-dose probenecid. 132 PDA plasma samples (66 per phase) and 108 PDA urine samples (54 per phase; intervals 0-6, 6-12, 12-24, then 24-h collections to 168 h), plus 36 probenecid plasma trough samples. The probenecid model was fit first; its individual empirical Bayes estimates were then fixed while PDA plasma and urine data from both phases were fit simultaneously in NONMEM (FOCE-I). Independently verified against the single-dose 1000 mg probenecid study of Shen et al. 2019 (n = 14)."
  )

  ini({
    # ---------------------------------------------------------------
    # PDA biomarker -- Ahmad 2021 Table 1, PDA block. Turnover model
    # of Eq. 2 (plasma) and Eq. 3 (urine), structure adopted from
    # Barnett et al. 2018 (coproporphyrin I).
    # ---------------------------------------------------------------
    lksyn <- log(58.6)
    label("PDA zero-order synthesis rate ksyn (ug/h)")
    # Table 1, PDA row 'k syn, ug/h' = 58.6 (SE 26%).

    lvc <- log(9.34)
    label("PDA volume of distribution V (L)")
    # Table 1, PDA row 'V, L' = 9.34 (SE 56%).

    lcl_renal <- log(15.2)
    label("PDA renal clearance CLr (L/h)")
    # Table 1, PDA row 'CL r, L/h' = 15.2 (SE 15%).

    lcl_nonren <- log(3.23)
    label("PDA nonrenal clearance CLnr (L/h)")
    # Table 1, PDA row 'CL nr, L/h' = 3.23 (SE 41%). CLr / (CLr + CLnr) =
    # 15.2 / 18.43 = 82%, the renal fraction stated in the Results.

    lki_oat13 <- log(54.5)
    label("Probenecid total-plasma OAT1/3 inhibition constant Ki, PDA as probe (umol/L)")
    # Table 1, Probenecid row 'K i, uM (PDA data)' = 54.5 (SE 14%). The
    # Table 1 note states this is the TOTAL Ki; the unbound value 3.4 uM
    # (fu = 0.062) is reported only for comparison with in vitro data and
    # is not used by the model.

    # ---------------------------------------------------------------
    # Probenecid perpetrator PK -- Table 1, Probenecid block. One
    # compartment with first-order absorption and linear elimination
    # (Figure 1; Supplementary Material 'Probenecid population PK model').
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
    etalksyn ~ 0.0974856 # Table 1, PDA IIV ksyn = 32% (SE 34%); log(1 + 0.32^2)
    etalcl_renal ~ 0.0703190 # Table 1, PDA IIV CLr = 27% (SE 15%); log(1 + 0.27^2)
    etalvc_prob ~ 0.1484200 # Table 1, Probenecid IIV V = 40% (SE 29%); log(1 + 0.40^2)
    etalcl_prob ~ 0.0515460 # Table 1, Probenecid IIV CL = 23% (SE 23%); log(1 + 0.23^2)
    etaiov_ksyn_1 ~ 0.0472723 # Table 1, PDA IOV ksyn = 22% (SE 57%); log(1 + 0.22^2)
    etaiov_ksyn_2 ~ fixed(0.0472723) # Table 1, PDA IOV ksyn = 22%; shared variance, occasion 2

    # ---------------------------------------------------------------
    # Residual error -- Table 1. Combined proportional + additive on
    # every output (Methods, 'Structural and statistic models').
    # ---------------------------------------------------------------
    propSd <- 0.127
    label("Proportional residual error, PDA plasma (fraction)")
    # Table 1, PDA row 'sigma prop (%) - plasma' = 12.7 (SE 7%).

    addSd <- 0.213
    label("Additive residual error, PDA plasma (ug/L = ng/mL)")
    # Table 1, PDA row 'sigma add, ng/ml - plasma' = 0.213 (SE 18%).

    propSd_Upda <- 0.321
    label("Proportional residual error, PDA urine amount (fraction)")
    # Table 1, PDA row 'sigma prop (%) - urine' = 32.1 (SE 20%).

    addSd_Upda <- 58
    label("Additive residual error, PDA urine amount (ug)")
    # Table 1, PDA row 'sigma add, ug - urine' = 58 (SE 78%).

    propSd_prob <- 0.188
    label("Proportional residual error, probenecid plasma (fraction)")
    # Table 1, Probenecid row 'sigma prop, (%)' = 18.8 (SE 30%).

    addSd_prob <- fixed(0.001)
    label("Additive residual error, probenecid plasma (umol/L)")
    # Table 1, Probenecid row 'sigma add, uM' = 0.001 (Fixed).
  })

  model({
    # --- Occasion indicators (IOV on ksyn) --------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_ksyn <- oc1 * etaiov_ksyn_1 + oc2 * etaiov_ksyn_2

    # --- Individual parameters -------------------------------------
    ksyn <- exp(lksyn + etalksyn + iov_ksyn)
    vc <- exp(lvc)
    cl_renal <- exp(lcl_renal + etalcl_renal)
    cl_nonren <- exp(lcl_nonren)
    ki_oat13 <- exp(lki_oat13)

    ka_prob <- exp(lka_prob)
    vc_prob <- exp(lvc_prob + etalvc_prob)
    cl_prob <- exp(lcl_prob + etalcl_prob)

    # --- Probenecid PK (umol, umol/L) -------------------------------
    Cc_prob <- central_prob / vc_prob

    # --- PDA turnover ----------------------------------------------
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
    # Css = ksyn / (CLr + CLnr); typical 58.6 / 18.43 = 3.18 ug/L.
    central(0) <- ksyn / (cl_renal + cl_nonren) * vc

    # --- Outputs ----------------------------------------------------
    # Upda is the cumulative amount excreted since time 0; an observed
    # collection-interval amount is the difference of Upda between the
    # interval end points.
    Upda <- urine

    Cc ~ add(addSd) + prop(propSd)
    Upda ~ add(addSd_Upda) + prop(propSd_Upda)
    Cc_prob ~ add(addSd_prob) + prop(propSd_prob)
  })
}
