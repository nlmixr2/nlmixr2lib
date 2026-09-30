Sime_2019b_ceftolozane <- function() {
  description <- paste(
    "Four-compartment intravenous population PK model for UNBOUND ceftolozane in",
    "six critically ill adults undergoing continuous venovenous",
    "hemodiafiltration (CVVHDF) in a quaternary referral intensive care unit in",
    "Brisbane, Australia. Fitted non-parametrically in Pmetrics (NPAG) to",
    "unbound prefilter plasma, postfilter plasma and CVVHDF effluent",
    "concentrations simultaneously. A prefilter central compartment receives the",
    "dose and loses drug by residual non-CVVHDF clearance and by CVVHDF",
    "clearance into an effluent compartment, which drains first-order; it",
    "exchanges with a postfilter compartment by the rate constants K12 / K21 and",
    "with a peripheral compartment by an intercompartmental clearance, and the",
    "postfilter compartment also exchanges with the same peripheral compartment",
    "(paper Figure 1). No covariate was retained. Every parameter carries its",
    "own inter-individual variability taken from the mean and SD of the NPAG",
    "support-point distribution; residual error is fixed(0) because the selected",
    "Pmetrics error model was not published. Ceftolozane and tazobactam were",
    "fitted as separate models and are supplied as two files; see",
    "modellib('Sime_2019b_tazobactam') for the partner component of the fixed",
    "2:1 ceftolozane-tazobactam combination.",
    sep = " "
  )
  reference <- paste(
    "Sime FB, Lassig-Smith M, Starr T, Stuart J, Pandey S, Parker SL, Wallis SC,",
    "Lipman J, Roberts JA. A population pharmacokinetic model-guided evaluation",
    "of ceftolozane-tazobactam dosing in critically ill patients undergoing",
    "continuous venovenous hemodiafiltration. Antimicrob Agents Chemother.",
    "2019 (published 20 December 2019; issue January 2020);64(1):e01655-19.",
    "doi:10.1128/AAC.01655-19. PMCID: PMC7187594.",
    "Structure: Figure 1. All parameter values: Table 3, 'Ceftolozane' columns.",
    sep = " "
  )
  vignette <- "Sime_2019b_ceftolozane_tazobactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The postfilter compartment is a fitted, empirical compartment of the
  # paper's Figure 1 (volume ~18 L), not the physical blood-side volume of the
  # hemofilter, and it is NOT the phase drug is cleared from (CVVHDF clearance
  # leaves the PREFILTER compartment), so the canonical circuit_plasma role
  # does not fit it.
  paper_specific_compartments <- c("postfilter")

  # What each ODE state holds. The assay measured UNBOUND concentrations
  # directly (Methods, 'Ceftolozane and tazobactam assay': unbound fraction
  # isolated by ultracentrifugation), and the whole administered milligram
  # amount was dosed into the model, so every volume here relates total drug
  # amount to an UNBOUND concentration. The effluent compartment is fed by the
  # CVVHDF clearance and drains at Kd, the dialysate-side role of the canonical
  # circuit_dialysate state; what drains out of it (the effluent bag) is not a
  # state of the published model.
  compartmentData <- list(
    central = list(analyte = "ceftolozane", units = "mg", specimen = "plasma", verified = TRUE),
    postfilter = list(analyte = "ceftolozane", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ceftolozane", units = "mg", specimen = "plasma", verified = TRUE),
    circuit_dialysate = list(analyte = "ceftolozane", units = "mg", specimen = "dialysate", verified = TRUE)
  )

  covariateData <- list()

  # Methods, 'Pharmacokinetic analysis': 'Available covariates considered for
  # analysis included sex, height, weight, body max index, body surface area,
  # albumin concentration, serum creatinine, sequential organ failure
  # assessment (SOFA) score, acute physiology and chronic health evaluation
  # (APACHE) II score, dialysate flow rate, transmembrane pressure, filter type,
  # and blood flow rate.' Results: 'none of the available covariates improved
  # model fit.' Tested on residual clearance, the intercompartmental
  # clearances, and the pre- and postfilter volumes.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex, 1 = female",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained. Table 1: 1 of 6 patients female."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened and not retained. Not tabulated in Table 1."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened and not retained. Table 1 median 72.5 kg (Q1 65, Q3 95)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened and not retained. Not tabulated in Table 1."
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Screened and not retained. Not tabulated in Table 1."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = paste(
        "Screened and not retained. Table 1 median 26.5 g/L (Q1 25.25, Q3",
        "28.5). The analysis used directly-measured UNBOUND concentrations."
      )
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened and not retained. Table 1 median 138 umol/L (Q1 92, Q3 151.75)."
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score",
      units = "(score)",
      type = "continuous",
      notes = "Screened and not retained. Table 1 median 11.5 (Q1 7.75, Q3 13.75)."
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation II score",
      units = "(score)",
      type = "continuous",
      notes = "Screened and not retained. Table 1 median 32 (Q1 25.25, Q3 36.5)."
    ),
    DFR = list(
      description = "Dialysate flow rate",
      units = "mL/h",
      type = "continuous",
      notes = "Screened and not retained. Table 2: 1,000 or 1,500 mL/h."
    ),
    TMP = list(
      description = "Hemofilter transmembrane pressure",
      units = "mmHg",
      type = "continuous",
      notes = "Screened and not retained. Values not reported."
    ),
    FILT_SA = list(
      description = "Hemofilter membrane surface area (filter type)",
      units = "m^2",
      type = "continuous",
      notes = paste(
        "Screened as 'filter type' and not retained. Table 2: AN69 ST100 (3",
        "patients) or ST150 (3 patients); the Table 2 footnote gives 1 and 1.5",
        "m^2, the Methods 0.9 and 1.50 m^2."
      )
    ),
    BFR = list(
      description = "Blood flow rate through the extracorporeal circuit",
      units = "mL/min",
      type = "continuous",
      notes = "Screened and not retained. Table 2: 100 to 200 mL/min."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 6L,
    n_studies = 1L,
    age_range = "23-66 years; median 61.5 (Q1 58, Q3 65)",
    weight_range = "65-103 kg; median 72.5 (Q1 65, Q3 95)",
    sex_female_pct = 16.7,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Critically ill adults in the Royal Brisbane and Women's Hospital ICU with",
      "a systemic infection known or suspected to be caused by an organism",
      "susceptible to ceftolozane-tazobactam, prescribed continuous renal",
      "replacement therapy (all six received CVVHDF). APACHE II median 32, SOFA",
      "median 11.5 (Table 1).",
      sep = " "
    ),
    renal_function = paste(
      "All on CVVHDF (Prismaflex, AN69 ST100 or ST150 filter). Table 2 settings:",
      "blood flow 100-200 mL/min, dialysate 1,000-1,500 mL/h, prefilter",
      "replacement 200-1,800 mL/h, postfilter replacement 0-1,500 mL/h, target",
      "fluid removal 30-500 mL/h. Serum creatinine median 138 umol/L.",
      sep = " "
    ),
    dose_range = paste(
      "1.5 g ceftolozane-tazobactam (2:1, i.e. 1000 mg ceftolozane / 500 mg",
      "tazobactam) every 8 h as a 1-h intravenous infusion per protocol.",
      sep = " "
    ),
    regions = "Australia (single centre, Royal Brisbane and Women's Hospital)",
    notes = paste(
      "Demographics and CVVHDF settings from Tables 1 and 2. Sampling within one",
      "dosing interval: prefilter pre-dose, 0.25, 0.75, 1.25, 2-7 and 8 h;",
      "postfilter 0.75, 2 and 6 h; effluent line 1, 2, 4, 6 and 8 h after the",
      "start of infusion. Unbound concentrations by UHPLC-MS/MS (ceftolozane",
      "calibrated 1-100 mg/L).",
      "ESTIMATION was NONPARAMETRIC (Pmetrics NPAG). Table 3 summarises the",
      "marginal support-point distributions as mean and SD; parameter",
      "correlations are not recoverable from the publication.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # STRUCTURAL PARAMETERS -- Sime 2019b Table 3, 'Ceftolozane Mean'
    # column. Each value is the MEAN of the NPAG support-point
    # distribution; the SD column is the between-subject spread (IIV), not
    # uncertainty. As in the sibling Sime_2019_ceftolozane.R (same group,
    # same software), each tabulated mean is carried as the MEDIAN of a
    # log-normal marginal; see the vignette Assumptions and deviations.
    # =====================================================================
    lcl_crrt <- log(2.659)
    label("CVVHDF clearance from the prefilter compartment, CL_CVVHDF (L/h)")
    # Table 3: CL CVVHDF mean 2.659, SD 0.783
    lcl <- log(0.596)
    label("Residual non-CVVHDF clearance from the prefilter compartment, CL_residual (L/h)")
    # Table 3: CL residual mean 0.596, SD 0.504
    lvc <- log(25.184)
    label("Prefilter (central) volume, V_pre (L)")
    # Table 3: V pre mean 25.184, SD 7.499
    lvpostfilter <- log(17.578)
    label("Postfilter compartment volume, V_post (L)")
    # Table 3: V post mean 17.578, SD 10.871
    lk12 <- log(0.43)
    label("Transfer rate constant prefilter -> postfilter, K12 (1/h)")
    # Table 3: K12 mean 0.43, SD 0.718
    lk21 <- log(0.676)
    label("Transfer rate constant postfilter -> prefilter, K21 (1/h)")
    # Table 3: K21 mean 0.676, SD 0.908
    lkd <- log(1.596)
    label("Effluent drainage rate constant, Kd (1/h)")
    # Table 3: Kd mean 1.596, SD 0.495
    lvdialysate <- log(2.178)
    label("Effluent compartment volume, V_effluent (L)")
    # Table 3: V effluent mean 2.178, SD 0.801
    lq <- log(0.834)
    label("Intercompartmental clearance prefilter <-> peripheral, Q_pre (L/h)")
    # Table 3: Q pre mean 0.834, SD 1.863
    lqpostfilter <- log(2.42)
    label("Intercompartmental clearance postfilter <-> peripheral, Q_post (L/h)")
    # Table 3: Q post mean 2.42, SD 1.451
    lvp <- log(73.379)
    label("Peripheral volume, V_peripheral (L)")
    # Table 3: V peripheral mean 73.379, SD 39.042

    # =====================================================================
    # INTER-INDIVIDUAL VARIABILITY -- Table 3 'SD' column, converted to a
    # coefficient of variation CV = SD / mean and then to a log-normal
    # variance omega^2 = log(CV^2 + 1), the convention of the sibling
    # Sime_2019_ceftolozane.R and Hughes_2024_vancomycin_nonparametric.R.
    # The NPAG joint density is not published, so the marginals are
    # encoded as independent.
    # =====================================================================
    etalcl_crrt ~ 0.08316 # CV 0.783/2.659 = 0.2945
    etalcl ~ 0.53947 # CV 0.504/0.596 = 0.8456
    etalvc ~ 0.08495 # CV 7.499/25.184 = 0.2978
    etalvpostfilter ~ 0.32387 # CV 10.871/17.578 = 0.6184
    etalk12 ~ 1.33187 # CV 0.718/0.43 = 1.6698
    etalk21 ~ 1.03111 # CV 0.908/0.676 = 1.3432
    etalkd ~ 0.09184 # CV 0.495/1.596 = 0.3102
    etalvdialysate ~ 0.12686 # CV 0.801/2.178 = 0.3678
    etalq ~ 1.79008 # CV 1.863/0.834 = 2.2338
    etalqpostfilter ~ 0.30712 # CV 1.451/2.42 = 0.5996
    etalvp ~ 0.24927 # CV 39.042/73.379 = 0.5321

    # =====================================================================
    # RESIDUAL UNEXPLAINED VARIABILITY is NOT reported. Methods,
    # 'Pharmacokinetic analysis', lists only the menu of Pmetrics error
    # models tested -- additive (SD^2 + lambda^2)^0.5, multiplicative
    # SD * gamma, and an assay polynomial C0 + C1 * obs -- without the
    # selected form or any of lambda, gamma, C0, C1; the supplement holds
    # only Figure S1 and Table S1. Carried as fixed(0), one pair per fitted
    # matrix, rather than invented.
    # =====================================================================
    propSd <- fixed(0)
    label("Proportional residual SD, prefilter plasma (fraction; not reported in the source)")
    addSd <- fixed(0)
    label("Additive residual SD, prefilter plasma (mg/L; not reported in the source)")
    propSd_Cpostfilter <- fixed(0)
    label("Proportional residual SD, postfilter plasma (fraction; not reported in the source)")
    addSd_Cpostfilter <- fixed(0)
    label("Additive residual SD, postfilter plasma (mg/L; not reported in the source)")
    propSd_Ceffluent <- fixed(0)
    label("Proportional residual SD, CVVHDF effluent (fraction; not reported in the source)")
    addSd_Ceffluent <- fixed(0)
    label("Additive residual SD, CVVHDF effluent (mg/L; not reported in the source)")
  })

  model({
    # 1. Individual parameters (Table 3 means as log-normal medians).
    cl_crrt <- exp(lcl_crrt + etalcl_crrt)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vpostfilter <- exp(lvpostfilter + etalvpostfilter)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    kd <- exp(lkd + etalkd)
    vdialysate <- exp(lvdialysate + etalvdialysate)
    q <- exp(lq + etalq)
    qpostfilter <- exp(lqpostfilter + etalqpostfilter)
    vp <- exp(lvp + etalvp)

    # 2. Concentrations in each compartment.
    Cc <- central / vc
    Cpostfilter <- postfilter / vpostfilter
    Cperipheral <- peripheral1 / vp
    Ceffluent <- circuit_dialysate / vdialysate

    # 3. Figure 1 topology. The dose enters the prefilter (central)
    #    compartment, which loses drug by CL_residual (out of the system)
    #    and CL_CVVHDF (into the effluent compartment), transfers to and
    #    from the postfilter compartment by the first-order rate constants
    #    K12 / K21, and exchanges with the peripheral compartment by Q_pre.
    #    The postfilter compartment exchanges with the same peripheral
    #    compartment by Q_post. The effluent compartment drains by Kd.
    #    The 1-h infusion is expressed in the event table.
    d/dt(central) <- -(cl + cl_crrt) * Cc - k12 * central + k21 * postfilter -
      q * (Cc - Cperipheral)
    d/dt(postfilter) <- k12 * central - k21 * postfilter -
      qpostfilter * (Cpostfilter - Cperipheral)
    d/dt(peripheral1) <- q * (Cc - Cperipheral) +
      qpostfilter * (Cpostfilter - Cperipheral)
    d/dt(circuit_dialysate) <- cl_crrt * Cc - kd * circuit_dialysate

    # 4. Observations: unbound prefilter plasma (Cc), unbound postfilter
    #    plasma and CVVHDF effluent, fitted simultaneously (Methods).
    Cc ~ add(addSd) + prop(propSd)
    Cpostfilter ~ add(addSd_Cpostfilter) + prop(propSd_Cpostfilter)
    Ceffluent ~ add(addSd_Ceffluent) + prop(propSd_Ceffluent)
  })
}
