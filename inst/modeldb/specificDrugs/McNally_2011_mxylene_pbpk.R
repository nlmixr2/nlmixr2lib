McNally_2011_mxylene_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, flow-limited, inhalation; recoded from the acslX CSL",
    "listing in the Appendix). Human m-xylene inhalation model with a",
    "bladder compartment for the urinary metabolite 3-methylhippuric acid",
    "(MHA), used by McNally et al. 2011 to demonstrate a Morris-screening +",
    "eFAST global sensitivity analysis workflow. Parent m-xylene is",
    "distributed over liver, adipose, rapidly perfused and slowly perfused",
    "tissue as flow-limited compartments; arterial blood is an algebraic",
    "steady-state lung equation coupling cardiac output to allometric",
    "pulmonary ventilation through the",
    "blood:air partition coefficient, so inhalation enters as a square-wave",
    "air concentration (ppm) rather than as a dose record. Saturable",
    "Michaelis-Menten hepatic metabolism, with Vmax extrapolated from the",
    "in vitro microsomal rate through the microsomal protein yield and liver",
    "mass. All metabolised m-xylene is converted mole-for-mole to MHA in",
    "blood, which moves first-order (K1) into a bladder that fills with urine",
    "at a constant rate; voiding is an event-table reset of the bladder MHA",
    "amount and urine volume. Urinary MHA is reported in mmol/mol creatinine",
    "through a fixed urinary creatinine concentration. Deterministic",
    "typical-value model (no IIV, no residual error): the paper evaluates",
    "the model with Table 1 mean point values against volunteer data rather",
    "than fitting it. Outputs are mixed-venous blood m-xylene (Cc, mg/L, and",
    "Cc_umol, umol/L as plotted in Figure 2A), the bladder-urine MHA",
    "concentration (Curine, mmol/mol creatinine, Figure 2B), tissue",
    "concentrations and cumulative inhaled / exhaled / metabolised amounts."
  )
  reference <- paste(
    "McNally K, Cotton R, Loizou GD. (2011). A workflow for global",
    "sensitivity analysis of PBPK models. Front Pharmacol 2:31.",
    "doi:10.3389/fphar.2011.00031. PMCID PMC3128931.",
    "Model structure adapted from Franks et al. 2006 (bladder / micturition)",
    "with parameter values similar to Loizou et al. 1999 and MacDonald et al.",
    "2002 (Table 1). Model code: Appendix, 'PBPK model CSL file'.",
    sep = " "
  )
  vignette <- "McNally_2011_mxylene"

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L",
    amount = "mg",
    weight = "kg"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Table 1 gives BW = 75 kg as the point value (range 49-92 kg, PopGen",
        "5th-95th percentiles). Carried as a covariate so the model can be",
        "re-scaled: tissue volumes are linear in body weight and cardiac",
        "output and pulmonary ventilation are allometric at exponent 0.75.",
        "Hepatic Vmax scales with liver mass (linear in body weight) through",
        "the in vitro-in vivo extrapolation.",
        sep = " "
      ),
      source_name = "BW"
    )
  )

  compartmentData <- list(
    a_rapidly_perfused = list(analyte = "m-xylene", units = "mg", specimen = "tissue", verified = TRUE),
    a_slowly_perfused = list(analyte = "m-xylene", units = "mg", specimen = "tissue", verified = TRUE),
    a_fat = list(analyte = "m-xylene", units = "mg", specimen = "tissue", verified = TRUE),
    a_liver = list(analyte = "m-xylene", units = "mg", specimen = "tissue", verified = TRUE),
    a_metabolized = list(analyte = "m-xylene", units = "mg", specimen = "not applicable", verified = TRUE),
    a_inhaled = list(analyte = "m-xylene", units = "mg", specimen = "not applicable", verified = TRUE),
    a_exhaled = list(analyte = "m-xylene", units = "mg", specimen = "not applicable", verified = TRUE),
    a_mha_blood = list(
      analyte = "3-methylhippuric acid",
      units = "mg",
      specimen = "whole blood",
      verified = TRUE
    ),
    a_mha_bladder = list(
      analyte = "3-methylhippuric acid",
      units = "mg",
      specimen = "urine",
      verified = TRUE
    ),
    urine_vol = list(analyte = "not applicable", units = "L", specimen = "urine", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 0L,
    n_studies = 0L,
    age_range = "adult (single 75 kg reference adult; no age term in the model)",
    weight_range = "75 kg point value (Table 1; sensitivity-analysis range 49-92 kg)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state = "healthy adult volunteers (occupational / environmental solvent exposure)",
    dose_range = paste(
      "Inhalation of m-xylene vapour for 4 h at a target 40 ppm (the human",
      "volunteer study); the sensitivity analysis also simulated 200 and",
      "500 ppm.",
      sep = " "
    ),
    regions = "United Kingdom (Health and Safety Laboratory, Buxton)",
    notes = paste(
      "Forward toxicokinetic PBPK simulation, not a population fit, so",
      "n_subjects = 0. The model was evaluated against venous blood m-xylene",
      "and urinary MHA from a Health and Safety Laboratory volunteer study",
      "(groups of four volunteers exposed for 4 h on two occasions to 40 ppm",
      "in an 8 m3 controlled-atmosphere facility; Figure 2 shows two",
      "volunteers). Physiological ranges came from the PopGen virtual",
      "population generator (50% male, white, 16-65 years, BMI 18.5-30);",
      "KM and in vitro Vmax from Tassaneeyakul et al. 1996; K1 was estimated",
      "by the authors from digitised Engstrom et al. 1978 urinary-excretion",
      "data.",
      sep = " "
    )
  )

  # MHA (3-methylhippuric acid) blood and bladder pools are amounts with no
  # volume, specific to this paper's empirical metabolite submodel.
  paper_specific_compartments <- c("a_mha_blood", "a_mha_bladder")

  ini({
    # =====================================================================
    # Every value is a literature constant or a Table 1 point value, so
    # every parameter is fixed(); nothing was estimated in this paper
    # except K1, which the authors fitted to external (Engstrom 1978) data
    # before using it as a fixed constant here.
    # =====================================================================

    # --- Molecular masses and ppm conversion --------------------------------
    mw <- fixed(106.17); label("Molecular mass of m-xylene (g/mol)") # Table 1, MWxyl
    mw_mha <- fixed(193.2); label("Molecular mass of 3-methylhippuric acid (g/mol)") # Table 1, MWMHA
    molar_vol <- fixed(24450); label("Ideal-gas molar volume for the ppm to mg/L conversion (mL/mol)") # Appendix CSL, IH = CIppm*binexpinh*(RMM/24450)

    # --- Cardiac output and ventilation --------------------------------------
    qc_per_kg <- fixed(12); label("Cardiac output (L/h/kg^0.75)") # Table 1, QCC
    qp_per_kg <- fixed(12); label("Respiratory rate (L/h/kg^0.75)") # Table 1, QPC
    e_wt_flow <- fixed(0.75); label("Allometric exponent on cardiac output and ventilation (unitless)") # Table 1, units 'L h-1 BW0.75'; Appendix CSL, CAE = RAE = 0.75
    # Table 1 also lists a respiratory dead space DS = 0.3. The CSL listing
    # computes QPa = QP*(1 - DS) but never uses it (gas exchange runs on the
    # undiminished ventilation), and the paper's own sensitivity results
    # give DS a main effect of 9.18e-07 (Appendix Lowry-plot data) and rank
    # it last (Table 2), so DS does not enter this model. See vignette Errata.

    # --- Blood flows, fraction of cardiac output ------------------------------
    fq_rapidly_perfused <- fixed(0.48); label("Rapidly perfused tissue blood flow, fraction of cardiac output") # Table 1, QrpdC
    fq_slowly_perfused <- fixed(0.22); label("Slowly perfused tissue blood flow, fraction of cardiac output") # Table 1, QspdC
    fq_fat <- fixed(0.05); label("Adipose blood flow, fraction of cardiac output") # Table 1, QfaC
    fq_liver <- fixed(0.25); label("Liver blood flow, fraction of cardiac output") # Table 1, QliC

    # --- Tissue masses, fraction of body mass ---------------------------------
    fv_rapidly_perfused <- fixed(0.09); label("Rapidly perfused tissue mass, fraction of body mass") # Table 1, VrpdC
    fv_slowly_perfused <- fixed(0.604); label("Slowly perfused tissue mass, fraction of body mass") # Table 1, VspdC
    fv_fat <- fixed(0.19); label("Adipose mass, fraction of body mass") # Table 1, VfaC
    fv_liver <- fixed(0.0257); label("Liver mass, fraction of body mass") # Table 1, VliC

    # --- Partition coefficients (tissue:air and blood:air) -------------------
    pc_blood_air <- fixed(19); label("Blood:air partition coefficient (unitless)") # Table 1, Pba
    pc_rapidly_perfused_air <- fixed(117); label("Rapidly perfused tissue:air partition coefficient (unitless)") # Table 1, Prpda
    pc_slowly_perfused_air <- fixed(53); label("Slowly perfused tissue:air partition coefficient (unitless)") # Table 1, Pspda
    pc_fat_air <- fixed(1874); label("Adipose:air partition coefficient (unitless)") # Table 1, Pfaa
    pc_liver_air <- fixed(279); label("Liver:air partition coefficient (unitless)") # Table 1, Plia

    # --- Hepatic metabolism (in vitro-in vivo extrapolation) -----------------
    # Table 1 prints the KM unit as 'mMol L-1'; the Appendix CSL carries the
    # same constant as Km2E1 = 1.253 mg/L = 11.8 umol/L x 106.17 g/mol, so
    # the value is micromolar. See vignette Errata.
    km_umol <- fixed(11.8); label("Michaelis constant for hepatic metabolism (umol/L)") # Table 1, KM; Appendix CSL, Km2E1 = 1.253 mg/L
    vmax_invitro <- fixed(895); label("In vitro maximum rate of metabolism (pmol/min/mg microsomal protein)") # Table 1, Vmax; Appendix CSL, VmaxCivM2E1
    mppgl <- fixed(32); label("Microsomal protein yield per gram wet weight liver (mg/g)") # Table 1, MPY

    # --- MHA bladder submodel -------------------------------------------------
    kurine_mha <- fixed(10); label("First-order rate of MHA removal from blood to bladder urine, K1 (1/h)") # Table 1, K1
    urprod <- fixed(0.07); label("Rate of urine production (L/h)") # Table 1, Rurine
    c_creat_urine <- fixed(12.5); label("Urinary creatinine concentration (mmol/L)") # Table 1, CRE

    # --- Exposure (overridden per scenario via rxSolve(params = ...)) --------
    conc_air_ppm <- fixed(40); label("Inhaled m-xylene air concentration during the exposure window (ppm)") # Methods 'The data', target 40 ppm; Appendix CSL, CIppm = 40
    inhale_start <- fixed(0); label("Start time of the inhalation exposure window (h)") # Appendix CSL, tinhon01 = 0
    inhale_time <- fixed(4); label("Duration of the inhalation exposure window (h)") # Methods 'The data', 4 h; Appendix CSL, tinhoff01 = 4

    # No IIV and no residual error: deterministic literature PBPK model
    # evaluated with Table 1 mean point values (Figure 2), not fitted.
  })

  model({
    # -------------------------------------------------------------------
    # 1. Scaled flows (L/h) and volumes (L). The Appendix CSL listing
    #    re-parameterises the fractional flows and masses (Gelman et al.
    #    1996) so that mass balance holds when they are sampled: every
    #    fractional flow is divided by the sum of the fractional flows,
    #    and the rapidly perfused, slowly perfused and adipose fractional
    #    masses are divided by the sum of all fractional masses, while the
    #    liver mass is VliC * BW without renormalisation.
    # -------------------------------------------------------------------
    qc <- qc_per_kg * WT^e_wt_flow # CSL, QC = QCC*BWc
    qp <- qp_per_kg * WT^e_wt_flow # CSL, BWr = BW**RAE; ventilation entering gas exchange (CSL QPMC)

    fq_total <- fq_rapidly_perfused + fq_slowly_perfused + fq_fat + fq_liver # CSL, qcci
    q_rapidly_perfused <- fq_rapidly_perfused / fq_total * qc
    q_slowly_perfused <- fq_slowly_perfused / fq_total * qc
    q_fat <- fq_fat / fq_total * qc
    q_liver <- fq_liver / fq_total * qc

    fv_total <- fv_rapidly_perfused + fv_slowly_perfused + fv_fat + fv_liver # CSL, vti
    v_rapidly_perfused <- fv_rapidly_perfused / fv_total * WT
    v_slowly_perfused <- fv_slowly_perfused / fv_total * WT
    v_fat <- fv_fat / fv_total * WT
    v_liver <- fv_liver * WT

    # Tissue:blood partition coefficients from tissue:air / blood:air
    pc_rapidly_perfused <- pc_rapidly_perfused_air / pc_blood_air
    pc_slowly_perfused <- pc_slowly_perfused_air / pc_blood_air
    pc_fat <- pc_fat_air / pc_blood_air
    pc_liver <- pc_liver_air / pc_blood_air

    # In vitro-in vivo extrapolation: pmol/min/mg protein x mg protein/g
    # liver x g liver (1000 * v_liver) x 60 min/h x 1e-12 mol/pmol x
    # mw g/mol x 1000 mg/g = mg/h. CSL, Vmax2E1 = VmaxCivM2E1*1E-12*Vli*RMM*MPY*1000000*60
    vmax <- vmax_invitro * 1e-12 * v_liver * mw * mppgl * 1e6 * 60
    km <- km_umol * mw / 1000 # umol/L -> mg/L

    # -------------------------------------------------------------------
    # 2. Inhalation exposure: a single square wave from inhale_start for
    #    inhale_time hours. CSL, IH = CIppm * binexpinh * (RMM/24450).
    # -------------------------------------------------------------------
    air_on <- (time >= inhale_start) * (time < inhale_start + inhale_time)
    ci <- conc_air_ppm * mw / molar_vol * air_on

    # -------------------------------------------------------------------
    # 3. Venous-effluent concentrations (mg/L), mixed venous blood and the
    #    algebraic steady-state lung (arterial) equation.
    # -------------------------------------------------------------------
    c_rapidly_perfused <- a_rapidly_perfused / v_rapidly_perfused
    c_slowly_perfused <- a_slowly_perfused / v_slowly_perfused
    c_fat <- a_fat / v_fat
    c_liver <- a_liver / v_liver
    cv_rapidly_perfused <- c_rapidly_perfused / pc_rapidly_perfused
    cv_slowly_perfused <- c_slowly_perfused / pc_slowly_perfused
    cv_fat <- c_fat / pc_fat
    cv_liver <- c_liver / pc_liver

    c_venous <- (q_slowly_perfused * cv_slowly_perfused +
      q_rapidly_perfused * cv_rapidly_perfused +
      q_fat * cv_fat + q_liver * cv_liver) / qc # CSL, CV
    c_arterial <- (qc * c_venous + qp * ci) / (qc + qp / pc_blood_air) # CSL, CA = ((QC*CV) + (QPMC*IH))/(QC + (QPMC/Pba))
    c_alveolar <- c_arterial / pc_blood_air # CSL, CX = CA/Pba

    r_metab <- vmax * cv_liver / (km + cv_liver) # CSL, MRli; paper Methods, MRLi

    # -------------------------------------------------------------------
    # 4. Parent m-xylene mass balance (CSL DERIVATIVE block).
    # -------------------------------------------------------------------
    d/dt(a_rapidly_perfused) <- q_rapidly_perfused * (c_arterial - cv_rapidly_perfused)
    d/dt(a_slowly_perfused) <- q_slowly_perfused * (c_arterial - cv_slowly_perfused)
    d/dt(a_fat) <- q_fat * (c_arterial - cv_fat)
    d/dt(a_liver) <- q_liver * (c_arterial - cv_liver) - r_metab
    d/dt(a_metabolized) <- r_metab
    d/dt(a_inhaled) <- qp * ci # CSL, Rdose = QPMC*IH
    d/dt(a_exhaled) <- qp * c_alveolar # CSL, RAX = QPMC*CX

    # -------------------------------------------------------------------
    # 5. MHA submodel (paper Methods equations; CSL RAMHAB / RAMHAU).
    #    The bladder fills with urine at urprod; voiding is an event-table
    #    reset (evid = 5, amt = 0) of a_mha_bladder and urine_vol.
    # -------------------------------------------------------------------
    d/dt(a_mha_blood) <- r_metab * mw_mha / mw - kurine_mha * a_mha_blood
    d/dt(a_mha_bladder) <- kurine_mha * a_mha_blood
    d/dt(urine_vol) <- urprod

    # -------------------------------------------------------------------
    # 6. Outputs. Cc is the mixed-venous blood m-xylene concentration (CV in
    #    the paper). Curine is the MHA concentration in bladder urine in
    #    mmol/mol creatinine: mg / (g/mol) = mmol, x 1000 = umol, divided by
    #    L urine x mmol/L creatinine = umol/mmol = mmol/mol. Observed just
    #    before a void it equals the concentration in the voided urine
    #    (CSL, CUrine = AUCREmmol/(VBladder*CREmmol)).
    # -------------------------------------------------------------------
    Cc <- c_venous
    Cc_umol <- c_venous / mw * 1000 # CSL, CVumol
    Cexhaled <- c_alveolar
    Cliver <- c_liver
    Cfat <- c_fat
    Crapid <- c_rapidly_perfused
    Cslow <- c_slowly_perfused
    mha_umol_bladder <- a_mha_bladder * 1000 / mw_mha # CSL, AUCREmmol
    Curine <- 0
    if (urine_vol > 0) {
      Curine <- mha_umol_bladder / (urine_vol * c_creat_urine)
    }
  })
}
