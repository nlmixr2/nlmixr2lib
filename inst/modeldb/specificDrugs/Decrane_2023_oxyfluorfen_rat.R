Decrane_2023_oxyfluorfen_rat <- function() {
  description <- "QSP (PBPK + thyroid-hormone kinetics). Preclinical (rat, Sprague-Dawley). Whole-body PBPK model for the diphenyl-ether herbicide oxyfluorfen coupled to a hypothalamic-pituitary-thyroid (HPT) axis thyroid-hormone kinetics sub-model (Decrane et al. 2023, Curr Res Toxicol 5:100138). Oxyfluorfen is distributed across 13 flow-limited compartments (gut lumen, GI tissue, liver, kidney, muscle, skin, fat, brain, slowly perfused, rapidly perfused, a diffusion-limited thyroid split into thyroid blood and thyroid tissue, and a single well-mixed blood pool) with hepatic metabolic clearance and GFR-driven renal excretion. The predicted oxyfluorfen concentration in thyroid blood drives an empirical in-vitro-derived sodium-iodide-symporter (NIS) inhibition function that scales thyroidal T4 and T3 synthesis (the paper's CBthy; see the model file and vignette Errata for the evidence that CBthy is the thyroid blood rather than the thyroid tissue concentration, a reading that reproduces both species' published points of departure); five thyroid-hormone states (T4 and T3 in thyroid tissue, T4, T3 and TSH in serum) carry the HPT feedback loop in which TSH stimulates T4 synthesis (stim) and serum T4 regulates TSH production and turnover (feed1, feed2). The model is deterministic typical-value (no IIV, no residual error) -- the paper reports point estimates only, with no parameter uncertainty, no between-animal variability and no measurement-error model. Companion model: Decrane_2023_oxyfluorfen_human."
  reference <- paste(
    "Decrane R, Stoker T, Murr A, Ford J, El-Masri H. (2023).",
    "Cross species extrapolation of the disruption of thyroid hormone synthesis by oxyfluorfen",
    "using in vitro data, physiologically based pharmacokinetic (PBPK), and thyroid hormone kinetics models.",
    "Current Research in Toxicology 5:100138.",
    "doi:10.1016/j.crtox.2023.100138. PMCID PMC10697989.",
    sep = " "
  )
  vignette <- "Decrane_2023_oxyfluorfen"

  # Non-canonical states are all registered canonicals as of this PR (see
  # inst/references/compartment-names.md); no paper_specific_compartments needed.

  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L",
    amount        = "mg",
    weight        = "kg"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      notes              = paste(
        "Scales cardiac output (allometric exponent 0.75), every tissue volume (linear),",
        "hepatic clearance (exponent 0.75), GFR (linear), the thyroidal and serum",
        "thyroid-hormone production / transfer rate constants (exponent 0.66) and the",
        "thyroid-hormone volumes of distribution (linear).",
        "Table 1 of the paper reports rat body weight as 'Varies' and never states a",
        "numeric value: the 8-day gavage study used growing adolescent Sprague-Dawley rats",
        "and the body-weight trajectory lives in the co-submitted Stoker et al. manuscript,",
        "which is not part of this paper. The paper further distinguishes an instantaneous BW",
        "(used to adjust the administered dose as the animals grew) from an 'Avg.BW' /",
        "terminal BW (used for CLR, GFR and the hormone volumes of distribution) but reports",
        "neither. This single WT covariate therefore collapses both; supply it per subject",
        "via the event table or via rxSolve(params = c(WT = ...)). See the vignette Errata.",
        sep = " "
      ),
      source_name        = "BW"
    )
  )

  compartmentData <- list(
    # --- oxyfluorfen (chemical PBPK layer), amounts in mg ---
    stomach            = list(analyte = "Oxyfluorfen", units = "mg", specimen = "administration site", verified = TRUE),
    a_gut              = list(analyte = "Oxyfluorfen", units = "mg", specimen = "tissue", verified = TRUE),
    a_liver            = list(analyte = "Oxyfluorfen", units = "mg", specimen = "tissue", verified = TRUE),
    a_kidney           = list(analyte = "Oxyfluorfen", units = "mg", specimen = "tissue", verified = TRUE),
    a_muscle           = list(analyte = "Oxyfluorfen", units = "mg", specimen = "tissue", verified = TRUE),
    a_skin             = list(analyte = "Oxyfluorfen", units = "mg", specimen = "tissue", verified = TRUE),
    a_fat              = list(analyte = "Oxyfluorfen", units = "mg", specimen = "tissue", verified = TRUE),
    a_brain            = list(analyte = "Oxyfluorfen", units = "mg", specimen = "tissue", verified = TRUE),
    a_slowly_perfused  = list(analyte = "Oxyfluorfen", units = "mg", specimen = "tissue", verified = TRUE),
    a_rapidly_perfused = list(analyte = "Oxyfluorfen", units = "mg", specimen = "tissue", verified = TRUE),
    a_thyroid_blood    = list(analyte = "Oxyfluorfen", units = "mg", specimen = "whole blood", verified = TRUE),
    a_thyroid_tissue   = list(analyte = "Oxyfluorfen", units = "mg", specimen = "tissue", verified = TRUE),
    a_blood            = list(analyte = "Oxyfluorfen", units = "mg", specimen = "whole blood", verified = TRUE),
    # --- thyroid hormones (TH kinetics layer), amounts in mg ---
    t4_thyroid         = list(analyte = "Thyroxine (T4)", units = "mg", specimen = "tissue", verified = TRUE),
    t3_thyroid         = list(analyte = "Triiodothyronine (T3)", units = "mg", specimen = "tissue", verified = TRUE),
    t4_serum           = list(analyte = "Thyroxine (T4)", units = "mg", specimen = "serum", verified = TRUE),
    t3_serum           = list(analyte = "Triiodothyronine (T3)", units = "mg", specimen = "serum", verified = TRUE),
    tsh_serum          = list(analyte = "Thyrotropin (TSH)", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species        = "rat (Sprague-Dawley)",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    age_range      = "adolescent (exact age not reported)",
    weight_range   = "not reported; Table 1 lists body weight as 'Varies'",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = "healthy (untreated controls plus oxyfluorfen-exposed groups)",
    dose_range     = "0.8125, 1.625, 3.25, 7.5, 15, 31.25 and 62.5 mg/kg once daily by oral gavage for 8 days (1% methyl cellulose suspension)",
    regions        = "United States (US EPA ORD)",
    notes          = paste(
      "Calibration data (oxyfluorfen in thyroid tissue and serum; serum T4 and T3) were",
      "'provided via personal communication and illustrated in a co-submitted paper'",
      "(Stoker et al.), so the number of animals per dose group, their age and their body",
      "weights are not reported in this paper. Basal thyroid-hormone levels came from the",
      "experimental controls plus Handa et al. (2021) and Hassan et al. (2020);",
      "physiological parameters from Brown et al. (1997) (Table 1); tissue:blood partition",
      "coefficients from GastroPlus 9.8 (Table 3).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Rat physiological parameters -- Table 1 (all from Brown et al. 1997).
    # Every value here is a literature physiological constant, so all are
    # fixed().
    # =====================================================================

    # Cardiac output: Table 1 gives "0.235*BW.75*60" with units L/hr/kg^.75.
    # The 0.235 is L/min/kg^0.75 and the *60 converts to L/hr; both factors
    # are carried in model() so the printed form stays recognisable.
    qc_coef            <- fixed(0.235)     ; label("Cardiac-output allometric coefficient (L/min/kg^0.75)")   # Table 1, Cardiac Output = 0.235*BW^.75*60

    # Tissue volumes, fraction of body weight.
    # NB Table 1 labels the blood pool "Plasma" = 0.074, but 0.074 is Brown
    # et al. (1997)'s rat TOTAL BLOOD volume fraction and the PBPK equations
    # use it as the volume of the well-mixed blood pool that supplies Cart.
    # Carried verbatim; see the vignette Errata.
    fv_plasma          <- fixed(0.074)     ; label("Blood-pool volume, fraction of body weight (Table 1 'Plasma')")  # Table 1
    fv_thyroid_total   <- fixed(0.00005)   ; label("Total thyroid volume, fraction of body weight")            # Table 1, Thyroid Total
    fv_thyroid_blood   <- fixed(0.18)      ; label("Thyroid blood sub-compartment, fraction of TOTAL THYROID volume")  # Table 1, Thyroid Blood
    fv_fat             <- fixed(0.070)     ; label("Fat volume, fraction of body weight")                      # Table 1
    fv_muscle          <- fixed(0.4043)    ; label("Muscle volume, fraction of body weight")                   # Table 1
    fv_skin            <- fixed(0.1903)    ; label("Skin volume, fraction of body weight")                     # Table 1
    fv_liver           <- fixed(0.05)      ; label("Liver volume, fraction of body weight")                    # Table 1
    fv_gut             <- fixed(0.026925)  ; label("GI-tract tissue volume, fraction of body weight")          # Table 1, GI tract
    fv_kidney          <- fixed(0.0073)    ; label("Kidney volume, fraction of body weight")                   # Table 1
    fv_brain           <- fixed(0.0057)    ; label("Brain volume, fraction of body weight")                    # Table 1
    fv_slowly_perfused <- fixed(0.073)     ; label("Slowly-perfused tissue volume, fraction of body weight")   # Table 1 (footnote 2: assumed mainly bone)
    fv_rapidly_perfused <- fixed(0.122)    ; label("Rapidly-perfused tissue volume, fraction of body weight")  # Table 1 (footnote 3: remaining body weight)

    # Blood flows, fraction of cardiac output.
    # Table 1's "Liver 0.174" is the TOTAL hepatic venous outflow and equals
    # Portal Vein (0.153) + Hepatic Artery (0.021) exactly -- an internal
    # consistency check the table passes. model() builds the liver flows from
    # the portal-vein and hepatic-artery rows, so fq_liver is declared for
    # traceability and used only in that identity check.
    fq_thyroid         <- fixed(0.0027)    ; label("Thyroid blood flow, fraction of cardiac output")           # Table 1
    fq_fat             <- fixed(0.070)     ; label("Fat blood flow, fraction of cardiac output")               # Table 1
    fq_muscle          <- fixed(0.278)     ; label("Muscle blood flow, fraction of cardiac output")            # Table 1
    fq_skin            <- fixed(0.058)     ; label("Skin blood flow, fraction of cardiac output")              # Table 1
    fq_portal_vein     <- fixed(0.153)     ; label("Portal-vein blood flow, fraction of cardiac output")       # Table 1, Portal Vein
    fq_hepatic_artery  <- fixed(0.021)     ; label("Hepatic-artery blood flow, fraction of cardiac output")    # Table 1, Hepatic Artery
    fq_kidney          <- fixed(0.141)     ; label("Kidney blood flow, fraction of cardiac output")            # Table 1
    fq_brain           <- fixed(0.020)     ; label("Brain blood flow, fraction of cardiac output")             # Table 1
    fq_slowly_perfused <- fixed(0.0823)    ; label("Slowly-perfused tissue blood flow, fraction of cardiac output")   # Table 1
    fq_rapidly_perfused <- fixed(0.3453)   ; label("Rapidly-perfused tissue blood flow, fraction of cardiac output")  # Table 1

    # =====================================================================
    # Tissue:blood partition coefficients -- Table 3 (GastroPlus 9.8).
    # Identical for rat and human. All fixed (software-predicted).
    # =====================================================================
    pc_thyroid         <- fixed(8.5)       ; label("Thyroid tissue:blood partition coefficient (unitless)")    # Table 3, Thyroid
    pc_muscle          <- fixed(9.5)       ; label("Muscle tissue:blood partition coefficient (unitless)")     # Table 3
    pc_skin            <- fixed(10.86)     ; label("Skin tissue:blood partition coefficient (unitless)")       # Table 3
    pc_fat             <- fixed(50.4)      ; label("Fat tissue:blood partition coefficient (unitless)")        # Table 3
    pc_liver           <- fixed(9.5)       ; label("Liver tissue:blood partition coefficient (unitless)")      # Table 3
    pc_kidney          <- fixed(9.5)       ; label("Kidney tissue:blood partition coefficient (unitless)")     # Table 3
    pc_brain           <- fixed(20.31)     ; label("Brain tissue:blood partition coefficient (unitless)")      # Table 3
    pc_slowly_perfused <- fixed(4.71)      ; label("Slowly-perfused tissue:blood partition coefficient (unitless)")   # Table 3, Slowly Perfused
    pc_rapidly_perfused <- fixed(9.7)      ; label("Rapidly-perfused tissue:blood partition coefficient (unitless)")  # Table 3, Rapidly Perfused
    pc_gut             <- fixed(9.5)       ; label("GI tissue:blood partition coefficient (unitless)")         # Table 3, GI

    # =====================================================================
    # Oxyfluorfen biochemical parameters -- Table 4.
    # Table 4 footnote 5 = "Fit to rat data"; those three are the only
    # ESTIMATED parameters in the rat chemical model, so they are NOT
    # wrapped in fixed(). frac and GFR are literature/software values.
    # =====================================================================
    pa_thyroid         <- 0.005            ; label("Thyroid blood-to-tissue permeability-area product (L/h)")  # Table 4, PAThy rat (footnote 5: fit to rat data)
    clr_coef           <- 1.275            ; label("Hepatic metabolic clearance allometric coefficient (L/h/kg^0.75)")  # Table 4, CLR rat = 1.275*(Avg.BW)^0.75 (footnote 5: fit to rat data)
    lka                <- log(0.00375)     ; label("Log first-order oral absorption rate constant from gut lumen (1/h)")  # Table 4, Ka rat = 0.00375 (footnote 5: fit to rat data)
    frac               <- fixed(0.08)      ; label("Fraction of oxyfluorfen unbound in blood (unitless)")      # Table 4, frac (footnote 2: GastroPlus 9.8)
    gfr_coef           <- fixed(0.60)      ; label("Glomerular filtration rate per body weight (L/h/kg)")      # Table 4, GFR rat = 60e-2*BW (footnote 3: Carrara et al. 2016)
    e_wt_cl            <- fixed(0.75)      ; label("Allometric exponent on cardiac output and hepatic clearance (unitless)")  # Table 1 (BW^.75) and Table 4 (Avg.BW^0.75)

    # =====================================================================
    # Empirical in-vitro NIS-inhibition function -- paper p. 5:
    #   Inhib_rat = 0.788*exp(-2.006*CBthy) + 0.2218
    # fitted to the Buckalew et al. (2020) FRTL-5 rat radioactive-iodide-
    # uptake data shown in Fig. 3a. CBthy is the oxyfluorfen concentration in
    # thyroid tissue (mg/L). At CBthy = 0 the function returns 1.0098, i.e.
    # it is normalised to the drug-free synthesis rate (see vignette).
    # =====================================================================
    inhib_amp          <- fixed(0.788)     ; label("NIS-inhibition exponential amplitude (unitless)")          # p. 5, Inhib_rat
    inhib_rate         <- fixed(-2.006)    ; label("NIS-inhibition exponential rate on thyroid oxyfluorfen concentration (L/mg)")  # p. 5, Inhib_rat
    inhib_floor        <- fixed(0.2218)    ; label("NIS-inhibition asymptotic residual synthesis fraction (unitless)")             # p. 5, Inhib_rat

    # =====================================================================
    # Thyroid-hormone kinetics -- Table 5.
    # Table 5's footnote 1 is NOT defined anywhere in the paper. From the
    # Methods ("Whenever applicable, parameters were estimated by fitting
    # model simulations to data ...") and by analogy with Table 4 footnote 5,
    # footnote 1 marks parameters CALIBRATED to the rat data; those are left
    # un-fixed(). Footnotes 3, 5, 7 and 8 are literature citations -> fixed().
    # See the vignette Errata for this inference.
    # =====================================================================
    fr_t4_t3           <- fixed(0.22)      ; label("Fraction of T4 that converts to T3 (unitless)")            # Table 5, fr (footnote 8: Silva and Matthews 1984)

    # Thyroidal synthesis and interconversion
    prod_t4_coef       <- 3.68e-4          ; label("Thyroidal T4 zero-order production coefficient (mg/h/kg^0.66)")  # Table 5, T4prod rat = 3.68e-4*rat_BW^0.66 (footnote 1: calibrated)
    prod_t3_coef       <- 3.43e-5          ; label("Thyroidal T3 zero-order production coefficient (mg/h/kg^0.66)")  # Table 5, T3prod rat = 3.43e-5*rat_BW^0.66 (footnote 1)
    met_t4thy          <- 0.065            ; label("Thyroidal T4-to-T3 conversion rate constant (1/h)")        # Table 5, T4met rat (footnote 1)
    met_t3thy          <- 0.045            ; label("Thyroidal T3 metabolism rate constant (1/h)")              # Table 5, T3met rat (footnote 1)
    e_wt_th            <- fixed(0.66)      ; label("Allometric exponent on thyroid-hormone production and transfer rates (unitless)")  # Table 5, rat_BW^0.66 exponents

    # Basal (homeostatic) levels used to set the initial conditions and to
    # normalise the HPT feedback terms.
    t4_bl_thy          <- fixed(208)       ; label("Basal T4 concentration in thyroid tissue (mg/L)")          # Table 5, T4_BL rat (footnote 3: Handa et al. 2021)
    t3_bl_thy          <- fixed(18)        ; label("Basal T3 concentration in thyroid tissue (mg/L)")          # Table 5, T3_BL rat (footnote 3)
    tsh_bl             <- fixed(2.5e-3)    ; label("Basal TSH concentration in serum (mg/L)")                 # Table 5, TSH_BL rat (footnote 3)
    t4_bl_srm          <- fixed(58.4e-3)   ; label("Basal T4 concentration in serum (mg/L)")                   # Table 5, T4_BL_srm rat (no footnote; from the experimental controls)
    t3_bl_srm          <- fixed(1.26e-3)   ; label("Basal T3 concentration in serum (mg/L)")                   # Table 5, T3_BL_srm rat (no footnote; from the experimental controls)

    # TSH turnover and the HPT power-law feedback exponents
    k_tsh              <- fixed(5e-6)      ; label("First-order TSH loss rate constant (1/h)")                 # Table 5, k_TSH (footnote 3: Handa et al. 2021)
    nf1                <- fixed(2.53)      ; label("Feedback exponent NF1 on serum T4 driving TSH production (unitless)")  # Table 5, NF1 (footnote 5: Ekerot et al. 2013)
    nf2                <- fixed(1.9)       ; label("Feedback exponent NF2 on serum T4 modulating TSH turnover (unitless)")  # Table 5, NF2 (footnote 5)
    nf3                <- fixed(0.11)      ; label("Stimulation exponent NF3 on serum TSH driving thyroidal T4 synthesis (unitless)")  # Table 5, NF3 (footnote 5)

    # Thyroid-to-serum transfer and systemic hormone disposition
    f_t4srm_coef       <- 0.13             ; label("Thyroid-to-serum T4 transfer coefficient (1/h/kg^0.66)")   # Table 5, fT4 rat = 0.13*rat_BW^0.66 (footnote 1)
    f_t3srm_coef       <- 0.04             ; label("Thyroid-to-serum T3 transfer coefficient (1/h/kg^0.66)")   # Table 5, fT3 rat = 0.04*rat_BW^0.66 (footnote 1)
    met_t4srm          <- 0.085            ; label("Systemic T4-to-T3 conversion rate constant in serum (1/h)")  # Table 5, T4met_srm (footnote 1)
    met_t3srm          <- 0.045            ; label("Systemic T3 loss rate constant in serum (1/h)")            # Table 5, T3met_srm (footnote 1)
    loss_srm           <- 0.045            ; label("Systemic non-deiodinative T4 loss rate constant in serum (1/h)")  # Table 5, T4loss_srm (footnote 1)

    # Thyroid-hormone volumes of distribution (L), linear in body weight.
    f_vd_t4            <- fixed(0.149)     ; label("T4 volume of distribution per body weight (L/kg)")         # Table 5, VDT4 = 0.149*Avg.BW (footnote 7: Dubois and Dussault 1977)
    f_vd_t3            <- fixed(1.62)      ; label("T3 volume of distribution per body weight (L/kg)")         # Table 5, VDT3 = 1.62*Avg.BW (footnote 7)
    f_vd_tsh           <- fixed(0.149)     ; label("TSH volume of distribution per body weight (L/kg)")        # Table 5, VDTSH = 0.149*Avg.BW (footnote 7)

    # =====================================================================
    # NO IIV (eta) and NO residual-error block. The paper reports single
    # point estimates for every parameter with no standard errors, no
    # between-animal variance components and no measurement-error model, so
    # this file is a deterministic typical-value mechanistic simulator.
    # Inventing variances would be fabrication. See the vignette Errata.
    # =====================================================================
  })

  model({
    # -------------------------------------------------------------------
    # 1. Cardiac output and regional blood flows (L/h).
    #    Table 1: cardiac output = 0.235 * BW^0.75 * 60 (L/min/kg^0.75 -> L/h).
    #    Each regional flow is its Table 1 fraction of cardiac output.
    # -------------------------------------------------------------------
    qc                 <- qc_coef * WT^e_wt_cl * 60

    q_thyroid          <- fq_thyroid          * qc
    q_fat              <- fq_fat              * qc
    q_muscle           <- fq_muscle           * qc
    q_skin             <- fq_skin             * qc
    q_portal_vein      <- fq_portal_vein      * qc
    q_hepatic_artery   <- fq_hepatic_artery   * qc
    q_kidney           <- fq_kidney           * qc
    q_brain            <- fq_brain            * qc
    q_slowly_perfused  <- fq_slowly_perfused  * qc
    q_rapidly_perfused <- fq_rapidly_perfused * qc

    # Qvenliver: "the venous blood leaving the liver is equal to sum of
    # portal vein and arterial blood flow to the tissue" (paper p. 3).
    q_ven_liver        <- q_portal_vein + q_hepatic_artery

    # Total flow drawn from the blood pool. Implemented as the ACTUAL sum of
    # the Table 1 flows rather than as qc, because the rat Table 1 fractions
    # sum to 1.171 of cardiac output rather than 1.000 (the human Table 2
    # fractions sum to 0.956). Using the actual sum keeps the circulatory
    # loop exactly mass-conserving by construction; using qc would not.
    # See the vignette Errata.
    q_total            <- q_ven_liver + q_kidney + q_muscle + q_skin + q_fat +
                          q_brain + q_slowly_perfused + q_rapidly_perfused + q_thyroid

    # -------------------------------------------------------------------
    # 2. Compartment volumes (L). Table 1 fractions of body weight.
    #    The thyroid is split into a vascular and a tissue sub-compartment:
    #    Table 1's "Thyroid Blood 0.18" is a fraction of the TOTAL thyroid
    #    volume (0.00005 of body weight), not of body weight.
    # -------------------------------------------------------------------
    v_blood            <- fv_plasma           * WT
    v_thyroid_total    <- fv_thyroid_total    * WT
    v_thyroid_blood    <- fv_thyroid_blood    * v_thyroid_total
    v_thyroid_tissue   <- v_thyroid_total - v_thyroid_blood
    v_fat              <- fv_fat              * WT
    v_muscle           <- fv_muscle           * WT
    v_skin             <- fv_skin             * WT
    v_liver            <- fv_liver            * WT
    v_gut              <- fv_gut              * WT
    v_kidney           <- fv_kidney           * WT
    v_brain            <- fv_brain            * WT
    v_slowly_perfused  <- fv_slowly_perfused  * WT
    v_rapidly_perfused <- fv_rapidly_perfused * WT

    # -------------------------------------------------------------------
    # 3. Hepatic clearance and renal filtration (L/h).
    # -------------------------------------------------------------------
    ka                 <- exp(lka)
    clr                <- clr_coef * WT^e_wt_cl
    gfr                <- gfr_coef * WT

    # -------------------------------------------------------------------
    # 4. Tissue concentrations (mg/L). c_<tissue> is the tissue
    #    concentration; cv_<tissue> = c_<tissue> / pc_<tissue> is the
    #    concentration in the blood leaving that tissue (flow-limited).
    #    a_thyroid_blood is BLOOD, so it carries no partition coefficient.
    # -------------------------------------------------------------------
    c_art               <- a_blood / v_blood

    c_gut               <- a_gut / v_gut
    cv_gut              <- c_gut / pc_gut
    c_liver             <- a_liver / v_liver
    cv_liver            <- c_liver / pc_liver
    c_kidney            <- a_kidney / v_kidney
    cv_kidney           <- c_kidney / pc_kidney
    c_muscle            <- a_muscle / v_muscle
    cv_muscle           <- c_muscle / pc_muscle
    c_skin              <- a_skin / v_skin
    cv_skin             <- c_skin / pc_skin
    c_fat               <- a_fat / v_fat
    cv_fat              <- c_fat / pc_fat
    c_brain             <- a_brain / v_brain
    cv_brain            <- c_brain / pc_brain
    c_slowly_perfused   <- a_slowly_perfused / v_slowly_perfused
    cv_slowly_perfused  <- c_slowly_perfused / pc_slowly_perfused
    c_rapidly_perfused  <- a_rapidly_perfused / v_rapidly_perfused
    cv_rapidly_perfused <- c_rapidly_perfused / pc_rapidly_perfused

    c_thyroid_blood     <- a_thyroid_blood / v_thyroid_blood
    c_thyroid_tissue    <- a_thyroid_tissue / v_thyroid_tissue
    cv_thyroid_tissue   <- c_thyroid_tissue / pc_thyroid

    # -------------------------------------------------------------------
    # 5. Chemical PBPK ODE system (13 states), paper pp. 2-3.
    # -------------------------------------------------------------------
    # dAstomach/dt = -ka*Astomach
    d/dt(stomach)            <- -ka * stomach

    # dAGI/dt = Qportalvein*(Cart - CGI/pcGI) + ka*Astomach
    d/dt(a_gut)              <- q_portal_vein * (c_art - cv_gut) + ka * stomach

    # dAliver/dt = Qartliver*Cart - Qvenliver*Cliver/pcliver
    #              + Qportalvein*CGI/pcGI - CLR*Aliver/(Vliver*pcliver)
    # The last term equals CLR*cv_liver by definition of cv_liver.
    d/dt(a_liver)            <- q_hepatic_artery * c_art - q_ven_liver * cv_liver +
                                q_portal_vein * cv_gut - clr * cv_liver

    # dAkidney/dt = Qkidney*(Cart - Ckidney/pckidney) - GFR*frac*Ckidney
    # NB the renal sink uses the TISSUE concentration Ckidney, not the
    # blood-side cv_kidney -- carried verbatim from the printed equation.
    d/dt(a_kidney)           <- q_kidney * (c_art - cv_kidney) - gfr * frac * c_kidney

    # dAx/dt = Qx*(Cart - Cx/pcx) for the flow-limited tissues
    d/dt(a_muscle)           <- q_muscle            * (c_art - cv_muscle)
    d/dt(a_skin)             <- q_skin              * (c_art - cv_skin)
    d/dt(a_fat)              <- q_fat               * (c_art - cv_fat)
    d/dt(a_brain)            <- q_brain             * (c_art - cv_brain)
    d/dt(a_slowly_perfused)  <- q_slowly_perfused   * (c_art - cv_slowly_perfused)
    d/dt(a_rapidly_perfused) <- q_rapidly_perfused  * (c_art - cv_rapidly_perfused)

    # Diffusion-limited thyroid, split into blood and tissue sub-compartments.
    # dAthyblood/dt  = Qthy*(Cart - Cthyblood) + PAthy*(Cthytissue/pcthy - Cthyblood)
    # dAthytissue/dt = -PAthy*(Cthytissue/pcthy - Cthyblood)
    d/dt(a_thyroid_blood)    <- q_thyroid * (c_art - c_thyroid_blood) +
                                pa_thyroid * (cv_thyroid_tissue - c_thyroid_blood)
    d/dt(a_thyroid_tissue)   <- -pa_thyroid * (cv_thyroid_tissue - c_thyroid_blood)

    # Single well-mixed blood pool. This is the ONE chemical ODE the paper
    # does not print. It is reconstructed by mass-balance closure of the ten
    # printed tissue equations: inflow = every tissue's venous outflow,
    # outflow = the total arterial flow times Cart. No modelling choice is
    # involved -- summing all 13 chemical ODEs then leaves exactly the two
    # elimination terms (hepatic clr*cv_liver and renal gfr*frac*c_kidney),
    # which is the mass-balance gate exercised in the vignette.
    # See the vignette Errata.
    d/dt(a_blood)            <- q_ven_liver * cv_liver +
                                q_kidney            * cv_kidney +
                                q_muscle            * cv_muscle +
                                q_skin              * cv_skin +
                                q_fat               * cv_fat +
                                q_brain             * cv_brain +
                                q_slowly_perfused   * cv_slowly_perfused +
                                q_rapidly_perfused  * cv_rapidly_perfused +
                                q_thyroid           * c_thyroid_blood -
                                q_total * c_art

    # -------------------------------------------------------------------
    # 6. NIS inhibition (paper p. 5, Fig. 3a). Inhib multiplies thyroidal T4
    #    and T3 synthesis; it equals 1.0098 at CBthy = 0 (normalised to the
    #    drug-free rate) and decays to the residual floor 0.2218.
    #
    #    CBthy IS THE THYROID *BLOOD* CONCENTRATION, not the thyroid tissue
    #    concentration. The paper's prose says "CBthy is the concentration of
    #    chemical in the thyroid cells (in vitro) or tissue (in vivo)", but
    #    four independent lines of evidence say otherwise, and the choice
    #    changes the predicted dose-response by exactly pc_thyroid = 8.5x:
    #
    #     1. The symbol. The ODEs use Cthytissue and Cthyblood; CBthy is a
    #        distinct third symbol whose "B" denotes blood.
    #     2. The mechanism. The paper's own Introduction states NIS "resides
    #        in the basolateral membrane of thyroid epithelial cells and
    #        simultaneously transports two Na+ and one I- from extracellular
    #        fluid (plasma) into the thyroid epithelial cell" -- the
    #        inhibitor concentration NIS is exposed to is the EXTRACELLULAR
    #        one.
    #     3. In-vitro correspondence. The Buckalew et al. (2020) FRTL-5 assay
    #        dosed cells via the culture medium, so the in-vivo analogue of
    #        the in-vitro exposure metric is thyroid blood, not intracellular
    #        tissue.
    #     4. Quantitative. Both of the paper's independently reported points
    #        of departure reproduce under this reading and miss by ~8.5x
    #        under the other:
    #          rat   7.5 mg/kg/day x 8 d -> 10% serum-T4 drop (Discussion p.7)
    #                blood reading: 9.2%   tissue reading: 51.0%
    #          human 57 mg/L drinking water -> 10% serum-T4 drop (Fig. 8a)
    #                blood reading: 10.0%  tissue reading: 18.0%
    #
    #    Flagged as the primary review point in the PR and quantified
    #    side-by-side in the vignette Errata. To switch to the literal prose
    #    reading, replace c_thyroid_blood with c_thyroid_tissue below.
    # -------------------------------------------------------------------
    inhib <- inhib_amp * exp(inhib_rate * c_thyroid_blood) + inhib_floor

    # -------------------------------------------------------------------
    # 7. Thyroid-hormone kinetics (paper pp. 3-4).
    #    Body-weight-scaled production and transfer rate constants.
    # -------------------------------------------------------------------
    prod_t4 <- prod_t4_coef   * WT^e_wt_th
    prod_t3 <- prod_t3_coef   * WT^e_wt_th
    f_t4srm <- f_t4srm_coef   * WT^e_wt_th
    f_t3srm <- f_t3srm_coef   * WT^e_wt_th

    # Hormone volumes of distribution (L), linear in body weight (Table 5).
    vd_t4  <- f_vd_t4  * WT
    vd_t3  <- f_vd_t3  * WT
    vd_tsh <- f_vd_tsh * WT

    # HPT feedback terms.
    #   stim  = (CTSH/TSHBL)^NF3            -- TSH stimulates T4 synthesis
    #   feed1 = (T4BL*VDT4/AT4srm)^NF1      -- low serum T4 raises TSH production
    #   feed2 = (AT4srm/(T4BL*VDT4))^NF2    -- low serum T4 slows TSH turnover
    # The printed symbol T4BL in feed1/feed2 is the SERUM basal T4
    # (Table 5 T4_BL_srm = 58.4e-3 mg/L), NOT the thyroid-tissue T4_BL
    # (208 mg/L) that shares the printed name: only the serum reading is
    # dimensionally consistent (mg/L * L = mg, matching the amount AT4srm)
    # and only it makes feed1 = feed2 = 1 at homeostasis. See the Errata.
    c_tsh <- tsh_serum / vd_tsh
    t4_srm_bl_amt <- t4_bl_srm * vd_t4

    stim  <- (c_tsh / tsh_bl)^nf3
    feed1 <- (t4_srm_bl_amt / t4_serum)^nf1
    feed2 <- (t4_serum / t4_srm_bl_amt)^nf2

    # Zero-order TSH production. Table 5 prints "kin_TSH = TSH_BL * k_TSH",
    # which has units mg/L/h. That is self-consistent only if ATSH is carried
    # as a CONCENTRATION. This file carries all five hormone states as
    # AMOUNTS in mg (see compartmentData), for which the same physical system
    # requires the VDTSH factor below -- which also matches Table 5's own
    # Units column, where kin_TSH is mg/h.
    #
    # The two formulations are algebraically identical: the TSH state enters
    # the model only through c_tsh = tsh_serum / vd_tsh, so dividing this
    # ODE through by vd_tsh recovers the printed concentration form exactly.
    # Either way stim = 1 at homeostasis; taking the printed expression while
    # carrying ATSH as an amount is the one combination that does NOT work
    # (it would give stim = (1/vd_tsh)^NF3, a ~41% inflation of baseline T4
    # synthesis at WT = 0.3 kg). Documented in the vignette Errata.
    kin_tsh <- tsh_bl * vd_tsh * k_tsh

    # dAT4thy/dt = prodT4*Inhib*stim - metT4thy*AT4thy*fr - fT4srm*AT4thy
    d/dt(t4_thyroid) <- prod_t4 * inhib * stim -
                        met_t4thy * t4_thyroid * fr_t4_t3 -
                        f_t4srm * t4_thyroid

    # dAT3thy/dt = prodT3*Inhib + metT4thy*AT4thy*fr - metT3thy*AT3thy - fT3srm*AT3thy
    d/dt(t3_thyroid) <- prod_t3 * inhib +
                        met_t4thy * t4_thyroid * fr_t4_t3 -
                        met_t3thy * t3_thyroid -
                        f_t3srm * t3_thyroid

    # dATSH/dt = kinTSH*feed1 - kTSH*feed2*ATSH
    d/dt(tsh_serum) <- kin_tsh * feed1 - k_tsh * feed2 * tsh_serum

    # dAT4srm/dt = fT4srm*AT4thy - metT4srm*AT4srm*fr - losssrm*AT4srm
    d/dt(t4_serum) <- f_t4srm * t4_thyroid -
                      met_t4srm * t4_serum * fr_t4_t3 -
                      loss_srm * t4_serum

    # dAT3srm/dt = fT3srm*AT3thy + metT4srm*AT4srm*fr - metT3srm*AT3srm
    d/dt(t3_serum) <- f_t3srm * t3_thyroid +
                      met_t4srm * t4_serum * fr_t4_t3 -
                      met_t3srm * t3_serum

    # -------------------------------------------------------------------
    # 8. Initial conditions. The chemical states all start empty; the five
    #    thyroid-hormone states start at their reported basal levels
    #    (Table 5) converted from concentration to amount.
    # -------------------------------------------------------------------
    t4_thyroid(0) <- t4_bl_thy * v_thyroid_tissue
    t3_thyroid(0) <- t3_bl_thy * v_thyroid_tissue
    t4_serum(0)   <- t4_bl_srm * vd_t4
    t3_serum(0)   <- t3_bl_srm * vd_t3
    tsh_serum(0)  <- tsh_bl    * vd_tsh

    # -------------------------------------------------------------------
    # 9. Observation variables (all mg/L).
    #    Cc              -- oxyfluorfen in blood / serum      (Fig. 5)
    #    Cthyroid        -- oxyfluorfen in thyroid tissue     (Fig. 4)
    #    Cthyroid_blood  -- oxyfluorfen in thyroid blood; this is CBthy,
    #                       the driver of the NIS inhibition function
    #    T4serum         -- serum thyroxine                   (Fig. 6)
    #    T3serum         -- serum triiodothyronine            (Fig. 7)
    #    TSHserum        -- serum thyrotropin                 (Fig. 6 insert)
    # -------------------------------------------------------------------
    Cc             <- c_art
    Cthyroid       <- c_thyroid_tissue
    Cthyroid_blood <- c_thyroid_blood
    T4serum  <- t4_serum / vd_t4
    T3serum  <- t3_serum / vd_t3
    TSHserum <- tsh_serum / vd_tsh
  })
}
