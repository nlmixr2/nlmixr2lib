Decrane_2023_oxyfluorfen_human <- function() {
  description <- "QSP (PBPK + thyroid-hormone kinetics). Human cross-species extrapolation of the rat oxyfluorfen model of Decrane et al. 2023 (Curr Res Toxicol 5:100138), used by the authors to forecast the drop in serum T4 and T3 from long-term oral exposure to the herbicide oxyfluorfen in drinking water. Structurally identical to the rat companion model (13 flow-limited chemical compartments including a diffusion-limited thyroid split into thyroid blood and thyroid tissue, plus five thyroid-hormone states carrying the HPT feedback loop), but with human physiology (Table 2), an age-dependent cardiac output, a hepatic clearance derived from human in-vitro hepatocyte data rather than fitted, human basal thyroid-hormone levels, non-allometric human hormone production and transfer rates, and a DIFFERENT empirical sodium-iodide-symporter (NIS) inhibition function fitted to human hNIS-HEK293T-EPA in-vitro data. There is no human in-vivo calibration data for oxyfluorfen, so the authors caution that the predictions 'should be considered with caution'. Deterministic typical-value (no IIV, no residual error). Companion model: Decrane_2023_oxyfluorfen_rat."
  reference <- paste(
    "Decrane R, Stoker T, Murr A, Ford J, El-Masri H. (2023).",
    "Cross species extrapolation of the disruption of thyroid hormone synthesis by oxyfluorfen",
    "using in vitro data, physiologically based pharmacokinetic (PBPK), and thyroid hormone kinetics models.",
    "Current Research in Toxicology 5:100138.",
    "doi:10.1016/j.crtox.2023.100138. PMCID PMC10697989.",
    sep = " "
  )
  vignette <- "Decrane_2023_oxyfluorfen"

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
        "Table 2 fixes the reference human at 80 kg. Scales cardiac output (allometric",
        "exponent 0.75), every tissue volume (linear) and the three thyroid-hormone volumes",
        "of distribution (linear). Unlike the rat model, human hepatic clearance, GFR and the",
        "thyroid-hormone production and transfer rate constants are reported as absolute",
        "values and are NOT body-weight scaled.",
        sep = " "
      ),
      source_name        = "Body Weight (kg)"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      notes              = paste(
        "Table 2 fixes the reference human at 30 years. Enters ONLY through the cardiac-output",
        "relationship Cardiac Output = -6.846*log10(age) + 16.775 (L/h/kg^0.75). No other",
        "parameter is age-dependent. At age 30 this gives 6.662 L/h/kg^0.75; see the vignette",
        "Errata, which notes this is below the Brown et al. (1997) adult human value.",
        sep = " "
      ),
      source_name        = "Age(years)"
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
    species        = "human",
    n_subjects     = 0L,
    n_studies      = 0L,
    age_range      = "30 years (single reference adult, Table 2)",
    weight_range   = "80 kg (single reference adult, Table 2)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = "healthy reference adult",
    dose_range     = "long-term daily oral exposure via drinking water; the paper reports the water concentrations predicted to give a 10% drop in serum T4 (57 mg/L) and serum T3 (89 mg/L), and compares against 300 mg/L (the rat-derived equivalent), assuming 2 L/day water consumption at 80 kg",
    regions        = "United States (US EPA ORD)",
    notes          = paste(
      "This is a forward EXTRAPOLATION, not a fit: there were no human in-vivo oxyfluorfen",
      "data, so n_subjects = 0. Human physiology is from Brown et al. (1997) (Table 2);",
      "tissue:blood partition coefficients from GastroPlus 9.8 (Table 3, identical to the",
      "rat); hepatic clearance calculated from an in-vitro hepatic clearance of",
      "7.88 ul/min/10^6 cells (Wetmore et al. 2012); GFR from Delanaye et al. (2012);",
      "basal serum T4/T3 from Kinjo et al. (2002); thyroidal T4/T3 production from",
      "Wiersinga et al. (2012). Basal human thyroid-TISSUE hormone levels were scaled from",
      "the human serum levels using the rat thyroid:serum ratio, and basal human serum TSH",
      "was assumed similar to the rat. The Discussion states the model predictions",
      "'should be considered with caution'.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Human physiological parameters -- Table 2 (all from Brown et al. 1997).
    # =====================================================================

    # Cardiac output (L/h/kg^0.75) = -6.846*log10(age) + 16.775.
    # NB the leading minus sign is present in the PDF text layer; the
    # docling-converted markdown drops it. At age 30 this evaluates to
    # 6.662 L/h/kg^0.75 -> 178 L/h at 80 kg, which is below the Brown et al.
    # (1997) adult human cardiac output. Carried verbatim per the
    # trust-the-printed-equation rule; see the vignette Errata.
    qc_age_slope       <- fixed(-6.846)    ; label("Cardiac-output slope on log10(age) (L/h/kg^0.75)")   # Table 2, Cardiac Output = -6.846*log10(age) + 16.775
    qc_intercept       <- fixed(16.775)    ; label("Cardiac-output intercept (L/h/kg^0.75)")             # Table 2

    # Tissue volumes, fraction of body weight.
    # As in the rat model, Table 2's "Plasma" row (0.079) is used as the
    # volume of the well-mixed blood pool supplying Cart.
    fv_plasma          <- fixed(0.079)     ; label("Blood-pool volume, fraction of body weight (Table 2 'Plasma')")  # Table 2
    fv_thyroid_total   <- fixed(0.00005)   ; label("Total thyroid volume, fraction of body weight")            # Table 2, Thyroid Total
    fv_thyroid_blood   <- fixed(0.18)      ; label("Thyroid blood sub-compartment, fraction of TOTAL THYROID volume")  # Table 2, Thyroid Blood
    fv_fat             <- fixed(0.214)     ; label("Fat volume, fraction of body weight")                      # Table 2
    fv_muscle          <- fixed(0.40)      ; label("Muscle volume, fraction of body weight")                   # Table 2
    fv_skin            <- fixed(0.037)     ; label("Skin volume, fraction of body weight")                     # Table 2
    fv_liver           <- fixed(0.026)     ; label("Liver volume, fraction of body weight")                    # Table 2
    fv_gut             <- fixed(0.017)     ; label("GI-tract tissue volume, fraction of body weight")          # Table 2, GI tract
    fv_kidney          <- fixed(0.004)     ; label("Kidney volume, fraction of body weight")                   # Table 2
    fv_brain           <- fixed(0.02)      ; label("Brain volume, fraction of body weight")                    # Table 2
    fv_slowly_perfused <- fixed(0.143)     ; label("Slowly-perfused tissue volume, fraction of body weight")   # Table 2 (footnote 2: assumed mainly bone)
    fv_rapidly_perfused <- fixed(0.142)    ; label("Rapidly-perfused tissue volume, fraction of body weight")  # Table 2 (footnote 3: remaining body weight)

    # Blood flows, fraction of cardiac output. As in the rat table,
    # "Liver 0.227" = Portal Vein (0.181) + Hepatic Artery (0.046) exactly.
    fq_thyroid         <- fixed(0.016)     ; label("Thyroid blood flow, fraction of cardiac output")           # Table 2
    fq_fat             <- fixed(0.052)     ; label("Fat blood flow, fraction of cardiac output")               # Table 2
    fq_muscle          <- fixed(0.191)     ; label("Muscle blood flow, fraction of cardiac output")            # Table 2
    fq_skin            <- fixed(0.058)     ; label("Skin blood flow, fraction of cardiac output")              # Table 2
    fq_portal_vein     <- fixed(0.181)     ; label("Portal-vein blood flow, fraction of cardiac output")       # Table 2, Portal Vein
    fq_hepatic_artery  <- fixed(0.046)     ; label("Hepatic-artery blood flow, fraction of cardiac output")    # Table 2, Hepatic Artery
    fq_kidney          <- fixed(0.175)     ; label("Kidney blood flow, fraction of cardiac output")            # Table 2
    fq_brain           <- fixed(0.114)     ; label("Brain blood flow, fraction of cardiac output")             # Table 2
    fq_slowly_perfused <- fixed(0.042)     ; label("Slowly-perfused tissue blood flow, fraction of cardiac output")   # Table 2
    fq_rapidly_perfused <- fixed(0.081)    ; label("Rapidly-perfused tissue blood flow, fraction of cardiac output")  # Table 2

    # =====================================================================
    # Tissue:blood partition coefficients -- Table 3 (GastroPlus 9.8).
    # Table 3 gives IDENTICAL values for rat and human.
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
    # Oxyfluorfen biochemical parameters -- Table 4, Human Value column.
    # NONE of the human chemical parameters were estimated: PAThy and Ka are
    # "assumed like the rat value" (footnote 6), CLR is calculated from human
    # in-vitro hepatocyte data (footnote 1), frac is GastroPlus (footnote 2)
    # and GFR is literature (footnote 4). All therefore fixed().
    # =====================================================================
    pa_thyroid         <- fixed(0.005)     ; label("Thyroid blood-to-tissue permeability-area product (L/h)")  # Table 4, PAThy human (footnote 6: assumed like the rat value)
    clr                <- fixed(118)       ; label("Hepatic metabolic clearance (L/h)")                        # Table 4, CLR human (footnote 1: calculated from in vitro hepatic clearance of 7.88 ul/min/10^6 cells, Wetmore et al. 2012)
    lka                <- fixed(log(0.00375)) ; label("Log first-order oral absorption rate constant from gut lumen (1/h)")  # Table 4, Ka human = 0.00375 (footnote 6: assumed like the rat value)
    frac               <- fixed(0.08)      ; label("Fraction of oxyfluorfen unbound in blood (unitless)")      # Table 4, frac (footnote 2: GastroPlus 9.8)
    gfr                <- fixed(107)       ; label("Glomerular filtration rate (L/h as printed; see Errata)")   # Table 4, GFR human (footnote 4: Delanaye et al. 2012)
    e_wt_cl            <- fixed(0.75)      ; label("Allometric exponent on cardiac output (unitless)")          # Table 2, Cardiac Output units L/hr/kg^.75

    # =====================================================================
    # Empirical in-vitro NIS-inhibition function -- paper p. 5:
    #   Inhib_human = -0.2917*CBthy^0.2739 + 1.006
    # fitted to the Buckalew et al. (2020) hNIS-HEK293T-EPA human-cell
    # radioactive-iodide-uptake data in Fig. 3b. This is a POWER function,
    # structurally different from the rat exponential. At CBthy = 0 it
    # returns 1.006 (normalised to the drug-free synthesis rate) but, unlike
    # the rat form, it has NO lower asymptote and crosses zero at
    # CBthy = 91.8 mg/L, becoming negative above that. See the Errata.
    # =====================================================================
    inhib_coef         <- fixed(-0.2917)   ; label("NIS-inhibition power-function coefficient (unitless)")      # p. 5, Inhib_human
    inhib_exp          <- fixed(0.2739)    ; label("NIS-inhibition power-function exponent on thyroid oxyfluorfen concentration (unitless)")  # p. 5, Inhib_human
    inhib_intercept    <- fixed(1.006)     ; label("NIS-inhibition power-function intercept (unitless)")        # p. 5, Inhib_human

    # =====================================================================
    # Thyroid-hormone kinetics -- Table 5, Human Value column.
    # Every human TH value is a literature or assumed value (footnotes 2, 3,
    # 4, 5, 6, 7, 8) EXCEPT the four thyroidal / systemic metabolism rate
    # constants carrying footnote 1, which Table 5 reports as identical to
    # the rat calibrated values. Those four are left un-fixed() to match the
    # rat file's encoding of the same numbers.
    # =====================================================================
    fr_t4_t3           <- fixed(0.22)      ; label("Fraction of T4 that converts to T3 (unitless)")            # Table 5, fr (footnote 8: Silva and Matthews 1984)

    # Thyroidal synthesis and interconversion. NB the human production rates
    # are ABSOLUTE (mg/h), not body-weight scaled as in the rat model.
    prod_t4            <- fixed(0.0416)    ; label("Thyroidal T4 zero-order production rate (mg/h)")           # Table 5, T4prod human (footnote 2: Wiersinga et al. 2012)
    prod_t3            <- fixed(2.7e-3)    ; label("Thyroidal T3 zero-order production rate (mg/h)")           # Table 5, T3prod human (footnote 2)
    met_t4thy          <- 0.065            ; label("Thyroidal T4-to-T3 conversion rate constant (1/h)")        # Table 5, T4met human (footnote 1; same value as rat)
    met_t3thy          <- 0.045            ; label("Thyroidal T3 metabolism rate constant (1/h)")              # Table 5, T3met human (footnote 1)

    # Basal (homeostatic) levels.
    t4_bl_thy          <- fixed(288)       ; label("Basal T4 concentration in thyroid tissue (mg/L)")          # Table 5, T4_BL human (footnote 4: calculated from Kinjo et al. 2002)
    t3_bl_thy          <- fixed(33.4)      ; label("Basal T3 concentration in thyroid tissue (mg/L)")          # Table 5, T3_BL human (footnote 4)
    tsh_bl             <- fixed(5e-3)      ; label("Basal TSH concentration in serum (mg/L)")                 # Table 5, TSH_BL human (footnote 4)
    t4_bl_srm          <- fixed(85e-3)     ; label("Basal T4 concentration in serum (mg/L)")                   # Table 5, T4_BL_srm human (footnote 4)
    t3_bl_srm          <- fixed(1.50e-3)   ; label("Basal T3 concentration in serum (mg/L)")                   # Table 5, T3_BL_srm human (footnote 4)

    # TSH turnover and the HPT power-law feedback exponents (same as rat).
    k_tsh              <- fixed(5e-6)      ; label("First-order TSH loss rate constant (1/h)")                 # Table 5, k_TSH (footnote 3: Handa et al. 2021)
    nf1                <- fixed(2.53)      ; label("Feedback exponent NF1 on serum T4 driving TSH production (unitless)")  # Table 5, NF1 (footnote 5: Ekerot et al. 2013)
    nf2                <- fixed(1.9)       ; label("Feedback exponent NF2 on serum T4 modulating TSH turnover (unitless)")  # Table 5, NF2 (footnote 5)
    nf3                <- fixed(0.11)      ; label("Stimulation exponent NF3 on serum TSH driving thyroidal T4 synthesis (unitless)")  # Table 5, NF3 (footnote 5)

    # Thyroid-to-serum transfer (absolute, not body-weight scaled) and
    # systemic hormone disposition.
    f_t4srm            <- fixed(0.06)      ; label("Thyroid-to-serum T4 transfer rate constant (1/h)")         # Table 5, fT4 human (footnote 6: set similar to rat values from Handa et al. 2021)
    f_t3srm            <- fixed(0.02)      ; label("Thyroid-to-serum T3 transfer rate constant (1/h)")         # Table 5, fT3 human (footnote 6)
    met_t4srm          <- 0.085            ; label("Systemic T4-to-T3 conversion rate constant in serum (1/h)")  # Table 5, T4met_srm (footnote 1)
    met_t3srm          <- 0.045            ; label("Systemic T3 loss rate constant in serum (1/h)")            # Table 5, T3met_srm (footnote 1)
    loss_srm           <- 0.045            ; label("Systemic non-deiodinative T4 loss rate constant in serum (1/h)")  # Table 5, T4loss_srm (footnote 1)

    # Thyroid-hormone volumes of distribution (L), linear in body weight.
    f_vd_t4            <- fixed(0.149)     ; label("T4 volume of distribution per body weight (L/kg)")         # Table 5, VDT4 = 0.149*Avg.BW (footnote 7: Dubois and Dussault 1977)
    f_vd_t3            <- fixed(1.62)      ; label("T3 volume of distribution per body weight (L/kg)")         # Table 5, VDT3 = 1.62*Avg.BW (footnote 7)
    f_vd_tsh           <- fixed(0.149)     ; label("TSH volume of distribution per body weight (L/kg)")        # Table 5, VDTSH = 0.149*Avg.BW (footnote 7)

    # =====================================================================
    # NO IIV and NO residual error -- see the rat companion file and the
    # vignette Errata. This model is a forward extrapolation with no human
    # in-vivo data to fit, so there is nothing to estimate variability from.
    # =====================================================================
  })

  model({
    # -------------------------------------------------------------------
    # 1. Cardiac output and regional blood flows (L/h).
    #    Table 2: Cardiac Output = -6.846*log10(age) + 16.775 (L/h/kg^0.75).
    #    The age term is aliased to a local first so the sum carries only one
    #    bare population parameter (rxode2 mu-reference parser).
    # -------------------------------------------------------------------
    ka                 <- exp(lka)

    qc_age_term        <- qc_age_slope * log10(AGE)
    qc_per_kg          <- qc_age_term + qc_intercept
    qc                 <- qc_per_kg * WT^e_wt_cl

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

    q_ven_liver        <- q_portal_vein + q_hepatic_artery

    # Total flow drawn from the blood pool = the ACTUAL sum of the Table 2
    # flows (0.956 of cardiac output, versus 1.171 for the rat Table 1).
    # Using the actual sum keeps the circulatory loop exactly
    # mass-conserving. See the vignette Errata.
    q_total            <- q_ven_liver + q_kidney + q_muscle + q_skin + q_fat +
                          q_brain + q_slowly_perfused + q_rapidly_perfused + q_thyroid

    # -------------------------------------------------------------------
    # 2. Compartment volumes (L). Table 2 fractions of body weight;
    #    "Thyroid Blood 0.18" is a fraction of the TOTAL thyroid volume.
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
    # 3. Tissue concentrations (mg/L), as in the rat model.
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
    # 4. Chemical PBPK ODE system (13 states), paper pp. 2-3. Identical
    #    structure to the rat model; only the parameter values differ.
    # -------------------------------------------------------------------
    d/dt(stomach)            <- -ka * stomach

    d/dt(a_gut)              <- q_portal_vein * (c_art - cv_gut) + ka * stomach

    d/dt(a_liver)            <- q_hepatic_artery * c_art - q_ven_liver * cv_liver +
                                q_portal_vein * cv_gut - clr * cv_liver

    # NB the renal sink uses the TISSUE concentration c_kidney, per the
    # printed equation.
    d/dt(a_kidney)           <- q_kidney * (c_art - cv_kidney) - gfr * frac * c_kidney

    d/dt(a_muscle)           <- q_muscle            * (c_art - cv_muscle)
    d/dt(a_skin)             <- q_skin              * (c_art - cv_skin)
    d/dt(a_fat)              <- q_fat               * (c_art - cv_fat)
    d/dt(a_brain)            <- q_brain             * (c_art - cv_brain)
    d/dt(a_slowly_perfused)  <- q_slowly_perfused   * (c_art - cv_slowly_perfused)
    d/dt(a_rapidly_perfused) <- q_rapidly_perfused  * (c_art - cv_rapidly_perfused)

    d/dt(a_thyroid_blood)    <- q_thyroid * (c_art - c_thyroid_blood) +
                                pa_thyroid * (cv_thyroid_tissue - c_thyroid_blood)
    d/dt(a_thyroid_tissue)   <- -pa_thyroid * (cv_thyroid_tissue - c_thyroid_blood)

    # Blood pool: the one chemical ODE the paper does not print, reconstructed
    # by mass-balance closure of the ten printed tissue equations. See Errata.
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
    # 5. NIS inhibition, HUMAN form (paper p. 5, Fig. 3b):
    #      Inhib_human = -0.2917*CBthy^0.2739 + 1.006
    #    A power function, not the rat's exponential. It has no lower
    #    asymptote and turns negative above CBthy = 91.8 mg/L; the vignette
    #    checks that the simulated exposures stay well inside that domain.
    #
    #    CBthy IS THE THYROID *BLOOD* CONCENTRATION -- see the extended note
    #    in the rat companion file for the four lines of evidence (the
    #    symbol's "B"; the basolateral-NIS mechanism described in the paper's
    #    own Introduction; correspondence with the in-vitro medium
    #    concentration; and reproduction of BOTH species' reported points of
    #    departure). Under this reading 57 mg/L of oxyfluorfen in drinking
    #    water gives a 10.0% drop in serum T4, matching Fig. 8a exactly;
    #    reading CBthy as the thyroid TISSUE concentration gives 18.0%.
    #    Quantified side-by-side in the vignette Errata.
    #
    #    NUMERICAL GUARD (not a change to the published function): a
    #    non-integer power is undefined for a negative base, so the solver
    #    returns NaN for the whole system the moment the thyroid blood state
    #    takes a transient negative excursion of order -1e-20 while rising
    #    from its zero initial condition. The concentration is non-negative
    #    on the entire physical domain, so clamping the base at 0 leaves
    #    Inhib bit-identical everywhere the model is defined and only removes
    #    the non-physical excursion. Without it the model cannot be solved at
    #    ANY nonzero dose. Documented in the vignette Errata.
    # -------------------------------------------------------------------
    c_thyroid_blood_nonneg <- max(c_thyroid_blood, 0)
    inhib <- inhib_coef * c_thyroid_blood_nonneg^inhib_exp + inhib_intercept

    # -------------------------------------------------------------------
    # 6. Thyroid-hormone kinetics (paper pp. 3-4). Unlike the rat model,
    #    prod_t4, prod_t3, f_t4srm and f_t3srm are absolute rate constants
    #    (Table 5 human column) and are NOT body-weight scaled.
    # -------------------------------------------------------------------
    vd_t4  <- f_vd_t4  * WT
    vd_t3  <- f_vd_t3  * WT
    vd_tsh <- f_vd_tsh * WT

    c_tsh <- tsh_serum / vd_tsh
    t4_srm_bl_amt <- t4_bl_srm * vd_t4

    # T4BL in feed1/feed2 is the SERUM basal T4 (Table 5 T4_BL_srm), not the
    # thyroid-tissue T4_BL that shares the printed symbol -- see the rat file
    # and the vignette Errata.
    stim  <- (c_tsh / tsh_bl)^nf3
    feed1 <- (t4_srm_bl_amt / t4_serum)^nf1
    feed2 <- (t4_serum / t4_srm_bl_amt)^nf2

    # kin_TSH carries the VDTSH factor required because this file holds the
    # TSH state as an AMOUNT (mg) rather than the concentration that Table 5's
    # printed "kin_TSH = TSH_BL * k_TSH" implies. The two forms are
    # algebraically identical -- see the extended note in the rat companion
    # file and the vignette Errata.
    kin_tsh <- tsh_bl * vd_tsh * k_tsh

    d/dt(t4_thyroid) <- prod_t4 * inhib * stim -
                        met_t4thy * t4_thyroid * fr_t4_t3 -
                        f_t4srm * t4_thyroid

    d/dt(t3_thyroid) <- prod_t3 * inhib +
                        met_t4thy * t4_thyroid * fr_t4_t3 -
                        met_t3thy * t3_thyroid -
                        f_t3srm * t3_thyroid

    d/dt(tsh_serum) <- kin_tsh * feed1 - k_tsh * feed2 * tsh_serum

    d/dt(t4_serum) <- f_t4srm * t4_thyroid -
                      met_t4srm * t4_serum * fr_t4_t3 -
                      loss_srm * t4_serum

    d/dt(t3_serum) <- f_t3srm * t3_thyroid +
                      met_t4srm * t4_serum * fr_t4_t3 -
                      met_t3srm * t3_serum

    # -------------------------------------------------------------------
    # 7. Initial conditions: chemical states empty, hormone states at the
    #    Table 5 human basal levels converted to amounts.
    # -------------------------------------------------------------------
    t4_thyroid(0) <- t4_bl_thy * v_thyroid_tissue
    t3_thyroid(0) <- t3_bl_thy * v_thyroid_tissue
    t4_serum(0)   <- t4_bl_srm * vd_t4
    t3_serum(0)   <- t3_bl_srm * vd_t3
    tsh_serum(0)  <- tsh_bl    * vd_tsh

    # -------------------------------------------------------------------
    # 8. Observation variables (all mg/L), as in the rat model.
    #    Cthyroid_blood is CBthy, the driver of the NIS inhibition function.
    # -------------------------------------------------------------------
    Cc             <- c_art
    Cthyroid       <- c_thyroid_tissue
    Cthyroid_blood <- c_thyroid_blood
    T4serum  <- t4_serum / vd_t4
    T3serum  <- t3_serum / vd_t3
    TSHserum <- tsh_serum / vd_tsh
  })
}
