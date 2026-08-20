Trujillo_2025_proinsulin_qsp <- function() {
  description <- paste(
    "QSP. Diabetes quantitative systems pharmacology model of glucose",
    "homeostasis with a proinsulin sub-model (Trujillo 2025), used for",
    "in silico hypothesis testing of a hypothetical therapy that converts",
    "circulating proinsulin to insulin in type 2 diabetes mellitus.",
    "Fifteen ODEs (Table S1) driven by 42 flux / rule expressions",
    "(Table S2): gastric emptying of carbohydrate and non-carbohydrate",
    "calories, intestinal glucose absorption, hepatic glucose output",
    "regulated by a glucose / insulin / glucagon ratio (HRS), muscle",
    "uptake via insulin-independent GLUT1 and insulin- and",
    "proinsulin-responsive GLUT4, brain and other-tissue uptake, renal",
    "filtration / reabsorption / urinary excretion, glucose-stimulated",
    "insulin secretion (GSIS) and proinsulin secretion (GSPS) both",
    "amplified by active GLP-1 incretin, two-compartment insulin",
    "disposition, one-compartment proinsulin disposition, glucagon",
    "secretion with a T2DM dysregulation factor, and 30-day / 90-day",
    "glucose delay states that generate HbA1c. Meals enter as dosing",
    "events on the gastric carbohydrate and gastric other compartments",
    "(kcal); 250 mg of glucose reaches the intestine per kcal of",
    "carbohydrate. The hypothetical drug adds first-order conversion of",
    "plasma proinsulin to plasma insulin, gated by the ON_TREATMENT",
    "covariate (Eq 6). Deterministic: no IIV and no residual error are",
    "reported. Default ini() values are the VPHealthy parameter set;",
    "the three type 2 diabetic virtual patients (VPT2DM-1, VPT2DM-2,",
    "VPT2DM-3) differ in ten parameters that are given in the",
    "source-trace comments and reproduced in the validation vignette.",
    sep = " "
  )
  reference <- paste(
    "Trujillo ME, Han Y, Baillie RA, Weis MC, Chung D, Hayes S,",
    "Carrington PE, Reed M (2025).",
    "In Silico Hypothesis Testing in Drug Discovery: Using Quantitative",
    "Systems Pharmacology Modeling to Evaluate the Therapeutic Value of",
    "Proinsulin Conversion to Insulin Therapy for Type 2 Diabetes",
    "Mellitus. Pharmaceutics 17(12):1522.",
    "doi:10.3390/pharmaceutics17121522. PMCID: PMC12736522.",
    "Model equations, fluxes and parameters are in the Supplementary",
    "Materials (Tables S1, S2 and S3 respectively); the main text",
    "prints Equations 1-6 only.",
    sep = " "
  )
  vignette <- "Trujillo_2025_proinsulin_qsp"

  paper_specific_compartments <- c(
    "gastric_cho", "gastric_other", "gut_glucose", "glucose",
    "tubular_glucose", "insulin", "insulin_per", "glp1", "glp1_inactive",
    "glucagon", "glucose_delay30", "glucose_delay90", "proinsulin",
    "insulin_musc", "proinsulin_musc"
  )

  units <- list(
    time          = "min",
    dosing        = "kcal",
    concentration = "mg/dL"
  )

  # Meals enter as FoodCarbonhydrateIntake / FoodOtherIntake on the two gastric
  # states (Table S1). Declared explicitly because buildModelDb()'s dosing
  # heuristic recognises only "depot" and "central".
  dosing <- c("gastric_cho", "gastric_other")

  covariateData <- list(
    ON_TREATMENT = list(
      description        = "Proinsulin-to-insulin conversion therapy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (untreated)",
      notes              = paste(
        "1 = the hypothetical proinsulin-converting drug is active, adding the",
        "first-order ProinsulinConv flux (Eq 6, rate constant 1/K_CPro) that",
        "moves plasma proinsulin into the plasma insulin pool; 0 = untreated,",
        "ProinsulinConv is zero. Trujillo 2025 Section 2.2 states the",
        "conversion term was 'added to the model' to represent the drug, so it",
        "is absent from the untreated base model; K_CPro = 20 min in Table S3",
        "is the treated value chosen in Section 2.6. The untreated reading is",
        "confirmed numerically: with ProinsulinConv switched off the VPHealthy",
        "fasting steady state is 5.5 pM proinsulin against the 5 pM of Table 2,",
        "whereas leaving it on gives 1.3 pM. Time-varying use is intended",
        "(treatment starts partway through a simulation).",
        sep = " "
      ),
      source_name        = "ON_TREATMENT"
    )
  )

  compartmentData <- list(
    gastric_cho      = list(analyte = "dietary carbohydrate", units = "kcal", specimen = "administration site", verified = TRUE),
    gastric_other    = list(analyte = "dietary non-carbohydrate calories", units = "kcal", specimen = "administration site", verified = TRUE),
    gut_glucose      = list(analyte = "glucose", units = "mg", specimen = "administration site", verified = TRUE),
    glucose          = list(analyte = "glucose", units = "mg", specimen = "plasma", verified = TRUE),
    tubular_glucose  = list(analyte = "glucose", units = "mg", specimen = "urine", verified = TRUE),
    insulin          = list(analyte = "insulin", units = "nmol", specimen = "plasma", verified = TRUE),
    insulin_per      = list(analyte = "insulin", units = "nmol", specimen = "tissue", verified = TRUE),
    glp1             = list(analyte = "active GLP-1", units = "pmol", specimen = "plasma", verified = TRUE),
    glp1_inactive    = list(analyte = "inactive GLP-1", units = "pmol", specimen = "plasma", verified = TRUE),
    glucagon         = list(analyte = "glucagon", units = "pmol", specimen = "plasma", verified = TRUE),
    glucose_delay30  = list(analyte = "glucose", units = "mg/dL", specimen = "not applicable", verified = TRUE),
    glucose_delay90  = list(analyte = "glucose", units = "mg/dL", specimen = "not applicable", verified = TRUE),
    proinsulin       = list(analyte = "proinsulin", units = "nmol", specimen = "plasma", verified = TRUE),
    insulin_musc     = list(analyte = "insulin", units = "uU/mL", specimen = "not applicable", verified = TRUE),
    proinsulin_musc  = list(analyte = "proinsulin", units = "uU/mL equivalent", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species       = "human",
    n_subjects    = 4,
    n_studies     = 2,
    disease_state = paste(
      "Four virtual patients spanning healthy physiology and the type 2",
      "diabetes spectrum: VPHealthy, VPT2DM-1 and VPT2DM-2 (early / moderate",
      "T2DM with insulin resistance and pancreatic compensation), and",
      "VPT2DM-3 (late-stage T2DM with pancreatic beta-cell failure)",
      "(Trujillo 2025 Section 2.5, Table 1).",
      sep = " "
    ),
    dose_range    = paste(
      "Three meals per day for 52 weeks, plus a 75 g oral glucose tolerance",
      "test after a 16 h fast (Section 2.6). The hypothetical conversion",
      "therapy was explored at conversion rates up to ~7 pM/min.",
      sep = " "
    ),
    notes         = paste(
      "The virtual patients are model constructs, not fitted individuals.",
      "They were calibrated against fasting insulin, proinsulin and the",
      "proinsulin/insulin ratio from the literature and from baseline samples",
      "of 307 participants with T2DM in two completed Merck Phase 3 studies",
      "(Study 11 / MK-3102-011, NCT01717313, n = 181; Study 24 /",
      "MK-3102-024, NCT01755156, n = 127) (Sections 2.3 and 3.1, Table 2).",
      "The underlying Diabetes QSP model was validated against the insulin",
      "glargine arm of Zinman 2012 (BEGIN Once Long) (Figure S1). No",
      "individual-level demographics are reported for the virtual patients.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # ALL parameters are FIXED. Trujillo 2025 is a hand-built SimBiology QSP
    # model: nothing is estimated, no IIV and no residual error are reported.
    # Values are Trujillo 2025 Supplementary Table S3, VPHealthy column,
    # unless noted. Where a parameter differs between virtual patients the
    # source-trace comment lists all four columns in the table's order
    # (VPHealthy / VPT2DM-1 / VPT2DM-2 / VPT2DM-3).
    #
    # Time is in MINUTES. Glucose amounts are in mg and glucose
    # concentration in mg/dL; insulin and proinsulin amounts are in nmol,
    # giving nmol/L (= nM) concentrations that the model reports as pM after
    # the x1000 in the observation block. Gastric contents are in kcal.
    #
    # Two rows of Table S3 are deliberately absent from ini() because no
    # flux or rule in Table S2 references them -- they are virtual-patient
    # CONSTRUCTION factors, not runtime parameters (see the note below the
    # secretion block):
    #   bfx     = 1 / 0.8 / 0.7 / 0.5           (beta-cell function factor)
    #   Km_MIns = 37.6137 in all four VPs       (reference muscle-insulin Km)
    # =====================================================================

    # --- Glucose distribution -------------------------------------------
    lvd_glucose  <- fixed(log(140))    ; label("Plasma glucose distribution volume Vd_PlasmaGlucose (dL)")            # Table S3 Vd_PlasmaGlucose

    # --- Insulin and proinsulin disposition ------------------------------
    lv_insulin   <- fixed(log(5.11))   ; label("Insulin central volume V_Insulin (L)")                                 # Table S3 V_Insulin
    lvp_insulin  <- fixed(log(31.6))   ; label("Insulin peripheral volume VP_Insulin (L)")                             # Table S3 VP_Insulin
    lq_insulin   <- fixed(log(0.4083)) ; label("Insulin inter-compartmental clearance Q_Insulin (L/min)")              # Table S3 Q_Insulin
    lcl_insulin  <- fixed(log(0.7283)) ; label("Insulin clearance CL_Insulin (L/min)")                                 # Table S3 CL_Insulin
    lv_proins    <- fixed(log(5.11))   ; label("Proinsulin volume V_Proinsulin (L); assumed equal to V_Insulin")       # Table S3 V_Proinsulin; Section 2.2 assumption (3)
    lcl_proins   <- fixed(log(0.084))  ; label("Proinsulin clearance CL_Proinsulin (L/min)")                           # Table S3 CL_Proinsulin; Section 2.2 assumption (2), insulin t1/2 5-10x shorter

    # --- Gastric emptying and intestinal absorption ----------------------
    lvmax_ge     <- fixed(log(8))      ; label("Maximum gastric emptying rate Vmax_GE (kcal/min)")                     # Table S3 Vmax_GE
    lkm_ge       <- fixed(log(100))    ; label("Gastric emptying Km on total gastric contents Km_GE (kcal)")           # Table S3 Km_GE
    lkm_gei      <- fixed(log(4000))   ; label("Intestinal feedback inhibition Km Km_GEi (mg glucose)")                # Table S3 Km_GEi
    lk_absorp    <- fixed(log(0.1))    ; label("Intestinal glucose absorption rate constant K_Absorption (1/min)")     # Table S3 K_Absorption

    # --- Hepatic glucose handling ----------------------------------------
    lvmax_g6p    <- fixed(log(560))    ; label("Maximum hepatic glucose release Vmax_G6P (mg/min) -- DERIVED, not tabulated")
    # ^^ NOT REPORTED ANYWHERE ON DISK. Vmax_G6P is used by the Ra_Liver rule
    #    in Table S2 but has no row in Table S3. Derived here from the paper's
    #    own numbers, not substituted from outside the paper:
    #      (a) At the model's own normal reference point (Conc_Glucose_Pl =
    #          NormalGlucose, EffectivePlasmaInsulin = NormalInsulin,
    #          PlasmaGlucagon = NormalGlucagon) the HRS rule gives HRS = 1
    #          exactly, so Ra_Liver = Vmax_G6P / (1 * Km_Ra_Liver + 1) =
    #          Vmax_G6P / 4 for VPHealthy. Vmax_G6P = 560 makes that
    #          140 mg/min, i.e. 2.0 mg/kg/min for a 70 kg adult.
    #      (b) Inverting the whole fasting steady state so that VPHealthy
    #          reproduces its Table 1 fasting glucose of 93 mg/dL returns
    #          Vmax_G6P = 559.0, i.e. 560 to the two significant figures the
    #          fasting glucose is reported to.
    #      (c) Cross-check on virtual patients not used in the derivation:
    #          with Vmax_G6P = 560 the fasting steady state predicts
    #          168.4 mg/dL for VPT2DM-3 (reported 167) and 125.8 mg/dL for
    #          VPT2DM-2 (reported 131).
    #    The derivation is reproduced end to end in the validation vignette.

    lbasal_gk    <- fixed(log(5))      ; label("Basal hepatic glucose uptake Basal_GK (mg/min)")                       # Table S3 Basal_GK
    lka_gk       <- fixed(log(0.3))    ; label("Hepatic glucose uptake per unit absorption Ka_GK (unitless)")          # Table S3 Ka_GK
    lkm_ra_liver <- fixed(log(3))      ; label("Hepatic glucose release suppression Km_Ra_Liver (unitless)")           # Table S3 Km_Ra_Liver = 3 / 1.5 / 2.4 / 2.1

    # --- Reference values entering the HRS regulation ratio --------------
    lnormal_glu  <- fixed(log(90))     ; label("Reference plasma glucose NormalGlucose (mg/dL)")                       # Table S3 NormalGlucose
    lnormal_ins  <- fixed(log(5))      ; label("Reference effective plasma insulin NormalInsulin (uU/mL)")             # Table S3 NormalInsulin
    lnormal_ggn  <- fixed(log(75))     ; label("Reference plasma glucagon NormalGlucagon (pmol)")                      # Table S3 NormalGlucagon

    # --- Muscle glucose uptake (GLUT1 / GLUT4) ---------------------------
    lvmax_mg1    <- fixed(log(24))     ; label("GLUT1 maximum muscle glucose uptake Vmax_MG1 (mg/min)")                # Table S3 Vmax_MG1
    lkm_mg1      <- fixed(log(18))     ; label("GLUT1 glucose Km Km_MG1 (mg/dL)")                                      # Table S3 Km_MG1
    lvmax_mg4    <- fixed(log(1800))   ; label("GLUT4 maximum muscle glucose uptake Vmax_MG4 (mg/min)")                # Table S3 Vmax_MG4
    lkm_mg4      <- fixed(log(180))    ; label("GLUT4 glucose Km Km_MG4 (mg/dL)")                                      # Table S3 Km_MG4
    lkm_ins_g4   <- fixed(log(37.6137)); label("Insulin Km for GLUT4 recruitment Km_Ins_GLUT4 (uU/mL)")                # Table S3 Km_Ins_GLUT4 = 37.6137 / 75.2274 / 47.017125 / 53.73386
    lkm_pro_g4   <- fixed(log(5187.19)); label("Proinsulin Km for GLUT4 recruitment Km_Pro_GLUT4 (uU/mL equivalent)")  # Table S3 Km_Pro_GLUT4 = 5187.19 / 10374.38 / 6483.9875 / 7410.271
    lh_mins      <- fixed(log(3.0612)) ; label("Hill exponent for GLUT4 recruitment h_MIns (unitless)")                # Table S3 h_MIns
    lk_mins      <- fixed(log(30))     ; label("Plasma-to-muscle insulin equilibration time K_MIns (min)")             # Table S3 K_MIns
    lk_mpro      <- fixed(log(30))     ; label("Plasma-to-muscle proinsulin equilibration time K_MPro (min)")          # Table S3 K_MPro

    # --- Brain and other-tissue glucose uptake ---------------------------
    lvmax_brain  <- fixed(log(96))     ; label("Maximum brain glucose uptake Vmax_Brain (mg/min)")                     # Table S3 Vmax_Brain
    lkm_brain    <- fixed(log(18))     ; label("Brain glucose uptake Km Km_Brain (mg/dL)")                             # Table S3 Km_Brain
    lvmax_ot     <- fixed(log(90))     ; label("Maximum other-tissue glucose uptake Vmax_OT (mg/min)")                 # Table S3 Vmax_OT
    lkm_ot       <- fixed(log(180))    ; label("Other-tissue glucose uptake Km Km_OT (mg/dL)")                         # Table S3 Km_OT

    # --- Renal handling ---------------------------------------------------
    lka_tubglu   <- fixed(log(1.13))   ; label("Glomerular filtration of glucose Ka_TubGlu (dL/min)")                  # Table S3 Ka_TubGlu
    ltmg         <- fixed(log(180))    ; label("Maximum tubular glucose reabsorption TmG (mg/min)")                    # Table S3 TmG
    lkm_ugr      <- fixed(log(100))    ; label("Tubular glucose reabsorption Km Km_UGR (mg)")                          # Table S3 Km_UGR
    lka_uge      <- fixed(log(0.01))   ; label("Urinary glucose excretion rate constant Ka_UGE (1/min)")               # Table S3 Ka_UGE

    # --- Glucose-stimulated insulin and proinsulin secretion --------------
    # Ten rows of Table S3 vary across the virtual patients. SIX of them are
    # generated by exactly two factors: an insulin-resistance factor
    # IR = 1 / 2 / 1.25 / 1.428571 and the tabulated beta-cell factor
    # bfx = 1 / 0.8 / 0.7 / 0.5, via
    #   Km_Ra_Liver  = 3 / IR            Km_Ins_GLUT4 = 37.6137 * IR
    #   Km_Pro_GLUT4 = 5187.19 * IR      Km_GSIS      = 155 / bfx
    #   Km_GSPS      = 350 / bfx         Vmax_GSIS = Vmax_GSPS = 0.2 * IR * bfx
    # These six relations reproduce all 24 tabulated values (6 parameters x 4
    # columns) to six significant figures, which is why bfx and Km_MIns appear
    # in Table S3 but in no flux: they are the generators, not runtime inputs.
    # The other four VP-varying rows (n_GSPS, afx, and bfx / Km_MIns) are set
    # independently and are NOT reproduced by these relations.
    lvmax_gsis   <- fixed(log(0.2))    ; label("Maximum glucose-stimulated insulin secretion Vmax_GSIS (nmol/min)")    # Table S3 Vmax_GSIS = 0.2 / 0.32 / 0.175 / 0.142857
    lkm_gsis     <- fixed(log(155))    ; label("Glucose Km for insulin secretion Km_GSIS (mg/dL)")                     # Table S3 Km_GSIS = 155 / 193.75 / 221.42857 / 310
    lvmax_gsps   <- fixed(log(0.2))    ; label("Maximum glucose-stimulated proinsulin secretion Vmax_GSPS (nmol/min)") # Table S3 Vmax_GSPS = 0.2 / 0.32 / 0.175 / 0.142857
    lkm_gsps     <- fixed(log(350))    ; label("Glucose Km for proinsulin secretion Km_GSPS (mg/dL)")                  # Table S3 Km_GSPS = 350 / 437.5 / 500 / 700
    ln_gsps      <- fixed(log(4.5))    ; label("Hill exponent for proinsulin secretion n_GSPS (unitless)")             # Table S3 n_GSPS = 4.5 / 4 / 4.3 / 4.2

    # --- Incretin (GLP-1) axis -------------------------------------------
    lvmax_incr   <- fixed(log(5))      ; label("Maximum incretin amplification of secretion Vmax_Incretin (unitless)") # Table S3 Vmax_Incretin
    lkm_incr     <- fixed(log(23))     ; label("Active GLP-1 Km for incretin amplification Km_Incretin (pmol)")        # Table S3 Km_Incretin
    lbasal_glp   <- fixed(log(1))      ; label("Basal active GLP-1 production Basal_GLP (pmol/min)")                   # Table S3 Basal_GLP
    lka_aglp     <- fixed(log(0.00025)); label("Intestinal glucose-driven GLP-1 production Ka_aGLP (pmol/min/mg)")     # Table S3 Ka_aGLP
    lkd_dpp4     <- fixed(log(0.1))    ; label("Active GLP-1 inactivation by DPP-4 Kd_DPP4 (1/min)")                   # Table S3 Kd_DPP4
    lkd_aglp     <- fixed(log(0.1))    ; label("Active GLP-1 inactivation by other enzymes Kd_aGLP (1/min)")           # Table S3 Kd_aGLP
    lkd_iglp     <- fixed(log(0.05))   ; label("Inactive GLP-1 clearance Kd_iGLP (1/min)")                             # Table S3 Kd_iGLP

    # --- Glucagon secretion ------------------------------------------------
    lnormal_vggn <- fixed(log(20))     ; label("Maximum glucose-suppressible glucagon secretion normal_Vm_GGN (pmol/min)")  # Table S3 normal_Vm_GGN
    lnormal_bggn <- fixed(log(3.25))   ; label("Basal glucagon secretion normal_BasalGGN (pmol/min)")                       # Table S3 normal_BasalGGN
    lkm_ggn      <- fixed(log(0.02))   ; label("Glucose Km for glucagon suppression Km_GGN (dL/mg)")                        # Table S3 Km_GGN
    ldefect_bggn <- fixed(log(4))      ; label("Glucagon dysregulation basal secretion defective_BasalGGN (pmol/min)")      # Table S3 defective_BasalGGN
    lkm_def_bggn <- fixed(log(0.75))   ; label("Km of the glucagon dysregulation term Km_defective_BasalGGN (unitless)")    # Table S3 Km_defective_BasalGGN
    afx          <- fixed(1)           ; label("Glucagon regulation intactness factor afx (1 = intact, unitless)")          # Table S3 afx = 1 / 0.8 / 0.9 / 0.8

    # --- Hypothetical proinsulin-to-insulin conversion therapy ------------
    lk_cpro      <- fixed(log(20))     ; label("Proinsulin-to-insulin conversion time constant K_CPro (min); gated by ON_TREATMENT")  # Table S3 K_CPro; Eq 6; Section 2.6 chose the ~7 pM/min rate this produces
  })

  model({
    # ===================================================================
    # Parameters back-transformed from the log scale.
    # ===================================================================
    vd_glucose  <- exp(lvd_glucose)
    v_insulin   <- exp(lv_insulin)
    vp_insulin  <- exp(lvp_insulin)
    q_insulin   <- exp(lq_insulin)
    cl_insulin  <- exp(lcl_insulin)
    v_proins    <- exp(lv_proins)
    cl_proins   <- exp(lcl_proins)
    vmax_ge     <- exp(lvmax_ge)
    km_ge       <- exp(lkm_ge)
    km_gei      <- exp(lkm_gei)
    k_absorp    <- exp(lk_absorp)
    vmax_g6p    <- exp(lvmax_g6p)
    basal_gk    <- exp(lbasal_gk)
    ka_gk       <- exp(lka_gk)
    km_ra_liver <- exp(lkm_ra_liver)
    normal_glu  <- exp(lnormal_glu)
    normal_ins  <- exp(lnormal_ins)
    normal_ggn  <- exp(lnormal_ggn)
    vmax_mg1    <- exp(lvmax_mg1)
    km_mg1      <- exp(lkm_mg1)
    vmax_mg4    <- exp(lvmax_mg4)
    km_mg4      <- exp(lkm_mg4)
    km_ins_g4   <- exp(lkm_ins_g4)
    km_pro_g4   <- exp(lkm_pro_g4)
    h_mins      <- exp(lh_mins)
    k_mins      <- exp(lk_mins)
    k_mpro      <- exp(lk_mpro)
    vmax_brain  <- exp(lvmax_brain)
    km_brain    <- exp(lkm_brain)
    vmax_ot     <- exp(lvmax_ot)
    km_ot       <- exp(lkm_ot)
    ka_tubglu   <- exp(lka_tubglu)
    tmg         <- exp(ltmg)
    km_ugr      <- exp(lkm_ugr)
    ka_uge      <- exp(lka_uge)
    vmax_gsis   <- exp(lvmax_gsis)
    km_gsis     <- exp(lkm_gsis)
    vmax_gsps   <- exp(lvmax_gsps)
    km_gsps     <- exp(lkm_gsps)
    n_gsps      <- exp(ln_gsps)
    vmax_incr   <- exp(lvmax_incr)
    km_incr     <- exp(lkm_incr)
    basal_glp   <- exp(lbasal_glp)
    ka_aglp     <- exp(lka_aglp)
    kd_dpp4     <- exp(lkd_dpp4)
    kd_aglp     <- exp(lkd_aglp)
    kd_iglp     <- exp(lkd_iglp)
    normal_vggn <- exp(lnormal_vggn)
    normal_bggn <- exp(lnormal_bggn)
    km_ggn      <- exp(lkm_ggn)
    defect_bggn <- exp(ldefect_bggn)
    km_def_bggn <- exp(lkm_def_bggn)
    k_cpro      <- exp(lk_cpro)

    # ===================================================================
    # Concentrations and effective (assay-scale) hormone levels.
    # Table S2 rules 2-6. The x1000 converts nmol/L to pmol/L and the /6
    # converts pmol/L to the uU/mL scale on which NormalInsulin = 5 and the
    # GLUT4 Km values are expressed.
    # ===================================================================
    conc_glucose    <- glucose / vd_glucose                    # Table S2 Conc_Glucose_Pl
    conc_insulin    <- insulin / v_insulin                     # Table S2 Conc_Insulin_Pl
    conc_proinsulin <- proinsulin / v_proins                   # Table S2 Conc_Proinsulin_Pl
    eff_plasma_ins  <- 1000 * insulin / (6 * v_insulin)        # Table S2 EffectivePlasmaInsulin
    eff_plasma_pro  <- 1000 * proinsulin / (6 * v_proins)      # Table S2 EffectivePlasmaProinsulin

    # ===================================================================
    # Gastric emptying and intestinal absorption. Table S2 rules 7-12.
    # ===================================================================
    gastric_total   <- gastric_cho + gastric_other                                        # Table S2 GastricTotal
    intest_inhib    <- 1 - gut_glucose / (km_gei + gut_glucose)                           # Table S2 IntestinalInhibition
    gastric_cho_emp <- intest_inhib * vmax_ge * gastric_cho / (km_ge + gastric_total)     # Table S2 GastricCarbohydrateEmptying
    gastric_oth_emp <- intest_inhib * vmax_ge * gastric_other / (km_ge + gastric_total)   # Table S2 GastricOtherEmptying
    glucose_absorp  <- k_absorp * gut_glucose                                             # Table S2 GlucoseAbsorption

    # ===================================================================
    # Hepatic glucose output. Table S2 rules 13-16.
    # HRS is the glucose x insulin / glucagon regulation ratio; it equals 1
    # when all three are at their reference values.
    # ===================================================================
    rd_liver  <- basal_gk + ka_gk * glucose_absorp                                        # Table S2 Rd_Liver
    hrs       <- ((conc_glucose / normal_glu) * (eff_plasma_ins / normal_ins)) /
      (glucagon / normal_ggn)                                                             # Table S2 HRS
    ra_liver  <- vmax_g6p / (hrs * km_ra_liver + 1)                                       # Table S2 Ra_Liver
    hep_gluc_out <- ra_liver - rd_liver                                                   # Table S2 HepaticGlucoseOutput

    # ===================================================================
    # Muscle glucose uptake. Table S2 rules 17-21.
    # GLUT4 recruitment is driven by effective MUSCLE insulin and proinsulin
    # (the delayed states), each through its own Hill term with a shared
    # exponent; proinsulin's Km is ~138-fold higher, giving it the weaker
    # effect noted in Section 2.2.
    # ===================================================================
    bound_g4  <- insulin_musc^h_mins / (km_ins_g4^h_mins + insulin_musc^h_mins) +
      proinsulin_musc^h_mins / (km_pro_g4^h_mins + proinsulin_musc^h_mins)                # Table S2 BoundG4
    glut1     <- vmax_mg1 * conc_glucose / (km_mg1 + conc_glucose)                        # Table S2 GLUT1
    glut4     <- bound_g4 * vmax_mg4 * conc_glucose / (km_mg4 + conc_glucose)             # Table S2 GLUT4
    musc_gluc_uptake <- glut1 + glut4                                                     # Table S2 MuscleGlucoseUptake
    ins_equil <- (1 / k_mins) * (eff_plasma_ins - insulin_musc)                           # Table S2 InsulinEquilibrium
    pro_equil <- (1 / k_mpro) * (eff_plasma_pro - proinsulin_musc)                        # Table S2 ProinsulinEqulibrium (spelled without the second "i" in the source)

    # ===================================================================
    # Brain, other tissue and renal handling. Table S2 rules 22-27.
    # ===================================================================
    brain_gluc_uptake <- vmax_brain * conc_glucose / (km_brain + conc_glucose)            # Table S2 BrainGlucoseUptake
    other_gluc_uptake <- vmax_ot * conc_glucose / (km_ot + conc_glucose)                  # Table S2 OtherTissueGlucoseUptake
    filtration    <- ka_tubglu * conc_glucose                                             # Table S2 Filtration
    reabsorption  <- tmg * tubular_glucose / (km_ugr + tubular_glucose)                   # Table S2 Reabsorption
    gluc_pl_to_tub <- filtration - reabsorption                                           # Table S2 GlucosePlasmaToTubular
    tub_gluc_excr <- ka_uge * tubular_glucose                                             # Table S2 TubularGlucoseExcretion

    # ===================================================================
    # Insulin and proinsulin secretion, distribution and clearance.
    # Table S2 rules 28-33 and 37-40.
    # ===================================================================
    incretin  <- vmax_incr * glp1 / (km_incr + glp1)                                      # Table S2 Incretin
    gsis      <- vmax_gsis * conc_glucose^4 / (km_gsis^4 + conc_glucose^4)                # Table S2 GSIS (Hill exponent hardcoded to 4 in the source)
    gsps      <- vmax_gsps * conc_glucose^n_gsps / (km_gsps^n_gsps + conc_glucose^n_gsps) # Table S2 GSPS
    insulin_secr    <- gsis * incretin                                                    # Table S2 InsulinSecretion
    proinsulin_secr <- gsps * incretin                                                    # Table S2 ProinsulinSecretion
    ins_per_to_pl   <- insulin_per * q_insulin / vp_insulin -
      insulin * q_insulin / v_insulin                                                     # Table S2 InsulinPeriphToPlasma
    insulin_clear    <- insulin * cl_insulin / v_insulin                                  # Table S2 InsulinClearance
    proinsulin_clear <- proinsulin * cl_proins / v_proins                                 # Table S2 ProinsulinClearance

    # Hypothetical proinsulin-converting drug (Eq 6). Section 2.2 states the
    # term was ADDED to the model to represent the therapy, so it is absent
    # from the untreated base model; ON_TREATMENT carries that switch.
    proinsulin_conv <- ON_TREATMENT * proinsulin / k_cpro                                 # Table S2 ProinsulinConv; Trujillo 2025 Eq 6

    # ===================================================================
    # GLP-1 and glucagon. Table S2 rules 34-36 and 41-42.
    # ===================================================================
    intest_glp_prod <- ka_aglp * gut_glucose                                              # Table S2 IntestinalGLPProduction
    glp_inact_dpp4  <- kd_dpp4 * glp1                                                     # Table S2 GLPInactivationByDPP4
    glp_inact_other <- kd_aglp * glp1                                                     # Table S2 GLPInactivationByGLP
    inactive_glp_cl <- kd_iglp * glp1_inactive                                            # Table S2 InactiveGLPClearance
    glucagon_secr   <- normal_bggn +
      defect_bggn * (1 - afx) / (km_def_bggn + (1 - afx)) +
      normal_vggn * afx / ((km_ggn^4) * (conc_glucose^4) + 1)                             # Table S2 PlasmaGlucagonSecretion
    glucagon_clear  <- glucagon / 15                                                      # Table S2 PlasmaGlucagonClearance (15 min time constant hardcoded in the source)

    # ===================================================================
    # ODE system -- Supplementary Table S1, in the order printed there.
    # Meals are dosing events on gastric_cho (FoodCarbonhydrateIntake) and
    # gastric_other (FoodOtherIntake), both in kcal. The factor 250 in the
    # intestinal glucose balance is mg of glucose per kcal of carbohydrate
    # (4 kcal/g x 1000 mg/g gives 250 mg/kcal), so a 75 g / 300 kcal oral
    # glucose load delivers exactly 75000 mg of glucose to the intestine.
    # ===================================================================
    d/dt(gastric_cho)     <- -gastric_cho_emp
    d/dt(gastric_other)   <- -gastric_oth_emp
    d/dt(gut_glucose)     <- 250 * gastric_cho_emp - glucose_absorp
    d/dt(glucose)         <- glucose_absorp + hep_gluc_out - musc_gluc_uptake -
      brain_gluc_uptake - other_gluc_uptake - gluc_pl_to_tub
    d/dt(tubular_glucose) <- gluc_pl_to_tub - tub_gluc_excr
    d/dt(insulin)         <- insulin_secr + ins_per_to_pl - insulin_clear + proinsulin_conv
    d/dt(insulin_per)     <- -ins_per_to_pl
    d/dt(glp1)            <- basal_glp + intest_glp_prod - glp_inact_dpp4 - glp_inact_other
    d/dt(glp1_inactive)   <- glp_inact_dpp4 + glp_inact_other - inactive_glp_cl
    d/dt(glucagon)        <- glucagon_secr - glucagon_clear
    d/dt(glucose_delay30) <- (conc_glucose - glucose_delay30) / 43200
    d/dt(glucose_delay90) <- (conc_glucose - glucose_delay90) / 129600
    d/dt(proinsulin)      <- proinsulin_secr - proinsulin_conv - proinsulin_clear
    d/dt(insulin_musc)    <- ins_equil
    d/dt(proinsulin_musc) <- pro_equil

    # ===================================================================
    # Observations. No residual error is reported: Trujillo 2025 is a
    # deterministic hand-built QSP model with no estimation step.
    # ===================================================================
    glucosePlasma    <- conc_glucose                                                      # mg/dL
    insulinPlasma    <- 1000 * conc_insulin                                               # pM
    proinsulinPlasma <- 1000 * conc_proinsulin                                            # pM
    glucagonPlasma   <- glucagon                                                          # pmol
    hba1c <- 1.4 + (0.7 * glucose_delay30 + 0.3 * glucose_delay90) / 28                   # Table S2 A1c (percent)
  })
}
