Landersdorfer_2012_vildagliptin_qsp <- function() {
  description <- paste(
    "QSP. Mechanism-based glucose-insulin-incretin systems model for",
    "vildagliptin in patients with type 2 diabetes (Landersdorfer 2012,",
    "BJCP 73:373-390). Couples the vildagliptin TMDD PK / DPP-4 layer of the",
    "companion paper (BJCP 73:391-401, already in the library as",
    "Landersdorfer_2012_vildagliptin) to three downstream endpoints: active",
    "GLP-1, plasma glucose and plasma insulin. Meals enter as zero-order",
    "input into four meal-specific gut compartments (breakfast, lunch,",
    "dinner, snack), each with its own absorption rate constant and",
    "glucose-bioavailability factor. Nutrient in the gut stimulates active",
    "GLP-1 secretion; active GLP-1 is cleared by free plasma DPP-4 (which",
    "vildagliptin inhibits) plus a non-saturable pathway. Glucose and insulin",
    "form a reciprocal feedback loop -- glucose stimulates insulin secretion",
    "and insulin stimulates glucose utilization -- and active GLP-1 amplifies",
    "both limbs. All PK and PD parameters were estimated simultaneously in one",
    "S-ADAPT MC-PEM run, so the whole system is one model. Five outputs:",
    "vildagliptin plasma concentration (Cc, ng/mL), DPP-4 activity (DPP4,",
    "mU/mL/min), active GLP-1 (Cglp, pmol/L), glucose (Cglc, mg/dL) and",
    "insulin (Cins, mIU/L). Two values are NOT published and were recovered",
    "from the paper's own figures -- see the vignette Errata before use.",
    sep = " "
  )
  reference <- paste(
    "Landersdorfer CB, He YL, Jusko WJ.",
    "Mechanism-based population modelling of the effects of vildagliptin on",
    "GLP-1, glucose and insulin in patients with type 2 diabetes.",
    "Br J Clin Pharmacol. 2012;73(3):373-390.",
    "doi:10.1111/j.1365-2125.2011.04109.x. PMCID: PMC3370342.",
    "The vildagliptin PK and DPP-4 layer is the companion report:",
    "Landersdorfer CB, He YL, Jusko WJ. Mechanism-based population",
    "pharmacokinetic modelling in diabetes: vildagliptin as a tight binding",
    "inhibitor and substrate of dipeptidyl peptidase IV.",
    "Br J Clin Pharmacol. 2012;73(3):391-401.",
    "doi:10.1111/j.1365-2125.2011.04108.x; see",
    "modellib('Landersdorfer_2012_vildagliptin').",
    sep = " "
  )
  vignette <- "Landersdorfer_2012_vildagliptin_qsp"

  units <- list(
    time = "h",
    dosing = paste(
      "Two kinds of dose record. (1) Vildagliptin: amt in mg into the depot",
      "compartment (converted internally to nmol via MW 303.40 g/mol and",
      "multiplied by F = 0.772). (2) Meals: amt in GRAMS of glucose into one",
      "of glucose_gut_breakfast / glucose_gut_lunch / glucose_gut_dinner /",
      "glucose_gut_snack, given as a zero-order infusion whose duration is the",
      "recorded duration of food intake (dur = tk0 on the event record). The",
      "paper's arbitrary meal anchor is 75000 mg = 75 g per meal, so use",
      "amt = 75; the meal-specific bioavailability factor F is applied inside",
      "the model via f(<gut compartment>).",
      sep = " "
    ),
    concentration = paste(
      "Cc = vildagliptin in plasma (ng/mL); DPP4 = DPP-4 activity in plasma",
      "(mU/mL/min); Cglp = active GLP-1 (pmol/L); Cglc = plasma glucose",
      "(mg/dL); Cins = plasma insulin (mIU/L).",
      sep = " "
    )
  )

  # buildModelDb() infers the registry's `dosing` column from a two-name
  # heuristic that recognises only `depot` and `central`. This model also takes
  # meal doses into four paper-specific gut compartments, so without this
  # declaration the registry would report a truthy-but-incomplete
  # "depot,central" and hide the meal-input route entirely. `central` is
  # retained because the TMDD structure accepts an i.v. vildagliptin dose (the
  # companion report's F was fixed from a model that included i.v. data).
  dosing <- c(
    "depot", "central",
    "glucose_gut_breakfast", "glucose_gut_lunch",
    "glucose_gut_dinner", "glucose_gut_snack"
  )

  paper_specific_compartments <- c(
    "complex_peripheral",
    "glucose_gut_breakfast", "glucose_gut_lunch",
    "glucose_gut_dinner", "glucose_gut_snack",
    "glp1", "glucose", "insulin"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The four glucose_gut_* states hold GRAMS of glucose
  # (the paper's Input = 75000 mg / tk0 anchor, carried in grams so that the
  # published S1 in g^-1 applies directly). The glp1 / glucose / insulin
  # states hold CONCENTRATIONS, not amounts, exactly as the paper writes them
  # (dCglp/dt, dCglc/dt, dCins/dt).
  compartmentData <- list(
    depot                 = list(analyte = "vildagliptin", units = "nmol", specimen = "administration site", verified = TRUE),
    transit1              = list(analyte = "vildagliptin", units = "nmol", specimen = "administration site", verified = TRUE),
    central               = list(analyte = "vildagliptin", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1           = list(analyte = "vildagliptin", units = "nmol", specimen = "tissue", verified = TRUE),
    complex               = list(analyte = "vildagliptin-DPP-4 complex", units = "nmol", specimen = "plasma", verified = TRUE),
    complex_peripheral    = list(analyte = "vildagliptin-DPP-4 complex", units = "nmol", specimen = "tissue", verified = TRUE),
    glucose_gut_breakfast = list(analyte = "glucose", units = "g", specimen = "administration site", verified = TRUE),
    glucose_gut_lunch     = list(analyte = "glucose", units = "g", specimen = "administration site", verified = TRUE),
    glucose_gut_dinner    = list(analyte = "glucose", units = "g", specimen = "administration site", verified = TRUE),
    glucose_gut_snack     = list(analyte = "glucose", units = "g", specimen = "administration site", verified = TRUE),
    glp1                  = list(analyte = "active GLP-1", units = "pmol/L", specimen = "plasma", verified = TRUE),
    glucose               = list(analyte = "glucose", units = "mg/dL", specimen = "plasma", verified = TRUE),
    insulin               = list(analyte = "insulin", units = "mIU/L", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 13L,
    n_studies      = 1L,
    age_range      = "37-64 years",
    age_median     = "53.5 years (mean)",
    weight_range   = "65-116 kg",
    weight_median  = "91 kg (mean)",
    height_range   = "148-183 cm (mean 166 cm)",
    sex_female_pct = 53.8,
    disease_state  = paste(
      "Type 2 diabetes mellitus diagnosed at least 3 months before screening;",
      "washout from hypoglycaemic drugs of up to 4 weeks before dosing.",
      sep = " "
    ),
    dose_range     = paste(
      "Oral vildagliptin 10, 25 or 100 mg twice daily for 28 days plus placebo,",
      "in a randomized, double-blind, placebo-controlled four-way crossover.",
      sep = " "
    ),
    co_medication  = "None (hypoglycaemic drugs washed out before each period).",
    notes          = paste(
      "Demographics from the Results section, first paragraph. Twelve subjects",
      "completed all four periods; one completed only the 10 and 25 mg",
      "treatments (6 male, 7 female). Patients were confined from the evening",
      "of day 26 to the morning of day 29 of each period and received a",
      "standard diet with identical meals across all four treatments (55%",
      "carbohydrate, 25% fat, 20% protein). GLP-1, glucose and insulin were",
      "sampled on day 28. Breakfast and dinner were consumed approximately",
      "30 min after the morning and evening doses; the evening dose was at",
      "10.5 h (from the sampling schedule). Lunch and snack times and the",
      "individual durations of food intake (tk0) are NOT reported -- see the",
      "vignette Assumptions and deviations.",
      sep = " "
    )
  )

  ini({
    # =================================================================
    # PK / DPP-4 layer -- companion article Table 1 (BJCP 73:391-401).
    # These were estimated jointly with the PD layer below ("In the full
    # PK/PD model all PK (shown in the companion report [7]) and PD
    # parameters were estimated at the same time", Methods p. 376), so
    # they belong in this file rather than being imported at run time.
    # Identical to inst/modeldb/specificDrugs/Landersdorfer_2012_vildagliptin.R.
    # =================================================================
    lka1    <- log(1.26)         ; label("Gut1 -> Gut2 absorption rate constant k_a1 (1/h)")           # Companion Table 1: k_a1 = 1.26 (BSV 46%, SE 15%)
    lka2    <- log(1.05)         ; label("Gut2 -> central absorption rate constant k_a2 (1/h)")        # Companion Table 1: k_a2 = 1.05 (BSV 14%, SE 4%)
    lcl     <- log(36.4)         ; label("Non-saturable vildagliptin clearance CL (L/h)")              # Companion Table 1: CL = 36.4 (BSV 25%, SE 9%)
    lvc     <- log(22.2)         ; label("Central volume of distribution V_C (L)")                     # Companion Table 1: V_C = 22.2 (BSV 37%, SE 11%)
    lvp     <- log(97.3)         ; label("Peripheral (tissue) volume V_P (L)")                         # Companion Table 1: V_P = 97.3 (BSV 37%, SE 13%)
    lq      <- log(40.1)         ; label("Inter-compartmental clearance CL_ic (L/h)")                  # Companion Table 1: CL_ic = 40.1 (BSV 34%, SE 11%)
    lfdepot <- fixed(log(0.772)) ; label("Oral bioavailability F (fraction)")                          # Companion Table 1 footnote *: F = 77.2%, fixed from the model including i.v. data
    lkd     <- log(71.9)         ; label("Equilibrium dissociation constant K_d (nmol/L)")             # Companion Table 1: K_d = 71.9 (BSV 54%, SE 16%)
    lk2     <- log(23.4)         ; label("Conversion rate weak -> high-affinity complex k_2 (1/h)")    # Companion Table 1: k_2 = 23.4 (BSV 70%, SE 22%)
    lkoff   <- log(0.612)        ; label("Dissociation rate of vildagliptin from DPP-4 k_off (1/h)")   # Companion Table 1: k_off = 0.612 (BSV 94%, SE 27%)
    lkdeg   <- log(0.110)        ; label("Hydrolysis rate of vildagliptin by DPP-4 k_deg (1/h)")       # Companion Table 1: k_deg = 0.110 (BSV 81%, SE 26%)
    lrmaxc  <- log(5.0)          ; label("Total DPP-4 in the central compartment R_maxC (nmol)")       # Companion Table 1: R_maxC = 5.0 (BSV 12%, SE 4%)
    lrmaxp  <- log(13000)        ; label("Total DPP-4 in the peripheral compartment R_maxP (nmol)")    # Companion Table 1 prints "R_maxP (mmol) 13.0"; the SI prefix is a typeset loss. 13 umol / 5 nmol = 2600-fold matches that paper's own ">2000-fold" statement, whereas 13 mmol would be 2.6e6-fold (and ~1.4 kg of a 110 kDa enzyme). Encoded as 13000 nmol = 13 umol.
    lcf1    <- log(2.80)         ; label("DPP-4 amount to activity conversion cf1 (mU/mL/min per nmol)") # Companion Table 1: cf1 = 2.80 (BSV 17%, SE 5%)

    # =================================================================
    # PD layer -- this paper, Table 1 (p. 380). BSV is reported by
    # S-ADAPT as sqrt(variance) ("S-ADAPT estimates the BSV as variance.
    # The square root of the variance is reported for BSV", Methods
    # p. 376), so omega^2 = (BSV/100)^2 directly -- no CV-to-variance
    # conversion.
    # =================================================================

    # --- Active GLP-1 ---
    lkout_glp_lin <- log(2.07)   ; label("Non-saturable elimination rate constant of active GLP-1 k_out_glp_lin (1/h)")  # Table 1: k_out_glp_lin = 2.07 (BSV 100%, SE 25%)
    ls1           <- log(0.049)  ; label("Stimulation of GLP-1 production by gut nutrient S1 (1/g)")                      # Table 1: S1 = 0.049 g^-1 (BSV 25%, SE 8%). The body text (p. 376) prints "S1 (mg-1)"; Table 1 is correct -- see vignette Errata.
    lcf2          <- log(0.641)  ; label("Free-DPP-4 amount to GLP-1 elimination rate conversion cf2 (1/h/nmol)")         # Table 1: cf2 = 0.641 (BSV 89%, SE 25%)
    lbl_glp       <- log(1.68)   ; label("Baseline active GLP-1 B_glp (pmol/L)")                                          # Table 1: B_glp = 1.68 (BSV 67%, SE 17%)

    # --- Meal glucose input: bioavailability factors and absorption rates ---
    lfdepot_breakfast <- log(0.796) ; label("Glucose bioavailability factor for breakfast F_B (fraction of the 75 g anchor)")  # Table 1: F_B = 0.796 (BSV 27%, SE 8%)
    lfdepot_lunch     <- log(0.865) ; label("Glucose bioavailability factor for lunch F_L (fraction of the 75 g anchor)")      # Table 1: F_L = 0.865 (BSV 16%, SE 4%)
    lfdepot_dinner    <- log(0.817) ; label("Glucose bioavailability factor for dinner F_D (fraction of the 75 g anchor)")     # Table 1: F_D = 0.817 (BSV 16%, SE 4%)
    lfdepot_snack     <- log(0.342) ; label("Glucose bioavailability factor for snack F_S (fraction of the 75 g anchor)")      # Table 1: F_S = 0.342 (BSV 56%, SE 13%)
    lka_breakfast     <- log(0.732) ; label("Gut absorption rate constant after breakfast k_aB (1/h)")                         # Table 1: k_aB = 0.732 (BSV 43%, SE 14%)
    lka_lunch         <- log(0.520) ; label("Gut absorption rate constant after lunch k_aL (1/h)")                             # Table 1: k_aL = 0.520 (BSV 44%, SE 13%)
    lka_dinner        <- log(0.252) ; label("Gut absorption rate constant after dinner k_aD (1/h)")                            # Table 1: k_aD = 0.252 (BSV 31%, SE 14%)
    lka_snack         <- log(0.169) ; label("Gut absorption rate constant after snack k_aS (1/h)")                             # Table 1: k_aS = 0.169 (BSV 89%, SE 8%)

    # --- Glucose ---
    lbl_glc   <- log(133)    ; label("Baseline plasma glucose B_glc (mg/dL)")                                # Table 1: B_glc = 133 (BSV 21%, SE 6%)
    lkout_glc <- log(0.334)  ; label("First-order glucose elimination rate constant k_out_glc (1/h)")         # Table 1: k_out_glc = 0.334 (BSV 85%, SE 21%)
    ls5       <- log(0.584)  ; label("Stimulation of glucose utilization by insulin S5 (mL/mIU; printed as l mIU-1)")  # Table 1: S5 = 0.584 (BSV 150%, SE 43%). Table 1 prints "l mIU-1" but the prefix is a typeset loss: Figure 7B plots ST_ins = S5*[1+S4*(Cglp-Bglp)] on a 0-0.009 axis, and at baseline GLP-1 ST_ins reduces to S5 exactly, so S5 = 5.84e-4 L/mIU = 0.584 mL/mIU. The printed 0.584 L/mIU is 640x above the top of the paper's own axis. Converted explicitly in model(). See vignette Errata.
    ls4       <- log(1.90)   ; label("Amplification of insulin-dependent glucose utilization by GLP-1 S4 (L/pmol)")    # Table 1: S4 = 1.90 (BSV 46%, SE 9%)

    # V_glc is required by the glucose ODE but is reported NOWHERE in the
    # paper, the companion article, or any supplement -- the Methods define
    # it ("where Vglc (dl) is the volume of distribution of glucose", p. 376)
    # and Table 1 omits it. Back-solved from the Figure 3 placebo-panel median
    # (operator ruling 2026-09-02). Digitising that median and fitting V_glc
    # gives 250 dL (95 combinations of the unreported meal times and food-
    # intake durations bracket it at 246-255 dL); the resulting endogenous
    # glucose production B_glc * k_out_glc * V_glc = 11.1 g/h = 2.0 mg/kg/min
    # at the cohort's mean 91 kg, which is the textbook value and was not used
    # in the fit. See vignette Errata.
    lv_glc <- fixed(log(250)) ; label("Apparent volume of distribution of glucose V_glc (dL; not published, back-solved from Figure 3)")

    # --- Insulin ---
    lbl_ins   <- log(9.75)   ; label("Baseline plasma insulin B_ins (mIU/L)")                                 # Table 1: B_ins = 9.75 (BSV 80%, SE 24%)
    lkout_ins <- log(14.0)   ; label("First-order insulin elimination rate constant k_out_ins (1/h)")          # Table 1: k_out_ins = 14.0 (BSV 98%, SE 31%)
    ls3       <- log(0.0185) ; label("Stimulation of insulin secretion by glucose S3 (dL/mg)")                 # Table 1: S3 = 0.0185 (BSV 135%, SE 41%)
    ls2       <- log(0.0701) ; label("Amplification of glucose-dependent insulin secretion by GLP-1 S2 (L/pmol)")  # Table 1: S2 = 0.0701 (BSV 188%, SE 55%)

    # =================================================================
    # Between-subject variability. Diagonal only.
    #
    # Both papers state that a FULL variance-covariance matrix was
    # estimated ("a full variance-covariance matrix for the PD parameters
    # was included. A full variance-covariance matrix was also implemented
    # for the PK parameters. No covariance was included between PK and PD
    # parameters", Methods p. 376), but only the diagonal is published as
    # BSV (%) in either Table 1. The off-diagonals are therefore omitted;
    # see vignette Assumptions and deviations. Between-occasion
    # variability was explicitly NOT included by the authors.
    # =================================================================
    etalka1   ~ 0.2116     # Companion Table 1: BSV k_a1 46%   -> 0.46^2
    etalka2   ~ 0.0196     # Companion Table 1: BSV k_a2 14%
    etalcl    ~ 0.0625     # Companion Table 1: BSV CL 25%
    etalvc    ~ 0.1369     # Companion Table 1: BSV V_C 37%
    etalvp    ~ 0.1369     # Companion Table 1: BSV V_P 37%
    etalq     ~ 0.1156     # Companion Table 1: BSV CL_ic 34%
    etalkd    ~ 0.2916     # Companion Table 1: BSV K_d 54%
    etalk2    ~ 0.4900     # Companion Table 1: BSV k_2 70%
    etalkoff  ~ 0.8836     # Companion Table 1: BSV k_off 94%
    etalkdeg  ~ 0.6561     # Companion Table 1: BSV k_deg 81%
    etalrmaxc ~ 0.0144     # Companion Table 1: BSV R_maxC 12%
    etalrmaxp ~ 0.4096     # Companion Table 1: BSV R_maxP 64%
    etalcf1   ~ 0.0289     # Companion Table 1: BSV cf1 17%

    etalkout_glp_lin ~ 1.0000     # Table 1: BSV k_out_glp_lin 100%
    etals1           ~ 0.0625     # Table 1: BSV S1 25%
    etalcf2          ~ 0.7921     # Table 1: BSV cf2 89%
    etalbl_glp       ~ 0.4489     # Table 1: BSV B_glp 67%

    etalfdepot_breakfast ~ 0.0729 # Table 1: BSV F_B 27%
    etalfdepot_lunch     ~ 0.0256 # Table 1: BSV F_L 16%
    etalfdepot_dinner    ~ 0.0256 # Table 1: BSV F_D 16%
    etalfdepot_snack     ~ 0.3136 # Table 1: BSV F_S 56%
    etalka_breakfast     ~ 0.1849 # Table 1: BSV k_aB 43%
    etalka_lunch         ~ 0.1936 # Table 1: BSV k_aL 44%
    etalka_dinner        ~ 0.0961 # Table 1: BSV k_aD 31%
    etalka_snack         ~ 0.7921 # Table 1: BSV k_aS 89%

    etalbl_glc   ~ 0.0441         # Table 1: BSV B_glc 21%
    etalkout_glc ~ 0.7225         # Table 1: BSV k_out_glc 85%
    etals5       ~ 2.2500         # Table 1: BSV S5 150%
    etals4       ~ 0.2116         # Table 1: BSV S4 46%

    etalbl_ins   ~ 0.6400         # Table 1: BSV B_ins 80%
    etalkout_ins ~ 0.9604         # Table 1: BSV k_out_ins 98%
    etals3       ~ 1.8225         # Table 1: BSV S3 135%
    etals2       ~ 3.5344         # Table 1: BSV S2 188%

    # =================================================================
    # Residual error. Both papers used a combined additive plus
    # proportional error model on every output ("The residual
    # unidentified variability was described by a combined additive and
    # proportional error model for active GLP-1, glucose and insulin
    # concentrations", Methods p. 376).
    # =================================================================
    propSd      <- 0.487  ; label("Vildagliptin proportional residual SD (fraction)")           # Companion Table 1: CV_Vilda = 48.7%
    addSd       <- 0.99   ; label("Vildagliptin additive residual SD (ng/mL)")                  # Companion Table 1: SD_Vilda = 0.99
    propSd_DPP4 <- 0.196  ; label("DPP-4 activity proportional residual SD (fraction)")         # Companion Table 1: CV_DPP-4 = 19.6%
    addSd_DPP4  <- 0.061  ; label("DPP-4 activity additive residual SD (mU/mL/min)")            # Companion Table 1: SD_DPP-4 = 0.061
    propSd_Cglp <- 0.344  ; label("Active GLP-1 proportional residual SD (fraction)")           # Table 1: CV_glp = 34.4%
    addSd_Cglp  <- 2.01   ; label("Active GLP-1 additive residual SD (pmol/L)")                 # Table 1: SD_glp = 2.01
    propSd_Cglc <- 0.169  ; label("Glucose proportional residual SD (fraction)")                # Table 1: CV_glc = 16.9%
    addSd_Cglc  <- 4.06   ; label("Glucose additive residual SD (mg/dL)")                       # Table 1: SD_glc = 4.06
    propSd_Cins <- 0.298  ; label("Insulin proportional residual SD (fraction)")                # Table 1: CV_ins = 29.8%
    addSd_Cins  <- 1.05   ; label("Insulin additive residual SD (mIU/L)")                       # Table 1: SD_ins = 1.05
  })

  model({
    # -----------------------------------------------------------------
    # Constants
    # -----------------------------------------------------------------
    # Vildagliptin molecular weight (PubChem CID 6918537, C17H25N3O2).
    # Converts the user-facing mg dose to internal nmol amounts and the
    # internal nM concentration back to ng/mL for the bioanalytical readout.
    mw_vilda <- 303.40
    # Meal gut amounts are carried in grams so that the published S1 (1/g)
    # applies directly to Glc_Gut; the glucose ODE needs mg/dL/h, hence this
    # g -> mg conversion on the absorption rate.
    mg_per_g <- 1000
    # S5 is printed in Table 1 as "l mIU-1" but reproduces the paper's own
    # Figures 3 and 4 only when it acts as mL/mIU. Applied here as an explicit
    # unit conversion so the printed 0.584 stays visible in ini().
    ml_per_l <- 1000

    # -----------------------------------------------------------------
    # Individual parameters
    # -----------------------------------------------------------------
    ka1   <- exp(lka1   + etalka1)
    ka2   <- exp(lka2   + etalka2)
    cl    <- exp(lcl    + etalcl)
    vc    <- exp(lvc    + etalvc)
    vp    <- exp(lvp    + etalvp)
    q     <- exp(lq     + etalq)
    fdep  <- exp(lfdepot)
    kd    <- exp(lkd    + etalkd)
    k2    <- exp(lk2    + etalk2)
    koff  <- exp(lkoff  + etalkoff)
    kdegv <- exp(lkdeg  + etalkdeg)
    rmaxc <- exp(lrmaxc + etalrmaxc)
    rmaxp <- exp(lrmaxp + etalrmaxp)
    cf1   <- exp(lcf1   + etalcf1)

    kout_glp_lin <- exp(lkout_glp_lin + etalkout_glp_lin)
    s1           <- exp(ls1           + etals1)
    cf2          <- exp(lcf2          + etalcf2)
    bl_glp       <- exp(lbl_glp       + etalbl_glp)

    f_breakfast  <- exp(lfdepot_breakfast + etalfdepot_breakfast)
    f_lunch      <- exp(lfdepot_lunch     + etalfdepot_lunch)
    f_dinner     <- exp(lfdepot_dinner    + etalfdepot_dinner)
    f_snack      <- exp(lfdepot_snack     + etalfdepot_snack)
    ka_breakfast <- exp(lka_breakfast     + etalka_breakfast)
    ka_lunch     <- exp(lka_lunch         + etalka_lunch)
    ka_dinner    <- exp(lka_dinner        + etalka_dinner)
    ka_snack     <- exp(lka_snack         + etalka_snack)

    bl_glc   <- exp(lbl_glc   + etalbl_glc)
    kout_glc <- exp(lkout_glc + etalkout_glc)
    s5       <- exp(ls5       + etals5)
    s4       <- exp(ls4       + etals4)
    v_glc    <- exp(lv_glc)

    bl_ins   <- exp(lbl_ins   + etalbl_ins)
    kout_ins <- exp(lkout_ins + etalkout_ins)
    s3       <- exp(ls3       + etals3)
    s2       <- exp(ls2       + etals2)

    # -----------------------------------------------------------------
    # Vildagliptin TMDD PK and DPP-4 activity (companion article).
    # Concentrations in nmol/L = nM; amounts in nmol; volumes in L.
    # -----------------------------------------------------------------
    Cc_nM <- central     / vc
    Cp_nM <- peripheral1 / vp

    # Free DPP-4 amount (nmol) in each compartment: (R_max - DR).
    freeR_c <- rmaxc - complex
    freeR_p <- rmaxp - complex_peripheral

    # Slow tight-binding rate, V_max = k_2 * C * (R_max - DR) / (K_d + C),
    # in nmol/h of complex formed.
    v_maxc <- k2 * Cc_nM * freeR_c / (kd + Cc_nM)
    v_maxp <- k2 * Cp_nM * freeR_p / (kd + Cp_nM)

    d/dt(depot)              <- -ka1 * depot
    d/dt(transit1)           <-  ka1 * depot - ka2 * transit1
    d/dt(central)            <-  ka2 * transit1 - (cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1 - v_maxc + koff * complex
    d/dt(peripheral1)        <-  (q / vc) * central - (q / vp) * peripheral1 - v_maxp + koff * complex_peripheral
    d/dt(complex)            <-  v_maxc - koff * complex - kdegv * complex
    d/dt(complex_peripheral) <-  v_maxp - koff * complex_peripheral - kdegv * complex_peripheral

    # Dose in mg -> bioavailable nmol. 100 mg * 0.772 * 1e6 / 303.40 = 254452 nmol.
    f(depot) <- fdep * 1e6 / mw_vilda

    # -----------------------------------------------------------------
    # Meal glucose in the gut (this paper, Methods p. 375).
    #   dA_GB/dt = Input * F_B - k_aB * A_GB    (and likewise L, D, S)
    # Input = 75000 mg glucose / (individual duration of food intake), an
    # arbitrary anchor per the paper. The zero-order input is supplied by the
    # event record (amt = 75 g with dur = tk0); f() applies the meal-specific
    # bioavailability factor, which scales the infusion RATE and leaves the
    # duration alone -- exactly "Input x F". All initial conditions are zero.
    # -----------------------------------------------------------------
    d/dt(glucose_gut_breakfast) <- -ka_breakfast * glucose_gut_breakfast
    d/dt(glucose_gut_lunch)     <- -ka_lunch     * glucose_gut_lunch
    d/dt(glucose_gut_dinner)    <- -ka_dinner    * glucose_gut_dinner
    d/dt(glucose_gut_snack)     <- -ka_snack     * glucose_gut_snack
    f(glucose_gut_breakfast) <- f_breakfast
    f(glucose_gut_lunch)     <- f_lunch
    f(glucose_gut_dinner)    <- f_dinner
    f(glucose_gut_snack)     <- f_snack

    # Total nutrient in the gut (g), Glc_Gut = A_GB + A_GL + A_GD + A_GS.
    glc_gut <- glucose_gut_breakfast + glucose_gut_lunch +
      glucose_gut_dinner + glucose_gut_snack

    # Glucose absorption rate (mg/dL/h):
    #   Glc_GutAb = (k_aB*A_GB + k_aL*A_GL + k_aD*A_GD + k_aS*A_GS) / V_glc
    glc_gut_ab <- mg_per_g *
      (ka_breakfast * glucose_gut_breakfast + ka_lunch  * glucose_gut_lunch +
         ka_dinner  * glucose_gut_dinner    + ka_snack  * glucose_gut_snack) / v_glc

    # -----------------------------------------------------------------
    # Active GLP-1 (pmol/L). Secretion is stimulated by gut nutrient;
    # elimination is the sum of a non-saturable pathway and a term
    # proportional to the free plasma DPP-4 amount, which vildagliptin
    # depletes. Steady state: k_in_glp = B_glp * (k_out_glp_lin + R_maxC*cf2).
    # -----------------------------------------------------------------
    kin_glp <- bl_glp * (kout_glp_lin + rmaxc * cf2)
    d/dt(glp1) <- kin_glp * (1 + s1 * glc_gut) -
      (kout_glp_lin + freeR_c * cf2) * glp1
    glp1(0) <- bl_glp

    # -----------------------------------------------------------------
    # Glucose (mg/dL) and insulin (mIU/L) with reciprocal feedback, both
    # limbs amplified by active GLP-1 above its baseline.
    #   ST_ins = S5 * [1 + S4 * (C_glp - B_glp)]     peripheral insulin sensitivity
    #   ST_glc = S3 * [1 + S2 * (C_glp - B_glp)]     pancreatic glucose sensitivity
    # ST_ins carries the mL/mIU -> L/mIU conversion discussed above.
    # -----------------------------------------------------------------
    st_ins <- (s5 / ml_per_l) * (1 + s4 * (glp1 - bl_glp))
    st_glc <- s3 * (1 + s2 * (glp1 - bl_glp))

    kin_glc <- bl_glc * kout_glc
    kin_ins <- bl_ins * kout_ins

    d/dt(glucose) <- kin_glc + glc_gut_ab -
      kout_glc * (1 + st_ins * (insulin - bl_ins)) * glucose
    d/dt(insulin) <- kin_ins * (1 + st_glc * (glucose - bl_glc)) -
      kout_ins * insulin
    glucose(0) <- bl_glc
    insulin(0) <- bl_ins

    # -----------------------------------------------------------------
    # Observations
    # -----------------------------------------------------------------
    Cc   <- Cc_nM * mw_vilda / 1000   # vildagliptin plasma concentration (ng/mL)
    DPP4 <- cf1 * freeR_c             # DPP-4 activity in plasma (mU/mL/min)
    Cglp <- glp1                      # active GLP-1 (pmol/L)
    Cglc <- glucose                   # plasma glucose (mg/dL)
    Cins <- insulin                   # plasma insulin (mIU/L)

    Cc   ~ prop(propSd)      + add(addSd)
    DPP4 ~ prop(propSd_DPP4) + add(addSd_DPP4)
    Cglp ~ prop(propSd_Cglp) + add(addSd_Cglp)
    Cglc ~ prop(propSd_Cglc) + add(addSd_Cglc)
    Cins ~ prop(propSd_Cins) + add(addSd_Cins)
  })
}
