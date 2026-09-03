Riggs_2012_ckd_mbd_qsp <- function() {
  description <- paste(
    "QSP. Multiscale physiology-based systems model of chronic kidney",
    "disease-mineral bone disorder (CKD-MBD) in a hypothetical adult",
    "patient. 31 ODE states spanning four biological scales: organ-level",
    "renal handling of calcium and phosphate driven by a progressively",
    "declining glomerular filtration rate; endocrine control by",
    "parathyroid hormone, 1-alpha-hydroxylase and calcitriol, including",
    "calcium-sensing-receptor feedback on PTH secretion and both pool and",
    "hypertrophy limbs of parathyroid gland capacity; the cellular",
    "RANK / RANK-ligand / osteoprotegerin axis with responding and",
    "fast/slow osteoblast populations, active osteoclasts, latent and",
    "active TGF-beta, and the intracellular RUNX2 / CREB / BCL-2 cascade",
    "governing osteoblast apoptosis; and a clinical-outcome layer linking",
    "the bone-formation and bone-resorption markers to lumbar spine bone",
    "mineral density through an indirect-response relationship. Two",
    "therapeutic interventions are built in and dosable: a hypothetical",
    "calcium-sensing-receptor agonist (calcimimetic), entered as constant",
    "millimolar calcium equivalents added to the CaSR feedback",
    "expression, and exogenous calcitriol, which is assumed to interact",
    "with the system identically to endogenous calcitriol. The model is",
    "deterministic: the authors state that no inter- or intra-individual",
    "variance terms were estimated, so it carries no IIV and no",
    "residual-error model.",
    sep = " "
  )
  reference <- paste(
    "Riggs MM, Peterson MC, Gastonguay MR. Multiscale physiology-based",
    "modeling of mineral bone disorder in patients with impaired kidney",
    "function. J Clin Pharmacol. 2012;52(1 Suppl):45S-53S.",
    "doi:10.1177/0091270011412967.",
    "The GFR forcing function (equation 1), the reiterated PTH-production",
    "and parathyroid-gland equations (equations 2-4), the lumbar spine BMD",
    "indirect-response extension (equation 5), the CaSR agonist term",
    "(equation 6) and the three fitted BMD parameters are from that paper.",
    "The remainder of the calcium-homeostasis / bone-remodeling backbone",
    "originates in Peterson MC, Riggs MM. A physiologically based",
    "mathematical model of integrated calcium homeostasis and bone",
    "remodeling. Bone. 2010;46(1):49-63. doi:10.1016/j.bone.2009.08.053,",
    "which the 2012 paper cites as reference 9 and does not reproduce.",
    "That paper is not open access and could not be obtained, so the",
    "backbone equations and parameter values here were transcribed from",
    "the authors' own open-source release of the same model,",
    "https://github.com/metrumresearchgroup/OpenBoneMin, file",
    "inst/OpenBoneMin.cpp at commit b3d59bfc (2026-01-26). Riggs and Gastonguay",
    "are authors of both that code and this paper, and the code's",
    "constants reproduce every value equation 2-6 prints. Where the two",
    "sources disagree the published paper governs; the differences are",
    "listed in the validation vignette Errata.",
    sep = " "
  )
  vignette <- "Riggs_2012_ckd_mbd_qsp"

  # State names. The source implementation uses single-letter names for many
  # states (P, B, T, R, S, Q, O, N, M, L), several of which are unusable in
  # rxode2 (`T` collides with TRUE, `t` is reserved for time) and none of
  # which is self-describing. Each has been renamed to a descriptive symbol;
  # the mapping to the OpenBoneMin.cpp name is given on every compartmentData
  # entry and reproduced as a table in the validation vignette.
  paper_specific_compartments <- c(
    "PTH", "PTpool", "PThypertrophy", "calcitriol", "alphaOHase",
    "plasmaCa", "ECCPhos", "PhosGut", "IntraPO", "gutCa", "gutCaAbsorp",
    "HAp", "boneCaExch", "boneCaNonExch", "OBfast", "OBslow", "OC", "ROB",
    "TGFB", "TGFBact", "RANKL", "RANK", "RANK_RANKL", "OPG", "OPG_RANKL",
    "RX2", "CREB", "BCL2", "BMDls", "urineCa", "casrAgonist"
  )

  # buildModelDb() infers the registry's dosing column only from states
  # literally named "depot" or "central"; neither exists here. Exogenous
  # calcitriol is dosed straight into the endogenous calcitriol pool (the
  # paper assumes the two are indistinguishable to the system), and the
  # hypothetical calcimimetic is dosed into a zero-derivative state so a
  # bolus sets, and holds, the constant millimolar calcium equivalent the
  # paper describes.
  dosing <- c("calcitriol", "casrAgonist")

  units <- list(
    time          = "h",
    dosing        = "pmol (calcitriol) or mmol/L calcium equivalents (CaSR agonist)",
    concentration = "pmol/L (PTH, calcitriol) or mmol/L (calcium, phosphate)"
  )

  # analyte / units / specimen are `verified = TRUE` where the source states
  # the moiety and its units explicitly -- for the six states the paper
  # itself writes out (equations 2-5) that is the paper; for the rest it is
  # the $INIT and [ODE] blocks of OpenBoneMin.cpp, whose comments name the
  # moiety and unit of every state. `verified = FALSE` marks the states whose
  # biological matrix had to be inferred because neither source names one.
  compartmentData <- list(
    PTH            = list(analyte = "parathyroid hormone", units = "pmol", specimen = "plasma", verified = TRUE),                      # obm PTH
    PTpool         = list(analyte = "parathyroid gland PTH-production capacity", units = "fraction of maximal", specimen = "tissue", verified = TRUE),  # obm S
    PThypertrophy  = list(analyte = "parathyroid gland hypertrophy factor", units = "ratio to baseline", specimen = "tissue", verified = TRUE),         # obm PTmax
    calcitriol     = list(analyte = "calcitriol (1,25-dihydroxyvitamin D3)", units = "pmol", specimen = "plasma", verified = TRUE),    # obm B
    alphaOHase     = list(analyte = "renal 1-alpha-hydroxylase", units = "pmol/h equivalent", specimen = "tissue", verified = TRUE),   # obm AOH
    plasmaCa       = list(analyte = "calcium", units = "mmol", specimen = "plasma", verified = TRUE),                                  # obm P
    ECCPhos        = list(analyte = "inorganic phosphate", units = "mmol", specimen = "plasma", verified = TRUE),                      # obm ECCPhos
    PhosGut        = list(analyte = "dietary inorganic phosphate", units = "mmol", specimen = "administration site", verified = TRUE), # obm PhosGut
    IntraPO        = list(analyte = "intracellular inorganic phosphate", units = "mmol", specimen = "tissue", verified = TRUE),        # obm IntraPO
    gutCa          = list(analyte = "dietary calcium", units = "mmol", specimen = "administration site", verified = TRUE),             # obm T
    gutCaAbsorp    = list(analyte = "calcitriol-dependent intestinal calcium-absorption capacity", units = "fraction of maximal", specimen = "not applicable", verified = TRUE), # obm R
    HAp            = list(analyte = "bone hydroxyapatite deposition capacity", units = "ratio to baseline", specimen = "tissue", verified = TRUE),      # obm HAp
    boneCaExch     = list(analyte = "immediately exchangeable bone calcium", units = "mmol", specimen = "tissue", verified = TRUE),    # obm Q
    boneCaNonExch  = list(analyte = "non-immediately-exchangeable bone calcium", units = "mmol", specimen = "tissue", verified = TRUE),# obm Qbone
    OBfast         = list(analyte = "fast-turnover osteoblasts", units = "pmol/L", specimen = "tissue", verified = FALSE),             # obm OBfast
    OBslow         = list(analyte = "slow-turnover osteoblasts", units = "pmol/L", specimen = "tissue", verified = FALSE),             # obm OBslow
    OC             = list(analyte = "active osteoclasts", units = "pmol/L", specimen = "tissue", verified = FALSE),                    # obm OC
    ROB            = list(analyte = "responding osteoblasts", units = "pmol/L", specimen = "tissue", verified = FALSE),                # obm ROB1
    TGFB           = list(analyte = "latent transforming growth factor beta", units = "pmol/L", specimen = "tissue", verified = FALSE),# obm TGFB
    TGFBact        = list(analyte = "active transforming growth factor beta", units = "pmol/L", specimen = "tissue", verified = FALSE),# obm TGFBact
    RANKL          = list(analyte = "free RANK ligand", units = "pmol/L", specimen = "tissue", verified = TRUE),                       # obm L
    RANK           = list(analyte = "free RANK receptor", units = "pmol/L", specimen = "tissue", verified = TRUE),                     # obm RNK
    RANK_RANKL     = list(analyte = "RANK-RANK-ligand complex", units = "pmol/L", specimen = "tissue", verified = TRUE),               # obm M
    OPG            = list(analyte = "free osteoprotegerin", units = "pmol/L", specimen = "tissue", verified = TRUE),                   # obm O
    OPG_RANKL      = list(analyte = "osteoprotegerin-RANK-ligand complex", units = "pmol/L", specimen = "tissue", verified = TRUE),    # obm N
    RX2            = list(analyte = "RUNX2 transcription factor", units = "arbitrary", specimen = "tissue", verified = TRUE),          # obm RX2
    CREB           = list(analyte = "cAMP response element-binding protein", units = "arbitrary", specimen = "tissue", verified = TRUE), # obm CREB
    BCL2           = list(analyte = "B-cell lymphoma 2 protein", units = "arbitrary", specimen = "tissue", verified = TRUE),           # obm BCL2
    BMDls          = list(analyte = "lumbar spine bone mineral density", units = "ratio to baseline", specimen = "not applicable", verified = TRUE),  # obm BMDlsDEN
    urineCa        = list(analyte = "cumulative excreted calcium", units = "mmol", specimen = "urine", verified = TRUE),               # obm UCA
    casrAgonist    = list(analyte = "hypothetical calcium-sensing-receptor agonist (calcimimetic)", units = "mmol/L calcium equivalents", specimen = "administration site", verified = TRUE)  # equation (6)
  )

  covariateData <- list()

  population <- list(
    species       = "human",
    n_subjects    = 1,
    n_studies     = 2,
    age_range     = "not reported",
    disease_state = paste(
      "Chronic kidney disease-mineral bone disorder (CKD-MBD). A single",
      "deterministic hypothetical patient starts from normal renal",
      "function (glomerular filtration rate 100 mL/min) and is carried",
      "through a typical 10-year course of progressive renal impairment",
      "by the exponential GFR decline of equation (1), reaching",
      "58 mL/min at month 28, 39 mL/min at month 50 and 16 mL/min at",
      "month 120."
    ),
    dose_range    = paste(
      "Hypothetical calcimimetic: constant 0.25, 0.5, 0.75 and 1 mM",
      "calcium equivalents from year 8.5 (Methods, equation 6); the",
      "Figure 4 legend labels the same simulations 0.33, 0.67 and",
      "1.0 mM Ca Eq. Calcitriol: 1.25 and 2.5 ug every other day from",
      "year 8.5 (Methods; Figure 5)."
    ),
    notes         = paste(
      "There is no fitted cohort. The model is deterministic and was run",
      "as a single hypothetical patient; the authors state explicitly",
      "that variance terms for inter- and intra-individual differences",
      "were not estimated, and that Berkeley Madonna, the fitting",
      "platform, does not provide standard errors. The n_studies count of 2",
      "refers to the two literature data sets the BMD extension was built",
      "against, not to a fitted cohort: the three BMD parameters were fit to",
      "the digitized clinical data of the paper's reference 12 and evaluated",
      "against the CKD bone mineral density data of Rix et al. (reference 13),",
      "which supply the observed symbols in Figures 2 and 3."
    )
  )

  ini({
    # =====================================================================
    # Time unit is hours. One year is taken as 8766 h, the convention used
    # throughout OpenBoneMin.cpp (line 765).
    #
    # PROVENANCE. Source tags used in the trailing comments:
    #   "eq (N)"  -- printed in equation N of Riggs 2012 (pp 47S-48S)
    #   "Results" -- printed in the Results text of Riggs 2012 (p 49S)
    #   "obm:NNN" -- line NNN of inst/OpenBoneMin.cpp in
    #                github.com/metrumresearchgroup/OpenBoneMin, the
    #                authors' open-source release of the upstream
    #                Peterson & Riggs 2010 backbone (see `reference`).
    #
    # FIXED vs ESTIMATED. Only gamOB, gamOCls and koutBMDls were fit in
    # this paper ("The parameters gamma_OB, gamma_OC and k_out,BMD were
    # fit (unpublished data) to reported clinical data", p 49S); they are
    # left estimable. Every other value is a fixed system constant
    # inherited from the upstream model and is wrapped in fixed().
    # =====================================================================

    # ---- Renal function forcing function, equation (1) -------------------
    gfrAsymptote <- fixed(10)   ; label("Asymptotic glomerular filtration rate at t = infinity (mL/min)")  # eq (1)
    gfrDecline   <- fixed(90)   ; label("Decline in glomerular filtration rate from baseline to asymptote (mL/min)")  # eq (1)
    gfrKdecay    <- fixed(0.27) ; label("First-order rate constant of glomerular filtration rate decline (1/year)")   # eq (1)
    hoursPerYear <- fixed(8766) ; label("Hours per year used to convert the equation (1) timescale (h/year)")         # obm:765
    mLminPerLh   <- fixed(16.667) ; label("Conversion from mL/min to L/h (mL/min per L/h)")                            # obm:263, obm:764

    # ---- Distribution volume ---------------------------------------------
    V1 <- fixed(14.0) ; label("Extracellular fluid volume, shared by PTH, calcium, phosphate and calcitriol (L)")  # eq (2) (PTHconc = PTH / 14 L); obm:54

    # ---- Baseline (t = 0) state values, the upstream steady state --------
    # These double as the reference values that many of the model's
    # normalised feedback terms divide by, so they are parameters rather
    # than literals in model().
    PTHinit        <- fixed(53.90)      ; label("Baseline plasma PTH amount (pmol)")                       # eq (2); obm:225
    PTpoolInit     <- fixed(0.50)       ; label("Baseline parathyroid gland pool (fraction of maximal)")   # eq (3); obm:226
    PThyperInit    <- fixed(1.00)       ; label("Baseline parathyroid gland hypertrophy (ratio)")          # eq (4); obm:227
    CtriolInit     <- fixed(1260.0)     ; label("Baseline plasma calcitriol amount (pmol)")                # obm:228
    AOHinitFactor  <- fixed(10.0)       ; label("Divisor giving baseline 1-alpha-hydroxylase from baseline calcitriol (unitless)")  # obm:42
    CaInit         <- fixed(32.90)      ; label("Baseline plasma calcium amount (mmol)")                   # obm:229
    ECCPhosInit    <- fixed(16.8)       ; label("Baseline extracellular phosphate amount (mmol)")          # obm:230
    GutCaInit      <- fixed(1.58471)    ; label("Baseline gut calcium amount (mmol)")                      # obm:231
    GutCaAbsInit   <- fixed(0.50)       ; label("Baseline calcitriol-dependent gut calcium absorption capacity (fraction)")  # obm:232
    HApInit        <- fixed(1.00)       ; label("Baseline hydroxyapatite deposition capacity (ratio)")     # obm:233
    PhosGutInit    <- fixed(0.839)      ; label("Baseline dietary phosphate amount in gut (mmol)")         # obm:234
    IntraPOInit    <- fixed(3226.0)     ; label("Baseline intracellular phosphate amount (mmol)")          # obm:235
    OCinit         <- fixed(0.00115398) ; label("Baseline active osteoclasts (pmol/L)")                    # obm:236
    ROBinit        <- fixed(0.00104122) ; label("Baseline responding osteoblasts (pmol/L)")                # obm:237
    RANKLinit      <- fixed(0.40)       ; label("Baseline free RANK ligand (pmol/L)")                      # obm:238
    RANKinit       <- fixed(10.0)       ; label("Baseline free RANK (pmol/L)")                             # obm:239
    OPGinit        <- fixed(4.0)        ; label("Baseline free osteoprotegerin (pmol/L)")                  # obm:240
    BoneCaExchInit <- fixed(100.0)      ; label("Baseline immediately exchangeable bone calcium (mmol)")   # obm:241
    BoneCaNonExchInit <- fixed(24900.0) ; label("Baseline non-exchangeable bone calcium (mmol)")           # obm:242
    RX2init        <- fixed(10.0)       ; label("Baseline RUNX2 (arbitrary units)")                        # obm:243
    CREBinit       <- fixed(10.0)       ; label("Baseline CREB (arbitrary units)")                         # obm:244
    BCL2init       <- fixed(100.0)      ; label("Baseline BCL-2, the product of baseline RUNX2 and CREB (arbitrary units)")  # obm:245
    BMDlsInit      <- fixed(1.0)        ; label("Baseline lumbar spine bone mineral density (ratio to baseline)")            # obm:260
    OBtot0         <- fixed(0.00501324) ; label("Baseline total osteoblasts (pmol/L)")                     # obm:49
    FracOBfast     <- fixed(0.797629)   ; label("Fraction of baseline osteoblasts in the fast-turnover pool (unitless)")     # obm:154
    Pic0           <- fixed(0.228142)   ; label("Baseline active TGF-beta concentration, also the baseline Pic feedback value (pmol/L)")  # obm:148
    TGFBscale      <- fixed(1000.0)     ; label("Ratio of baseline latent to baseline active TGF-beta (unitless)")           # obm:36, obm:308

    # ---- Calcium fluxes: bone, gut, kidney -------------------------------
    CaDay      <- fixed(88.0)     ; label("Daily calcium turnover between plasma and the exchangeable bone pool (mmol/day)")  # obm:55
    FracJ14    <- fixed(0.107763) ; label("Osteoclast-dependent fraction of calcium flux from bone to plasma (unitless)")     # obm:57
    J14OCmax   <- fixed(0.543488) ; label("Maximal osteoclast-driven bone-to-plasma calcium flux term (1/h)")                 # obm:58
    J14OCgam   <- fixed(1.6971)   ; label("Hill exponent of the osteoclast effect on bone-to-plasma calcium flux (unitless)") # obm:59
    FracJ15    <- fixed(0.114376) ; label("Hydroxyapatite-dependent fraction of calcium flux from plasma to bone (unitless)") # obm:60
    k14a       <- fixed(0.0000244437) ; label("Rate constant for calcium transfer from non-exchangeable to exchangeable bone (1/h)")  # obm:87
    HApMRT     <- fixed(3.60609)  ; label("Mean residence time of the hydroxyapatite deposition-capacity signal (h)")         # obm:88
    CaPratio   <- fixed(0.464)    ; label("Molar ratio of phosphate to calcium in hydroxyapatite (unitless)")                 # obm:300, obm:377
    OralCa     <- fixed(1.002291667) ; label("Oral calcium intake rate, 24.055 mmol/day (mmol/h)")                            # obm:134 (24.055/24)
    T28        <- fixed(0.9)      ; label("Maximal saturable gut-to-plasma calcium transfer rate (1/h)")                      # obm:133
    T310       <- fixed(0.105929) ; label("Fraction of the maximal gut-to-plasma calcium rate attained at baseline gut calcium (unitless)")  # obm:135
    T81        <- fixed(0.75)     ; label("Half-saturation gut calcium amount for calcitriol-dependent absorption (mmol)")    # obm:131
    T87        <- fixed(0.0495)   ; label("First-order passive gut-to-plasma calcium transfer rate constant (1/h)")           # obm:132
    T77        <- fixed(0.909359) ; label("Maximal fraction of an oral calcium dose absorbed (unitless)")                     # obm:136
    T80        <- fixed(4.0)      ; label("Hill exponent of gut calcium-absorption capacity on the absorbed fraction (unitless)")  # obm:137
    T33        <- fixed(0.003)    ; label("Minimal turnover rate constant of gut calcium-absorption capacity (1/h)")          # obm:110
    T34        <- fixed(0.037)    ; label("Maximal turnover rate constant of gut calcium-absorption capacity (1/h)")          # obm:111
    T35        <- fixed(90.0)     ; label("Calcitriol concentration at half-maximal effect on gut calcium-absorption capacity (pmol/L)")  # obm:112
    CaPOgam    <- fixed(1.0)      ; label("Hill exponent of calcitriol on gut calcium-absorption capacity (unitless)")        # obm:113
    Reabs50    <- fixed(1.57322)  ; label("Plasma calcium concentration at half-maximal renal tubular reabsorption (mmol/L)") # obm:105
    T16        <- fixed(1.06147)  ; label("Maximal fold-effect of PTH on renal calcium reabsorption (unitless)")              # obm:82
    T7         <- fixed(2.0)      ; label("Maximal fold-effect of calcitriol on renal calcium excretion (unitless)")          # obm:106
    T9         <- fixed(90.0)     ; label("Calcitriol concentration at half-maximal effect on renal calcium excretion (pmol/L)")  # obm:107
    CaFiltFu   <- fixed(0.6)      ; label("Fraction of plasma calcium unbound and available for glomerular filtration (unitless)")  # obm:447
    CaFiltPTHdep <- fixed(0.5)    ; label("Fraction of filtered calcium reabsorbed by the PTH-dependent route (unitless)")    # obm:447
    ReabsMaxSlope <- fixed(0.3)   ; label("Slope of maximal tubular calcium reabsorption on filtration rate (unitless)")      # obm:453
    ReabsMaxInt <- fixed(0.149997) ; label("Intercept of maximal tubular calcium reabsorption (mmol/h)")                      # obm:453

    # ---- Phosphate --------------------------------------------------------
    OralPhos     <- fixed(0.4375)   ; label("Oral phosphate intake rate, 10.5 mmol/day (mmol/h)")                            # obm:116 (10.5/24)
    FracPhosAbs  <- fixed(0.7)      ; label("Fraction of dietary phosphate absorbed (unitless)")                             # obm:117 (F12)
    T52          <- fixed(0.365)    ; label("First-order rate constant of phosphate absorption from gut (1/h)")              # obm:115
    T49          <- fixed(51.8)     ; label("Rate constant for phosphate movement from extracellular to intracellular pool (L/h)")  # obm:118
    T55          <- fixed(0.019268) ; label("Rate constant for phosphate movement from intracellular to extracellular pool (1/h)")  # obm:119
    T46          <- fixed(1.142)    ; label("Renal tubular phosphate reabsorption threshold concentration (mmol/L)")         # obm:114
    PhosFiltFrac <- fixed(0.88)     ; label("Fraction of plasma phosphate available for glomerular filtration (unitless)")   # obm:487
    PhosBase     <- fixed(1.2)      ; label("Baseline extracellular phosphate concentration (mmol/L)")                       # obm:419, obm:429
    PhosEff0     <- fixed(1.52493)  ; label("Effect of phosphate on 1-alpha-hydroxylase at zero phosphate (unitless)")       # obm:100
    PhosEff50    <- fixed(1.3021)   ; label("Phosphate concentration at half-maximal inhibition of 1-alpha-hydroxylase (mmol/L)")  # obm:101
    PhosEffGam   <- fixed(8.25229)  ; label("Hill exponent of phosphate inhibition of 1-alpha-hydroxylase (unitless)")       # obm:102
    PO4inhPTHgam <- fixed(0.0)      ; label("Hill exponent of the direct phosphate effect on 1-alpha-hydroxylase, not activated in this model (unitless)")  # obm:103

    # ---- PTH, parathyroid gland, calcitriol -------------------------------
    kout        <- fixed(7.142857143) ; label("First-order PTH elimination rate constant, printed as 7.14 in equation (2) (1/h)")  # eq (2); obm:142 (100/14)
    T58         <- fixed(6249.09)  ; label("CaSR feedback on PTH secretion at zero calcium, printed as 6249 in equation (2) (pmol/h)")  # eq (2); obm:143
    T59         <- fixed(11.7387)  ; label("Hill exponent of the CaSR feedback on PTH secretion, printed as 11.7 in equation (2) (unitless)")  # eq (2); obm:144
    T61         <- fixed(96.25)    ; label("CaSR feedback on PTH secretion at infinite calcium; 6249 - 96.25 is the 6153 of equation (2) (pmol/h)")  # eq (2); obm:145
    SPTHbase    <- fixed(385.0)    ; label("PTH secretion rate at baseline, used to back-solve the CaSR half-saturation constant (pmol/h)")  # eq (2) text; obm:568
    T70         <- fixed(0.01)     ; label("Scale of the calcitriol effect on parathyroid gland pool turnover (1/h)")        # eq (3); obm:108
    T71         <- fixed(0.03)     ; label("Scale factor inside the tanh of the parathyroid pool feedback (L/pmol)")         # eq (3); obm:109
    ScaEffGam   <- fixed(0.9)      ; label("Hill exponent of calcium on the parathyroid pool calcitriol set-point (unitless)")  # eq (3); obm:99
    CtriolSetPt <- fixed(90.0)     ; label("Calcitriol set-point of the parathyroid pool feedback at baseline calcium (pmol/L)")  # eq (3); obm:476
    PTout       <- fixed(0.0001604) ; label("First-order turnover rate constant of parathyroid gland hypertrophy (1/h)")     # eq (4); obm:141
    CtriolMax   <- fixed(4.1029)   ; label("Parathyroid hypertrophy drive at zero calcitriol, printed as 4.1 in equation (4) (unitless)")  # eq (4); obm:139
    CtriolMin   <- fixed(0.9)      ; label("Parathyroid hypertrophy drive at infinite calcitriol; 4.1029 - 0.9 is the 3.2 of equation (4) (unitless)")  # eq (4); obm:140
    CtriolPTgam <- fixed(12.5033)  ; label("Hill exponent of calcitriol on parathyroid hypertrophy, printed as 12.5 in equation (4) (unitless)")  # eq (4); obm:138
    T64         <- fixed(0.05)     ; label("First-order elimination rate constant of 1-alpha-hydroxylase (1/h)")             # obm:83
    T65         <- fixed(6.3)      ; label("Maximal 1-alpha-hydroxylase production rate (pmol/h)")                           # obm:84
    T67         <- fixed(1.54865)  ; label("PTH concentration at half-maximal 1-alpha-hydroxylase production (pmol/L)")      # obm:85
    AlphOHgam   <- fixed(0.111241) ; label("Hill exponent of PTH on 1-alpha-hydroxylase production (unitless)")              # obm:86
    T69         <- fixed(0.10)     ; label("First-order calcitriol elimination rate constant (1/h)")                         # obm:104

    # ---- Bone cells: TGF-beta, osteoblast and osteoclast control ----------
    Da          <- fixed(0.02916666667) ; label("Osteoclast and precursor differentiation rate constant, 0.7/day (1/h)")     # obm:63 (0.7/24)
    kb          <- fixed(0.000605516) ; label("Baseline osteoblast elimination rate constant (1/h)")                         # obm:152
    koutTGF0    <- fixed(0.0000298449) ; label("Baseline first-order elimination rate constant of latent TGF-beta (1/h)")    # obm:65
    koutTGFGam  <- fixed(0.919131) ; label("Hill exponent of latent TGF-beta on its own elimination (unitless)")             # obm:66
    OBtgfGAM    <- fixed(0.0111319) ; label("Hill exponent of osteoblasts on latent TGF-beta production (unitless)")         # obm:64
    OCtgfGAM    <- fixed(0.593891) ; label("Hill exponent of osteoclasts on TGF-beta activation (unitless)")                 # obm:67
    EmaxPicROB  <- fixed(3.9745)   ; label("Maximal active-TGF-beta effect on responding-osteoblast production (unitless)")  # obm:68
    PicROBgam   <- fixed(1.80968)  ; label("Hill exponent of active TGF-beta on responding-osteoblast production (unitless)") # obm:69
    FracPicROB  <- fixed(0.883824) ; label("Baseline fraction of the responding-osteoblast TGF-beta effect that is ligand independent (unitless)")  # obm:70
    EmaxPicOB   <- fixed(0.251636) ; label("Maximal active-TGF-beta effect on osteoblast differentiation (unitless)")        # obm:73
    PicOBgam    <- fixed(0.122313) ; label("Hill exponent of active TGF-beta on osteoblast differentiation (unitless)")      # obm:71
    FracPicOB   <- fixed(0.000244818) ; label("Baseline fraction of the osteoblast-differentiation TGF-beta effect that is ligand independent (unitless)")  # obm:72
    EmaxPicOC   <- fixed(1.9746)   ; label("Maximal active-TGF-beta effect on osteoclast differentiation (unitless)")        # obm:77
    PicOCgam    <- fixed(1.0168)   ; label("Hill exponent of active TGF-beta on osteoclast differentiation (unitless)")      # obm:79
    FracPicOC   <- fixed(0.878215) ; label("Baseline fraction of the osteoclast-differentiation TGF-beta effect that is ligand independent (unitless)")  # obm:78
    E0Meff      <- fixed(0.388267) ; label("RANK-RANKL-complex effect on osteoclast formation at zero complex (unitless)")   # obm:74
    EmaxMeffOC  <- fixed(3.15667)  ; label("Maximal RANK-RANKL-complex effect on osteoclast formation (unitless)")           # obm:75
    kinOCgam    <- fixed(8.53065)  ; label("Hill exponent of the RANK-RANKL complex on osteoclast formation (unitless)")     # obm:76
    MOCratioGam <- fixed(0.603754) ; label("Hill exponent of the RANK-RANKL-complex-to-osteoclast ratio on bone calcium release (unitless)")  # obm:62
    PicOBgamkb  <- fixed(2.92375)  ; label("Hill exponent of active TGF-beta on osteoblast apoptosis (unitless)")            # obm:120
    MultPicOBkb <- fixed(3.11842)  ; label("Osteoblast-apoptosis TGF-beta effect at zero active TGF-beta, as a multiple of Pic0 (unitless)")  # obm:121
    FracPic0kb  <- fixed(0.764028) ; label("Osteoblast-apoptosis TGF-beta effect at infinite active TGF-beta, as a fraction of Pic0 (unitless)")  # obm:122
    Frackb      <- fixed(0.313186) ; label("Fraction of the osteoblast apoptosis rate applying to the slow-turnover pool (unitless)")  # obm:130

    # ---- RANK / RANK-ligand / osteoprotegerin axis ------------------------
    k1          <- fixed(0.00000624) ; label("Association rate constant of osteoprotegerin with RANK ligand (L/pmol/h)")     # obm:50
    k2          <- fixed(0.112013)   ; label("Dissociation rate constant of the osteoprotegerin-RANK-ligand complex (1/h)")  # obm:51
    k3          <- fixed(0.00000624) ; label("Association rate constant of RANK with RANK ligand (L/pmol/h)")                # obm:52
    k4          <- fixed(0.112013)   ; label("Dissociation rate constant of the RANK-RANK-ligand complex (1/h)")             # obm:53
    koutL       <- fixed(0.00293273) ; label("First-order elimination rate constant of free RANK ligand (1/h)")              # obm:89
    koutRNK     <- fixed(0.00323667) ; label("First-order elimination rate constant of free RANK (1/h)")                     # obm:61
    kinRNKgam   <- fixed(0.151825)   ; label("Hill exponent of active TGF-beta on RANK production (unitless)")               # obm:60
    koutOPG     <- fixed(15.8885)    ; label("First-order elimination rate constant of free osteoprotegerin (1/h)")          # obm:151 (kO)
    opgPTH50    <- fixed(3.85)       ; label("PTH concentration at half-maximal suppression of osteoprotegerin production (pmol/L)")  # obm:91
    OsteoEffectGam <- fixed(0.173833) ; label("Hill exponent of osteoblasts on RANK-ligand production (unitless)")           # obm:90
    EmaxLpth    <- fixed(1.30721)    ; label("Maximal PTH effect on RANK-ligand production (unitless)")                      # obm:150
    E0RANKL     <- fixed(3.80338)    ; label("Osteoclast survival factor at zero RANK-ligand occupancy (unitless)")          # obm:80
    EmaxL       <- fixed(0.469779)   ; label("Osteoclast survival factor at maximal RANK-ligand occupancy (unitless)")       # obm:81
    LsurvOCgam  <- fixed(3.09023)    ; label("Hill exponent of RANK-ligand occupancy on osteoclast survival (unitless)")     # obm:153
    PiLscale    <- fixed(10.0)       ; label("Scale factor converting the RANK-RANK-ligand complex to the osteoclast-survival driver (unitless)")  # obm:350

    # ---- Intracellular RUNX2 / CREB / BCL-2 cascade -----------------------
    RX2Kout0    <- fixed(0.693)   ; label("Baseline first-order elimination rate constant of RUNX2 (1/h)")                   # obm:92
    E0rx2Kout   <- fixed(0.125)   ; label("RUNX2 elimination rate constant at zero PTH (1/h)")                               # obm:93
    EmaxPTHRX2x <- fixed(5.0)     ; label("Maximal PTH effect on the RUNX2 elimination rate constant (1/h)")                 # obm:94
    E0crebKin   <- fixed(0.5)     ; label("CREB production at zero PTH, as a fraction of baseline (unitless)")               # obm:95
    EmaxPTHcreb <- fixed(3.39745) ; label("Maximal PTH effect on CREB production (unitless)")                                # obm:96
    crebKout    <- fixed(0.00279513) ; label("First-order elimination rate constant of CREB (1/h)")                          # obm:97
    bcl2Kout    <- fixed(0.693)   ; label("First-order elimination rate constant of BCL-2, a 1-h half-life (1/h)")           # obm:98
    RUNX2thresh <- fixed(105.0)   ; label("BCL-2 threshold above which RUNX2 tracks BCL-2 (arbitrary units)")                # obm:517
    RUNX2offset <- fixed(90.0)    ; label("Offset subtracted from BCL-2 to give RUNX2 above the threshold (arbitrary units)") # obm:517
    RUNX20      <- fixed(10.0)    ; label("RUNX2 value below the BCL-2 threshold, and its reference value (arbitrary units)") # obm:129
    RUNkbGAM    <- fixed(3.81644) ; label("Hill exponent of RUNX2 on the osteoblast apoptosis rate constant (unitless)")     # obm:124
    RUNkbMaxFact <- fixed(0.638114) ; label("Maximal RUNX2-driven reduction in the osteoblast apoptosis rate constant, as a fraction (unitless)")  # obm:128
    E0RUNX2kbEffFACT <- fixed(1.01) ; label("Osteoblast apoptosis rate constant at zero RUNX2, as a multiple of kb (unitless)")  # obm:123

    # ---- Lumbar spine BMD, equation (5) ----------------------------------
    # The only three parameters this paper estimated. In OpenBoneMin these
    # are the "DEN" parameter set (koutBMDlsDEN, gamOClsDEN), which carries
    # exactly the values this paper prints; gamOB differs (see vignette
    # Errata) and the published 0.0739 governs.
    koutBMDls <- 0.000145 ; label("First-order lumbar spine BMD elimination rate constant (1/h)")           # Results p 49S; obm:171
    gamOB     <- 0.0739   ; label("Exponent of the relative change in bone-specific alkaline phosphatase on BMD formation (unitless)")  # Results p 49S
    gamOCls   <- 0.0679   ; label("Exponent of the relative change in serum CTx on BMD elimination (unitless)")  # Results p 49S; obm:175
  })

  model({
    # =====================================================================
    # Equation-by-equation transcription. Local names follow
    # OpenBoneMin.cpp so that a reader can audit this block line-by-line
    # against inst/OpenBoneMin.cpp; state names are the descriptive ones
    # documented in compartmentData.
    #
    # Reductions relative to OpenBoneMin.cpp, each behaviourally inert for
    # the scenarios of this paper and justified in the vignette Errata:
    #   - the estrogen / menopause limb is dropped (its master switch ESTON
    #     is 0 by default, which pins EST at 1 and makes every EST term the
    #     identity),
    #   - the denosumab, teriparatide and generic-drug PK compartments are
    #     dropped (all have zero rate constants unless dosed, and none is
    #     used by this paper),
    #   - the femoral-neck BMD state is dropped (this paper reports lumbar
    #     spine only).
    # =====================================================================

    # ---- Renal function, equation (1) ------------------------------------
    # Printed in mL/min; the flux equations below consume it as L/h.
    gfrMLmin <- gfrAsymptote + gfrDecline * exp(-gfrKdecay * t / hoursPerYear)
    gfr <- gfrMLmin / mLminPerLh

    # ---- Baselines derived from the ini() baseline amounts ----------------
    PTHconc0    <- PTHinit / V1
    CaConc0     <- CaInit / V1
    CtriolConc0 <- CtriolInit / V1
    AOHinit     <- CtriolInit / AOHinitFactor          # obm:42
    OBfast0     <- OBtot0 * FracOBfast                 # obm:38
    OBslow0     <- OBtot0 * (1 - FracOBfast)           # obm:39
    OB0         <- OBtot0
    TGFB0       <- Pic0 * TGFBscale                    # obm:36
    TGFBact0    <- Pic0                                # obm:37
    M0          <- k3 * RANKinit * RANKLinit / k4      # obm:40
    N0          <- k1 * OPGinit * RANKLinit / k2       # obm:41

    # ---- Concentrations ---------------------------------------------------
    PTHconc <- PTH / V1
    CaConc  <- plasmaCa / V1
    C2      <- ECCPhos / V1
    C8      <- calcitriol / V1
    Osteoblast <- OBfast + OBslow

    # ---- Calcium flux from bone to plasma (J14) ---------------------------
    T13 <- (CaDay / 24) / BoneCaExchInit                                       # obm:278
    T15 <- CaDay / (CaConc0 * V1 * 24)                                         # obm:280
    J14OC50 <- exp(log((J14OCmax * OCinit^J14OCgam / T13) - OCinit^J14OCgam) / J14OCgam)  # obm:282
    OCeqn <- (J14OCmax * OC^J14OCgam) / (OC^J14OCgam + J14OC50^J14OCgam)       # obm:284
    MOCratio <- RANK_RANKL / OC                                                # obm:288
    MOCratio0 <- M0 / OCinit                                                   # obm:290
    MOCratioEff <- (MOCratio / MOCratio0)^MOCratioGam                          # obm:292
    J14OCdepend <- OCeqn * BoneCaExchInit * FracJ14 * MOCratioEff              # obm:294
    J14 <- T13 * BoneCaExchInit * (1 - FracJ14) + J14OCdepend                  # obm:296
    J41 <- CaPratio * J14                                                      # obm:300

    # ---- Calcium flux from plasma to bone (J15), and slow bone exchange ----
    J15 <- T15 * plasmaCa * (1 - FracJ15) + T15 * plasmaCa * FracJ15 * HAp     # obm:374
    J42 <- CaPratio * J15                                                      # obm:377
    k15a <- k14a * BoneCaNonExchInit / BoneCaExchInit                          # obm:362
    J14a <- k14a * boneCaNonExch                                               # obm:364
    J15a <- k15a * boneCaExch                                                  # obm:366
    kLShap <- 1 / HApMRT                                                       # obm:369
    kHApIn <- kLShap / OB0                                                     # obm:371

    # ---- TGF-beta ---------------------------------------------------------
    kinTGF <- koutTGF0 * TGFB0                                                 # obm:304
    koutTGF <- koutTGF0 * (TGFB / TGFB0)^koutTGFGam                            # obm:306
    koutTGFact <- koutTGF0 * TGFBscale                                         # obm:308
    koutTGFeqn <- koutTGF * TGFB * (OC / OCinit)^OCtgfGAM                      # obm:310

    # ---- TGF-beta effects on osteoblast lineage (Pic terms) ---------------
    bigDb <- kb * OB0 * Pic0 / ROBinit                                         # obm:302
    Dr <- kb * OB0 / Pic0                                                      # obm:318
    E0PicROB <- FracPicROB * Pic0                                              # obm:312
    EC50PicROBparen <- (EmaxPicROB * TGFBact0^PicROBgam / (Pic0 - E0PicROB)) - TGFBact0^PicROBgam  # obm:314
    EC50PicROB <- exp(log(EC50PicROBparen) / PicROBgam)                        # obm:316
    PicROB <- E0PicROB + EmaxPicROB * TGFBact^PicROBgam / (TGFBact^PicROBgam + EC50PicROB^PicROBgam)  # obm:320
    ROBin <- Dr * PicROB                                                       # obm:322
    E0PicOB <- FracPicOB * Pic0                                                # obm:324
    EC50PicOBparen <- (EmaxPicOB * TGFBact0^PicOBgam / (Pic0 - E0PicOB)) - TGFBact0^PicOBgam  # obm:326
    EC50PicOB <- exp(log(EC50PicOBparen) / PicOBgam)                           # obm:328
    PicOB <- E0PicOB + EmaxPicOB * TGFBact^PicOBgam / (TGFBact^PicOBgam + EC50PicOB^PicOBgam)  # obm:330
    KPT <- bigDb / PicOB                                                       # obm:332

    # ---- Osteoclast formation and survival --------------------------------
    EC50MeffOC <- exp(log(M0^kinOCgam * EmaxMeffOC / (1 - E0Meff) - M0^kinOCgam) / kinOCgam)  # obm:334
    MeffOC <- E0Meff + EmaxMeffOC * RANK_RANKL^kinOCgam / (RANK_RANKL^kinOCgam + EC50MeffOC^kinOCgam)  # obm:336
    kinOC2 <- Da * Pic0 * MeffOC * OCinit                                      # obm:338 (PicOCkin is #defined as Pic0)
    E0PicOC <- FracPicOC * Pic0                                                # obm:340
    EC50PicOCparen <- (EmaxPicOC * TGFBact0^PicOCgam / (Pic0 - E0PicOC)) - TGFBact0^PicOCgam  # obm:342
    EC50PicOC <- exp(log(EC50PicOCparen) / PicOCgam)                           # obm:344
    PicOC <- E0PicOC + EmaxPicOC * TGFBact^PicOCgam / (TGFBact^PicOCgam + EC50PicOC^PicOCgam)  # obm:346
    PiL0 <- (k3 / k4) * RANKLinit                                              # obm:348
    PiL <- RANK_RANKL / PiLscale                                               # obm:350
    EC50survInPar <- (E0RANKL - EmaxL) * (PiL0^LsurvOCgam / (E0RANKL - 1)) - PiL0^LsurvOCgam  # obm:352
    EC50surv <- exp(log(EC50survInPar) / LsurvOCgam)                           # obm:354
    LsurvOC <- E0RANKL - (E0RANKL - EmaxL) * (PiL^LsurvOCgam / (PiL^LsurvOCgam + EC50surv^LsurvOCgam))  # obm:356
    KLSoc <- Da * PicOC * LsurvOC                                              # obm:358

    # ---- RANK-ligand and osteoprotegerin production -----------------------
    kinRNK <- (koutRNK * RANKinit + k3 * RANKinit * RANKLinit - k4 * M0) / TGFBact0^kinRNKgam  # obm:286
    kinLbase <- koutL * RANKLinit                                              # obm:381
    OsteoEffect <- (Osteoblast / OB0)^OsteoEffectGam                           # obm:383
    PTH50 <- EmaxLpth * PTHconc0 - PTHconc0                                    # obm:385
    LpthEff <- EmaxLpth * PTHconc / ((PTH50 * OsteoEffect) + PTHconc)          # obm:387
    kinL <- kinLbase * OsteoEffect * LpthEff                                   # obm:389
    pObase <- koutOPG * OPGinit                                                # obm:391
    pO <- pObase * (ROB / ROBinit) * ((PTHconc + (opgPTH50 * (ROB / ROBinit))) / (2 * PTHconc))  # obm:393

    # ---- RUNX2 / CREB / BCL-2 ---------------------------------------------
    RX2Kin <- RX2Kout0 * RX2init                                               # obm:395
    EC50PTHRX2x <- ((EmaxPTHRX2x * PTHconc0) / (RX2Kout0 - E0rx2Kout)) - PTHconc0  # obm:397
    RX2Kout <- E0rx2Kout + EmaxPTHRX2x * PTHconc / (PTHconc + EC50PTHRX2x)     # obm:399
    EC50PTHcreb <- ((EmaxPTHcreb * PTHconc0) / (1 - E0crebKin)) - PTHconc0     # obm:405
    crebKin0 <- crebKout * CREBinit                                            # obm:407
    crebKin <- crebKin0 * (E0crebKin + EmaxPTHcreb * PTHconc / (PTHconc + EC50PTHcreb))  # obm:409

    # ---- Phosphate effect on 1-alpha-hydroxylase, and phosphate fluxes ----
    PO4inhPTH <- (C2 / PhosBase)^PO4inhPTHgam                                  # obm:419
    PhosEffTop <- (PhosEff0 - 1) * (PhosBase^PhosEffGam + PhosEff50^PhosEffGam)  # obm:421
    PhosEffBot <- PhosEff0 * PhosBase^PhosEffGam                               # obm:423
    PhosEffMax <- PhosEffTop / PhosEffBot                                      # obm:425
    PhosEff <- PhosEff0 - (PhosEffMax * PhosEff0 * C2^PhosEffGam / (C2^PhosEffGam + PhosEff50^PhosEffGam))  # obm:427
    if (C2 > PhosBase) {                                                       # obm:429
      PhosEffect <- PhosEff
    } else {
      PhosEffect <- 1
    }
    T66 <- (T67^AlphOHgam + PTHconc0^AlphOHgam) / PTHconc0^AlphOHgam           # obm:360
    T68 <- T66 * PTHconc^AlphOHgam / (T67^AlphOHgam * PO4inhPTH + PTHconc^AlphOHgam)  # obm:431
    SE <- T65 * T68 * PhosEffect                                               # obm:433
    T47 <- T46 * PhosFiltFrac * gfr                                            # obm:487
    J48a <- PhosFiltFrac * gfr * C2 - T47                                      # obm:489
    if (J48a < 0) {                                                            # obm:491
      J48 <- 0
    } else {
      J48 <- J48a
    }
    J53 <- T52 * PhosGut                                                       # obm:494
    J54 <- T49 * C2                                                            # obm:496
    J56 <- T55 * IntraPO                                                       # obm:498

    # ---- Calcitriol-dependent gut calcium absorption ----------------------
    T36 <- T33 + (T34 - T33) * (C8^CaPOgam / (T35^CaPOgam + C8^CaPOgam))       # obm:436
    T37 <- T34 - (T34 - T33) * (C8^CaPOgam / (T35^CaPOgam + C8^CaPOgam))       # obm:437
    T29 <- (T28 * GutCaInit - T310 * GutCaInit) / T310                         # obm:538
    T31 <- T28 * gutCa / (gutCa + T29)                                         # obm:540
    T83 <- gutCaAbsorp / GutCaAbsInit                                          # obm:543
    J40 <- T31 * gutCa * T83 / (gutCa + T81) + T87 * gutCa                     # obm:546
    T85Rpart <- gutCaAbsorp^T80 / (gutCaAbsorp^T80 + T81^T80)                  # obm:550
    T85 <- T77 * T85Rpart                                                      # obm:551

    # ---- Renal calcium handling -------------------------------------------
    CaFilt <- CaFiltFu * CaFiltPTHdep * gfr * CaConc                           # obm:447
    ReabsMax <- (ReabsMaxSlope * gfr * CaConc0 - ReabsMaxInt) * (Reabs50 + CaConc0) / CaConc0  # obm:453 (tmEST = 1 with the estrogen limb off)
    T17 <- PTHconc0 * T16 - PTHconc0                                           # obm:456
    ReabsPTHeff <- (T16 * PTHconc) / (PTHconc + T17)                           # obm:458
    CaReabsActive <- (ReabsMax * CaConc / (Reabs50 + CaConc)) * ReabsPTHeff    # obm:462
    T20 <- CaFilt - CaReabsActive                                              # obm:464
    T10 <- T7 * C8 / (C8 + T9)                                                 # obm:466
    J27a <- (2 - T10) * T20                                                    # obm:469
    if (J27a < 0) {                                                            # obm:472
      J27 <- 0
    } else {
      J27 <- J27a
    }

    # ---- Parathyroid gland pool, equation (3) -----------------------------
    ScaEff <- (CaConc0 / CaConc)^ScaEffGam                                     # obm:474
    T72 <- CtriolSetPt * ScaEff                                                # obm:476
    T73 <- T71 * (C8 - T72)                                                    # obm:478
    T74 <- (exp(T73) - exp(-T73)) / (exp(T73) + exp(-T73))                     # obm:480 (tanh)
    T75 <- T70 * (0.85 * (1 + T74) + 0.15)                                     # obm:482
    T76 <- T70 * (0.85 * (1 - T74) + 0.15)                                     # obm:484

    # ---- Osteoblast apoptosis --------------------------------------------
    E0PicOBkb <- MultPicOBkb * Pic0                                            # obm:501
    EmaxPicOBkb <- FracPic0kb * Pic0                                           # obm:503
    EC50PicOBparenKb <- ((E0PicOBkb - EmaxPicOBkb) * TGFBact0^PicOBgamkb) / (E0PicOBkb - Pic0) - TGFBact0^PicOBgamkb  # obm:505
    EC50PicOBkb <- exp(log(EC50PicOBparenKb) / PicOBgamkb)                     # obm:507
    PicOBkb <- E0PicOBkb - (E0PicOBkb - EmaxPicOBkb) * TGFBact^PicOBgamkb / (TGFBact^PicOBgamkb + EC50PicOBkb^PicOBgamkb)  # obm:509
    PicOBkbEff <- PicOBkb / Pic0                                               # obm:512 (EST = 1)
    E0RUNX2kbEff <- E0RUNX2kbEffFACT * kb                                      # obm:515
    if (BCL2 > RUNX2thresh) {                                                  # obm:517
      RUNX2 <- BCL2 - RUNX2offset
    } else {
      RUNX2 <- RUNX20
    }
    RUNkbMax <- E0RUNX2kbEff * RUNkbMaxFact                                    # obm:519
    INparen <- (RUNkbMax * RUNX20^RUNkbGAM) / (E0RUNX2kbEff - kb) - RUNX20^RUNkbGAM  # obm:521
    RUNkb50 <- exp(log(INparen) / RUNkbGAM)                                    # obm:523
    RUNX2kbPrimeEff <- RUNkbMax * RUNX2^RUNkbGAM / (RUNX2^RUNkbGAM + RUNkb50^RUNkbGAM)  # obm:525
    kbprime <- E0RUNX2kbEff * PicOBkbEff - RUNX2kbPrimeEff                     # obm:527
    kbslow <- kbprime * Frackb                                                 # obm:529
    kbfast <- (kb * OB0 + kbslow * OBfast0 - kbslow * OB0) / OBfast0           # obm:531
    Frackb2 <- kbfast / kbprime                                                # obm:533

    # ---- Parathyroid hypertrophy, equation (4) ----------------------------
    INparenCtriol <- ((CtriolMax - CtriolMin) * CtriolConc0^CtriolPTgam) / (CtriolMax - 1) - CtriolConc0^CtriolPTgam  # obm:557
    Ctriol50 <- exp(log(INparenCtriol) / CtriolPTgam)                          # obm:559
    CtriolPTeff <- CtriolMax - (CtriolMax - CtriolMin) * C8^CtriolPTgam / (C8^CtriolPTgam + Ctriol50^CtriolPTgam)  # obm:561
    PTin <- PTout * CtriolPTeff                                                # obm:563

    # ---- PTH secretion, equations (2) and (6) -----------------------------
    FCTD <- (PTpool / PTpoolInit) * PThypertrophy                              # eq (2); obm:566
    # casrAgonist carries the constant millimolar calcium equivalents of
    # equation (6); it is zero unless dosed, in which case equation (6)
    # reduces to equation (2).
    CaSRdriver <- CaConc + casrAgonist                                         # eq (6)
    INparenCa <- (T58 - T61) * CaConc0^T59 / (T58 - SPTHbase) - CaConc0^T59    # obm:568
    T60 <- exp(log(INparenCa) / T59)                                           # obm:569
    T63 <- T58 - (T58 - T61) * CaSRdriver^T59 / (CaSRdriver^T59 + T60^T59)     # eq (6); obm:570
    SPTH <- T63 * FCTD                                                         # eq (2); obm:573

    # ---- Lumbar spine BMD, equation (5) -----------------------------------
    # BSAP / BSAP_baseline is represented by total osteoblasts relative to
    # baseline and CTx / CTx_baseline by active osteoclasts relative to
    # baseline; OpenBoneMin's [TABLE] block labels exactly these two ratios
    # "bone-specific alkaline phosphatase, %change from baseline" (obm:781)
    # and "serum CTx, %change from baseline" (obm:782). kin,BMD is not
    # printed in the paper; it is the baseline steady-state product
    # kout,BMD * BMD_LS(0) (obm:672).
    kinBMDls <- koutBMDls * BMDlsInit                                          # obm:672

    # =====================================================================
    # Differential equations
    # =====================================================================
    d/dt(PTH) <- SPTH - kout * PTH                                             # eq (2); obm:583
    d/dt(PTpool) <- (1 - PTpool) * T76 - PTpool * T75                          # eq (3); obm:586
    d/dt(PThypertrophy) <- PTin - PTout * PThypertrophy                        # eq (4); obm:589
    d/dt(calcitriol) <- alphaOHase - T69 * calcitriol                          # obm:591
    d/dt(alphaOHase) <- SE - T64 * alphaOHase                                  # obm:594
    d/dt(plasmaCa) <- J14 - J15 - J27 + J40                                    # obm:601
    d/dt(PhosGut) <- OralPhos * FracPhosAbs - J53                              # obm:613
    d/dt(IntraPO) <- J54 - J56                                                 # obm:614
    d/dt(ECCPhos) <- J41 - J42 - J48 + J53 - J54 + J56                         # obm:615
    d/dt(gutCa) <- OralCa * T85 - J40                                          # obm:621
    d/dt(gutCaAbsorp) <- T36 * (1 - gutCaAbsorp) - T37 * gutCaAbsorp           # obm:625
    d/dt(HAp) <- kHApIn * Osteoblast - kLShap * HAp                            # obm:628
    d/dt(OBfast) <- (bigDb / PicOB) * ROB * FracOBfast * Frackb2 - kbfast * OBfast  # obm:648
    d/dt(OBslow) <- (bigDb / PicOB) * ROB * (1 - FracOBfast) * Frackb - kbslow * OBslow  # obm:650
    d/dt(OC) <- kinOC2 - KLSoc * OC                                            # obm:653
    d/dt(ROB) <- ROBin - KPT * ROB                                             # obm:656
    d/dt(TGFB) <- kinTGF * (Osteoblast / OB0)^OBtgfGAM - koutTGFeqn            # obm:659
    d/dt(TGFBact) <- koutTGFeqn - koutTGFact * TGFBact                         # obm:664
    d/dt(BMDls) <- kinBMDls * (Osteoblast / OB0)^gamOB - koutBMDls * (OC / OCinit)^gamOCls * BMDls  # eq (5); obm:679
    d/dt(RANKL) <- kinL - koutL * RANKL - k1 * OPG * RANKL + k2 * OPG_RANKL - k3 * RANK * RANKL + k4 * RANK_RANKL  # obm:694
    d/dt(RANK) <- kinRNK * TGFBact^kinRNKgam - koutRNK * RANK - k3 * RANK * RANKL + k4 * RANK_RANKL  # obm:697
    d/dt(RANK_RANKL) <- k3 * RANK * RANKL - k4 * RANK_RANKL                    # obm:700
    d/dt(OPG_RANKL) <- k1 * OPG * RANKL - k2 * OPG_RANKL                       # obm:703
    d/dt(OPG) <- pO - k1 * OPG * RANKL + k2 * OPG_RANKL - koutOPG * OPG        # obm:706
    d/dt(boneCaExch) <- J15 - J14 + J14a - J15a                                # obm:713
    d/dt(boneCaNonExch) <- J15a - J14a                                         # obm:715
    d/dt(RX2) <- RX2Kin - RX2Kout * RX2                                        # obm:725
    d/dt(CREB) <- crebKin - crebKout * CREB                                    # obm:727
    d/dt(BCL2) <- bcl2Kout * CREB * RX2 - bcl2Kout * BCL2                      # obm:729
    d/dt(urineCa) <- J27                                                       # obm:774
    d/dt(casrAgonist) <- 0                                                     # eq (6); a bolus sets and holds the constant Ca equivalent

    # =====================================================================
    # Initial conditions: the upstream steady state (obm $INIT, obm:222-263,
    # with the four values that obm [MAIN] overrides computed above).
    # =====================================================================
    PTH(0)           <- PTHinit
    PTpool(0)        <- PTpoolInit
    PThypertrophy(0) <- PThyperInit
    calcitriol(0)    <- CtriolInit
    alphaOHase(0)    <- AOHinit
    plasmaCa(0)      <- CaInit
    ECCPhos(0)       <- ECCPhosInit
    PhosGut(0)       <- PhosGutInit
    IntraPO(0)       <- IntraPOInit
    gutCa(0)         <- GutCaInit
    gutCaAbsorp(0)   <- GutCaAbsInit
    HAp(0)           <- HApInit
    boneCaExch(0)    <- BoneCaExchInit
    boneCaNonExch(0) <- BoneCaNonExchInit
    OBfast(0)        <- OBfast0
    OBslow(0)        <- OBslow0
    OC(0)            <- OCinit
    ROB(0)           <- ROBinit
    TGFB(0)          <- TGFB0
    TGFBact(0)       <- TGFBact0
    RANKL(0)         <- RANKLinit
    RANK(0)          <- RANKinit
    RANK_RANKL(0)    <- M0
    OPG(0)           <- OPGinit
    OPG_RANKL(0)     <- N0
    RX2(0)           <- RX2init
    CREB(0)          <- CREBinit
    BCL2(0)          <- BCL2init
    BMDls(0)         <- BMDlsInit
    urineCa(0)       <- 0
    casrAgonist(0)   <- 0

    # =====================================================================
    # Reported outputs. The paper plots PTH, calcium, BSAP, sCTx and BMD as
    # a percent of their baseline values (Figures 2-5); obm [TABLE]
    # (obm:778-787) defines the same quantities.
    # =====================================================================
    PTHpm <- PTHconc                    # PTH concentration (pmol/L); obm:778
    CaC <- CaConc                       # total serum calcium (mmol/L); obm:783
    PhosC <- C2                         # serum phosphate (mmol/L)
    CtriolC <- C8                       # serum calcitriol (pmol/L)
    PTHpct <- 100 * PTH / PTHinit       # PTH, percent of baseline
    CApct <- 100 * plasmaCa / CaInit    # calcium, percent of baseline; obm:780
    BSAPpct <- 100 * Osteoblast / OB0   # BSAP surrogate, percent of baseline; obm:781
    CTXpct <- 100 * OC / OCinit         # sCTx surrogate, percent of baseline; obm:782
    BMDlspct <- (BMDls - 1) * 100       # lumbar spine BMD, percent change from baseline; obm:784
  })
}
