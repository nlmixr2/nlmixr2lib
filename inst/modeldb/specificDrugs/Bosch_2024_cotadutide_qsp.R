Bosch_2024_cotadutide_qsp <- function() {
  description <- paste(
    "QSP. 4GI quantitative systems pharmacology model (glucose, insulin,",
    "GLP-1, glucagon, GIP) coupled to a one-compartment first-order-",
    "absorption cotadutide PK model in adults with type 2 diabetes",
    "mellitus (Bosch 2024). Cotadutide is a dual GLP-1/glucagon receptor",
    "agonist; in vivo EC50s for cotadutide on each receptor are derived",
    "from the in vitro EC50 ratio vs the endogenous ligand (Eq 1). The",
    "drug's free-fraction-corrected central concentration drives four",
    "saturable Emax effects on the system: (1) stimulation of glucose-",
    "dependent insulin secretion via GLP-1R, (2) inhibition of meal-",
    "glucose absorption via GLP-1R, (3) inhibition of glucagon",
    "production via GLP-1R, and (4) stimulation of glucose production",
    "via GCGR. A fifth Emax inhibits endogenous active GLP-1 production",
    "(Eq 3). The placebo arm's lifestyle-change effect on fasting plasma",
    "glucose is modelled as an inverse Bateman attenuation of endogenous",
    "glucose production (Eq 2). Cotadutide PK structure and typical",
    "values are fixed from the upstream popPK analysis of Guan et al.",
    "2022 (KA=0.343 1/h, CL=1.04 L/h, V=18.7 L). All 4GI system-",
    "specific disposition and effect parameters are fixed from the",
    "upstream 4GI model of Bosch et al. 2022; meal-effect, baseline,",
    "lifestyle and EMAX_5/EC50_5S parameters were re-estimated against",
    "the cotadutide MAD/Ph2a dataset (NCT02548585; n=51, T2DM). Five",
    "outputs: plasma glucose (mmol/L), insulin (pmol/L), GLP-1 (pmol/L),",
    "glucagon (pmol/L) and GIP (pmol/L), each with proportional residual",
    "error. Individual fasting plasma glucose enters via the FPG",
    "covariate; meal glucose enters as dosing events on the glucose-gut",
    "compartment. Defaults are T2DM; healthy-volunteer parameter set",
    "from Bosch 2022 is given in source-trace comments. No IIV is",
    "encoded (sequential model fit with individual PK / glucose-",
    "baseline inputs from Guan 2022 and the observed dataset).",
    sep = " "
  )
  reference <- paste(
    "Bosch R, Petrone M, Arends R, Vicini P, Sijbrands EJG, Hoefman S,",
    "Snelder N (2024).",
    "Characterisation of cotadutide's dual GLP-1/glucagon receptor",
    "agonistic effects on glycaemic control using an in vivo human",
    "glucose regulation quantitative systems pharmacology model.",
    "British Journal of Pharmacology 181(12):1874-1885.",
    "doi:10.1111/bph.16336. PMID: 38403793.",
    "Cotadutide PK structure fixed from Guan H et al. (2022),",
    "population pharmacokinetics of cotadutide (cited in Bosch 2024",
    "Methods Section 2.3); 4GI system parameters fixed from Bosch R",
    "et al. (2022), the original 4GI model (cited in Bosch 2024",
    "Methods Section 2.3).",
    sep = " "
  )
  vignette <- "Bosch_2024_cotadutide_qsp"

  paper_specific_compartments <- c(
    "glucose", "glucose_per", "glucose_gut", "glucose_buffer",
    "glucose_tr1", "glucose_tr2", "glucose_tr3",
    "insulin", "insulin_eff",
    "glp1", "glucagon", "gip", "gip_per"
  )

  units <- list(
    time          = "h",
    dosing        = paste(
      "cotadutide dose into the depot compartment must be in pmol",
      "(convert from ug via amt_pmol = amt_ug * 1e6 / MW_cotadutide_g_per_mol;",
      "cotadutide MW is approximately 4549 g/mol so 100 ug ~= 22000 pmol).",
      "Meal glucose dose into the glucose_gut compartment must be in mmol",
      "(convert from grams via amt_mmol = amt_g * 1000 / 180.16).",
      sep = " "
    ),
    concentration = paste(
      "Cotadutide total plasma concentration Cmedi is in pmol/L",
      "(central / vc); free concentration Cmedif = Cmedi * fumedi.",
      "Endogenous outputs: plasma glucose Cglc in mmol/L, insulin Cins",
      "in pmol/L, GLP-1 Cglp in pmol/L, glucagon Cglg in pmol/L,",
      "GIP Cgip in pmol/L.",
      sep = " "
    )
  )

  covariateData <- list(
    FPG = list(
      description        = "Individual fasting plasma glucose at study baseline (BSLglc in the paper notation). Used as the per-subject initial condition for the central glucose state and as the reference value for the power-function feedback of glucose on glucagon production. Time-fixed per subject.",
      units              = "mmol/L (paper reports baseline glucose in mmol/L for the cotadutide MAD/Ph2a cohort; multiply mg/dL by 0.0555 to convert)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required per-subject covariate. Set to the observed individual baseline value; for typical-subject simulations use the cohort median (approximately 7.67 mmol/L per supplement THETA(64) median MAD/Ph2a baseline). Bosch 2024 Methods Section 2.5: 'Individual glucose baseline values were used as input'. The paper estimates baseline insulin / GLP-1 / glucagon / GIP as typical values (Table 3) so only the glucose baseline needs a per-subject covariate.",
      source_name        = "BSLglc (IBGLC in the supplement mrgSolve code)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 45L,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = paste(
      "Overweight or obese adults with type 2 diabetes mellitus",
      "(NCT02548585; the MAD/Ph2a cotadutide study reported in Ambery",
      "et al. 2018). N = 51 randomised (25 cotadutide, 26 placebo);",
      "6 subjects (2 placebo, 4 cotadutide) excluded as outliers per",
      "Bosch 2024 Methods Section 2.4 (CWRES > 6 or extreme baseline",
      "insulin / GLP-1), leaving n = 45 in the population PD",
      "analysis. MAD part: 5 cohorts (A-E) with up-titrated SC daily",
      "doses to top doses of 100, 150, 200, 300, 300 ug. Phase 2a",
      "part: 25 subjects to a 200 ug top dose with 4-day 100 ug + 4-",
      "day 150 ug up-titration followed by 33 days at 200 ug; 26",
      "placebo. Detailed cohort table is Bosch 2024 Table 1.",
      sep = " "
    ),
    dose_range     = paste(
      "Cotadutide 100-300 ug subcutaneous once daily before",
      "breakfast, with up-titration over 7-41 days depending on",
      "cohort (Bosch 2024 Table 1). MMTT (mixed meal tolerance test):",
      "standardised Ensure Plus(R) breakfast on selected study days.",
      sep = " "
    ),
    regions        = NA_character_,
    notes          = paste(
      "Cohort sizes from Bosch 2024 Table 1; the paper does not",
      "tabulate detailed demographics (age / weight / sex / race",
      "distribution) for the Bosch 2024 analysis cohort and refers",
      "the reader to Ambery et al. 2018 for the parent study design.",
      "PK samples taken pre-dose and up to 72 h post-dose; on MMTT",
      "days, plasma glucose, insulin, GLP-1, glucagon and GIP were",
      "measured before and up to 240 min after the MMTT meal.",
      "Individual cotadutide PK time-courses came from the empirical",
      "Bayes estimates (CL, VC, KA) of the upstream Guan 2022 popPK",
      "model; in this packaged model file the typical values of",
      "those EBEs (KA = 0.343 1/h, CL = 1.04 L/h, V = 18.7 L) are",
      "fixed structural defaults and IIV is not encoded.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Cotadutide PK -- FIXED at the Guan 2022 popPK typical values per
    # the supplement mrgSolve $PK block defaults (`KA = 0.343`, `CL = 1.04`,
    # `V = 18.7` when individual EBE inputs are not supplied). Bosch 2024
    # Section 2.3 fixes the cotadutide PK structure and individual EBEs
    # from Guan et al. 2022; this packaged model uses the typical values.
    # =====================================================================
    lka <- fixed(log(0.343))       ; label("Cotadutide absorption rate ka (1/h; FIXED from Guan 2022 typical value)")  # Bosch 2024 supplement $PK block default `if(IKA <= 0) KA = 0.343`
    lcl <- fixed(log(1.04))        ; label("Cotadutide clearance CL (L/h; FIXED from Guan 2022 typical value)")        # Bosch 2024 supplement $PK block default `if(ICL <= 0) CL = 1.04`
    lvc <- fixed(log(18.7))        ; label("Cotadutide central volume Vc (L; FIXED from Guan 2022 typical value)")     # Bosch 2024 supplement $PK block default `if(IV <= 0) V = 18.7`

    fumedi <- fixed(0.0023)        ; label("Cotadutide plasma free fraction (unitless; FIXED at 0.23 percent)")        # Bosch 2024 Section 2.2 and supplement `fumedi = 0.0023`

    # =====================================================================
    # 4GI glucose disposition (T2DM) -- all FIXED from Bosch 2022 (Table S1).
    # Healthy-volunteer alternatives: CLglc = 5.36 L/h (supplement THETA46),
    # CLglci = 0.072 (L/h)/(pmol/L) (THETA47); not used here because the
    # MAD/Ph2a cohort is entirely T2DM. Glucose amounts are in mmol;
    # concentrations in mmol/L.
    # =====================================================================
    lclglc   <- fixed(log(1.72))   ; label("Insulin-independent glucose clearance CLglc T2DM (L/h; FIXED from Bosch 2022 Table S1)")           # Bosch 2024 supplement Table S1; THETA(1)
    lclgi    <- fixed(log(0.0256)) ; label("Insulin-dependent glucose clearance CLglci T2DM ((L/h)/(pmol/L); FIXED from Bosch 2022 Table S1)") # Bosch 2024 supplement Table S1; THETA(2)
    lqglc    <- fixed(log(26.5))   ; label("Glucose inter-compartmental clearance Qglc (L/h; FIXED from Bosch 2022 Table S1)")                  # Bosch 2024 supplement Table S1; THETA(3)
    lvcglc   <- fixed(log(9.33))   ; label("Glucose central volume VCglc (L; FIXED from Bosch 2022 Table S1)")                                  # Bosch 2024 supplement Table S1; THETA(4)
    lvpglc   <- fixed(log(8.56))   ; label("Glucose peripheral volume VPglc (L; FIXED from Bosch 2022 Table S1)")                               # Bosch 2024 supplement Table S1; THETA(5)
    lkeglc   <- fixed(log(0.281))  ; label("Glucose buffer-to-transit rate constant Keglc (1/h; FIXED from Bosch 2022 Table S1)")               # Bosch 2024 supplement Table S1; THETA(8)
    lkelglc  <- fixed(log(1.93))   ; label("Glucose transit chain rate constant Kelglc (1/h; FIXED from Bosch 2022 Table S1)")                  # Bosch 2024 supplement Table S1; THETA(9)

    # MMTT-meal glucose absorption -- estimated in Bosch 2024 Table 3.
    lkaglc    <- log(3.58)                  ; label("Glucose absorption rate constant for MMTT meal Kaglc (1/h)")        # Bosch 2024 Table 3 estimate (RSE 7.58%)
    logitfglc <- log(0.334 / (1 - 0.334))   ; label("Logit of glucose bioavailability Fglc for MMTT meal (unitless)")    # Bosch 2024 Table 3 estimate 0.334 (RSE 8.89%)

    # =====================================================================
    # 4GI insulin disposition -- all FIXED from Bosch 2022 (Table S1).
    # Insulin amounts are in pmol; concentrations in pmol/L.
    # =====================================================================
    lclins   <- fixed(log(73.2))      ; label("Insulin clearance CLins (L/h; FIXED from Bosch 2022 Table S1)")                                  # Bosch 2024 supplement Table S1; THETA(10)
    lvcins   <- fixed(log(6.09))      ; label("Insulin central volume VCins (L; FIXED from Bosch 2022 Table S1)")                               # Bosch 2024 supplement Table S1; THETA(11)
    lke0ins  <- fixed(-0.158729)      ; label("Log effect-compartment insulin equilibration rate KE0ins (log-1/h; FIXED, exp = 0.853)")          # Bosch 2024 supplement THETA(12) (log-scale: exp(-0.159) = 0.853 1/h)

    # =====================================================================
    # 4GI GLP-1 disposition -- all FIXED from Bosch 2022 (Table S1).
    # GLP-1 amounts are in pmol; concentrations in pmol/L.
    # Michaelis-Menten elimination: the source $THETA stores VM and KM on
    # the log scale (`VM = exp(THETA(14))`, `KM = exp(THETA(15))`).
    # =====================================================================
    lvcglp     <- fixed(log(16.0))    ; label("GLP-1 central volume VCglp (L; FIXED from Bosch 2022 Table S1)")                                 # Bosch 2024 supplement Table S1; THETA(13)
    lvmax_glp1 <- fixed(7.96952)      ; label("Log GLP-1 maximum MM elimination rate VM (log-pmol/(L*h); FIXED, exp = 2893)")                   # Bosch 2024 supplement THETA(14); VM = exp(7.97) = 2893 pmol/(L*h)
    lkm_glp1   <- fixed(4.90603)      ; label("Log GLP-1 MM half-saturation KM (log-pmol/L; FIXED, exp = 135)")                                  # Bosch 2024 supplement THETA(15); KM = exp(4.91) = 135 pmol/L

    # =====================================================================
    # 4GI glucagon disposition -- all FIXED from Bosch 2022 (Table S1).
    # Glucagon amounts are in pmol; concentrations in pmol/L.
    # =====================================================================
    lclglg <- fixed(log(453))         ; label("Glucagon clearance CLglg (L/h; FIXED from Bosch 2022 Table S1)")                                  # Bosch 2024 supplement Table S1; THETA(16)
    lvcglg <- fixed(log(64.6))        ; label("Glucagon central volume VCglg (L; FIXED from Bosch 2022 Table S1)")                               # Bosch 2024 supplement Table S1; THETA(17)

    # =====================================================================
    # 4GI GIP disposition -- all FIXED from Bosch 2022 (Table S1).
    # GIP amounts are in pmol; concentrations in pmol/L.
    # =====================================================================
    lclgip <- fixed(log(86.8))        ; label("GIP clearance CLgip (L/h; FIXED from Bosch 2022 Table S1)")                                       # Bosch 2024 supplement Table S1; THETA(18)
    lvcgip <- fixed(log(9.21))        ; label("GIP central volume VCgip (L; FIXED from Bosch 2022 Table S1)")                                    # Bosch 2024 supplement Table S1; THETA(19)
    lqgip  <- fixed(log(49.4))        ; label("GIP inter-compartmental clearance Qgip (L/h; FIXED from Bosch 2022 Table S1)")                    # Bosch 2024 supplement Table S1; THETA(20)
    lvpgip <- fixed(log(22.8))        ; label("GIP peripheral volume VPgip (L; FIXED from Bosch 2022 Table S1)")                                 # Bosch 2024 supplement Table S1; THETA(21)

    # =====================================================================
    # MMTT-specific meal-effect parameters -- estimated in Bosch 2024
    # (Table 3). These multiply the buffer-compartment glucose amount
    # (or the end-of-transit amount for FDGLP_2) to drive food-induced
    # increases in GLP-1 / GIP / glucagon secretion. Units are 1/mmol so
    # that (parameter * amount_in_mmol) is unitless.
    # =====================================================================
    lfdglp   <- log(0.0150)           ; label("Food effect on GLP-1 via glucose buffer FDGLP (1/mmol)")                                       # Bosch 2024 Table 3 estimate 0.0150 (RSE 23.2%)
    lfdglp_2 <- log(0.113)            ; label("Food effect on GLP-1 via end-of-transit glucose FDGLP_2 (1/mmol)")                              # Bosch 2024 Table 3 estimate 0.113 (RSE 25.7%)
    lfdgip   <- log(0.107)            ; label("Food effect on GIP via glucose buffer FDGIP (1/mmol)")                                          # Bosch 2024 Table 3 estimate 0.107 (RSE 39.2%)
    lfdglg   <- log(0.0201)           ; label("Food effect on glucagon via glucose buffer FDGLG (1/mmol)")                                     # Bosch 2024 Table 3 estimate 0.0201 (RSE 64.4%)

    # =====================================================================
    # Glucose-on-other-species feedback exponents (Bosch 2022 Table S2).
    # GLCINS_S is the glucose-on-insulin secretion power exponent (used
    # in `1 + STglc * Cglc^GLCINS_S`). GLCGLG_POWH is the high-glucose
    # branch of the glucose-on-glucagon production feedback; the T2DM
    # low-glucose branch GLCGLG_POWL is fixed at 0 in the source code
    # logic (no hypoglycaemic glucagon counter-regulation in T2DM).
    # =====================================================================
    lglcins_s    <- fixed(log(2.46))   ; label("Glucose-on-insulin secretion power coefficient GLCINS_S (1/mM; FIXED)")             # Bosch 2024 supplement THETA(26); Table S2 GLCINS_S = 2.46
    lglcglg_powh <- fixed(log(0.925))  ; label("Glucose-on-glucagon production power exponent high-glc branch GLCGLG_POWH (unitless; FIXED)")  # Bosch 2024 supplement THETA(27); Table S2 GLCGLG_POWH = 0.925

    # =====================================================================
    # GLP-1 receptor effect parameters on (1) insulin secretion, (2)
    # glucose absorption (gastric emptying), (3) glucagon production --
    # all FIXED from Bosch 2022 (Table S2). The source $THETA stores
    # EC50_1, EC50_2, EC50_3 on the log scale; back-transformed values
    # match Table S2 (26.6, 144, 99.5 pmol/L). The in vivo cotadutide
    # EC50 for each effect is derived from these by the in-vitro EC50
    # ratio (Bosch 2024 Eq 1), computed in model() as ECGLP1 / 2 / 3.
    # =====================================================================
    lemax_1 <- fixed(log(10.7))       ; label("GLP-1 Emax on insulin secretion EMAX_1 (unitless; FIXED)")                                       # Bosch 2024 supplement THETA(28); Table S2 EMAX_1 = 10.7
    lec50_1 <- fixed(log(26.6))       ; label("GLP-1 EC50 on insulin secretion EC50_1 (pmol/L; FIXED)")                                          # Bosch 2024 supplement THETA(29); back-transformed exp(3.29) = 26.6; Table S2 26.6
    lhill_1 <- fixed(log(1.79))       ; label("GLP-1 Hill on insulin secretion HILL_1 (unitless; FIXED)")                                        # Bosch 2024 supplement THETA(30); Table S2 HILL_1 = 1.79

    lemax_2 <- fixed(log(1))          ; label("GLP-1 Emax on glucose absorption EMAX_2 (unitless; FIXED at 1)")                                  # Bosch 2024 supplement THETA(31); Table S2 EMAX_2 = 1
    lec50_2 <- fixed(log(144))        ; label("GLP-1 EC50 on glucose absorption EC50_2 (pmol/L; FIXED)")                                         # Bosch 2024 supplement THETA(32); back-transformed exp(4.97) = 144; Table S2 144
    lhill_2 <- fixed(log(1))          ; label("GLP-1 Hill on glucose absorption HILL_2 (unitless; FIXED at 1)")                                  # Bosch 2024 supplement THETA(33); Table S2 HILL_2 = 1

    lemax_3 <- fixed(log(1))          ; label("GLP-1 Emax on glucagon production EMAX_3 (unitless; FIXED at 1)")                                 # Bosch 2024 supplement THETA(34); Table S2 EMAX_3 = 1
    lec50_3 <- fixed(log(99.5))       ; label("GLP-1 EC50 on glucagon production EC50_3 (pmol/L; FIXED)")                                        # Bosch 2024 supplement THETA(35); back-transformed exp(4.60) = 99.5; Table S2 99.5
    lhill_3 <- fixed(log(1))          ; label("GLP-1 Hill on glucagon production HILL_3 (unitless; FIXED at 1)")                                 # Bosch 2024 supplement THETA(36); Table S2 HILL_3 = 1

    # Cotadutide effect on endogenous GLP-1 (Bosch 2024 Eq 3 and Table 3).
    # EC50_5S is a scaling factor on the cotadutide in vitro EC50 for GLP-1
    # (free): EC50_5 = ECmGLP * EC50_5S (~ 0.076 * 10.9 = 0.83 pmol/L free
    # cotadutide). HILL_5 fixed at 5 per Section 3.1.
    lemax_5  <- log(0.321)            ; label("Cotadutide Emax inhibition of endogenous GLP-1 EMAX_5 (unitless)")                                # Bosch 2024 Table 3 estimate 0.321 (RSE 25.7%)
    lec50_5s <- log(10.9)             ; label("Cotadutide in-vitro-to-in-vivo EC50 scaling factor EC50_5S (unitless)")                            # Bosch 2024 Table 3 estimate 10.9 (RSE 37.9%)
    lhill_5  <- fixed(log(5))         ; label("Cotadutide Hill on endogenous GLP-1 inhibition HILL_5 (unitless; FIXED at 5)")                     # Bosch 2024 Table 3 footnote and Section 3.1

    # =====================================================================
    # Glucagon receptor effect on glucose production (Bosch 2022 Table S2).
    # EC50_4 stored log-scale in source $THETA. In vivo cotadutide GCGR
    # EC50 is derived from the in vitro ratio (Eq 1), computed in model()
    # as ECGLG1.
    # =====================================================================
    lemax_4 <- fixed(log(6.73))       ; label("Glucagon Emax on glucose production EMAX_4 (unitless; FIXED)")                                    # Bosch 2024 supplement THETA(40); Table S2 EMAX_4 = 6.73
    lec50_4 <- fixed(log(98.5))       ; label("Glucagon EC50 on glucose production EC50_4 (pmol/L; FIXED)")                                       # Bosch 2024 supplement THETA(41); back-transformed exp(4.59) = 98.5; Table S2 98.5
    lhill_4 <- fixed(log(1))          ; label("Glucagon Hill on glucose production HILL_4 (unitless; FIXED at 1)")                                # Bosch 2024 supplement THETA(42); Table S2 HILL_4 = 1

    # =====================================================================
    # GIP-on-glucagon production power exponent (Bosch 2022 Table S2). The
    # T2DM GIP-on-insulin scaler GIPINS and exponent POW_3 are both fixed
    # at 0 in the source code (HV alternatives: GIPINS = 1, POW_3 = 0.286
    # per Table S2); they appear hard-coded as 0 inside model().
    # =====================================================================
    lpow_4  <- fixed(log(0.109))      ; label("GIP-on-glucagon production power POW_4 (unitless; FIXED)")                                          # Bosch 2024 supplement THETA(45); Table S2 POW_4 = 0.109

    # In vitro EC50s (free cotadutide and endogenous ligands) -- Table 2.
    # These FIXED inputs feed the in-vivo cotadutide EC50 derivation (Eq 1).
    lecmglp <- fixed(log(0.076))      ; label("Cotadutide in vitro EC50 on GLP-1R (pmol/L; FIXED at Table 2 value)")                              # Bosch 2024 Table 2 in vitro GLP-1R EC50 cotadutide = 0.076 pmol/L (free, unpublished)
    lecmglg <- fixed(log(0.088))      ; label("Cotadutide in vitro EC50 on GCGR (pmol/L; FIXED at Table 2 value)")                                # Bosch 2024 Table 2 in vitro GCGR EC50 cotadutide = 0.088 pmol/L (free, unpublished)
    lecglp  <- fixed(log(1.92))       ; label("Endogenous GLP-1 in vitro EC50 on GLP-1R (pmol/L; FIXED at Table 2 value)")                        # Bosch 2024 Table 2 in vitro GLP-1R EC50 endogenous GLP-1 = 1.92 pmol/L (unpublished)
    lecglg  <- fixed(log(1.54))       ; label("Endogenous glucagon in vitro EC50 on GCGR (pmol/L; FIXED at Table 2 value)")                       # Bosch 2024 Table 2 in vitro GCGR EC50 endogenous glucagon = 1.54 pmol/L (unpublished)

    # =====================================================================
    # Lifestyle-change effect on endogenous glucose production (Bosch 2024
    # Eq 2 -- inverse Bateman attenuation). LSCI = 0.566 maximum fractional
    # reduction in KINglc; Klsc = 10 1/day FIXED (fast onset); Kred is the
    # decay rate of the lifestyle effect. NOTE: rate units are 1/day; the
    # source paper text gives Kred in 1/day (the supplement Table 3 column
    # header `(h-1)` is a typo confirmed by the Section 3.1 narrative
    # 'Kred was 0.0164 d^-1, consistent with a half-life ... of 45 days').
    # =====================================================================
    llsci <- log(0.566)               ; label("Lifestyle-change amplitude LSCI (fraction)")                                                       # Bosch 2024 Table 3 estimate 0.566 (RSE 13.1%)
    lklsc <- fixed(log(10))           ; label("Lifestyle-change onset rate Klsc (1/day; FIXED)")                                                    # Bosch 2024 Table 3 fixed at 10 per Section 3.1
    lkred <- log(0.0164)              ; label("Lifestyle-change reduction rate Kred (1/day)")                                                       # Bosch 2024 Table 3 estimate 0.0164 (RSE 55.0%)

    # =====================================================================
    # Baseline typical values for the endogenous species other than glucose
    # (which uses the per-subject FPG covariate). Bosch 2024 Table 3.
    # =====================================================================
    lbslins <- log(138)               ; label("Typical baseline insulin BSLins (pmol/L)")                                                          # Bosch 2024 Table 3 estimate 138 (RSE 8.44%)
    lbslglp <- log(18.0)              ; label("Typical baseline GLP-1 BSLglp (pmol/L)")                                                            # Bosch 2024 Table 3 estimate 18.0 (RSE 6.53%)
    lbslglg <- log(49.1)              ; label("Typical baseline glucagon BSLglg (pmol/L)")                                                         # Bosch 2024 Table 3 estimate 49.1 (RSE 11.5%)
    lbslgip <- log(19.8)              ; label("Typical baseline GIP BSLgip (pmol/L)")                                                              # Bosch 2024 Table 3 estimate 19.8 (RSE 18.7%)

    # =====================================================================
    # Residual error -- Bosch 2024 Table 3 reports sigma^2 for proportional
    # error Y = IPRED * (1 + EPS) with Var(EPS) = sigma^2; the nlmixr2
    # `prop()` parameter is sqrt(sigma^2).
    # =====================================================================
    propSd_Cglc <- 0.198              ; label("Proportional residual SD on plasma glucose (fraction)")                                            # Bosch 2024 Table 3 sigma^2_glc = 0.0393; SD = sqrt(0.0393) = 0.198 (RSE 8.72%)
    propSd_Cins <- 0.804              ; label("Proportional residual SD on plasma insulin (fraction)")                                            # Bosch 2024 Table 3 sigma^2_ins = 0.647; SD = sqrt(0.647) = 0.804 (RSE 14.1%)
    propSd_Cglp <- 0.599              ; label("Proportional residual SD on plasma GLP-1 (fraction)")                                              # Bosch 2024 Table 3 sigma^2_glp = 0.358; SD = sqrt(0.358) = 0.599 (RSE 27.5%)
    propSd_Cglg <- 0.633              ; label("Proportional residual SD on plasma glucagon (fraction)")                                           # Bosch 2024 Table 3 sigma^2_glg = 0.401; SD = sqrt(0.401) = 0.633 (RSE 27.7%)
    propSd_Cgip <- 0.811              ; label("Proportional residual SD on plasma GIP (fraction)")                                                # Bosch 2024 Table 3 sigma^2_gip = 0.658; SD = sqrt(0.658) = 0.811 (RSE 10.5%)
  })

  model({
    # ---------------------------------------------------------------------
    # 1. Individual structural parameters (no IIV encoded -- sequential
    #    model fit with individual PK and glucose-baseline inputs).
    # ---------------------------------------------------------------------
    ka  <- exp(lka)
    cl  <- exp(lcl)
    vc  <- exp(lvc)
    kel <- cl / vc

    clglc  <- exp(lclglc)
    clgi   <- exp(lclgi)
    qglc   <- exp(lqglc)
    vcglc  <- exp(lvcglc)
    vpglc  <- exp(lvpglc)
    keglc  <- exp(lkeglc)
    kelglc <- exp(lkelglc)
    kaglc  <- exp(lkaglc)
    fglc   <- 1 / (1 + exp(-logitfglc))

    clins  <- exp(lclins)
    vcins  <- exp(lvcins)
    ke0ins <- exp(lke0ins)

    vcglp     <- exp(lvcglp)
    vmax_glp1 <- exp(lvmax_glp1)
    km_glp1   <- exp(lkm_glp1)

    clglg <- exp(lclglg)
    vcglg <- exp(lvcglg)

    clgip <- exp(lclgip)
    vcgip <- exp(lvcgip)
    qgip  <- exp(lqgip)
    vpgip <- exp(lvpgip)

    fdglp   <- exp(lfdglp)
    fdglp_2 <- exp(lfdglp_2)
    fdgip   <- exp(lfdgip)
    fdglg   <- exp(lfdglg)
    fdins   <- 0  # FIXED at 0 per supplement code -- food effect on insulin is wholly mediated through GLP-1 / GIP

    glcins_s    <- exp(lglcins_s)
    glcglg_powh <- exp(lglcglg_powh)

    emax_1 <- exp(lemax_1) ; ec50_1 <- exp(lec50_1) ; hill_1 <- exp(lhill_1)
    emax_2 <- exp(lemax_2) ; ec50_2 <- exp(lec50_2) ; hill_2 <- exp(lhill_2)
    emax_3 <- exp(lemax_3) ; ec50_3 <- exp(lec50_3) ; hill_3 <- exp(lhill_3)
    emax_4 <- exp(lemax_4) ; ec50_4 <- exp(lec50_4) ; hill_4 <- exp(lhill_4)
    emax_5  <- exp(lemax_5)
    ec50_5s <- exp(lec50_5s)
    hill_5  <- exp(lhill_5)

    pow_4 <- exp(lpow_4)

    ecmglp <- exp(lecmglp)
    ecmglg <- exp(lecmglg)
    ecglp  <- exp(lecglp)
    ecglg  <- exp(lecglg)
    fu     <- fumedi

    lsci <- exp(llsci)
    klsc <- exp(lklsc)
    kred <- exp(lkred)

    bslglc <- FPG
    bslins <- exp(lbslins)
    bslglp <- exp(lbslglp)
    bslglg <- exp(lbslglg)
    bslgip <- exp(lbslgip)

    # ---------------------------------------------------------------------
    # 2. In-vivo cotadutide EC50 derivation (Bosch 2024 Eq 1):
    #      EC50_cota_in_vivo = (EC50_cota_in_vitro / EC50_endogenous_in_vitro)
    #                          * EC50_endogenous_in_vivo
    # ---------------------------------------------------------------------
    ecglp1 <- ecmglp / ecglp * ec50_1
    ecglp2 <- ecmglp / ecglp * ec50_2
    ecglp3 <- ecmglp / ecglp * ec50_3
    ecglg1 <- ecmglg / ecglg * ec50_4
    ec50_5 <- ecmglp * ec50_5s

    # ---------------------------------------------------------------------
    # 3. Cotadutide plasma concentration (free) and endogenous species
    #    concentrations.
    # ---------------------------------------------------------------------
    Cmedi  <- central / vc
    Cmedif <- Cmedi * fu

    Cglc <- glucose  / vcglc
    Cins <- insulin  / vcins
    Cglp <- glp1     / vcglp
    Cglg <- glucagon / vcglg
    Cgip <- gip      / vcgip

    # ---------------------------------------------------------------------
    # 4. Glucose-on-glucagon feedback (T2DM): high-glc branch when
    #    Cglc >= bslglc, low-glc branch (exponent = 0, i.e. no feedback)
    #    when Cglc < bslglc. The T2DM cohort loses hypoglycaemic glucagon
    #    counter-regulation (Bosch 2024 supplement `if(PAT == 0 & Cglc > 0
    #    & Cglc < BSLglc) POW_2 = 0`).
    # ---------------------------------------------------------------------
    glcglg_pow_eff <- glcglg_powh * (Cglc >= bslglc)
    glcEFFglg <- (bslglc / Cglc)^glcglg_pow_eff

    # ---------------------------------------------------------------------
    # 5. GLP-1 receptor effects -- competitive-additive Emax over endogenous
    #    GLP-1 and free cotadutide at their respective in-vivo EC50s
    #    (Bosch 2024 supplement NUM/DENOM construct).
    # ---------------------------------------------------------------------
    glpins_s0 <- emax_1 * (bslglp / ec50_1)^hill_1 / (1 + (bslglp / ec50_1)^hill_1)
    num1 <- (Cglp / ec50_1)^hill_1 + (Cmedif / ecglp1)^hill_1
    den1 <- 1 + (Cglp / ec50_1)^hill_1 + (Cmedif / ecglp1)^hill_1
    glpins_s <- emax_1 * num1 / den1

    num2 <- (Cglp / ec50_2)^hill_2 + (Cmedif / ecglp2)^hill_2
    den2 <- 1 + (Cglp / ec50_2)^hill_2 + (Cmedif / ecglp2)^hill_2
    glpglu_ai <- emax_2 * num2 / den2

    glpglg_i0 <- emax_3 * (bslglp / ec50_3)^hill_3 / (1 + (bslglp / ec50_3)^hill_3)
    num3 <- (Cglp / ec50_3)^hill_3 + (Cmedif / ecglp3)^hill_3
    den3 <- 1 + (Cglp / ec50_3)^hill_3 + (Cmedif / ecglp3)^hill_3
    glpglg_i <- emax_3 * num3 / den3
    glpEFFglg <- (1 - glpglg_i) / (1 - glpglg_i0)

    # Cotadutide self-inhibition of endogenous GLP-1 production (Eq 3).
    glpglp_i <- emax_5 * (Cmedif / ec50_5)^hill_5 / (1 + (Cmedif / ec50_5)^hill_5)
    glpEFFglp <- 1 - glpglp_i

    # ---------------------------------------------------------------------
    # 6. Glucagon-on-glucose-production effect (Emax; cotadutide free
    #    concentration adds competitively at the GCGR).
    # ---------------------------------------------------------------------
    glgglc_s0 <- emax_4 * (bslglg / ec50_4)^hill_4 / (1 + (bslglg / ec50_4)^hill_4)
    num4 <- (Cglg / ec50_4)^hill_4 + (Cmedif / ecglg1)^hill_4
    den4 <- 1 + (Cglg / ec50_4)^hill_4 + (Cmedif / ecglg1)^hill_4
    glgglc_s <- emax_4 * num4 / den4
    glgEFFglc <- (1 + glgglc_s) / (1 + glgglc_s0)

    # ---------------------------------------------------------------------
    # 7. GIP feedback effects (T2DM: GIPINS = 0, POW_3 = 0; GIP does not
    #    drive insulin secretion). GIP-on-glucagon term uses POW_4 = 0.109.
    # ---------------------------------------------------------------------
    gipins_s   <- 0
    gipins_s0  <- 0
    gipEFFglg  <- (Cgip / bslgip)^pow_4

    # ---------------------------------------------------------------------
    # 8. Combined incretin stimulation of insulin secretion (in T2DM only
    #    the GLP-1 arm contributes).
    # ---------------------------------------------------------------------
    stglc  <- glpins_s + gipins_s
    stglc0 <- glpins_s0 + gipins_s0

    # ---------------------------------------------------------------------
    # 9. Food-effect terms driven by buffer-compartment and end-of-transit
    #    glucose amounts (units 1/mmol * mmol -> unitless).
    # ---------------------------------------------------------------------
    fdglp_s  <- fdglp   * glucose_buffer
    fdglp_s2 <- fdglp_2 * glucose_tr3
    fdgip_s  <- fdgip   * glucose_buffer
    fdglg_s  <- fdglg   * glucose_buffer
    fdins_s  <- fdins   * glucose_buffer  # fdins is FIXED at 0 in this paper

    # ---------------------------------------------------------------------
    # 10. Steady-state baseline kinetics. Each KIN is derived so the
    #     steady-state ODE balance reproduces baseline at t = 0 with no
    #     cotadutide on board.
    # ---------------------------------------------------------------------
    kinglc <- bslglc * (clglc + clgi * bslins)
    kinins <- bslins * clins / (1 + stglc0 * bslglc^glcins_s)
    kinglp <- vmax_glp1 * (bslglp * vcglp) / (km_glp1 + bslglp)
    kinglg <- bslglg * clglg
    kingip <- bslgip * clgip

    # ---------------------------------------------------------------------
    # 11. Derived rate constants for compartmental transfers.
    # ---------------------------------------------------------------------
    kaglc2 <- kaglc * (1 - glpglu_ai)
    k27    <- qglc / vcglc
    k72    <- qglc / vpglc
    k612   <- qgip / vcgip
    k126   <- qgip / vpgip

    # ---------------------------------------------------------------------
    # 12. Lifestyle-change attenuation of endogenous glucose production
    #     (Bosch 2024 Eq 2; TIMEd in days; LSCeff = 1 at t = 0).
    # ---------------------------------------------------------------------
    TIMEd <- t / 24
    LSCeff <- 1 - lsci * klsc * (exp(-klsc * TIMEd) - exp(-kred * TIMEd)) / (kred - klsc)

    # ---------------------------------------------------------------------
    # 13. ODE system (Bosch 2024 supplement S12 DADT block, omitting the
    #     liraglutide arm A(17-18), the AUC accumulator A(8), and the
    #     unused A(11) / A(16)).
    # ---------------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    d/dt(glucose_gut)    <- -kaglc2 * glucose_gut
    d/dt(glucose_buffer) <-  kaglc2 * glucose_gut - kaglc * glucose_buffer - keglc * glucose_buffer
    d/dt(glucose_tr1)    <-  keglc  * glucose_buffer - kelglc * glucose_tr1
    d/dt(glucose_tr2)    <-  kelglc * (glucose_tr1 - glucose_tr2)
    d/dt(glucose_tr3)    <-  kelglc * (glucose_tr2 - glucose_tr3)

    d/dt(glucose)     <-  kaglc * glucose_buffer + kinglc * glgEFFglc * LSCeff -
                          k27 * glucose - clglc * Cglc - clgi * insulin_eff * Cglc +
                          k72 * glucose_per
    d/dt(glucose_per) <-  k27 * glucose - k72 * glucose_per

    d/dt(insulin)      <-  kinins * (1 + fdins_s) * (1 + stglc * Cglc^glcins_s) - clins * Cins
    d/dt(insulin_eff)  <-  ke0ins * (Cins - insulin_eff)

    d/dt(glp1)     <-  kinglp * (1 + fdglp_s + fdglp_s2) * glpEFFglp - vmax_glp1 * glp1 / (km_glp1 + Cglp)
    d/dt(glucagon) <-  kinglg * (1 + fdglg_s) * glpEFFglg * glcEFFglg * gipEFFglg - clglg * Cglg

    d/dt(gip)     <-  kingip * (1 + fdgip_s) - clgip * Cgip - k612 * gip + k126 * gip_per
    d/dt(gip_per) <-  k612 * gip - k126 * gip_per

    # ---------------------------------------------------------------------
    # 14. Steady-state initial conditions (amounts). The glucose state is
    #     anchored on the per-subject FPG covariate; the other four
    #     endogenous species use estimated typical-value baselines.
    # ---------------------------------------------------------------------
    glucose(0)     <- bslglc * vcglc
    glucose_per(0) <- bslglc * vpglc
    insulin(0)     <- bslins * vcins
    insulin_eff(0) <- bslins
    glp1(0)        <- bslglp * vcglp
    glucagon(0)    <- bslglg * vcglg
    gip(0)         <- bslgip * vcgip
    gip_per(0)     <- bslgip * vpgip

    # ---------------------------------------------------------------------
    # 15. Bioavailability of the meal-glucose dose (cotadutide depot is
    #     F = 1 implicitly). MMTT-meal Fglc = 0.334 was estimated in
    #     Bosch 2024 (Table 3).
    # ---------------------------------------------------------------------
    f(glucose_gut) <- fglc

    # ---------------------------------------------------------------------
    # 16. Observed concentrations and residual error.
    # ---------------------------------------------------------------------
    Cglc ~ prop(propSd_Cglc)
    Cins ~ prop(propSd_Cins)
    Cglp ~ prop(propSd_Cglp)
    Cglg ~ prop(propSd_Cglg)
    Cgip ~ prop(propSd_Cgip)
  })
}
