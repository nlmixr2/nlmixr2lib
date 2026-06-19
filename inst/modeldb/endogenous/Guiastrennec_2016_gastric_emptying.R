Guiastrennec_2016_gastric_emptying <- function() {
  description <- paste(
    "Mechanism-based model of postprandial gastric emptying (GE),",
    "cholecystokinin (CCK) plasma kinetics, and gallbladder emptying",
    "(GBE) in response to caloric intake (Guiastrennec 2016). Acetaminophen",
    "(1.5 g oral) is used as a gastric-emptying marker: it shares the",
    "stomach -> upper-small-intestine transit rate KG with the nutrients,",
    "and is absorbed by KA into a two-compartment paracetamol disposition",
    "(CL, VC, Q, VP) where it is observed as plasma Cc. KG is gated by a",
    "Hill-type onset function in time-after-dose (sigmoidicity SIG;",
    "T50OGTT 15.7 min for glucose-only drinks, T50Fat 23.1 min for fat-",
    "containing drinks; T50 = 0 for water) and inhibited by a linear feedback",
    "of caloric content in the upper SI (SLPCAL = -0.0173 1/kcal, 40.7%",
    "stronger in females via SEXF). Per-nutrient amounts (g of fat,",
    "protein, carbohydrate) are tracked through three parallel signal",
    "tracks: an 'upper SI' track that drives the GE caloric feedback",
    "(stomach -> upper_si via KG, drains via KUL with a fixed saturable",
    "MM absorption RAMAX/KM scaled from glucose), a 'duodenum' track that",
    "drives both the CCK-fast (CCKF) Emax stimulus and the GBE Emax",
    "stimulus (stomach -> duodenum via KG, drains via fixed KDJ), and an",
    "'upper jejunum' track that drives the CCK-late (CCKL) linear stimulus",
    "(duodenum -> upper_jejunum via KDJ, drains via fixed KJI). Both",
    "downstream tracks run in parallel from a single stomach (literal",
    "reading of Figure 1; the joint-topology choice is documented in the",
    "validation vignette's Errata). CCKF and CCKL are encoded as paired",
    "precursor + plasma indirect-response models (KRF stimulated by",
    "duodenal-signal Emax; KRL stimulated by jejunal-signal linear; DIS_DIAB",
    "depresses carbohydrate potency POTcarbC by -81.1%). GBE is an",
    "indirect-response model where the duodenal-nutrient signal increases",
    "the gallbladder release rate KRB (Emax with S50_BILE_eff that",
    "scales positively with AGE at +2.15%/yr above reference 58 y, and a",
    "WT-driven baseline gallbladder volume BASEBILE_eff at +1.19%/kg",
    "above reference 88 kg). No recirculation of emptied bile. Observed",
    "outputs are paracetamol plasma concentration (Cc, uM), total CCK",
    "plasma concentration (TCCK = CCKF + CCKL, pM), and gallbladder",
    "volume (GVol, mL); each carries its own residual-error component as",
    "reported in the paper's Table 3."
  )
  reference <- paste(
    "Guiastrennec B, Sonne DP, Hansen M, Bagger JI, Lund A, Rehfeld JF,",
    "Alskar O, Karlsson MO, Vilsboll T, Knop FK, Bergstrand M (2016).",
    "Mechanism-Based Modeling of Gastric Emptying Rate and Gallbladder",
    "Emptying in Response to Caloric Intake.",
    "CPT Pharmacometrics Syst Pharmacol 5(12):692-700.",
    "doi:10.1002/psp4.12152. PMID 28028939; PMCID PMC5192972.",
    sep = " "
  )
  vignette <- "Guiastrennec_2016_gastric_emptying"
  paper_specific_compartments <- c("upper_si", "stom_fat", "stom_prot", "stom_carb", "fat_usi", "prot_usi", "carb_usi", "fat_duod", "prot_duod", "carb_duod", "fat_jej", "prot_jej", "carb_jej", "pool_f", "cckf", "pool_l", "cckl")

  units <- list(
    time          = "min",
    dosing        = "mg (acetaminophen) and g (fats / proteins / carbohydrates) into the per-substrate stomach compartments",
    concentration = "umol/L (acetaminophen Cc); secondary outputs in pmol/L (TCCK) and mL (gallbladder volume GVol)"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline. Linear-deviation effect on the baseline gallbladder volume (BASEBILE_eff = BASEBILE * (1 + 0.0119 * (WT - 88))).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference weight 88 kg, set to the pooled study-population median (Studies A-C: Hansen 89.9, Bagger 89.6, Sonne 86.2 kg). The supplementary control stream that would unambiguously fix the reference is not on disk; see vignette Errata.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age at baseline. Linear-deviation effect on the gallbladder S50_BILE (S50_BILE_eff = S50_BILE * (1 + 0.0215 * (AGE - 58))).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference age 58 years, set to the pooled study-population median (Studies A-C: Hansen 62, Bagger 57, Sonne 60). The supplementary control stream that would unambiguously fix the reference is not on disk; see vignette Errata.",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male). Multiplicative effect on the caloric-feedback slope SLPCAL: SLPCAL_eff = SLPCAL * (1 + SEX_SLPCAL * SEXF) with SEX_SLPCAL = 0.407, meaning the caloric inhibition of gastric emptying is 40.7% stronger in females.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "SEXF = 0 (male)",
      notes              = "Paper uses an FEM flag with FEM = 1 for females; the canonical alias is SEXF.",
      source_name        = "FEM"
    ),
    DIS_DIAB = list(
      description        = "Type-2-diabetes-mellitus indicator (1 = DIS_DIAB patient, 0 = nondiabetic control). Multiplicative -81.1% depression of POTcarbC (the carbohydrate potency on CCK release): POTcarbC_eff = POTcarbC * (1 - 0.811 * DIS_DIAB). All other parameters are common between cohorts.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "DIS_DIAB = 0 (nondiabetic control)",
      notes              = "Same canonical name as NA_NA_paracetamol.R. Paper uses a flag set to 1 for T2D and 0 for matched controls.",
      source_name        = "T2D"
    ),
    DRINK_OGTT = list(
      description        = "Binary indicator that the test drink is glucose-only (oral glucose tolerance test, OGTT). Selects T50OGTT (15.7 min) as the half-onset time of the gastric-emptying delay Hill function. Mutually exclusive with DRINK_FAT.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "DRINK_OGTT = 0 (not a glucose-only drink)",
      notes              = "Set to 1 for Study B's 25 / 75 / 125 g OGTT drinks and the Study C 75 g OGTT arm. Set to 0 for water (Study A) and all fat-containing drinks (Study C low / medium / high fat, Study D medium-high fat).",
      source_name        = "DRINK_OGTT"
    ),
    DRINK_FAT = list(
      description        = "Binary indicator that the test drink contains fat (any nonzero fat content). Selects T50Fat (23.1 min) as the half-onset time of the gastric-emptying delay Hill function. Mutually exclusive with DRINK_OGTT.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "DRINK_FAT = 0 (no fat in the drink)",
      notes              = "Set to 1 for Study C low / medium / high fat and Study D medium-high fat drinks. Set to 0 for water and pure OGTT drinks.",
      source_name        = "DRINK_FAT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 66L,
    n_studies      = 3L,
    age_range      = "Pooled across Studies A-C: 38-75 years (Bagger 38-75, Sonne 42-71, Hansen 44-74); medians ~57-62 y per study",
    age_median     = "58 years (pooled study-population median; reference for AGE-S50_BILE effect)",
    weight_range   = "63.6-122 kg pooled across Studies A-C",
    weight_median  = "88 kg (pooled study-population median; reference for WT-BASEBILE effect)",
    sex_female_pct = 39,
    disease_state  = "33 patients with type 2 diabetes and 33 matched (gender, age, BMI) nondiabetic controls. Cross-over test-drink challenge (water, glucose only at 25/75/125 g, or isocaloric fat-mixed drinks at low/medium/high fat) after an overnight (10 h) fast. An external validation cohort (Study D, 10 nondiabetic controls, medium-high fat drink) was used for the GE / CCK predictive checks but not pooled into the model-development fit.",
    dose_range     = "Acetaminophen 1.5 g oral (1-min infusion into stomach in Study A via a nasogastric tube; bolus in Studies B / C / D). Test-drink macronutrient composition spans 0-125 g carbohydrate, 0-13 g protein, 0-40 g fat (Table 2 of the paper); caloric range 0-506 kcal.",
    regions        = "Denmark (Studies A-D conducted in Copenhagen / Hellerup); fits performed in Sweden (Uppsala).",
    notes          = "Subject demographics, study design, and assessment times come from Table 1 of Guiastrennec 2016. Subject counts and study counts reflect the model-development pool (Studies A + B + C); the Study D external-validation cohort is not pooled in. Pooled medians for WT and AGE were derived from the per-study median ranges; the NONMEM control stream that would unambiguously fix the covariate-effect reference values is in Supplementary Material S1 (not on disk; see vignette Errata)."
  )

  ini({
    # ---------------------------------------------------------------------
    # Acetaminophen (paracetamol) disposition parameters. Paper Table 3
    # column "Gastric emptying model". F1 is fixed at 1 with estimated
    # BSV (22% CV) per the "1 fixed" entry in the table. HL-KA carries a
    # frequentist prior from Hens 2008 (HL-KA_prior = 6.8 +/- 0.9 min);
    # the final estimate is 8.19 min (per Table 3 footnote a). The
    # absorption-rate constant is then KA = log(2) / HL-KA in min^-1.
    # ---------------------------------------------------------------------
    lcl       <- log(0.441)         ; label("Acetaminophen apparent clearance CL/F (L/min)")               # Table 3 GE column: CL/F = 4.41E-01
    lvc       <- log(19.5)          ; label("Acetaminophen central volume VC/F (L)")                       # Table 3 GE column: V_c/F = 1.95E+01
    lq        <- log(1.56)          ; label("Acetaminophen inter-compartmental clearance Q/F (L/min)")     # Table 3 GE column: Q/F = 1.56E+00
    lvp       <- log(48.1)          ; label("Acetaminophen peripheral volume VP/F (L)")                    # Table 3 GE column: V_p/F = 4.81E+01
    lfdepot   <- fixed(log(1))      ; label("Relative acetaminophen bioavailability F1 (FIXED at 1)")      # Table 3 GE column: F_1 = 1 fixed, BSV 22% CV
    lhl_ka    <- log(8.19)          ; label("Acetaminophen absorption half-life HL-KA (min, prior Hens 2008)") # Table 3 GE column: HL-K_A = 8.19 (prior 6.8 +/- 0.9)

    # ---------------------------------------------------------------------
    # Gastric emptying structural parameters. Paper Table 3 "Gastric
    # emptying model" column. KG0 is the baseline first-order emptying
    # rate constant (min^-1); SLPCAL is the linear caloric-feedback
    # slope on KG (units 1/kcal; reported as the unicode minus sign
    # 0.0173 in the trimmed-markdown rendering of the paper table).
    # SEX-SLPCAL is the relative increase in SLPCAL magnitude for females
    # (+40.7%). The GE onset is a Hill function in time-after-dose
    # (sigmoidicity SIG); T50OGTT applies to glucose-only drinks and
    # T50Fat to fat-containing drinks. For water (DRINK_OGTT = DRINK_FAT
    # = 0) the onset gating defaults to 1 throughout, matching the paper's
    # observation that no onset delay was needed for the Study A water
    # arm.
    # ---------------------------------------------------------------------
    lkg0      <- log(1.06)          ; label("Baseline gastric-emptying rate constant KG0 (1/min)")          # Table 3 GE column: K_G0 = 1.06; BSV 155% CV
    lkul      <- log(0.0266)        ; label("Upper- to lower-SI calorie transfer rate KUL (1/min)")         # Table 3 GE column: K_UL = 2.66E-02; BSV 56% CV
    ramax     <- fixed(2.292)       ; label("Maximal nonlinear absorption rate of nutrients RA_MAX (g/min; FIXED from glucose, Hens 2008)")   # Table 3 GE column: RA_MAX = 2.292 fixed
    km        <- fixed(25.12)       ; label("Michaelis constant for nutrient absorption K_M (g equivalent; FIXED from glucose, Hens 2008)")   # Table 3 GE column: K_M = 25.12 fixed
    lslpcal_mag <- log(0.0173)      ; label("Log of the caloric-feedback slope magnitude |SLPCAL| (1/kcal); the sign is fixed negative in model()") # Table 3 GE column: SLP_CAL = 2 0.0173; BSV 19% CV
    sig       <- 0.285              ; label("Sigmoidicity of the gastric-emptying onset Hill function SIG (unitless)") # Table 3 GE column: SIG = 0.285
    lt50ogtt  <- log(15.7)          ; label("Half-onset time of GE for glucose-only (OGTT) drinks T50OGTT (min)") # Table 3 GE column: T_50OGTT = 15.7
    lt50fat   <- log(23.1)          ; label("Half-onset time of GE for fat-containing drinks T50Fat (min)") # Table 3 GE column: T_50Fat = 23.1
    e_sexf_slpcal <- 0.407          ; label("Relative increase in caloric-feedback magnitude SLPCAL for females (+40.7%)")  # Table 3 GE column: SEX-SLP_CAL = 1 40.7%

    # ---------------------------------------------------------------------
    # Cholecystokinin (CCK) sub-model. Paper Table 3 "Cholecystokinin
    # model" column. Two precursor-pool indirect-response models for
    # fast (CCKF) and late (CCKL) plasma CCK; KRF and KRL (derived from
    # RprodF / POOLCCKF and RprodL / POOLCCKL) are released into plasma
    # under stimulation by the duodenal-nutrient and jejunal-nutrient
    # signals respectively. KDJ and KJI are fixed to the literature
    # values cited in the paper (Hens et al. 24, the glucose-absorption
    # reference). DIS_DIAB enters only on POTcarbC via a -81.1% depression.
    # POTfatC is fixed to 100% as the per-paper reference; POTprotC is
    # fixed to 0% (the paper found no statistical support for a protein
    # potency on the fast CCK release).
    # ---------------------------------------------------------------------
    lbase_cckf  <- log(0.506)       ; label("Baseline plasma CCKF concentration BASE_CCKF (pM)")            # Table 3 CCK column: BASE_CCKF = 5.06E-01; BSV 58% CV
    lbase_cckl  <- log(0.0882)      ; label("Baseline plasma CCKL concentration BASE_CCKL (pM)")            # Table 3 CCK column: BASE_CCKL = 8.82E-02
    lpool_cckf  <- log(17.1)        ; label("CCKF precursor-pool size POOLCCKF (pM)")                       # Table 3 CCK column: POOL_CCKF = 1.71E+01
    lpool_cckl  <- log(1.95)        ; label("CCKL precursor-pool size POOLCCKL (pM)")                       # Table 3 CCK column: POOL_CCKL = 1.95; BSV 42% CV
    lkoutf      <- log(0.355)       ; label("CCKF plasma disappearance rate KoutF (1/min)")                 # Table 3 CCK column: K_outF = 3.55E-01; BSV 72% CV
    lkoutl      <- log(0.00546)     ; label("CCKL plasma disappearance rate KoutL (1/min)")                 # Table 3 CCK column: K_outL = 5.46E-03
    smax_cckf   <- 0.193            ; label("Maximal stimulatory effect on CCKF release SMAX-CCKF (unitless)")  # Table 3 CCK column: S_MAX-CCKF = 1.93E-01
    ls50_cckf   <- log(0.592)       ; label("Duodenal-nutrient-signal at half-maximum CCKF stimulation S50-CCKF (g-equivalent)")    # Table 3 CCK column: S_50-CCKF = 5.92E-01
    slp_cckl    <- 0.0377           ; label("Slope of jejunal-nutrient linear stimulation on CCKL release SLP_CCKL (1/g)") # Table 3 CCK column: SLP_CCKL = 3.77E-02
    lkdj        <- fixed(log(0.0833)) ; label("Duodenum-to-jejunum nutrient transfer rate KDJ (1/min; FIXED, Hens 2008)") # Table 3 CCK column: K_DJ = 8.33E-02 fixed
    lkji        <- fixed(log(0.0111)) ; label("Upper- to lower-jejunum nutrient transfer rate KJI (1/min; FIXED, Hens 2008)") # Table 3 CCK column: K_JI = 1.11E-02 fixed
    potfatc     <- fixed(100)       ; label("Fat potency on CCK release POTfatC (%; FIXED reference)")      # Table 3 CCK column: POT_fatC = 100 fixed
    potprotc    <- fixed(0)         ; label("Protein potency on CCK release POTprotC (%; FIXED at 0, not supported)")  # Table 3 CCK column: POT_protC = 0 fixed
    potcarbc    <- 10.1             ; label("Carbohydrate potency on CCK release POTcarbC (%)")             # Table 3 CCK column: POT_carbC = 1.01E+01
    e_t2dm_potcarbc <- -0.811       ; label("DIS_DIAB effect on POTcarbC (multiplicative -81.1%)")              # Table 3 CCK column: T2D-POT_carbC = 2 8.11E-01

    # ---------------------------------------------------------------------
    # Gallbladder emptying (GBE) sub-model. Paper Table 3 "Gallbladder
    # emptying model" column. Indirect-response on gallbladder volume
    # with constant production RprodB and a duodenal-nutrient-driven
    # Emax stimulation of the release rate KRB. WT scales the baseline
    # volume BASEBILE; AGE scales the half-effective signal S50_BILE.
    # POTfatB, POTprotB, POTcarbB are the relative potencies of fats,
    # proteins, and carbohydrates on the gallbladder signal (percent of
    # the fat reference, which is fixed to 100).
    # ---------------------------------------------------------------------
    lbase_bile   <- log(36.3)       ; label("Baseline gallbladder volume BASE_BILE (mL)")                   # Table 3 GBE column: BASE_BILE = 3.63E+01; BSV 27% CV
    lkrb         <- log(0.0618)     ; label("Gallbladder release rate constant K_RB (1/min)")               # Table 3 GBE column: K_RB = 6.18E-02; BSV 90% CV
    smax_bile    <- 6.62            ; label("Maximal gallbladder-release stimulatory effect SMAX-BILE (unitless)")  # Table 3 GBE column: S_MAX-BILE = 6.62
    ls50_bile    <- log(5.52)       ; label("Duodenal-signal at half-maximum gallbladder release S50-BILE (g-equivalent)") # Table 3 GBE column: S_50-BILE = 5.52; BSV 80% CV
    potfatb      <- fixed(100)      ; label("Fat potency on GBE POTfatB (%; FIXED reference)")              # Table 3 GBE column: POT_fatB = 100 fixed
    potprotb     <- 67.9            ; label("Protein potency on GBE POTprotB (%)")                          # Table 3 GBE column: POT_protB = 6.79E+01
    potcarbb     <- 2.25            ; label("Carbohydrate potency on GBE POTcarbB (%)")                     # Table 3 GBE column: POT_carbB = 2.25
    e_wt_base_bile <- 0.0119        ; label("WT-deviation effect on BASEBILE (+1.19%/kg above reference 88 kg)")  # Table 3 GBE column: WT-BASE_BILE = 1 1.19%/kg
    e_age_s50_bile <- 0.0215        ; label("AGE-deviation effect on S50-BILE (+2.15%/yr above reference 58 y)")  # Table 3 GBE column: AGE-S_50-BILE = 1 2.15%/yr

    # ---------------------------------------------------------------------
    # Inter-individual variability. Variances on the log-transformed
    # parameters (omega^2 = log(CV^2 + 1)). Only the parameters that
    # Table 3 reports BSV for are eta-bearing; the rest carry no IIV.
    # All etas are on log-transformed structural parameters (etalvc on
    # lvc, etc.); the un-prefixed SLPCAL also carries BSV per the paper
    # (additive on the linear-scale slope -- Table 3 reports the BSV in
    # CV% and the conversion log(CV^2 + 1) is computed in line as the
    # variance for nlmixr2).
    # ---------------------------------------------------------------------
    etalvc        ~ log(0.68^2 + 1)             # Table 3 GE column: V_c/F BSV = 68% CV -> var = log(1.4624) = 0.380
    etalfdepot    ~ log(0.22^2 + 1)             # Table 3 GE column: F_1   BSV = 22% CV -> var = log(1.0484) = 0.0473
    etalkg0       ~ log(1.55^2 + 1)             # Table 3 GE column: K_G0  BSV = 155% CV -> var = log(3.4025) = 1.224
    etalkul       ~ log(0.56^2 + 1)             # Table 3 GE column: K_UL  BSV = 56% CV  -> var = log(1.3136) = 0.273
    etalslpcal_mag ~ log(0.19^2 + 1)            # Table 3 GE column: SLP_CAL BSV = 19% CV -> var = log(1.0361) = 0.0355
    etalbase_cckf ~ log(0.58^2 + 1)             # Table 3 CCK column: BASE_CCKF BSV = 58% CV -> var = log(1.3364) = 0.290
    etalpool_cckl ~ log(0.42^2 + 1)             # Table 3 CCK column: POOL_CCKL BSV = 42% CV -> var = log(1.1764) = 0.163
    etalkoutf     ~ log(0.72^2 + 1)             # Table 3 CCK column: K_outF   BSV = 72% CV -> var = log(1.5184) = 0.418
    etalbase_bile ~ log(0.27^2 + 1)             # Table 3 GBE column: BASE_BILE BSV = 27% CV -> var = log(1.0729) = 0.0703
    etalkrb       ~ log(0.90^2 + 1)             # Table 3 GBE column: K_RB    BSV = 90% CV -> var = log(1.81)   = 0.593
    etals50_bile  ~ log(0.80^2 + 1)             # Table 3 GBE column: S_50-BILE BSV = 80% CV -> var = log(1.64)  = 0.495

    # ---------------------------------------------------------------------
    # Residual error. Each observation has its own residual-error model
    # per Table 3. Acetaminophen Cc uses a proportional error
    # (14.8% CV); total CCK uses a proportional error (29.5% CV);
    # gallbladder volume uses a combined additive (2.33 mL) +
    # proportional (7.66%) error per the "Add. Err. (mL)" and
    # "Prop. Err. (%)" entries in the GBE column.
    # ---------------------------------------------------------------------
    propSd          <- 0.148            ; label("Proportional residual SD for acetaminophen Cc (fraction)") # Table 3 GE column: Prop. Err = 14.8%
    propSd_TCCK     <- 0.295            ; label("Proportional residual SD for total CCK (fraction)")        # Table 3 CCK column: Prop. Err = 29.5%
    addSd_GVol      <- 2.33             ; label("Additive residual SD for gallbladder volume (mL)")         # Table 3 GBE column: Add. Err. = 2.33 mL
    propSd_GVol     <- 0.0766           ; label("Proportional residual SD for gallbladder volume (fraction)") # Table 3 GBE column: Prop. Err = 7.66%
  })

  model({
    # Individual structural parameters. eta on log-transformed parameters
    # is mu-referenced on the log scale; SLPCAL is on the linear scale
    # so its eta is added on the linear scale (matching Table 3's BSV
    # of 19% CV interpreted as the CV of the linear-scale slope).
    cl       <- exp(lcl) * exp(0)                                    # no IIV on CL/F per Table 3 (no BSV reported)
    vc       <- exp(lvc + etalvc)
    q        <- exp(lq)
    vp       <- exp(lvp)
    fdepot   <- exp(lfdepot + etalfdepot)
    hl_ka    <- exp(lhl_ka)
    ka       <- log(2) / hl_ka                                       # acetaminophen absorption rate (1/min) from prior half-life
    kg0      <- exp(lkg0 + etalkg0)
    kul      <- exp(lkul + etalkul)
    slpcal_mag <- exp(lslpcal_mag + etalslpcal_mag)
    slpcal_i <- -slpcal_mag * (1 + e_sexf_slpcal * SEXF)             # SLPCAL is sign-fixed negative; covariate effect is on its magnitude
    t50ogtt  <- exp(lt50ogtt)
    t50fat   <- exp(lt50fat)

    # CCK sub-model individual parameters.
    base_cckf  <- exp(lbase_cckf + etalbase_cckf)
    base_cckl  <- exp(lbase_cckl)
    pool_cckf_typ <- exp(lpool_cckf)
    pool_cckl_typ <- exp(lpool_cckl + etalpool_cckl)
    koutf      <- exp(lkoutf + etalkoutf)
    koutl      <- exp(lkoutl)
    s50_cckf   <- exp(ls50_cckf)
    kdj        <- exp(lkdj)
    kji        <- exp(lkji)
    potcarbc_i <- potcarbc * (1 + e_t2dm_potcarbc * DIS_DIAB)

    # Derived CCK rates (paper Table 3 footnote b: "derived from other
    # estimated parameters" -- RprodF = BASE_CCKF * KoutF, KRF = RprodF
    # / POOLCCKF; same for the late branch).
    rprodf <- base_cckf * koutf
    rprodl <- base_cckl * koutl
    krf    <- rprodf / pool_cckf_typ
    krl    <- rprodl / pool_cckl_typ

    # GBE sub-model individual parameters with covariate effects.
    base_bile  <- exp(lbase_bile + etalbase_bile) *
                  (1 + e_wt_base_bile * (WT - 88))
    krb        <- exp(lkrb + etalkrb)
    s50_bile   <- exp(ls50_bile + etals50_bile) *
                  (1 + e_age_s50_bile * (AGE - 58))
    rprodb     <- base_bile * krb                                    # baseline-anchored production

    # Acetaminophen disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---------------------------------------------------------------------
    # Gastric-emptying onset Hill function (paper Equation 3). Because
    # the test drink and acetaminophen are co-administered at t = 0 in
    # the published study design, the model uses the raw time variable
    # t as a stand-in for time-after-dose; downstream users dosing at a
    # nonzero time should shift the data accordingly. The branch
    # (1 - DRINK_OGTT - DRINK_FAT) keeps the onset factor at 1 for
    # water (no GE delay was needed for the Study A water arm). Small
    # epsilons (1e-9) prevent division by zero when both binary
    # indicators are 0 or at t = 0.
    # ---------------------------------------------------------------------
    tad_eff <- t
    t50_eff <- t50ogtt * DRINK_OGTT + t50fat * DRINK_FAT + 1e-9
    onset   <- (1 - DRINK_OGTT - DRINK_FAT) +
               (DRINK_OGTT + DRINK_FAT) *
               (tad_eff^sig / (t50_eff^sig + tad_eff^sig + 1e-9))

    # ---------------------------------------------------------------------
    # Caloric content in the upper SI (CALUSI, kcal). The macronutrient
    # energy densities are 9 kcal/g for fats and 4 kcal/g for both
    # proteins and carbohydrates (paper Methods, "GE model").
    # ---------------------------------------------------------------------
    calusi <- 9 * fat_usi + 4 * prot_usi + 4 * carb_usi

    # Effective gastric-emptying rate. Per the paper's Equation 2, the
    # feedback is implemented as a linear modifier on KG0; SLPCAL is
    # signed negative so increasing CALUSI lowers KG. The slpcal_i is
    # clamped to keep the feedback factor non-negative under high
    # caloric inputs (a protective bound for the simulator; the original
    # NONMEM fit operated within data ranges where 1 + slpcal_i * CALUSI
    # remained positive).
    feedback <- 1 + slpcal_i * calusi
    kg_pos   <- (feedback > 0) * feedback
    kg       <- kg0 * kg_pos * onset

    # ---------------------------------------------------------------------
    # Duodenal and jejunal nutrient signals. The signal is the
    # per-macronutrient amount weighted by the macronutrient's relative
    # potency (percent of fat); divided by 100 because the potencies are
    # reported in percent. CCKF and GBE see the duodenal signal; CCKL
    # sees the jejunal signal. POTprotC is fixed at 0 so protein does
    # not enter the CCK signal at all.
    # ---------------------------------------------------------------------
    signal_cckf <- (fat_duod * potfatc + prot_duod * potprotc +
                    carb_duod * potcarbc_i) / 100
    signal_cckl <- (fat_jej  * potfatc + prot_jej  * potprotc +
                    carb_jej  * potcarbc_i) / 100
    signal_bile <- (fat_duod * potfatb + prot_duod * potprotb +
                    carb_duod * potcarbb) / 100

    # Dynamic release rates (KRF and KRL stimulated as in Methods CCK
    # section; KRB stimulated as in GBE section).
    krf_dyn <- krf * (1 + smax_cckf * signal_cckf /
                                  (s50_cckf + signal_cckf + 1e-12))
    krl_dyn <- krl * (1 + slp_cckl * signal_cckl)
    krb_dyn <- krb * (1 + smax_bile * signal_bile /
                                  (s50_bile + signal_bile + 1e-12))

    # Saturable nonlinear absorption of nutrients (paper Methods; fixed
    # to glucose published values for all macronutrients). Applied to
    # each per-substrate compartment in the upper-SI track and in the
    # duodenum / upper-jejunum tracks.
    abs_fat_usi  <- ramax * fat_usi  / (km + fat_usi  + 1e-12)
    abs_prot_usi <- ramax * prot_usi / (km + prot_usi + 1e-12)
    abs_carb_usi <- ramax * carb_usi / (km + carb_usi + 1e-12)
    abs_fat_duod <- ramax * fat_duod / (km + fat_duod + 1e-12)
    abs_prot_duod <- ramax * prot_duod / (km + prot_duod + 1e-12)
    abs_carb_duod <- ramax * carb_duod / (km + carb_duod + 1e-12)
    abs_fat_jej  <- ramax * fat_jej  / (km + fat_jej  + 1e-12)
    abs_prot_jej <- ramax * prot_jej / (km + prot_jej + 1e-12)
    abs_carb_jej <- ramax * carb_jej / (km + carb_jej + 1e-12)

    # ---------------------------------------------------------------------
    # ODE system. Acetaminophen track (stomach -> upper_si -> central +
    # peripheral), three parallel per-macronutrient signal tracks
    # (per-substrate stomach -> upper_si -> drain via KUL, per-substrate
    # stomach -> duodenum via KG -> upper_jejunum via KDJ -> drain via
    # KJI), and the CCK / GBE indirect-response states. Each
    # per-substrate stomach compartment is a separate state because the
    # data carries per-macronutrient gram doses; the joint-model
    # parallel-tracks topology is documented in the validation vignette
    # Errata.
    # ---------------------------------------------------------------------
    # Acetaminophen
    d/dt(stomach)     <- -kg * stomach
    d/dt(upper_si)    <-  kg * stomach - ka * upper_si
    d/dt(central)     <-  ka * upper_si - kel * central -
                          k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Per-macronutrient stomach states (shared KG with acetaminophen).
    d/dt(stom_fat)    <- -kg * stom_fat
    d/dt(stom_prot)   <- -kg * stom_prot
    d/dt(stom_carb)   <- -kg * stom_carb

    # Upper-SI track (drives the GE feedback; drains via KUL plus a
    # saturable MM absorption fixed to the glucose-literature values).
    d/dt(fat_usi)     <-  kg * stom_fat  - kul * fat_usi  - abs_fat_usi
    d/dt(prot_usi)    <-  kg * stom_prot - kul * prot_usi - abs_prot_usi
    d/dt(carb_usi)    <-  kg * stom_carb - kul * carb_usi - abs_carb_usi

    # Duodenum + upper-jejunum track (drives the CCK / GBE signals;
    # transit via KDJ then KJI; saturable MM absorption at each level).
    d/dt(fat_duod)    <-  kg * stom_fat  - kdj * fat_duod  - abs_fat_duod
    d/dt(prot_duod)   <-  kg * stom_prot - kdj * prot_duod - abs_prot_duod
    d/dt(carb_duod)   <-  kg * stom_carb - kdj * carb_duod - abs_carb_duod
    d/dt(fat_jej)     <-  kdj * fat_duod  - kji * fat_jej  - abs_fat_jej
    d/dt(prot_jej)    <-  kdj * prot_duod - kji * prot_jej - abs_prot_jej
    d/dt(carb_jej)    <-  kdj * carb_duod - kji * carb_jej - abs_carb_jej

    # CCK precursor + plasma states. Pool baseline matches the
    # population-typical pool size (POOLCCKF / POOLCCKL); plasma baseline
    # matches BASE_CCKF / BASE_CCKL.
    d/dt(pool_f)      <-  rprodf - krf_dyn * pool_f
    d/dt(cckf)        <-  krf_dyn * pool_f - koutf * cckf
    d/dt(pool_l)      <-  rprodl - krl_dyn * pool_l
    d/dt(cckl)        <-  krl_dyn * pool_l - koutl * cckl
    pool_f(0)         <-  pool_cckf_typ
    cckf(0)           <-  base_cckf
    pool_l(0)         <-  pool_cckl_typ
    cckl(0)           <-  base_cckl

    # Gallbladder volume (open-loop indirect-response model; no
    # recirculation per the paper's GBE Methods).
    d/dt(gallbladder) <-  rprodb - krb_dyn * gallbladder
    gallbladder(0)    <-  base_bile

    # Bioavailability on acetaminophen dose enters the stomach
    # compartment.
    f(stomach) <- fdepot

    # ---------------------------------------------------------------------
    # Observation variables. Acetaminophen plasma concentration in uM
    # (paper LLOQ 0.001 mM = 1 uM; convert from mg amount via molecular
    # weight 151.17 g/mol -> umol = mg / 151.17e-3; concentration uM =
    # umol / vc_L). Total plasma CCK is the sum of CCKF and CCKL
    # ("Assuming similar volume of distribution, the total plasma CCK
    # concentrations were obtained by summing the predicted
    # concentrations of CCKF and CCKL").  Gallbladder volume is the
    # state itself in mL.
    # ---------------------------------------------------------------------
    Cc   <- central / 151.17 / vc * 1000
    TCCK <- cckf + cckl
    GVol <- gallbladder

    Cc   ~ prop(propSd)
    TCCK ~ prop(propSd_TCCK)
    GVol ~ add(addSd_GVol) + prop(propSd_GVol)
  })
}
