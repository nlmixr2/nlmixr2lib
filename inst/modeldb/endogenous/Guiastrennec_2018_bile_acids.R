Guiastrennec_2018_bile_acids <- function() {
  description <- paste(
    "Mechanism-based model of gastric emptying, gallbladder emptying and",
    "refilling, and the enterohepatic circulation (EHC) of endogenous total",
    "bile acids (TBA) after a test drink in healthy adults and adults with",
    "type 2 diabetes (Guiastrennec 2018). The test drink is dosed in kcal",
    "into the stomach; gastric emptying (KG, inhibited by calories in the",
    "upper small intestine, with a logistic onset in time after the drink)",
    "feeds parallel upper-small-intestine and duodenal calorie tracks. The",
    "fat-equivalent amount in the duodenum stimulates gallbladder emptying",
    "linearly, interrupted by a logistic postprandial refilling window. TBA",
    "released from the gallbladder (113 mM) is absorbed in the upper small",
    "intestine or transits a three-compartment lower small intestine to the",
    "liver, where a constant hepatic extraction EH returns it to the",
    "gallbladder and the escaping fraction spills over into a",
    "two-compartment plasma distribution. Outputs are plasma TBA (Cc, umol/L)",
    "and gallbladder volume (GVol, mL). Gastric-emptying parameters are the",
    "typical values of the upstream Guiastrennec 2016 model."
  )
  reference <- paste(
    "Guiastrennec B, Sonne DP, Bergstrand M, Vilsboll T, Knop FK, Karlsson MO",
    "(2018). Model-Based Prediction of Plasma Concentration and",
    "Enterohepatic Circulation of Total Bile Acids in Humans.",
    "CPT Pharmacometrics Syst Pharmacol 7(9):603-612. doi:10.1002/psp4.12325.",
    "PMID 30070437; PMCID PMC6157686. Gastric-emptying parameters from",
    "Guiastrennec B et al. (2016) CPT Pharmacometrics Syst Pharmacol",
    "5(12):692-700, doi:10.1002/psp4.12152."
  )
  vignette <- "Guiastrennec_2018_bile_acids"
  units <- list(
    time = "min",
    dosing = "kcal (energy content of the test drink, dosed into the stomach)",
    concentration = "umol/L (plasma total bile acids); gallbladder volume in mL"
  )

  paper_specific_compartments <- c(
    "upper_si",
    "bile_upper_si",
    "lower_si1",
    "lower_si2",
    "lower_si3"
  )
  paper_specific_etas <- c("etatrefill_off")

  covariateData <- list(
    WT = list(
      description = "Body weight. Allometric scaling of the baseline gallbladder volume and the TBA volumes (exponent 1) and of the hepatic blood flow and TBA distribution clearance (exponent 0.75); also enters the DuBois body-surface area used for the liver volume.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference weight 70 kg (Table 2 footnote b; supplementary NONMEM code Figure S2).",
      source_name = "WT"
    ),
    HT = list(
      description = "Body height. Enters only the DuBois body-surface area, BSA = 0.007184 * WT^0.425 * HT^0.725, from which the liver volume is computed with the Small et al. equation.",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "The source data column HT is in metres and the NONMEM code multiplies it by 100; supply centimetres here.",
      source_name = "HT"
    ),
    AGE = list(
      description = "Age. Linear effect on the slope of the fat-equivalent stimulation of gallbladder emptying, Slopeefat * (1 - 0.025 * (AGE - 64)); also enters the Small et al. liver-volume equation.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference age 64 years (Eq. 5 and supplementary code Figure S2).",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male). Shortens the TBA lower-small-intestine mean transit time by 27.9% (Eq. 6) and strengthens the caloric inhibition of gastric emptying by 40.7% (Guiastrennec 2016 Eq. 4).",
      units = "(binary)",
      type = "binary",
      reference_category = "SEXF = 0 (male)",
      notes = "Source column SEX coded 1 = male, 2 = female; SEXF = SEX - 1. Eq. 6 uses an FEM flag with the same orientation as SEXF.",
      source_name = "SEX"
    ),
    DIS_DIAB = list(
      description = "Type 2 diabetes indicator (1 = T2D, 0 = matched healthy control). Shortens the basal gallbladder-emptying half-life by 48.1% (Eq. 4) and lowers Slopeefat by 32.9% (Eq. 5).",
      units = "(binary)",
      type = "binary",
      reference_category = "DIS_DIAB = 0 (healthy)",
      notes = "Source column T2D (0 = healthy, 1 = T2D). All patients in the cohort had type 2 diabetes.",
      source_name = "T2D"
    ),
    DRINK_OGTT = list(
      description = "Glucose-only test-drink indicator (1 = 75 g oral glucose tolerance test drink, 0 = a fat-containing drink). Shifts the logit of the fraction of TBA absorbed in the upper small intestine by +1.25, and selects the gastric-emptying onset time T50 (15.7 min for OGTT, 23.1 min for fat-containing drinks).",
      units = "(binary)",
      type = "binary",
      reference_category = "DRINK_OGTT = 0 (fat-containing drink)",
      notes = "Source column DRINK == 2. Per-occasion covariate in the crossover design. Water (no calories) was not part of this study, so every DRINK_OGTT = 0 drink is treated as fat-containing for the gastric-emptying onset.",
      source_name = "DRINK"
    ),
    FED_HIGHFAT = list(
      description = "High-fat test-drink indicator (1 = the 40 g fat / 500 kcal drink, 0 = any other test drink). Shifts the logit of the hepatic extraction EH by -0.295.",
      units = "(binary)",
      type = "binary",
      reference_category = "FED_HIGHFAT = 0 (OGTT, low-fat or medium-fat drink)",
      notes = "Source column DRINK == 6. The high-fat drink carries 72% of its energy as fat (Table 1). Per-occasion covariate.",
      source_name = "DRINK"
    ),
    DOSE_FAT_G = list(
      description = "Fat content of the ingested test drink. Enters the fat-equivalent amount of the drink (potency 1 per g) and its energy content (9 kcal/g).",
      units = "g",
      type = "continuous",
      reference_category = NULL,
      notes = "Source column FAT. Table 1: OGTT 0, low fat 2.5, medium fat 10, high fat 40 g.",
      source_name = "FAT"
    ),
    DOSE_CARB_G = list(
      description = "Carbohydrate content of the ingested test drink. Enters the fat-equivalent amount (potency Pcarb per g), the energy content (4 kcal/g) and the TBA lower-small-intestine mean transit time, MTT * (1 + 0.00614 * (DOSE_CARB_G - 93)).",
      units = "g",
      type = "continuous",
      reference_category = NULL,
      notes = "Source column CARB. Reference 93 g (Eq. 6). Table 1: OGTT 75, low fat 107, medium fat 93, high fat 32 g.",
      source_name = "CARB"
    ),
    DOSE_PROT_G = list(
      description = "Protein content of the ingested test drink. Enters the fat-equivalent amount (potency Pprot per g) and the energy content (4 kcal/g).",
      units = "g",
      type = "continuous",
      reference_category = NULL,
      notes = "Hard-coded per drink in the NONMEM code (PROT). Table 1: OGTT 0, low fat 13, medium fat 11, high fat 3 g.",
      source_name = "PROT"
    )
  )

  compartmentData <- list(
    stomach = list(analyte = "test-drink calories", units = "kcal", specimen = "administration site", verified = TRUE),
    upper_si = list(analyte = "test-drink calories", units = "kcal", specimen = "administration site", verified = TRUE),
    duodenum = list(analyte = "test-drink calories", units = "kcal", specimen = "administration site", verified = TRUE),
    gallbladder = list(analyte = "bile (volume)", units = "mL", specimen = "bile", verified = TRUE),
    bile_upper_si = list(
      analyte = "total bile acids",
      units = "umol",
      specimen = "administration site",
      verified = TRUE
    ),
    lower_si1 = list(analyte = "total bile acids", units = "umol", specimen = "administration site", verified = TRUE),
    lower_si2 = list(analyte = "total bile acids", units = "umol", specimen = "administration site", verified = TRUE),
    lower_si3 = list(analyte = "total bile acids", units = "umol", specimen = "administration site", verified = TRUE),
    liver = list(analyte = "total bile acids", units = "umol", specimen = "tissue", verified = TRUE),
    central = list(analyte = "total bile acids", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "total bile acids", units = "umol", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 30L,
    n_studies = 1L,
    age_range = "42-71 years",
    age_mean = "60 years",
    weight_range = "63.6-122 kg",
    weight_mean = "86.2 kg",
    height_range = "1.58-1.97 m (mean 1.76 m)",
    bmi_range = "24.0-31.9 kg/m^2 (mean 27.9)",
    sex_female_pct = 40,
    disease_state = "15 patients with type 2 diabetes and 15 healthy controls matched on sex, age and body-mass index.",
    dose_range = "Four test drinks in a crossover after a 10-hour fast: 75 g OGTT (300 mL, 300 kcal) and three isocaloric (500 kcal), isovolemic (350 mL) drinks with low (2.5 g), medium (10 g) and high (40 g) fat content.",
    regions = "Denmark (Capital Region).",
    notes = "Demographics from Supplementary Table S1 (18 male / 12 female). 920 gallbladder volumes (-10 to 60 min) and 1,236 plasma TBA concentrations (-20 to 240 min). TBA distribution (VC, VP, QD) was fitted separately to digitised radiolabelled cholic-acid data in 8 healthy subjects (Engelking et al.) and fixed."
  )

  ini({
    # --- Gastric emptying (Guiastrennec 2016 typical values) ---------------
    # In the 2018 fit these entered as each subject's empirical Bayes
    # estimates from the 2016 model (supplementary code IT50, IKG0, ISLP,
    # IKUL, ISTEEP). For a new population the 2016 typical values and BSV
    # are used and held fixed.
    lkg0 <- fixed(log(1.06))
    label("Baseline gastric-emptying rate constant KG0 (1/min)") # Guiastrennec 2016 Table 3: KG0 = 1.06
    lkul <- fixed(log(0.0266))
    label("Calorie transfer rate constant from upper to lower small intestine KUL (1/min)") # Guiastrennec 2016 Table 3: KUL = 0.0266
    lslpcal <- fixed(log(0.0173))
    label("Magnitude of the caloric-feedback slope on KG, SLPCAL (1/kcal; sign negative in model)") # Guiastrennec 2016 Table 3: SLPCAL = -0.0173
    e_sexf_slpcal <- fixed(0.407)
    label("Female effect on SLPCAL (fraction)") # Guiastrennec 2016 Table 3: SEX-SLPCAL = +40.7%
    steep <- fixed(0.285)
    label("Steepness of the logistic gastric-emptying onset SIG (1/min)") # Guiastrennec 2016 Table 3 and Eq. 3: SIG = 0.285
    lt50ogtt <- fixed(log(15.7))
    label("Gastric-emptying onset T50 for the OGTT drink (min)") # Guiastrennec 2016 Table 3: T50OGTT = 15.7
    lt50fat <- fixed(log(23.1))
    label("Gastric-emptying onset T50 for fat-containing drinks (min)") # Guiastrennec 2016 Table 3: T50Fat = 23.1
    ramax <- fixed(2.292)
    label("Maximal calorie absorption rate RAMAX (kcal/min)") # Table 2: RAMAX 2.29 fixed; Figure S2 code RAMAXD = 0.573*4
    km <- fixed(25.12)
    label("Calorie amount at half-maximal absorption KM (kcal)") # Table 2: KM 25.1 fixed; Figure S2 code KM = 6.28*4
    lhldj <- log(2.52)
    label("Half-life of calorie transfer from duodenum to jejunum HLDJ (min)") # Table 2: HLDJ = 2.52 min

    # --- Gallbladder emptying ---------------------------------------------
    lbase_gb <- log(29.4)
    label("Baseline (fasting) gallbladder volume for 70 kg (mL)") # Table 2: Base GB = 29.4 mL
    lhlrb <- log(1250)
    label("Half-life of basal gallbladder emptying HLRB (min)") # Table 2: HLRB = 1,250 min
    e_dis_diab_hlrb <- -0.481
    label("Type 2 diabetes effect on HLRB (fraction)") # Table 2: T2D-HLRB = -48.1%; Eq. 4
    lslope_efat <- log(101)
    label("Slope of fat-equivalent stimulation of gallbladder emptying Slopeefat (1/g)") # Table 2: Slopeefat = 101
    e_age_slope_efat <- -0.0250
    label("Age effect on Slopeefat (fraction per year from 64 y)") # Table 2: Age-Slopeefat = -2.50 %/year; Eq. 5
    e_dis_diab_slope_efat <- -0.329
    label("Type 2 diabetes effect on Slopeefat (fraction)") # Table 2: T2D-Slopeefat = -32.9%; Eq. 5
    pfat <- fixed(1)
    label("Potency of fat on gallbladder emptying (fat equivalents per g)") # Table 2: Pfat = 100% fixed
    pcarb <- 0.0574
    label("Potency of carbohydrate on gallbladder emptying (fat equivalents per g)") # Table 2: Pcarb = 5.74%
    pprot <- 0.823
    label("Potency of protein on gallbladder emptying (fat equivalents per g)") # Table 2: Pprot = 82.3%
    slope_refill <- 0.155
    label("Steepness of gallbladder refilling onset and offset Sloperefill (1/min)") # Table 2: Sloperefill = 0.155
    ltrefill_onset <- log(25.6)
    label("Time to 50% of gallbladder refilling onset (min)") # Table 2: Trefill_onset = 25.6 min
    ltrefill_dur <- log(156)
    label("Duration of gallbladder refilling (min)") # Table 2: Trefill_duration = 156 min
    gb_cc <- fixed(113)
    label("Total bile acid concentration in gallbladder bile (mM = umol/mL)") # Table 2: TBA conc. GB = 113 mM fixed

    # --- Enterohepatic circulation ----------------------------------------
    lhltul <- log(59.7)
    label("Half-life of TBA transfer from upper to lower small intestine HLTUL (min)") # Table 2: HLTUL = 59.7 min
    logitfau <- log(0.216 / (1 - 0.216))
    label("Logit of the fraction of TBA absorbed in the upper small intestine FAU (unitless)") # Table 2: FAU = 21.6%
    e_drink_ogtt_fau <- 1.25
    label("OGTT drink shift on the logit of FAU (unitless)") # Table 2: OGTT-FAU = +125; additive on the logit in Figure S2 code (THETA(16))
    lmtt_al <- log(109)
    label("Mean transit time of TBA through the lower small intestine MTTA,L (min)") # Table 2: MTTA,L = 109 min
    e_carb_mtt_al <- 0.00614
    label("Carbohydrate effect on MTTA,L (fraction per g from 93 g)") # Table 2: Carbs-MTTA,L = +0.614 %/g; Eq. 6
    e_sexf_mtt_al <- -0.279
    label("Female effect on MTTA,L (fraction)") # Table 2: Female-MTTA,L = -27.9%; Eq. 6
    logiteh <- log(0.638 / (1 - 0.638))
    label("Logit of the hepatic extraction of TBA EH (unitless)") # Table 2: EH = 63.8%
    e_fed_highfat_eh <- -0.294
    label("High-fat drink shift on the logit of EH (unitless)") # Table 2: High fat-EH = -29.4; additive on the logit in Figure S2 code (THETA(18))
    lqh <- fixed(log(1.5))
    label("Fasting hepatic blood flow for 70 kg QH (L/min)") # Table 2: QH = 1.5 fixed (fasting); Figure S2 code TVQH = 1.5
    e_fed_qh <- fixed(0.3)
    label("Fractional increase in hepatic blood flow after the drink") # Figure S2 code: TVQH * 1.3 (Table 2: 2.0 fixed, fed)
    lvc <- fixed(log(2.839))
    label("TBA central volume for 70 kg VC (L)") # Table 2 and Table in Figure S3: VC = 2.84 fixed; Figure S2 code 2.839
    lvp <- fixed(log(4.072))
    label("TBA peripheral volume for 70 kg VP (L)") # Table 2 and Figure S3: VP = 4.07 fixed; Figure S2 code 4.072
    lq <- fixed(log(0.1732))
    label("TBA intercompartmental clearance for 70 kg QD (L/min)") # Table 2 and Figure S3: QD = 0.173 fixed; Figure S2 code 0.1732

    # --- Between-subject variability --------------------------------------
    # Total variability (BSV + BOV, one eta per parameter), Table S3 final
    # variance-covariance. Order: Base GB, HLRB, Slopeefat, FAU (logit), MTTA,L,
    # EH (logit), Trefill_onset, refill offset time.
    etalbase_gb + etalhlrb + etalslope_efat + etalogitfau + etalmtt_al + etalogiteh + etaltrefill_onset + etatrefill_off ~ c(
      0.0718,
      0.0343, 0.894,
      0.00748, 0.688, 1.04,
      -0.0667, -0.0986, -0.102, 0.390,
      0.00578, -0.252, -0.150, 0.158, 0.222,
      0.0546, -0.216, -0.091, -0.0403, 0.0966, 0.347,
      -0.00433, 0.0599, -0.319, -0.00497, -0.0727, 0.00667, 0.491,
      -0.0144, 0.0338, 0.0711, 0.078, 0.0377, -0.0315, -0.0448, 0.123
    ) # Table S3 (variance-covariance of the etas)
    # Gastric-emptying BSV from Guiastrennec 2016 Table 3 (CV% = sqrt(omega),
    # as in that paper's supplementary code and in Table 2 here).
    etalkg0 ~ fixed(2.4025) # Guiastrennec 2016 Table 3: KG0 BSV 155% -> 1.55^2
    etalkul ~ fixed(0.3136) # Guiastrennec 2016 Table 3: KUL BSV 56% -> 0.56^2
    etalslpcal ~ fixed(0.0361) # Guiastrennec 2016 Table 3: SLPCAL BSV 19% -> 0.19^2

    # --- Residual error ----------------------------------------------------
    propSd <- 0.318
    label("Proportional residual error, plasma TBA (fraction)") # Table 2: Prop. err. 31.8%; Figure S2 SIGMA(3,3) = 0.101 variance
    propSd_GVol <- 0.102
    label("Proportional residual error, gallbladder volume (fraction)") # Table 2: Prop. err. 10.2%; Figure S2 SIGMA(1,1) = 0.0104 variance
    addSd_GVol <- 2.066
    label("Additive residual error, gallbladder volume (mL)") # Table 2: Add. err. 4.27 is the variance (Figure S2 SIGMA(2,2)); SD = sqrt(4.27)
  })

  model({
    # Test-drink composition: energy (Table 1 footnote: 4 kcal/g carbohydrate
    # and protein, 9 kcal/g fat) and fat-equivalent amount (Eq. 1).
    kcal_drink <- 4 * DOSE_CARB_G + 4 * DOSE_PROT_G + 9 * DOSE_FAT_G
    efat_drink <- pfat * DOSE_FAT_G + pprot * DOSE_PROT_G + pcarb * DOSE_CARB_G

    # Gastric emptying (Guiastrennec 2016)
    kg0 <- exp(lkg0 + etalkg0)
    kul <- exp(lkul + etalkul)
    slpcal <- -exp(lslpcal + etalslpcal) * (1 + e_sexf_slpcal * SEXF)
    t50 <- exp(lt50ogtt) * DRINK_OGTT + exp(lt50fat) * (1 - DRINK_OGTT)
    kdj <- log(2) / exp(lhldj)

    # Gallbladder emptying (Eqs. 3-5)
    base_gb <- exp(lbase_gb + etalbase_gb) * (WT / 70)
    hlrb <- exp(lhlrb + etalhlrb) * (1 + e_dis_diab_hlrb * DIS_DIAB)
    krb <- log(2) / hlrb
    slope_efat <- exp(lslope_efat + etalslope_efat) *
      (1 + e_age_slope_efat * (AGE - 64)) *
      (1 + e_dis_diab_slope_efat * DIS_DIAB)
    trefill_onset <- exp(ltrefill_onset + etaltrefill_onset)
    # The eta sits on the offset time (onset + duration), as in the
    # supplementary code (TOFF = (TON + THETA(14)) * EXP(ETA(11))).
    trefill_off <- (trefill_onset + exp(ltrefill_dur)) * exp(etatrefill_off)

    # Enterohepatic circulation (Eq. 6; logit-scale fractions)
    ktul <- log(2) / exp(lhltul)
    logitfau_i <- logitfau + e_drink_ogtt_fau * DRINK_OGTT + etalogitfau
    fau <- expit(logitfau_i)
    kau <- ktul * fau / (1 - fau)
    mtt_al <- exp(lmtt_al + etalmtt_al) *
      (1 + e_sexf_mtt_al * SEXF) *
      (1 + e_carb_mtt_al * (DOSE_CARB_G - 93))
    kal <- 3 / mtt_al
    logiteh_i <- logiteh + e_fed_highfat_eh * FED_HIGHFAT + etalogiteh
    eh <- expit(logiteh_i)
    fh <- 1 - eh

    # Physiology: DuBois BSA (HT in cm), Small et al. liver volume (L),
    # allometric TBA distribution and hepatic blood flow.
    bsa <- 0.007184 * WT^0.425 * HT^0.725
    vh <- 0.822 * bsa^(1.108 * AGE^(-0.022))
    vc <- exp(lvc) * (WT / 70)
    vp <- exp(lvp) * (WT / 70)
    q <- exp(lq) * (WT / 70)^0.75
    qh_fast <- exp(lqh) * (WT / 70)^0.75

    # Time after the test drink. Before the first drink the system is in its
    # fasting steady state: no gastric emptying, no refilling window, fasting
    # hepatic blood flow.
    tdrink <- tad()
    fed <- 1
    lag_ge <- 0
    relax <- 1
    if (is.na(tdrink)) {
      fed <- 0
    } else {
      lag_ge <- expit(steep * (tdrink - t50))
      relax <- (1 - expit(slope_refill * (tdrink - trefill_onset))) +
        expit(slope_refill * (tdrink - trefill_off))
    }
    qh <- qh_fast * (1 + e_fed_qh * fed)

    # Calorie amounts floored at zero for the rate expressions
    # ("compartment safety" in the supplementary code).
    x2 <- upper_si
    if (x2 < 0) x2 <- 0
    x3 <- duodenum
    if (x3 < 0) x3 <- 0
    efat_duod <- 0
    if (kcal_drink > 0) efat_duod <- x3 / kcal_drink * efat_drink

    kg <- kg0 * (1 + slpcal * x2) * lag_ge
    rad <- ramax * x2 / (km + x2)
    rad2 <- ramax * x3 / (km + x3)
    b_rel <- (1 + slope_efat * efat_duod) * relax

    d/dt(stomach) <- -kg * stomach
    d/dt(upper_si) <- kg * stomach - rad - kul * x2
    d/dt(duodenum) <- kg * stomach - rad2 - kdj * x3
    d/dt(gallbladder) <- liver * qh * eh / (vh * gb_cc) - gallbladder * krb * b_rel
    d/dt(bile_upper_si) <- gallbladder * krb * b_rel * gb_cc - bile_upper_si * (kau + ktul)
    d/dt(lower_si1) <- bile_upper_si * ktul - lower_si1 * kal
    d/dt(lower_si2) <- (lower_si1 - lower_si2) * kal
    d/dt(lower_si3) <- (lower_si2 - lower_si3) * kal
    d/dt(liver) <- bile_upper_si * kau + lower_si3 * kal - liver * qh / vh + central * qh / vc
    d/dt(central) <- liver * qh * fh / vh - central * (qh + q) / vc + peripheral1 * q / vp
    d/dt(peripheral1) <- central * q / vc - peripheral1 * q / vp

    # Fasting steady state (supplementary code A_0 block).
    gallbladder(0) <- base_gb
    bile_upper_si(0) <- base_gb * krb * gb_cc / (ktul + kau)
    lower_si1(0) <- base_gb * krb * gb_cc * ktul / ((ktul + kau) * kal)
    lower_si2(0) <- base_gb * krb * gb_cc * ktul / ((ktul + kau) * kal)
    lower_si3(0) <- base_gb * krb * gb_cc * ktul / ((ktul + kau) * kal)
    liver(0) <- base_gb * krb * gb_cc * vh / (qh_fast * eh)
    central(0) <- base_gb * krb * gb_cc * fh * vc / (qh_fast * eh)
    peripheral1(0) <- base_gb * krb * gb_cc * fh * vp / (qh_fast * eh)

    Cc <- central / vc
    GVol <- gallbladder
    Cc ~ prop(propSd)
    GVol ~ add(addSd_GVol) + prop(propSd_GVol)
  })
}
