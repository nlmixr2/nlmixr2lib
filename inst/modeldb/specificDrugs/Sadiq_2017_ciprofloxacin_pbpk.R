Sadiq_2017_ciprofloxacin_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, 13 perfusion-limited tissues, NONMEM).",
    "Ciprofloxacin disposition in 102 adult intensive-care-unit patients,",
    "fitted to plasma concentrations alone by non-linear mixed effects with",
    "frequentist priors (NWPRI) on the eleven tissue-to-plasma partition",
    "coefficients. Lung, brain, heart, skin, muscle, adipose, spleen, gut,",
    "liver, kidney and a lumped rest-of-body are strung between arterial and",
    "venous blood; spleen and gut drain into the liver alongside the hepatic",
    "artery, so the whole splanchnic bed leaves through the hepatic vein.",
    "Every tissue volume and blood flow is an individual function of body",
    "weight and sex rather than a 70 kg reference, and cardiac output is",
    "allometric in weight. Elimination is split into a renal arm driven by",
    "the individual's measured creatinine clearance (glomerular filtration",
    "of unbound drug, augmented by a fitted tubular-secretion factor) and a",
    "non-renal arm driven by the unbound liver concentration; a single",
    "between-subject random effect scales both arms together and a second,",
    "common random effect scales all eleven partition coefficients. The",
    "observation is venous plasma concentration with an additive residual on",
    "the natural-log scale.",
    sep = " "
  )
  reference <- paste(
    "Sadiq MW, Nielsen EI, Khachman D, Conil JM, Georges B, Houin G,",
    "Laffont CM, Karlsson MO, Friberg LE. A whole-body physiologically based",
    "pharmacokinetic (WB-PBPK) model of ciprofloxacin: a step towards",
    "predicting bacterial killing at sites of infection.",
    "J Pharmacokinet Pharmacodyn. 2017;44(2):69-79.",
    "doi:10.1007/s10928-016-9486-9. PMID 27578330. PMCID PMC5376394.",
    "Parameter estimates are Table 1. The structural model - tissue volumes,",
    "blood flows, the mass-balance ODEs and the clearance parameterisation -",
    "is transcribed from the Electronic Supplementary Material",
    "(10928_2016_9486_MOESM1_ESM.pdf, 'NONMEM code WB-PBPK-PD",
    "ciprofloxacin.mod'), blocks $MODEL, $PK, $DES and $ERROR. The renal",
    "clearance form is the article's Equation 3 and the total-clearance split",
    "is Equation 2. See the vignette Errata for the two places where the",
    "deposited control stream and Table 1 disagree.",
    sep = " "
  )
  vignette <- "Sadiq_2017_ciprofloxacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Ciprofloxacin was given as an intravenous infusion into a peripheral
  # vein; the dose therefore enters the venous blood pool. `venous` is
  # outside the auto-detected depot/central set, so it is declared here.
  dosing <- c("venous")

  compartmentData <- list(
    arterial = list(analyte = "ciprofloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    venous = list(analyte = "ciprofloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    lung = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    heart = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    gut = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    liver = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    other = list(analyte = "ciprofloxacin", units = "mg", specimen = "tissue", verified = TRUE),
    a_urine = list(analyte = "ciprofloxacin", units = "mg", specimen = "urine", verified = TRUE),
    a_metabolized = list(analyte = "ciprofloxacin", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives every tissue volume (as a fraction of body weight) and cardiac",
        "output (allometrically, 15 * WT^0.74 L/h). The authors used individual",
        "weights rather than a 70 kg reference patient. Cohort mean 77 +/- 16 kg.",
        sep = " "
      ),
      source_name = "WT"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Selects between two complete sets of fractional tissue volumes and",
        "fractional blood flows (ICRP body composition). VALUES ARE INVERTED",
        "relative to the source: the deposited control stream codes SEX = 0 for",
        "female and SEX = 1 for male, so SEXF = 1 - SEX. The direction is fixed",
        "independently by two blocks of the $PK code - the SEX = 0 branch",
        "carries the higher adipose fraction (0.300 vs 0.171 of body weight) and",
        "lower muscle fraction (0.2916 vs 0.3973), and the matching flow branch",
        "carries higher adipose flow (8.5% vs 5% of cardiac output) and lower",
        "muscle flow (12% vs 17%). Cohort 27 women and 75 men.",
        sep = " "
      ),
      source_name = "SEX"
    ),
    CRCL = list(
      description = "Measured creatinine clearance, raw (NOT BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Measured creatinine clearance in raw mL/min, not BSA-normalized. Enters",
        "the renal clearance arm through the article's Equation 3 with no",
        "centering. TRUNCATED AT 150 mL/min by the deposited control stream",
        "('Cap high values where Cockcroft-Gault less reliable'), so a user",
        "supplying real data must apply the same cap or the term extrapolates",
        "past what was fitted; compare Mouksassi_2015_thrombomodulinAlfa.R,",
        "which truncates at the same value. Cohort 82 +/- 51 mL/min.",
        sep = " "
      ),
      source_name = "CRCL"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 102L,
    n_studies = 1L,
    age_range = "mean 60 +/- 17 years",
    weight_range = "mean 77 +/- 16 kg",
    sex_female_pct = 26.5,
    disease_state = paste(
      "Adults admitted to the intensive care unit for a range of indications,",
      "all mechanically ventilated and receiving ciprofloxacin infusion therapy",
      "during their ICU stay. Treatment duration 3-21 days (average 12).",
      sep = " "
    ),
    renal_function = "Measured creatinine clearance 82 +/- 51 mL/min",
    dose_range = paste(
      "86 patients received 400 mg twice daily as a 1 h infusion, 9 received",
      "400 mg three times daily, 6 received 200 mg twice daily as a 30 min",
      "infusion, and 1 received 600 mg twice daily.",
      sep = " "
    ),
    regions = "France (Toulouse: Hopital Purpan / Hopital Rangueil)",
    notes = paste(
      "588 plasma concentrations, on average 5.8 samples per patient over an",
      "average of 3.1 dosing-interval occasions. Plasma concentration was the",
      "only dependent variable: no tissue was sampled in this study, so every",
      "tissue profile the model produces is a prediction informed by the",
      "literature Kp priors rather than by data. The same dataset had",
      "previously been described by a two-compartment population PK model",
      "(reference 23 of the paper), which is not part of this extraction.",
      sep = " "
    )
  )

  ini({
    # ================= Estimated parameters (Table 1) =================
    # Table 1 reports each estimated parameter twice: a log-scale column
    # with standard error, and a normal-scale column. The two agree
    # throughout (e.g. exp(2.60) = 13.5 for CL_NR), which is what fixes
    # the scale of each row below.
    lcl_nonren <- 2.60
    label("Log non-renal clearance (L/h)") # Table 1 row 'CL NR (l h-1)', log-scale 2.60 +/- 0.14 -> normal 13.5; ESM $THETA 1 '~CLH'

    # Fractional augmentation of glomerular filtration by tubular
    # secretion, entering Equation 3 as (1 + fsec). NATURAL scale, not
    # log: Table 1 places 0.674 in the 'Model estimate normal scale'
    # column and leaves this row's log-scale cell EMPTY, and Equation 3
    # is written with a bare (1 + f_Secretion). See vignette Errata -
    # the deposited $PK block computes RSEC = EXP(THETA(2)), which would
    # make the factor exp(0.674) = 1.96 and is inconsistent with both.
    fsec <- 0.674
    label("Fraction of renal clearance attributable to tubular secretion") # Table 1 row 'f secretion', normal scale 0.674 (26% RSE)

    # Log tissue-to-plasma partition coefficients. Estimated with the
    # literature values as NWPRI informative priors (25% uncertainty on
    # each), which is what made an 11-tissue model identifiable from
    # plasma data alone.
    lkp_lung <- 1.20
    label("Log tissue-to-plasma partition coefficient, lung (dimensionless)") # Table 1 row 'Kp,lung', log-scale 1.20 +/- 0.25 -> normal 3.32 (prior 3.3)
    lkp_brain <- -0.257
    label("Log tissue-to-plasma partition coefficient, brain (dimensionless)") # Table 1 row 'Kp,brain', log-scale -0.257 +/- 0.25 -> normal 0.773 (prior 0.771)
    lkp_heart <- 1.30
    label("Log tissue-to-plasma partition coefficient, heart (dimensionless)") # Table 1 row 'Kp,heart', log-scale 1.30 +/- 0.25 -> normal 3.67 (prior 3.67)
    lkp_skin <- -0.335
    label("Log tissue-to-plasma partition coefficient, skin (dimensionless)") # Table 1 row 'Kp,skin', log-scale -0.335 +/- 0.24 -> normal 0.715 (prior 0.718)
    # Table 1 gives -0.0229 in the log-scale column AND 0.977 in the
    # normal-scale column; exp(-0.0229) = 0.9774, so the two printed
    # columns agree with each other. The deposited $THETA 7 prints
    # +0.0229 without the minus sign. Table 1 is taken as authoritative
    # because it is self-consistent across two columns; see Errata.
    lkp_muscle <- -0.0229
    label("Log tissue-to-plasma partition coefficient, muscle (dimensionless)") # Table 1 row 'Kp,muscle', log-scale -0.0229 +/- 0.16 -> normal 0.977 (prior 1.6)
    lkp_adipose <- -0.885
    label("Log tissue-to-plasma partition coefficient, adipose (dimensionless)") # Table 1 row 'Kp,adipose', log-scale -0.885 +/- 0.23 -> normal 0.413 (prior 0.449)
    lkp_spleen <- 0.668
    label("Log tissue-to-plasma partition coefficient, spleen (dimensionless)") # Table 1 row 'Kp,spleen', log-scale 0.668 +/- 0.25 -> normal 1.95 (prior 1.954)
    lkp_gut <- 1.21
    label("Log tissue-to-plasma partition coefficient, gastrointestinal tract (dimensionless)") # Table 1 row 'Kp,GIT', log-scale 1.21 +/- 0.23 -> normal 3.35 (prior 3.39)
    lkp_liver <- 1.27
    label("Log tissue-to-plasma partition coefficient, liver (dimensionless)") # Table 1 row 'Kp,liver', log-scale 1.27 +/- 0.23 -> normal 3.56 (prior 3.67)
    lkp_kidney <- 2.09
    label("Log tissue-to-plasma partition coefficient, kidney (dimensionless)") # Table 1 row 'Kp,kidney', log-scale 2.09 +/- 0.25 -> normal 8.09 (prior 8.2)
    lkp_other <- 1.35
    label("Log tissue-to-plasma partition coefficient, rest of body (dimensionless)") # Table 1 row 'Kp,rest', log-scale 1.35 +/- 0.11 -> normal 3.86 (prior 2.77)

    # ================= System constants =================
    fu <- fixed(0.65)
    label("Fraction unbound in plasma") # Methods: 'fu,plasma is the fraction unbound of ciprofloxacin in plasma (fu = 0.65)'; ESM $PK 'FUP = 0.65'

    # Cardiac output, allometric in body weight. ESM $PK
    # 'CO = (15*(WT)**(0.74)) ; Cardiac output in L/h'.
    co_coef <- fixed(15)
    label("Cardiac output coefficient (L/h per kg^exponent)") # ESM $PK 'CO = (15*(WT)**(0.74))'
    co_exp <- fixed(0.74)
    label("Cardiac output allometric exponent on body weight") # ESM $PK 'CO = (15*(WT)**(0.74))'

    # Upper bound applied to creatinine clearance before it enters the
    # renal arm. ESM $PK 'IF (CRCL.GT.150) CRCL2=150'.
    crcl_cap <- fixed(150)
    label("Creatinine clearance cap (mL/min)") # ESM $PK 'IF (CRCL.GT.150) CRCL2=150 ; Cap high values where Cockcroft-Gault less reliable'

    # Tissue densities used to turn a mass fraction of body weight into a
    # volume; only skin and adipose carry one in the deposited code.
    dens_skin <- fixed(1.18)
    label("Skin density (kg/L)") # ESM $PK 'VSKN = 0.0383*(WT/1.18)' (female) and '0.045205*(WT/1.18)' (male)
    dens_adipose <- fixed(0.916)
    label("Adipose density (kg/L)") # ESM $PK 'VADI = 0.3*(WT/0.916)' (female) and '0.171233*(WT/0.916)' (male)

    # ---- Fractional tissue volumes, female (ESM $PK, IF (SEX.EQ.0)) ----
    fvol_arterial_f <- fixed(0.017083)
    label("Fractional volume, arterial blood, female (L/kg)") # ESM $PK SEX.EQ.0 'VART = 0.017083*(WT)'
    fvol_venous_f <- fixed(0.051248)
    label("Fractional volume, venous blood, female (L/kg)") # ESM $PK SEX.EQ.0 'VVEN = 0.051248*(WT)'
    fvol_lung_f <- fixed(0.00643836)
    label("Fractional volume, lung, female (L/kg)") # ESM $PK SEX.EQ.0 'VLUN = 0.00643836*(WT)'
    fvol_brain_f <- fixed(0.0191781)
    label("Fractional volume, brain, female (L/kg)") # ESM $PK SEX.EQ.0 'VBRA = 0.0191781*(WT)'
    fvol_heart_f <- fixed(0.004167)
    label("Fractional volume, heart, female (L/kg)") # ESM $PK SEX.EQ.0 'VHRT = 0.004167*(WT)'
    fvol_skin_f <- fixed(0.0383)
    label("Fractional mass, skin, female (kg/kg; divided by skin density)") # ESM $PK SEX.EQ.0 'VSKN = 0.0383*(WT/1.18)'
    fvol_muscle_f <- fixed(0.2916)
    label("Fractional volume, muscle, female (L/kg)") # ESM $PK SEX.EQ.0 'VMUS = 0.2916*(WT)'
    fvol_adipose_f <- fixed(0.3)
    label("Fractional mass, adipose, female (kg/kg; divided by adipose density)") # ESM $PK SEX.EQ.0 'VADI = 0.3*(WT/0.916)'
    fvol_spleen_f <- fixed(0.00247)
    label("Fractional volume, spleen, female (L/kg)") # ESM $PK SEX.EQ.0 'VSPL = 0.00247*(WT)'
    fvol_gut_f <- fixed(0.01644)
    label("Fractional volume, gastrointestinal tract, female (L/kg)") # ESM $PK SEX.EQ.0 'VGIO = 0.01644*(WT) ; stomach+SI colon ICRP'
    fvol_liver_f <- fixed(0.020724)
    label("Fractional volume, liver, female (L/kg)") # ESM $PK SEX.EQ.0 'VHEP = 0.020724*(WT)'
    fvol_kidney_f <- fixed(0.0042466)
    label("Fractional volume, kidney, female (L/kg)") # ESM $PK SEX.EQ.0 'VKID = 0.0042466*(WT)'
    fvol_other_f <- fixed(0.23220194)
    label("Fractional volume, rest of body, female (L/kg)") # ESM $PK SEX.EQ.0 'VRES = 0.23220194*(WT) ; Rest of body volume bones ICRP'

    # ---- Fractional tissue volumes, male (ESM $PK, IF (SEX.EQ.1)) ----
    fvol_arterial_m <- fixed(0.01918)
    label("Fractional volume, arterial blood, male (L/kg)") # ESM $PK SEX.EQ.1 'VART = 0.01918*(WT)'
    fvol_venous_m <- fixed(0.05753)
    label("Fractional volume, venous blood, male (L/kg)") # ESM $PK SEX.EQ.1 'VVEN = 0.05753*(WT)'
    fvol_lung_m <- fixed(0.00643836)
    label("Fractional volume, lung, male (L/kg)") # ESM $PK SEX.EQ.1 'VLUN = 0.00643836*(WT)'
    fvol_brain_m <- fixed(0.01918)
    label("Fractional volume, brain, male (L/kg)") # ESM $PK SEX.EQ.1 'VBRA = 0.01918*(WT)'
    fvol_heart_m <- fixed(0.004521)
    label("Fractional volume, heart, male (L/kg)") # ESM $PK SEX.EQ.1 'VHRT = 0.004521*(WT)'
    fvol_skin_m <- fixed(0.045205)
    label("Fractional mass, skin, male (kg/kg; divided by skin density)") # ESM $PK SEX.EQ.1 'VSKN = 0.045205*(WT/1.18)'
    fvol_muscle_m <- fixed(0.3973)
    label("Fractional volume, muscle, male (L/kg)") # ESM $PK SEX.EQ.1 'VMUS = 0.3973*(WT)'
    fvol_adipose_m <- fixed(0.171233)
    label("Fractional mass, adipose, male (kg/kg; divided by adipose density)") # ESM $PK SEX.EQ.1 'VADI = 0.171233*(WT/0.916)'
    fvol_spleen_m <- fixed(0.00247)
    label("Fractional volume, spleen, male (L/kg)") # ESM $PK SEX.EQ.1 'VSPL = 0.00247*(WT)'
    fvol_gut_m <- fixed(0.01644)
    label("Fractional volume, gastrointestinal tract, male (L/kg)") # ESM $PK SEX.EQ.1 'VGIO = 0.01644*(WT) ; stomach+SI colon ICRP'
    fvol_liver_m <- fixed(0.02466)
    label("Fractional volume, liver, male (L/kg)") # ESM $PK SEX.EQ.1 'VHEP = 0.02466*(WT)'
    fvol_kidney_m <- fixed(0.004247)
    label("Fractional volume, kidney, male (L/kg)") # ESM $PK SEX.EQ.1 'VKID = 0.004247*(WT)'
    fvol_other_m <- fixed(0.23612164)
    label("Fractional volume, rest of body, male (L/kg)") # ESM $PK SEX.EQ.1 'VRES = 0.23612164*(WT) ; Rest of body volume bones ICRP'

    # ---- Fractional blood flows, female (ESM $PK, IF (SEX.EQ.0)) ----
    # The ten arterial-side shares below sum to exactly 1.000000.
    fq_brain_f <- fixed(0.12)
    label("Fractional blood flow, brain, female (fraction of cardiac output)") # ESM $PK SEX.EQ.0 'QBRA = CO*(0.12)'
    fq_heart_f <- fixed(0.05)
    label("Fractional blood flow, heart, female (fraction of cardiac output)") # ESM $PK SEX.EQ.0 'QHRT = CO*(0.05)'
    fq_skin_f <- fixed(0.05)
    label("Fractional blood flow, skin, female (fraction of cardiac output)") # ESM $PK SEX.EQ.0 'QSKN = CO*(0.05)'
    fq_muscle_f <- fixed(0.12)
    label("Fractional blood flow, muscle, female (fraction of cardiac output)") # ESM $PK SEX.EQ.0 'QMUS = CO*(0.12)'
    fq_adipose_f <- fixed(0.085)
    label("Fractional blood flow, adipose, female (fraction of cardiac output)") # ESM $PK SEX.EQ.0 'QADI = CO*(0.085)'
    fq_spleen_f <- fixed(0.03)
    label("Fractional blood flow, spleen, female (fraction of cardiac output)") # ESM $PK SEX.EQ.0 'QSPL = CO*(0.03)'
    fq_gut_f <- fixed(0.16)
    label("Fractional blood flow, gastrointestinal tract, female (fraction of cardiac output)") # ESM $PK SEX.EQ.0 'QGIO = CO*(0.16)'
    fq_kidney_f <- fixed(0.17)
    label("Fractional blood flow, kidney, female (fraction of cardiac output)") # ESM $PK SEX.EQ.0 'QKID = CO*(0.17)'
    fq_hepatic_artery_f <- fixed(0.065)
    label("Fractional blood flow, hepatic artery, female (fraction of cardiac output)") # ESM $PK SEX.EQ.0 'QHEPA= CO*(0.065)'
    fq_other_f <- fixed(0.15)
    label("Fractional blood flow, rest of body, female (fraction of cardiac output)") # ESM $PK SEX.EQ.0 'QRES = CO*(0.15)'

    # ---- Fractional blood flows, male (ESM $PK, IF (SEX.EQ.1)) ----
    # The ten arterial-side shares below also sum to exactly 1.000000.
    fq_brain_m <- fixed(0.12)
    label("Fractional blood flow, brain, male (fraction of cardiac output)") # ESM $PK SEX.EQ.1 'QBRA = CO*(0.12)'
    fq_heart_m <- fixed(0.04)
    label("Fractional blood flow, heart, male (fraction of cardiac output)") # ESM $PK SEX.EQ.1 'QHRT = CO*(0.04)'
    fq_skin_m <- fixed(0.05)
    label("Fractional blood flow, skin, male (fraction of cardiac output)") # ESM $PK SEX.EQ.1 'QSKN = CO*(0.05)'
    fq_muscle_m <- fixed(0.17)
    label("Fractional blood flow, muscle, male (fraction of cardiac output)") # ESM $PK SEX.EQ.1 'QMUS = CO*(0.17)'
    fq_adipose_m <- fixed(0.05)
    label("Fractional blood flow, adipose, male (fraction of cardiac output)") # ESM $PK SEX.EQ.1 'QADI = CO*(0.05)'
    fq_spleen_m <- fixed(0.03)
    label("Fractional blood flow, spleen, male (fraction of cardiac output)") # ESM $PK SEX.EQ.1 'QSPL = CO*(0.03)'
    fq_gut_m <- fixed(0.15)
    label("Fractional blood flow, gastrointestinal tract, male (fraction of cardiac output)") # ESM $PK SEX.EQ.1 'QGIO = CO*(0.15)'
    fq_kidney_m <- fixed(0.19)
    label("Fractional blood flow, kidney, male (fraction of cardiac output)") # ESM $PK SEX.EQ.1 'QKID = CO*(0.19)'
    fq_hepatic_artery_m <- fixed(0.065)
    label("Fractional blood flow, hepatic artery, male (fraction of cardiac output)") # ESM $PK SEX.EQ.1 'QHEPA= CO*(0.065)'
    fq_other_m <- fixed(0.135)
    label("Fractional blood flow, rest of body, male (fraction of cardiac output)") # ESM $PK SEX.EQ.1 'QRES = CO*(0.135)'

    # ================= Random effects =================
    # One random effect on clearance: the deposited $PK applies EXP(ETA(1))
    # to the non-renal arm AND to the renal arm, so it scales total
    # clearance rather than either component alone.
    etalcl ~ 0.316 # ESM $OMEGA 1 '~IIV_CLH' = 0.316; Table 1 row 'IIV CL (CV %)' = 56 (9.3% RSE), and sqrt(0.316) = 0.562
    # One random effect shared by all eleven partition coefficients
    # ('a common value for all', Results).
    etalkp ~ 0.306 # ESM $OMEGA 2 '~IIV_KP' = 0.306; Table 1 row 'IIV Kp (CV %)' = 55 (15% RSE), and sqrt(0.306) = 0.553

    # Additive residual on the natural-log scale. Methods: 'All plasma
    # concentrations from patients were transformed into natural
    # logarithms before the data analysis. An additive error model was
    # used on log-transformed data.' The deposited $ERROR confirms it:
    # IPRED = LOG(IPRD2) and Y = IPRED + EPS(1).
    expSd <- 0.334664
    label("Additive residual standard deviation on the natural-log scale") # ESM $SIGMA 1 '~ Prop RES_ERR' = 0.112; sqrt(0.112) = 0.334664; Table 1 row 'Proportional residual error (%)' = 33 (7.1% RSE)
  })

  model({
    # ================= Individual physiology =================
    # Every volume and flow is an individual function of body weight and
    # sex. SEXF = 1 selects the female branch of the deposited $PK code
    # (its SEX.EQ.0 branch); SEXF = 0 selects the male branch.
    co <- co_coef * WT^co_exp # ESM $PK 'CO = (15*(WT)**(0.74)) ; Cardiac output in L/h'

    # Tissue volumes (L). Skin and adipose fractions are masses and are
    # divided by a tissue density to become volumes.
    v_arterial <- WT * (SEXF * fvol_arterial_f + (1 - SEXF) * fvol_arterial_m)
    v_venous <- WT * (SEXF * fvol_venous_f + (1 - SEXF) * fvol_venous_m)
    v_lung <- WT * (SEXF * fvol_lung_f + (1 - SEXF) * fvol_lung_m)
    v_brain <- WT * (SEXF * fvol_brain_f + (1 - SEXF) * fvol_brain_m)
    v_heart <- WT * (SEXF * fvol_heart_f + (1 - SEXF) * fvol_heart_m)
    v_skin <- WT / dens_skin * (SEXF * fvol_skin_f + (1 - SEXF) * fvol_skin_m)
    v_muscle <- WT * (SEXF * fvol_muscle_f + (1 - SEXF) * fvol_muscle_m)
    v_adipose <- WT / dens_adipose * (SEXF * fvol_adipose_f + (1 - SEXF) * fvol_adipose_m)
    v_spleen <- WT * (SEXF * fvol_spleen_f + (1 - SEXF) * fvol_spleen_m)
    v_gut <- WT * (SEXF * fvol_gut_f + (1 - SEXF) * fvol_gut_m)
    v_liver <- WT * (SEXF * fvol_liver_f + (1 - SEXF) * fvol_liver_m)
    v_kidney <- WT * (SEXF * fvol_kidney_f + (1 - SEXF) * fvol_kidney_m)
    v_other <- WT * (SEXF * fvol_other_f + (1 - SEXF) * fvol_other_m)

    # Tissue blood flows (L/h). Lung, arterial and venous all carry the
    # whole cardiac output.
    q_lung <- co # ESM $PK 'QLUN = CO'
    q_arterial <- co # ESM $PK 'QART = CO'
    q_venous <- co # ESM $PK 'QVEN = CO'
    q_brain <- co * (SEXF * fq_brain_f + (1 - SEXF) * fq_brain_m)
    q_heart <- co * (SEXF * fq_heart_f + (1 - SEXF) * fq_heart_m)
    q_skin <- co * (SEXF * fq_skin_f + (1 - SEXF) * fq_skin_m)
    q_muscle <- co * (SEXF * fq_muscle_f + (1 - SEXF) * fq_muscle_m)
    q_adipose <- co * (SEXF * fq_adipose_f + (1 - SEXF) * fq_adipose_m)
    q_spleen <- co * (SEXF * fq_spleen_f + (1 - SEXF) * fq_spleen_m)
    q_gut <- co * (SEXF * fq_gut_f + (1 - SEXF) * fq_gut_m)
    q_kidney <- co * (SEXF * fq_kidney_f + (1 - SEXF) * fq_kidney_m)
    q_hepatic_artery <- co * (SEXF * fq_hepatic_artery_f + (1 - SEXF) * fq_hepatic_artery_m)
    q_other <- co * (SEXF * fq_other_f + (1 - SEXF) * fq_other_m)

    # Total hepatic outflow is the hepatic artery plus the two portal
    # tributaries. ESM $PK 'QHEPT= QHEPA + QSPL + QGIO'.
    q_hepatic <- q_hepatic_artery + q_spleen + q_gut

    # ================= Individual parameters =================
    # Partition coefficients share one random effect.
    kp_lung <- exp(lkp_lung + etalkp)
    kp_brain <- exp(lkp_brain + etalkp)
    kp_heart <- exp(lkp_heart + etalkp)
    kp_skin <- exp(lkp_skin + etalkp)
    kp_muscle <- exp(lkp_muscle + etalkp)
    kp_adipose <- exp(lkp_adipose + etalkp)
    kp_spleen <- exp(lkp_spleen + etalkp)
    kp_gut <- exp(lkp_gut + etalkp)
    kp_liver <- exp(lkp_liver + etalkp)
    kp_kidney <- exp(lkp_kidney + etalkp)
    kp_other <- exp(lkp_other + etalkp)

    # Renal clearance, article Equation 3: CL_R = CRCL x fu,plasma x
    # (1 + f_Secretion). The creatinine clearance is capped first and
    # converted from mL/min to L/h by 60/1000. The same random effect
    # that scales non-renal clearance also scales this arm.
    crcl_capped <- min(CRCL, crcl_cap) # ESM $PK 'CRCL2=CRCL' then 'IF (CRCL.GT.150) CRCL2=150'
    cl_renal <- crcl_capped * 60 / 1000 * fu * (1 + fsec) * exp(etalcl) # Eq 3; ESM $PK 'CLR = (((CRCL2*60/1000)*FUP)*(1+RSEC))*EXP(ETA(1))'

    # Non-renal (hepatic) clearance, acting on the unbound liver
    # concentration in the mass balance below.
    cl_nonren <- exp(lcl_nonren + etalcl) # ESM $PK 'TVCLH= EXP(THETA(1))' and 'CLH = TVCLH*EXP(ETA(1))'

    # Total clearance, article Equation 2: CL = CL_R + CL_NR. Reported
    # for diagnostics; the ODEs use the two arms separately.
    cl <- cl_renal + cl_nonren # Eq 2; ESM $PK 'CL = CLR+CLH'

    # ================= Tissue concentrations (mg/L) =================
    # ESM $DES 'C1 = A(1)/VART' ... 'C13 = A(13)/VRES'.
    c_arterial <- arterial / v_arterial
    c_venous <- venous / v_venous
    c_lung <- lung / v_lung
    c_brain <- brain / v_brain
    c_heart <- heart / v_heart
    c_skin <- skin / v_skin
    c_muscle <- muscle / v_muscle
    c_adipose <- adipose / v_adipose
    c_spleen <- spleen / v_spleen
    c_gut <- gut / v_gut
    c_liver <- liver / v_liver
    c_kidney <- kidney / v_kidney
    c_other <- other / v_other

    # ================= Mass-balance abbreviations =================
    # Blood returning to the venous pool. Spleen and gut are absent
    # because they drain into the liver, which returns the whole
    # splanchnic bed at the total hepatic flow.
    venin1 <- (c_brain * q_brain / kp_brain) + (c_heart * q_heart / kp_heart) +
      (c_skin * q_skin / kp_skin) + (c_muscle * q_muscle / kp_muscle) # ESM $DES 'VENIN1 = ...'
    venin2 <- (c_adipose * q_adipose / kp_adipose) + (c_liver * q_hepatic / kp_liver) +
      (c_kidney * q_kidney / kp_kidney) + (c_other * q_other / kp_other) # ESM $DES 'VENIN2 = ...'
    venout <- q_lung * c_venous # ESM $DES 'VENOUT = (QLUN*C2)'

    # Liver inflow is the hepatic artery plus the two portal tributaries;
    # outflow is the hepatic vein plus non-renal elimination of the
    # unbound liver concentration.
    hepin <- (c_arterial * q_hepatic_artery) + (c_spleen * q_spleen / kp_spleen) +
      (c_gut * q_gut / kp_gut) # ESM $DES 'HEPIN = (C1*QHEPA)+(C9*QSPL/KSPL)+(C10*QGIO/KGIO)'
    # Non-renal elimination rate (mg/h), carried separately so it can
    # feed both the liver ODE and the cumulative-metabolism integrator.
    r_nonren <- c_liver * cl_nonren * fu / kp_liver # ESM $DES 'HEPOUT' second term '(C11*CLH*FUP/KHEP)'
    hepout <- (c_liver * q_hepatic / kp_liver) + r_nonren # ESM $DES 'HEPOUT = (C11*QHEPT/KHEP)+(C11*CLH*FUP/KHEP)'

    # Renal elimination rate (mg/h). Note that the deposited code drives
    # it with the ARTERIAL concentration, not the kidney concentration.
    r_renal <- c_arterial * cl_renal # ESM $DES 'DADT(12)' third term '-C1*CLR'

    # ================= Mass balance (mg/h) =================
    d/dt(arterial) <- (q_lung * c_lung / kp_lung) - (q_arterial * c_arterial) # ESM $DES 'DADT(1) = (QLUN*C3/KLUN)-(QART*C1)'
    d/dt(venous) <- venin1 + venin2 - venout # ESM $DES 'DADT(2) = VENIN1+VENIN2-VENOUT'
    d/dt(lung) <- (q_venous * c_venous) - (q_lung * c_lung / kp_lung) # ESM $DES 'DADT(3) = (QVEN*C2)-(QLUN*C3/KLUN)'
    d/dt(brain) <- q_brain * (c_arterial - c_brain / kp_brain) # ESM $DES 'DADT(4) = (QBRA*C1)-(QBRA*C4/KBRA)'
    d/dt(heart) <- q_heart * (c_arterial - c_heart / kp_heart) # ESM $DES 'DADT(5) = (QHRT*C1)-(QHRT*C5/KHRT)'
    d/dt(skin) <- q_skin * (c_arterial - c_skin / kp_skin) # ESM $DES 'DADT(6) = (QSKN*C1)-(QSKN*C6/KSKN)'
    d/dt(muscle) <- q_muscle * (c_arterial - c_muscle / kp_muscle) # ESM $DES 'DADT(7) = (QMUS*C1)-(QMUS*C7/KMUS)'
    d/dt(adipose) <- q_adipose * (c_arterial - c_adipose / kp_adipose) # ESM $DES 'DADT(8) = (QADI*C1)-(QADI*C8/KADI)'
    d/dt(spleen) <- q_spleen * (c_arterial - c_spleen / kp_spleen) # ESM $DES 'DADT(9) = (QSPL*C1)-(QSPL*C9/KSPL)'
    d/dt(gut) <- q_gut * (c_arterial - c_gut / kp_gut) # ESM $DES 'DADT(10)= (QGIO*C1)-(QGIO*C10/KGIO)'
    d/dt(liver) <- hepin - hepout # ESM $DES 'DADT(11)= HEPIN-HEPOUT'
    d/dt(kidney) <- q_kidney * (c_arterial - c_kidney / kp_kidney) - r_renal # ESM $DES 'DADT(12)= (QKID*C1)-(QKID*C12/KKID)-C1*CLR'
    d/dt(other) <- q_other * (c_arterial - c_other / kp_other) # ESM $DES 'DADT(13)= (QRES*C1)-(QRES*C13/KRES)'

    # The deposited code carries a single eliminated-drug integrator,
    # 'DADT(14)= C1*CLR+(C11*CLH*FUP/KHEP)'. It is split here into its
    # two printed terms so the renal and non-renal routes can be read
    # separately; their sum reproduces DADT(14) exactly.
    d/dt(a_urine) <- r_renal # ESM $DES 'DADT(14)' first term 'C1*CLR'
    d/dt(a_metabolized) <- r_nonren # ESM $DES 'DADT(14)' second term '(C11*CLH*FUP/KHEP)'

    # ================= Observation =================
    # Venous plasma concentration. The deposited $ERROR observes
    # compartment 2 as IPRD2 = A(2)/VVEN, with the residual additive on
    # the natural-log scale, which is lnorm() in nlmixr2.
    Cc <- c_venous
    Cc ~ lnorm(expSd)
  })
}
