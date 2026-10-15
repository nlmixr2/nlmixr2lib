Sier_2017_estradiol_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, COPASI). Model of endogenous and exogenous",
    "17-beta-estradiol (E2)",
    "in women (base PBPK variant of Sier 2017). Venous and arterial blood with",
    "explicit reversible E2 binding to SHBG and albumin, 16 tissue",
    "compartments (well-mixed except permeability-limited liver and",
    "uterus+ovaries 'gonads'), estrogen-receptor binding in liver and gonad",
    "cells, zero-order E2 biosynthesis in the gonads, first-order hepatic",
    "intrinsic clearance scaled to liver weight and first-order extra-hepatic",
    "clearance from the kidney. Organ volumes and blood flows follow the Bosgra",
    "2012 anthropometric equations in body weight, height and age. Oral doses",
    "enter a dose compartment that empties to the intestine; i.v. boluses enter",
    "a dose compartment that empties to venous blood at 250 1/h. States start",
    "at the deposited endogenous steady state of the 60 kg, 1.63 m, 35-year",
    "reference woman. Observation Cc is FREE (unbound) venous E2. Sibling",
    "model with a mechanistic liver-metabolism network in place of the",
    "intrinsic clearance: modellib('Sier_2017_estradiol_liverode_pbpk')."
  )
  reference <- paste(
    "Sier JH, Thumser AE, Plant NJ. Linking physiologically-based",
    "pharmacokinetic and genome-scale metabolic networks to understand",
    "estradiol biology. BMC Syst Biol. 2017;11:141.",
    "doi:10.1186/s12918-017-0520-3. Executable model deposited by the authors",
    "as BioModels MODEL1711210001 (E2_PBPK_v7.3)."
  )
  vignette <- "Sier_2017_estradiol"
  # Every ODE state is a concentration (nM). The deposited COPASI model sets
  # every compartment size to 1, so a dose in nmol entered into a dose
  # compartment is numerically its starting concentration; the paper enters
  # 0.3 mg i.v. as 1101 nmol and 2 / 4 mg oral as 7342 / 14685 nmol.
  units <- list(time = "h", dosing = "nmol", concentration = "nM")

  # Two dose routes: oral into `depot` (emptying to the intestine) and i.v.
  # bolus into `depot_iv` (emptying to venous blood).
  dosing <- c("depot", "depot_iv")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives body surface area (DuBois), body mass index, adipose volume",
        "(Deurenberg) and, through the 93% lean-mass balance, muscle volume.",
        "The deposited reference woman weighs 60 kg; the paper's virtual",
        "population draws weight from NHANES 2013-14 women aged 18-45."
      ),
      source_name = "Body weight (kg)"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The source equations take height in metres; model() converts",
        "ht_m = HT / 100. Drives most organ volumes (Bosgra 2012",
        "ln(height) regressions) and body surface area. Deposited reference",
        "height 1.63 m (163 cm)."
      ),
      source_name = "Height (m)"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives adipose, brain, gonad and thymus volumes. Deposited",
        "reference age 35 years; virtual population ages 18-45 years."
      ),
      source_name = "Age (years)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "estradiol", units = "nmol", specimen = "administration site", verified = TRUE),
    depot_iv = list(analyte = "estradiol", units = "nmol", specimen = "administration site", verified = TRUE),
    venous = list(analyte = "estradiol (free)", units = "nM", specimen = "plasma", verified = TRUE),
    shbg_venous = list(analyte = "SHBG (free)", units = "nM", specimen = "plasma", verified = TRUE),
    complex_shbg_venous = list(analyte = "estradiol-SHBG complex", units = "nM", specimen = "plasma", verified = TRUE),
    albumin_venous = list(analyte = "albumin (free)", units = "nM", specimen = "plasma", verified = TRUE),
    complex_albumin_venous = list(
      analyte = "estradiol-albumin complex",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    arterial = list(analyte = "estradiol (free)", units = "nM", specimen = "plasma", verified = TRUE),
    shbg_arterial = list(analyte = "SHBG (free)", units = "nM", specimen = "plasma", verified = TRUE),
    complex_shbg_arterial = list(
      analyte = "estradiol-SHBG complex",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    albumin_arterial = list(analyte = "albumin (free)", units = "nM", specimen = "plasma", verified = TRUE),
    complex_albumin_arterial = list(
      analyte = "estradiol-albumin complex",
      units = "nM",
      specimen = "plasma",
      verified = TRUE
    ),
    lung = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    bone = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    brain = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    heart = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    intestine = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    spleen = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    stomach = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    pancreas = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    thymus = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    kidney = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    other = list(analyte = "estradiol", units = "nM", specimen = "tissue", verified = TRUE),
    vp_liver = list(analyte = "estradiol", units = "nM", specimen = "whole blood", verified = TRUE),
    liver = list(analyte = "estradiol (free)", units = "nM", specimen = "tissue", verified = TRUE),
    target_er_liver = list(analyte = "estrogen receptor (free)", units = "nM", specimen = "tissue", verified = TRUE),
    complex_er_liver = list(
      analyte = "estradiol-estrogen receptor complex",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    ),
    vp_gonads = list(analyte = "estradiol", units = "nM", specimen = "whole blood", verified = TRUE),
    gonads = list(analyte = "estradiol (free)", units = "nM", specimen = "tissue", verified = TRUE),
    target_er_gonads = list(analyte = "estrogen receptor (free)", units = "nM", specimen = "tissue", verified = TRUE),
    complex_er_gonads = list(
      analyte = "estradiol-estrogen receptor complex",
      units = "nM",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1495L,
    n_studies = 2L,
    age_range = "18-45 years (virtual population)",
    sex_female_pct = 100,
    disease_state = "healthy women (endogenous estradiol plus single exogenous doses)",
    dose_range = "Single 0.3 mg i.v. bolus; single 2 mg or 4 mg oral dose",
    regions = "United States (NHANES anthropometrics)",
    notes = paste(
      "Not a population fit: a deterministic mechanistic model with no",
      "inter-individual variability or residual error. The 1495 subjects are",
      "a virtual population of women aged 18-45 years whose weight and height",
      "come from NHANES 2013-14 (Methods, 'Virtual patient population';",
      "Additional file 1 Figure S7). Model behaviour was compared against",
      "published clinical means: Kuhnz 1993 (0.3 mg i.v. bolus in young",
      "women) and Lyrenas 1981 (2 and 4 mg oral, plasma followed for 50 h).",
      "Physiology follows Bosgra 2012; estradiol-specific constants follow",
      "Plowchalk and Teeguarden 2002 except where the paper states otherwise."
    )
  )

  ini({
    # All values are fixed point values from the paper / Additional file 1 /
    # the authors' BioModels deposit (MODEL1711210001). Nothing was estimated
    # against data in the sense of a population fit.

    # E2 input -----------------------------------------------------------
    ksyn_gonads <- fixed(2)
    label("E2 biosynthesis rate into gonad cells (nM/h)") # BioModels MODEL1711210001 reaction 'Gonads: E2 Biosynthesis' k1 = 2 times the fixed 1 nM Biosynthesis species; curator note 'E2 biosynthesis rate constant = 2'
    kiv <- fixed(250)
    label("First-order delivery rate constant from the i.v. dose compartment to venous blood (1/h)") # Main text Results: rate constant of 250 h-1; Additional file 1 v12 K1 = 250 h-1
    lka <- fixed(log(0.025))
    label("First-order absorption rate constant from the oral dose compartment to the intestine (1/h)") # Main text Results: rate constant of 0.025 h-1; Additional file 1 v14 K1 = 0.025 h-1
    kfeces <- fixed(0.01)
    label("First-order loss of unabsorbed oral dose (faecal excretion) (1/h)") # Main text Results: first order rate constant of 0.01 h-1; Additional file 1 v37 k1 = 0.01 h-1

    # Elimination ----------------------------------------------------------
    clint_per_g <- fixed(0.018466)
    label("Hepatic intrinsic clearance per gram of liver (1/h per g)") # Additional file 1 v11: k1 = (0.018466*liver weight) h-1
    kcleh <- fixed(5)
    label("Extra-hepatic clearance rate constant acting on kidney E2 (1/h)") # Additional file 1 v10 k1 = 5 h-1

    # Permeability-limited tissues -------------------------------------------
    kperm_in_liver <- fixed(1000)
    label("Liver uptake rate constant, vascular to hepatocyte (1/h)") # Main text Results: uptake 1000 L.h-1; Additional file 1 v32 k1 = 1000
    kperm_out_liver <- fixed(277.8)
    label("Liver efflux rate constant, hepatocyte to vascular (1/h)") # Main text Results: efflux 277.8 L.h-1; Additional file 1 v32 k2 = 277.8
    kperm_in_gonads <- fixed(1000)
    label("Gonad uptake rate constant, vascular to cells (1/h)") # Additional file 1 v19 k1 = 1000; BioModels 'Gonads: Permeability-limited access' k1 = 1000
    kperm_out_gonads <- fixed(100)
    label("Gonad efflux rate constant, cells to vascular (1/h)") # Additional file 1 v19 k2 = 100; BioModels k2 = 100

    # Plasma protein binding ----------------------------------------------------
    kon_shbg <- fixed(3600)
    label("E2-SHBG association rate constant (1/(nM*h))") # Additional file 1 v5 / v48 k1 = 3600 h-1
    koff_shbg <- fixed(5400)
    label("E2-SHBG dissociation rate constant (1/h)") # Additional file 1 v5 / v48 k2 = 5400 h-1 (Kd 1.5 nM)
    kon_albumin <- fixed(3600)
    label("E2-albumin association rate constant (1/(nM*h))") # Additional file 1 v4 / v47 k1 = 3600 h-1
    koff_albumin <- fixed(6.12e7)
    label("E2-albumin dissociation rate constant (1/h)") # Additional file 1 v4 / v47 k2 = 6.12e+7 h-1 (Kd 17 uM)
    bl_shbg <- fixed(40)
    label("Total (free + E2-bound) SHBG concentration in venous and arterial plasma (nM)") # BioModels MODEL1711210001 initial concentrations: free SHBG + E2.SHBG = 40 nM in both blood compartments
    bl_albumin <- fixed(500000)
    label("Total (free + E2-bound) albumin concentration in venous and arterial plasma (nM)") # BioModels MODEL1711210001 initial concentrations: free albumin + E2.albumin = 500000 nM in both blood compartments

    # Estrogen-receptor binding ----------------------------------------------
    kon_er_liver <- fixed(3600)
    label("E2-estrogen receptor association rate constant in liver cells (1/(nM*h))") # Additional file 1 v25 k1 = 3600 h-1
    koff_er_liver <- fixed(900)
    label("E2-estrogen receptor dissociation rate constant in liver cells (1/h)") # Additional file 1 v25 k2 = 900 h-1
    kon_er_gonads <- fixed(3600)
    label("E2-estrogen receptor association rate constant in gonad cells (1/(nM*h))") # BioModels MODEL1711210001 'Gonads: ER binding' k1 = 3600 (Additional file 1 v18 prints 1000; see vignette)
    koff_er_gonads <- fixed(900)
    label("E2-estrogen receptor dissociation rate constant in gonad cells (1/h)") # BioModels MODEL1711210001 'Gonads: ER binding' k2 = 900 (Additional file 1 v18 prints 100; see vignette)
    bl_er_liver <- fixed(12.5)
    label("Total (free + E2-bound) estrogen receptor in liver cells (nM)") # BioModels MODEL1711210001 initial concentrations: ER + E2.ER (liver cells) = 9.565 + 2.935 = 12.5 nM
    bl_er_gonads <- fixed(140.1)
    label("Total (free + E2-bound) estrogen receptor in gonad cells (nM)") # BioModels MODEL1711210001 initial concentrations: ER + E2.ER (gonad cells) = 0.7848 + 139.3152 = 140.1 nM
  })
  model({
    # ---- Anthropometrics (Additional file 1, 'Other Measurements') ----
    ht_m <- HT / 100 # source equations use height in metres
    bsa <- 0.007184 * WT^0.425 * HT^0.725 # DuBois BSA (m^2), height in cm
    bmi <- WT / ht_m^2
    # The deposited equations use the literals 2.718 and 2.71 in place of
    # e; both are needed to reproduce the deposited volumes exactly.
    e_a <- 2.718
    e_b <- 2.71

    # ---- Organ volumes, L (Additional file 1 Global Quantities; Bosgra 2012) ----
    v_adipose <- WT * ((1.2 * bmi + 0.23 * AGE - 5.4) / 100) / 0.916 # Deurenberg; the printed '-10.8*0' sex term is 0 for women
    v_plasma_arterial <- ((2.66 * bsa - 0.46) * 0.33 * 0.55) / 1.025
    v_plasma_venous <- ((2.66 * bsa - 0.46) * 0.67 * 0.55) / 1.025
    v_bone <- e_a^(2.67 * log(ht_m) + 0.0689) / 1.99
    v_brain <- (0.373 * (3.68 - 2.68 * e_b^(-AGE / 0.89)) * (e_b^(-AGE / 629))) / 1.035
    v_gonads <- ((3.3 + 90 * (1 - e_b^(-(AGE / 16.8)^6.7))) / 1000) / 1.03
    v_heart <- e_a^(2.13 * log(ht_m) - 2.502) / 1.03
    v_intestine <- e_a^(2.47 * log(ht_m) - 1.351) / 1.042
    v_kidney <- e_a^(1.93 * log(ht_m) - 2.306) / 1.05
    v_liver <- e_a^(1.98 * log(ht_m) - 0.6786) / 1.03
    v_liver_blood <- v_liver * 0.075 # computed by the source; not used by any rate law
    v_lung <- e_a^(2.1 * log(ht_m) - 2.092) / 1.05
    v_pancreas <- e_a^(2.43 * log(ht_m) - 3.431) / 1.045
    v_other <- e_a^(1.95 * log(ht_m) - 0.072) / 1.03
    v_skin <- e_b^(1.63 * bsa - 1.93) / 1.5 # deposited form; the supplement prints e^(2.63*ln(Height)-1.93)
    v_spleen <- e_a^(2.16 * log(ht_m) - 3.123) / 1.054
    v_stomach <- e_a^(2.45 * log(ht_m) - 3.266) / 1.044 # deposited density 1.044; the supplement prints 1.44
    v_thymus <- ((14 * (7.1 - 6.1 * e_b^(-AGE / 11.9)) * (0.14 + 0.86 * e_b^(-AGE / 10.3))) / 1000) / 1.03
    v_muscle <- (WT * 0.93 - (v_adipose + v_plasma_arterial + v_plasma_venous +
      v_bone + v_brain + v_gonads + v_heart + v_intestine + v_kidney + v_liver +
      v_lung + v_pancreas + v_other + v_skin + v_spleen + v_stomach + v_thymus)) / 1.041

    # ---- Blood flows, L/h: volume * specific flow * density * 0.55 ----
    # (Additional file 1, 'Organ Blood Flow', Q1-Q17)
    q_adipose <- v_adipose * 1.4 * 0.916 * 0.55
    q_bone <- v_bone * 1.8 * 1.99 * 0.55
    q_brain <- v_brain * 32.7 * 1.035 * 0.55
    q_gonads <- v_gonads * 0.8 * 1.03 * 0.55
    q_heart <- v_heart * 70.8 * 1.03 * 0.55
    q_intestine <- v_intestine * 59.0 * 1.042 * 0.55
    q_kidney <- v_kidney * 215 * 1.05 * 0.55
    q_liver_artery <- v_liver * 16.4 * 1.03 * 0.55
    q_muscle <- v_muscle * 2.4 * 1.041 * 0.55
    q_pancreas <- v_pancreas * 29.5 * 1.045 * 0.55
    q_other <- v_other * 18.8 * 1.03 * 0.55
    q_skin <- v_skin * 7.7 * 1.5 * 0.55
    q_spleen <- v_spleen * 81.5 * 1.054 * 0.55
    q_stomach <- v_stomach * 25.3 * 1.044 * 0.55
    q_thymus <- v_thymus * 18.5 * 1.03 * 0.55
    q_liver_vein <- q_liver_artery + q_intestine + q_pancreas + q_spleen + q_stomach # Q9
    # Cardiac output = sum of the 15 arterial-side flows (deposited rule);
    # it is the lung flow (Q10).
    q_lung <- q_adipose + q_bone + q_brain + q_gonads + q_heart + q_intestine +
      q_kidney + q_liver_artery + q_muscle + q_pancreas + q_other + q_skin +
      q_spleen + q_stomach + q_thymus

    clint <- clint_per_g * v_liver * 1000 * 1.03 # v11; liver weight (g) = Liver_cells * 1000 * 1.03
    ka <- exp(lka)

    # ---- Fluxes (nM/h) ----
    # Every flow term is divided by the volume of the SOURCE compartment of
    # the flux (BioModels functions 'PBPK: blood -> tissue' = Qt/Vt * Ca and
    # 'PBPK: tissue -> blood' = Qt/Vt * Ct/(Kp/BP), with Kp = BP = 1 for every
    # tissue), exactly as deposited. This is not volume-consistent mass
    # balance for the blood pools; see the vignette.
    f_lung_in <- q_lung / v_lung * venous # v34
    f_lung_out <- q_lung / v_lung * lung # v33
    f_adipose_in <- q_adipose / v_adipose * arterial
    f_adipose_out <- q_adipose / v_adipose * adipose
    f_bone_in <- q_bone / v_bone * arterial
    f_bone_out <- q_bone / v_bone * bone
    f_brain_in <- q_brain / v_brain * arterial
    f_brain_out <- q_brain / v_brain * brain
    f_heart_in <- q_heart / v_heart * arterial
    f_heart_out <- q_heart / v_heart * heart
    f_muscle_in <- q_muscle / v_muscle * arterial
    f_muscle_out <- q_muscle / v_muscle * muscle
    f_skin_in <- q_skin / v_skin * arterial
    f_skin_out <- q_skin / v_skin * skin
    f_kidney_in <- q_kidney / v_kidney * arterial
    f_kidney_out <- q_kidney / v_kidney * kidney
    f_thymus_in <- q_thymus / v_thymus * arterial
    f_thymus_out <- q_thymus / v_thymus * thymus
    f_other_in <- q_other / v_other * arterial
    f_other_out <- q_other / v_other * other
    f_gonads_in <- q_gonads / v_gonads * arterial
    f_gonads_out <- q_gonads / v_gonads * vp_gonads
    f_intestine_in <- q_intestine / v_intestine * arterial
    f_intestine_liver <- q_intestine / v_intestine * intestine
    f_spleen_in <- q_spleen / v_spleen * arterial
    f_spleen_liver <- q_spleen / v_spleen * spleen
    f_stomach_in <- q_stomach / v_stomach * arterial
    f_stomach_liver <- q_stomach / v_stomach * stomach
    f_pancreas_in <- q_pancreas / v_pancreas * arterial
    f_pancreas_liver <- q_pancreas / v_pancreas * pancreas
    f_liver_artery <- q_liver_artery / v_liver * arterial # v27 (deposit divides by Liver_cells volume)
    f_liver_vein <- q_liver_vein / v_liver * vp_liver # v26 (deposit divides by Liver_cells volume)

    f_perm_liver <- kperm_in_liver * vp_liver - kperm_out_liver * liver # v32
    f_perm_gonads <- kperm_in_gonads * vp_gonads - kperm_out_gonads * gonads # v19
    f_er_liver <- kon_er_liver * target_er_liver * liver - koff_er_liver * complex_er_liver # v25
    f_er_gonads <- kon_er_gonads * target_er_gonads * gonads - koff_er_gonads * complex_er_gonads # v18
    f_shbg_venous <- kon_shbg * venous * shbg_venous - koff_shbg * complex_shbg_venous # v48
    f_albumin_venous <- kon_albumin * venous * albumin_venous - koff_albumin * complex_albumin_venous # v47
    f_shbg_arterial <- kon_shbg * arterial * shbg_arterial - koff_shbg * complex_shbg_arterial # v5
    f_albumin_arterial <- kon_albumin * arterial * albumin_arterial - koff_albumin * complex_albumin_arterial # v4
    f_clint <- clint * liver # v11
    f_cleh <- kcleh * kidney # v10
    f_iv <- kiv * depot_iv # v12
    f_oral <- ka * depot # v14

    # ---- Initial conditions: deposited endogenous steady state ----
    # (BioModels MODEL1711210001 initial concentrations, reference woman
    # 60 kg / 1.63 m / 35 y). The arterial-side tissues share one value
    # because Qt/Vt cancels for each of them at steady state.
    venous(0) <- 0.012880324439466
    complex_shbg_venous(0) <- 0.340551046408471
    shbg_venous(0) <- bl_shbg - 0.340551046408471
    complex_albumin_venous(0) <- 0.378832784791379
    albumin_venous(0) <- bl_albumin - 0.378832784791379
    arterial(0) <- 0.0227598152748276
    complex_shbg_arterial(0) <- 0.597856997340431
    shbg_arterial(0) <- bl_shbg - 0.597856997340431
    complex_albumin_arterial(0) <- 0.669405435577776
    albumin_arterial(0) <- bl_albumin - 0.669405435577776
    lung(0) <- 0.012880324439466
    adipose(0) <- 0.0227598152748276
    bone(0) <- 0.0227598152748276
    brain(0) <- 0.0227598152748276
    heart(0) <- 0.0227598152748276
    muscle(0) <- 0.0227598152748276
    skin(0) <- 0.0227598152748276
    intestine(0) <- 0.0227598152748276
    spleen(0) <- 0.0227598152748276
    stomach(0) <- 0.0227598152748276
    pancreas(0) <- 0.0227598152748276
    thymus(0) <- 0.0227598152748276
    other(0) <- 0.0227598152748276
    kidney(0) <- 0.0218787617463336
    vp_liver(0) <- 0.023199174381171
    liver(0) <- 0.0767047091066076
    complex_er_liver(0) <- 2.93478740007918
    target_er_liver(0) <- bl_er_liver - 2.93478740007918
    vp_gonads(0) <- 4.43582248094885
    gonads(0) <- 44.3782248094892
    complex_er_gonads(0) <- 139.315182
    target_er_gonads(0) <- bl_er_gonads - 139.315182

    # ---- Balance equations (Additional file 1 'Balance Equations') ----
    d/dt(depot) <- -f_oral - kfeces * depot
    d/dt(depot_iv) <- -f_iv
    d/dt(venous) <- f_adipose_out + f_bone_out + f_brain_out + f_heart_out +
      f_muscle_out + f_skin_out + f_liver_vein + f_kidney_out + f_thymus_out +
      f_gonads_out + f_other_out + f_iv - f_lung_in - f_shbg_venous - f_albumin_venous
    d/dt(shbg_venous) <- -f_shbg_venous
    d/dt(complex_shbg_venous) <- f_shbg_venous
    d/dt(albumin_venous) <- -f_albumin_venous
    d/dt(complex_albumin_venous) <- f_albumin_venous
    d/dt(arterial) <- f_lung_out - f_adipose_in - f_bone_in - f_brain_in -
      f_heart_in - f_muscle_in - f_skin_in - f_intestine_in - f_liver_artery -
      f_spleen_in - f_kidney_in - f_thymus_in - f_gonads_in - f_stomach_in -
      f_pancreas_in - f_other_in - f_shbg_arterial - f_albumin_arterial
    d/dt(shbg_arterial) <- -f_shbg_arterial
    d/dt(complex_shbg_arterial) <- f_shbg_arterial
    d/dt(albumin_arterial) <- -f_albumin_arterial
    d/dt(complex_albumin_arterial) <- f_albumin_arterial
    d/dt(lung) <- f_lung_in - f_lung_out
    d/dt(adipose) <- f_adipose_in - f_adipose_out
    d/dt(bone) <- f_bone_in - f_bone_out
    d/dt(brain) <- f_brain_in - f_brain_out
    d/dt(heart) <- f_heart_in - f_heart_out
    d/dt(muscle) <- f_muscle_in - f_muscle_out
    d/dt(skin) <- f_skin_in - f_skin_out
    d/dt(intestine) <- f_intestine_in + f_oral - f_intestine_liver
    d/dt(spleen) <- f_spleen_in - f_spleen_liver
    d/dt(stomach) <- f_stomach_in - f_stomach_liver
    d/dt(pancreas) <- f_pancreas_in - f_pancreas_liver
    d/dt(thymus) <- f_thymus_in - f_thymus_out
    d/dt(kidney) <- f_kidney_in - f_kidney_out - f_cleh
    d/dt(other) <- f_other_in - f_other_out
    d/dt(vp_liver) <- f_liver_artery + f_intestine_liver + f_spleen_liver +
      f_stomach_liver + f_pancreas_liver - f_liver_vein - f_perm_liver
    d/dt(liver) <- f_perm_liver - f_er_liver - f_clint
    d/dt(target_er_liver) <- -f_er_liver
    d/dt(complex_er_liver) <- f_er_liver
    d/dt(vp_gonads) <- f_gonads_in - f_gonads_out - f_perm_gonads
    d/dt(gonads) <- f_perm_gonads - f_er_gonads + ksyn_gonads
    d/dt(target_er_gonads) <- -f_er_gonads
    d/dt(complex_er_gonads) <- f_er_gonads

    # Free (unbound) venous E2, nM: the quantity that reproduces Table 1
    # (total venous E2 is venous + complex_shbg_venous + complex_albumin_venous).
    Cc <- venous
  })
}
