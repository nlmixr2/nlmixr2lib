Chou_2022_penicillinG_swine_pbpk <- function() {
  description <- paste(
    "Veterinary (swine). PBPK (whole-body, flow-limited, mrgsolve) for",
    "penicillin g in swine, one of the six drug-species calibrations of",
    "the Chou 2022 interactive generic PBPK (igPBPK) platform used to",
    "predict withdrawal intervals in edible tissues. Plasma, liver,",
    "kidney, muscle, fat and a lumped rest of body are perfused in",
    "parallel; only the unbound plasma fraction exchanges with tissues.",
    "Intravenous doses go into `plasma`; oral doses into `stomach`, which",
    "empties into the small intestine, from which drug is absorbed into",
    "the liver, reabsorbed into the liver (enterohepatic circulation) or",
    "lost in faeces. Intramuscular and subcutaneous injections use a",
    "two-compartment injection site: dose the full amount into both",
    "`depot` and `depot2` (IM) or `depot3` and `depot4` (SC);",
    "bioavailability terms place Frac at the fast-absorption site and 1 -",
    "Frac in a slow-release depot that dissolves into it. Elimination is",
    "first-order hepatic metabolism, biliary excretion into the small",
    "intestine and urinary excretion from the kidney. The paper does not",
    "track penicillin G metabolites, but the deposited code runs its",
    "metabolite submodel with default parameters and returns biliary",
    "metabolite to the gut as parent (Supplementary Eq S13); it is kept",
    "because it roughly doubles late liver residues in cattle. Fixed",
    "between-animal variability reproduces the authors' Monte Carlo",
    "population model; propSd is a placeholder. Amounts in mg,",
    "concentrations in mg/L (ug/g in tissue)."
  )
  reference <- paste(
    "Chou WC, Tell LA, Baynes RE, Davis JL, Maunsell FP, Riviere JE, Lin Z.",
    "An interactive generic physiologically based pharmacokinetic (igPBPK)",
    "modeling platform to predict drug withdrawal intervals in cattle and",
    "swine: a case study on flunixin, florfenicol, and penicillin G.",
    "Toxicol Sci. 2022;188(2):180-197. doi:10.1093/toxsci/kfac056.",
    "Model equations from Supplementary Equations S1-S16 and the authors'",
    "deposited mrgsolve code (https://github.com/UFPBPK/FARAD-igPBPK,",
    "GenPBPK.R); physiology from Table 2; chemical-specific values from",
    "Tables 3 and 4 with fitted values at the precision of the deposited",
    "Fit_*.rds files; Monte Carlo distributions from Materials and Methods.",
    sep = " "
  )
  vignette <- "Chou_2022_igpbpk"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  # Intravenous -> plasma, oral -> stomach, intramuscular -> depot
  # (+ depot2), subcutaneous -> depot3 (+ depot4).
  dosing <- c("plasma", "stomach", "depot", "depot2", "depot3", "depot4")

  compartmentData <- list(
    depot = list(analyte = "penicillin g", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "penicillin g", units = "mg", specimen = "administration site", verified = TRUE),
    depot3 = list(analyte = "penicillin g", units = "mg", specimen = "administration site", verified = TRUE),
    depot4 = list(analyte = "penicillin g", units = "mg", specimen = "administration site", verified = TRUE),
    stomach = list(analyte = "penicillin g", units = "mg", specimen = "administration site", verified = TRUE),
    a_small_intestine = list(analyte = "penicillin g", units = "mg", specimen = "administration site", verified = TRUE),
    plasma = list(analyte = "penicillin g", units = "mg", specimen = "plasma", verified = TRUE),
    a_liver = list(analyte = "penicillin g", units = "mg", specimen = "tissue", verified = TRUE),
    a_kidney = list(analyte = "penicillin g", units = "mg", specimen = "tissue", verified = TRUE),
    a_muscle = list(analyte = "penicillin g", units = "mg", specimen = "tissue", verified = TRUE),
    a_fat = list(analyte = "penicillin g", units = "mg", specimen = "tissue", verified = TRUE),
    a_remainder = list(analyte = "penicillin g", units = "mg", specimen = "tissue", verified = TRUE),
    plasma_metab = list(
      analyte = "penicillin G metabolites (unnamed)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    a_liver_metab = list(
      analyte = "penicillin G metabolites (unnamed)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    a_kidney_metab = list(
      analyte = "penicillin G metabolites (unnamed)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    a_muscle_metab = list(
      analyte = "penicillin G metabolites (unnamed)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    a_remainder_metab = list(
      analyte = "penicillin G metabolites (unnamed)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    a_urine = list(analyte = "penicillin g", units = "mg", specimen = "urine", verified = TRUE),
    a_feces = list(analyte = "penicillin g", units = "mg", specimen = "faeces", verified = TRUE),
    a_metabolized = list(analyte = "penicillin g", units = "mg", specimen = "not applicable", verified = TRUE),
    a_urine_metab = list(
      analyte = "penicillin G metabolites (unnamed)",
      units = "mg",
      specimen = "urine",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      reference_value = 92,
      notes = paste(
        "Every flow and volume is a fraction of body weight and the metabolic,",
        "enterohepatic, biliary and urinary terms are per-kg constants scaled",
        "by BW, exactly as in the code (`QC = QCC*BW`, `Kmet = KmetC*BW`,",
        "...). Because the first-order rate constants KmetC and KehcC are",
        "themselves multiplied by BW, metabolism and reabsorption speed up in",
        "heavier animals. Calibration used each study's own body weight;",
        "reference_value is the body weight of the Monte Carlo",
        "withdrawal-interval simulations (92 kg, sampled normally with CV 20%",
        "in the code; Table 2 does not list BW)."
      ),
      source_name = "BW"
    )
  )

  population <- list(
    species = "swine",
    n_subjects = NA_integer_,
    n_studies = 4L,
    age_range = NA_character_,
    weight_range = "piglets to heavy sows across the calibration studies; the Monte Carlo simulations use 92 kg (CV 20%)",
    sex_female_pct = NA_real_,
    disease_state = "healthy food-producing animals; lactating animals excluded",
    dose_range = "IM 6.5-99 mg/kg and SC 99 mg/kg procaine penicillin G, single or once daily for 3-5 days",
    regions = "not reported (published studies collected from the FARAD Comparative Pharmacokinetic Database)",
    notes = paste(
      "Four literature datasets (Chou 2022 Table 1): Ranheim 2002, Korsrud",
      "1998 and Li 2019b for calibration and Lupton 2014 for evaluation; IM",
      "and SC procaine penicillin G 6.5-99 mg/kg, single or daily for 3-5",
      "days, in piglets, market-age swine and heavy sows; plasma, liver,",
      "kidney, muscle and fat. Mean concentrations were taken from tables or",
      "digitised from figures; individual animal counts are not reported.",
      "Parameters were calibrated by weighted least squares in FME."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Swine physiology (Table 2, mean (SD); Lin 2020b except QFC and QLC,
    # which are from Li 2017). Flows are fractions of
    # cardiac output and volumes fractions of body weight. The rest-of-body
    # fractions are the complements printed in Table 2. All fixed: the
    # model was calibrated with physiology held at these values.
    # ---------------------------------------------------------------
    lqcc <- fixed(log(8.7))
    label("Cardiac output per kg body weight QCC (L/h/kg)") # Table 2 QCC 8.7 (1.62)
    lfq_liver <- fixed(log(0.273))
    label("Fraction of cardiac output to liver QLC (unitless)") # Table 2 QLC 0.273 (0.082)
    lfq_kidney <- fixed(log(0.114))
    label("Fraction of cardiac output to kidney QKC (unitless)") # Table 2 QKC 0.114 (0.032)
    lfq_muscle <- fixed(log(0.342))
    label("Fraction of cardiac output to muscle QMC (unitless)") # Table 2 QMC 0.342 (0.306)
    lfq_fat <- fixed(log(0.128))
    label("Fraction of cardiac output to fat QFC (unitless)") # Table 2 QFC 0.128 (0.038)
    lfq_remainder <- fixed(log(0.143))
    label("Fraction of cardiac output to the rest of body QRestC (unitless)") # Table 2 QRestC = 1-QLC-QKC-QMC-QFC = 0.143
    hct <- fixed(0.412)
    label("Hematocrit Htc (unitless)") # Table 2 Htc 0.412 (0.05)
    fv_blood <- fixed(0.0412)
    label("Blood volume as a fraction of body weight, Table 2 'Plasma' VPC (unitless)") # Table 2 VPC 0.0412 (0.0046)
    fv_liver <- fixed(0.0204)
    label("Liver volume as a fraction of body weight VLC (unitless)") # Table 2 VLC 0.0204 (0.0033)
    fv_kidney <- fixed(0.0037)
    label("Kidney volume as a fraction of body weight VKC (unitless)") # Table 2 VKC 0.0037 (0.0011)
    fv_muscle <- fixed(0.3632)
    label("Muscle volume as a fraction of body weight VMC (unitless)") # Table 2 VMC 0.3632 (0.0266)
    fv_fat <- fixed(0.1544)
    label("Fat volume as a fraction of body weight VFC (unitless)") # Table 2 VFC 0.1544 (0.0265)
    fv_remainder <- fixed(0.4171)
    label("Rest-of-body volume as a fraction of body weight VRestC (unitless)") # Table 2 VRestC = 1-VPC-VLC-VKC-VMC-VFC = 0.4171

    # ---------------------------------------------------------------
    # Chemical-specific parameters (Table 4 and the authors' deposited
    # code). Values the paper marks as fitted are unfixed and carried at
    # the precision of the deposited Fit_Swine_PG.rds (the table rounds
    # or truncates them); literature, in silico and code-only values are
    # fixed. Per-kg rate constants and clearances are scaled by WT in model().
    # ---------------------------------------------------------------
    lka_im <- log(0.10364)
    label("Absorption rate constant Kim from the intramuscular injection site to plasma (1/h)") # Table 4 Kim 0.10 (fitted)
    lka_sc <- fixed(log(0.25))
    label("Absorption rate constant Ksc from the subcutaneous injection site to plasma (1/h)") # Table 4 Ksc 0.25 (Li 2017)
    frac_im <- 0.550286
    label("Fraction Fracim of an intramuscular dose immediately available at the injection site (unitless)") # Table 4 Fracim 0.55 (fitted)
    frac_sc <- fixed(0.5)
    label("Fraction Fracsc of a subcutaneous dose immediately available at the injection site (unitless)") # Table 4 Fracsc 0.5 (Li 2017)
    lkdiss_im <- fixed(log(0.007))
    label("Dissolution rate constant Kdissim from the intramuscular slow-release depot (1/h)") # Table 4 Kdissim 0.007 (Li 2017)
    lkdiss_sc <- fixed(log(0.005))
    label("Dissolution rate constant Kdisssc from the subcutaneous slow-release depot (1/h)") # Table 4 Kdisssc 0.005 (Li 2017)
    lka <- fixed(log(1.91))
    label("Intestinal absorption rate constant Kabs, small intestine to liver (1/h)") # Table 4 Kabs 1.9 (assumed equal to florfenicol, Yang 2019); code 1.91
    lkfec <- log(0.813264)
    label("Faecal elimination rate constant Kunabs from the small intestine (1/h)") # Table 4 Kunabs 0.81 (fitted)
    lkreab <- log(0.0106406)
    label("Enterohepatic reabsorption rate constant per kg KehcC, small intestine to liver (1/h/kg)") # Table 4 KehcC 0.01 (fitted)
    lkmet <- log(0.594519)
    label("Hepatic metabolic rate constant per kg KmetC, applied to the liver amount (1/h/kg)") # Table 4 KmetC 0.59 (fitted)
    lfu <- log(0.956455)
    label("Unbound fraction of parent drug in plasma fR (unitless)") # Table 4 fR 0.96 (fitted)
    lfu_metab <- fixed(log(0.634))
    label("Unbound fraction of penicillin G metabolites (unnamed) in plasma fR1 (unitless)") # not in Table 4; figure and app code fR1 = 0.634
    lcl_bile <- log(0.521208)
    label("Biliary clearance per kg KbileC of parent, applied to the liver venous concentration (L/h/kg)") # Table 4 KbileC 0.52 (fitted)
    lcl_bile_metab <- fixed(log(0.1))
    label("Biliary clearance per kg KbileC1 of penicillin G metabolites (unnamed), applied to the liver venous concentration (L/h/kg)") # not in Table 4; mrgsolve default KbileC1 = 0.1
    lcl_renal <- log(1.54088)
    label("Urinary clearance per kg KurineC of parent, applied to the kidney venous concentration (L/h/kg)") # Table 4 KurineC 1.54 (fitted)
    lcl_renal_metab <- fixed(log(0.1))
    label("Urinary clearance per kg KurineC1 of penicillin G metabolites (unnamed), applied to the kidney venous concentration (L/h/kg)") # not in Table 4; mrgsolve default KurineC1 = 0.1
    lkp_liver <- log(0.07405)
    label("Liver:plasma partition coefficient PL of parent (unitless)") # Table 4 PL 0.07 (fitted)
    lkp_kidney <- log(1.42595)
    label("Kidney:plasma partition coefficient PK of parent (unitless)") # Table 4 PK 1.43 (fitted)
    lkp_muscle <- log(0.0788637)
    label("Muscle:plasma partition coefficient PM of parent (unitless)") # Table 4 PM 0.08 (fitted)
    lkp_fat <- fixed(log(0.249))
    label("Fat:plasma partition coefficient PF of parent (unitless)") # Table 4 PF 0.25 (in silico prediction); code 0.249
    lkp_remainder <- log(0.462613)
    label("Rest-of-body:plasma partition coefficient PRest of parent (unitless)") # Table 4 PRest 0.46 (fitted)
    lkp_liver_metab <- fixed(log(10.52))
    label("Liver:plasma partition coefficient PL1 of penicillin G metabolites (unnamed) (unitless)") # not in Table 4; mrgsolve default PL1 = 10.52
    lkp_kidney_metab <- fixed(log(4))
    label("Kidney:plasma partition coefficient PK1 of penicillin G metabolites (unnamed) (unitless)") # not in Table 4; mrgsolve default PK1 = 4
    lkp_muscle_metab <- fixed(log(0.5))
    label("Muscle:plasma partition coefficient PM1 of penicillin G metabolites (unnamed) (unitless)") # not in Table 4; mrgsolve default PM1 = 0.5
    lkp_remainder_metab <- fixed(log(8))
    label("Rest-of-body:plasma partition coefficient PRest1 of penicillin G metabolites (unnamed) (unitless)") # not in Table 4; mrgsolve default PRest1 = 8
    lkst <- fixed(log(0.182))
    label("Gastric emptying rate constant GE, stomach to small intestine (1/h)") # not in the paper's tables; code `GE = 0.182` (cited there to Yang 2013), never overridden or sampled

    # ---------------------------------------------------------------
    # Between-animal variability of the Monte Carlo population PBPK model
    # (Materials and Methods, 'Establishment of a population PBPK model').
    # These are the authors' assumed coefficients of variation, not
    # estimates, so all are fixed. Log-normal parameters (blood flows,
    # cardiac output and every chemical-specific rate, clearance,
    # partition coefficient and unbound fraction) carry
    # omega^2 = log(1 + CV^2) on the log scale; normally distributed
    # parameters (hematocrit, tissue volume fractions, Frac) enter as
    # value + eta with omega^2 = SD^2. CVs: Table 2 SD / mean for
    # physiology, the default 30% for the rest-of-body flow and volume, 20%
    # for partition coefficients, 30% for rate constants and clearances, 10%
    # for unbound fractions and Frac (the code's default for parameters
    # whose names start with F). Only the parameters the authors' Monte
    # Carlo code samples carry an eta. The authors truncated every draw at
    # its 2.5th / 97.5th percentiles; see the vignette.
    # ---------------------------------------------------------------
    etalqcc ~ fixed(0.0340854)
    # Table 2 SD 1.62 on mean 8.7, log-normal: log(1 + (1.62/8.7)^2)
    etalfq_liver ~ fixed(0.0863794)
    # Table 2 SD 0.082 on mean 0.273, log-normal: log(1 + (0.082/0.273)^2)
    etalfq_kidney ~ fixed(0.0758433)
    # Table 2 SD 0.032 on mean 0.114, log-normal: log(1 + (0.032/0.114)^2)
    etalfq_muscle ~ fixed(0.588094)
    # Table 2 SD 0.306 on mean 0.342, log-normal: log(1 + (0.306/0.342)^2)
    etalfq_fat ~ fixed(0.084465)
    # Table 2 SD 0.038 on mean 0.128, log-normal: log(1 + (0.038/0.128)^2)
    etalfq_remainder ~ fixed(0.0861777)
    # default CV 30% (SD 0.0429), log-normal: log(1 + (0.0429/0.143)^2)
    etahct ~ fixed(0.0025)
    # Table 2 SD 0.05 on mean 0.412, normal: SD^2 = 0.05^2
    etafv_blood ~ fixed(2.116e-05)
    # Table 2 SD 0.0046 on mean 0.0412, normal: SD^2 = 0.0046^2
    etafv_liver ~ fixed(1.089e-05)
    # Table 2 SD 0.0033 on mean 0.0204, normal: SD^2 = 0.0033^2
    etafv_kidney ~ fixed(1.21e-06)
    # Table 2 SD 0.0011 on mean 0.0037, normal: SD^2 = 0.0011^2
    etafv_muscle ~ fixed(0.00070756)
    # Table 2 SD 0.0266 on mean 0.3632, normal: SD^2 = 0.0266^2
    etafv_fat ~ fixed(0.00070225)
    # Table 2 SD 0.0265 on mean 0.1544, normal: SD^2 = 0.0265^2
    etafv_remainder ~ fixed(0.0156575)
    # default CV 30% (SD 0.1251), normal: SD^2 = 0.1251^2
    etalka_im ~ fixed(0.0861777)
    # CV 30%, log-normal: log(1 + 0.3^2)
    etalka_sc ~ fixed(0.0861777)
    # CV 30%, log-normal: log(1 + 0.3^2)
    etafrac_im ~ fixed(0.00302815)
    # CV 10%, normal: SD^2 = (0.1 * 0.550286)^2
    etafrac_sc ~ fixed(0.0025)
    # CV 10%, normal: SD^2 = (0.1 * 0.5)^2
    etalkdiss_im ~ fixed(0.0861777)
    # CV 30%, log-normal: log(1 + 0.3^2)
    etalkdiss_sc ~ fixed(0.0861777)
    # CV 30%, log-normal: log(1 + 0.3^2)
    etalka ~ fixed(0.0861777)
    # CV 30%, log-normal: log(1 + 0.3^2)
    etalkfec ~ fixed(0.0861777)
    # CV 30%, log-normal: log(1 + 0.3^2)
    etalkreab ~ fixed(0.0861777)
    # CV 30%, log-normal: log(1 + 0.3^2)
    etalkmet ~ fixed(0.0861777)
    # CV 30%, log-normal: log(1 + 0.3^2)
    etalfu ~ fixed(0.00995033)
    # CV 10%, log-normal: log(1 + 0.1^2)
    etalfu_metab ~ fixed(0.00995033)
    # CV 10%, log-normal: log(1 + 0.1^2)
    etalcl_bile ~ fixed(0.0861777)
    # CV 30%, log-normal: log(1 + 0.3^2)
    etalcl_renal ~ fixed(0.0861777)
    # CV 30%, log-normal: log(1 + 0.3^2)
    etalkp_liver ~ fixed(0.0392207)
    # CV 20%, log-normal: log(1 + 0.2^2)
    etalkp_kidney ~ fixed(0.0392207)
    # CV 20%, log-normal: log(1 + 0.2^2)
    etalkp_muscle ~ fixed(0.0392207)
    # CV 20%, log-normal: log(1 + 0.2^2)
    etalkp_fat ~ fixed(0.0392207)
    # CV 20%, log-normal: log(1 + 0.2^2)
    etalkp_remainder ~ fixed(0.0392207)
    # CV 20%, log-normal: log(1 + 0.2^2)

    # ---------------------------------------------------------------
    # Chou 2022 calibrated by weighted least squares (FME, Levenberg-
    # Marquardt) and reports no residual-error model. nlmixr2 requires one,
    # so propSd is a fixed placeholder and must NOT be read as an estimate.
    # ---------------------------------------------------------------
    propSd <- fixed(0.1)
    label("Proportional residual error placeholder, plasma (fraction)") # not reported in Chou 2022; placeholder only
  })

  model({
    # =================================================================
    # Physiology. Individual flow fractions and volume fractions are
    # renormalised so that each set sums to one, as in the code's $MAIN
    # block (sumQ, sumV); at the typical values the sums are already one.
    # The normally distributed etas are clamped at +/-1.96 SD, the
    # 2.5th / 97.5th percentile bounds the authors sampled within, which
    # also keeps the volume fractions positive.
    # =================================================================
    q_co <- exp(lqcc + etalqcc) * WT # code `QC = QCC*BW`
    fq_liver_i <- exp(lfq_liver + etalfq_liver)
    fq_kidney_i <- exp(lfq_kidney + etalfq_kidney)
    fq_muscle_i <- exp(lfq_muscle + etalfq_muscle)
    fq_fat_i <- exp(lfq_fat + etalfq_fat)
    fq_remainder_i <- exp(lfq_remainder + etalfq_remainder)
    fq_sum <- fq_liver_i + fq_kidney_i + fq_muscle_i + fq_fat_i + fq_remainder_i
    q_liver <- q_co * fq_liver_i / fq_sum
    q_kidney <- q_co * fq_kidney_i / fq_sum
    q_muscle <- q_co * fq_muscle_i / fq_sum
    q_fat <- q_co * fq_fat_i / fq_sum
    q_remainder <- q_co - q_liver - q_kidney - q_muscle - q_fat # code `QRestC = 1-QLC-QKC-QMC-QFC`
    # The metabolite submodel has no fat compartment; its rest of body
    # takes the fat flow and volume too (code `QRest1C`, `VRest1C`).
    q_remainder_metab <- q_co - q_liver - q_kidney - q_muscle

    fv_blood_e <- fv_blood + etafv_blood
    fv_blood_i <- min(max(fv_blood_e, fv_blood - 0.009016), fv_blood + 0.009016)
    fv_liver_e <- fv_liver + etafv_liver
    fv_liver_i <- min(max(fv_liver_e, fv_liver - 0.006468), fv_liver + 0.006468)
    fv_kidney_e <- fv_kidney + etafv_kidney
    fv_kidney_i <- min(max(fv_kidney_e, fv_kidney - 0.002156), fv_kidney + 0.002156)
    fv_muscle_e <- fv_muscle + etafv_muscle
    fv_muscle_i <- min(max(fv_muscle_e, fv_muscle - 0.05214), fv_muscle + 0.05214)
    fv_fat_e <- fv_fat + etafv_fat
    fv_fat_i <- min(max(fv_fat_e, fv_fat - 0.05194), fv_fat + 0.05194)
    fv_remainder_e <- fv_remainder + etafv_remainder
    fv_remainder_i <- min(max(fv_remainder_e, fv_remainder - 0.2453), fv_remainder + 0.2453)
    fv_sum <- fv_blood_i + fv_liver_i + fv_kidney_i + fv_muscle_i + fv_fat_i + fv_remainder_i
    v_blood <- WT * fv_blood_i / fv_sum
    v_liver <- WT * fv_liver_i / fv_sum
    v_kidney <- WT * fv_kidney_i / fv_sum
    v_muscle <- WT * fv_muscle_i / fv_sum
    v_fat <- WT * fv_fat_i / fv_sum
    v_remainder <- WT - v_blood - v_liver - v_kidney - v_muscle - v_fat
    v_remainder_metab <- v_remainder + v_fat
    hct_e <- hct + etahct
    hct_i <- min(max(hct_e, hct - 0.098), hct + 0.098)
    v_plasma <- v_blood * (1 - hct_i) # code `VPlas = Vblood*(1-Htc)`

    # =================================================================
    # Individual chemical-specific parameters
    # =================================================================
    ka_im <- exp(lka_im + etalka_im)
    ka_sc <- exp(lka_sc + etalka_sc)
    frac_im_e <- frac_im + etafrac_im
    frac_im_i <- min(max(frac_im_e, frac_im - 0.1079), frac_im + 0.1079)
    frac_sc_e <- frac_sc + etafrac_sc
    frac_sc_i <- min(max(frac_sc_e, frac_sc - 0.098), frac_sc + 0.098)
    kdiss_im <- exp(lkdiss_im + etalkdiss_im)
    kdiss_sc <- exp(lkdiss_sc + etalkdiss_sc)
    ka <- exp(lka + etalka)
    kfec <- exp(lkfec + etalkfec)
    kreab <- exp(lkreab + etalkreab) * WT
    kmet <- exp(lkmet + etalkmet) * WT
    fu <- exp(lfu + etalfu)
    fu_metab <- exp(lfu_metab + etalfu_metab)
    cl_bile <- exp(lcl_bile + etalcl_bile) * WT
    cl_bile_metab <- exp(lcl_bile_metab) * WT
    cl_renal <- exp(lcl_renal + etalcl_renal) * WT
    cl_renal_metab <- exp(lcl_renal_metab) * WT
    kp_liver <- exp(lkp_liver + etalkp_liver)
    kp_kidney <- exp(lkp_kidney + etalkp_kidney)
    kp_muscle <- exp(lkp_muscle + etalkp_muscle)
    kp_fat <- exp(lkp_fat + etalkp_fat)
    kp_remainder <- exp(lkp_remainder + etalkp_remainder)
    kp_liver_metab <- exp(lkp_liver_metab)
    kp_kidney_metab <- exp(lkp_kidney_metab)
    kp_muscle_metab <- exp(lkp_muscle_metab)
    kp_remainder_metab <- exp(lkp_remainder_metab)
    kst <- exp(lkst)

    # Molecular weights (Table 3 (PubChem) 334; the code uses 334.4 for parent and metabolite). The code
    # tracks amounts in mmol so that metabolism is mole-for-mole; here
    # amounts are in mg and the molar ratio converts between species.
    mw <- 334
    mw_metab <- 334
    mw_ratio <- mw_metab / mw

    # =================================================================
    # Concentrations (mg/L). `plasma` is the code's `APlas_free`: the IV
    # dose and the injection-site absorption enter it, and it carries the
    # whole plasma mass balance, but the code reports the TOTAL plasma
    # concentration as `CPlas = CPlas_free/Free`, i.e. plasma / (V * fu).
    # Tissue inflow is `Q*CPlas*Free` and outflow `Q*CV*Free`, where
    # CV = C_tissue / P is the venous-equilibrium concentration.
    # =================================================================
    c_plasma <- plasma / (v_plasma * fu)
    cv_liver <- a_liver / v_liver / kp_liver
    cv_kidney <- a_kidney / v_kidney / kp_kidney
    cv_muscle <- a_muscle / v_muscle / kp_muscle
    cv_fat <- a_fat / v_fat / kp_fat
    cv_remainder <- a_remainder / v_remainder / kp_remainder
    c_venous <- (q_liver * cv_liver + q_kidney * cv_kidney + q_fat * cv_fat +
      q_muscle * cv_muscle + q_remainder * cv_remainder) / q_co

    c_plasma_metab <- plasma_metab / (v_plasma * fu_metab)
    cv_liver_metab <- a_liver_metab / v_liver / kp_liver_metab
    cv_kidney_metab <- a_kidney_metab / v_kidney / kp_kidney_metab
    cv_muscle_metab <- a_muscle_metab / v_muscle / kp_muscle_metab
    cv_remainder_metab <- a_remainder_metab / v_remainder_metab / kp_remainder_metab
    c_venous_metab <- (q_liver * cv_liver_metab + q_kidney * cv_kidney_metab +
      q_muscle * cv_muscle_metab + q_remainder_metab * cv_remainder_metab) / q_co

    # =================================================================
    # Fluxes (mg/h)
    # =================================================================
    r_im <- ka_im * depot # code `Rim = Kim*Amtsiteim`
    r_sc <- ka_sc * depot3 # code `Rsc = Ksc*Amtsitesc`
    r_met <- kmet * a_liver # code `Rmet = Kmet*AL` (Supplementary Eq S10)
    r_bile <- cl_bile * cv_liver # code `Rbile = Kbile*CVL`
    r_bile_metab <- cl_bile_metab * cv_liver_metab # code `Rbile1 = Kbile1*CVL1`
    r_urine <- cl_renal * cv_kidney # code `Rurine = Kurine*CVK` (Eq S14)
    r_urine_metab <- cl_renal_metab * cv_kidney_metab # code `Rurine1`
    r_abs <- ka * a_small_intestine # code `RabsSI = Kabs*ASI` (Eq S9)
    r_reab <- kreab * a_small_intestine # code `Rehc = Kehc*ASI` (Eq S11)
    r_fec <- kfec * a_small_intestine # code `Rfeces = Kunabs*ASI`

    # =================================================================
    # ODEs (deposited GenPBPK.R $ODE; Supplementary Equations S1-S16)
    # =================================================================
    # Two-compartment injection sites (Eqs S1-S6): `depot` / `depot3` are
    # the fast-absorption sites and `depot2` / `depot4` the slow-release
    # depots that dissolve into them.
    d/dt(depot) <- kdiss_im * depot2 - r_im
    d/dt(depot2) <- -kdiss_im * depot2
    d/dt(depot3) <- kdiss_sc * depot4 - r_sc
    d/dt(depot4) <- -kdiss_sc * depot4
    # Oral: stomach emptying (code `GE`) into the small intestine (Eqs S7-S9).
    d/dt(stomach) <- -kst * stomach
    # Bile of parent AND metabolite empties into the small intestine; the
    # metabolite is assumed to revert instantly to parent there (Eq S13).
    d/dt(a_small_intestine) <- kst * stomach + r_bile + r_bile_metab / mw_ratio -
      r_abs - r_reab - r_fec
    d/dt(plasma) <- q_co * fu * (c_venous - c_plasma) + r_im + r_sc # code `RPlas_free`
    d/dt(a_liver) <- q_liver * fu * (c_plasma - cv_liver) - r_met - r_bile + r_abs + r_reab # Eq S12
    d/dt(a_kidney) <- q_kidney * fu * (c_plasma - cv_kidney) - r_urine # Eq S15
    d/dt(a_muscle) <- q_muscle * fu * (c_plasma - cv_muscle)
    d/dt(a_fat) <- q_fat * fu * (c_plasma - cv_fat)
    d/dt(a_remainder) <- q_remainder * fu * (c_plasma - cv_remainder)

    # Metabolite submodel (no fat compartment), formed mole-for-mole in the liver.
    d/dt(plasma_metab) <- q_co * fu_metab * (c_venous_metab - c_plasma_metab)
    d/dt(a_liver_metab) <- q_liver * fu_metab * (c_plasma_metab - cv_liver_metab) +
      r_met * mw_ratio - r_bile_metab
    d/dt(a_kidney_metab) <- q_kidney * fu_metab * (c_plasma_metab - cv_kidney_metab) - r_urine_metab
    d/dt(a_muscle_metab) <- q_muscle * fu_metab * (c_plasma_metab - cv_muscle_metab)
    d/dt(a_remainder_metab) <- q_remainder_metab * fu_metab * (c_plasma_metab - cv_remainder_metab)

    # Cumulative elimination (mg), for mass-balance checks.
    d/dt(a_urine) <- r_urine
    d/dt(a_feces) <- r_fec
    d/dt(a_metabolized) <- r_met
    d/dt(a_urine_metab) <- r_urine_metab

    # Split each injection: dose the full amount into BOTH the fast site
    # and its slow depot; these fractions apportion it (code `DOSEfast =
    # Dose*Frac`, `DOSEslow = Dose*(1-Frac)`).
    f(depot) <- frac_im_i
    f(depot2) <- 1 - frac_im_i
    f(depot3) <- frac_sc_i
    f(depot4) <- 1 - frac_sc_i

    # =================================================================
    # Outputs (mg/L = ug/mL; ug/g for tissues at density 1)
    # =================================================================
    Cc <- c_plasma # code `Plasma = CPlas`
    Cliver <- a_liver / v_liver
    Ckidney <- a_kidney / v_kidney
    Cmuscle <- a_muscle / v_muscle
    Cfat <- a_fat / v_fat

    Cc ~ prop(propSd)
  })
}
