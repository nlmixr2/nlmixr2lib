Chou_2022_flunixin_swine_pbpk <- function() {
  description <- paste(
    "Veterinary (swine). PBPK (whole-body, flow-limited, mrgsolve) for",
    "flunixin and its metabolite 5-hydroxyflunixin in swine, one of the",
    "six drug-species calibrations of the Chou 2022 interactive generic",
    "PBPK (igPBPK) platform used to predict withdrawal intervals in edible",
    "tissues. Plasma, liver, kidney, muscle, fat and a lumped rest of body",
    "are perfused in parallel; only the unbound plasma fraction exchanges",
    "with tissues. Intravenous doses go into `plasma`; oral doses into",
    "`stomach`, which empties into the small intestine, from which drug is",
    "absorbed into the liver, reabsorbed into the liver (enterohepatic",
    "circulation) or lost in faeces. Intramuscular doses go into `depot`",
    "and subcutaneous doses into `depot3`, each absorbed first-order into",
    "plasma. Elimination is first-order hepatic metabolism, biliary",
    "excretion into the small intestine and urinary excretion from the",
    "kidney. A parallel 5-hydroxyflunixin submodel (no fat compartment) is",
    "formed mole-for-mole in the liver; its biliary excretion re-enters",
    "the gut as parent (Supplementary Eq S13). Fixed between-animal",
    "variability reproduces the authors' Monte Carlo population model;",
    "propSd is a placeholder. Amounts in mg, concentrations in mg/L (ug/g",
    "in tissue)."
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
  # and subcutaneous -> depot3.
  dosing <- c("plasma", "stomach", "depot", "depot3")

  compartmentData <- list(
    depot = list(analyte = "flunixin", units = "mg", specimen = "administration site", verified = TRUE),
    depot3 = list(analyte = "flunixin", units = "mg", specimen = "administration site", verified = TRUE),
    stomach = list(analyte = "flunixin", units = "mg", specimen = "administration site", verified = TRUE),
    a_small_intestine = list(analyte = "flunixin", units = "mg", specimen = "administration site", verified = TRUE),
    plasma = list(analyte = "flunixin", units = "mg", specimen = "plasma", verified = TRUE),
    a_liver = list(analyte = "flunixin", units = "mg", specimen = "tissue", verified = TRUE),
    a_kidney = list(analyte = "flunixin", units = "mg", specimen = "tissue", verified = TRUE),
    a_muscle = list(analyte = "flunixin", units = "mg", specimen = "tissue", verified = TRUE),
    a_fat = list(analyte = "flunixin", units = "mg", specimen = "tissue", verified = TRUE),
    a_remainder = list(analyte = "flunixin", units = "mg", specimen = "tissue", verified = TRUE),
    plasma_5oh = list(analyte = "5-hydroxyflunixin", units = "mg", specimen = "plasma", verified = TRUE),
    a_liver_5oh = list(analyte = "5-hydroxyflunixin", units = "mg", specimen = "tissue", verified = TRUE),
    a_kidney_5oh = list(analyte = "5-hydroxyflunixin", units = "mg", specimen = "tissue", verified = TRUE),
    a_muscle_5oh = list(analyte = "5-hydroxyflunixin", units = "mg", specimen = "tissue", verified = TRUE),
    a_remainder_5oh = list(analyte = "5-hydroxyflunixin", units = "mg", specimen = "tissue", verified = TRUE),
    a_urine = list(analyte = "flunixin", units = "mg", specimen = "urine", verified = TRUE),
    a_feces = list(analyte = "flunixin", units = "mg", specimen = "faeces", verified = TRUE),
    a_metabolized = list(analyte = "flunixin", units = "mg", specimen = "not applicable", verified = TRUE),
    a_urine_5oh = list(analyte = "5-hydroxyflunixin", units = "mg", specimen = "urine", verified = TRUE)
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
    n_studies = 7L,
    age_range = NA_character_,
    weight_range = "3.4-168 kg across the calibration studies (code); the Monte Carlo simulations use 92 kg (CV 20%)",
    sex_female_pct = NA_real_,
    disease_state = "healthy food-producing animals; lactating animals excluded",
    dose_range = "IV 2-3 mg/kg, IM 2.2-2.4 mg/kg and oral 2.2-3.3 mg/kg, single or once daily for 3 days",
    regions = "not reported (published studies collected from the FARAD Comparative Pharmacokinetic Database)",
    notes = paste(
      "Seven literature datasets (Chou 2022 Table 1): Howard 2014,",
      "Pairis-Garcia 2013, FDA (2005), Buur 2006, EMA (1999) and Bates 2020",
      "for calibration and Kittrell 2020 for evaluation; IV, IM and oral",
      "2-3.3 mg/kg, single or daily for 3 days; plasma, liver, muscle,",
      "kidney, fat and urine; flunixin and 5-hydroxyflunixin. Mean",
      "concentrations were taken from tables or digitised from figures;",
      "individual animal counts are not reported. Parameters were calibrated",
      "by weighted least squares in FME."
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
    # the precision of the deposited Fit_Swine_FLU.rds (the table rounds
    # or truncates them); literature, in silico and code-only values are
    # fixed. Per-kg rate constants and clearances are scaled by WT in model().
    # ---------------------------------------------------------------
    lka_im <- fixed(log(1))
    label("Absorption rate constant Kim from the intramuscular injection site to plasma (1/h)") # Table 4 Kim 1 (Li 2019a) and the calibration script; the figure and app code use 0.5
    lka_sc <- fixed(log(0.4))
    label("Absorption rate constant Ksc from the subcutaneous injection site to plasma (1/h)") # Table 4 Ksc 0.4 (Li 2019a)
    frac_im <- fixed(1)
    label("Fraction Fracim of an intramuscular dose immediately available at the injection site (unitless)") # Table 4 Fracim blank; code Fracim = 1
    frac_sc <- fixed(1)
    label("Fraction Fracsc of a subcutaneous dose immediately available at the injection site (unitless)") # Table 4 Fracsc blank; code Fracsc = 1
    lka <- fixed(log(0.4))
    label("Intestinal absorption rate constant Kabs, small intestine to liver (1/h)") # Table 4 Kabs 0.4 (Li 2019a) and the calibration script; the figure and app code use 1.9
    lkfec <- fixed(log(0.5))
    label("Faecal elimination rate constant Kunabs from the small intestine (1/h)") # Table 4 Kunabs 0.5 (Li 2019a)
    lkreab <- log(0.00406284)
    label("Enterohepatic reabsorption rate constant per kg KehcC, small intestine to liver (1/h/kg)") # Table 4 KehcC 0.004 (fitted)
    lkmet <- log(0.00425709)
    label("Hepatic metabolic rate constant per kg KmetC, applied to the liver amount (1/h/kg)") # Table 4 KmetC 0.004 (fitted)
    lfu <- log(0.12747)
    label("Unbound fraction of parent drug in plasma fR (unitless)") # Table 4 fR 0.13 (fitted)
    lfu_5oh <- fixed(log(0.01))
    label("Unbound fraction of 5-hydroxyflunixin in plasma fR1 (unitless)") # not in Table 4; code fR1 = 0.01
    lcl_bile <- log(0.401223)
    label("Biliary clearance per kg KbileC of parent, applied to the liver venous concentration (L/h/kg)") # Table 4 KbileC 0.401 (fitted)
    lcl_bile_5oh <- log(0.0338498)
    label("Biliary clearance per kg KbileC1 of 5-hydroxyflunixin, applied to the liver venous concentration (L/h/kg)") # Table 4 KbileC1 0.034 (fitted)
    lcl_renal <- log(0.253339)
    label("Urinary clearance per kg KurineC of parent, applied to the kidney venous concentration (L/h/kg)") # Table 4 KurineC 0.253 (fitted)
    lcl_renal_5oh <- log(0.0071342)
    label("Urinary clearance per kg KurineC1 of 5-hydroxyflunixin, applied to the kidney venous concentration (L/h/kg)") # Table 4 KurineC1 0.007 (fitted)
    lkp_liver <- log(2.30087)
    label("Liver:plasma partition coefficient PL of parent (unitless)") # Table 4 PL 2.30 (fitted)
    lkp_kidney <- log(6.99774)
    label("Kidney:plasma partition coefficient PK of parent (unitless)") # Table 4 PK 6.99 (fitted)
    lkp_muscle <- log(0.314037)
    label("Muscle:plasma partition coefficient PM of parent (unitless)") # Table 4 PM 0.31 (fitted)
    lkp_fat <- fixed(log(0.6))
    label("Fat:plasma partition coefficient PF of parent (unitless)") # Table 4 PF 0.6 (in silico prediction)
    lkp_remainder <- log(20.7167)
    label("Rest-of-body:plasma partition coefficient PRest of parent (unitless)") # Table 4 PRest 20.7 (fitted)
    lkp_liver_5oh <- fixed(log(9.26))
    label("Liver:plasma partition coefficient PL1 of 5-hydroxyflunixin (unitless)") # not in Table 4; code PL1 = 9.26
    lkp_kidney_5oh <- fixed(log(4))
    label("Kidney:plasma partition coefficient PK1 of 5-hydroxyflunixin (unitless)") # not in Table 4; code PK1 = 4
    lkp_muscle_5oh <- fixed(log(1.57))
    label("Muscle:plasma partition coefficient PM1 of 5-hydroxyflunixin (unitless)") # not in Table 4; figure and app code PM1 = 1.57 (calibration used the in silico 1.5747)
    lkp_remainder_5oh <- fixed(log(5))
    label("Rest-of-body:plasma partition coefficient PRest1 of 5-hydroxyflunixin (unitless)") # not in Table 4; code PRest1 = 5
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
    etafrac_im ~ fixed(0.01)
    # CV 10%, normal: SD^2 = (0.1 * 1)^2
    etafrac_sc ~ fixed(0.01)
    # CV 10%, normal: SD^2 = (0.1 * 1)^2
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
    etalfu_5oh ~ fixed(0.00995033)
    # CV 10%, log-normal: log(1 + 0.1^2)
    etalcl_bile ~ fixed(0.0861777)
    # CV 30%, log-normal: log(1 + 0.3^2)
    etalcl_bile_5oh ~ fixed(0.0861777)
    # CV 30%, log-normal: log(1 + 0.3^2)
    etalcl_renal ~ fixed(0.0861777)
    # CV 30%, log-normal: log(1 + 0.3^2)
    etalcl_renal_5oh ~ fixed(0.0861777)
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
    etalkp_liver_5oh ~ fixed(0.0392207)
    # CV 20%, log-normal: log(1 + 0.2^2)
    etalkp_kidney_5oh ~ fixed(0.0392207)
    # CV 20%, log-normal: log(1 + 0.2^2)
    etalkp_muscle_5oh ~ fixed(0.0392207)
    # CV 20%, log-normal: log(1 + 0.2^2)
    etalkp_remainder_5oh ~ fixed(0.0392207)
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
    q_remainder_5oh <- q_co - q_liver - q_kidney - q_muscle

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
    v_remainder_5oh <- v_remainder + v_fat
    hct_e <- hct + etahct
    hct_i <- min(max(hct_e, hct - 0.098), hct + 0.098)
    v_plasma <- v_blood * (1 - hct_i) # code `VPlas = Vblood*(1-Htc)`

    # =================================================================
    # Individual chemical-specific parameters
    # =================================================================
    ka_im <- exp(lka_im + etalka_im)
    ka_sc <- exp(lka_sc + etalka_sc)
    frac_im_e <- frac_im + etafrac_im
    frac_im_i <- min(max(frac_im_e, frac_im - 0.196), frac_im + 0.196)
    frac_sc_e <- frac_sc + etafrac_sc
    frac_sc_i <- min(max(frac_sc_e, frac_sc - 0.196), frac_sc + 0.196)
    ka <- exp(lka + etalka)
    kfec <- exp(lkfec + etalkfec)
    kreab <- exp(lkreab + etalkreab) * WT
    kmet <- exp(lkmet + etalkmet) * WT
    fu <- exp(lfu + etalfu)
    fu_5oh <- exp(lfu_5oh + etalfu_5oh)
    cl_bile <- exp(lcl_bile + etalcl_bile) * WT
    cl_bile_5oh <- exp(lcl_bile_5oh + etalcl_bile_5oh) * WT
    cl_renal <- exp(lcl_renal + etalcl_renal) * WT
    cl_renal_5oh <- exp(lcl_renal_5oh + etalcl_renal_5oh) * WT
    kp_liver <- exp(lkp_liver + etalkp_liver)
    kp_kidney <- exp(lkp_kidney + etalkp_kidney)
    kp_muscle <- exp(lkp_muscle + etalkp_muscle)
    kp_fat <- exp(lkp_fat + etalkp_fat)
    kp_remainder <- exp(lkp_remainder + etalkp_remainder)
    kp_liver_5oh <- exp(lkp_liver_5oh + etalkp_liver_5oh)
    kp_kidney_5oh <- exp(lkp_kidney_5oh + etalkp_kidney_5oh)
    kp_muscle_5oh <- exp(lkp_muscle_5oh + etalkp_muscle_5oh)
    kp_remainder_5oh <- exp(lkp_remainder_5oh + etalkp_remainder_5oh)
    kst <- exp(lkst)

    # Molecular weights (Table 3 (PubChem) 296.4 / 312.4; the code uses 296.24 / 312.24). The code
    # tracks amounts in mmol so that metabolism is mole-for-mole; here
    # amounts are in mg and the molar ratio converts between species.
    mw <- 296.4
    mw_5oh <- 312.4
    mw_ratio <- mw_5oh / mw

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

    c_plasma_5oh <- plasma_5oh / (v_plasma * fu_5oh)
    cv_liver_5oh <- a_liver_5oh / v_liver / kp_liver_5oh
    cv_kidney_5oh <- a_kidney_5oh / v_kidney / kp_kidney_5oh
    cv_muscle_5oh <- a_muscle_5oh / v_muscle / kp_muscle_5oh
    cv_remainder_5oh <- a_remainder_5oh / v_remainder_5oh / kp_remainder_5oh
    c_venous_5oh <- (q_liver * cv_liver_5oh + q_kidney * cv_kidney_5oh +
      q_muscle * cv_muscle_5oh + q_remainder_5oh * cv_remainder_5oh) / q_co

    # =================================================================
    # Fluxes (mg/h)
    # =================================================================
    r_im <- ka_im * depot # code `Rim = Kim*Amtsiteim`
    r_sc <- ka_sc * depot3 # code `Rsc = Ksc*Amtsitesc`
    r_met <- kmet * a_liver # code `Rmet = Kmet*AL` (Supplementary Eq S10)
    r_bile <- cl_bile * cv_liver # code `Rbile = Kbile*CVL`
    r_bile_5oh <- cl_bile_5oh * cv_liver_5oh # code `Rbile1 = Kbile1*CVL1`
    r_urine <- cl_renal * cv_kidney # code `Rurine = Kurine*CVK` (Eq S14)
    r_urine_5oh <- cl_renal_5oh * cv_kidney_5oh # code `Rurine1`
    r_abs <- ka * a_small_intestine # code `RabsSI = Kabs*ASI` (Eq S9)
    r_reab <- kreab * a_small_intestine # code `Rehc = Kehc*ASI` (Eq S11)
    r_fec <- kfec * a_small_intestine # code `Rfeces = Kunabs*ASI`

    # =================================================================
    # ODEs (deposited GenPBPK.R $ODE; Supplementary Equations S1-S16)
    # =================================================================
    # Flunixin has no dissolution step (Table 4 gives no Frac or
    # Kdiss): the whole injected dose sits at the absorption site.
    d/dt(depot) <- -r_im
    d/dt(depot3) <- -r_sc
    # Oral: stomach emptying (code `GE`) into the small intestine (Eqs S7-S9).
    d/dt(stomach) <- -kst * stomach
    # Bile of parent AND metabolite empties into the small intestine; the
    # metabolite is assumed to revert instantly to parent there (Eq S13).
    d/dt(a_small_intestine) <- kst * stomach + r_bile + r_bile_5oh / mw_ratio -
      r_abs - r_reab - r_fec
    d/dt(plasma) <- q_co * fu * (c_venous - c_plasma) + r_im + r_sc # code `RPlas_free`
    d/dt(a_liver) <- q_liver * fu * (c_plasma - cv_liver) - r_met - r_bile + r_abs + r_reab # Eq S12
    d/dt(a_kidney) <- q_kidney * fu * (c_plasma - cv_kidney) - r_urine # Eq S15
    d/dt(a_muscle) <- q_muscle * fu * (c_plasma - cv_muscle)
    d/dt(a_fat) <- q_fat * fu * (c_plasma - cv_fat)
    d/dt(a_remainder) <- q_remainder * fu * (c_plasma - cv_remainder)

    # Metabolite submodel (no fat compartment), formed mole-for-mole in the liver.
    d/dt(plasma_5oh) <- q_co * fu_5oh * (c_venous_5oh - c_plasma_5oh)
    d/dt(a_liver_5oh) <- q_liver * fu_5oh * (c_plasma_5oh - cv_liver_5oh) +
      r_met * mw_ratio - r_bile_5oh
    d/dt(a_kidney_5oh) <- q_kidney * fu_5oh * (c_plasma_5oh - cv_kidney_5oh) - r_urine_5oh
    d/dt(a_muscle_5oh) <- q_muscle * fu_5oh * (c_plasma_5oh - cv_muscle_5oh)
    d/dt(a_remainder_5oh) <- q_remainder_5oh * fu_5oh * (c_plasma_5oh - cv_remainder_5oh)

    # Cumulative elimination (mg), for mass-balance checks.
    d/dt(a_urine) <- r_urine
    d/dt(a_feces) <- r_fec
    d/dt(a_metabolized) <- r_met
    d/dt(a_urine_5oh) <- r_urine_5oh

    # Frac is 1 for flunixin; the authors' Monte Carlo still samples it
    # (normal, CV 10%), which acts as a bioavailability factor here.
    f(depot) <- frac_im_i
    f(depot3) <- frac_sc_i

    # =================================================================
    # Outputs (mg/L = ug/mL; ug/g for tissues at density 1)
    # =================================================================
    Cc <- c_plasma # code `Plasma = CPlas`
    Cliver <- a_liver / v_liver
    Ckidney <- a_kidney / v_kidney
    Cmuscle <- a_muscle / v_muscle
    Cfat <- a_fat / v_fat
    Cc_5oh <- c_plasma_5oh
    Cliver_5oh <- a_liver_5oh / v_liver
    Ckidney_5oh <- a_kidney_5oh / v_kidney
    Cmuscle_5oh <- a_muscle_5oh / v_muscle

    Cc ~ prop(propSd)
  })
}
