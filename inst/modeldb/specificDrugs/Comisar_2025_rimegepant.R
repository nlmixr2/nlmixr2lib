Comisar_2025_rimegepant <- function() {
  description <- "Two-compartment oral population PK model for the small-molecule CGRP receptor antagonist rimegepant with a six-compartment pre-systemic absorption chain (dosing/transit-0 plus four transit compartments draining at ktr, then a terminal absorption compartment draining at ka), empirical allometric body-weight scaling (exponent 0.75 on CL/F and Q/F, 1.0 on Vc/F and Vp/F, reference 70 kg), moderate and severe hepatic impairment plus fluconazole and itraconazole coadministration on CL/F, fed status and a dose power term on relative bioavailability, and fed status, itraconazole, capsule and orally-disintegrating-tablet formulation, and low (10/25 mg) dose level on ktr. Pooled from 11 phase 1 studies in healthy adults, elderly adults with stable chronic illness, and adults with renal or hepatic impairment (n = 349)."
  reference   <- paste(
    "Comisar CM, Hughes JH, Francis J, Chinda Y, Sano Y, Muto C, Neumar C,",
    "Bhardwaj R, Bertz R, Liu J (2025).",
    "Population pharmacokinetic modeling of the oral calcitonin gene-related",
    "peptide receptor antagonist rimegepant in adults.",
    "CPT: Pharmacometrics & Systems Pharmacology 14(8):1332-1342.",
    "doi:10.1002/psp4.70051.",
    "Structure and covariate equation forms taken from Supplementary Listing 1",
    "(the final-model NONMEM control stream) and the Supplemental Methods",
    "covariate-model equations; all parameter VALUES are the Stage 3 final",
    "estimates in Table 3 of the main paper (the control stream lists initial",
    "estimates only)."
  )
  vignette <- "Comisar_2025_rimegepant"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Compartment roles verified against the $MODEL block of
  # the final-model control stream (Supplementary Listing 1), which names
  # COMP(1) "DOSING AND TRANSIT 0", COMP(4)-COMP(7) "TRANSIT 1"-"TRANSIT 4",
  # COMP(8) "ABSORPTION", COMP(2) "CENTRAL", and COMP(3) "PERIPHERAL".
  compartmentData <- list(
    depot       = list(analyte = "rimegepant", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "rimegepant", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "rimegepant", units = "mg", specimen = "administration site", verified = TRUE),
    transit3    = list(analyte = "rimegepant", units = "mg", specimen = "administration site", verified = TRUE),
    transit4    = list(analyte = "rimegepant", units = "mg", specimen = "administration site", verified = TRUE),
    transit5    = list(analyte = "rimegepant", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "rimegepant", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "rimegepant", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight; drives empirical allometric scaling.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference 70 kg (control-stream terms `(BWGT/70)**THETA(22)` and",
        "`(BWGT/70)**THETA(23)`; Comisar 2025 Discussion paragraph 3 refers to",
        "'the typical clearance for a 70 kg individual estimated by the model').",
        "Exponents are FIXED, not estimated: 0.75 on CL/F and Q/F, 1.0 on Vc/F",
        "and Vp/F (Table 3 note; Methods 2.2 paragraph 3 'empirical allometric",
        "scaling with fixed allometric coefficients of 0.75 for clearances and",
        "1 for volumes of distribution'). Baseline (time-fixed) weight, not",
        "time-varying. Analysis-population median 75 kg, range 46-134 kg."
      ),
      source_name        = "BWGT"
    ),
    HEPIMP_MOD = list(
      description        = "Moderate hepatic impairment indicator (1 = moderate, 0 = normal or mild).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function or mild impairment; mild was screened and found not statistically significant, so mild subjects sit in the reference group).",
      notes              = paste(
        "Multiplicative fractional effect on CL/F: CL/F * (1 + (-0.229)), a",
        "22.9% decrease (Comisar 2025 Table 3 'Moderate hepatic impairment on",
        "CL/F' = -0.229, %RSE 27.8). Corresponds to control-stream",
        "`IF(HEPATC.EQ.2) CLHEPATC = (1 + THETA(24))`. Mutually exclusive with",
        "HEPIMP_SEV. The paper does not state whether the mild/moderate/severe",
        "classification follows Child-Pugh or NCI ODWG criteria; the source",
        "hepatic-impairment study BHV3000-107 enrolled parallel normal / mild /",
        "moderate / severe groups (n = 18/6/6/6). Simulated median AUCtau rises",
        "1.3-fold vs normal hepatic function (Table 4)."
      ),
      source_name        = "HEPATC (value 2)"
    ),
    HEPIMP_SEV = list(
      description        = "Severe hepatic impairment indicator (1 = severe, 0 = normal, mild, or moderate).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function or mild impairment; moderate impairment is carried by the separate HEPIMP_MOD indicator).",
      notes              = paste(
        "Multiplicative fractional effect on CL/F: CL/F * (1 + (-0.423)), a",
        "42.3% decrease (Comisar 2025 Table 3 'Severe hepatic impairment on",
        "CL/F' = -0.423, %RSE 24.9). Corresponds to control-stream",
        "`IF(HEPATC.EQ.3) CLHEPATC = (1 + THETA(11))`. Mutually exclusive with",
        "HEPIMP_MOD. Classification scheme not stated in the paper (see",
        "HEPIMP_MOD notes). Simulated median AUCtau rises 1.7-fold vs normal",
        "hepatic function (Table 4); rimegepant is to be avoided in severe",
        "hepatic impairment (Discussion paragraph 4)."
      ),
      source_name        = "HEPATC (value 3)"
    ),
    CONMED_FLUCONAZOLE = list(
      description        = "Concomitant fluconazole (moderate CYP3A4 / strong CYP2C9 inhibitor) coadministration indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant fluconazole).",
      notes              = paste(
        "Multiplicative fractional effect on CL/F: CL/F * (1 + (-0.427)), a",
        "42.7% decrease (Comisar 2025 Table 3 'Fluconazole use on CL/F' =",
        "-0.427, %RSE 2.88). Corresponds to control-stream",
        "`IF(FLCZ.EQ.1) CLFLCZ = (1 + THETA(8))`. Source study BHV3000-105 gave",
        "fluconazole 400 mg daily on days 5-12 in a fixed-sequence design, so",
        "the flag is time-varying within subject. Fluconazole was also screened",
        "on ktr (control-stream KAFLCZ / THETA(21)) but that effect was fixed to",
        "0 in the final model. Simulated median AUCtau rises 1.8-fold (Table 4);",
        "no dose adjustment is recommended for moderate CYP3A4 inhibitors."
      ),
      source_name        = "FLCZ"
    ),
    CONMED_ITRACONAZOLE = list(
      description        = "Concomitant itraconazole (strong CYP3A4 / P-gp inhibitor) coadministration indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant itraconazole).",
      notes              = paste(
        "Two multiplicative fractional effects. On CL/F: CL/F * (1 + (-0.743)),",
        "a 74.3% decrease (Table 3 'Itraconazole use on CL/F' = -0.743, %RSE",
        "1.25; control-stream `IF(ITCZ.EQ.1) CLITCZ = (1 + THETA(12))`). On ktr:",
        "ktr * (1 + (-0.351)), a 35.1% slowing (Table 3 'Itraconazole use on",
        "ktr' = -0.351, %RSE 28.9; control-stream `IF(ITCZ.EQ.1) KAITCZ = (1 +",
        "THETA(19))`, which multiplies KTR rather than KA despite the KA-prefixed",
        "variable name). Source study BHV3000-103 gave itraconazole 200 mg oral",
        "solution daily on days 5-11 in a fixed-sequence design, so the flag is",
        "time-varying within subject. Simulated median AUCtau rises 3.9-fold",
        "(Table 4); coadministration is not recommended."
      ),
      source_name        = "ITCZ"
    ),
    FED = list(
      description        = "Dose-record fed-vs-fasted indicator (1 = fed, 0 = fasted).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted; all studies dosed fasting except the food-effect study BHV3000-112).",
      notes              = paste(
        "Two multiplicative fractional effects. On relative bioavailability:",
        "F1 * (1 + (-0.315)), a 31.5% decrease (Table 3 'Fed-status on relative",
        "bioavailability (F1)' = -0.315, %RSE 9.74; control-stream",
        "`IF(FED.EQ.1) F1BIOFED = (1 + THETA(13))`). On ktr: ktr * (1 +",
        "(-0.706)), a 70.6% slowing (Table 3 'Fed-status on ktr' = -0.706, %RSE",
        "3.34; control-stream `IF(FED.EQ.1) KAFED = (1 + THETA(16))`). The only",
        "fed data come from BHV3000-112, where the fed arm received a HIGH-FAT",
        "meal, so the effect size is calibrated on a high-fat challenge even",
        "though the paper labels the covariate generically as 'Food effect",
        "(fasting/fed)' (Table 2) and applies it as a general fed flag. FED is",
        "used rather than FED_HIGHFAT to match that general labelling; a",
        "reviewer wanting the high-fat refinement can set FED_HIGHFAT = FED.",
        "Crossover design, so a subject contributes both FED levels."
      ),
      source_name        = "FED"
    ),
    FORM_ODT = list(
      description        = "Orally disintegrating tablet formulation indicator (1 = ODT, 0 = the reference immediate-release tablet).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (immediate-release tablet, the most common formulation in the analysis dataset at 73.1% of subjects; control-stream `IF(FORM.EQ.1) KAFORM = 1 ; Most common`).",
      notes              = paste(
        "Multiplicative fractional effect on ktr only: ktr * (1 + 0.355), a",
        "35.5% faster transit (Table 3 'Oral disintegrating tablet on ktr' =",
        "0.355, %RSE 29.6; control-stream `IF(FORM.EQ.3) KAFORM = (1 +",
        "THETA(17))`). Formulation was also screened on relative bioavailability",
        "(control-stream F1BIOFORM / THETA(14)-(15)) and on CL/F (CLFORM /",
        "THETA(9)-(10)) but both were fixed to 0 in the final model, so ODT and",
        "tablet are bioequivalent in exposure and differ only in absorption rate.",
        "Mutually exclusive with FORM_CAPSULE. The ODT is the marketed Nurtec /",
        "Vydura presentation and is administered sublingually or on the tongue."
      ),
      source_name        = "FORM (value 3)"
    ),
    FORM_CAPSULE = list(
      description        = "Capsule formulation indicator (1 = capsule, 0 = the reference immediate-release tablet).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (immediate-release tablet; the comparator for this model is the tablet, as in Lacy 2018 cabozantinib and Gupta 2016 lenvatinib, not the oral solution used by Hennig 2006 itraconazole).",
      notes              = paste(
        "Multiplicative fractional effect on ktr only: ktr * (1 + 2.03), i.e.",
        "3.03-fold faster transit (Table 3 'Capsule formulation on ktr' = 2.03,",
        "%RSE 24.1; control-stream `IF(FORM.EQ.2) KAFORM = (1 + THETA(18))`).",
        "Formulation on relative bioavailability and on CL/F was screened but",
        "fixed to 0 in the final model, so the capsule is bioequivalent to the",
        "tablet in exposure. Mutually exclusive with FORM_ODT. Only 17 subjects",
        "(4.9%) received the capsule (3 x 25 mg capsules in BHV3000-102), which",
        "is consistent with the wide confidence interval on this effect; the",
        "paper judges the faster ktr not clinically relevant because ka (3.86",
        "1/h) is the slower of the two absorption steps (Discussion paragraph 8)."
      ),
      source_name        = "FORM (value 2)"
    ),
    DOSE = list(
      description        = "Administered rimegepant dose level per dose record, in mg.",
      units              = "mg",
      type               = "continuous",
      reference_category = "10 mg is the normalising reference for the relative-bioavailability power term ((DOSE/10)^0.191); 75 mg is the reference level for the categorical low-dose effect on ktr (control-stream `IF(DOSE.EQ.75) KTRDOSE = 1 ; Most common`).",
      notes              = paste(
        "Use case (a) of the DOSE canonical (per-record assigned dose level",
        "driving a dose-dependent PK parameter), with BOTH a continuous and a",
        "categorical use in the same model.",
        "(1) Continuous power term on relative bioavailability:",
        "`F1BIODOSE = (DOSE/10)**THETA(25)` with THETA(25) = 0.191, so F1 = 1 at",
        "10 mg and 1.47 at 75 mg. This is the paper's encoding of greater-than-",
        "dose-proportional exposure over 10-150 mg, attributed to dose-dependent",
        "autoinhibition of CYP3A-mediated first-pass metabolism (Discussion",
        "paragraph 7). Note that F1 is a RELATIVE bioavailability anchored at 1",
        "for a 10 mg fasted dose, so CL/F and the volumes in Table 3 are the",
        "apparent values on that anchor; absolute fasted oral bioavailability is",
        "about 64% (Discussion paragraph 7).",
        "(2) Categorical low-dose indicator on ktr: the control stream sets",
        "`KTRDOSE = (1 + THETA(29))` for DOSE = 10 and DOSE = 25 and",
        "`KTRDOSE = 1` for DOSE = 75, with the 150 mg level's THETA(28) fixed to",
        "0 (i.e. also 1). The model file reproduces this as `(DOSE < 50)`; 50 mg",
        "is the midpoint of the 25-75 mg gap and is exact over the four dose",
        "levels actually studied (10, 25, 75, 150 mg). Doses outside those four",
        "levels are an extrapolation - the control stream leaves KTRDOSE",
        "undefined there.",
        "Studied levels and subject counts (Table 2, Stage 3): 10 mg n = 41,",
        "25 mg n = 54, 75 mg n = 322, 150 mg n = 14."
      ),
      source_name        = "DOSE"
    )
  )

  # Covariates the paper screened but did NOT retain in the final model, either
  # because they were not statistically significant or because the effect was
  # not clinically relevant. Documented for provenance; not referenced in
  # model(). Full screening matrix in Comisar 2025 Table S1; the retained set
  # after forward inclusion / backward elimination is Table S2.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on volume and clearance parameters (Table S1); not statistically significant after allometric body-weight scaling was included a priori (Results 3.3.3 paragraph 2). Population median 43 years, range 18-77; study BHV3000-108 specifically compared elderly vs nonelderly adults."
    ),
    SEXF = list(
      description = "Female-sex indicator (1 = female, 0 = male).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on volume and clearance parameters (Table S1). Sex on Vc/F did reach the forward-inclusion threshold in Stage 1 (Table S2, dOFV = -22.9) but did not survive backward elimination; the control stream retains it only as the commented-out `; antiquated sex testing on V2` term with THETA(20) fixed to 0. Population 73.6% male."
    ),
    HEPIMP_MILD = list(
      description = "Mild hepatic impairment indicator (1 = mild, 0 = otherwise).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as one level of the categorical hepatic-function covariate (Table S1); not statistically significant (Results 3.3.3 paragraph 2), so mild-impairment subjects fall in the reference group of HEPIMP_MOD / HEPIMP_SEV. n = 6 (1.7%). No dose adjustment is required for mild or moderate hepatic impairment (Discussion paragraph 4)."
    ),
    CRCL = list(
      description = "Renal function as eGFR by the 4-variable MDRD equation.",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened on clearance parameters (Table S1; Methods 2.2 paragraph 3, MDRD 4-variable eGFR). Not statistically significant; mild, moderate, and severe renal impairment (n = 6 each) require no dose adjustment (Conclusions). Population median 92 mL/min, range 9-150. Reported by the paper in mL/min rather than the BSA-normalised mL/min/1.73 m^2 of the CRCL canonical."
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on volume and clearance parameters (Table S1); race was not a statistically significant covariate (Results 3.3.3, Conclusions). n = 31 (8.9%)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (pooling the Japanese and Chinese cohorts of studies 111 and 118 with nonspecified Asian subjects).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in Stage 3 as the overall Asian vs non-Asian effect on rimegepant clearance and volume of distribution; not statistically significant and not retained (Results 3.3.3 paragraph 2; pcVPC by Asian race in Figure S10)."
    ),
    RACE_JAPANESE = list(
      description = "Japanese-ethnicity indicator (study BHV3000-111).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in Stage 2 on clearance and volume parameters; CL/F and Vd/F were not significantly different between Japanese and non-Japanese groups (Results 3.3.2; pcVPC in Figure S8). n = 18 (5.2%)."
    ),
    RACE_CHINESE = list(
      description = "Chinese-ethnicity indicator (study BHV3000-118).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in Stage 3 on clearance and volume parameters; not statistically significant and not retained (Results 3.3.3 paragraph 2; pcVPC in Figure S9). n = 12 (3.4%)."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 349L,
    n_studies        = 11L,
    age_range        = "18-77 years",
    age_median       = "43 years",
    weight_range     = "46-134 kg",
    weight_median    = "75 kg",
    sex_female_pct   = 26.4,
    race_ethnicity   = c(
      Caucasian               = 81.9,
      `Black/African American` = 8.9,
      `Japanese ethnicity`     = 5.2,
      `Chinese ethnicity`      = 3.4,
      `Nonspecified Asian`     = 0.6
    ),
    disease_state    = "Predominantly healthy adults, plus elderly adults with stable chronic illness(es) and adults with renal or hepatic impairment. No migraine patients were included.",
    hepatic_function = "94.8% normal (or not listed); 1.7% each mild, moderate, and severe impairment (n = 6 each, all from BHV3000-107).",
    renal_function   = "94.8% normal (or not listed); 1.7% each mild, moderate, and severe impairment (n = 6 each, all from BHV3000-106). Median eGFR 92 mL/min (range 9-150).",
    co_medication    = "Fixed-sequence DDI arms: itraconazole 200 mg daily (n = 22, 6.3%) and fluconazole 400 mg daily (n = 23, 6.6%).",
    dose_range       = "Single and multiple oral doses of 10, 25, 75, and 150 mg rimegepant as immediate-release tablet, capsule, or orally disintegrating tablet.",
    regions          = "United States and Japan (phase 1 studies; Caucasian, Black/African American, Japanese, and Chinese cohorts).",
    notes            = paste(
      "Stage 3 (final) analysis population, Comisar 2025 Table 2 rightmost",
      "column and Table 1 (per-study designs). 12,587 concentration records were",
      "considered and 1842 (14.6%) excluded as below the limit of quantitation,",
      "leaving 10,745 observations from 349 subjects (Results 3.1 and 3.3.3).",
      "LLOQ was 0.5 ng/mL for studies 106, 107, 108, 110, 111, 112, 117, 118 and",
      "10 ng/mL for studies 102, 103, 105. Percentages for formulation, dose,",
      "fed status, and DDI in Table 2 sum to more than 100% because the",
      "crossover designs let one subject contribute several categorical states.",
      "The model was developed in three sequential stages (8, then 10, then 11",
      "studies); this file encodes the Stage 3 final model, whose parameters",
      "differ from Stage 2 by less than 2%."
    )
  )

  ini({
    # =========================================================================
    # Structural parameters: Comisar 2025 Table 3 'NONMEM Model Estimates'
    # (Stage 3 final model), referenced to a 70 kg adult. These are APPARENT
    # (oral) parameters anchored on a relative bioavailability of 1 for a 10 mg
    # fasted dose -- see the DOSE covariate notes.
    #
    # NOTE ON SOURCING: the supplement's Supplementary Listing 1 is the final
    # model's NONMEM control stream, but its $THETA / $OMEGA / $SIGMA blocks
    # carry INITIAL estimates (e.g. TH1 CL = 17.6 vs the final 24.1 L/h). The
    # control stream is used here ONLY for structure and covariate equation
    # forms; every value below is the final estimate from Table 3.
    # =========================================================================
    lcl  <- log(24.1) ; label("Apparent clearance CL/F (L/h) at 70 kg")                        # Comisar 2025 Table 3 'Apparent clearance, CL/F (L/h)' = 24.1 (%RSE 4.86)
    lvc  <- log(114)  ; label("Apparent central volume of distribution Vc/F (L) at 70 kg")     # Comisar 2025 Table 3 'Apparent central volume of distribution, Vc/F (L)' = 114 (%RSE 5.36)
    lvp  <- log(46.0) ; label("Apparent peripheral volume of distribution Vp/F (L) at 70 kg")  # Comisar 2025 Table 3 'Apparent peripheral volume of distribution, Vp/F (L)' = 46.0 (%RSE 5.3)
    lq   <- log(3.94) ; label("Apparent inter-compartmental clearance Q/F (L/h) at 70 kg")     # Comisar 2025 Table 3 'Apparent inter-compartmental clearance, Q/F (L/h)' = 3.94 (%RSE 6.37)
    lktr <- log(8.23) ; label("Transit absorption rate constant ktr (1/h)")                    # Comisar 2025 Table 3 'Transit absorption rate constant, ktr (h-1)' = 8.23 (%RSE 8.24)
    lka  <- log(3.86) ; label("First-order absorption rate constant ka, terminal absorption compartment to central (1/h)") # Comisar 2025 Table 3 'First-order absorption rate constant, ka (h-1)' = 3.86 (%RSE 28.4)

    # =========================================================================
    # Empirical allometric exponents, FIXED (not estimated). Table 3 note:
    # "CL/F and Q/F parameters have body weight exponents fixed to 0.75, and
    # Vc/F and Vp/F have body weight exponents fixed to 1.0 (empirical
    # allometry)." Control stream: `$THETA 0.75 FIX ; TH22:CLBWGT1` and
    # `$THETA 1 FIX ; TH23:V2BWGT1`.
    # =========================================================================
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent on CL/F and Q/F (unitless)")   # Comisar 2025 Table 3 note; Methods 2.2 paragraph 3; control stream TH22 = 0.75 FIX
    allo_v  <- fixed(1.0)  ; label("Allometric exponent on Vc/F and Vp/F (unitless)")  # Comisar 2025 Table 3 note; Methods 2.2 paragraph 3; control stream TH23 = 1 FIX

    # =========================================================================
    # Covariate effects, Comisar 2025 Table 3 'Covariates' block. Every
    # categorical effect uses the proportional-constant form given in the
    # Supplemental Methods: theta_i = theta_TV,REF * (1 + theta_x)^X_i with
    # X_i the 0/1 indicator -- i.e. multiply by (1 + theta) when the indicator
    # fires. Reference states: normal-or-mild hepatic function, no azole
    # coadministration, fasted, immediate-release tablet, 75 mg dose.
    #
    # The (1 + theta) reading is corroborated four ways: the Supplemental
    # Methods equation, the control stream's `(1 + THETA(n))` literals, the
    # paper's own percentage restatements (e.g. Discussion paragraph 5 "itra-
    # conazole ... decreased rimegepant clearance by 74.3%" for -0.743), and
    # the Table 4 simulated exposure ratios (e.g. AUCtau 16,200 / 4160 = 3.89
    # against 1 / (1 - 0.743) = 3.89).
    # =========================================================================
    e_hepimp_mod_cl <- -0.229 ; label("Fractional effect of moderate hepatic impairment on CL/F")  # Comisar 2025 Table 3 'Moderate hepatic impairment on CL/F' = -0.229 (%RSE 27.8); control stream TH24
    e_hepimp_sev_cl <- -0.423 ; label("Fractional effect of severe hepatic impairment on CL/F")    # Comisar 2025 Table 3 'Severe hepatic impairment on CL/F' = -0.423 (%RSE 24.9); control stream TH11
    e_flcz_cl       <- -0.427 ; label("Fractional effect of fluconazole coadministration on CL/F") # Comisar 2025 Table 3 'Fluconazole use on CL/F' = -0.427 (%RSE 2.88); control stream TH8
    e_itcz_cl       <- -0.743 ; label("Fractional effect of itraconazole coadministration on CL/F") # Comisar 2025 Table 3 'Itraconazole use on CL/F' = -0.743 (%RSE 1.25); control stream TH12

    e_itcz_ktr      <- -0.351 ; label("Fractional effect of itraconazole coadministration on ktr")  # Comisar 2025 Table 3 'Itraconazole use on ktr' = -0.351 (%RSE 28.9); control stream TH19 (KAITCZ, which multiplies KTR)
    e_fed_ktr       <- -0.706 ; label("Fractional effect of fed status on ktr")                     # Comisar 2025 Table 3 'Fed-status on ktr' = -0.706 (%RSE 3.34); control stream TH16
    e_form_odt_ktr  <-  0.355 ; label("Fractional effect of the orally disintegrating tablet on ktr") # Comisar 2025 Table 3 'Oral disintegrating tablet on ktr' = 0.355 (%RSE 29.6); control stream TH17 (FORM = 3)
    e_form_cap_ktr  <-  2.03  ; label("Fractional effect of the capsule formulation on ktr")        # Comisar 2025 Table 3 'Capsule formulation on ktr' = 2.03 (%RSE 24.1); control stream TH18 (FORM = 2)
    e_dose_low_ktr  <-  0.536 ; label("Fractional effect of a low (10 or 25 mg) dose level on ktr") # Comisar 2025 Table 3 '10/25mg dose effect on ktr' = 0.536 (%RSE 48.7); control stream KTRDOSE2and3, applied as (1 + THETA(29)) for DOSE = 10 and 25

    e_fed_f         <- -0.315 ; label("Fractional effect of fed status on relative bioavailability F1") # Comisar 2025 Table 3 'Fed-status on relative bioavailability (F1)' = -0.315 (%RSE 9.74); control stream TH13
    e_dose_f        <-  0.191 ; label("Power exponent for dose on relative bioavailability F1 (unitless)") # Comisar 2025 Table 3 'Dose effect on F1' = 0.191 (%RSE 12.6); control stream `F1BIODOSE = (DOSE/10)**THETA(25)`; Methods 2.2 paragraph 4 'F1 with dose impact = F1 without dose impact x (dose/10 mg)^thetaF1'

    # =========================================================================
    # Inter-individual variability, Comisar 2025 Table 3 'IIV% (%RSE)' column
    # and 'Correlation matrix' block. IIV on CL/F, Vc/F, and Vp/F is a
    # correlated BLOCK(3); IIV on ktr is diagonal (control stream: `$OMEGA
    # BLOCK(3)` for ETA1-ETA3 then a separate `$OMEGA` for ETA4).
    #
    # OMEGA SCALE: the reported IIV% is the approximate CV on the SD scale,
    # i.e. omega^2 = (IIV%/100)^2, NOT sqrt(exp(omega^2) - 1). Two lines of
    # evidence: (a) round-tripping the control stream's initial OMEGA
    # diagonals row by row under this reading reproduces all four reported
    # finals within 2.5% (0.0979 -> 31.3% vs 31.2; 0.1577 -> 39.7% vs 40.6;
    # 0.0833 -> 28.9% vs 28.2; 0.2972 -> 54.5% vs 53.1), whereas the
    # exp-based reading misses ktr by 10% (58.8% vs 53.1); and (b) the same
    # table's residual-error row is unambiguously on the SD scale -- the
    # initial $SIGMA proportional variance 0.1581 gives sqrt = 39.8% against
    # the reported final 37.7%, and the additive variance 0.1744 gives
    # sqrt = 0.418 against the reported final 0.43.
    #
    # Off-diagonals are reconstructed from the reported CORRELATIONS as
    # cov(i,j) = corr(i,j) * omega_i * omega_j:
    #   CL/F:Vc/F 0.758 * 0.312 * 0.406 = 0.096017
    #   CL/F:Vp/F 0.540 * 0.312 * 0.282 = 0.047511
    #   Vc/F:Vp/F 0.667 * 0.406 * 0.282 = 0.076366
    # The resulting correlation matrix is positive definite (eigenvalues all
    # > 0; determinant 0.235).
    # =========================================================================
    etalcl + etalvc + etalvp ~ c(
      0.097344,
      0.096017, 0.164836,
      0.047511, 0.076366, 0.079524
    ) # Comisar 2025 Table 3: IIV 31.2% on CL/F (%RSE 5.1), 40.6% on Vc/F (%RSE 6.1), 28.2% on Vp/F (%RSE 8.2); correlations w1,2 CL/F:Vc/F = 75.8% (%RSE 6.8), w1,3 CL/F:Vp/F = 54.0% (%RSE 9.2), w2,3 Vc/F:Vp/F = 66.7% (%RSE 8.4)
    etalktr ~ 0.281961 # Comisar 2025 Table 3: IIV 53.1% on ktr (%RSE 6.9); 0.531^2 = 0.281961

    # No IIV was estimated on Q/F or ka (Table 3 shows '-' in the IIV column
    # for both; the control stream sets Q = TVQ and KA = TVKA with no ETA).

    # =========================================================================
    # Residual unexplained variability: combined proportional + additive on the
    # linear concentration scale (control stream $ERROR `Y = F + F*EP + EA`).
    # =========================================================================
    propSd <- 0.377 ; label("Proportional residual error (fraction)")  # Comisar 2025 Table 3 'Proportional error (%)' = 37.7% (%RSE 1.9)
    addSd  <- 0.43  ; label("Additive residual error (ng/mL)")         # Comisar 2025 Table 3 'Additive error' = 0.43 (%RSE 12.0). Table 3 prints the unit as 'ng/L'; that is a typo for ng/mL -- see the vignette Errata. The control stream scales the prediction to ng/mL (`S2 = V2/1000 ; ng/ml to mg in final volume units of L`), its initial additive $SIGMA variance 0.1744 has sqrt = 0.418 on that same ng/mL scale, and 0.43 ng/mL sits just below the 0.5 ng/mL LLOQ as expected for an additive term, whereas 0.43 ng/L would be 0.00043 ng/mL and numerically inert.
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Derived covariate multipliers. Each categorical effect enters as
    #    (1 + theta) when its indicator fires, per the Supplemental Methods
    #    categorical covariate model and the control stream's `(1 + THETA(n))`
    #    literals. The products below mirror the control stream's
    #    `CLCOV = CLFLCZ*CLFORM*CLHEPATC*CLITCZ*CLBWGT`,
    #    `KTRCOV = KAFED*KAFORM*KAITCZ*KAFLCZ*KTRDOSE`, and
    #    `F1BIOCOV = F1BIOFED*F1BIOFORM*F1BIODOSE`, with the terms the paper
    #    fixed to 0 (formulation on CL/F, formulation on F1, fluconazole on
    #    ktr) omitted because they evaluate to exactly 1.
    # -----------------------------------------------------------------------
    clCov <- (1 + e_flcz_cl * CONMED_FLUCONAZOLE) *
      (1 + e_itcz_cl * CONMED_ITRACONAZOLE) *
      (1 + e_hepimp_mod_cl * HEPIMP_MOD) *
      (1 + e_hepimp_sev_cl * HEPIMP_SEV)

    # Low-dose indicator for the ktr effect. The control stream fires
    # `KTRDOSE = (1 + THETA(29))` for DOSE = 10 and DOSE = 25 only, with
    # 75 mg as the reference and the 150 mg theta fixed to 0. Over the four
    # dose levels studied (10, 25, 75, 150 mg) the 50 mg cut is exact.
    doseLow <- (DOSE < 50)

    ktrCov <- (1 + e_fed_ktr * FED) *
      (1 + e_form_odt_ktr * FORM_ODT) *
      (1 + e_form_cap_ktr * FORM_CAPSULE) *
      (1 + e_itcz_ktr * CONMED_ITRACONAZOLE) *
      (1 + e_dose_low_ktr * doseLow)

    # Relative bioavailability: a continuous power term on dose (reference
    # 10 mg, where F1 = 1) times the fed-status fractional effect.
    f1Cov <- (DOSE / 10)^e_dose_f * (1 + e_fed_f * FED)

    # -----------------------------------------------------------------------
    # 2. Individual parameters with empirical allometric body-weight scaling
    #    (reference 70 kg). Q/F and ka carry no IIV.
    # -----------------------------------------------------------------------
    cl  <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * clCov
    vc  <- exp(lvc + etalvc) * (WT / 70)^allo_v
    vp  <- exp(lvp + etalvp) * (WT / 70)^allo_v
    q   <- exp(lq)           * (WT / 70)^e_wt_cl
    ktr <- exp(lktr + etalktr) * ktrCov
    ka  <- exp(lka)

    # -----------------------------------------------------------------------
    # 3. Two-compartment disposition with a six-compartment pre-systemic
    #    chain, transcribed from the control stream's $DES block:
    #
    #      DADT(1) = -KTR*A(1)                       -> depot
    #      DADT(4) = KTR*A(1) - KTR*A(4)             -> transit1
    #      DADT(5) = KTR*A(4) - KTR*A(5)             -> transit2
    #      DADT(6) = KTR*A(5) - KTR*A(6)             -> transit3
    #      DADT(7) = KTR*A(6) - KTR*A(7)             -> transit4
    #      DADT(8) = KTR*A(7) - KA*A(8)              -> transit5
    #      DADT(2) = KA*A(8) - CL*C2 - Q*C2 + Q*C3   -> central
    #      DADT(3) = Q*C2 - Q*C3                     -> peripheral1
    #
    #    Canonical naming: `depot` is the control stream's COMP(1) "DOSING AND
    #    TRANSIT 0" (it receives the dose and drains at ktr); `transit1` to
    #    `transit4` are its COMP(4)-COMP(7) "TRANSIT 1"-"TRANSIT 4"; and
    #    `transit5` is its COMP(8) "ABSORPTION", the terminal pre-systemic
    #    compartment, which drains at ka rather than ktr. Naming that last
    #    compartment `transit5` follows the house idiom established by
    #    Bukkems_2021_raltegravir, where the terminal transit compartment
    #    likewise empties into central at ka. So there are five sequential
    #    ktr transfers and one ka transfer, even though the paper counts only
    #    the four explicitly labelled transit compartments ("n ~ 4" in
    #    Figure 2).
    # -----------------------------------------------------------------------
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3
    d/dt(transit4)    <-  ktr * transit3 - ktr * transit4
    d/dt(transit5)    <-  ktr * transit4 - ka  * transit5
    d/dt(central)     <-  ka  * transit5 - (cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    # -----------------------------------------------------------------------
    # 4. Relative bioavailability on the dosing compartment (control stream
    #    `F1 = TVF1BIO`, which in NONMEM applies to COMP(1) = depot here).
    # -----------------------------------------------------------------------
    f(depot) <- f1Cov

    # -----------------------------------------------------------------------
    # 5. Observation in ng/mL. Amounts are in mg and vc in L, so central/vc is
    #    mg/L and the factor 1000 converts to ng/mL (= ug/L), matching the
    #    control stream's `S2 = V2/1000`.
    # -----------------------------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
