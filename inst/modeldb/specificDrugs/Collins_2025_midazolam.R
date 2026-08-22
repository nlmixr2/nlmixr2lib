Collins_2025_midazolam <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral midazolam and its",
    "1-OH-midazolam metabolite in 72 healthy volunteers given a 1 mg oral",
    "midazolam CYP3A probe on two occasions - once with a single 600 mg dose",
    "of efavirenz and again after 17 days of efavirenz 600 mg/day (Collins",
    "2025). Both analytes have two-compartment disposition; midazolam is",
    "absorbed first-order from a depot and a fixed fraction (0.7) of the",
    "midazolam elimination flux forms 1-OH-midazolam. Chronic (multiple-dose)",
    "efavirenz is a covariate on midazolam clearance, on the absorption rate",
    "constant, and on bioavailability; CYP3A5 expresser status, CYP3A4*22",
    "carriage, and female sex are covariates on midazolam clearance. All",
    "clearances and volumes carry fixed allometric scaling on body weight",
    "(0.75 and 1) referenced to 73 kg. Amounts are nmol and concentrations",
    "nM, as fitted. Parameter values are taken from the deposited final",
    "NONMEM control stream (Supporting Information File S5), which for six",
    "parameters disagrees with the back-transformed values printed in Table 2",
    "- see the vignette Assumptions and deviations section."
  )
  reference <- paste(
    "Collins KS, Aruldhas BW, Metzger IF, Lu JBL, Heathman MA, Quinney SK,",
    "Desta Z. A Population Pharmacokinetic Approach to Understand the Effect",
    "of Efavirenz on CYP3A Activity in Healthy Volunteers Using Midazolam as",
    "a Probe. CPT Pharmacometrics Syst Pharmacol. 2025;14:2095-2106.",
    "doi:10.1002/psp4.70116.",
    "Parameter values from Supporting Information File S5 (final NONMEM",
    "control stream).",
    sep = " "
  )
  vignette <- "Collins_2025_midazolam"
  units <- list(
    time          = "h",
    dosing        = "nmol",
    concentration = "nM"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling with fixed exponents, referenced to 73 kg (the",
        "Table 1 median). File S5 codes the scaling on the log scale inside",
        "the MU-referenced typical value, e.g.",
        "TVCLP = THETA(2) + LOG(WT/73)*(3/4) and",
        "TVV2 = THETA(3) + LOG(WT/73), which is algebraically",
        "(WT/73)^0.75 on every clearance and (WT/73)^1 on every volume,",
        "for both midazolam and 1-OH-midazolam. Methods 'Model Development':",
        "'Weight influence was accounted for using allometric scaling (fixed",
        "exponents of 0.75 for clearance and 1 for volume) across all",
        "clearance and volume parameter estimates.'",
        "Table 1 weight median (range): 73 (53-104) kg."
      ),
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Self-reported female sex",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "1 = female, 0 = male. Multiplicative effect on midazolam clearance.",
        "The source dataset codes the opposite polarity: File S5 sets",
        "K0 = 1 IF(SEX.EQ.0) with the trailing comment ';female', so the",
        "paper's SEX column is 0 = female / 1 = male and",
        "SEXF = 1 - SEX. The effect enters as (1 + THETA(14)*K0), i.e. it is",
        "applied to females with males as the reference, which matches the",
        "Discussion ('clearance of midazolam was slightly but significantly",
        "higher in females compared to males') and the Table 2 row label",
        "'Female'. Table 1: 27/72 (37.5%) female in the single-dose cohort."
      ),
      source_name        = "SEX (0 = female, 1 = male)"
    ),
    CYP3A5_EXPR = list(
      description        = "CYP3A5 expresser status",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5 non-expresser)",
      notes              = paste(
        "1 = expresser, 0 = non-expresser. Methods 'DNA Genotyping':",
        "'CYP3A5 genotypes were classified as non-expressors (two copies of",
        "*3, *6, or *7 alleles) and expressors (one or no copies of *3, *6,",
        "and *7 alleles).' This pools *1/*1 and *1/*3 into the expresser",
        "group, exactly the dichotomy the CYP3A5_EXPR canonical encodes.",
        "Note that this study genotyped *6 (rs10264272) and *7",
        "(rs41303343) in addition to *3 (rs776746), so a non-expresser here",
        "may carry a loss-of-function allele other than *3.",
        "File S5 column A5, with N1 = 1 IF(A5.EQ.1) ';CYP3A5 Expressors'.",
        "Table 1: 18/72 (25%) expressers in the single-dose cohort."
      ),
      source_name        = "A5"
    ),
    SNP_CYP3A4_RS35599367 = list(
      description        = "CYP3A4*22 (rs35599367) reduced-function allele carrier",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A4 *1/*1 normal metabolizer)",
      notes              = paste(
        "1 = carries at least one *22 allele, 0 = *1/*1. Methods 'DNA",
        "Genotyping': 'CYP3A4 genotypes were categorized as normal",
        "metabolizers (*1/*1) and intermediate metabolizers (*1/*22)' and",
        "genotyping was 'for CYP3A4*22 (rs35599367)'. The cohort contained",
        "no *22/*22 subjects, so the paper's 'intermediate metabolizer'",
        "stratum is exactly the *22 carrier stratum this canonical encodes.",
        "File S5 column A4, with M1 = 1 IF(A4.EQ.1) ';CYP3A4 IM'.",
        "Table 1: 4/72 (6%) intermediate metabolizers in the single-dose",
        "cohort; the Discussion notes the effect 'did not reach statistical",
        "significance, likely due to the small number of intermediate",
        "metabolizers enrolled in the study (n = 4)'."
      ),
      source_name        = "A4"
    ),
    CONMED_EFV_MD = list(
      description        = "Multiple-dose (steady-state) efavirenz co-administration",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (single 600 mg dose of efavirenz)",
      notes              = paste(
        "1 = midazolam probe dose given on study Day 24, after efavirenz",
        "600 mg once daily on Days 7-23 plus the Day 24 morning dose;",
        "0 = midazolam probe dose given on Day 1, one hour after a single",
        "600 mg efavirenz dose. Note that BOTH levels are efavirenz-exposed:",
        "the reference is single-dose, not efavirenz-free, so this is not",
        "the same contrast as the CONMED_EFV canonical. The paper",
        "acknowledges this in the Discussion ('Our study design did not",
        "include a midazolam-alone arm ... the true magnitude of this change",
        "may be underestimated'). Occasion-level, not subject-level: 58 of",
        "the 72 subjects contributed both levels.",
        "File S5 column EFV with reversed polarity - L0 = 1 IF(EFV.EQ.0)",
        "';multiple dose efavirenz' - so CONMED_EFV_MD = 1 - EFV.",
        "The Results text states this explicitly: 'In the NONMEM control",
        "stream, EFV = 0 represents the multiple-dose efavirenz condition,",
        "and EFV = 1 represents the single-dose condition ... The reference",
        "(baseline) state is thus the single-dose efavirenz condition.'"
      ),
      source_name        = "EFV (0 = multiple dose, 1 = single dose)"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "midazolam", units = "nmol",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "midazolam", units = "nmol", specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "midazolam", units = "nmol", specimen = "plasma", verified = TRUE
    ),
    central_1ohm = list(
      analyte = "1-OH-midazolam", units = "nmol", specimen = "plasma", verified = TRUE
    ),
    peripheral1_1ohm = list(
      analyte = "1-OH-midazolam", units = "nmol", specimen = "plasma", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 72,
    n_studies      = 1,
    age_range      = "18-50 years",
    age_median     = "25 years",
    weight_range   = "53-104 kg",
    weight_median  = "73 kg",
    sex_female_pct = 37.5,
    race_ethnicity = c(White = 72, Black = 21, Other = 7),
    disease_state  = "healthy volunteers",
    dose_range     = paste(
      "Midazolam 1 mg oral syrup as a single CYP3A probe dose on each of two",
      "occasions, given one hour after efavirenz. Efavirenz 600 mg oral:",
      "a single dose on Day 1 and 600 mg once daily on Days 7-23 with a",
      "final dose on Day 24."
    ),
    regions        = "United States (Indiana CTSI Clinical Research Center, Indianapolis)",
    notes          = paste(
      "Demographics from Table 1 (single-dose efavirenz column, n = 72).",
      "All 72 subjects completed the Day 1 (single-dose efavirenz) session;",
      "58 completed both sessions, so the multiple-dose efavirenz occasion",
      "has n = 58. Body mass index 20-32 kg/m^2 was an eligibility criterion",
      "(Table 1 median 24, range 18-32). The probe cocktail also contained",
      "150 mg caffeine, 250 mg tolbutamide and 20 mg omeprazole; only",
      "midazolam and 1-OH-midazolam are modelled here. Plasma was assayed",
      "after beta-glucuronidase deconjugation, so the 1-OH-midazolam",
      "measurements are total (free plus glucuronide) metabolite. Dataset:",
      "2309 plasma concentrations (1153 midazolam, 1156 1-OH-midazolam);",
      "33% of midazolam and 0.2% of 1-OH-midazolam concentrations were below",
      "the 0.3 ng/mL limit of quantification and were retained uncensored",
      "under the Keizer 'all data' approach. Estimation used SAEM followed",
      "by importance sampling in NONMEM 7.4 with PsN 4.9.01. Base model OFV",
      "9166.31. CYP2B6 genotype was assayed but deliberately not carried",
      "into the model (Methods 'Model Development'). ClinicalTrials.gov",
      "NCT00668395."
    )
  )

  ini({
    # =================================================================
    # SOURCING NOTE (read before changing any value below).
    #
    # Every value in this ini() block is taken from the deposited final
    # NONMEM control stream, Supporting Information File S5
    # ('$PROBLEM Final'), whose $THETA / $OMEGA / $SIGMA blocks hold the
    # final estimates (PsN update_inits). That file agrees exactly with
    # Table 2 'Final model estimate' for all eight structural typical
    # values, all ten omegas and both sigmas.
    #
    # It DISAGREES with Table 2 for the five covariate coefficients and
    # for the multiple-dose bioavailability, because Table 2 reports
    # exp(THETA) for every THETA while File S5 codes those six as
    # proportional (1 + THETA) multipliers and as a logit-scale
    # bioavailability. Figure 4 (the paper's own simulated single-dose
    # vs multiple-dose profiles) resolves the conflict in favour of
    # File S5 for both analytes; see the vignette's 'Assumptions and
    # deviations' section for the digitised comparison. Table 2's Ka
    # (0.77) is the one value that reconciles with neither reading of
    # THETA(1) = -0.303 (exp = 0.7386); File S5 is used.
    # =================================================================

    # -----------------------------------------------------------------
    # Midazolam absorption and disposition. Typical values are for a
    # 73 kg male CYP3A4 *1/*1, CYP3A5 non-expresser on single-dose
    # efavirenz.
    # -----------------------------------------------------------------
    lka <- -0.303
    label("Midazolam absorption rate constant, single-dose efavirenz (1/h)")   # File S5 $THETA(1) = -0.303 ';KA SD'; exp = 0.7386 1/h. Table 2 prints 0.77 (bootstrap 0.77, 95% CI 0.71-0.84), which is not exp(-0.303); see SOURCING NOTE.
    lcl <- 3.38
    label("Midazolam clearance, reference subject (L/h)")                      # File S5 $THETA(2) = 3.38 ';CLP MD'; exp = 29.370 L/h. Matches Table 2 'CL MDZ, SD EFV, NM, Nonexp, 73kg, Males' = 29.37 and the Results text 'estimated at 29.37 L/h'. This is a true CL, not CL/F: F is coded explicitly in the model.
    lvc <- 1.67
    label("Midazolam central volume (L)")                                      # File S5 $THETA(3) = 1.67 ';V2'; exp = 5.3122 L. Table 2 'Vc MDZ' = 5.31.
    lvp <- 4.42
    label("Midazolam peripheral volume (L)")                                   # File S5 $THETA(4) = 4.42 ';V4'; exp = 83.096 L. Table 2 'Vp MDZ' = 83.10.
    lq <- 2.73
    label("Midazolam intercompartmental clearance (L/h)")                      # File S5 $THETA(5) = 2.73 ';QP'; exp = 15.333 L/h. Table 2 'Q MDZ' = 15.33.

    # -----------------------------------------------------------------
    # 1-OH-midazolam disposition (two-compartment; the metabolite model
    # was built on the parent predictions assuming linear formation,
    # Methods 'Model Development').
    # -----------------------------------------------------------------
    lcl_1ohm <- 0.451
    label("1-OH-midazolam clearance (L/h)")                                    # File S5 $THETA(6) = 0.451 ';CLM'; exp = 1.5699 L/h. Table 2 "CL 1'-OH-MDZ" = 1.57.
    lvc_1ohm <- -0.416
    label("1-OH-midazolam central volume (L)")                                 # File S5 $THETA(7) = -0.416 ';V3'; exp = 0.65965 L. Table 2 "Vc 1'-OH-MDZ" = 0.66.
    lvp_1ohm <- 3.02
    label("1-OH-midazolam peripheral volume (L)")                              # File S5 $THETA(8) = 3.02 ';V5'; exp = 20.492 L. Table 2 "Vp 1'-OH-MDZ" = 20.49.
    lq_1ohm <- -0.494
    label("1-OH-midazolam intercompartmental clearance (L/h)")                 # File S5 $THETA(9) = -0.494 ';QM'; exp = 0.61019 L/h. Table 2 "Q 1'-OH-MDZ" = 0.61.

    # -----------------------------------------------------------------
    # Fraction of the midazolam elimination flux forming
    # 1-OH-midazolam. Held constant, so the metabolite CL and volumes
    # above are identifiable only up to this value.
    # -----------------------------------------------------------------
    lfm <- fixed(log(0.7))
    label("Fraction of midazolam metabolised to 1-OH-midazolam (unitless)")    # File S5 $PK 'FMET=0.7; fraction of parent metabolized to 1-OH metabolite'. Results: 'A fixed fraction of 0.7 was used for the metabolism of midazolam into 1'-hydroxy-midazolam [40, 41].' Both amounts are molar (nmol), so no molecular-weight correction is needed.

    # -----------------------------------------------------------------
    # Bioavailability. The two efavirenz occasions carry structurally
    # different bioavailability terms in File S5, so each stratum gets
    # its own parameter (parameter-names.md 'Stratum-suffixed
    # parameters'): the single-dose value is a fixed constant with no
    # between-subject variability, and the multiple-dose value is
    # estimated on the logit scale and does carry IIV.
    # -----------------------------------------------------------------
    lfdepot_sdefv <- fixed(log(0.5))
    label("Midazolam oral bioavailability, single-dose efavirenz (unitless)")  # File S5 $PK 'F1=0.5' outside the IF block; not estimated. Methods 'Model Development': 'For midazolam, bioavailability was fixed at 0.5 based on available literature after a single dose of efavirenz, with relative bioavailability estimated after multiple doses [39].' Table 2 'F SD-EFV (fixed)' = 0.50.
    logitfdepot_mdefv <- -0.511
    label("Midazolam oral bioavailability, multiple-dose efavirenz (logit scale)")  # File S5 $THETA(10) = -0.511 ';F1 MD', entering as F1 = EXP(LF1)/(1+EXP(LF1)) so F = expit(-0.511) = 0.3749. Table 2 prints 0.60 = exp(-0.511), i.e. the generic exp() back-transform applied to a logit-scale THETA; Figure 4 falsifies 0.60 and confirms 0.375. See SOURCING NOTE and the vignette.

    # -----------------------------------------------------------------
    # Allometric exponents - fixed at the canonical 0.75 / 1 and shared
    # across midazolam and 1-OH-midazolam.
    # -----------------------------------------------------------------
    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent on body weight for all clearances (unitless)")  # File S5 $PK: LOG(WT/73)*(3/4) added to TVCLP, TVQP, TVCLM and TVQM. Methods 'Model Development': 'fixed exponents of 0.75 for clearance'.
    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent on body weight for all volumes (unitless)")     # File S5 $PK: LOG(WT/73) added to TVV2, TVV4, TVV3 and TVV5. Methods 'Model Development': 'fixed exponents of ... 1 for volume'.

    # -----------------------------------------------------------------
    # Covariate effects. All four clearance covariates and the single
    # absorption covariate enter as proportional (1 + coefficient)
    # multipliers in File S5, i.e. the coefficient is the fractional
    # change relative to the reference subject.
    # -----------------------------------------------------------------
    e_conmed_efv_md_cl <- 0.652
    label("Fractional change in midazolam CL with multiple-dose efavirenz (unitless)")  # File S5 $THETA(11) = 0.652 ';CL MD' in CLP = ... *(1+L0*THETA(11)); multiplier 1.652. The Study Highlights give '1.64-fold'. Table 2 / abstract / Figure 5 report 1.92 = exp(0.652); Figure 4 supports 1.652 (see SOURCING NOTE).
    e_snp_cyp3a4_rs35599367_cl <- -0.117
    label("Fractional change in midazolam CL for CYP3A4*22 carriers (unitless)")        # File S5 $THETA(12) = -0.117 ';CL IM' in CLP = ... *(1+M1*THETA(12)); multiplier 0.883. Table 2 prints 0.89 = exp(-0.117); bootstrap median 0.94.
    e_cyp3a5_expr_cl <- 0.254
    label("Fractional change in midazolam CL for CYP3A5 expressers (unitless)")         # File S5 $THETA(13) = 0.254 ';CL Exp' in CLP = ... *(1+N1*THETA(13)); multiplier 1.254. Table 2 prints 1.29 = exp(0.254); bootstrap median 1.27.
    e_sexf_cl <- 0.199
    label("Fractional change in midazolam CL for females (unitless)")                   # File S5 $THETA(14) = 0.199 ';CL F' in CLP = ... *(1+K0*THETA(14)); multiplier 1.199. Table 2 prints 1.22 = exp(0.199); bootstrap median 1.30.
    e_conmed_efv_md_ka <- 0.255
    label("Fractional change in midazolam Ka with multiple-dose efavirenz (unitless)")  # File S5 $THETA(15) = 0.255 ';KA MD' in KA = ... *(1+L0*THETA(15)); multiplier 1.255. Table 2 prints 1.29 = exp(0.255).

    # -----------------------------------------------------------------
    # Between-subject variability. File S5 $OMEGA is a diagonal matrix
    # of variances on the (log or logit) transformed scale, in the same
    # order as ETA(1)..ETA(10). Table 2 reports these as %CV via
    # sqrt(exp(omega) - 1) * 100; each cross-check is given below.
    # -----------------------------------------------------------------
    etalka ~ 0.0122                    # File S5 $OMEGA ';KA'      -> sqrt(exp(0.0122)-1)*100 = 11.1% CV; Table 2 'Ka' 11.1%
    etalcl ~ 0.132                     # File S5 $OMEGA ';CLP'     -> 37.6% CV; Table 2 'CL MDZ' 37.6%
    etalvc ~ 0.865                     # File S5 $OMEGA ';V2'      -> 117.3% CV; Table 2 'Vc MDZ' 117.3%
    etalvp ~ 0.359                     # File S5 $OMEGA ';V4'      -> 65.7% CV; Table 2 'Vp MDZ' 65.7%
    etalq ~ 0.111                      # File S5 $OMEGA ';QP'      -> 34.3% CV; Table 2 'Q MDZ' 34.3%
    # NOTE: eta lines carry no label(), so rxode2 turns these trailing comments
    # into label() calls. Double quotes here make the generated code unparseable
    # -- keep the metabolite row names unquoted.
    etalcl_1ohm ~ 0.0521               # File S5 $OMEGA ';CLM'     -> 23.1% CV; Table 2 row CL 1-OH-MDZ 23.1%
    etalvc_1ohm ~ 0.0488               # File S5 $OMEGA ';V3'      -> 22.4% CV; Table 2 row Vc 1-OH-MDZ 22.4%
    etalvp_1ohm ~ 0.595                # File S5 $OMEGA ';V5'      -> 90.2% CV; Table 2 row Vp 1-OH-MDZ 90.2%
    etalq_1ohm ~ 0.229                 # File S5 $OMEGA ';QM'      -> 50.7% CV; Table 2 row Q 1-OH-MDZ 50.7%
    etalogitfdepot_mdefv ~ 0.134       # File S5 $OMEGA ';F1'      -> 37.9% CV; Table 2 'F' 37.9%. On the LOGIT scale, and only active on the multiple-dose efavirenz occasion (File S5 puts ETA(10) inside the IF(EFV.EQ.0) branch), so single-dose bioavailability has no IIV.

    # -----------------------------------------------------------------
    # Residual error - proportional for both analytes. File S5 $SIGMA
    # holds variances; Table 2 reports them as %CV = sqrt(sigma)*100,
    # so the proportional SD is sqrt(sigma). A combined additive plus
    # proportional model was tried and rejected (Results).
    # -----------------------------------------------------------------
    propSd <- sqrt(0.0779)
    label("Midazolam proportional residual SD (fraction)")                     # File S5 $SIGMA ';Parent Prop' = 0.0779 -> sqrt = 0.2791; Table 2 'Proportional error MDZ (%CV)' 27.9%
    propSd_1ohm <- sqrt(0.0688)
    label("1-OH-midazolam proportional residual SD (fraction)")                # File S5 $SIGMA ';Metabolite Prop' = 0.0688 -> sqrt = 0.2623; Table 2 "Proportional error 1'-OH-MDZ (%CV)" 26.2%
  })

  model({
    # ---------------------------------------------------------------
    # 1. Derived covariate terms.
    # ---------------------------------------------------------------
    wt_ref <- 73                       # Table 1 median body weight, kg; the File S5 allometric denominator
    wt_cl  <- (WT / wt_ref)^e_wt_cl_q  # applies to every clearance (parent and metabolite)
    wt_v   <- (WT / wt_ref)^e_wt_vc_vp # applies to every volume (parent and metabolite)

    # Multiplicative covariate factor on midazolam clearance, in the
    # File S5 $PK order:
    #   CLP = EXP(MU_2+ETA(2)) * (1+L0*THETA(11)) * (1+M1*THETA(12))
    #                          * (1+N1*THETA(13)) * (1+K0*THETA(14))
    cl_cov <- (1 + e_conmed_efv_md_cl * CONMED_EFV_MD) *
      (1 + e_snp_cyp3a4_rs35599367_cl * SNP_CYP3A4_RS35599367) *
      (1 + e_cyp3a5_expr_cl * CYP3A5_EXPR) *
      (1 + e_sexf_cl * SEXF)

    # ---------------------------------------------------------------
    # 2. Individual parameters.
    # ---------------------------------------------------------------
    ka <- exp(lka + etalka) * (1 + e_conmed_efv_md_ka * CONMED_EFV_MD)
    cl <- exp(lcl + etalcl) * wt_cl * cl_cov
    vc <- exp(lvc + etalvc) * wt_v
    vp <- exp(lvp + etalvp) * wt_v
    q  <- exp(lq + etalq) * wt_cl

    cl_1ohm <- exp(lcl_1ohm + etalcl_1ohm) * wt_cl
    vc_1ohm <- exp(lvc_1ohm + etalvc_1ohm) * wt_v
    vp_1ohm <- exp(lvp_1ohm + etalvp_1ohm) * wt_v
    q_1ohm  <- exp(lq_1ohm + etalq_1ohm) * wt_cl

    fm <- exp(lfm)

    # Bioavailability switches wholesale between the two efavirenz
    # occasions: a fixed 0.5 with no IIV on the single-dose occasion,
    # and a logit-scale estimate with IIV on the multiple-dose
    # occasion. This reproduces the File S5 branch
    #   F1 = 0.5
    #   IF(EFV.EQ.0) THEN ... F1 = EXP(LF1)/(1+EXP(LF1))
    # without an if(), so etalogitfdepot_mdefv has no effect when
    # CONMED_EFV_MD is 0.
    fdepot_sdefv <- exp(lfdepot_sdefv)
    fdepot_mdefv <- expit(logitfdepot_mdefv + etalogitfdepot_mdefv)
    fdepot <- fdepot_sdefv * (1 - CONMED_EFV_MD) + fdepot_mdefv * CONMED_EFV_MD

    # ---------------------------------------------------------------
    # 3. Micro-constants (File S5 $PK KEP / K24 / K42 / KEM / K35 / K53).
    # ---------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    kel_1ohm <- cl_1ohm / vc_1ohm
    k12_1ohm <- q_1ohm / vc_1ohm
    k21_1ohm <- q_1ohm / vp_1ohm

    # ---------------------------------------------------------------
    # 4. ODE system - File S5 $DES, with COMP 1/2/3/4/5 =
    #    DEPOT/PARENT/METABOLITE/PPERIH/MPERIH. A fixed fraction fm of
    #    the midazolam elimination flux appears as 1-OH-midazolam;
    #    the remaining (1 - fm) leaves the system.
    # ---------------------------------------------------------------
    d/dt(depot)            <- -ka * depot
    d/dt(central)          <-  ka * depot - kel * central -
                               k12 * central + k21 * peripheral1
    d/dt(peripheral1)      <-  k12 * central - k21 * peripheral1
    d/dt(central_1ohm)     <-  fm * kel * central - kel_1ohm * central_1ohm -
                               k12_1ohm * central_1ohm + k21_1ohm * peripheral1_1ohm
    d/dt(peripheral1_1ohm) <-  k12_1ohm * central_1ohm - k21_1ohm * peripheral1_1ohm

    # ---------------------------------------------------------------
    # 5. Bioavailability on the oral depot.
    # ---------------------------------------------------------------
    f(depot) <- fdepot

    # ---------------------------------------------------------------
    # 6. Observations. States are nmol and volumes L, so both
    #    concentrations are nM - the scale File S5 fits in
    #    ('S2=V2 ;Dose and concentrations in nM').
    # ---------------------------------------------------------------
    Cc      <- central / vc
    Cc_1ohm <- central_1ohm / vc_1ohm

    Cc      ~ prop(propSd)
    Cc_1ohm ~ prop(propSd_1ohm)
  })
}
