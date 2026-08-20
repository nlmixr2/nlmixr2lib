Pei_2023_tacrolimus_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, SimBiology 5.8.2). Perfusion-limited 15-compartment",
    "physiologically based model for oral tacrolimus in adult heart",
    "transplant recipients (Pei 2023 Pharmaceutics). Thirteen well-stirred",
    "flow-limited tissues (gut, spleen, pancreas, liver, muscle, kidney,",
    "brain, heart, lung, skin, tendon, other, adipose) plus arterial and",
    "venous blood, on the Levitt PKQuest standard-human physiology (organ",
    "weights and per-kg perfusions from Table S1, 70 kg / 25 percent body",
    "fat reference). Gut, spleen and pancreas drain into the liver via the",
    "portal vein; the liver is the only eliminating organ. Rodgers and",
    "Rowland tissue-to-plasma partition coefficients (Table S2) are",
    "converted to tissue-to-blood by the blood-to-plasma ratio BPR and",
    "multiplied by a Kp scaling factor. BPR is hematocrit-dependent",
    "through the red-cell binding capacity Bmax and affinity constant KD",
    "(Eq 5). Hepatic blood clearance follows the well-stirred extraction",
    "ratio (Eq 3) scaled by the CYP3A5 / CYP3A4 metabolic fractions and",
    "genotype-specific activity levels (Eq 6), with an optional reversible",
    "voriconazole inhibition term (supplement Eq 1-2) driven by a",
    "user-supplied voriconazole whole-blood concentration. Absorption is",
    "first-order from a depot into the gut compartment with a fixed",
    "absorbed fraction Fg = 0.2 (Eq 2). Ka, KD, Bmax and CLint are the",
    "four fitted parameters (Table 3, heart-transplant model-building",
    "column); no IIV or residual-error variance was reported. Observation",
    "is venous whole blood. NOTE: the paper's printed Kp scaling factor of",
    "350 is refuted by the paper's own Tables S5 and 6 by 4-11 fold; the",
    "model uses 9.15 = 11.9 / 1.3, the mouse-to-human average-Kp ratio the",
    "Methods sentence motivates. See the vignette Errata."
  )
  reference <- paste(
    "Pei L, Li R, Zhou H, Du W, Gu Y, Jiang Y, Wang Y, Chen X, Sun J,",
    "Zhu J (2023). A Physiologically Based Pharmacokinetic Approach to",
    "Recommend an Individual Dose of Tacrolimus in Adult Heart Transplant",
    "Recipients. Pharmaceutics 15(11):2580.",
    "doi:10.3390/pharmaceutics15112580.",
    "Physiology (organ weights and perfusions, Table S1) from Levitt DG",
    "(2002) BMC Clin Pharmacol 2:5 and Levitt DG, Schnider TW (2005) BMC",
    "Anesthesiol 5:4; see modellib('Levitt_2005_propofol_pbpk')."
  )
  vignette <- "Pei_2023_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Tissue states hold total tacrolimus amount in the
  # well-stirred organ; arterial / venous hold the amount in blood.
  compartmentData <- list(
    depot     = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    gut       = list(analyte = "tacrolimus", units = "ug", specimen = "tissue", verified = TRUE),
    spleen    = list(analyte = "tacrolimus", units = "ug", specimen = "tissue", verified = TRUE),
    pancreas  = list(analyte = "tacrolimus", units = "ug", specimen = "tissue", verified = TRUE),
    liver     = list(analyte = "tacrolimus", units = "ug", specimen = "tissue", verified = TRUE),
    muscle    = list(analyte = "tacrolimus", units = "ug", specimen = "tissue", verified = TRUE),
    kidney    = list(analyte = "tacrolimus", units = "ug", specimen = "tissue", verified = TRUE),
    brain     = list(analyte = "tacrolimus", units = "ug", specimen = "tissue", verified = TRUE),
    heart     = list(analyte = "tacrolimus", units = "ug", specimen = "tissue", verified = TRUE),
    lung      = list(analyte = "tacrolimus", units = "ug", specimen = "tissue", verified = TRUE),
    skin      = list(analyte = "tacrolimus", units = "ug", specimen = "tissue", verified = TRUE),
    tendon    = list(analyte = "tacrolimus", units = "ug", specimen = "tissue", verified = TRUE),
    other     = list(analyte = "tacrolimus", units = "ug", specimen = "tissue", verified = TRUE),
    adipose   = list(analyte = "tacrolimus", units = "ug", specimen = "tissue", verified = TRUE),
    arterial  = list(analyte = "tacrolimus", units = "ug", specimen = "whole blood", verified = TRUE),
    venous    = list(analyte = "tacrolimus", units = "ug", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Pei 2023 Methods 2.5: 'The volume of each organ, V, was adjusted",
        "in relation to the bodyweight (BW) and the proportion of adipose",
        "tissue [24]', citing Levitt 2002 PKQuest. No explicit scaling",
        "formula is printed, so this model follows the scaling used by the",
        "cited physiology source (and by the sibling extraction",
        "Levitt_2005_propofol_pbpk.R): lean organ weights and their blood",
        "flows scale by lean_mass / 52.5 kg, and the adipose compartment",
        "weight is set directly to WT * BODYFAT_PCT / 100. The reference",
        "70 kg / 25 percent body fat human reproduces Table S1 exactly",
        "(organ weights sum to 70.00 kg). Study cohort median (IQR) weight",
        "67.50 (57.50, 75.00) kg, Table 2."
      ),
      source_name        = "Weight (kg)"
    ),
    BODYFAT_PCT = list(
      description        = "Percent total body fat",
      units              = "% (percent, 0-100)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Pei 2023 Methods 2.5 names 'the proportion of adipose tissue' as",
        "the second determinant of organ volumes, and Table 5 reports its",
        "local sensitivity (0.02 on AUC0-last, 0.17 on Cmax, 0.05 on",
        "Ctrough). The paper does not tabulate body-fat percentage per",
        "subject; the reference value 25 percent is the Levitt PKQuest",
        "standard human implied by Table S1 (adipose 17.5 kg of 70 kg",
        "total). Not measured in the Pei 2023 cohort - supply from the",
        "user's own body-composition data, or leave at the 25 percent",
        "reference."
      ),
      source_name        = "proportion of adipose tissue"
    ),
    HCT = list(
      description        = "Hematocrit",
      units              = "% (volume fraction times 100)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters the blood-to-plasma ratio through Eq 5:",
        "BPR = 1 + Bmax * HCT / (KD * HCTm), where HCTm is the median",
        "hematocrit of the heart transplant population. Pei 2023 Table 2",
        "gives the cohort median (IQR) as 31.98 (29.72, 34.96) percent,",
        "which is the value used for HCTm in ini(). Table 5 ranks",
        "hematocrit as the second most sensitive input after fraction",
        "unbound in plasma (0.91 on AUC0-last, 0.67 on Cmax, 1.06 on",
        "Ctrough). Table 7 stratifies the dosing recommendation on",
        "hematocrit bands 0.2-0.3 and 0.3-0.4 (i.e. 20-30 and 30-40",
        "percent on this column's scale). Time-varying in principle;",
        "baseline in the Pei 2023 analysis."
      ),
      source_name        = "HCT (%)"
    ),
    CYP3A5_EXPR = list(
      description        = "CYP3A5 expresser status (at least one CYP3A5*1 allele)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5*3/*3 nonexpresser)",
      notes              = paste(
        "Pei 2023 Results 3.2: 'we identified those carrying CYP3A5*1 and",
        "CYP3A4*18B with extensive metabolizers (EM), and those carrying",
        "CYP3A5*3/*3 and CYP3A4*1/*1 with poor metabolizers (PM). The",
        "values of EM were set at 1 for FACYP3A5 and FACYP3A4. The values",
        "of PM were set at 0.3 and 0.5 for FACYP3A5 and FACYP3A4,",
        "respectively.' The PBPK model therefore pools *1/*1 and *1/*3",
        "into a single expresser stratum, which is exactly the pooling",
        "CYP3A5_EXPR encodes (the companion popPK model",
        "Pei_2023_tacrolimus.R keeps the three levels apart via",
        "CYP3A5_STAR1_HET / CYP3A5_STAR1_HOM). Genotype counts at rs776746",
        "in the 86 genotyped subjects (Table 2): CC 38 (44.2 percent), CT",
        "45 (52.3 percent), TT 3 (3.5 percent); the C allele is CYP3A5*3",
        "(allele frequency 0.703, consistent with the reported *3",
        "frequency in Chinese Han), so CC = *3/*3 = nonexpresser and",
        "CYP3A5_EXPR = 1 for CT and TT."
      ),
      source_name        = "CYP3A5*3 (rs776746)"
    ),
    SNP_CYP3A4_RS2242480_VAR_COUNT = list(
      description        = "CYP3A4 rs2242480 (CYP3A4*18B) variant-allele count",
      units              = "(count, 0/1/2 alleles per subject)",
      type               = "continuous",
      reference_category = "n/a (0 = CYP3A4*1/*1 wild-type homozygote)",
      notes              = paste(
        "Pei 2023 pools CYP3A4*18B carriers (*18B/*18B and *1/*18B) into",
        "the EM stratum with FACYP3A4 = 1 and assigns FACYP3A4 = 0.5 to",
        "the CYP3A4*1/*1 wild-type homozygotes (Results 3.2). Inside",
        "model() the count is collapsed to a carrier indicator",
        "(VAR_COUNT >= 1) to reproduce that pooling exactly. Genotype",
        "counts at rs2242480 in the 86 genotyped subjects (Table 2): CC 45",
        "(52.3 percent), CT 40 (46.5 percent), TT 1 (1.2 percent). The",
        "Discussion notes linkage disequilibrium between CYP3A5*3 and",
        "CYP3A4*18B in this cohort."
      ),
      source_name        = "CYP3A4*18B (rs2242480)"
    ),
    CONC_VORI_NGML = list(
      description        = "Voriconazole whole-blood concentration driving the reversible CYP3A inhibition term",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = "0 (no voriconazole coadministration; DDIRE = 1, no inhibition)",
      notes              = paste(
        "Supplement Eq 1 defines the reversible-inhibition ratio",
        "DDIRE = 1 / (1 + Cb,vor * fub,vor / KI), where Cb,vor is the",
        "voriconazole whole-blood concentration, fub,vor its unbound",
        "fraction in blood (fup 0.42 / BPR 1.00 = 0.42, Table S6) and",
        "KI = 8.70 ng/mL (Table S6). Supplement Eq 2 applies DDIRE to the",
        "CYP3A-mediated fraction of hepatic clearance only. Pei 2023",
        "generated Cb,vor from a separate voriconazole PBPK model whose",
        "tissue-to-plasma partition coefficients, ka and Fg are NOT",
        "tabulated anywhere on disk (Table S6 gives only MW, pKa, LogP,",
        "fup, BPR, KI and CLint), so the voriconazole disposition is not",
        "reproducible from the source. Rather than substitute partition",
        "coefficients from another paper, this model exposes the",
        "voriconazole whole-blood concentration as a (time-varying) input",
        "column so the published inhibition equation and its constants are",
        "encoded faithfully. Set to 0 for tacrolimus alone. See the",
        "vignette for the Cb,vor implied by the paper's own reported",
        "5.80-fold AUC ratio."
      ),
      source_name        = "Cb,vor"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex, female indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Pei 2023 Discussion: 'the mean tacrolimus concentration was",
        "higher in female patients (3.77 +/- 2.82 ng/mL) compared to male",
        "patients (2.75 +/- 3.31 ng/mL) ... We reduced the number of male",
        "patients and performed parameter fitting in the PBPK model after",
        "balancing the gender ratio. No statistical differences were found",
        "in the fitting results.' Sex is therefore screened but not",
        "retained in the PBPK model; its effect is carried indirectly",
        "through WT and BODYFAT_PCT."
      )
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = paste(
        "A significant covariate on CL/F in the companion popPK model",
        "(Table S4, exponent -0.19) but NOT part of the PBPK model, which",
        "has no hepatic-function scaling on CLint. Recorded here so the",
        "difference between the two modelling approaches - which Pei 2023",
        "discusses explicitly - is visible from the model file."
      )
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 115L,
    n_studies        = 1L,
    n_concentrations = 443L,
    age_range        = "adults >= 18 years; median (IQR) 52.00 (46.00, 61.00) years (Table 2)",
    age_median       = "52.00 years",
    weight_range     = "median (IQR) 67.50 (57.50, 75.00) kg (Table 2)",
    weight_median    = "67.50 kg",
    sex_female_pct   = 19.1,
    disease_state    = paste(
      "Adult heart transplant recipients at Nanjing First Hospital",
      "(November 2012 - January 2023), on tacrolimus (Prograf) with",
      "mycophenolate and corticosteroids. Multi-organ transplant",
      "recipients and patients with incomplete clinical data were",
      "excluded. Data collected from post-transplant day 1 to day 30;",
      "median (IQR) postoperative day 23.00 (19.00, 29.00). Median (IQR)",
      "hematocrit 31.98 (29.72, 34.96) percent, total bilirubin 15.90",
      "(12.13, 21.10) umol/L, albumin 36.74 (34.48, 38.68) g/L."
    ),
    dose_range       = paste(
      "Oral tacrolimus titrated to a whole-blood trough target of 10-15",
      "ng/mL during the early postoperative period; median (IQR) daily",
      "dose 5.00 (4.00, 6.00) mg. The PBPK simulations in the paper use",
      "single 1.0-4.5 mg doses (Table 4), 2.5 mg single dose and 2.5 mg",
      "q12h (Table 6), and 2, 3 and 5 mg single doses in healthy adults",
      "(Figures 2 and 4, Table S5)."
    ),
    regions          = "China (Nanjing First Hospital; external evaluation cohort at Wuhan Union Hospital)",
    genotyping       = paste(
      "20 SNPs genotyped in 86 of the 115 subjects. rs35599367, rs1135840,",
      "rs150461093, rs2229109 and rs4253728 failed Hardy-Weinberg",
      "equilibrium (p < 0.05) and were excluded. CYP3A5*3 (rs776746),",
      "CYP3A4*18B (rs2242480) and IL-10 G-1082A (rs1800896) were the",
      "significant genetic covariates."
    ),
    external_evaluation = paste(
      "100 heart transplant recipients at Wuhan Union Hospital; the four",
      "fitted parameters were re-estimated (Table 3, third column: Ka 1.4",
      "1/h, KD 5.9 ng/mL, Bmax 176.8 ng/mL, CLint 10,256.2 L/h) with no",
      "statistical difference from the model-building estimates."
    ),
    notes            = paste(
      "Software: SimBiology 5.8.2 (MathWorks) for the PBPK model; the",
      "companion popPK model (Pei_2023_tacrolimus.R) used Phoenix NLME",
      "8.3. Whole-blood tacrolimus by CMIA, LLOQ 2 ng/mL, quantitative",
      "range 2-30 ng/mL. Ethics: Nanjing Medical University College Ethics",
      "Committee KY20190404-03-KS-01. The PBPK model was built on the",
      "healthy-adult single-5-mg-dose data of reference [37] and then",
      "re-fitted to the 115 heart transplant recipients; Table 3 reports",
      "all three parameter columns."
    )
  )

  ini({
    # ---- Absorption (Eq 2) -------------------------------------------------
    # Ka for the heart-transplant model-building dataset. Table 3 also
    # reports 4.4 1/h (healthy adults, n = 5) and 1.4 1/h (external
    # evaluation, n = 100); those columns are exercised in the vignette via
    # rxSolve(params = ...) rather than as separate model files.
    lka <- log(1.9); label("First-order absorption rate constant into the gut compartment (1/h)")  # Pei 2023 Table 3, model building dataset (SE 0.2)

    # Fraction of the dose absorbed. Pei 2023 Methods 2.5: "Fg was assumed
    # to be constant (Fg = 0.2) [26]." Assumed, not estimated.
    fg <- fixed(0.2); label("Fraction of the oral dose absorbed into the gut compartment (unitless)")  # Pei 2023 Methods 2.5 (Eq 2), assumed from ref [26]

    # ---- Red-cell binding / blood-to-plasma ratio (Eq 5) -------------------
    lbmax <- log(145.9); label("Red-cell binding capacity Bmax (ng/mL)")           # Pei 2023 Table 3, model building dataset (SE 18.4)
    lkd   <- log(7.2);   label("Red-cell binding affinity constant KD (ng/mL)")    # Pei 2023 Table 3, model building dataset (SE 0.7)

    # Median hematocrit of the heart transplant population, the HCTm
    # normaliser in Eq 5. Not fitted - it is the cohort median.
    hct_median <- fixed(31.98); label("Median hematocrit of the heart transplant population, HCTm (%)")  # Pei 2023 Table 2, HCT median

    # ---- Hepatic elimination (Eqs 3, 4, 6) ---------------------------------
    lclint <- log(11535); label("Hepatic intrinsic clearance CLint (L/h)")  # Pei 2023 Table 3, model building dataset, "11,535" (SE 506.3)

    fup <- fixed(0.013); label("Fraction unbound in plasma (unitless)")  # Pei 2023 Table 1, from ref [17]

    # Fractions of hepatic clearance by each route. Table 1 reports
    # fmCYP3A4 = 0.35 and fmCYP3A5 = 0.55 as "Calculated"; the Discussion
    # states "The CYP3A metabolic fraction fmCYP3A of tacrolimus was 0.9",
    # so the residual non-CYP3A route carries fm_other = 1 - 0.9 = 0.10.
    fm_cyp3a5 <- fixed(0.55); label("Fraction of hepatic clearance mediated by CYP3A5 (unitless)")  # Pei 2023 Table 1
    fm_cyp3a4 <- fixed(0.35); label("Fraction of hepatic clearance mediated by CYP3A4 (unitless)")  # Pei 2023 Table 1
    fm_other  <- fixed(0.10); label("Fraction of hepatic clearance mediated by non-CYP3A routes (unitless)")  # Pei 2023 Discussion, fmCYP3A = 0.9 so fm_other = 1 - 0.9

    # Genotype-specific enzyme activity levels (Eq 6). Pei 2023 Results 3.2:
    # "The values of EM were set at 1 for FACYP3A5 and FACYP3A4. The values
    # of PM were set at 0.3 and 0.5 for FACYP3A5 and FACYP3A4,
    # respectively." Assigned, not estimated.
    fa_cyp3a5_em <- fixed(1.0); label("CYP3A5 activity level in extensive metabolizers, CYP3A5*1 carriers (unitless)")  # Pei 2023 Results 3.2
    fa_cyp3a5_pm <- fixed(0.3); label("CYP3A5 activity level in poor metabolizers, CYP3A5*3/*3 (unitless)")             # Pei 2023 Results 3.2
    fa_cyp3a4_em <- fixed(1.0); label("CYP3A4 activity level in extensive metabolizers, CYP3A4*18B carriers (unitless)")  # Pei 2023 Results 3.2
    fa_cyp3a4_pm <- fixed(0.5); label("CYP3A4 activity level in poor metabolizers, CYP3A4*1/*1 (unitless)")               # Pei 2023 Results 3.2

    # ---- Voriconazole reversible inhibition (supplement Eqs 1-2) ----------
    ki_vori  <- fixed(8.70); label("Voriconazole inhibition constant KI for CYP3A (ng/mL)")            # Pei 2023 Table S6, "Predicted"
    fub_vori <- fixed(0.42); label("Voriconazole unbound fraction in blood, fup 0.42 / BPR 1.00 (unitless)")  # Pei 2023 Table S6 (fup 0.42, BPR 1.00)

    # ---- Tissue distribution ----------------------------------------------
    # Kp scaling factor `a`. Pei 2023 Methods 2.5 introduces it as: "The
    # average KT:p of tacrolimus obtained using Rodgers and Rowland's
    # formulas was 1.3. The average KT:p of mice reported in the literature
    # was 11.9, which was nearly 10-fold higher [30]. Therefore, the human
    # KT:p was applied by a scaling factor a." Results 3.2 then states "The
    # scaling factor a was fitted to be 350."
    #
    # The printed value 350 is REFUTED by the paper's own validation tables:
    # applied as (a * KT:p / BPR) - the only placement Eq 1's nested
    # fraction allows - it under-predicts every Table S5 healthy-adult AUC
    # and Cmax by 4-10 fold and every Table 6 heart-transplant exposure by
    # up to 16 fold. The value encoded here, 9.15 = 11.9 / 1.3, is the
    # mouse-to-human average-Kp ratio that the Methods sentence above
    # motivates directly; it is read off the paper rather than fitted by
    # us. Operator-ratified 2026-08-05 (sidecar oare_PMC10675244 q1 = A).
    # The companion popPK model corroborates the 9-30 range independently:
    # Vd/F 656.8 L at Fg 0.2 implies a true Vd of about 131 L, which needs
    # a of roughly 29 rather than 350. See the vignette Errata.
    kp_scale <- fixed(9.15); label("Kp scaling factor a applied to the Rodgers-Rowland tissue-to-plasma partition coefficients (unitless)")  # Pei 2023 Methods 2.5 (11.9 mouse average / 1.3 human average); printed fitted value 350 refuted - see vignette Errata

    # Rodgers and Rowland tissue-to-plasma partition coefficients, Table S2.
    # The "Additional Organ" row (1.9) covers tendon and other.
    kp_adipose  <- fixed(1.90); label("Adipose tissue-to-plasma partition coefficient (unitless)")   # Pei 2023 Table S2, Adipose
    kp_brain    <- fixed(1.91); label("Brain tissue-to-plasma partition coefficient (unitless)")     # Pei 2023 Table S2, Brain
    kp_gut      <- fixed(1.56); label("Gut tissue-to-plasma partition coefficient (unitless)")       # Pei 2023 Table S2, Gut
    kp_heart    <- fixed(0.74); label("Heart tissue-to-plasma partition coefficient (unitless)")     # Pei 2023 Table S2, Heart
    kp_kidney   <- fixed(0.96); label("Kidney tissue-to-plasma partition coefficient (unitless)")    # Pei 2023 Table S2, Kidney
    kp_liver    <- fixed(1.33); label("Liver tissue-to-plasma partition coefficient (unitless)")     # Pei 2023 Table S2, Liver
    kp_lung     <- fixed(0.52); label("Lung tissue-to-plasma partition coefficient (unitless)")      # Pei 2023 Table S2, Lung
    kp_pancreas <- fixed(1.37); label("Pancreas tissue-to-plasma partition coefficient (unitless)")  # Pei 2023 Table S2, Pancreas
    kp_muscle   <- fixed(0.96); label("Muscle tissue-to-plasma partition coefficient (unitless)")    # Pei 2023 Table S2, Muscle
    kp_skin     <- fixed(1.07); label("Skin tissue-to-plasma partition coefficient (unitless)")      # Pei 2023 Table S2, Skin
    kp_spleen   <- fixed(0.97); label("Spleen tissue-to-plasma partition coefficient (unitless)")    # Pei 2023 Table S2, Spleen
    kp_addl     <- fixed(1.90); label("Tendon and 'other' tissue-to-plasma partition coefficient (unitless)")  # Pei 2023 Table S2, Additional Organ

    # ---- Residual error ----------------------------------------------------
    # Pei 2023 Results 3.2: "the exponential error model with the minimum
    # Akaike information criterion (AIC) and Bayesian information criterion
    # (BIC) was selected as an error model." No residual-error magnitude is
    # reported anywhere in the paper or supplement (Table 3 lists only the
    # four fitted parameters and their standard errors, despite its caption
    # promising "corresponding residual error"). Fixed at 0 rather than
    # invented; see the vignette Errata.
    propSd <- fixed(0); label("Proportional residual error (fraction); magnitude not reported by Pei 2023")  # Pei 2023 Results 3.2 - error model stated (exponential) but no SD reported
  })

  model({
    # ======================================================================
    # 1. Reference physiology - Levitt PKQuest standard human, Table S1.
    #    Weights in kg, perfusions in L/h per kg of organ. The 16 tabulated
    #    organ weights sum to exactly 70.00 kg. Bone (4 kg, perfusion 0) is
    #    omitted as a state: it is absent from the Figure S1 flow diagram
    #    and, with zero blood flow, is kinetically inert. The remaining 13
    #    tissues plus arterial and venous blood are the paper's "15
    #    compartments".
    # ======================================================================
    ref_wt          <- 70      # kg, total body weight of the Table S1 human
    ref_bodyfat_pct <- 25      # %, 17.5 kg adipose of 70 kg
    ref_lean_mass   <- ref_wt * (1 - ref_bodyfat_pct / 100)   # 52.5 kg

    # Reference organ weights (kg) - Table S1 "Weight (kg)"
    w_ref_arterial <- 1.1
    w_ref_venous   <- 4.4
    w_ref_gut      <- 1.17
    w_ref_spleen   <- 0.15
    w_ref_pancreas <- 0.14
    w_ref_liver    <- 1.8
    w_ref_muscle   <- 26
    w_ref_kidney   <- 0.31
    w_ref_brain    <- 1.4
    w_ref_heart    <- 0.33
    w_ref_lung     <- 0.54
    w_ref_skin     <- 2.6
    w_ref_tendon   <- 3
    w_ref_other    <- 5.56
    w_ref_adipose  <- 17.5

    # Organ perfusions (L/h per kg of organ) - Table S1 "Blood Flow (L/h/kg)"
    perf_gut      <- 50
    perf_spleen   <- 78
    perf_pancreas <- 27.9
    perf_liver    <- 15      # hepatic artery; the portal contribution arrives from gut + spleen + pancreas
    perf_muscle   <- 1.4
    perf_kidney   <- 240
    perf_brain    <- 33.6
    perf_heart    <- 48
    perf_skin     <- 6
    perf_tendon   <- 0.6
    perf_other    <- 1.2
    perf_adipose  <- 2.6

    # ======================================================================
    # 2. Per-subject body composition. Pei 2023 Methods 2.5 states only that
    #    organ volumes are "adjusted in relation to the bodyweight (BW) and
    #    the proportion of adipose tissue [24]" (ref [24] = Levitt 2002
    #    PKQuest) without printing a formula, so the scaling of the cited
    #    physiology source is used: lean organ weights scale by lean mass,
    #    the adipose compartment weight is set from body-fat percent, and
    #    each organ's blood flow follows its weight through the Table S1
    #    per-kg perfusion. See the vignette Errata.
    # ======================================================================
    v_adipose  <- WT * BODYFAT_PCT / 100
    lean_scale <- (WT - v_adipose) / ref_lean_mass

    v_arterial <- w_ref_arterial * lean_scale
    v_venous   <- w_ref_venous   * lean_scale
    v_gut      <- w_ref_gut      * lean_scale
    v_spleen   <- w_ref_spleen   * lean_scale
    v_pancreas <- w_ref_pancreas * lean_scale
    v_liver    <- w_ref_liver    * lean_scale
    v_muscle   <- w_ref_muscle   * lean_scale
    v_kidney   <- w_ref_kidney   * lean_scale
    v_brain    <- w_ref_brain    * lean_scale
    v_heart    <- w_ref_heart    * lean_scale
    v_lung     <- w_ref_lung     * lean_scale
    v_skin     <- w_ref_skin     * lean_scale
    v_tendon   <- w_ref_tendon   * lean_scale
    v_other    <- w_ref_other    * lean_scale

    # Blood flows (L/h) = organ weight * per-kg perfusion (Table S1).
    q_gut      <- v_gut      * perf_gut
    q_spleen   <- v_spleen   * perf_spleen
    q_pancreas <- v_pancreas * perf_pancreas
    q_liverart <- v_liver    * perf_liver
    q_muscle   <- v_muscle   * perf_muscle
    q_kidney   <- v_kidney   * perf_kidney
    q_brain    <- v_brain    * perf_brain
    q_heart    <- v_heart    * perf_heart
    q_skin     <- v_skin     * perf_skin
    q_tendon   <- v_tendon   * perf_tendon
    q_other    <- v_other    * perf_other
    q_adipose  <- v_adipose  * perf_adipose

    # Total hepatic blood flow = hepatic artery + portal vein (gut + spleen
    # + pancreas). At the 70 kg reference this is
    # 27 + 58.5 + 11.7 + 3.906 = 101.106 L/h, i.e. 29 percent of the
    # 344.358 L/h cardiac output summed below - a physiologically standard
    # hepatic fraction, which confirms that the Table S1 liver perfusion of
    # 15 L/h/kg is the hepatic-artery flow rather than total liver flow.
    q_livertot <- q_liverart + q_gut + q_spleen + q_pancreas

    # The lung receives the whole cardiac output. Summing every organ's
    # flow at the reference gives 344.358 L/h, within 1.4 percent of the
    # Table S1 lung entry (0.54 kg * 629 L/h/kg = 339.66 L/h); the summed
    # value is used so that mass balance is exact.
    q_lung <-
      q_livertot + q_muscle + q_kidney + q_brain + q_heart +
      q_skin + q_tendon + q_other + q_adipose

    # ======================================================================
    # 3. Blood binding and hepatic clearance (Eqs 3-6)
    # ======================================================================
    # Back-transform the four fitted parameters (all constrained positive).
    bmax  <- exp(lbmax)
    kd    <- exp(lkd)
    clint <- exp(lclint)

    # Eq 5: BPR = 1 + Bmax * HCT / (KD * HCTm)
    bpr <- 1 + bmax * HCT / (kd * hct_median)

    # Eq 4: fu_b = fu_p / BPR
    fub <- fup / bpr

    # Eq 3: hepatic extraction ratio E = fu_b * CLint / (Q_liver + fu_b * CLint)
    e_liver <- fub * clint / (q_livertot + fub * clint)

    # Genotype-specific enzyme activity levels (Pei 2023 Results 3.2). The
    # CYP3A4 variant-allele count is collapsed to a carrier indicator
    # because the paper pools *18B/*18B with *1/*18B.
    cyp3a4_carrier <- (SNP_CYP3A4_RS2242480_VAR_COUNT >= 1)
    fa_cyp3a5 <- fa_cyp3a5_pm + (fa_cyp3a5_em - fa_cyp3a5_pm) * CYP3A5_EXPR
    fa_cyp3a4 <- fa_cyp3a4_pm + (fa_cyp3a4_em - fa_cyp3a4_pm) * cyp3a4_carrier

    # Supplement Eq 1: reversible voriconazole inhibition of CYP3A.
    # CONC_VORI_NGML = 0 gives ddi_re = 1, i.e. no inhibition.
    ddi_re <- 1 / (1 + CONC_VORI_NGML * fub_vori / ki_vori)

    # Eq 6 (with supplement Eq 2 for the DDI): hepatic blood clearance.
    cl_liver <-
      ((fm_cyp3a5 * fa_cyp3a5 + fm_cyp3a4 * fa_cyp3a4) * ddi_re + fm_other) *
      q_livertot * e_liver

    # ODE elimination coefficient (L/h) applied to the liver's effluent
    # blood concentration. For a well-stirred organ with elimination
    # k * C_out the blood clearance is Q * k / (Q + k); solving for k at
    # CL_H = cl_liver gives the expression below, so the ODE reproduces
    # Eq 6's blood clearance exactly. cl_liver < q_livertot always (the
    # Eq 6 bracket is at most 1 and E < 1), so the denominator is positive.
    cl_effluent <- q_livertot * cl_liver / (q_livertot - cl_liver)

    # ======================================================================
    # 4. Tissue-to-blood partition coefficients. Eq 1's nested fraction
    #    KT:p / BPR converts the Rodgers-Rowland tissue-to-PLASMA
    #    coefficients of Table S2 to tissue-to-BLOOD; the scaling factor a
    #    (Methods 2.5) multiplies the human KT:p.
    # ======================================================================
    kpb_gut      <- kp_scale * kp_gut      / bpr
    kpb_spleen   <- kp_scale * kp_spleen   / bpr
    kpb_pancreas <- kp_scale * kp_pancreas / bpr
    kpb_liver    <- kp_scale * kp_liver    / bpr
    kpb_muscle   <- kp_scale * kp_muscle   / bpr
    kpb_kidney   <- kp_scale * kp_kidney   / bpr
    kpb_brain    <- kp_scale * kp_brain    / bpr
    kpb_heart    <- kp_scale * kp_heart    / bpr
    kpb_lung     <- kp_scale * kp_lung     / bpr
    kpb_skin     <- kp_scale * kp_skin     / bpr
    kpb_tendon   <- kp_scale * kp_addl     / bpr
    kpb_other    <- kp_scale * kp_addl     / bpr
    kpb_adipose  <- kp_scale * kp_adipose  / bpr

    # ======================================================================
    # 5. Concentrations. Tissue states hold ug, volumes are L, so every
    #    concentration below is ug/L = ng/mL (the paper's unit). c_*_out is
    #    the blood concentration leaving the organ in equilibrium with the
    #    well-stirred tissue.
    # ======================================================================
    c_arterial <- arterial / v_arterial
    c_venous   <- venous   / v_venous

    c_gut_out      <- gut      / v_gut      / kpb_gut
    c_spleen_out   <- spleen   / v_spleen   / kpb_spleen
    c_pancreas_out <- pancreas / v_pancreas / kpb_pancreas
    c_liver_out    <- liver    / v_liver    / kpb_liver
    c_muscle_out   <- muscle   / v_muscle   / kpb_muscle
    c_kidney_out   <- kidney   / v_kidney   / kpb_kidney
    c_brain_out    <- brain    / v_brain    / kpb_brain
    c_heart_out    <- heart    / v_heart    / kpb_heart
    c_lung_out     <- lung     / v_lung     / kpb_lung
    c_skin_out     <- skin     / v_skin     / kpb_skin
    c_tendon_out   <- tendon   / v_tendon   / kpb_tendon
    c_other_out    <- other    / v_other    / kpb_other
    c_adipose_out  <- adipose  / v_adipose  / kpb_adipose

    ka <- exp(lka)

    # ======================================================================
    # 6. ODE system. Eq 1 for every non-eliminating tissue, Eq 2 for the
    #    gut (which also receives the absorbed dose), and the liver as the
    #    single eliminating organ fed by the hepatic artery and the portal
    #    vein. The depot holds mg (the dosing unit); the factor 1000
    #    converts the absorbed mg/h into the ug/h that the tissue states
    #    use.
    # ======================================================================
    d/dt(depot) <- -ka * depot

    d/dt(gut)      <- q_gut      * (c_arterial - c_gut_out) + fg * ka * depot * 1000
    d/dt(spleen)   <- q_spleen   * (c_arterial - c_spleen_out)
    d/dt(pancreas) <- q_pancreas * (c_arterial - c_pancreas_out)
    d/dt(muscle)   <- q_muscle   * (c_arterial - c_muscle_out)
    d/dt(kidney)   <- q_kidney   * (c_arterial - c_kidney_out)
    d/dt(brain)    <- q_brain    * (c_arterial - c_brain_out)
    d/dt(heart)    <- q_heart    * (c_arterial - c_heart_out)
    d/dt(skin)     <- q_skin     * (c_arterial - c_skin_out)
    d/dt(tendon)   <- q_tendon   * (c_arterial - c_tendon_out)
    d/dt(other)    <- q_other    * (c_arterial - c_other_out)
    d/dt(adipose)  <- q_adipose  * (c_arterial - c_adipose_out)

    d/dt(liver) <-
      q_liverart * c_arterial +
      q_gut      * c_gut_out +
      q_spleen   * c_spleen_out +
      q_pancreas * c_pancreas_out -
      q_livertot * c_liver_out -
      cl_effluent * c_liver_out

    d/dt(lung)     <- q_lung * (c_venous - c_lung_out)
    d/dt(arterial) <- q_lung * c_lung_out - q_lung * c_arterial

    d/dt(venous) <-
      q_muscle   * c_muscle_out +
      q_kidney   * c_kidney_out +
      q_brain    * c_brain_out +
      q_heart    * c_heart_out +
      q_skin     * c_skin_out +
      q_tendon   * c_tendon_out +
      q_other    * c_other_out +
      q_adipose  * c_adipose_out +
      q_livertot * c_liver_out -
      q_lung     * c_venous

    # ======================================================================
    # 7. Observation: venous whole-blood tacrolimus (ng/mL), the matrix and
    #    sampling site of the CMIA assay used throughout Pei 2023.
    # ======================================================================
    Cc <- c_venous
    Cc ~ prop(propSd)
  })
}
