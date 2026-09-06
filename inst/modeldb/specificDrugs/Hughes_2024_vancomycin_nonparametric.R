Hughes_2024_vancomycin_nonparametric <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for vancomycin in hospitalized adults with class",
    "3 obesity (BMI >= 40 kg/m2), developed as a NONPARAMETRIC fit with Pmetrics using the",
    "nonparametric adaptive grid (NPAG) algorithm. Structurally identical to its parametric sibling:",
    "fat-free mass is computed inside the model by the Janmahasatian equations from body weight,",
    "height and sex, creatinine clearance is then computed by a Cockcroft-Gault form that substitutes",
    "fat-free mass for total body weight, clearance scales as (CRCL/100)^1.05, and both volumes scale",
    "linearly on (FFM/70). Unlike the parametric fit, between-subject variability is carried on all",
    "four disposition parameters including intercompartmental clearance, and the residual error is a",
    "fixed Pmetrics assay-error polynomial rather than an estimated proportional term. NPAG estimates",
    "a discrete joint distribution of individual parameters that has no closed form; it is",
    "approximated here by independent lognormal marginals, so the shape of the joint density is not",
    "recoverable from this encoding. The parametric counterpart is",
    "modellib('Hughes_2024_vancomycin_parametric').",
    sep = " "
  )
  reference <- paste(
    "Hughes MSA, Hughes JH, Endicott J, Langton M, Ahern JW, Keizer RJ.",
    "Developing parametric and nonparametric models for model-informed precision dosing:",
    "a quality improvement effort in vancomycin for patients with obesity.",
    "Ther Drug Monit 2024;46(5):575-583. doi:10.1097/FTD.0000000000001214.",
    "Parameter estimates from Table 2, Pmetrics column; structural model from the Table 2 footnote,",
    "Table 3 and Supplemental Digital Content 1 section S5 (Pmetrics model file),",
    "http://links.lww.com/TDM/A753.",
    "The fat-free-mass equations are Janmahasatian S, Duffull SB, Ash S, Ward LC, Byrne NM, Green B.",
    "Quantification of lean bodyweight. Clin Pharmacokinet 2005;44(10):1051-1065.",
    sep = " "
  )
  vignette <- "Hughes_2024_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central     = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Table 1, development dataset: median 134 kg (range 91.6-218). Listed as WEIGHT in the",
        "Supplemental Digital Content 1 S5 Pmetrics #Cov block. Total body weight does NOT enter the",
        "disposition parameters directly -- Results/Nonparametric: 'FFM implemented as an allometric",
        "power model was a better predictor than TBW'. WT enters only through the two internally",
        "derived quantities that do: body mass index (BMI = WT / (HT/100)^2) and, through BMI,",
        "Janmahasatian fat-free mass.",
        sep = " "
      ),
      source_name        = "WEIGHT"
    ),
    HT = list(
      description        = "Body height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Table 1, development dataset: median 170 cm (range 122-190). The Supplemental Digital",
        "Content 1 S4 control stream (the same derived covariates, written out in NONMEM) divides by",
        "100 to obtain metres before squaring (BMI = WT / ((HT/100)**2)), which fixes the unit as cm.",
        "Enters only through BMI and hence fat-free mass. The Pmetrics S5 listing takes FFM as a",
        "precomputed data column rather than deriving it, so HT does not appear in its #Cov block.",
        sep = " "
      ),
      source_name        = "HT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Table 1, development dataset: median 56.3 years (range 24.2-89.3). Used only inside the",
        "fat-free-mass-based Cockcroft-Gault creatinine-clearance calculation, which the Pmetrics S5",
        "listing takes as a precomputed data column (CRCLf) rather than deriving; it is derived",
        "inside the model here so the model is self-contained from raw demographics.",
        sep = " "
      ),
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Table 1, development dataset: 41 male / 42 female (50.6% female). The Supplemental Digital",
        "Content 1 S4 control stream uses the OPPOSITE polarity, SEX with 1 = male: the default",
        "fat-free-mass line uses the male Janmahasatian coefficients (6680, 216) and IF(SEX.EQ.0)",
        "overrides them with the female coefficients (8780, 244), while the Cockcroft-Gault term is",
        "written 0.85**(1-SEX) so that females receive the 0.85 factor. Converted to the canonical",
        "SEXF (1 = female) via SEXF = 1 - SEX, so SEXF selects the female fat-free-mass coefficients",
        "and the Cockcroft-Gault factor becomes 0.85^SEXF.",
        sep = " "
      ),
      source_name        = "SEX"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Table 1, development dataset: median 0.84 mg/dL (range 0.40-2.51); listed as CREAT in the",
        "Supplemental Digital Content 1 S5 Pmetrics #Cov block. Used as the denominator of the",
        "Cockcroft-Gault equation with the 72 constant, which fixes the unit as mg/dL. NOTE: the S5",
        "#Cov block declares both CRCL and a CAPCRCL column, but the #Sec block uses the UNCAPPED",
        "fat-free-mass form CRCLf -- consistent with Results/Nonparametric, 'Capping with CrCl did",
        "not improve the fit'. No cap is applied here.",
        sep = " "
      ),
      source_name        = "CREAT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 83L,
    n_studies      = 1L,
    age_range      = "24.2-89.3 years",
    age_median     = "56.3 years",
    weight_range   = "91.6-218 kg",
    weight_median  = "134 kg",
    sex_female_pct = 50.6,
    race_ethnicity = "Not reported",
    disease_state  = paste(
      "Hospitalized adults with class 3 obesity (body mass index at least 40 kg/m2 at any point",
      "during treatment) receiving intravenous vancomycin under routine model-informed precision",
      "dosing, with at least one vancomycin therapeutic-drug-monitoring concentration. Patient status,",
      "level of care and indication for vancomycin were not collected (stated as a limitation in the",
      "Discussion), so no sepsis / critical-illness stratification is available.",
      sep = " "
    ),
    dose_range     = paste(
      "Intravenous vancomycin per the institutional protocol: 20 mg/kg loading dose (maximum 4000 mg),",
      "then an intermittent-infusion maintenance regimen selected to reach an AUC target of",
      "400-600 mg*h/L",
      sep = " "
    ),
    regions        = "United States (single centre: University of Vermont Medical Center, Burlington, Vermont)",
    bmi_range      = "40-70.3 kg/m2",
    bmi_median     = "46.3 kg/m2",
    renal_function = paste(
      "Serum creatinine median 0.84 mg/dL (range 0.40-2.51); Cockcroft-Gault creatinine clearance",
      "median 174 mL/min (range 30.2-579) on total body weight, or 84.1 mL/min (range 13.0-255) on",
      "fat-free mass",
      sep = " "
    ),
    n_concentrations = 272L,
    notes          = paste(
      "Table 1, development dataset column -- the SAME 83 patients and 272 levels the parametric",
      "sibling was fitted to; the two models differ only in estimation method. Median 2 levels per",
      "patient (range 1-29), of which 31 were peaks (within 2 h after end of infusion), 115 troughs",
      "(within 1 h of the next administration) and 126 random. Data were collected 1 November 2021 to",
      "14 February 2023 from the InsightRX Nova precision-dosing software database. Three levels were",
      "removed during data cleaning. All levels were above the lower limit of quantification. The",
      "Discussion records that no subpopulations could be identified in the nonparametric analysis,",
      "which the authors attribute to routine-care data being too sparse for NPAG to resolve",
      "multimodality.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Structural parameters: Table 2, "Pmetrics" column. NPAG estimates a
    # discrete joint density over the support points rather than a point
    # estimate plus a parametric variance; the tabulated "Estimate" is the
    # summary value of that density and is encoded here as the median of a
    # lognormal marginal (log(x) is the mu of a lognormal whose median is x).
    # The Table 2 "Range" column spans the support points and is NOT an
    # uncertainty interval -- it is recorded in each label below but is not
    # encoded as a variance. See vignette Errata.
    # ------------------------------------------------------------------------
    lcl <- log(6.74); label("Clearance at CRCL = 100 mL/min (L/h); support 3.57-10.2")
    # Table 2, CL row, Pmetrics Estimate = 6.74 L/h (RSE 4.5%, range 3.57-10.2).
    lvc <- log(104); label("Central volume at FFM = 70 kg (L); support 61.5-161.4")
    # Table 2, Vc row, Pmetrics Estimate = 104 L (RSE 5.3%, range 61.5-161.4).
    lq <- log(7.5); label("Intercompartmental clearance Q (L/h); support 2.12-19.2")
    # Table 2, Q row, Pmetrics Estimate = 7.5 L/h (RSE 42%, range 2.12-19.2).
    lvp <- log(142); label("Peripheral volume at FFM = 70 kg (L); support 47.9-277")
    # Table 2, Vp row, Pmetrics Estimate = 142 L (RSE 34%, range 47.9-277).

    # ------------------------------------------------------------------------
    # Covariate effects. Same final structure as the parametric sibling --
    # Results/Nonparametric: "once CrCl was included as a predictor of CL, the
    # inclusion of the allometric exponent on CL and Q did not improve the model
    # and was, therefore, removed, resulting in the same final structure as the
    # parametric models." Confirmed by the Supplemental Digital Content 1 S5
    # #Sec block, which writes Q = Q_0 with no covariate term:
    #   CL = CL_0*(CRCLf/100)**TH1
    #   V  = V_0 * (FFM/70)
    #   Q  = Q_0
    #   V2 = V2_0 * (FFM/70)
    # ------------------------------------------------------------------------
    e_crcl_cl <- 1.05; label("Power exponent on (CRCL/100) for CL (unitless)")
    # Table 2, u3 row ("Effect of CrCl on CL"), Pmetrics Estimate = 1.05
    # (RSE 10%), i.e. TH1 in the S5 #Sec block. CRCLf is the Cockcroft-Gault
    # value computed on FAT-FREE MASS -- Results/Nonparametric: "Using FFM as an
    # input to the Cockcroft-Gault equation improved the fit significantly
    # compared with using TBW". Table 2 reports no BSV for this coefficient
    # ("-"), so it carries no eta even though NPAG places every model parameter
    # in the joint density.
    e_ffm_vc_vp <- fixed(1); label("Allometric exponent on (FFM/70) for Vc and Vp (unitless)")
    # Table 2, u2 row ("Effect of FFM on Vc and Vp"), Pmetrics = 1.0 FIX. A
    # linear exponent, matching the S5 #Sec lines "V = V_0 * (FFM/70)" and
    # "V2 = V2_0 * (FFM/70)", which are written as plain proportionalities.
    # Results/Nonparametric: the exponents "could not be estimated reliably, but
    # when fixed to 0.75 and 1.0 for Q/CL and Vc/Vp, respectively, led to
    # significant improvement in fit" -- the 0.75 CL/Q term was then dropped
    # once CrCl entered, leaving only this fixed linear volume exponent.

    # ------------------------------------------------------------------------
    # Inter-individual variability on ALL FOUR disposition parameters, including
    # Q. Table 3, "Between-subject variability" row: nonparametric =
    # "CL, Vc, Q, Vp". Results/Nonparametric: "The discrete distribution of
    # individual parameters was estimated for all 4 basic PK parameters, even
    # though fixing Q and Vp to a single population value did not result in a
    # worse AIC", and Discussion: "because nonparametric model development
    # workflows do not conventionally assess the statistical relevance of BSV to
    # basic PK parameters, BSV was incorporated for all 4 basic PK parameters in
    # the base model." This is the one structural difference from the parametric
    # sibling, which has no eta on Q.
    #
    # SCALE CONVENTION: unlike the parametric columns (where the paper's own
    # NONMEM OMEGA listing shows the quoted %CV is 100*sqrt(omega^2)), the
    # Pmetrics %CV is a genuine descriptive coefficient of variation of the
    # discrete distribution, computed by Pmetrics as SD/mean over the support
    # points. Matching it with a lognormal marginal therefore requires the exact
    # conversion omega^2 = log(CV^2 + 1), which reproduces each reported %CV
    # exactly. See vignette Errata for the side-by-side of the two conventions.
    #
    # Eta correlations are not reported for the nonparametric fit -- Supplemental
    # Digital Content 1 S2 shows a correlation figure for the posterior
    # individual estimates but prints no coefficients -- so the marginals are
    # encoded as independent. See vignette Errata.
    # ------------------------------------------------------------------------
    etalcl ~ 0.04810803
    # Table 2, CL row, Pmetrics BSV = 22.2 %CV; log(0.222^2 + 1) = 0.04810803.
    etalvc ~ 0.06156912
    # Table 2, Vc row, Pmetrics BSV = 25.2 %CV; log(0.252^2 + 1) = 0.06156912.
    etalq ~ 0.2558818
    # Table 2, Q row, Pmetrics BSV = 54 %CV; log(0.54^2 + 1) = 0.2558818.
    etalvp ~ 0.1155583
    # Table 2, Vp row, Pmetrics BSV = 35.0 %CV; log(0.350^2 + 1) = 0.1155583.

    # ------------------------------------------------------------------------
    # Residual error: the Pmetrics assay-error polynomial, FIXED. Methods:
    # "For the nonparametric model, a combined proportional (10%) and additive
    # (0.5 mg/L) residual error model was implemented using the polynomial assay
    # error, with values fixed to 0.1 (C1) and 0.5 (C0), respectively, and a
    # lambda of 0.1". Table 2 gives RUV = 0.1 FIX and RUVadd = 0.5 FIX; Table 3
    # gives "10% + 0.5 mg/L (fixed), L = 0.1"; and the Supplemental Digital
    # Content 1 S5 #Err block is
    #   L=0.1
    #   0.5,0.1,0,0
    # i.e. lambda = 0.1 with polynomial coefficients C0 = 0.5, C1 = 0.1,
    # C2 = C3 = 0. All four sources agree exactly.
    #
    # Pmetrics weights each observation by its assay SD = C0 + C1*C, a LINEAR
    # SUM in the observed concentration. nlmixr2's default add()+prop()
    # combines in quadrature instead, so combined1() is required to reproduce
    # the linear-sum form. Both coefficients are stated assay constants held
    # fixed during estimation, hence fixed().
    # ------------------------------------------------------------------------
    addSd <- fixed(0.5); label("Assay-error polynomial intercept C0 (mg/L)")
    # Table 2, RUVadd row, Pmetrics = 0.5 FIX; S5 #Err polynomial first
    # coefficient.
    propSd <- fixed(0.1); label("Assay-error polynomial slope C1 (fraction)")
    # Table 2, RUV row, Pmetrics = 0.1 FIX; S5 #Err polynomial second
    # coefficient.
  })
  model({
    # 1. Derived covariate terms. The Supplemental Digital Content 1 S5 Pmetrics
    #    listing consumes FFM and CRCLf as precomputed data columns (its #Cov
    #    block declares WEIGHT, FFM, CREAT, CRCL, CAPCRCL, CRCLf); they are
    #    derived here from raw demographics so the model is self-contained,
    #    using the same equations the paper writes out in Data Collection and in
    #    the S4 $PK block for the parametric fit:
    #      BMI  = WT / ((HT/100)**2)
    #      FFM  = 9270 * WT / (6680 + 216*BMI)      [male,   SEX = 1]
    #      FFM  = 9270 * WT / (8780 + 244*BMI)      [female, SEX = 0]
    #      CRCL = (140-AGE) * FFM * 0.85**(1-SEX) / (72*CR)
    #    with SEXF = 1 - SEX, so SEXF selects the female branch and the
    #    Cockcroft-Gault sex factor becomes 0.85^SEXF. CRCLf is the uncapped
    #    fat-free-mass form.
    bmi_i    <- WT / (HT / 100)^2
    ffm_male <- 9270 * WT / (6680 + 216 * bmi_i)
    ffm_fem  <- 9270 * WT / (8780 + 244 * bmi_i)
    ffm_i    <- ffm_male + SEXF * (ffm_fem - ffm_male)
    crcl_i   <- (140 - AGE) * ffm_i * 0.85^SEXF / (72 * CREAT)

    # 2. Individual PK parameters. CL scales on renal function only and the two
    #    volumes on fat-free mass only; Q carries no covariate but, unlike the
    #    parametric sibling, does carry an eta.
    cl <- exp(lcl + etalcl) * (crcl_i / 100)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (ffm_i / 70)^e_ffm_vc_vp
    vp <- exp(lvp + etalvp) * (ffm_i / 70)^e_ffm_vc_vp
    q  <- exp(lq + etalq)

    # 3. Micro-constants. The S5 #Sec block defines Ke = CL/V explicitly and
    #    relies on the Pmetrics analytic two-compartment library for the
    #    distribution rate constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. Standard two-compartment IV disposition, written out
    #    explicitly for house consistency with the parametric sibling.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # 5. Observation and error. The S5 #Out block is Y(1) = X(1)/V, so dose in
    #    mg over volume in L gives mg/L. combined1() reproduces the Pmetrics
    #    assay SD = C0 + C1*C linear sum rather than nlmixr2's default
    #    quadrature combination.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
