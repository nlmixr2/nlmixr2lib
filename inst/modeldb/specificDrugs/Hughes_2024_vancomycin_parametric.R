Hughes_2024_vancomycin_parametric <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for vancomycin in hospitalized adults with class",
    "3 obesity (BMI >= 40 kg/m2), developed as a parametric (maximum-likelihood) fit with nlmixr2",
    "SAEM. Fat-free mass is computed inside the model by the Janmahasatian equations from body weight,",
    "height and sex; creatinine clearance is then computed by a Cockcroft-Gault form that substitutes",
    "fat-free mass for total body weight. Clearance scales as (CRCL/100)^0.887 with no allometric term",
    "of its own; central and peripheral volumes scale linearly on (FFM/70). Intercompartmental",
    "clearance carries neither a covariate nor an eta. This is the parametric half of a",
    "parametric-versus-nonparametric comparison fitted to a single institution's routine",
    "therapeutic-drug-monitoring data; the nonparametric counterpart is",
    "modellib('Hughes_2024_vancomycin_nonparametric').",
    sep = " "
  )
  reference <- paste(
    "Hughes MSA, Hughes JH, Endicott J, Langton M, Ahern JW, Keizer RJ.",
    "Developing parametric and nonparametric models for model-informed precision dosing:",
    "a quality improvement effort in vancomycin for patients with obesity.",
    "Ther Drug Monit 2024;46(5):575-583. doi:10.1097/FTD.0000000000001214.",
    "Parameter estimates from Table 2, nlmixr2 (SAEM) column; structural model from the Table 2",
    "footnote and from Supplemental Digital Content 1 sections S3 (nlmixr2 model code) and S4",
    "(NONMEM control stream), http://links.lww.com/TDM/A753.",
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
        "Table 1, development dataset: median 134 kg (range 91.6-218). Total body weight does NOT",
        "enter the disposition parameters directly -- the paper's central finding is that fat-free",
        "mass outperformed total body weight as the size descriptor in this obese cohort. WT enters",
        "only through the two internally derived quantities that do: body mass index",
        "(BMI = WT / (HT/100)^2) and, through BMI, Janmahasatian fat-free mass.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    HT = list(
      description        = "Body height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Table 1, development dataset: median 170 cm (range 122-190). The Supplemental Digital",
        "Content 1 S4 control stream divides by 100 to obtain metres before squaring",
        "(BMI = WT / ((HT/100)**2)), which fixes the unit as cm. Enters only through BMI and hence",
        "fat-free mass.",
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
        "internal fat-free-mass-based Cockcroft-Gault creatinine-clearance calculation.",
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
        "Table 1, development dataset: median 0.84 mg/dL (range 0.40-2.51). Used as the denominator",
        "of the Cockcroft-Gault equation with the 72 constant, which fixes the unit as mg/dL. NOTE:",
        "no cap is applied to the resulting creatinine clearance -- the Results state that capping",
        "CrCl at 150 or 200 mL/min worsened the fit, so the (CRCL/100)^0.887 power term is applied to",
        "the raw computed value. Median CrCl computed this way is 84.1 mL/min (Table 1, 'Creatinine",
        "clearance (based on FFM)'), against 174 mL/min when total body weight is used instead.",
        sep = " "
      ),
      source_name        = "CR"
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
      "Table 1, development dataset column. 83 patients contributed 272 therapeutic-drug-monitoring",
      "levels (median 2 per patient, range 1-29), of which 31 were peaks (within 2 h after end of",
      "infusion), 115 troughs (within 1 h of the next administration) and 126 random. Data were",
      "collected 1 November 2021 to 14 February 2023 from the InsightRX Nova precision-dosing",
      "software database. Three levels were removed during data cleaning. All levels were above the",
      "lower limit of quantification. A separate external validation dataset of 576 patients across",
      "74 US organizations (1191 levels) was used to assess predictive performance but not to",
      "estimate any parameter.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Structural parameters: Table 2, "nlmixr2 (SAEM)" column. The paper fitted
    # the same model in NONMEM as a software cross-check and reports that column
    # too; the NONMEM parameter set is separately available in the library as
    # modellib('Tong_2026_vancomycin_hughes'). Results/Parametric states "we
    # present only the results for nlmixr2", so the nlmixr2 column is the
    # paper's own headline parametric model and is what is encoded here.
    # ------------------------------------------------------------------------
    lcl <- log(5.09); label("Clearance at CRCL = 100 mL/min (L/h)")
    # Table 2, CL row, nlmixr2 (SAEM) Estimate = 5.09 L/h (RSE 4.1%).
    lvc <- log(64.9); label("Central volume at FFM = 70 kg (L)")
    # Table 2, Vc row, nlmixr2 (SAEM) Estimate = 64.9 L (RSE 8.4%).
    lq <- log(6.36); label("Intercompartmental clearance Q (L/h)")
    # Table 2, Q row, nlmixr2 (SAEM) Estimate = 6.36 L/h (RSE 18%).
    lvp <- log(66.4); label("Peripheral volume at FFM = 70 kg (L)")
    # Table 2, Vp row, nlmixr2 (SAEM) Estimate = 66.4 L (RSE 13%).

    # ------------------------------------------------------------------------
    # Covariate effects. Table 2 footnote gives the general form
    #   CLi = CL*(FFM/70)^u1*(CRCL/100)^u3,  Vi  = V*(FFM/70)^u2,
    #   Qi  = Q *(FFM/70)^u1,                Vp,i = Vp*(FFM/70)^u2
    # but u1 is absent from the Table 2 parameter rows because it was dropped
    # from the final model: Results/Parametric states that "the inclusion of the
    # allometric exponent on clearance (CL) and intercompartmental clearance (Q)
    # considerably worsened the model fit when fixed and was estimated to be
    # close to zero when estimated; therefore, this exponent was excluded", and
    # Table 3 lists the final parametric covariates as CL~CrCl, Vc~FFM, Vp~FFM
    # only. Confirmed by the Supplemental Digital Content 1 S4 control stream,
    # whose $PK writes TVQ = THETA(3) with no covariate term. So u1 = 0 here:
    # CL carries only the renal term and Q carries no covariate at all.
    # ------------------------------------------------------------------------
    e_crcl_cl <- 0.887; label("Power exponent on (CRCL/100) for CL (unitless)")
    # Table 2, u3 row ("Effect of CrCl on CL"), nlmixr2 (SAEM) Estimate = 0.887
    # (RSE 6.1%). CRCL is the Cockcroft-Gault value computed on FAT-FREE MASS,
    # not total body weight -- Results/Parametric: "using FFM as input to the
    # Cockcroft-Gault equation instead of TBW improved the fit considerably".
    e_ffm_vc_vp <- fixed(1); label("Allometric exponent on (FFM/70) for Vc and Vp (unitless)")
    # Table 2, u2 row ("Effect of FFM on Vc and Vp"), nlmixr2 (SAEM) = 1.0 FIX.
    # Results/Parametric: "Estimation of the exponent of Vc/Vp also did not
    # provide a better fit than fixing to 1.0." A linear exponent, matching the
    # Supplemental Digital Content 1 S3 nlmixr2 code line
    # "tFFM_V <- fixed(1.0)".

    # ------------------------------------------------------------------------
    # Inter-individual variability on CL, Vc and Vp; none on Q. Table 3,
    # "Between-subject variability" row: parametric = "CL, Vc, Vp", and
    # Discussion: "The final version of the parametric model included BSV in
    # only 3 parameters because including it in the peripheral clearance (Q) did
    # not produce a statistically significantly better fit." Table 2's Q row
    # correspondingly shows "-" in the nlmixr2 BSV column.
    #
    # SCALE CONVENTION: the paper reports BSV as "%CV", but its own NONMEM
    # listings show that the quoted number is 100*sqrt(omega^2) rather than the
    # exact lognormal 100*sqrt(exp(omega^2)-1). The Table 2 NONMEM column
    # (21.8 / 18.0 / 79.7 %CV) is reproduced by the same authors' published
    # OMEGA BLOCK(3) diagonal (0.0473154, 0.0322236, 0.634656) as
    # sqrt -> 21.75 / 17.95 / 79.67, an exact three-parameter match, whereas the
    # lognormal conversion would give 22.0 / 18.1 / 94.1. Variances below are
    # therefore (%CV/100)^2. See vignette Errata.
    #
    # The Supplemental Digital Content 1 S3 nlmixr2 code fits these three etas
    # as a correlated BLOCK(3); the final off-diagonal estimates are not
    # published in either the paper or the supplement, so the etas are encoded
    # as independent here rather than inventing a correlation. See vignette
    # Errata.
    # ------------------------------------------------------------------------
    etalcl ~ 0.062001
    # Table 2, CL row, nlmixr2 BSV = 24.9 %CV (RSE_BSV 13%); 0.249^2 = 0.062001.
    etalvc ~ 0.031684
    # Table 2, Vc row, nlmixr2 BSV = 17.8 %CV (RSE_BSV 21%); 0.178^2 = 0.031684.
    etalvp ~ 0.644809
    # Table 2, Vp row, nlmixr2 BSV = 80.3 %CV (RSE_BSV 32%); 0.803^2 = 0.644809.

    # ------------------------------------------------------------------------
    # Residual error. Structure comes from the supplement code, which is
    # unambiguous: Supplemental Digital Content 1 S3 (nlmixr2) declares
    # "prop_sd <- 0.2" (estimated) and "add_sd <- fixed(1e-3)", and S4 (NONMEM)
    # builds W = SQRT(IPRED**2 * PROP**2 + ADD**2) with $THETA "0.001 ; add
    # error". That is a proportional model with a nominal, fixed, effectively
    # nil additive term for numerical stability -- consistent with Table 2,
    # whose "RUVadd / Additive error" row is "-" for both parametric columns.
    # nlmixr2's add()+prop() default combines in quadrature, matching the
    # NONMEM W above.
    #
    # ERRATUM: Table 3 instead reports the parametric residual error as
    # "15.6% + 1.2 mg/L". That 1.2 mg/L additive term appears nowhere else --
    # not in Table 2, not in S3, not in S4 -- and the S4 listing's own
    # proportional value (0.157812) does not match Table 2's nlmixr2 column
    # either. Table 2 is the paper's dedicated parameter-estimate table and its
    # nlmixr2 column is the model encoded here, so Table 2 is used. See vignette
    # Errata for the full three-way comparison.
    # ------------------------------------------------------------------------
    propSd <- 0.168; label("Proportional residual error (fraction)")
    # Table 2, RUV row ("Proportional error"), nlmixr2 (SAEM) = 16.8 (percent).
    addSd <- fixed(0.001); label("Additive residual error (mg/L), nominal value for numerical stability")
    # Supplemental Digital Content 1 S3, "add_sd <- fixed(1e-3)"; equivalently
    # S4 $THETA(6) "0.001 ; add error".
  })
  model({
    # 1. Derived covariate terms, transcribed from the Supplemental Digital
    #    Content 1 S4 $PK block. Body mass index, then Janmahasatian fat-free
    #    mass with sex-specific coefficients, then a Cockcroft-Gault creatinine
    #    clearance that substitutes fat-free mass for total body weight:
    #      BMI  = WT / ((HT/100)**2)
    #      FFM  = 9270 * WT / (6680 + 216*BMI)      [male,   SEX = 1]
    #      FFM  = 9270 * WT / (8780 + 244*BMI)      [female, SEX = 0]
    #      CRCL = (140-AGE) * FFM * 0.85**(1-SEX) / (72*CR)
    #    with SEXF = 1 - SEX, so SEXF selects the female branch and the
    #    Cockcroft-Gault sex factor becomes 0.85^SEXF. The FFM equations are
    #    also written out in the paper's Data Collection section.
    bmi_i    <- WT / (HT / 100)^2
    ffm_male <- 9270 * WT / (6680 + 216 * bmi_i)
    ffm_fem  <- 9270 * WT / (8780 + 244 * bmi_i)
    ffm_i    <- ffm_male + SEXF * (ffm_fem - ffm_male)
    crcl_i   <- (140 - AGE) * ffm_i * 0.85^SEXF / (72 * CREAT)

    # 2. Individual PK parameters. CL scales on renal function only, the two
    #    volumes on fat-free mass only, and Q carries neither a covariate nor
    #    an eta.
    cl <- exp(lcl + etalcl) * (crcl_i / 100)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (ffm_i / 70)^e_ffm_vc_vp
    vp <- exp(lvp + etalvp) * (ffm_i / 70)^e_ffm_vc_vp
    q  <- exp(lq)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. The source used the closed-form $SUBROUTINE ADVAN3 TRANS4
    #    (two-compartment, IV) in NONMEM and linCmt() in nlmixr2; the standard
    #    two-compartment system is written out explicitly here for house
    #    consistency with the sibling nonparametric model.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # 5. Observation and error. S1 = V1 in the control stream, so dose in mg
    #    over volume in L gives mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
